"""Binding Graph Model and Utils"""
import copy
import random
import sys
import datetime
import functools
import os
import logging
import re
import tempfile
from collections import defaultdict
import itertools
import types
import uuid

import networkx as nx
from xmlbuilder import XMLBuilder
from pbcommand.models import ResourceTypes
from pbcommand.pb_io.common import write_pipeline_chunks

from pbsmrtpipe.exceptions import (TaskIdNotFound, MalformedBindingError,
                                   InvalidEntryPointError,
                                   MalformedBindingStrError,
                                   TaskExecutionError, ChunkGatheringError,
                                   TaskChunkingError)
from pbsmrtpipe.opts_graph import resolve_di
from pbsmrtpipe.exceptions import (MalformedBindingGraphError,
                                   BindingFileTypeIncompatiblyError)
from pbsmrtpipe.models import MetaScatterTask, TaskStates
import pbsmrtpipe.constants as GlobalConstants
from pbsmrtpipe.pb_io import strip_entry_prefix, binding_str_to_task_id_and_instance_id
from pbsmrtpipe.utils import validate_type_or_raise, nfs_refresh
from pbsmrtpipe.graph.models import (TaskBindingNode,
                                     ConstantsNodes,
                                     BindingInFileNode,
                                     BindingOutFileNode,
                                     _TaskLike,
                                     EntryPointNode,
                                     EntryOutBindingFileNode,
                                     TaskScatterBindingNode,
                                     VALID_FILE_NODE_CLASSES,
                                     VALID_TASK_NODE_CLASSES,
                                     TaskChunkedBindingNode,
                                     BindingChunkInFileNode,
                                     BindingChunkOutFileNode,
                                     TaskGatherBindingNode,
                                     VALID_ALL_TASK_NODE_CLASSES)

log = logging.getLogger(__name__)
slog = logging.getLogger(GlobalConstants.SLOG_PREFIX + "__name__")
#logging.basicConfig(level=logging.DEBUG)


class BindingsGraph(nx.DiGraph):

    # This is the new model. This will replace the Abstract Graph

    def _validate_type(self, n):
        _allowed_types = tuple(itertools.chain(VALID_TASK_NODE_CLASSES, VALID_FILE_NODE_CLASSES))

        if not isinstance(n, _allowed_types):
            msg = "Got type {t} for {k}. Allowed types {a}.".format(a=_allowed_types, t=type(n), k=self.__class__.__name__)
            log.error(msg)
            raise TypeError(msg)
        return n

    def add_file_in_to_out(self, in_node, out_node):
        validate_type_or_raise(in_node, BindingInFileNode)
        validate_type_or_raise(out_node, BindingOutFileNode)
        self.add_edge(in_node, out_node)

    def add_edge(self, u, v, attr_dict=None, **attr):
        for n in (u, v):
            self._validate_type(n)
        super(BindingsGraph, self).add_edge(u, v, attr_dict=attr_dict, **attr)

    def add_node(self, n, attr_dict=None, **attr):
        self._validate_type(n)
        super(BindingsGraph, self).add_node(n, attr_dict=attr_dict, **attr)

    def _get_nodes_by_klasses(self, klasses, data=False):
        return [n for n in list(self.nodes_iter(data=data)) if isinstance(n, klasses)]

    def _get_sorted_nodes_by_klass(self, klasses):
        nodes = nx.topological_sort(self)
        return [n for n in nodes if isinstance(n, klasses)]

    def _get_next_instance_id(self, meta_task):
        xd = defaultdict(lambda : 0)
        for tnode in self.all_task_type_nodes():
            if not isinstance(tnode, EntryPointNode):
                i = max(tnode.instance_id, xd[tnode.meta_task.task_id])
                xd[tnode.meta_task.task_id] = i
        return xd[meta_task.task_id] + 1

    def add_meta_task(self, meta_task):
        """Generate a TaskBindingNode and assign an instance-id

        Returns the TaskBindingNode instance
        """
        # After graph initialization, every access point should be
        # done here
        instance_id = self._get_next_instance_id(meta_task)
        t = TaskBindingNode(meta_task, instance_id)
        add_node_by_type(self, t)
        return t

    def add_chunked_meta_task(self, meta_task, chunk_id, chunk_group_id, operator_id):
        instance_id = self._get_next_instance_id(meta_task)
        t = TaskChunkedBindingNode(meta_task, instance_id, chunk_id, chunk_group_id, operator_id)
        add_node_by_type(self, t)
        return t

    def add_scatter_meta_task(self, meta_task, original_nid, original_task_type_id, chunk_group_id):
        instance_id = self._get_next_instance_id(meta_task)
        t = TaskScatterBindingNode(meta_task, original_nid, original_task_type_id, instance_id, chunk_group_id)
        add_node_by_type(self, t)
        return t

    def add_gather_meta_task(self, meta_task, chunk_key):
        instance_id = self._get_next_instance_id(meta_task)
        t = TaskGatherBindingNode(meta_task, instance_id, chunk_key)
        add_node_by_type(self, t)
        return t

    def _get_next_file_instance_id(self, file_node_class, file_type):
        xd = defaultdict(lambda : 0)
        for fnode in self.file_nodes():
            if not isinstance(fnode, file_node_class):
                i = max(fnode.instance_id, xd[file_type.file_type_id])
                xd[file_type.file_type_id] = i
        return xd[file_type.file_type_id] + 1

    def add_binding_in(self, meta_task, index, file_type):
        i = self._get_next_file_instance_id((BindingInFileNode, BindingChunkInFileNode), file_type)
        n = BindingInFileNode(meta_task, i, index, file_type)
        add_node_by_type(self, n)
        return n

    def add_binding_out(self, meta_task, index, file_type_id):
        i = self._get_next_file_instance_id(BindingOutFileNode, file_type_id)
        n = BindingOutFileNode(meta_task, i, index, file_type_id)
        add_node_by_type(self, n)
        return n

    def add_binding_file_chunk_in(self, meta_task, index, file_type, chunk_id, chunk_group_id):
        i = self._get_next_file_instance_id(
            (BindingInFileNode, BindingChunkInFileNode), file_type)
        n = BindingChunkInFileNode(meta_task, i, index, file_type, chunk_id,
                                   chunk_group_id)
        add_node_by_type(self, n)
        return n

    def add_binding_file_chunk_out(self, meta_task, index, file_type, chunk_id, chunk_group_id):
        i = self._get_next_file_instance_id(
            (BindingOutFileNode, BindingChunkOutFileNode), file_type)
        n = BindingChunkOutFileNode(meta_task, i, index, file_type, chunk_id,
                                    chunk_group_id)
        add_node_by_type(self, n)
        return n

    def task_nodes(self, data=False):
        """
        This always returns a t-sorted list of task nodes
        """
        # FIXME this only returns non-chunked tasks types
        return self._get_sorted_nodes_by_klass(VALID_TASK_NODE_CLASSES)

    def all_task_type_nodes(self):
        """Returns Task of very type"""
        return self._get_sorted_nodes_by_klass(_TaskLike)

    def chunked_task_nodes(self):
        return self._get_nodes_by_klasses((TaskChunkedBindingNode, ))

    def scattered_task_nodes(self):
        return self._get_nodes_by_klasses((TaskScatterBindingNode, ))

    def gathered_task_nodes(self):
        return self._get_nodes_by_klasses((TaskGatherBindingNode, ))

    def file_nodes(self, data=False):
        return self._get_sorted_nodes_by_klass(VALID_FILE_NODE_CLASSES)

    def entry_binding_nodes(self, data=False):
        # these are files-esque
        nodes = nx.topological_sort(self)
        return [n for n in nodes if isinstance(n, EntryOutBindingFileNode)]

    def entry_point_nodes(self, data=False):
        # these are task-esque
        nodes = nx.topological_sort(self)
        return [n for n in nodes if isinstance(n, EntryPointNode)]

    def get_tasks_by_state(self, state):
        return get_tasks_by_state(self, state)

    def is_workflow_complete(self):
        return _is_workflow_complete(self)

    def __repr__(self):
        d = dict(k=self.__class__.__name__,
                 n=len(self.nodes()),
                 t=len(self.all_task_type_nodes()),
                 f=len(self.file_nodes()),
                 e=len(self.edges()),
                 p=len(self.entry_point_nodes()))
        return "<{k} Tasks:{t} Files:{f} EntryPoints:{p} node:{n} edges:{e} >".format(**d)


def validate_binding_graph_integrity(bg):
    """
    Check for malformed graphs with dangling input file nodes.


    :raises: MalformedBindingGraphError
    :param bg: Binding Graph

    :type bg: BindingsGraph
    :return:
    """
    for n in bg.all_task_type_nodes():
        for i in bg.predecessors(n):
            # the in degree should be 1,
            # or 0 if the node is an Entry Point
            i_d, o_d = bg.in_degree(i), bg.out_degree(i)
            d = o_d - i_d
            # print i_d, o_d, d, n, i

            emsg = "Invalid In-degree ({x}) of task id {d} with file node {f}.".format(x=i_d, d=n, f=i)

            if i_d == 0:
                if not isinstance(i, EntryPointNode):
                    raise MalformedBindingGraphError(emsg)
            elif i_d == 1:
                # this is the expected
                pass
            else:
                # Gather Tasks are allowed to have an in-degree > 1
                if not isinstance(n, TaskGatherBindingNode):
                    raise MalformedBindingGraphError(emsg)

    return True


def validate_compatible_binding_file_types(bg):
    """
    Validate that File bindings are of the same Type. This should warn if
    FileTypes.Fasta -> FileTypes.GFF are incompatible.

    :param bg:
    :type bg: BindingsGraph
    :return:
    """

    for n in bg.all_task_type_nodes():
        if isinstance(n, TaskBindingNode):
            meta_task = n.meta_task
            input_types = meta_task.input_types

            for file_binding_node in bg.predecessors(n):
                i = file_binding_node.index

                expected_type = input_types[i]
                # should be try because of how the graph is built, but
                # adding case for safety/debugging.
                if expected_type == file_binding_node.file_klass:
                    for x in bg.predecessors(file_binding_node):
                        if expected_type != x.file_klass:
                            msg = "Binding Type Incompatibly for task {t}. Expected type {e}, got type {g}".format(t=n, e=expected_type, g=x.file_klass)
                            log.error(msg)
                            raise BindingFileTypeIncompatiblyError(msg)
                else:
                    msg = "Binding Type Incompatibly for task {t}. Expected type {e}, got type {g}".format(t=n, e=expected_type, g=file_binding_node.file_klass)
                    log.error(msg)
                    raise BindingFileTypeIncompatiblyError(msg)

    return True


def add_node_by_type(g, node):
    """
    Adds a node to the graph and initializes the node "properties" defined
    in NODE_ATTRS class var of the node type.
    """
    kw = {}
    for k, v in node.NODE_ATTRS.iteritems():
        value = v() if isinstance(v, types.FunctionType) else v
        kw[k] = value

    g.add_node(node, **kw)
    return g


def _validate_binding_format(b):
    """
    Validates that a task binding is well-formed.

    Two cases

    - taskA -> taskB
    - entry -> taskA

    """
    rxs = (GlobalConstants.RX_BINDING_TASK_ADVANCED,
           GlobalConstants.RX_BINDING_TASK,
           GlobalConstants.RX_BINDING_ENTRY)

    for rx in rxs:
        if rx.search(b) is not None:
            return True

    _d = zip(("advanced task", "task", "entry"), rxs)
    msg = ", ".join(["'{a} pattern' {p}".format(a=a, p=x.pattern) for a, x in _d])
    raise MalformedBindingError("Binding '{b}' expected to match {m}".format(b=b, m=msg))


def _has_validate_binding_format(b):
    try:
        _validate_binding_format(b)
        return True
    except MalformedBindingError:
        return False


def _to_objs_from_binding_str(registered_tasks_d, b_out, b_in):
    """
    Convert bindings to Task Node and File Node objs

    This needs to have more validation and raise better exceptions
    """

    for x in (b_out, b_in):
        _validate_binding_format(x)

    def _log_error(msg):
        emsg = "Failed to process bindings {o} -> {i}. {e}".format(o=b_out, i=b_in, e=msg)
        log.error(emsg)
        return emsg

    def _get_index_or_raise(desc, task_id_, in_out_types_, index_):
        emsg = "MetaTask '{x}' Invalid index. Unable to get index {i} from max index {n} of {d}_types. {t}".format(i=index_, t=in_out_types_, x=task_id_, d=desc, n=(len(in_out_types_) - 1))

        if index_ < len(in_out_types_):
            return list(in_out_types_)[index_]

        raise MalformedBindingStrError(emsg)

    def _get_input_index_or_raise(task_id_, input_types_, index_):
        return _get_index_or_raise("input", task_id_, input_types_, index_)

    def _get_output_index_or_raise(task_id_, output_types_, index_):
        return _get_index_or_raise("output", task_id_, output_types_, index_)

    def _get_meta_task_or_raise(task_id_):
        meta_task_ = registered_tasks_d.get(task_id_, None)
        if meta_task_ is None:
            raise TaskIdNotFound("Unable to find task id '{i}'".format(i=task_id_))

        return meta_task_

    try:
        # Task Id, Instance id, in out position index
        ti_id, ti_in, ti_index = binding_str_to_task_id_and_instance_id(b_in)

        log.debug("Binding parsed binding '{}' => task id:{} instance:{} index:{}".format(b_in, ti_id, ti_in, ti_index))

        # meta task instance
        ti_meta_task = _get_meta_task_or_raise(ti_id)

        if ti_meta_task is None:
            raise TaskIdNotFound("task {i}".format(i=ti_meta_task))

        # file_type_ids = list(ti_meta_task.input_types)
        # file_type_instance = file_type_ids[ti_index]

        file_type_instance = _get_input_index_or_raise(ti_meta_task.task_id, ti_meta_task.input_types, ti_index)

        ti_node = TaskBindingNode(ti_meta_task, ti_in)
        fo_node = BindingInFileNode(ti_meta_task, ti_in, ti_index, file_type_instance)

        # in_node can be an entry point of file
        if b_out.startswith(GlobalConstants.ENTRY_PREFIX):
            eid = b_out.split(GlobalConstants.ENTRY_PREFIX)[-1]

            file_types = list(ti_meta_task.input_types)
            file_type_instance = file_types[ti_index]
            # use the file klass from the input as the type
            in_node = EntryPointNode(eid, file_type_instance)
            fi_node = EntryOutBindingFileNode(eid, file_type_instance)
        else:
            # Out Task
            to_id, to_in, to_i = binding_str_to_task_id_and_instance_id(b_out)

            to_meta_task = _get_meta_task_or_raise(to_id)

            #fi_file_type_instance = to_meta_task.output_types[to_i]
            fi_file_type_instance = _get_output_index_or_raise(to_meta_task.task_id, to_meta_task.output_types, to_i)

            in_node = TaskBindingNode(to_meta_task, instance_id=to_in)
            fi_node = BindingOutFileNode(to_meta_task, to_in, to_i, fi_file_type_instance)
    except (IndexError, TaskIdNotFound) as e:
        _log_error(e.message)
        raise
    except Exception as e:
        _log_error(e.message)
        raise

    result = (in_node, fi_node, ti_node, fo_node)
    log.debug("Binding parsed Raw {o} -> {i} successfully to {r}".format(o=b_out, i=b_in, r=result))
    # Task/EP node, File out, TaskIn, File in/None
    return result


def add_nodes_to_binding_graph(g, to_node, fo_node, ti_node, fi_node, add_node_func):
    """
    Add nodes to graph

    to_node can be a EntryPoint node
    fo_node can be None

    """

    ns = [to_node, fo_node, ti_node, fi_node]

    nodes = [n for n in ns if n is not None]
    tnodes = [n for n in ns if isinstance(n, TaskBindingNode)]

    for node in nodes:
        add_node_func(g, node)

    # binding file self, meta_task, instance_id, index, in_or_out, file_klass
    in_f = lambda a, b: g.add_edge(a, b)
    out_f = lambda a, b: g.add_edge(b, a)

    # add inputs and outputs file nodes that might not be defined in the bindings
    for node in tnodes:
        for binding_klass, add_edge_func in [(BindingInFileNode, in_f), (BindingOutFileNode, out_f)]:
            for i, in_out_file_types in enumerate(getattr(node.meta_task, binding_klass.ATTR_NAME)):
                fnode = binding_klass(node.meta_task, node.instance_id, i, in_out_file_types)
                add_node_func(g, fnode)
                add_edge_func(fnode, node)

    # add edges
    g.add_edge(to_node, fo_node)

    if fi_node is not None:
        g.add_edge(fo_node, fi_node)
        g.add_edge(fi_node, ti_node)


def binding_strs_to_binding_graph(registered_tasks, bindings):
    """Convert a list of binding tuples to BindingGraph

    This should be the ONLY way a Binding Graph is created

    A List of binding strings

    [("pbsmrtpipe.tasks.dev_01", "pbsmrtpipe.tasks.dev_02"), ("$entry:e_id", "pbsmrtpipe.tasks.dev_04")]

    :param registered_tasks: A dict of Registered MetaTasks
    :param bindings: A list of binding strings

    :type registered_tasks: {task_id:MetaTask}

    :rtype: BindingsGraph
    """
    objs = [_to_objs_from_binding_str(registered_tasks, b_in, b_out) for b_in, b_out in {x for x in bindings}]

    bg = BindingsGraph()
    _ = [add_nodes_to_binding_graph(bg, a, b, c, d, add_node_by_type) for a, b, c, d in objs]

    initialize_file_node_attrs(bg)
    initialize_task_node_attrs(bg)

    validate_binding_graph_integrity(bg)

    return bg


def binding_strs_to_xml(bindings):
    import pbsmrtpipe.pb_io as IO

    eps = []
    bs = []
    for e, b in bindings:
        if e.startswith(GlobalConstants.ENTRY_PREFIX):
            eps.append((e, b))
        else:
            bs.append((e, b))

    b = XMLBuilder(IO.Constants.WORKFLOW_ROOT)

    # Have to do this odd getattr to get around how the builder API works
    with getattr(b, 'entry-points'):
        for ep, x in eps:
            getattr(b, 'entry-point')(id=ep, out=x)

    with b.bindings:
        for out_b, in_b in bs:
            d = {'out': out_b, 'in': in_b}
            b.binding(**d)

    return b


def get_next_task_instance_id(g, task_node_klasses, meta_task_id):
    i = 0
    for node in g.nodes():
        if isinstance(g, task_node_klasses):
            if node.meta_task.task_id == meta_task_id:
                if node.instance_id > i:
                    i = node.instance_id
    return i + 1


def _get_next_in_out_file_instance_id(g, file_node_klass, file_type_id):
    i = 0
    for node in g.nodes():
        if isinstance(node, file_node_klass):
            if node.file_klass.file_type_id == file_type_id:
                if node.instance_id > i:
                    i = node.instance_id
    return i + 1


def get_next_in_file_instance_id(g, file_type_id):
    return _get_next_in_out_file_instance_id(g, BindingInFileNode, file_type_id)


def get_next_out_file_instance_id(g, file_type_type_id):
    return _get_next_in_out_file_instance_id(g, BindingOutFileNode, file_type_type_id)


def initialize_file_node_attrs(g):
    """
    Add is_resolved and path to file nodes

    """
    default_attrs = [(ConstantsNodes.FILE_ATTR_IS_RESOLVED, False),
                     (ConstantsNodes.FILE_ATTR_PATH, None)]
    for attr_name, value in default_attrs:
        for n in g.nodes():
            if isinstance(n, VALID_FILE_NODE_CLASSES):
                g.node[n][attr_name] = value


def initialize_task_node_attrs(g):
    # add stderr, stdout, log path
    default_attrs = [(ConstantsNodes.TASK_ATTR_STATE, 'created'),
                     (ConstantsNodes.TASK_ATTR_ROPTS, {}),
                     (ConstantsNodes.TASK_ATTR_NPROC, 1),
                     (ConstantsNodes.TASK_ATTR_CMDS, []),
                     (ConstantsNodes.TASK_ATTR_RUN_TIME, None),
                     (ConstantsNodes.TASK_ATTR_EMESSAGE, None),
                     (ConstantsNodes.TASK_ATTR_IS_CHUNKABLE, False)]

    for attr_name, value in default_attrs:
        for n in g.all_task_type_nodes():
            g.node[n][attr_name] = value


def get_node_attributes(g, name):
    # for backward compatibility with networkx version in SMRTAnalysis
    values = {}
    for n in g.nodes():
        try:
            v = g.node[n][name]
            values[n] = v
        except KeyError:
            pass
    return values


def _is_workflow_complete(g):
    """
    For the workflow to be complete, the task states must all be completed
    and all the files are resolved

    1. All the task nodes must be in the 'finished' state

    :rtype: bool
    """

    for task_node in g.all_task_type_nodes():
        states = get_node_attributes(g, ConstantsNodes.TASK_ATTR_STATE)
        state = states[task_node]
        if state not in TaskStates.COMPLETED_STATES():
            return False

    states = get_node_attributes(g, 'is_resolved')
    all_resolved = all(states.values())
    if not all_resolved:
        return False

    # made it here all the files are resolved and the tasks are all in
    # a completed state
    return True


def was_task_successful(bg, task_like_node):
    if not isinstance(task_like_node, _TaskLike):
        raise TypeError("Not a TaskLike instance {c}".format(c=task_like_node))
    return bg.node[task_like_node][ConstantsNodes.TASK_ATTR_STATE] == TaskStates.SUCCESSFUL


def was_task_successful_with_resolve_outputs(bg, task_like_node):
    if not was_task_successful(bg, task_like_node):
        return False
    return all(bg.node[output_file_node][ConstantsNodes.FILE_ATTR_IS_RESOLVED] for output_file_node in bg.successors(task_like_node))


def was_workflow_successful(bg):
    return all(was_task_successful(bg, t) for t in bg.all_task_type_nodes())


def _are_all_inputs_resolved(bg, tnode):
    _ns = {}
    for fnode in bg.predecessors(tnode):
        if isinstance(fnode, VALID_FILE_NODE_CLASSES):
            is_resolved = bg.node[fnode][ConstantsNodes.FILE_ATTR_IS_RESOLVED]
            _ns[fnode] = is_resolved

    # must have at least one value and all are resolved
    if _ns:
        return all(_ns.values())

    return False


def get_next_runnable_task(g):

    if g.is_workflow_complete():
        return None

    # this should probably do a top sort, then return
    for tnode in g.all_task_type_nodes():
        if isinstance(tnode, TaskBindingNode):
            state = g.node[tnode][ConstantsNodes.TASK_ATTR_STATE]
            is_chunkable = g.node[tnode][ConstantsNodes.TASK_ATTR_IS_CHUNKABLE]
            if state == TaskStates.SCATTERED:
                # these tasks are labeled as on-hold and will be deleted
                # once the gather step is successful
                continue
            elif is_chunkable is True:
                # Skip original 'unchunked' tasks.
                continue
            elif state in TaskStates.RUNNABLE_STATES():
                if _are_all_inputs_resolved(g, tnode):
                    return tnode
            else:
                # log.debug("Skipping chunkable tasks {n}".format(n=tnode))
                pass

    # log.debug("Unable to find runnable task")
    return None


def has_task_in_states(g, task_states):
    # All tasks are running or completed
    return any((g.node[t][ConstantsNodes.TASK_ATTR_STATE] not in task_states for t in g.all_task_type_nodes()))


def are_all_tasks_running(g):
    return all((g.node[t][ConstantsNodes.TASK_ATTR_STATE] not in TaskStates.RUNNABLE_STATES() for t in g.all_task_type_nodes()))


def has_running_task(g):
    return any(g.node[x][ConstantsNodes.TASK_ATTR_STATE] == TaskStates.RUNNING for x in g.all_task_type_nodes())


def has_next_runnable_task(g):
    """
    If there is a valid runnable task

    All tasks are running

    :type g: BindingsGraph
    :rtype: Boolean
    """

    if g.is_workflow_complete():
        return False

    # All tasks are running or completed
    if all((g.node[t][ConstantsNodes.TASK_ATTR_STATE] not in TaskStates.RUNNABLE_STATES() for t in g.all_task_type_nodes())):
        return False

    for tnode in get_tasks_by_state(g, TaskStates.RUNNABLE_STATES()):
        _ns = {}
        for fnode in g.predecessors(tnode):
            if isinstance(fnode, VALID_FILE_NODE_CLASSES):
                is_resolved = g.node[fnode][ConstantsNodes.FILE_ATTR_IS_RESOLVED]
                if is_resolved:
                    _ns[fnode] = is_resolved
                else:
                    # try to find task that is running
                    _ns[fnode] = is_resolved

        if all(_ns.values()):
            return True

    return False


def get_task_input_files(g, tnode):
    """

    :type tnode: TaskBindingNode
    """
    ninput_files = len(tnode.meta_task.input_types)

    # [(index, path)]
    files = {}
    for fnode in g.predecessors(tnode):
        path = g.node[fnode][ConstantsNodes.FILE_ATTR_PATH]
        files[fnode.index] = path

    if ninput_files != len(files):
        log.error("Expected {n} files. Got {f}".format(n=ninput_files, f=files))
        raise ValueError("Expected inputs to have {n}. Got {x}".format(n=ninput_files, x=len(files)))

    return [p for _, p in sorted([(k, v) for k, v in files.iteritems()])]


def get_task_output_files(g, tnode):

    # [(index, path)]
    files = []
    for fnode in g.successors(tnode):
        path = g.node[fnode][ConstantsNodes.FILE_ATTR_PATH]
        files.append((fnode.index, path))

    return [p for _, p in sorted(files)]


def get_tasks_by_state(g, state_or_states):
    if isinstance(state_or_states, (list, type)):
        states = state_or_states
    else:
        states = [state_or_states]

    node_states = {}
    for n in g.nodes():
        if isinstance(n, VALID_ALL_TASK_NODE_CLASSES):
            state = g.node[n][ConstantsNodes.TASK_ATTR_STATE]
            if state in states:
                node_states[n] = state

    return node_states


def update_task_state(g, tnode, state):
    if state not in TaskStates.ALL_STATES():
        raise ValueError("Invalid task state '{s}'".format(s=state))
    g.node[tnode][ConstantsNodes.TASK_ATTR_STATE] = state
    return g


def update_task_state_with_runtime(g, tnode, state, run_time):
    update_task_state(g, tnode, state)
    g.node[tnode][ConstantsNodes.TASK_ATTR_RUN_TIME] = run_time
    return g


def update_task_state_to_failed(g, tnode, run_time, error_message):
    update_task_state_with_runtime(g, tnode, TaskStates.FAILED, run_time)
    g.node[tnode][ConstantsNodes.TASK_ATTR_EMESSAGE] = error_message
    return g


def validate_outputs_and_update_task_to_success(g, tnode, run_time, output_files):
    """

    :param g:
    :param tnode:
    :param run_time:
    :return:

    :type tnode: TaskBindingNode
    """
    slog.debug("Validating task {t} output files {o}".format(o=output_files, t=tnode.idx))
    for output_file in output_files:
        nfs_refresh(output_file)
        if os.path.exists(output_file):
            log.debug("Successfully validated output of task {t} -> {p}".format(t=tnode.idx, p=output_file))
        else:
            e_msg = "Task {n} Failed to find required OUTPUT file '{c}'".format(c=output_file, n=tnode)
            raise TaskExecutionError(e_msg)

    # if we got here everything is fine
    return update_task_state_to_success(g, tnode, run_time)


def update_task_state_to_success(g, tnode, run_time):
    update_task_state(g, tnode, TaskStates.SUCCESSFUL)
    g.node[tnode][ConstantsNodes.TASK_ATTR_RUN_TIME] = run_time
    return True


def update_file_state_to_resolved(g, file_node, path):
    # this should do a node type check
    if not isinstance(file_node, VALID_FILE_NODE_CLASSES):
        raise TypeError("Unable to update state on {f}".format(f=file_node))

    if g.node[file_node][ConstantsNodes.FILE_ATTR_PATH] is None:
        g.node[file_node][ConstantsNodes.FILE_ATTR_PATH] = path
        g.node[file_node][ConstantsNodes.FILE_ATTR_RESOLVED_AT] = datetime.datetime.now()
        g.node[file_node][ConstantsNodes.FILE_ATTR_IS_RESOLVED] = True

    return True


def update_or_set_node_attrs(g, attrs_tuple, nodes):
    for attr_name, value in attrs_tuple:
        for node in nodes:
            g.node[node][attr_name] = value


def update_task_output_file_nodes(bg, tnode, task):
    """


    :type bg: BindingsGraph
    :type tnode: TaskBindingNode
    :type task: pbsmrtpipe.pb_tasks.core.Task

    :param bg:
    :param tnode:
    :param task:
    :return:
    """

    for fnode in bg.successors(tnode):
        i = fnode.index
        update_file_state_to_resolved(bg, fnode, task.output_files[i])

    return True


def resolve_successor_binding_file_path(g):
    """update linked bound files


    :type g: BindingsGraph
    """

    for fnode in g.file_nodes():
        attrs = g.node[fnode]
        is_resolved = attrs.get(ConstantsNodes.FILE_ATTR_IS_RESOLVED, False)
        path = attrs.get(ConstantsNodes.FILE_ATTR_PATH, None)

        if is_resolved and path is None:
            log.debug("Incompatible attrs. Resolved files, must have path defined. File {f}".format(f=fnode))

        if is_resolved and path is not None:
            snodes = g.successors(fnode)
            # log.debug("Updating {n} nodes".format(n=len(snodes)))
            for s in snodes:
                if isinstance(s, (BindingInFileNode, BindingOutFileNode)):
                    update_file_state_to_resolved(g, s, path)

    return True


def resolve_entry_point(g, entry_id, path):
    """
    Update the path and state of path of an entry point based on entry_id.

    An EntryPoint Node is like a task.
    """
    # FIXME
    eid = strip_entry_prefix(entry_id)

    eps = g.entry_point_nodes()

    ep_ids = [e.idx for e in eps]
    ep_ids.sort()

    # It might be better to raise an Exception?
    if eid not in ep_ids:
        raise InvalidEntryPointError("Unable to resolve required entry point id '{i}'. Valid entry point ids {e}".format(i=eid, e=ep_ids))

    attrs_t = [('style', 'solid'), ('path', path), ('is_resolved', True)]
    for ep_node in eps:
        if ep_node.idx == eid:
            update_or_set_node_attrs(g, attrs_t, [ep_node])

            # It's like a task
            update_or_set_node_attrs(g, [(ConstantsNodes.TASK_ATTR_STATE, TaskStates.SUCCESSFUL)], [ep_node])

            # EntryBindingNodes
            fnodes = g.successors(ep_node)
            log.info("{e} updating entry point successors {i}".format(i=fnodes, e=ep_node))
            update_or_set_node_attrs(g, attrs_t, fnodes)

            for xnode in fnodes:
                xnodes = g.successors(xnode)
                update_or_set_node_attrs(g, attrs_t, xnodes)

    return True


def resolve_entry_points(g, ep_d):
    for entry_id, path in ep_d.iteritems():
        # Allowing a bit of slop here. The "$entry:X" can be optionally given
        eid = strip_entry_prefix(entry_id)
        resolve_entry_point(g, eid, path)
        resolve_successor_binding_file_path(g)


def resolve_entry_binding_points(g):
    for n in g.nodes():
        # Task-esque node
        if isinstance(n, EntryOutBindingFileNode):
            # this is pretty awkward
            g.node[n][ConstantsNodes.TASK_ATTR_STATE] = TaskStates.SUCCESSFUL
            g.node[n][ConstantsNodes.TASK_ATTR_RUN_TIME] = 1.0
            g.node[n][ConstantsNodes.FILE_ATTR_IS_RESOLVED] = True
        # File-esque node
        elif isinstance(n, EntryPointNode):
            g.node[n][ConstantsNodes.FILE_ATTR_IS_RESOLVED] = True
            g.node[n][ConstantsNodes.TASK_ATTR_RUN_TIME] = 1.0


def add_scatter_task(g, scatterable_task_node, scatter_meta_task):
    """
    Add a scatter task to the graph, re-map the inputs of the chunkable task
    to the new Scatter Task


    :type g: BindingsGraph
    :param scatterable_task_node: Original task to scatter
    :type scatterable_task_node: TaskBindingNode

    :param scatter_meta_task: Companion scatter-able task, this must have the same input signature as the original task
    :type scatter_meta_task: pbsmrtpipe.models.MetaTask


    F1 -> T1 -> F2

    Get's mapped to

    F1 -> T1 -> F2
    F1 -> T1* -> FC

    Where T1 is the companion Chunk task and emits a Chunk JSON file (FC)
    """
    validate_type_or_raise(scatterable_task_node, TaskBindingNode)
    validate_type_or_raise(scatter_meta_task, MetaScatterTask)

    chunk_group_id = g.node[scatterable_task_node][ConstantsNodes.TASK_ATTR_CHUNK_GROUP_ID]
    if chunk_group_id is None:
        raise ValueError("Invalid chunk group id for node {n}".format(n=scatterable_task_node))

    # this will generate the chunk.json file
    scatter_task_node = g.add_scatter_meta_task(scatter_meta_task, scatterable_task_node.nid, scatterable_task_node.idx, chunk_group_id)
    slog.debug("Added scattered task {t} to graph.".format(t=scatter_task_node))

    # Create new InFile nodes for the Scattered Task and Re-map the inputs
    # of the original task to the inputs of the scatter task
    # note, this does not create new File Binding nodes.
    for input_node in g.predecessors(scatterable_task_node):
        log.debug(("Adding edge ", input_node, scatter_task_node))
        c_in_node = g.add_binding_in(input_node.meta_task, input_node.index, input_node.file_klass)

        g.add_edge(c_in_node, scatter_task_node)
        # duplicate input bindings from
        for i_node in g.predecessors(input_node):
            g.add_edge(i_node, c_in_node)

        g.add_edge(input_node, scatter_task_node)

    # Add the Chunk.json Output file to the Graph
    for positional_index, output_file_type in enumerate(scatter_meta_task.output_types):
        out_file_node = g.add_binding_out(scatter_meta_task, positional_index, output_file_type)
        g.add_edge(scatter_task_node, out_file_node)

    # update the scatter task-node
    attrs = [(ConstantsNodes.TASK_ATTR_WAS_CHUNKED, False),
             (ConstantsNodes.TASK_ATTR_CHUNK_GROUP_ID, chunk_group_id)]
    update_or_set_node_attrs(g, attrs, [scatter_task_node])

    resolve_successor_binding_file_path(g)
    validate_binding_graph_integrity(g)
    return g


def _get_scatterable_task_id(chunk_operators_d, task_id):
    """Get the companion scatterable task id from the original meta task id"""
    ids = []
    for operator_id, chunk_operator in chunk_operators_d.iteritems():
        if task_id == chunk_operator.scatter.task_id:
            ids.append(chunk_operator.scatter.scatter_task_id)

    return ids if ids else None


def _is_task_id_scatterable(chunk_operators_d, task_id):
    return _get_scatterable_task_id(chunk_operators_d, task_id) is not None


def apply_scatterable(bg, chunk_operators_d, task_registry_d):
    """Add All companion scatterable tasks to task from Chunk Operator

    For Operators that have the same scatter-task-id **and** chunk-keys, these
    are considered identical operations.

    The Chunking will be:

    a -> Tx -> b
    a -> Ty -> c

    To

    a -> TS -> chunk.json -> {Tx_i ... Tx_n} -> chunk.gather.json -> Txg -> b
                          -> {Ty_i ... Ty_n} -> chunk.gather.json -> Tyg -> c

    A single scatter task will be applied and the same mappings
    of the chunk keys in the chunk.json will be used to map
    the keys to the corresponding inputs.

    If the chunk keys are different, then a NEW scatter task will be created.

    """

    def is_task_id_scatterable(task_id_):
        return _is_task_id_scatterable(chunk_operators_d, task_id_)

    def get_scatter_task_ids_from_task_id(task_id_):
        return _get_scatterable_task_id(chunk_operators_d, task_id_)

    for tnode in bg.all_task_type_nodes():
        if not isinstance(tnode, (TaskChunkedBindingNode, EntryPointNode)):

            if is_task_id_scatterable(tnode.meta_task.task_id):
                scatterable_task_ids = get_scatter_task_ids_from_task_id(tnode.meta_task.task_id)

                was_chunked = bg.node[tnode][ConstantsNodes.TASK_ATTR_WAS_CHUNKED]
                # Only process tasks that have not been 'scattered' to create the companion task
                if not was_chunked and scatterable_task_ids is not None:
                    for scatterable_task_id in scatterable_task_ids:
                        log.debug("Resolved scatter task {i} from task {x}".format(i=scatterable_task_id, x=tnode.meta_task.task_id))
                        # add scatterable task to generate the chunk.json file
                        # and assign the chunk group id
                        add_scatter_task(bg, tnode, task_registry_d[scatterable_task_id])
                        # mark original Task so the operator is not applied again
                        bg.node[tnode][ConstantsNodes.TASK_ATTR_WAS_CHUNKED] = True

    resolve_successor_binding_file_path(bg)
    validate_binding_graph_integrity(bg)
    return bg


def add_chunkable_task_nodes_to_bgraph(bg, scatter_task_node, pipeline_chunks, chunk_operators, registered_tasks_d):
    """

    Add N TaskChunkedBindingNode(s) to graph and maps the inputs from PipelineChunks


    This must update the scattered-task to WAS_CHUNKED = True at the end
    so the process is only applied once.

    :type bg: BindingsGraph
    :type scatter_task_node: TaskScatterBindingNode
    :type pipeline_chunks: list[PipelineChunk]
    :type chunk_operator: list[ChunkOperator]

    :return:
    """
    validate_type_or_raise(scatter_task_node, TaskScatterBindingNode)

    if bg.node[scatter_task_node][ConstantsNodes.TASK_ATTR_WAS_CHUNKED]:
        return []

    _to_i = lambda x: int(x.split(":")[-1])

    # Original Task Type To Scatter
    task_type_to_scatter = registered_tasks_d[scatter_task_node.original_task_id]
    slog.debug("Chunking by {m}".format(m=scatter_task_node.meta_task))

    # Chunked file from the scatter task, Scattered Task should have
    # only one output of type FileTypes.CHUNK
    chunk_file_node = bg.successors(scatter_task_node)[0]

    chunk_group_id = scatter_task_node.chunk_group_id

    total_chunked_nodes = []
    for chunk_operator in chunk_operators:
        slog.debug("Starting to chunk task type {i} with chunk-group {g} for operator {o}".format(i=task_type_to_scatter.task_id, g=chunk_group_id, o=chunk_operator.idx))

        # {chunk_key -> task in index}
        scatter_ckey_in_index_d = {c.chunk_key: _to_i(c.task_input) for c in chunk_operator.scatter.chunks}

        chunked_task_nodes = []
        for chunk_number, pipeline_chunk in enumerate(pipeline_chunks):
            chunked_task_node = bg.add_chunked_meta_task(task_type_to_scatter, pipeline_chunk.chunk_id, str(chunk_group_id), chunk_operator.idx)

            chunked_task_nodes.append(chunked_task_node)

            for chunk_key, in_index in scatter_ckey_in_index_d.iteritems():

                if chunk_key not in pipeline_chunk.chunk_d:
                    raise KeyError("Unable to find required chunk key '{i}' in chunk {c}. Chunk keys found {k}.".format(i=chunk_key, k=pipeline_chunk.chunk_keys, c=pipeline_chunk))

                # create new InBindingFile(s) from chunk keys
                file_type_instance = task_type_to_scatter.input_types[in_index]

                datum = pipeline_chunk.chunk_d[chunk_key]
                slog.debug("Mapping chunk key {k} -> index {i} {f} with datum {x}".format(k=chunk_key, i=in_index, f=file_type_instance, x=datum))

                in_node = bg.add_binding_file_chunk_in(task_type_to_scatter, in_index, file_type_instance, pipeline_chunk.chunk_id, chunk_group_id)

                # Update the state to resolved
                bg.node[in_node][ConstantsNodes.FILE_ATTR_PATH] = datum
                bg.node[in_node][ConstantsNodes.FILE_ATTR_RESOLVED_AT] = datetime.datetime.now()
                bg.node[in_node][ConstantsNodes.FILE_ATTR_IS_RESOLVED] = True

                bg.add_edge(chunk_file_node, in_node)
                bg.add_edge(in_node, chunked_task_node)

            # Create new outputs of the chunked tasks
            out_nodes = []
            for out_index, out_file_type in enumerate(task_type_to_scatter.output_types):
                out_node = bg.add_binding_file_chunk_out(task_type_to_scatter, out_index, out_file_type, pipeline_chunk.chunk_id, chunk_group_id)
                bg.add_edge(chunked_task_node, out_node)
                out_nodes.append(out_node)

        # If NO chunked tasks were added, there was a serious problem
        if not chunked_task_nodes:
            raise TaskChunkingError("Starting to chunk task-id {i} with chunk-group {g}".format(i=task_type_to_scatter.task_id, g=chunk_group_id))
        else:
            total_chunked_nodes.extend(chunked_task_nodes)

    update_or_set_node_attrs(bg, [(ConstantsNodes.TASK_ATTR_WAS_CHUNKED, True)], [scatter_task_node])

    # log.debug(to_binding_graph_summary(bg))
    slog.info("Chunked Tasks added {n} from task-id {i}".format(n=len(total_chunked_nodes), i=task_type_to_scatter.task_id))
    resolve_successor_binding_file_path(bg)
    validate_binding_graph_integrity(bg)
    return total_chunked_nodes


def label_chunkable_tasks(g, operators_d):
    """
    Adds is_chunkable to Task nodes that have a companion Chunkable task in
    Chunk Operators

    :type operators_d: dict[str, pbsmrtpipe.models.ChunkOperator]
    :type g: BindingsGraph

    :return:
    """
    found_chunkable_task = False

    for task_node in g.task_nodes():
        if isinstance(task_node, TaskBindingNode):

            task_id = task_node.meta_task.task_id

            for operator_id, chunk_operator in operators_d.iteritems():
                if chunk_operator.scatter.task_id == task_id:
                    scatter_task_id = chunk_operator.scatter.scatter_task_id

                    # This Id should be used for ChunkTaskBinding Nodes
                    chunk_group_id = uuid.uuid4()

                    slog.info("Found chunkable task '{i}' assigning chunk-group {s}".format(i=task_id, s=str(chunk_group_id)))

                    # Setting the is-chunkable is trigger that a chunk/scatter task
                    # can be created
                    g.node[task_node][ConstantsNodes.TASK_ATTR_IS_CHUNKABLE] = True
                    # Set the operator id
                    g.node[task_node][ConstantsNodes.TASK_ATTR_OPERATOR_ID] = chunk_operator.idx
                    # Keep track of the original task_id for book-keeping
                    g.node[task_node][ConstantsNodes.TASK_ATTR_COMPANION_CHUNK_TASK_TYPE_ID] = scatter_task_id
                    # This is not great. In the TaskChunkedBindingNode this is a attribute of the node
                    g.node[task_node][ConstantsNodes.TASK_ATTR_CHUNK_GROUP_ID] = str(chunk_group_id)

                    found_chunkable_task = True

    if not found_chunkable_task:
        slog.warn("Unable to find any chunkable tasks from {n} chunk operators.".format(n=len(operators_d)))

    resolve_successor_binding_file_path(g)
    validate_binding_graph_integrity(g)
    return g


def _get_chunk_operators_by_scatter_task_id(scatter_task_id, chunk_operators_d):
    ops = []
    for operator_id, chunk_operator in chunk_operators_d.iteritems():
        if scatter_task_id == chunk_operator.scatter.scatter_task_id:
            ops.append(chunk_operator)

    if ops:
        return ops

    raise KeyError("Unable to find chunk operator for scatter task id {i}".format(i=scatter_task_id))


def apply_chunk_operator(bg, chunk_operators_d, registered_tasks_d, max_nchunks):
    """
    Look for all successfully completed Tasks that were chunked (to
    generate chunk.json), then add corresponding TaskChunkBindingNode(s)
    from the PipelineChunks.

    :param bg:
    :type bg: BindingsGraph
    :param chunk_operators_d:
    :param registered_tasks_d:

    """

    import pbsmrtpipe.pb_io as IO

    # Add chunkabled tasks if necessary
    for tnode_ in bg.scattered_task_nodes():
        if bg.node[tnode_][ConstantsNodes.TASK_ATTR_STATE] == TaskStates.SUCCESSFUL:
            was_chunked = bg.node[tnode_][ConstantsNodes.TASK_ATTR_WAS_CHUNKED]

            # Companion Chunked Task was successful and created chunk.json
            if not was_chunked:
                pipeline_chunks = IO.load_pipeline_chunks_from_json(bg.node[tnode_]['task'].output_files[0])

                if len(pipeline_chunks) > max_nchunks:
                    raise TaskChunkingError("Task {i} created too many {n} chunks. Max chunks must be <= {m}".format(n=len(pipeline_chunks), m=max_nchunks, i=tnode_))
                if not pipeline_chunks:
                    raise TaskChunkingError("Task {i} generated 0 chunks. Chunks must be >0 and <= {m} file: {f}".format(i=tnode_, m=max_nchunks, f=bg.node[tnode_]['task'].output_files[0]))

                # This applies the chunk operator once (or more) if task is chunked by multiple
                # This needs to be revisited fixed.
                chunk_operators = _get_chunk_operators_by_scatter_task_id(bg.node[tnode_]['task'].task_id, chunk_operators_d)
                chunked_nodes = add_chunkable_task_nodes_to_bgraph(bg, tnode_, pipeline_chunks, chunk_operators, registered_tasks_d)
                if chunked_nodes:
                    bg.node[tnode_][ConstantsNodes.TASK_ATTR_WAS_CHUNKED] = True
                    bg.node[tnode_][ConstantsNodes.TASK_ATTR_OPERATOR_ID] = [chunk_operator.idx for chunk_operator in chunk_operators]
                    # bg.node[tnode_][ConstantsNodes.TASK_ATTR_CHUNK_KEYS] = [c.chunk_key for c in chunk_operator.scatter.chunks]
                    slog.info("Successfully applying chunk operator {i} to {x} to generate {n} tasks.".format(x=tnode_, i=[c.idx for c in chunk_operators], n=len(chunked_nodes)))
                    slog.info(chunked_nodes)
            else:
                log.debug("Skipping {t}. Node was already chunked".format(t=tnode_))

    resolve_successor_binding_file_path(bg)
    validate_binding_graph_integrity(bg)
    return bg


def get_companion_unscattered_task_node(bg, chunk_group_id):
    for tnode in bg.task_nodes():
        if tnode.__class__ == TaskBindingNode:
            if bg.node[tnode][ConstantsNodes.TASK_ATTR_CHUNK_GROUP_ID] == chunk_group_id:
                return tnode
    return KeyError("Unable to find scattered companion task for chunk-group {g}".format(g=chunk_group_id))


def add_gather_to_completed_task_chunks(bg, chunk_operators_d, registered_tasks_d, tasks_root_dir):
    """Create the gathered.chunk.json by gathering the Chunked Task Instances.

    1. Find all scattered task nodes
    2. Find all the chunked tasks by chunk-group-id
    3. Load operator
    4. load pipeline gathered chunks from JSON
    5. load

    gather tasks to chunked instances of tasks nodes if the all the chunked tasks are
    completed and were successful

    :type bg: BindingsGraph
    """

    import pbsmrtpipe.pb_io as IO

    for node in bg.scattered_task_nodes():
        # skip if not already gathered
        _was_gathered = bg.node[node][ConstantsNodes.TASK_ATTR_WAS_GATHERED]
        if _was_gathered:
            log.debug("Skipping {n} was-gathered? {g}".format(n=node, g=_was_gathered))
            continue

        # The Task was already scattered and Chunk.json file is being
        # generated
        _state = bg.node[node][ConstantsNodes.TASK_ATTR_STATE]
        log.debug("Scattered task {n} {k} in state '{s}'".format(n=node, k=node.__class__, s=_state))
        # look for a completed Chunk.json file
        if _state in TaskStates.FAILURE_STATES():
            raise TaskExecutionError("Scattered task {n} failed.".format(i=node))
        elif _state != TaskStates.SUCCESSFUL:
            # wait till it's finished to process.
            continue
        else:
            # everything was good
            pass

        # Get Chunk group id
        chunk_group_id = bg.node[node][ConstantsNodes.TASK_ATTR_CHUNK_GROUP_ID]
        # slog.info("Scattered task was successful. {n} {k} in state '{s}'".format(n=node, k=node.__class__, s=_state))

        # Chunk Groups now have a fuzzy definition. A single scattering can
        # used by multiple tasks downstream tasks
        # Use the (chunk group id, operator id) as unique key
        operator_ids = bg.node[node][ConstantsNodes.TASK_ATTR_OPERATOR_ID]

        for operator_id in operator_ids:
            chunk_operator = chunk_operators_d[operator_id]

            # used with the chunk operator to find the gather task(s)
            # Un-scattered task type id
            original_task_id = node.original_task_id
            # original meta task
            original_task = registered_tasks_d[original_task_id]
            # Add check to see if
            required_chunk_keys = node.meta_task.tool_contract.task.chunk_keys

            _to_i = lambda x: int(x.split(":")[-1])

            # {index -> (chunk_key, gather task id, X)}
            gs = {_to_i(g.task_input): (g.chunk_key, g.gather_task_id, g.task_input) for g in chunk_operator.gather.chunks}

            # Get all the output chunk.json files. This should only find one
            chunk_file_node = bg.successors(node)[0]
            # load pipeline chunks generated from the task
            scattered_chunked_json_path = bg.node[chunk_file_node][ConstantsNodes.FILE_ATTR_PATH]

            if scattered_chunked_json_path is None:
                # this should raise. The task was successful, thee
                # chunk.json file should be present
                raise ChunkGatheringError("Unable to find output chunked JSON file from task {t}".format(t=chunk_file_node))

            scattered_pipeline_chunks = IO.load_pipeline_chunks_from_json(scattered_chunked_json_path)
            log.debug("Loaded {n} pipeline scattered chunks from {i} {p}".format(n=len(scattered_pipeline_chunks), p=scattered_chunked_json_path, i=node))

            # The Scattering task should have written keys that are
            # defined in the ToolContract
            for required_chunk_key in required_chunk_keys:
                found_keys = set(itertools.chain(*[p.chunk_keys for p in scattered_pipeline_chunks]))
                if required_chunk_key not in found_keys:
                    raise ChunkGatheringError("Unable to find required chunk keys {k} (chunk-operator {o}) from {f}".format(k=required_chunk_key, f=scattered_chunked_json_path, o=operator_id))

            pipeline_chunks_d = {c.chunk_id: c for c in scattered_pipeline_chunks}

            def _filter_chunks_by_group_and_operator(cnode, chunk_group_id_):
                return isinstance(cnode, TaskChunkedBindingNode) and cnode.chunk_group_id == chunk_group_id_ and cnode.operator_id == operator_id

            # Look for all TaskChunkedBindingNodes with the same chunk group
            chunked_task_states = [(chunk_node_, bg.node[chunk_node_][ConstantsNodes.TASK_ATTR_STATE]) for chunk_node_ in bg.chunked_task_nodes() if _filter_chunks_by_group_and_operator(chunk_node_, node.chunk_group_id)]

            # Check for failure in chunked tasks
            if any(s in TaskStates.FAILURE_STATES() for s, c in chunked_task_states):
                raise TaskExecutionError("Chunked task failure. {c}".format(c=chunked_task_states))

            if not chunked_task_states:
                raise ChunkGatheringError("No chunked tasks found for {t} chunk-group {g}".format(t=node, g=chunk_group_id))

            if all(s == TaskStates.SUCCESSFUL for cnode, s in chunked_task_states):
                # Check if all chunked tasks have completed and output files have been resolved
                if all(was_task_successful_with_resolve_outputs(bg, cnode) for cnode, s in chunked_task_states):

                    slog.info("Starting chunking gathering process for task {n} chunk-group {g} with operator {i}".format(n=node, g=chunk_group_id, i=operator_id))
                    gathered_pipeline_chunks_d = copy.deepcopy(pipeline_chunks_d)

                    # Found completed chunked files. Now:
                    # 1. create Gathered JSON File and GatheredFileNode
                    # 2. Create Gather tasks
                    # 3. Map output of first chunked task to input of Gathered File Node (this is a hack)

                    # bind the first output of the first chunked task. this is a hack to get around the degree constraints

                    all_chunked_out_files_nodes = []

                    for chunked_task_node, state_ in chunked_task_states:
                        log.debug(("Chunked task node", chunked_task_node, type(chunked_task_node), len(chunked_task_states)))

                        chunk_id = chunked_task_node.chunk_id

                        for output_node in bg.successors(chunked_task_node):
                            all_chunked_out_files_nodes.append(output_node)

                            # map to $chunk_key. Add better error message
                            if output_node.index not in gs:
                                log.error(("Chunk Operator {i} gather index ".format(i=operator_id), gs))
                                raise ChunkGatheringError("Chunk Operator {i} Failed to map {n}".format(i=operator_id, n=output_node))

                            output_chunk_key, _, _ = gs[output_node.index]
                            cpath = bg.node[output_node][ConstantsNodes.FILE_ATTR_PATH]

                            slog.debug("Mapping outputs of chunked task {n} Chunk {i} key: {k} {p}".format(n=output_node.idx, i=chunk_id, k=output_chunk_key, p=cpath))
                            gathered_pipeline_chunks_d[chunk_id]._datum[output_chunk_key] = bg.node[output_node][ConstantsNodes.FILE_ATTR_PATH]

                    comment = "Gathered pipeline chunks {t}. Scattered {f}".format(t=node, f=scattered_chunked_json_path)
                    gathered_json = os.path.join(tasks_root_dir, ".{t}-{u}-gathered-pipeline.chunks.json".format(t=node.meta_task.task_id, u=uuid.uuid4()))
                    write_pipeline_chunks(gathered_pipeline_chunks_d.values(), gathered_json, comment)

                    # Create New Gathered InFile Node
                    # Create all Gathered Tasks
                    # This will map the output gchunk.gather_task_id -> outs of
                    # the original unchunked tasks (that shouldn't have run)
                    # list of [(GatherChunk, out-file-node), ...]
                    g_out_gchunk_fnodes = []
                    for gi, gchunk in enumerate(chunk_operator.gather.chunks):
                        # need to make a copy here, since individual gather
                        # tasks may be applied to multiple outputs and the file
                        # labels will get mixed up otherwise

                        g_meta_task = copy.deepcopy(registered_tasks_d[gchunk.gather_task_id])
                        g_node = bg.add_gather_meta_task(g_meta_task, gchunk.chunk_key)
                        # the name, description, and sourceId fields in the
                        # datastore should be populated from the original task,
                        # not the gather task
                        g_meta_task.output_file_display_names[0] = original_task.output_file_display_names[gi]
                        g_meta_task.output_file_descriptions[0] = original_task.output_file_descriptions[gi]
                        g_meta_task.datastore_source_id = "{t}-out-{i}".format(
                            t=original_task.task_id, i=gi)

                        # Both In/Out Gather Binding Types only have one input and
                        # one output, hence, the positional in/out index is ALWAYS 0

                        # this is still using the mixed form of [(FileType, label, desc), ]
                        g_in_file = bg.add_binding_in(g_meta_task, 0, g_meta_task.input_types[0])
                        g_out_file = bg.add_binding_out(g_meta_task, 0, g_meta_task.output_types[0])

                        # update the state, path of the resolved gathered file
                        update_file_state_to_resolved(bg, g_in_file, gathered_json)

                        bg.add_edge(g_in_file, g_node)
                        bg.add_edge(g_node, g_out_file)
                        if all_chunked_out_files_nodes:
                            for out_node in all_chunked_out_files_nodes:
                                bg.add_edge(out_node, g_in_file)

                        g_out_gchunk_fnodes.append((gchunk, g_out_file, binding_str_to_task_id_and_instance_id(gchunk.task_input)))

                    # Finally remap the outputs of the original unchunked task to
                    # the outputs of the gathered chunks and delete the original task node
                    # for origin_out_node in bg.successors()

                    original_unchunked_tnode = get_companion_unscattered_task_node(bg, chunk_group_id)
                    # # {Positional-index: OutFileNode}
                    g_lookup = {x[-1][-1]: x[1] for x in g_out_gchunk_fnodes}

                    for original_out_file in bg.successors(original_unchunked_tnode):
                        new_mapped_input_node = g_lookup[original_out_file.index]
                        for mapped_in_node in bg.successors(original_out_file):
                            slog.debug("Mapping new input {i} to {o}".format(i=new_mapped_input_node, o=mapped_in_node))
                            bg.add_edge(new_mapped_input_node, mapped_in_node)

                    # Delete original task, since new gather'ed mappings have created
                    bg.remove_nodes_from(bg.successors(original_unchunked_tnode))
                    bg.remove_node(original_unchunked_tnode)

                    # update chunked task properties
                    attrs = [(ConstantsNodes.TASK_ATTR_WAS_CHUNKED, True),
                             (ConstantsNodes.TASK_ATTR_WAS_GATHERED, True),
                             (ConstantsNodes.TASK_ATTR_CHUNK_GROUP_ID, chunk_group_id),
                             (ConstantsNodes.TASK_ATTR_COMPANION_CHUNK_TASK_TYPE_ID, original_task.task_id),
                             (ConstantsNodes.TASK_ATTR_OPERATOR_ID, chunk_operator.idx)]

                    update_or_set_node_attrs(bg, attrs, [node])
                    log.debug("Updated chunked task is_running {t}".format(t=node))

                    slog.info("complete chunking task {n} chunk-group {g}".format(n=node, g=chunk_group_id))
                else:
                    log.debug("Chunked tasks not completed, or resolved for {n} group: {g}".format(n=node, g=chunk_group_id))
            else:
                log.debug("Chunked tasks for are not completed, or successful {n} group: {g}".format(n=node, g=chunk_group_id))

    resolve_successor_binding_file_path(bg)
    validate_binding_graph_integrity(bg)
    return bg


def _get_tmpdir():
    return tempfile.mkdtemp()


def _get_tmpfile(suffix=".file"):
    t = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    t.close()
    return t.name


def _get_logfile(output_dir):
    suffix = ".log"
    t = tempfile.NamedTemporaryFile(suffix=suffix, delete=False, dir=output_dir)
    t.close()
    return t.name


def resolve_di_resources(task_output_dir, specials_di):
    """Convert the ResourceTypes to Resource instances"""

    to_log = functools.partial(_get_logfile, task_output_dir)
    to_o = lambda: task_output_dir

    r = {ResourceTypes.LOG_FILE: to_log,
         ResourceTypes.TMP_DIR: _get_tmpdir,
         ResourceTypes.TMP_FILE: _get_tmpfile,
         ResourceTypes.OUTPUT_DIR: to_o}

    resolved_specials = resolve_di(r, specials_di)
    return resolved_specials


def to_resolve_di_resources(task_output_dir, root_tmp_dir=None):
    def wrapper_(resource_di):
        return resolve_di_resources(task_output_dir, resource_di)
    return wrapper_


def resolve_io_files(id_to_count, output_dir, input_files, output_file_type_list, output_files_names, mutable_files):
    """

    output_dir = str

    id_to_count {ft_id: number of times the file instances were created}

    file_di =

    :param output_files_names: Optional value which will be a list of [(base_name, ext), ] for each file

    (This allows unique files to be created.)

    """

    # this is some ugly logic
    if output_files_names is None:
        override_names = [None for _ in output_file_type_list] if output_files_names is None else output_files_names
    else:
        if len(output_files_names) != len(output_file_type_list):
            log.debug("IGNORING override file names. Incompatible file type list ({i}) and override names ({f})".format(i=len(output_file_type_list), f=len(output_files_names)))
            override_names = [None for _ in output_file_type_list]
        else:
            override_names = output_files_names

    def _to_mutable_index(s):
        """Mutable file format"""
        rx = re.compile(r'\$(inputs|outputs).([0-9]*)')
        try:
            return int(rx.match(s).groups()[1])
        except (TypeError, IndexError, AttributeError) as e:
            msg_ = "Format error in mutable file {s}. Format should match pattern {p}".format(s=s, p=rx.pattern)
            log.error(msg_)
            log.error(e)
            raise ValueError(msg_)

    # to make the book keeping easier
    mutable_indices = []
    mutable_outputfile_indices = []
    # Overwrite output default file types
    if mutable_files is not None:
        for in_mfile, out_mfile in mutable_files:
            in_i = _to_mutable_index(in_mfile)
            out_i = _to_mutable_index(out_mfile)
            mutable_indices.append((in_i, out_i))
            mutable_outputfile_indices.append(out_i)

    # this should be done via global id?
    to_p = lambda d, n: os.path.join(d, n)
    to_f = functools.partial(to_p, output_dir)

    # Output Paths
    paths = []
    # TODO FIXME. This needs to be more robust
    ftypes_overrides = zip(output_file_type_list, override_names)
    for i, x in enumerate(ftypes_overrides):

        # Mutable files take priority over overridden names
        if i in mutable_outputfile_indices:
            paths.append(input_files[i])
        else:
            ftype, override_name = x

            if override_name is None:
                base_name, ext = ('file', 'txt') if ftype is None else (ftype.base_name, ftype.ext)
            else:
                log.info(override_name)
                base_name, ext = override_name

            if ftype.file_type_id not in id_to_count:
                id_to_count[ftype.file_type_id] = 0

            instance_id = id_to_count[ftype.file_type_id] = 0

            if instance_id == 0:
                p = to_f('.'.join([base_name, ext]))
            else:
                p = to_f('.'.join([base_name + '-' + str(instance_id), ext]))

            id_to_count[ftype.file_type_id] = instance_id + 1
            paths.append(p)

    log.debug("ID Counter")
    log.debug(id_to_count)
    return paths


def to_resolve_files(file_type_id_to_count):
    # Need to create a closure over the counter which will act
    # as a global-esque var when files names are created.
    def f2(task_dir, input_files, output_types, output_files_names, mutable_files):
        return resolve_io_files(file_type_id_to_count, task_dir, input_files, output_types, output_files_names, mutable_files)
    return f2
