import logging
import json
import datetime
import os
import functools

import networkx as nx

import pbsmrtpipe.external_tools as ET

from pbsmrtpipe.models import TaskStates
from pbsmrtpipe.graph.models import (ConstantsNodes,
                                     DotColorConstants,
                                     DotStyleConstants,
                                     DotShapeConstants,
                                     TaskBindingNode,
                                     TaskScatterBindingNode,
                                     TaskChunkedBindingNode,
                                     EntryOutBindingFileNode,
                                     BindingChunkInFileNode, BindingInFileNode,
                                     BindingOutFileNode, EntryPointNode,
                                     BindingChunkOutFileNode,
                                     VALID_FILE_NODE_CLASSES,
                                     VALID_TASK_NODE_CLASSES)
log = logging.getLogger(__name__)


def bindings_graph_to_dict(bg):

    _d = dict(_comment="Updated at {x}".format(x=datetime.datetime.now()))

    def _to_a(n, name):
        return bg.node[n][name]

    def _to_d(n_):
        _x = {a: _to_a(n_, a) for a in n_.NODE_ATTRS.keys()}
        _x['klass'] = n_.__class__.__name__
        _x['node_id'] = n.idx
        return _x

    nodes = []
    edges = set([])

    for n in bg.nodes():
        nodes.append(_to_d(n))

    _d['nodes'] = nodes
    _d['nnodes'] = len(nodes)
    _d['edges'] = edges
    _d['nedges'] = len(edges)

    return _d


class DateTimeEncoder(json.JSONEncoder):

    def default(self, o):
        if isinstance(o, datetime.datetime):
            return o.isoformat()


def write_bindings_graph_to_json(bg, path):

    d = bindings_graph_to_dict(bg)

    with open(path, 'w+') as w:
        w.write(json.dumps(d, indent=4, sort_keys=True, cls=DateTimeEncoder))


def write_binding_graph_images(g, root_dir):

    dot_file = os.path.join(root_dir, 'workflow.dot')
    s = binding_graph_to_dot(g)
    with open(dot_file, 'w') as f:
        f.write(s)

    formats = ('png', 'svg')
    for f in formats:
        p = os.path.join(root_dir, '.'.join(['workflow', f]))
        ET.dot_to_image(f, dot_file, p)

    workflow_json = os.path.join(root_dir, 'workflow-graph.json')
    write_bindings_graph_to_json(g, workflow_json)


def binding_graph_to_dot(g):
    """
    Custom generation of dot file format

    :rtype: str
    """

    outs = ["strict digraph G {"]
    _add = outs.append

    def _to_s(n):
        return "".join(['"', str(n), '"'])

    def _to_l(n):
        return n + " ;"

    def _to_attr_s(d):
        # Convert to a dot-friendly format
        # [ style=filled shape=rectangle fillcolor=grey]
        vs = " ".join(['='.join([k, v]) for k, v in d.iteritems()])
        attrs_str = " ".join(['[', vs, ']'])
        return attrs_str

    def _get_attr_or_default(g_, n_, attr_name_, default_value):
        try:
            return g_.node[n_][attr_name_]
        except KeyError:
            return default_value

    def _node_to_view_d(g_, n_):
        _d = {}
        _d['fillcolor'] = n_.DOT_COLOR
        _d['color'] = n_.DOT_COLOR
        _d['style'] = DotStyleConstants.FILLED
        _d['shape'] = n_.DOT_SHAPE
        return _d

    def _node_to_dot(g_, n_):
        s = _to_s(n_)
        _d = _node_to_view_d(g_, n_)
        attrs_str = _to_attr_s(_d)
        return ' '.join([s, attrs_str])

    def _task_node_to_dot(g_, n_):
        _d = _node_to_view_d(g_, n_)
        state = g_.node[n_]['state']
        state_color = DotColorConstants.RED if state == TaskStates.FAILED else node.DOT_COLOR
        _d['fillcolor'] = state_color
        _d['color'] = state_color

        # Chunk Operator Id
        operator_ids = _get_attr_or_default(g_, n_, 'operator_id', None)
        if operator_ids is not None:
            #_d['operator_ids'] = " ".join([operator_id.replace('.', '_') for operator_id in operator_ids])
            pass

        attrs_str = _to_attr_s(_d)
        return ' '.join([_to_s(n_), attrs_str])

    def _binding_file_to_dot(g_, n_):
        s = _to_s(n_)
        _d = _node_to_view_d(g_, n_)
        is_resolved = g_.node[n_][ConstantsNodes.FILE_ATTR_IS_RESOLVED]
        if not is_resolved:
            _d['style'] = DotStyleConstants.DOTTED
        attrs_str = _to_attr_s(_d)
        return ' '.join([s, attrs_str])

    # write the node metadata
    for node in g.nodes():
        funcs = {TaskBindingNode: _task_node_to_dot,
                 TaskChunkedBindingNode: _task_node_to_dot,
                 TaskScatterBindingNode: _task_node_to_dot,
                 EntryOutBindingFileNode: _node_to_dot,
                 BindingInFileNode: _binding_file_to_dot,
                 BindingChunkInFileNode: _binding_file_to_dot,
                 BindingOutFileNode: _binding_file_to_dot,
                 BindingChunkOutFileNode: _binding_file_to_dot,
                 EntryPointNode: _node_to_dot}

        f = funcs.get(node.__class__, _node_to_dot)
        x = f(g, node)
        _add(_to_l(x))

    for i, f in g.edges():
        s = ' -> ' .join([_to_s(i), _to_s(f)])
        _add(_to_l(s))

    _add("}")
    return "\n".join(outs)


def to_binding_graph_summary(bg):
    """
    General func for getting a summary of BindingGraph instance
    """

    header = "Binding Graph Status Summary"
    _n = 80
    sp = "-" * _n
    ssp = "*" * _n
    outs = []
    _add = outs.append
    _add_sp = functools.partial(_add, sp)
    _add_ssp = functools.partial(_add, ssp)

    _add_ssp()
    _add(header)
    _add_sp()
    _add("Workflow complete? {c}".format(c=bg.is_workflow_complete()))

    tn_states = {s: bg.get_tasks_by_state(s) for s in TaskStates.ALL_STATES()}
    tn_s = " ".join([":".join([k, str(v)]) for k, v in tn_states.iteritems()])

    _add_sp()
    _add("Task Summary {n} tasks ({s})".format(n=len(bg.task_nodes()), s=str(tn_s)))
    _add_sp()

    sorted_nodes = nx.topological_sort(bg)

    _add(" ".join(["resolved inputs".ljust(20), "resolved outputs".ljust(20), "state".ljust(12), "NodeType".ljust(30), "N inputs".ljust(12), "N outputs".ljust(12), "run time".ljust(12), "Id".ljust(60), ]))
    _add_sp()
    for tnode in sorted_nodes:
        if isinstance(tnode, VALID_TASK_NODE_CLASSES):
            state = bg.node[tnode]['state']
            _is_resolved = lambda it: all(bg.node[n][ConstantsNodes.FILE_ATTR_IS_RESOLVED] for n in it)
            inputs_resolved = _is_resolved(bg.predecessors(tnode))
            outputs_resolved = _is_resolved(bg.successors(tnode))
            ninputs = len(bg.predecessors(tnode))
            noutputs = len(bg.successors(tnode))

            run_time = bg.node[tnode]['run_time']

            s = str(inputs_resolved).ljust(20), str(outputs_resolved).ljust(20), state.ljust(12), tnode.__class__.__name__.ljust(30), str(ninputs).ljust(12), str(noutputs).ljust(12), str(run_time).ljust(12), str(tnode).ljust(60)
            _add(" ".join(s))

    # File-esque summary
    def _to_summary(fnode_):
        # print type(fnode_), fnode_
        s = bg.node[fnode_][ConstantsNodes.FILE_ATTR_IS_RESOLVED]
        p = bg.node[fnode_][ConstantsNodes.FILE_ATTR_PATH]
        if p is None:
            ppath = str(None)
        else:
            ppath = '... ' + str(p)[-35:]
        _add(" ".join([str(s).ljust(10), fnode.__class__.__name__.ljust(18), ppath.ljust(40), str(fnode)]))

    _add("")
    _add_sp()
    _add("File Summary ({n}) files".format(n=len(bg.file_nodes())))
    _add_sp()
    _add(" ".join(["resolved".ljust(10), "NodeType".ljust(18), "Path".ljust(40), "Id"]))
    _add_sp()

    for fnode in sorted_nodes:
        if isinstance(fnode, VALID_FILE_NODE_CLASSES):
            _to_summary(fnode)

    _add("")
    _add_sp()
    _add("Entry Point Node Summary ({n})".format(n=len(bg.entry_point_nodes())))
    _add_sp()
    _add(" ".join(["resolved".ljust(10), "NodeType".ljust(18), "Path".ljust(40), "Id"]))
    _add_sp()

    for x in bg.entry_point_nodes():
        _to_summary(x)

    _add_sp()

    return "\n".join(outs)


def to_binding_graph_task_summary(bg):
    """

    :type bg: BindingsGraph
    :param bg:
    :return:
    """
    def _to_p(p_or_none):
        if p_or_none is None:
            return None
        else:
            return '... ' + str(p_or_none)[-35:]

    for x, tnode in enumerate(bg.all_task_type_nodes()):
        print "Task {x} ".format(x=x), bg.node[tnode]['state'], tnode
        inodes = bg.predecessors(tnode)
        print "Inputs:"
        for i in inodes:
            print "Path ", _to_p(bg.node[i]['path']), bg.node[i][ConstantsNodes.FILE_ATTR_IS_RESOLVED], i

        print "Outputs:"
        onodes = bg.successors(tnode)
        for i in onodes:
            print "Path ", _to_p(bg.node[i]['path']), bg.node[i][ConstantsNodes.FILE_ATTR_IS_RESOLVED], i

        task = bg.node[tnode].get('task', None)
        if task is not None:
            for i, output in enumerate(task.output_files):
                print "Output ", i, output

        print ""

    return "Summary"
