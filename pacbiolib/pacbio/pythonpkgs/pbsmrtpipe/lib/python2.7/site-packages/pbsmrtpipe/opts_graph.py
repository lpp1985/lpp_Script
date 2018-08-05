"""Utils for Creating and Resolving the DI options graph"""
import sys
import inspect
import json
import logging
import functools
import re
import types
import pprint
import jsonschema

from pbcommand.models import FileTypes, TaskTypes, SymbolTypes
from pbsmrtpipe.exceptions import (InvalidDependencyInjectError,
                                   MalformedDependencyInjectionFileMetadataError)

from pbsmrtpipe.models import (MetaScatterTask,
                               MetaTask, ScatterTask, Task, MetaGatherTask,
                               GatherTask, ScatterToolContractMetaTask,
                               GatherToolContractMetaTask,
                               ToolContractMetaTask)

from pbsmrtpipe.dataset_io import (dispatch_metadata_resolver,
                                   has_metadata_resolver,
                                   DatasetMetadata)


import networkx as nx

log = logging.getLogger(__name__)
# logging.basicConfig(level=logging.DEBUG)


def get_report_json_attribute(report_file, attribute_id):

    try:
        with open(report_file, 'r') as f:
            s = f.read()

        d = json.loads(s)
        for attr in d['attributes']:
            if attr['id'] == attribute_id:
                value = attr['value']
                log.debug("Extracted report attribute {a} from {r}".format(a=attribute_id, r=report_file))
                return value
    except (ValueError, IOError, KeyError) as e:
        msg = "Unable to load report attribute '{a}' from {p}".format(a=attribute_id, p=report_file)
        log.error(msg)
        raise

    raise InvalidDependencyInjectError("Unable to find attribute '{a}' in report {r}".format(a=attribute_id, r=report_file))


def get_dataset_metadata_from_file(file_type, path, attribute_name):

    if attribute_name not in DatasetMetadata.SUPPORTED_ATTRS:
        raise InvalidDependencyInjectError("File metadata Attribute '{a}' is not support in {p}".format(a=attribute_name, p=path))

    dataset_metadata = dispatch_metadata_resolver(file_type, path)

    return getattr(dataset_metadata, attribute_name)


def get_metadata_from_file(file_type, path, attribute_name):

    if file_type == FileTypes.REPORT:
        return get_report_json_attribute(path, attribute_name)
    else:
        return get_dataset_metadata_from_file(file_type, path, attribute_name)


def resolve_di(resolved_keys_d, di_list_model):
    """
    Resolve DI values in a list

    :param resolved_keys_d: {di key : func} The Func has empty signature
    :param di_list_model: DI list model (e.g., ['$tmpfile', '$tmpdir','$nproc'])
    :return: a list of resolved DI values
    """

    resolved_specials = []
    for s in di_list_model:
        if s in resolved_keys_d:
            func = resolved_keys_d[s]
            v = func()
            resolved_specials.append(v)
        else:
            _d = dict(s=s, t=resolved_keys_d.keys())
            raise ValueError("Unable to resolve special '{s}'. Valid specials are '{t}'".format(**_d))

    return resolved_specials


def get_tail_func_or_raise(alist):
    f = alist[-1]

    if isinstance(f, types.FunctionType):
        return f
    else:
        _d = dict(a=alist, f=f, t=type(f))
        raise TypeError("Invalid DI definition for '{a}'. Expected func type {t} for {f}".format(**_d))


def to_di_graph_with_funcs(meta_task):
    """

    :type meta_task: MetaTask

    """
    # backward compatible for task type
    TASK_TYPES = {True: "pbsmrtpipe.task_types.distributed", False: "pbsmrtpipe.task_types.local"}

    to_funcs = {'to_nproc': None,
                'to_task_type': None,
                'to_nchunks': None,
                'to_ropts': None}

    g = nx.DiGraph()

    def _add_di_nodes(alist, to_func_node_id):
        for i in alist:
            if isinstance(i, str):
                if i.startswith('$'):
                    g.add_node(i)
                    g.add_edge(i, to_func_node_id)
            elif isinstance(i, dict):
                # FIXME skip these for now.
                continue
            else:
                # Just add stuff?
                x = str(i)
                g.add_node(x)
                g.add_edge(x, to_func_node_id)

    # add nodes for task type
    g.add_node('to_task_type', shape="rectangle")
    g.add_node(SymbolTypes.TASK_TYPE, style="filled", fillcolor="aquamarine")
    g.add_edge('to_task_type', SymbolTypes.TASK_TYPE)

    if isinstance(meta_task.is_distributed, bool):
        task_type_label = TASK_TYPES[meta_task.is_distributed]
        g.add_node(task_type_label)
        g.add_edge(task_type_label, 'to_task_type')
    else:
        raise ValueError("Unsupported type")

    # add option dependencies
    to_ropts_fname = "to_ropts"
    if isinstance(meta_task.option_schemas, (list, tuple)):
        f = get_tail_func_or_raise(meta_task.option_schemas)
        to_funcs['to_ropts'] = f
        to_ropts_fname += "--{f}".format(f=f.__name__)
        # the raw opts schema is given as first arg
        xs = meta_task.option_schemas[:-1]
        xs.pop(0)
        xs.insert(0, SymbolTypes.OPTS)

        _add_di_nodes(meta_task.option_schemas[:-1], to_ropts_fname)

    elif isinstance(meta_task.option_schemas, dict):
        # no DI
        pass
    else:
        raise ValueError("Malformed option DI. '{x}' for task {t}".format(x=meta_task.option_schemas, t=meta_task.task_id))

    # add option dependencies
    g.add_node(SymbolTypes.OPTS)
    g.add_node(SymbolTypes.SCHEMA_OPTS)
    g.add_node(SymbolTypes.RESOLVED_OPTS, style="filled", fillcolor="aquamarine")
    g.add_node(to_ropts_fname, shape="rectangle")
    g.add_edge(to_ropts_fname, SymbolTypes.RESOLVED_OPTS)
    g.add_edge(SymbolTypes.OPTS, to_ropts_fname)
    g.add_edge(SymbolTypes.SCHEMA_OPTS, to_ropts_fname)

    # add nproc
    to_nproc_fname = 'to_nproc'

    if isinstance(meta_task.nproc, (list, tuple)):
        f = get_tail_func_or_raise(meta_task.nproc)
        to_funcs['to_nproc'] = f
        to_nproc_fname += "--{f}".format(f=f.__name__)
        _add_di_nodes(meta_task.nproc[:-1], to_nproc_fname)
    else:
        if meta_task.nproc == SymbolTypes.MAX_NPROC:
            g.add_node(SymbolTypes.MAX_NPROC)
            g.add_edge(SymbolTypes.MAX_NPROC, to_nproc_fname)
        elif isinstance(meta_task.nproc, int):
            g.add_node(meta_task.nproc)
            g.add_edge(meta_task.nproc, to_nproc_fname)
        else:
            raise TypeError("expected primitive value for nproc. Got type {t} from '{v}'".format(t=type(meta_task.nproc), v=meta_task.nproc))

    # add outputs of func to $nproc
    g.add_node(SymbolTypes.NPROC, style="filled", fillcolor="aquamarine")
    g.add_node(to_nproc_fname, shape="rectangle")
    g.add_edge(to_nproc_fname, SymbolTypes.NPROC)

    to_nchunks_fname = "to_nchunks"

    if isinstance(meta_task, MetaScatterTask):
        # Handle DI case
        if isinstance(meta_task.chunk_di, (list, tuple)):
            f = get_tail_func_or_raise(meta_task.chunk_di)
            to_funcs['to_nchunks'] = f
            to_nchunks_fname += "--{f}".format(f=f.__name__)
            _add_di_nodes(meta_task.chunk_di[:-1], to_nchunks_fname)
        else:
            # primitive value was given
            if meta_task.chunk_di == SymbolTypes.MAX_NCHUNKS:
                g.add_node(SymbolTypes.MAX_NCHUNKS)
                g.add_edge(SymbolTypes.MAX_NCHUNKS, to_nchunks_fname)
            elif isinstance(meta_task.chunk_di, int):
                g.add_node(meta_task.chunk_di)
                g.add_edge(meta_task.chunk_di, to_nchunks_fname)
            else:
                raise TypeError("Expected primitive value or {x} for nchunks. Got type {t} for {v}".format(t=type(meta_task.chunk_di), v=meta_task.chunk_di, x=SymbolTypes.MAX_NCHUNKS))

        g.add_node(SymbolTypes.NCHUNKS, style="filled", fillcolor="aquamarine")
        g.add_node(to_nchunks_fname, shape="rectangle")
        g.add_edge(to_nchunks_fname, SymbolTypes.NCHUNKS)

    log.debug(to_funcs)
    return g, to_funcs


def to_di_graph(meta_task):
    g, _ = to_di_graph_with_funcs(meta_task)
    return g


def get_report_di(meta_task):
    """
    Get all DI models from meta task and extract the '$inputs' that have
    report JSON

    :param meta_task: MetaTask

    Returns a list of (inputs index, attribute_id)
    """
    report_dis = []

    attr_names = ['option_schemas', 'nproc', 'is_distributed']

    if isinstance(meta_task, MetaScatterTask):
        attr_names.append('chunk_di')

    # DI $inputs.[0-9].attribute_id
    rx = re.compile(r'\$inputs.(\d).([A-z0-9_\.]*)')

    for attr_name in attr_names:
        value_or_list = getattr(meta_task, attr_name)
        if is_di_list(value_or_list):
            for v in value_or_list:
                if is_dollar_value(v):
                    if rx.match(v) is not None:

                        m = rx.match(v)
                        if m is not None:
                            gs = m.groups()

                            input_file_index = int(gs[0])
                            attr_id = gs[1]

                            file_type = meta_task.input_types[input_file_index]

                            if has_metadata_resolver(file_type):
                                if attr_id in DatasetMetadata.SUPPORTED_ATTRS:
                                    report_dis.append((input_file_index, attr_id))
                                else:
                                    MalformedDependencyInjectionFileMetadataError("Unsupported metadata attribute '{a}' for File type {f}. Supported values {v}".format(a=attr_id, f=file_type, v=DatasetMetadata.SUPPORTED_ATTRS))
                            elif file_type == FileTypes.REPORT:
                                report_dis.append((input_file_index, attr_id))
                            else:
                                _d = dict(i=input_file_index, f=file_type, p=meta_task.input_types, r=FileTypes.REPORT, t=meta_task.task_id)
                                raise MalformedDependencyInjectionFileMetadataError("Incompatible file type for file index {i} {f} for task id '{t}'. Expected '{r}' type. Inputs types {p}".format(**_d))
                        else:
                            raise MalformedDependencyInjectionFileMetadataError("Malformed report input DI '{i} for task '{t}'".format(i=v, t=meta_task.task_id))

    return report_dis


def is_di_list(alist_or_item):
    """
    Validate a list or item is a properly constructed DI list

    A valid DI is

    1. Last item is a function
    2. the function has the same signature as len(di_list) - 1
    3. All DI items are Symbol types, ints, strings
    """
    if isinstance(alist_or_item, (list, tuple)):
        if isinstance(list(alist_or_item)[-1], types.FunctionType):
            return True
    return False


def is_dollar_value(x_):
    if isinstance(x_, str):
        if x_.startswith('$'):
            return True
    return False


def is_valid(schema, v):
    """Returns a bool if the schema is valid"""
    try:
        jsonschema.validate(v, schema)
        return True
    except jsonschema.ValidationError:
        pass
    return False


def default_to_ropts(user_opts, opts_schemas):
    """
    'Resolves' the options and returns a {id:value}

    or raises jsonschema.ValidationError
    """
    ropts = {}
    for opt_id, schema in opts_schemas.iteritems():
        if opt_id in user_opts:
            v = user_opts[opt_id]
            jsonschema.validate({opt_id: v}, schema)
            ropts[opt_id] = v
        else:
            # must have a default value or null?
            v = schema['pb_option']['default']
            ropts[opt_id] = v

    return ropts


def meta_task_to_task(meta_task,
                      input_files,
                      all_task_options,
                      output_dir,
                      max_nproc, max_nchunks,
                      to_resources_func,
                      to_resolve_files_func):
    """
    Converts a MetaTask to a Task instance

    :param input_files: List of resolved input files
    :param all_task_options: Raw Unresolved task options provided by the user {id:value}
    :param output_dir: Root path to write output to (or Task-id dir?)
    :param to_resources_func: Function that has signature (resources_di) -> resources
    :param to_resolve_files_func: Function (output_types, override_file_names, mutable_files=None) -> output file paths

    :type max_nchunks: int
    :type max_nproc: int
    :type output_dir: str
    :type to_resources_func: types.FunctionType

    :rtype: Task or ScatterTask or GatherTask
    """
    # Type checking

    if isinstance(to_resolve_files_func, functools.partial):
        f1 = to_resolve_files_func.func
        args_obj = inspect.getargspec(f1)
        if (len(args_obj.args) - len(to_resolve_files_func.args)) != 5:
            TypeError("Incorrect Function {f} args {a}".format(f=str(to_resolve_files_func), a=args_obj.args))
    elif isinstance(to_resolve_files_func, types.FunctionType):
        args_obj = inspect.getargspec(to_resolve_files_func)
        if len(args_obj.args) != 5:
            raise TypeError("Incorrect Function {f} args {a}".format(f=str(to_resolve_files_func), a=args_obj.args))
    else:
        raise TypeError("Expected function type. Got {f}".format(f=type(to_resolve_files_func)))

    v_to_resolve = [SymbolTypes.OPTS, SymbolTypes.SCHEMA_OPTS,
                    SymbolTypes.RESOLVED_OPTS, SymbolTypes.TASK_TYPE,
                    SymbolTypes.NPROC]

    # filter the opts with only the options that in the opts schema
    opts = {k: v for k, v in all_task_options.iteritems() if k in meta_task.option_schemas}

    # these are already resolved values
    resolved_values = {SymbolTypes.OPTS: opts,
                       SymbolTypes.SCHEMA_OPTS: meta_task.option_schemas,
                       SymbolTypes.MAX_NPROC: max_nproc,
                       SymbolTypes.MAX_NCHUNKS: max_nchunks}

    log.debug("Initial resolved DI values")
    log.debug(pprint.pformat(resolved_values.keys()))

    # [(input file index, report attribute id), ]
    report_dis = get_report_di(meta_task)

    if report_dis:
        log.debug("Report Input DIs")
        log.debug(report_dis)
        for f_index, attr_name in report_dis:
            k = '.'.join(['$inputs', str(f_index), attr_name])
            report_file_name = input_files[f_index]
            file_type = meta_task.input_types[f_index]
            v = get_metadata_from_file(file_type, report_file_name, attr_name)
            resolved_values[k] = v

    def _default_resolve_nproc():
        if meta_task.nproc == SymbolTypes.MAX_NPROC:
            return max_nproc
        elif isinstance(meta_task.nproc, int):
            return min(meta_task.nproc, max_nproc)
        else:
            return 1

    def _default_task_type():
        return meta_task.is_distributed

    def _default_nchunks():
        # the old model should be removed.
        if isinstance(meta_task, (MetaScatterTask, )):
            if meta_task.chunk_di == SymbolTypes.MAX_NCHUNKS:
                return max_nchunks
            elif isinstance(meta_task.chunk_di, int):
                return min(meta_task.chunk_di, SymbolTypes.MAX_NCHUNKS)

        # default to the max number of chunks. At the tool level, it
        # should sort out how to the max number of chunks.
        log.info("Default nchunks {x}".format(x=max_nchunks))
        return max_nchunks

    default_resolve_ropts = functools.partial(default_to_ropts, resolved_values[SymbolTypes.OPTS], resolved_values[SymbolTypes.SCHEMA_OPTS])

    # Default Resolution Functions. They have no args
    default_funcs = {'to_ropts': default_resolve_ropts,
                     'to_nproc': _default_resolve_nproc,
                     'to_task_type': _default_task_type,
                     'to_nchunks': _default_nchunks}

    nchunks_ = _default_nchunks()

    # the resolved func with update the resolved values with the value of dict
    # {str:(Symbol, func() -> value)}
    method_id_to_dollar_t = {'to_ropts': (SymbolTypes.RESOLVED_OPTS, meta_task.option_schemas),
                             'to_nproc': (SymbolTypes.NPROC, meta_task.nproc),
                             'to_nchunks': (SymbolTypes.NCHUNKS, nchunks_),
                             'to_task_type': (SymbolTypes.TASK_TYPE, meta_task.is_distributed)}

    if isinstance(meta_task, MetaScatterTask):
        method_id_to_dollar_t['to_nchunks'] = (SymbolTypes.NCHUNKS, meta_task.chunk_di)

    def get_method_id_prefix(s):
        if isinstance(s, str):
            for to_x in method_id_to_dollar_t.keys():
                if s.startswith(to_x):
                    return to_x
        return False

    def begins_with_method_prefix(s):
        return get_method_id_prefix(s) is not False

    g = to_di_graph(meta_task)

    nodes = nx.topological_sort(g)
    log.debug(pprint.pformat(nodes))

    for node in nodes:
        # log.debug("Trying to resolve '{r}'".format(r=node))
        if node in resolved_values:
            # nothing to do here
            # log.debug("Value {n} has been resolved. '{v}'".format(n=node, v=resolved_values[node]))
            pass
        else:
            # Are we using default funcs to resolve values
            if begins_with_method_prefix(node):

                method_prefix = get_method_id_prefix(node)

                dollar_key, di_values = method_id_to_dollar_t[method_prefix]

                if is_di_list(di_values):
                    f = get_tail_func_or_raise(di_values)

                    # this will be resolved args to func
                    # For example, [$a, $b, func(a, b)]
                    # then, the resolved value will computed via
                    # value = func(a, b)
                    injectable = []
                    for x in di_values[:-1]:
                        if is_dollar_value(x):
                            if x in resolved_values:
                                v = resolved_values[x]
                                injectable.append(v)
                            else:
                                raise ValueError("$ value '{x}' not resolved.".format(x=x))
                        else:
                            injectable.append(x)
                    # log.debug(f.__name__)
                    # log.debug(injectable)

                    # this should do an argsinpsect. Putting a TypeError try/catch
                    # could be misleading
                    value = f(*injectable)

                    log.debug("resolved '{k}' -> '{v}'".format(k=dollar_key, v=value))
                    resolved_values[dollar_key] = value
                else:
                    # use default value, it was supplied as a primitive.
                    # this still needs to validate the final value.
                    # nchunks can still be computed to be > $MAX_NCHUNKS

                    f = default_funcs[method_prefix]
                    value = f()
                    log.debug("resolved '{k}' -> '{v}' (default resolution)".format(k=dollar_key, v=value))
                    resolved_values[dollar_key] = value
            else:
                msg = "potentially unsupported value '{n}'".format(n=node)
                #log.warn(msg)
                #raise ValueError(msg)

    # Sanity Check to make sure required values are resolved
    for x in v_to_resolve:
        if x not in resolved_values:
            raise ValueError("{x} was never resolved. Resolved keys {k}".format(x=x, k=resolved_values.keys()))

    # Resolve Resources
    rfiles = to_resources_func(meta_task.resource_types)

    # override output file names (if necessary)
    ofiles = to_resolve_files_func(output_dir, input_files, meta_task.output_types, meta_task.output_file_names, meta_task.mutable_files)
    log.debug(("Resolved output files", ofiles))

    if isinstance(meta_task, (MetaScatterTask, ScatterToolContractMetaTask)):
        cmd_str = meta_task.to_cmd(input_files, ofiles, resolved_values[SymbolTypes.RESOLVED_OPTS], resolved_values[SymbolTypes.NPROC], rfiles, nchunks_)
        t = ScatterTask(meta_task.task_id, meta_task.is_distributed, input_files, ofiles,
                        resolved_values[SymbolTypes.RESOLVED_OPTS], resolved_values[SymbolTypes.NPROC], rfiles, cmd_str, resolved_values[SymbolTypes.NCHUNKS], output_dir, meta_task.chunk_keys)
    elif isinstance(meta_task, (MetaGatherTask, GatherToolContractMetaTask)):
        cmd_str = meta_task.to_cmd(input_files, ofiles, resolved_values[SymbolTypes.RESOLVED_OPTS], resolved_values[SymbolTypes.NPROC], rfiles)
        t = GatherTask(meta_task.task_id, meta_task.is_distributed, input_files, ofiles, resolved_values[SymbolTypes.RESOLVED_OPTS], resolved_values[SymbolTypes.NPROC], rfiles, cmd_str, output_dir)
    elif isinstance(meta_task, (MetaTask, ToolContractMetaTask)):
        cmd_str = meta_task.to_cmd(input_files, ofiles, resolved_values[SymbolTypes.RESOLVED_OPTS], resolved_values[SymbolTypes.NPROC], rfiles)
        t = Task(meta_task.task_id, meta_task.is_distributed, input_files, ofiles, resolved_values[SymbolTypes.RESOLVED_OPTS], resolved_values[SymbolTypes.NPROC], rfiles, cmd_str, output_dir)
    else:
        raise TypeError("Unsupported meta task type {m}".format(m=meta_task))

    # log.debug(t.__dict__)

    return t
