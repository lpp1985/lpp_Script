"""General validation functions"""

import logging
import types
import inspect

from pbcommand.models import (FileTypes, TaskTypes, SymbolTypes, ResourceTypes, FileType)
from pbsmrtpipe.exceptions import (MalformedMetaTaskError, MalformedPipelineError)

from pbsmrtpipe.constants import RX_VERSION, RX_CHUNK_KEY, RX_TASK_ID


log = logging.getLogger(__name__)


def to_list_if_necessary(tuple_or_s):
    if isinstance(tuple_or_s, tuple):
        return list(tuple_or_s)
    return tuple_or_s


def __validate_output_file_names(output_types, output_file_names):
    errors = []
    if isinstance(output_file_names, (list, tuple)):
        for x in output_file_names:
            if isinstance(x, (list, tuple)):
                if len(x) == 2:
                    pass
                else:
                    errors.append("Malformed output file name {x}. Expected 2-tuple (str, str).".format(x=x))
        if len(output_file_names) == len(output_types):
            # this is only branch where the outputs and override outputs file names are valid
            return True
        else:
            errors.append("Malformed output file names. Expected {n}. Got {m}".format(n=output_types, m=output_file_names))
    else:
        errors.append("Malformed output file name. Expected [(str, str)]. Got {m}".format(m=output_file_names))

    log.error("\n".join(errors))
    raise ValueError("\n".join(errors))


def _validate_output_file_names(output_types, output_files_names_or_none):
    if output_files_names_or_none is not None:
        __validate_output_file_names(output_types, output_files_names_or_none)
        out_names = output_files_names_or_none
    else:
        out_names = []

    return out_names


def _validate_output_file_names_or_raise(output_types, output_file_names_or_none, task_id):
    try:
        return _validate_output_file_names(output_types, output_file_names_or_none)
    except Exception:
        msg = "task id {i} Output file types ({n}) and override output file names ({m}) supplied are NOT compatible ".format(n=len(output_types), m=len(output_file_names_or_none), i=task_id)
        log.error(msg)
        raise


def __validate_provided_file_types(file_types):
    def _is_valid(x):
        if not isinstance(x, FileType):
            raise TypeError("Invalid FileType. Got {t}".format(t=type(x)))

    for i in file_types:
        _is_valid(i)

    return file_types


def validate_provided_file_types(file_type_or_file_types):
    if isinstance(file_type_or_file_types, FileType):
        return [file_type_or_file_types]

    return __validate_provided_file_types(file_type_or_file_types)


def __is_schema(d):
    try:
        p = d['property']
        keys = p.keys()
        oid = keys[0]
        if len(keys) != 1:
            raise IndexError
        _ = d['property'][oid]['default']
        return True
    except (KeyError, IndexError):
        return False


def _is_schema_list(value):
    if isinstance(value, dict):
        if value:
            __is_schema(value)
            return value
        else:
            # empty
            return value

    if isinstance(value, (list, tuple)):
        return all(__is_schema(x) for x in value)

    return False


def _validate_mutable_files(mutable_files_or_none, input_types, output_types):
    """
    1. Valid well-formed string
    2. Valid that indices are pointed
    3. Valid the input/output types are the same

    :param mutable_files_or_none: [("$input.0, "$output.0"), ]
    :param input_types:
    :param output_types:
    :return:

    """
    if mutable_files_or_none is None:
        return ()

    for mutable_file in mutable_files_or_none:

        if len(mutable_file) != 2:
            raise ValueError("Malformed mutable file name {n}. Expected tuple of (str, str)".format(n=mutable_file))

        in_f, out_f = mutable_file
        if not in_f.startswith("$inputs.") or not out_f.startswith("$outputs.0"):
            raise ValueError("Malformed mutable file ({i}, {f})".format(i=in_f, f=out_f))

        ix = int(in_f.split("$inputs.")[-1])
        ox = int(out_f.split("$outputs.")[-1])
        ixt = input_types[ix]
        oxt = output_types[ox]
        if ixt != oxt:
            log.warn("Mutable types are different. {i} {t}  -> {o} {f}".format(i=ix, t=ixt, o=ox, f=oxt))

    # we got here, everything is fine
    return mutable_files_or_none


def _validate_func_with_n_args(nargs, func):
    p = inspect.getargspec(func)
    if len(p.args) != nargs:
        raise ValueError("Expected func {x} to have {n} args. Got {y}".format(x=func.__name__, n=nargs, y=len(p.args)))
    return func


def _validate_chunk_only_input_type(in_out_types):
    if len(in_out_types) == 1:
        if in_out_types[0] == FileTypes.CHUNK:
            return True
    raise ValueError("Expected chunk file type in {i}".format(i=in_out_types))


def _get_class_attr_or_raise(class_name, attr, d):
    if attr not in d:
        raise MalformedMetaTaskError("MetaTask class '{c}' is missing required class var '{n}'".format(c=class_name, n=attr))
    else:
        return d[attr]


def _raise_malformed_task_attr(msg):
    def _wrapper(m=None):
        msg_ = msg + " " + m if m is not None else msg
        raise MalformedMetaTaskError(msg_)
    return _wrapper


def validate_task_type(x):

    _raise = _raise_malformed_task_attr("IS_DISTRIBUTED must be a DI List or primitive value {x}.".format(x=bool))

    if isinstance(x, bool):
        return x
    else:
        _raise("Incompatible type. Expected bool")

    return x


def _validate_in_out_types(x):

    _raise = _raise_malformed_task_attr("In/Out types must be defined as list of (FileTypes.FILE, label, description) or a single value FileTypes.FILE")
    processed_in_out_types = []

    if isinstance(x, (list, tuple)):
        for i in x:
            _raise_type = lambda: _raise("Expected FileType. Got {t}".format(t=type(i)))
            # Support the new and old format
            if isinstance(i, (list, tuple)):
                if len(i) == 3:
                    if isinstance(i[0], FileType):
                        processed_in_out_types.append(i[0])
                    else:
                        _raise_type()
                else:
                    _raise_type()
            elif isinstance(i, FileType):
                processed_in_out_types.append(i)
            else:
                _raise_type()

    else:
        _raise("Got type {t}".format(t=type(x)))

    return processed_in_out_types


def _validate_chunk_in_out_type(msg):
    def _f(x):
        _raise = _raise_malformed_task_attr("{m} Got {x}".format(x=x, m=msg))
        x = _validate_in_out_types(x)
        if isinstance(x, (list, tuple)):
            if len(x) == 1:
                if x[0] is FileTypes.CHUNK:
                    return x
        _raise()

    return _f


def _validate_scatter_output_types(x):
    _f = _validate_chunk_in_out_type("Scatter outputs type must be ONLY one chunk file type")
    return _f(x)


def _validate_gather_input_types(x):
    _f = _validate_chunk_in_out_type("Gather input type must be ONLY one chunk file type")
    return _f(x)


def _validate_gather_output_types(x):
    _raise = _raise_malformed_task_attr("Gather output types must be a single file type. Got {x}".format(x=x))
    if isinstance(x, (list, tuple)):
        # Only one output is allowed
        if len(x) == 1:
            # old format
            if isinstance(x, FileType):
                return [x]
            # New Format [(FileType, label, desc)]
            if isinstance(x[0], (list, tuple)):
                if isinstance(x[0][0], FileType):
                    return [x[0][0]]
    _raise()


def _validate_schema_options(x):
    _raise = _raise_malformed_task_attr("Schema options must be provided as DI list or dict. Got type {t}.".format(t=type(x)))

    # Standard form
    if isinstance(x, dict):
        if x:
            is_valid = _is_schema_list(x)
            if is_valid:
                return x
        else:
            # emtpy dict
            return x
    elif isinstance(x, (list, tuple)):
        x = to_list_if_necessary(x)
        if not __is_schema(x[0]):
            _raise("When task options are provided as DI model list, the first item must be the schema options for the Task.")
        return x

    # All other cases fail
    _raise()


def _validate_nproc(x):

    msg = "NPROC ('{x}') must be a DI LIST or primitive int value, or {s}".format(x=x, s=SymbolTypes.MAX_NPROC)
    _raise = _raise_malformed_task_attr(msg)
    if isinstance(x, int):
        if x > 0:
            return x
        _raise("NPROC must be > 0")
    elif isinstance(x, str):
        if x == SymbolTypes.MAX_NPROC:
            return x
        else:
            _raise("")
    elif isinstance(x, (tuple, list)):
        # Validate DI
        return to_list_if_necessary(x)
    else:
        _raise("Got type (t)").format(t=type(x))

    return x


def _validate_task_id(x):
    if isinstance(x, str):
        if RX_TASK_ID.match(x):
            return x
    else:
        raise MalformedMetaTaskError("Task id '{n}' must match {p}".format(p=RX_TASK_ID.pattern, n=x))


def _validate_resource_types(x):
    if x is None:
        return ()

    _raise = _raise_malformed_task_attr("Resource types must be a list with valid values {x}.".format(x=ResourceTypes.ALL()))

    if isinstance(x, (list, tuple)):
        for i in x:
            if i not in ResourceTypes.ALL():
                _raise("Invalid resource value '{x}'".format(x=i))
        return x
    else:
        _raise("Invalid type {x}".format(x=x))


def _validate_version(x):
    m = RX_VERSION.match(x)
    if m is None:
        _raise_malformed_task_attr("Version '{v}' should match {p}".format(v=x, p=RX_VERSION.pattern))
    return x


def _validate_nchunks(x):
    if isinstance(x, str):
        if x == SymbolTypes.MAX_NCHUNKS:
            return x
    if isinstance(x, int):
        return x
    if isinstance(x, (list, tuple)):
        if isinstance(list(x)[-1], types.FunctionType):
            # Add more validation here
            return x

    msg = "Chunk only supports int or {x}, or DI model list".format(x=SymbolTypes.MAX_NCHUNKS)
    _raise_malformed_task_attr(msg)


def _validate_chunk_keys(chunk_keys):
    if isinstance(chunk_keys, (list, tuple)):
        if not chunk_keys:
            _raise_malformed_task_attr("CHUNK_KEYS can NOT be empty.")
        for chunk_key in chunk_keys:
            if RX_CHUNK_KEY.match(chunk_key) is None:
                _raise_malformed_task_attr("CHUNK_KEYS '{x}' must match pattern {p}".format(x=chunk_key, p=RX_CHUNK_KEY))
        return chunk_keys

    _raise_malformed_task_attr("CHUNK_KEYS '{m}' is malformed".format(m=chunk_keys))

