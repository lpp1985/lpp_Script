"""DI Model utils used for Manipulating TaskOptions (in JSONSchema)"""
import logging

import jsonschema
from pbsmrtpipe.constants import to_task_option_ns


log = logging.getLogger(__name__)

__all__ = ['to_opt_id', 'is_valid', 'to_option_schema']


def to_opt_id(s):
    return to_task_option_ns(s)


def validate_value(schema, v):
    return jsonschema.validate(v, schema)


def is_valid(schema, v):
    """Returns a bool if the schema is valid"""
    try:
        validate_value(schema, v)
        return True
    except jsonschema.ValidationError:
        pass
    return False


def get_default_from_schema(schema):
    option_id = schema['properties'].keys()[0]
    return schema['properties'][option_id]['default']


def validate_schema(f):
    """Deco for validate the returned jsonschema against Draft 4 of the spec"""
    def w(*args, **kwargs):
        schema = f(*args, **kwargs)
        x = jsonschema.Draft4Validator(schema)
        return schema

    w.__doc__ = f.__doc__
    w.__name__ = f.__name__

    return w


@validate_schema
def to_option_schema(option_id, dtype_or_dtypes, display_name, description, default_value):
    """
    Simple util factory method


    :param dtype_or_dtypes: single data type or list of data types
    :param option_id: globally unique task option id. Must begin with
    'pbsmrtpipe.task_options.'
    :param display_name: display name of task options
    :param description: Short description of the task options
    :param required: Is the option required.
    """
    # annoying that you can't specify a tuple
    if isinstance(dtype_or_dtypes, tuple):
        dtype_or_dtypes = list(dtype_or_dtypes)

    d = {'$schema': "http://json-schema.org/draft-04/schema#",
         'type': 'object',
         'title': "JSON Schema for {o}".format(o=option_id),
         'properties': {option_id: {'description': description,
                                    'title': display_name,
                                    'type': dtype_or_dtypes}}
         }

    d['required'] = [option_id]
    d['properties'][option_id]['default'] = default_value
    return d


def crude_coerce_type_from_str(s, type_or_types):
    """
    Convert the type from a raw string to a specific type.

    :param s: raw string
    :param type_or_types: List or single value of option type (as a string) to coerce Example "bool"
    """

    def _to_bool(x):
        _d = {'true': True, 'false': False}
        if isinstance(x, str):
            if x.lower() in _d:
                return _d[x.lower()]
            else:
                raise TypeError
        raise TypeError

    # need to better handle null here
    _coerce = {'integer': int,
               'number': float,
               'null': lambda x: None,
               "string": str,
               "boolean": _to_bool}

    if not isinstance(type_or_types, (list, tuple)):
        type_or_types = [type_or_types]

    # Make sure Null is the last item processed
    if 'null' in type_or_types:
        i = type_or_types.index('null')
        type_or_types.pop(i)
        type_or_types.append('null')

    for type_ in type_or_types:
        if type_ in _coerce:
            try:
                f = _coerce[type_]
                return f(s)
            except TypeError:
                pass

    return s
