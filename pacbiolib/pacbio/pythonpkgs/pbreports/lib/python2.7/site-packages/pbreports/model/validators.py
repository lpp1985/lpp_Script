import logging
from pprint import pformat
import functools

from jsonschema import validate

log = logging.getLogger(__name__)


def _validate(schema, dct):
    log.debug("Validating\n" + pformat(dct))
    log.debug("with schema\n" + pformat(schema))
    # this will raise an ValidationError
    validate(dct, schema)
    log.debug("Successfully validated")


def _id_schema():
    return {'type': 'string'}


def _get_attribute_schema():
    # value can be anything?
    schema = {"$schema": "http://json-schema.org/draft-04/schema#",
              'type': 'object',
              'properties': {'id': _id_schema(),
                             'name': {'type': 'string', 'optional': True},
                             'value': {'type': ['string', 'number']}
                             },
              'required': ['id', 'value'],
              "additionalProperties": False}

    return schema


def _get_plot_schema():
    schema = {'$schema': 'http://json-schema.org/draft-04/schema#',
              'type': 'object',
              'properties': {'id': _id_schema(),
                             'name': {'type': 'string', 'optional': True},
                             },
              'required': ['id']}

    return schema


def _get_plot_group_schema():
    schema = {'$schema': 'http://json-schema.org/draft-04/schema#',
              'type': 'object',
              'properties': {'id': _id_schema(),
                             'title': {'type': 'string'},
                             'legend': {'type': 'string'},
                             'thumbnail': {'type': 'string'},
                             'plots': {'type': 'array',
                                       'items': _get_plot_schema()},
                             },
              'required': ['id', 'title']}

    return schema


def _get_column_schema():
    # the values can be anything?
    schema = {'$schema': 'http://json-schema.org/draft-04/schema#',
              'type': 'object',
              'properties': {'id': _id_schema(),
                             'header': {'type': 'string'}},
              'required': ['id', 'values']}
    return schema


def _get_table_schema():
    schema = {'$schema': 'http://json-schema.org/draft-04/schema#',
              'type': 'object',
              'properties': {'id': _id_schema(),
                             'title': {'type': 'string'},
                             'columns': {'type': 'array',
                                         'items': _get_column_schema()}},
              'required': ['id', 'columns']}

    return schema


def _get_report_schema():
    schema = {'$schema': 'http://json-schema.org/draft-04/schema#',
              'type': 'object',
              'properties': {'id': _id_schema(),
                             'tables': {'type': 'array',
                                        'items': _get_table_schema()},
                             'attributes': {'type': 'array',
                                            'items': _get_attribute_schema()},
                             'plotGroups': {'type': 'array',
                                            'items': _get_plot_group_schema()}},
              'required': ['id']}

    return schema


validate_attribute = functools.partial(_validate, _get_attribute_schema())
validate_plot = functools.partial(_validate, _get_plot_schema())
validate_plot_group = functools.partial(_validate, _get_plot_group_schema())
validate_column = functools.partial(_validate, _get_column_schema())
validate_table = functools.partial(_validate, _get_table_schema())
validate_report = functools.partial(_validate, _get_report_schema())
