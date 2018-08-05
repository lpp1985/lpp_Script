import os

import functools

import avro.schema
from avro.io import validate

SCHEMA_REGISTRY = {}

__all__ = ['validate_pipeline_template',
           'validate_pipeline_template_view_rules',
           'validate_tc',
           'validate_rtc',
           'SCHEMA_REGISTRY']


def _load_schema(idx, name):

    d = os.path.dirname(__file__)
    schema_path = os.path.join(d, name)
    with open(schema_path, 'r') as f:
        schema = avro.schema.parse(f.read())
    SCHEMA_REGISTRY[idx] = schema
    return schema

PT_SCHEMA = _load_schema("pipeline_template", "pipeline_template.avsc")
PTVR_SCHEMA = _load_schema("pipeline_template_view_rules", "pipeline_template_view_rules.avsc")


def _validate(schema, d):
    """Validate a python dict against a avro schema"""
    return validate(schema, d)

validate_pipeline_template = functools.partial(_validate, PT_SCHEMA)
validate_pipeline_template_view_rules = functools.partial(_validate, PTVR_SCHEMA)