"""DataStore Validators"""

import os
import functools
import logging
import json

from pbcore.io import FastaReader, FastqReader
from pbcommand.models import FileTypes

log = logging.getLogger(__name__)

_REQUIRED_JSON_REPORT_KEYS = 'attributes id plotGroups tables'.split()


def _validate_json(path):
    """Simple smoke test to validate the json is well formed"""
    with open(path, 'r') as f:
        s = json.loads(f.read())
    return True


def _validate_json_report(path):
    """Smoke Test to make sure pbreports JSON file has the required
    root level keys.
    """
    with open(path, 'r') as f:
        s = json.loads(f.read())

    for key in _REQUIRED_JSON_REPORT_KEYS:
        if key not in s:
            _d = dict(s=key, p=path)
            msg = "Unable to find {s} in {p}".format(**_d)
            log.error(msg)
            raise KeyError(msg)

    return True


class DataStoreFileValidator(object):

    def __init__(self, file_type_id, func, **kwargs):
        self.file_type_id = file_type_id
        # Func must have signature of (path, **kwargs) => True|False\Raise exception
        self.func = func
        self.validator_func_kwargs = kwargs

    def __repr__(self):
        name = getattr(self.func, "name") if hasattr(self.func, "name") else str(self.func)
        _d = dict(k=self.__class__.__name__,
                  n=name,
                  w=self.validator_func_kwargs)
        return "<{k} name:{n} kw:{k} >".format(**_d)


def _validate_fastx(reader_klass, path):
    n = 0
    with reader_klass(path) as reader:
        for _ in reader:
            n += 1

    log.info("Successfully validated {k} records from {p}".format(k=reader_klass.__name__, p=path))
    return True

validate_fasta = functools.partial(_validate_fastx, FastaReader)
validate_fastq = functools.partial(_validate_fastx, FastqReader)


ValidateFasta = DataStoreFileValidator(FileTypes.FASTA.file_type_id, validate_fasta)
ValidateFastq = DataStoreFileValidator(FileTypes.FASTQ.file_type_id, validate_fastq)

ValidateJsonReport = DataStoreFileValidator(FileTypes.REPORT.file_type_id, _validate_json_report)
ValidateJson = DataStoreFileValidator(FileTypes.JSON.file_type_id, _validate_json)
