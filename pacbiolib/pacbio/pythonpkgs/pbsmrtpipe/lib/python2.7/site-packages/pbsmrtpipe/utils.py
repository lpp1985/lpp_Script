"""
General func/tools used in pbsmrtpipe
"""
import os
import sys
import re
import time
import logging
import logging.config
import logging.handlers

import functools

from jinja2 import Environment, PackageLoader

from pbcore.util.Process import backticks

# for backward compatibility
from pbcommand.utils import setup_log, compose, nfs_exists_check, nfs_refresh
from pbcommand.validators import (validate_file, validate_dir, validate_fofn, validate_output_dir, fofn_to_files)

from pbsmrtpipe.decos import ignored
from pbsmrtpipe.constants import SLOG_PREFIX

HTML_TEMPLATE_ENV = Environment(loader=PackageLoader('pbsmrtpipe', 'html_templates'))


log = logging.getLogger(__name__)
slog = logging.getLogger(SLOG_PREFIX + __name__)


def validate_type_or_raise(obj, klasses, msg=None):
    if not isinstance(obj, klasses):
        emsg = "{o} Got type {x}, expected type {y}.".format(o=obj, x=type(obj), y=klasses)
        if msg is not None:
            emsg = " ".join([emsg, msg])
        raise TypeError(emsg)
    return obj


def log_timing(func):
    """Simple deco to log the runtime of func"""
    started_at = time.time()

    def wrapper(*args, **kw):
        return func(*args, **kw)

    run_time = time.time() - started_at
    name = func.__name__
    log.info("Func {f} took {s:.2f} sec ({m:.2f} min)".format(f=name, s=run_time, m=run_time / 60.0))

    return wrapper


class StdOutStatusLogFilter(logging.Filter):

    def filter(self, record):
        return record.name.startswith(SLOG_PREFIX)


def is_verified(path, max_nfs_refresh=3):
    """Validate that a file exists. Force NFS refresh if necessary"""
    for i in xrange(max_nfs_refresh):
        with ignored(OSError):
            # Try to force an NFS refresh
            os.listdir(os.path.dirname(path))
            if os.path.exists(path):
                return True
            # time.sleep(0.25)

    return False


def get_default_logging_config_dict(master_log, master_level, pb_log, stdout_level):
    """Returns a dict configuration of the logger. """
    d = {
        'version': 1,
        'disable_existing_loggers': False,  # this fixes the problem
        'formatters': {
            'console': {
                'format': '%(message)s'
            },
            'standard': {
                'format': '[%(levelname)s] %(asctime)-15sZ [%(name)s] %(message)s'
            },
            'full': {
                'format': '[%(levelname)s] %(asctime)-15sZ [%(name)s %(funcName)s %(lineno)d] %(message)s'
            }
        },
        'filters': {
            "slog_filter": {
                '()': StdOutStatusLogFilter,
            }
        },
        'handlers': {
            'console': {
                'level': logging.getLevelName(stdout_level),
                'class': 'logging.StreamHandler',
                'formatter': 'console',
                'stream': 'ext://sys.stdout',
                'filters': ['slog_filter']
            },
            "debug_file_handler": {
                "class": 'logging.handlers.RotatingFileHandler',
                "level": logging.getLevelName(master_level),
                "formatter": "full",
                "filename": master_log,
                "maxBytes": "10485760",
                "backupCount": "20",
                "encoding": "utf8"
            },
            "info_file_handler": {
                "class": 'logging.handlers.RotatingFileHandler',
                "level": "INFO",
                "formatter": "standard",
                "filename": pb_log,
                "maxBytes": "10485760",
                "backupCount": "20",
                "encoding": "utf8",
                "filters": ['slog_filter']
            }
        },
        'loggers': {
            '': {
                'handlers': ['console', 'info_file_handler', 'debug_file_handler'],
                'level': 'DEBUG',
                'propagate': True
            }
        },
        'root': {
            'level': 'DEBUG',
            'handlers': ['console', 'debug_file_handler', 'info_file_handler']
        }

    }
    return d


def setup_internal_logs(master_log, master_level, pb_log, stdout_level):
    d = get_default_logging_config_dict(master_log, master_level, pb_log, stdout_level)
    logging.config.dictConfig(d)
    logging.Formatter.converter = time.gmtime
    return d
