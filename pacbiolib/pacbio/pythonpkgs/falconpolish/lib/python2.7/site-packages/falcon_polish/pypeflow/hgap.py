from __future__ import absolute_import
from . import flow
from ..functional import stricter_json
import copy
import json
import logging
import logging.config
import os
import pprint
import sys
import time
import ConfigParser as configparser
import StringIO
log = logging.getLogger(__name__)


default_logging_config = """
[loggers]
keys=root,pypeflow,hgap

[handlers]
keys=stream,file_pypeflow,file_hgap

[formatters]
keys=form01,form02

[logger_root]
level=NOTSET
handlers=stream

[logger_pypeflow]
level=DEBUG
handlers=file_pypeflow
qualname=pypeflow
propagate=1

[logger_hgap]
level=NOTSET
handlers=file_hgap
qualname=.
propagate=1

[handler_stream]
class=StreamHandler
level=INFO
formatter=form02
args=(sys.stderr,)

[handler_file_pypeflow]
class=FileHandler
level=INFO
formatter=form01
args=('pypeflow.log',)

[handler_file_hgap]
class=FileHandler
level=DEBUG
formatter=form01
args=('hgap.log',)

[formatter_form01]
format=%(asctime)s - %(name)s - %(levelname)s - %(message)s

[formatter_form02]
format=[%(levelname)s]%(message)s
"""

dl = \
{
    'version': 1,
    'disable_existing_loggers': True,
    #'incremental': True,
    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
        },
    },
    'handlers': {
        'default': {
            'level': 'INFO',
            'class': 'logging.StreamHandler',
        },
    },
    'loggers': {
        '': {
            'handlers': ['default'],
            'level': 'INFO',
            'propagate': True
        },
        'django.request': {
            'handlers': ['default'],
            'level': 'WARN',
            'propagate': False
        },
    }
}
def setup_logger(logging_config_fn):
    """See https://docs.python.org/2/library/logging.config.html
    """
    logging.Formatter.converter = time.gmtime # cannot be done in .ini

    #logging.config.dictConfig(dl)
    #return
    if logging_config_fn:
        if logging_config_fn.endswith('.json'):
            logging.config.dictConfig(json.loads(stricter_json(open(logging_config_fn).read())))
            #print repr(logging.Logger.manager.loggerDict) # to debug
            return
        logger_fileobj = open(logging_config_fn)
    else:
        logger_fileobj = StringIO.StringIO(default_logging_config)
    defaults = {
    }
    logging.config.fileConfig(logger_fileobj, defaults=defaults, disable_existing_loggers=False)

def cfg2dict(ifp):
    cp = configparser.ConfigParser()
    cp.readfp(ifp)
    return {section: cp.items(section) for section in cp.sections()}

DEFAULT_OPTIONS = """
{
  "hgap": {
    "GenomeSize": 8000,
    "min_length_cutoff": 1,
    "job_type": "local",
    "job_queue": "DEFAULT_JOB_QUEUE",
    "use_tmpdir": "true",
    "~comment": "Overrides for full HGAP pipeline"
  },
  "falcon": {
    "falcon_sense_option": "--output_multi --min_idt 0.77 --min_cov 10 --max_n_read 2000 --n_core 6",
    "genome_size": "48500",
    "length_cutoff": "-1",
    "length_cutoff_pr": "50",
    "overlap_filtering_setting": "--max_diff 1000 --max_cov 100000 --min_cov 0 --bestn 1000 --n_core 4",
    "ovlp_DBsplit_option": "-s50 -a",
    "ovlp_HPCdaligner_option": "-v -k15 -h60 -w6 -e.95 -l40 -s100 -M16",
    "pa_DBsplit_option": "-x250 -s500 -a",
    "pa_HPCdaligner_option": "-v -k16 -h35 -w7 -e.70 -l40 -s100 -M16",
    "seed_coverage": "74",
    "~comment": "Overrides for FALCON"
  },
  "pbalign": {
    "options": "--hitPolicy randombest --minAccuracy 70.0 --minLength 50 --algorithm=blasr",
    "algorithmOptions": "--minMatch 12 --bestn 10 --minPctSimilarity 70.0",
    "~jdnotes": "--maxHits 1 --minAnchorSize 12 --maxDivergence=30 --minAccuracy=0.75 --minLength=50 --hitPolicy=random --seed=1",
    "~comment": "Overrides for blasr alignment (prior to polishing)"
  },
  "variantCaller": {
    "options": "--algorithm arrow --minConfidence 40 --minCoverage 5",
    "~comment": "Overrides for genomic consensus (polishing)"
  },
  "pbcoretools.tasks.filterdataset": {
    "other_filters": "rq >= 0.7, length lte 50000",
    "read_length": 1
  },
  "pbreports.tasks.summarize_coverage": {
    "options": "--num_regions 1000 --region_size 0",
    "~comment": "--force_num_regions"
  },
  "pbcoretools.task_options": {
    "scatter_subread_max_nchunks": "5",
    "scatter_alignments_reference_max_nchunks": "12",
    "~comment": "Overrides for pbcoretools task_options (mainly for scatter/gather)"
  },
  "pbsmrtpipe": {
    "~comment": "Overrides for pbsmrtpipe"
  },
  "~comment": "https://github.com/PacificBiosciences/ExperimentalPipelineOptionsDocs/HGAP",
}
"""
DEFAULT_OPTIONS = json.loads(stricter_json(DEFAULT_OPTIONS))

falcon_lambda_defaults = """
[General]
falcon_sense_option = --output_multi --min_idt 0.77 --min_cov 10 --max_n_read 2000 --n_core 6
length_cutoff = 1
length_cutoff_pr = 1
overlap_filtering_setting = --max_diff 1000 --max_cov 100000 --min_cov 0 --bestn 1000 --n_core 4
ovlp_DBsplit_option = -s50 -a
ovlp_HPCdaligner_option = -v -k15 -h60 -w5 -H1 -e.95 -l40 -s100 -M4
pa_DBsplit_option = -x250 -s500 -a
pa_HPCdaligner_option = -v -k15 -h35 -w5 -H1 -e.70 -l40 -s100 -M4
"""
# These are used in the passing siv test for lambda, same subread input! But GL=5M
falcon_test_options = """
[General]
#pbsmrtpipe.hgap_coresmax_str = 40
#pbsmrtpipe.hgap_genomelength_str = 5000000
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 1 --max_n_read 20000 --n_core 6
length_cutoff = 1
length_cutoff_pr = 1
overlap_filtering_setting = --max_diff 10000 --max_cov 100000 --min_cov 0 --bestn 1000 --n_core 4
ovlp_dbsplit_option = -x500 -s50
ovlp_hpcdaligner_option = -v -dal4 -t32 -h60 -e.96 -l500 -s1000
pa_dbsplit_option = -x500 -s50
pa_hpcdaligner_option = -v -dal4 -t16 -e.70 -l1000 -s1000
"""

def update2(start, updates):
    """Update 'start' dict with 'updates' to 2 levels.
    Skip '~comment' entries.
    """
    for key1, val1 in updates.iteritems():
        if key1.startswith('~'):
            continue
        if key1 not in start:
            start[key1] = copy.deepcopy(val1)
            continue
        assert isinstance(start[key1], dict), "%s %s" %(repr(key1), repr(start))
        for key2, val2 in val1.iteritems():
            if key2.startswith('~'):
                continue
            start[key1][key2] = copy.deepcopy(val2)

def run(input_config_fn, logging_config_fn):
    #global log
    #import logging_tree
    #logging_tree.printout()
    setup_logger(logging_config_fn)
    #logging_tree.printout()
    #log = logging.getLogger(__name__)
    log.info('Read logging config from {!r}'.format(logging_config_fn))
    log.info('Reading HGAP config from {!r}'.format(input_config_fn))
    if input_config_fn.endswith('.json'):
        config = json.loads(stricter_json(open(input_config_fn).read()))
    else:
        config = cfg2dict(open(input_config_fn))
    log.info('config=\n{}'.format(pprint.pformat(config)))
    log.info('defaults=\n{}'.format(pprint.pformat(DEFAULT_OPTIONS)))
    opts = dict()
    opts.update(DEFAULT_OPTIONS)
    update2(opts, config)
    #opts['falcon'].update(cfg2dict(StringIO.StringIO(falcon_test_options))['General']) # TODO: gen_config from GenomeLength!
    log.info('opts=\n{}'.format(pprint.pformat(opts)))
    flow.flow(opts)
