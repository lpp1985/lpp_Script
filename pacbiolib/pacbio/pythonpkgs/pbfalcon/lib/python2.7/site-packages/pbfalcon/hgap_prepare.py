"""Generate FALCON cfg (.ini file).

We plan to generate cfg in a complicated way.
But for now, we just use a look-up table,
based on ranges of the length of a genome.
"""
from __future__ import absolute_import
from falcon_polish.functional import stricter_json
from falcon_polish.sys import symlink
from falcon_polish.pypeflow.hgap import update2
from falcon_kit import run_support as support
from . import gen_config # for some option names
import collections
import commands
import ConfigParser as configparser
import copy
import json
import logging
import os
import pprint
import re
import StringIO
import sys


#logging.basicConfig()
log = logging.getLogger(__name__)
#log.setLevel(logging.DEBUG)

TASK_HGAP_OPTIONS = \
        'falcon_ns.task_options.' + 'HGAP_Options_JSON'
TASK_HGAP_GENOME_LENGTH = \
        'falcon_ns.task_options.' + gen_config.OPTION_GENOME_LENGTH
TASK_HGAP_SEED_LENGTH_CUTOFF = \
        'falcon_ns.task_options.' + gen_config.OPTION_SEED_LENGTH_CUTOFF
TASK_HGAP_SEED_COVERAGE = \
        'falcon_ns.task_options.' + gen_config.OPTION_SEED_COVERAGE

DEFAULT_LOGGING_CFG = {
    'version': 1,
    'formatters': {
        'format_full': {
            'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        },
        'format_brief': {
            'format': '%(levelname)s: %(message)s',
        }
    },
    'filters': {
    },
    'handlers': {
        'handler_file_all': {
            'class': 'logging.FileHandler',
            'level': 'DEBUG',
            'formatter': 'format_full',
            'filename': 'all.log',
            'mode': 'w',
        },
        'handler_stream': {
            'class': 'logging.StreamHandler',
            'level': 'INFO',
            'formatter': 'format_brief',
            'stream': 'ext://sys.stderr',
        },
    },
    'loggers': {
    },
    'root': {
        'handlers': ['handler_stream', 'handler_file_all'],
        'level': 'NOTSET',
    },
    'disable_existing_loggers': False
}
"""
    'handlers': {
        'handler_file_pypeflow': {
            'class': 'logging.FileHandler',
            'level': 'INFO',
            'formatter': 'format_full',
            'filename': 'pypeflow.log',
            'mode': 'w',
        },
        'handler_file_pbfalcon': {
            'class': 'logging.FileHandler',
            'level': 'DEBUG',
            'formatter': 'format_full',
            'filename': 'pbfalcon.log',
            'mode': 'w',
        }
    },
    'loggers': {
        'pypeflow': {
            'level': 'NOTSET',
            'propagate': 1,
            'handlers': ['handler_file_pypeflow'],
        },
        'pbfalcon': {
            'level': 'NOTSET',
            'propagate': 1,
            'handlers': ['handler_file_pbfalcon'],
        },
    },
"""
OPTION_SECTION_HGAP = 'hgap'
OPTION_SECTION_FALCON = 'falcon'
OPTION_SECTION_PBALIGN = 'pbalign'
OPTION_SECTION_VARIANTCALLER = 'variantcaller'
OPTION_SECTION_PBSMRTPIPE = 'pbsmrtpipe'

def say(msg):
    sys.stderr.write('Im just sayin:' + msg + '\n') # since pbsmrtpipe suppresses our logging
    log.info(msg)

def get_pbsmrtpipe_opts(d):
    with open(os.path.join(d, 'resolved-tool-contract.json')) as f:
        rtc = json.loads(f.read())
    opts = rtc['resolved_tool_contract']
    assert 'nproc' in opts
    assert 'is_distributed' in opts
    assert 'resources' in opts
    # TODO: Removed any unneeded rtc opts.
    return opts

def dump_as_json(data, ofs):
    as_json = json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))
    say('Dumping JSON:\n{}'.format(
        as_json))
    ofs.write(as_json)

def learn_submit_template(content):
    """Infer the submission template from the pbsmrtpipe runnable-task.json file.
    """
    table = json.loads(content)
    cluster_table = table['cluster']
    if cluster_table is None:
        raise RuntimeError('No cluster information in "{}"'.format(content))
    start_string = cluster_table['start']
    stop_string = cluster_table['stop']
    return start_string, stop_string

def update_for_grid(all_cfg, run_dir):
    """Update job_queue,
    based on pbsmrtpipe runnable-task.json file.
    """
    job_queue = all_cfg[OPTION_SECTION_HGAP].get('job_queue')
    # Set job_queue. (job_type means nothing now.)
    if not job_queue:
        runnable_task_fn = os.path.join(run_dir, 'runnable-task.json')
        try:
            # Get start_tmpl from pbsmrtpipe.
            with open(runnable_task_fn) as ifs:
                submit_template, kill_template = learn_submit_template(ifs.read())
            job_queue = submit_template
        except Exception:
            # For testing, this is not really an error. But customers should never see this code-path.
            log.exception('Running locally instead of via cluster.')
            job_queue = '/bin/bash -c "${CMD}"'
    all_cfg[OPTION_SECTION_HGAP]['job_queue'] = job_queue

def update_pwatcher(all_cfg):
    # Set defaults, in case of job_type=sge?
    #fc_cfg = all_cfg[OPTION_SECTION_FALCON]
    #sge_option_names = (
    #        'sge_option_da', 'sge_option_la',
    #        'sge_option_pda', 'sge_option_pla',
    #        'sge_option_fc', 'sge_option_cns',
    #)
    #for option_name in sge_option_names:
    #    fc_cfg[option_name] = sge_queue_option

    pwatcher_type = all_cfg[OPTION_SECTION_HGAP].get('pwatcher_type', 'blocking')
    all_cfg[OPTION_SECTION_HGAP]['pwatcher_type'] = pwatcher_type
    if pwatcher_type in ('fs_based', 'network_based'):
        log.error('pwatcher_type {!r} not currently supported within pbsmrtpipe. It may fail when run with an installed tarball.'.format(
            pwatcher_type))
        job_type = 'sge'
    else:
        job_type = all_cfg[OPTION_SECTION_HGAP].get('job_type', 'string') # for info only
    all_cfg[OPTION_SECTION_HGAP]['job_type'] = job_type
    all_cfg[OPTION_SECTION_FALCON]['job_type'] = job_type # I think this is used per-task, maybe.

def update_falcon(all_cfg):
    # Propagate use_tmpdir from hgap section to falcon section.
    use_tmpdir = all_cfg[OPTION_SECTION_HGAP].get('use_tmpdir')
    if use_tmpdir:
        all_cfg[OPTION_SECTION_FALCON]['use_tmpdir'] = use_tmpdir
    # Set concurrency per stage. Just use the global default for now.
    max_conc = 24
    all_cfg[OPTION_SECTION_FALCON].setdefault('cns_concurrent_jobs', max_conc)
    all_cfg[OPTION_SECTION_FALCON].setdefault('ovlp_concurrent_jobs', max_conc)
    all_cfg[OPTION_SECTION_FALCON].setdefault('pa_concurrent_jobs', max_conc)

def run_hgap_prepare(input_files, output_files, options):
    """Generate a config-file from options.
    """
    say('options to run_hgap_prepare:\n{}'.format(pprint.pformat(options)))
    i_subreadset_fn, = input_files
    o_hgap_cfg_fn, o_logging_cfg_fn, o_log_fn = output_files
    run_dir = os.path.dirname(o_hgap_cfg_fn)
    symlink(os.path.join(run_dir, 'stderr'), o_log_fn)

    # This will be the cfg we pass to hgap_run.
    all_cfg = collections.defaultdict(lambda: collections.defaultdict(str))

    # Get grid options, for job-distribution.
    update_for_grid(all_cfg, run_dir)

    # Set some other falcon options, based on hgap options.
    update_falcon(all_cfg)

    # Override from pbsmrtpipe config/preset.xml.
    all_cfg[OPTION_SECTION_FALCON]['genome_size'] = options[TASK_HGAP_GENOME_LENGTH].strip()
    all_cfg[OPTION_SECTION_FALCON]['length_cutoff'] = options[TASK_HGAP_SEED_LENGTH_CUTOFF].strip()
    all_cfg[OPTION_SECTION_FALCON]['seed_coverage'] = options[TASK_HGAP_SEED_COVERAGE].strip()
    cfg_json = options[TASK_HGAP_OPTIONS].strip()
    if not cfg_json:
        cfg_json = '{}'
    override_cfg = json.loads(stricter_json(cfg_json))
    update2(all_cfg, override_cfg)

    update_pwatcher(all_cfg)

    # Get options from pbsmrtpipe.
    pbsmrtpipe_opts = get_pbsmrtpipe_opts(run_dir)
    if OPTION_SECTION_PBSMRTPIPE not in all_cfg:
        all_cfg[OPTION_SECTION_PBSMRTPIPE] = dict()
    pbsmrtpipe_opts.update(all_cfg[OPTION_SECTION_PBSMRTPIPE])
    all_cfg[OPTION_SECTION_PBSMRTPIPE] = pbsmrtpipe_opts

    # Dump all_cfg.
    say('Dumping to {}'.format(repr(o_hgap_cfg_fn)))
    dump_as_json(all_cfg, open(o_hgap_cfg_fn, 'w'))

    # Get logging cfg.
    logging_cfg = DEFAULT_LOGGING_CFG

    # Dump logging cfg.
    say('Dumping to {}'.format(repr(o_logging_cfg_fn)))
    dump_as_json(logging_cfg, open(o_logging_cfg_fn, 'w'))
