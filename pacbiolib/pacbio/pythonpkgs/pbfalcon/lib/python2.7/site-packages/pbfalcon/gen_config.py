"""Generate FALCON cfg (.ini file).

We plan to generate cfg in a complicated way.
But for now, we just use a look-up table,
based on ranges of the length of a genome.
"""
from falcon_kit import run_support as support
from . import tusks
import ConfigParser as configparser
import logging
import os
import re
import StringIO


#logging.basicConfig()
log = logging.getLogger(__name__)
#log.setLevel(logging.DEBUG)
OPTION_GENOME_LENGTH = 'HGAP_GenomeLength_str'
OPTION_SEED_COVERAGE = 'HGAP_SeedCoverage_str'
OPTION_SEED_LENGTH_CUTOFF = 'HGAP_SeedLengthCutoff_str'
OPTION_CORES_MAX = 'HGAP_CoresMax_str'
OPTION_CFG = 'HGAP_FalconAdvanced_str'
OPTION_AGGRESSIVE_ASM = 'HGAP_AggressiveAsm_bool'

# Override pa_hpcdaligner_option if aggressive (greedy) mode is on
OPTION_HPC = 'pa_hpcdaligner_option'
AGGRESSIVE_HPC_OPTION_VALUE = "-v -dal24 -t14 -h70 -e.70 -l1000 -s100 -k14"

defaults_old = """\
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 1 --max_n_read 20000 --n_core 6
length_cutoff = 1000
length_cutoff_pr = 1
pa_DBsplit_option = -x5 -s50 -a
pa_HPCdaligner_option =  -v -k25 -h35 -w5 -e.95 -l40 -s1000 -t27
overlap_filtering_setting = --max_diff 10000 --max_cov 100000 --min_cov 0 --bestn 1000 --n_core 4
ovlp_HPCdaligner_option =  -v -k25 -h35 -w5 -e.99 -l40 -s1000 -t27
ovlp_DBsplit_option = -x5 -s50 -a
falcon_sense_greedy = False
"""
old_defaults_lambda = """\
genome_size = 48000
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 6
length_cutoff = -1
length_cutoff_pr = 12000
pa_DBsplit_option = -x500 -s50 -a
pa_HPCdaligner_option = -v -dal4 -t16 -e.70 -l1000 -s1000
overlap_filtering_setting = --max_diff 100 --max_cov 50 --min_cov 1 --bestn 10 --n_core 24
ovlp_HPCdaligner_option = -v -dal4 -t32 -h60 -e.96 -l500 -s1000
ovlp_DBsplit_option = -x500 -s50 -a
falcon_sense_greedy = False
"""
# These values will need more adjusting, but at least they worked on some dataset.
# The current values are from my latest experiment, producing 2 contigs and a total
# draft assembly of 49820b. I meant to experiment some more, but we are in a rush to
# to update defaults in 3.0.4 for ecoli.
defaults_lambda = """
genome_size = 48502
seed_coverage = 80
falcon_sense_option = --output_multi --min_idt 0.77 --min_cov 10 --max_n_read 2000 --n_core 6
length_cutoff = -1
length_cutoff_pr = 50
overlap_filtering_setting = --max_diff 1000 --max_cov 100000 --min_cov 0 --bestn 1000 --n_core 4
ovlp_DBsplit_option = -s50 -a
ovlp_hpcdaligner_option = -v -k15 -h60 -w6 -e.95 -l40 -s100 -M16
pa_DBsplit_option = -x250 -s500 -a
pa_HPCdaligner_option =   -v -k15 -h35 -w7 -e.70 -l40 -s100 -M16
falcon_sense_greedy = False
"""
# ecoli based on Jim's run: http://smrtlink-beta:8080/#/analysis-job/2437
defaults_ecoli = """
genome_size = 4500000
seed_coverage = 27
length_cutoff = -1
length_cutoff_pr = 500
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 6
overlap_filtering_setting = --max_diff 60 --max_cov 100 --min_cov 4 --bestn 10 --n_core 4
ovlp_dbsplit_option = -x500 -s200 -a
ovlp_hpcdaligner_option = -v -dal24 -t16 -h35 -e.93 -l1000 -s100 -k25
pa_dbsplit_option =   -x500 -s200 -a
pa_hpcdaligner_option =   -v -dal24 -t14 -h70 -e.75 -l1000 -s100 -k18
falcon_sense_greedy = False
"""
defaults_yeast = """
genome_size = 12000000
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 8
length_cutoff = -1
length_cutoff_pr = 500
overlap_filtering_setting = --max_diff 40 --max_cov 80 --min_cov 2 --n_core 12
ovlp_DBsplit_option = -x15000 -s40
ovlp_HPCdaligner_option =  -v -dal4 -k24 -e.96  -s200 -M16 -l2500 -h1024
pa_DBsplit_option = -a -x500 -s100
pa_HPCdaligner_option =    -v -dal4 -k18 -e0.70 -s200 -M16 -l4800 -h480 -w8
falcon_sense_greedy = False
"""
defaults_human = """
genome_size = 3000000000
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 8
length_cutoff = -1
length_cutoff_pr = 500
overlap_filtering_setting = --max_diff 40 --max_cov 80 --min_cov 2 --n_core 12
ovlp_DBsplit_option = -x15000 -s40
ovlp_HPCdaligner_option =  -v -dal4 -k24 -e.96  -s200 -M16 -l2500 -h1024
pa_DBsplit_option = -a -x500 -s500
pa_HPCdaligner_option =    -v -dal4 -k18 -e0.70 -s200 -M16 -l4800 -h480 -w8
falcon_sense_greedy = False
"""
defaults_human_in_falcon = """
genome_size = 3000000000
length_cutoff = -1
length_cutoff_pr = 7000
pa_HPCdaligner_option =  -v -dal128 -t16 -e.70 -l1000 -s1000
ovlp_HPCdaligner_option = -v -dal128 -t32 -h60 -e.96 -l500 -s1000

pa_DBsplit_option = -x500 -s400
ovlp_DBsplit_option = -x500 -s400

falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 16

overlap_filtering_setting = --max_diff 60 --max_cov 60 --min_cov 2 --n_core 24
falcon_sense_greedy = False
"""
# See /lustre/hpcprod/jchin/CHM1_P6_project/asm3
defaults_human_recent = """
genome_size = 3000000000
length_cutoff = -1
length_cutoff_pr = 15000
pa_DBsplit_option = -a -x500 -s400
pa_HPCdaligner_option =  -v -dal128 -t16 -e0.75 -M24 -l4800 -k18 -h480 -w8 -s100
ovlp_DBsplit_option = -s400
ovlp_HPCdaligner_option =  -v -dal128  -M24 -k24 -h1024 -e.96 -l2500 -s100
overlap_filtering_setting = --max_diff 40 --max_cov 80 --min_cov 2 --n_core 12
falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 400 --n_core 12
falcon_sense_skip_contained = False
falcon_sense_greedy = False
"""
# also see:
#   https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide/
#   https://dazzlerblog.wordpress.com/2014/06/01/the-dazzler-db/
#   https://github.com/PacificBiosciences/FALCON/wiki/Manual
#   http://bugzilla.nanofluidics.com/show_bug.cgi?id=29491

defaults = list(sorted([
    (       0, defaults_old),
    ( 8*10**3, old_defaults_lambda),
    (10*10**3, defaults_lambda),
    ( 1*10**6, defaults_ecoli),
    #(10*10**6, defaults_yeast), # For now, use ecoli as default for bacteria.
    ( 1*10**9, defaults_human_recent), # These are old; maybe not good with Sequel.
]))


def sorted_str(s):
    return '\n'.join(sorted(s.splitlines()))

def _populate_falcon_options(options):
    length = int(options[OPTION_GENOME_LENGTH]) # required!
    index = 0
    # We could binary-search, but we will just walk thru.
    while index < len(defaults) - 1:
        if defaults[index+1][0] <= length:
            index += 1
        else:
            break
    fc = ini2dict(sorted_str(defaults[index][1]))

    # Also keep everything except a few which could be mal-formatted,
    # but prefix the ones from special pbsmrtpipe options.
    excluded = [OPTION_CFG]
    for key in options:
        if key not in excluded:
            fc['pbsmrtpipe.' + key] = options[key]

    # Translate some to FALCON options. (Note that names are different in Falcon.)
    fc['genome_size'] = int(options[OPTION_GENOME_LENGTH])
    fc['length_cutoff'] = int(options.get(OPTION_SEED_LENGTH_CUTOFF, '-9'))
    fc['seed_coverage'] = float(options.get(OPTION_SEED_COVERAGE, '21'))
    fc['falcon_sense_greedy'] = bool(options.get(OPTION_AGGRESSIVE_ASM, False))

    # Greedy mode, override pa_HPCdaligner_option too
    if fc['falcon_sense_greedy'] is True:
       log.info('Agressive (greedy) mode is turned on, override %s=%s' % (OPTION_HPC, AGGRESSIVE_HPC_OPTION_VALUE))
       fc[OPTION_HPC] = AGGRESSIVE_HPC_OPTION_VALUE # "-v -dal24 -t14 -h70 -e.70 -l1000 -s100 -k14"
    return fc

def _options_dict_with_base_keys(options_dict, prefix='falcon_ns.task_options.'):
    """Remove leading namespaces from key names,
    in a copy of options_dict.

    prefix: should include trailing dot
    """
    new_dict = dict()
    for key, val in options_dict.items():
        if key.startswith(prefix):
            tail = key[len(prefix):]
            if '.' in tail:
                log.warning('prefix {!r} found on option {!r}'.format(
                    prefix, key))
            new_dict[tail] = val
    return new_dict

def _gen_config(options_dict):
    """Generate ConfigParser object from dict.
    """
    cfg = configparser.ConfigParser()
    sec = "General"
    cfg.add_section(sec)
    for key, val in options_dict.items():
        # Strip leading and trailing ws, b/c the pbsmrtpipe
        # misinterprets XML as a data-interchanges language.
        # (It is only mark-up, so ws is never meaningful.)
        # Also, we want only strings; hopefully, we can fix
        # the TC later to drop the type-info (e.g. integer).
        cfg.set(sec, key, str(val).strip())
    return cfg

def _write_config(config, config_fn):
    with open(config_fn, 'w') as ofh:
        # I wish ConfigParser would sort. Oh, well.
        config.write(ofh)

def ini2dict(ini_text):
    ifp = StringIO.StringIO('[General]\n' + ini_text)
    cp = configparser.ConfigParser()
    cp.readfp(ifp)
    return dict(cp.items('General'))

re_semicolon = re.compile(r'\s*;+\s*')

def option_text2ini(option_text):
    # Basically, just translate semicolons into linefeeds.
    return re_semicolon.sub('\n', option_text)

re_newline = re.compile(r'\s*\n\s*', re.MULTILINE)

def ini2option_text(ini):
    # Basically, just translate linefeeds into semicolons.
    return re_newline.sub(';', ini)

def get_falcon_overrides(cfg_content, OPTION_CFG=OPTION_CFG):
    """options keys are bare (no namespaces)
    """
    if '\n' in cfg_content and '' != cfg_content.strip():
        log.error('linefeed found in option "%s", which is ok here but should have been prevented earler' %(
            OPTION_CFG))
        cfg_content = ini2option_text(cfg_content)
    if cfg_content.strip().startswith('['):
        log.error('Option "%s" seems to have .ini-style [brackets], which is an error. It is not really a .ini file.' %(
            OPTION_CFG))
        # Try to strip the first line, and hope there are no others.
        cfg_content = cfg_content[cfg_content.index(']'):]
    ini = option_text2ini(cfg_content)
    log.info(ini)
    # Now, parse the overrides, but skip it on any error.
    try:
        overrides = ini2dict(ini)
    except Exception as exc:
        log.exception('For option "%s" (for overrides) we had a problem parsing its contents:\n%s' %(
            OPTION_CFG, cfg_content))
        overrides = dict()
    return overrides

def run_falcon_gen_config(input_files, output_files, options):
    """Generate a config-file from options.
    """
    i_fofn_fn, = input_files
    o_cfg_fn, = output_files
    import pprint
    log.info('options to run_falcon_gen_config:\n{}'.format(pprint.pformat(options)))
    options = _options_dict_with_base_keys(options)
    falcon_options = _populate_falcon_options(options)
    log.info('falcon_options to run_falcon_gen_config:\n{}'.format(pprint.pformat(falcon_options)))
    if OPTION_CFG in options:
        overrides = get_falcon_overrides(options[OPTION_CFG], OPTION_CFG)
        log.info('overrides:\n%s'% pprint.pformat(overrides))
        falcon_options.update(overrides)
    else:
        raise Exception("Could not find %s" %OPTION_CFG)
    config = _gen_config(falcon_options)
    with tusks.cd(os.path.dirname(i_fofn_fn)):
        return _write_config(config, o_cfg_fn) # Write lower-case keys, which is fine.

