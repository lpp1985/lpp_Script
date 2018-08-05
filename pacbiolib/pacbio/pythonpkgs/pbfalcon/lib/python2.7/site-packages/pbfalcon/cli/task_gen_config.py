from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcommand.models import (ResourceTypes, FileTypes)
from pbcommand.models.parser import get_pbparser
from .. import gen_config
import sys
import logging

__version__ = '1.0.0'
log = logging.getLogger(__name__)
TOOL_ID = 'falcon_ns.tasks.task_falcon_gen_config'


def add_args_and_options(p):
    # FileType, label, name, description
    p.add_input_file_type(FileTypes.FOFN, "fofn_in", "FileOfFileNames", "FOFN for fasta files")
    # File Type, label, name, description, default file name
    p.add_output_file_type(FileTypes.CFG, "cfg_out", "INI File", "FALCON cfg (aka 'ini')", 'fc_run')
    # Option id, label, default value, name, description
    p.add_str("falcon_ns.task_options." + gen_config.OPTION_GENOME_LENGTH, "genome-length", '5000000',
            "Genome length", "Approx. number of base pairs expected in the genome. We choose many hidden settings automatically, based on this. (To learn what we generate, see fc_*.cfg, currently called 'falcon_ns.tasks.task_falcon0_build_rdb-PacBio.FileTypes.txt' amongst output files.)")
    p.add_str("falcon_ns.task_options." + gen_config.OPTION_SEED_COVERAGE, "seed-coverage", '30',
            "Seed coverage", "A target for the total # of bases in the 'raw' (post primary) reads, divided by the total number in the 'seed' reads.")
    p.add_str("falcon_ns.task_options." + gen_config.OPTION_SEED_LENGTH_CUTOFF, "seed-length-cutoff", '-1',
            "Seed length cutoff", "Only reads as long as this will be used as 'seeds' for the draft assembly. (Shorter reads will be used for correction and polishing, if they pass the dataset filters.) If '-1', then this will be calculated automatically, such that the total number of seed bases nearly equals GenomeLength*SeedCoverage.")
    #p.add_str("falcon_ns.task_options." + gen_config.OPTION_CORES_MAX, "cores-max", '40',
    #        "Cores Max IGNORED.", "IGNORE - not currently used")
    p.add_str("falcon_ns.task_options." + gen_config.OPTION_CFG, "falcon-overrides", '',
            "FALCON cfg overrides", "This is intended to allow support engineers to override the cfg which we will generate from other options. It is a semicolon-separated list of key=val pairs. Newlines are allowed but ignored. For more details on the available options, see https://github.com/PacificBiosciences/FALCON/wiki/Manual")
    p.add_boolean("falcon_ns.task_options." + gen_config.OPTION_AGGRESSIVE_ASM, "falcon-greedy", False,
            "Aggressive option", "Greedy option for LA4Falcon. For more details on the available options, see https://github.com/PacificBiosciences/FALCON/wiki/Manual")
    return p

def get_contract_parser():
    # Number of processors to use, can also be SymbolTypes.MAX_NPROC
    nproc = 1
    # Log file, tmp dir, tmp file. See ResourceTypes in models, ResourceTypes.TMP_DIR
    resource_types = ()
    # Commandline exe to call "{exe}" /path/to/resolved-tool-contract.json
    driver_exe = "python -m pbfalcon.cli.task_gen_config --resolved-tool-contract "
    desc = "Generate FALCON cfg from pbcommand options."
    name = 'Tool FalconConfigGenerator'
    p = get_pbparser(TOOL_ID, __version__, name, desc, driver_exe,
            is_distributed=False, nproc=nproc, resource_types=resource_types)
    add_args_and_options(p)
    return p

def run_my_main(input_files, output_files, options):
    # do stuff. Main should return an int exit code
    rc = gen_config.run_falcon_gen_config(input_files, output_files, options)
    if rc:
        return rc
    else:
        return 0

def _args_runner(args):
    # this is the args from parser.parse_args()
    # the properties of args are defined as "labels" in the add_args_and_options func.
    # TODO: Convert 'args' to a dict somehow?
    return run_my_main([args.fasta_in], [args.fasta_out], args)

def _resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    # all options are referenced by globally namespaced id. This allows tools to use other tools options
    # e.g., pbalign to use blasr defined options.
    return run_my_main(rtc.task.input_files, rtc.task.output_files, rtc.task.options)

def main(argv=sys.argv):
    log.info("Starting {f} version {v} pbcommand example dev app".format(f=__file__, v=__version__))
    p = get_contract_parser()
    return pbparser_runner(argv[1:],
                           p,
                           _args_runner, # argparse runner func
                           _resolved_tool_contract_runner, # tool contract runner func
                           log, # log instance
                           setup_log # setup log func
                           )
if __name__ == '__main__':
    sys.exit(main())
