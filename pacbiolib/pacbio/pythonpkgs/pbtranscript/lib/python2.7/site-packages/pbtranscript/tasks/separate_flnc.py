# wraps 'pbtranscript.tasks.separate_flnc'

"""
Separate full-length-non-chimeric ccs reads into multiple
bins based on either primer or read length.
Inputs: isoseq_flnc.fasta|contigset.xml
Outputs:
        <out_pickle>
        <root_dir>/<bin>/isoseq_flnc.fasta|contigset.xml
"""

import logging
import os.path as op
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.utils import setup_log
from pbcommand.models import FileTypes

from pbtranscript.PBTranscriptOptions import get_base_contract_parser
from pbtranscript.separate_flnc import SeparateFLNCRunner


log = logging.getLogger(__name__)


class Constants(object):
    """Constants used in tool contract."""
    TOOL_ID = "pbtranscript.tasks.separate_flnc"
    DRIVER_EXE = "python -m pbtranscript.tasks.separate_flnc --resolved-tool-contract"
    PARSER_DESC = __doc__

    BIN_BY_PRIMER_ID = "pbtranscript.task_options.bin_by_primer"
    BIN_BY_PRIMER_DEFAULT = False
    BIN_BY_PRIMER_DESC = "Instead of binning by size, bin by primer " + \
                         "(overwrites --bin_size_kb and --bin_manual) " + \
                         "(default: %s)" % BIN_BY_PRIMER_DEFAULT

    BIN_SIZE_KB_ID = "pbtranscript.task_options.bin_size_kb"
    BIN_SIZE_KB_DEFAULT = 1
    BIN_SIZE_KB_DESC = "Bin size by kb (default: %s)" % BIN_SIZE_KB_DEFAULT

    BIN_MANUAL_ID = "pbtranscript.task_options.bin_manual"
    BIN_MANUAL_DEFAULT = "[]"
    BIN_MANUAL_DESC = "Bin manual (ex: (1,2,3,5)), overwrites bin_size_kb " + \
                      "(default: %s)" % BIN_MANUAL_DEFAULT

    MAX_BASE_LIMIT_MB_DEFAULT = 600
    MAX_BASE_LIMIT_MB_DESC = "Maximum number of bases per partitioned bin, in MB " + \
                             "(default: %s)" % MAX_BASE_LIMIT_MB_DEFAULT


def add_separate_flnc_io_arguments(arg_parser):
    """Add separate flnc io arguments."""
    helpstr = "Input full-length non-chimeric reads in FASTA or ContigSet format, " + \
              "used for clustering consensus isoforms, e.g., isoseq_flnc.fasta"
    arg_parser.add_argument("flnc_fa", type=str, help=helpstr)

    helpstr = "A directory to store separated reads."
    arg_parser.add_argument("root_dir", type=str, help=helpstr)

    helpstr = "Python pickle file of how flnc reads are separated"
    arg_parser.add_argument("out_pickle", type=str, help=helpstr)

    return arg_parser


def add_separate_flnc_arguments(arg_parser):
    """Add arg parser arguments."""
    sepa_group = arg_parser.add_argument_group("Separate FLNC arguments")
    sepa_group.add_argument("--bin_size_kb", default=Constants.BIN_SIZE_KB_DEFAULT,
                            type=int, help=Constants.BIN_SIZE_KB_DESC)

    sepa_group.add_argument("--bin_manual", default=Constants.BIN_MANUAL_DEFAULT,
                            help=Constants.BIN_MANUAL_DESC)

    sepa_group.add_argument("--bin_by_primer", default=Constants.BIN_BY_PRIMER_DEFAULT,
                            action="store_true", help=Constants.BIN_BY_PRIMER_DESC)

    sepa_group.add_argument("--max_base_limit_MB", default=Constants.MAX_BASE_LIMIT_MB_DEFAULT,
                            type=int, help=Constants.MAX_BASE_LIMIT_MB_DESC)
    return arg_parser


def add_separate_flnc_tcp_options(tcp):
    """Add tcp options."""
    tcp.add_int(Constants.BIN_SIZE_KB_ID, "bin_size_kb",
                default=Constants.BIN_SIZE_KB_DEFAULT,
                name="Bin by read length in KB", description=Constants.BIN_SIZE_KB_DESC)
    tcp.add_str(Constants.BIN_MANUAL_ID, "bin_manual",
                default=Constants.BIN_MANUAL_DEFAULT,
                name="Bin by read length manually", description=Constants.BIN_MANUAL_DESC)
    tcp.add_boolean(Constants.BIN_BY_PRIMER_ID, "bin_by_primer",
                    default=Constants.BIN_BY_PRIMER_DEFAULT,
                    name="Bin by primer", description=Constants.BIN_BY_PRIMER_DESC)
    return tcp


def args_runner(args):
    """Run given input args"""
    s = SeparateFLNCRunner(flnc_fa=args.flnc_fa, root_dir=args.root_dir, out_pickle=args.out_pickle,
                           bin_size_kb=args.bin_size_kb, bin_by_primer=args.bin_by_primer,
                           bin_manual=args.bin_manual, max_base_limit_MB=args.max_base_limit_MB)
    s.run()


def get_contract_parser():
    """Return a tool contract parser for separate_flnc.
    Input:
        idx 0 - flnc.fa|fq|ds
    Output:
        idx 0 - out.pickle
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    add_separate_flnc_io_arguments(p.arg_parser.parser)
    add_separate_flnc_arguments(p.arg_parser.parser)

    # tool contract parser
    tcp = p.tool_contract_parser
    tcp.add_input_file_type(FileTypes.DS_CONTIG, "flnc_fa",
                            name="FASTA or ContigSet file",
                            description="FLNC reads ContigSet")
    tcp.add_output_file_type(FileTypes.PICKLE, "out_pickle",
                             default_name="separate_flnc",
                             name="Bins of FLNC Reads", description="output bins in pickle.")
    add_separate_flnc_tcp_options(tcp)
    return p


def resolved_tool_contract_runner(rtc):
    """Given an input resolved tool contract, convert it to args and
       call args_runner."""
    # rtc has only one input, idx 0, isoseq_flnc.contigset.xml
    # rtc has only one output, idx 0, separate_flnc.pickle
    flnc_fa = rtc.task.input_files[0]
    out_pickle = rtc.task.output_files[0]
    root_dir = op.dirname(out_pickle)

    bin_by_primer = rtc.task.options[Constants.BIN_BY_PRIMER_ID] is True

    bin_manual = rtc.task.options[Constants.BIN_MANUAL_ID]
    if (isinstance(bin_manual, str) or isinstance(bin_manual, unicode)) and \
       len(str(bin_manual).translate(None, '[]() ')) > 0:
        log.info("Setting bin manually, bin_manual is %s", bin_manual)
        bin_manual = [int(x) for x in str(bin_manual).translate(None, '[]() ').split(',')]
    else: # bin_manual is None
        bin_manual = None

    bin_size_kb = rtc.task.options[Constants.BIN_SIZE_KB_ID]

    s = SeparateFLNCRunner(flnc_fa=flnc_fa, root_dir=root_dir, out_pickle=out_pickle,
                           bin_size_kb=bin_size_kb, bin_by_primer=bin_by_primer,
                           bin_manual=bin_manual,
                           max_base_limit_MB=Constants.MAX_BASE_LIMIT_MB_DEFAULT)
    s.run()
    return 0


def main(argv=sys.argv[1:]):
    """Main"""
    return pbparser_runner(
        argv=argv,
        parser=get_contract_parser(),
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
