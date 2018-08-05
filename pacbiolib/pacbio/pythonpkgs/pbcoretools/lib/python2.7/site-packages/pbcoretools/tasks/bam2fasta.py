
"""
bam2fasta export with tmp dir support
"""

import logging
import sys

from pbcommand.models import FileTypes, ResourceTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.tasks.converters import run_bam_to_fasta

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.bam2fasta"
    VERSION = "0.2.0"
    DRIVER = "python -m pbcoretools.tasks.bam2fasta --resolved-tool-contract"

def get_parser():
    p = get_pbparser(Constants.TOOL_ID,
                     Constants.VERSION,
                     "bam2fasta export",
                     __doc__,
                     Constants.DRIVER,
                     is_distributed=True,
                     resource_types=(ResourceTypes.TMP_DIR,))
    p.add_input_file_type(FileTypes.DS_SUBREADS, "subreads",
                          "Input Subreads",
                          "Input SubreadSet XML")
    p.add_output_file_type(FileTypes.FASTA,
                           "fasta_out",
                           "FASTA subreads",
                           description="Exported FASTA",
                           default_name="subreads")
    return p


def run_args(args):
    return run_bam_to_fasta(args.subreads, args.fasta_out)


def run_rtc(rtc):
    return run_bam_to_fasta(rtc.task.input_files[0], rtc.task.output_files[0],
                            tmp_dir=rtc.task.tmpdir_resources[0].path)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           run_args,
                           run_rtc,
                           log,
                           setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
