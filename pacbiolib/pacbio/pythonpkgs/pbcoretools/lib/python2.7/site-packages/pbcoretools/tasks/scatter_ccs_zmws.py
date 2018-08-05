
"""
Scatter consensus reads by ZMW range, used for input to IsoSeq processing.
"""

import logging
import os
import sys

from pbcore.io import FastaWriter, FastaReader
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

from pbcoretools.tasks.scatter_subread_zmws import get_contract_parser_impl
import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.ccsset_zmw_scatter"
    DEFAULT_NCHUNKS = 5
    DRIVER_EXE = "python -m pbcoretools.tasks.scatter_ccs_zmws --resolved-tool-contract "
    DATASET_TYPE = FileTypes.DS_CCS
    CHUNK_KEYS = ("$chunk.ccsset_id", )
    READ_TYPE = "ConsensusRead"
    READ_TYPE_ABBREV = "ccs"


def run_main(chunk_output_json, dataset_xml, max_nchunks, output_dir):
    return CU.write_ccsset_zmw_chunks_to_file(
        chunk_file=chunk_output_json,
        dataset_path=dataset_xml,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="chunk_subreadset",
        chunk_ext=FileTypes.DS_CCS.ext)


def _args_runner(args):
    return run_main(args.chunk_report_json, args.subreadset,
                    args.max_nchunks, os.path.dirname(args.chunk_report_json))


def _rtc_runner(rtc):
    output_dir = os.path.dirname(rtc.task.output_files[0])
    max_nchunks = rtc.task.max_nchunks
    return run_main(rtc.task.output_files[0], rtc.task.input_files[0],
                    max_nchunks, output_dir)


def main(argv=sys.argv):
    mp = get_contract_parser_impl(Constants)
    return pbparser_runner(argv[1:],
                           mp,
                           _args_runner,
                           _rtc_runner,
                           log,
                           setup_log)

if __name__ == '__main__':
    sys.exit(main())
