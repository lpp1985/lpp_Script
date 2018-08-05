import logging
import os
import sys

from pbcore.io import FastaWriter, FastaReader
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.ccsset_align_scatter"


class Constants(object):
    DEFAULT_NCHUNKS = 5
    CHUNK_KEYS = ('$chunk.ccsset_id', "$chunk.reference_id")

def get_contract_parser():
    driver = "python -m pbcoretools.tasks.scatter_ccs_reference --resolved-tool-contract "

    p = get_scatter_pbparser(TOOL_ID, "0.1.3", "ConsensusReadSet scatter",
                             "Scatter ConsensusRead DataSet", driver,
                             Constants.CHUNK_KEYS, is_distributed=True)

    p.add_input_file_type(FileTypes.DS_CCS,
                          "ccsset",
                          "ConsensusReadSet",
                          "Pac Bio Fasta format")

    p.add_input_file_type(FileTypes.DS_REF,
                          "ds_reference",
                          "ReferenceSet",
                          "Pac Bio Fasta format")

    p.add_output_file_type(FileTypes.CHUNK,
                           "cjson_out",
                           "Chunk ConsensusReadSet",
                           "PacBio Chunked JSON ConsensusReadSet",
                           "ccsset_chunked")

    # max nchunks for this specific task
    # FIXME using same option names as scatter_subread_reference.py - it would
    # be nice if these were more generic
    p.add_int("pbcoretools.task_options.scatter_subread_max_nchunks", "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")

    p.add_str("pbcoretools.task_options.scatter_subreadset_chunk_key", "chunk_key",
              "$chunk:fasta_id", "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def run_main(chunk_output_json, ccs_xml, reference_set_xml, max_nchunks, output_dir):
    return CU.write_ccsset_chunks_to_file(chunk_output_json,
                                          ccs_xml,
                                          reference_set_xml,
                                          max_nchunks,
                                          output_dir,
                                          "chunk_ccsset", FileTypes.DS_CCS.ext)


def _args_runner(args):
    return run_main(args.chunk_report_json, args.ccsset,
                    args.ds_reference, args.max_total_chunks, args.output_dir)


def _rtc_runner(rtc):
    output_dir = os.path.dirname(rtc.task.output_files[0])
    max_nchunks = rtc.task.max_nchunks
    return run_main(rtc.task.output_files[0], rtc.task.input_files[0], rtc.task.input_files[1], max_nchunks, output_dir)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           _args_runner,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
