import logging
import os
import sys

from pbcore.io import FastaWriter, FastaReader
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.subreadset_align_scatter"


class Constants(object):
    DEFAULT_NCHUNKS = 5


def get_contract_parser():
    driver = "python -m pbcoretools.tasks.scatter_subread_reference --resolved-tool-contract "

    # These Keys are expected to be PipelineChunks produced by this tool
    chunk_keys = ("$chunk.reference_id", "$chunk.subreadset_id")
    p = get_scatter_pbparser(TOOL_ID, "0.1.3", "SubreadSet scatter",
                             "Scatter Subread DataSet", driver, chunk_keys,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.DS_SUBREADS,
                          "subreads",
                          "SubreadSet",
                          "Pac Bio Fasta format")

    p.add_input_file_type(FileTypes.DS_REF,
                          "ds_reference",
                          "ReferenceSet",
                          "Pac Bio Fasta format")

    p.add_output_file_type(FileTypes.CHUNK,
                           "cjson_out",
                           "Chunk SubreadSet",
                           "PacBio Chunked JSON SubreadSet",
                           "subreadset_chunked")

    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.scatter_subread_max_nchunks", "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")

    # This should only be added at the argparse level.
    # Disabling for now.
    # FIXME. This should support --reference-chunk-key and --subread-key
    # p.arg_parser.add_str("pbcoretools.task_options.scatter_subreadset_chunk_key",
    #                      "chunk_key",
    #                      "$chunk:fasta_id", "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def run_main(chunk_output_json, subread_xml, reference_set_xml, max_nchunks, output_dir):
    return CU.write_subreadset_chunks_to_file(chunk_output_json,
                                              subread_xml,
                                              reference_set_xml,
                                              max_nchunks,
                                              output_dir,
                                              "chunk_subreadset",
                                              FileTypes.DS_SUBREADS.ext)


def _args_run_to_random_fasta_file(args):
    output_dir = os.path.dirname(args.cjson_out)
    return run_main(args.cjson_out, args.subreads,
                    args.ds_reference, args.max_nchunks, output_dir)


def _rtc_runner(rtc):
    output_dir = os.path.dirname(rtc.task.output_files[0])
    max_nchunks = rtc.task.max_nchunks
    return run_main(rtc.task.output_files[0], rtc.task.input_files[0], rtc.task.input_files[1], max_nchunks, output_dir)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           _args_run_to_random_fasta_file,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
