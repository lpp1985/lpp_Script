
"""
Scatter subreads by ZMW, used for bam2bam barcoding.
"""

import functools
import logging
import os
import sys

from pbcore.io import FastaWriter, FastaReader
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes

import pbcoretools.chunking.chunk_utils as CU

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbcoretools.tasks.scatter_subreads_bam2bam"
    DEFAULT_NCHUNKS = 5
    DRIVER_EXE = "python -m pbcoretools.tasks.scatter_subreads_bam2bam --resolved-tool-contract "
    DATASET_TYPE = FileTypes.DS_SUBREADS
    CHUNK_KEYS = ("$chunk.subreadset_id", )
    READ_TYPE = "Subread"
    READ_TYPE_ABBREV = "subread"

def get_contract_parser_impl(C):
    p = get_scatter_pbparser(C.TOOL_ID, "0.1.3",
        "%sSet ZMW scatter" % C.READ_TYPE,
        "Scatter %s DataSet for barcoding" % C.READ_TYPE, C.DRIVER_EXE,
        C.CHUNK_KEYS, is_distributed=True)

    p.add_input_file_type(C.DATASET_TYPE,
                          "dataset",
                          "%sSet" % C.READ_TYPE,
                          "Pac Bio Fasta format")
    p.add_input_file_type(FileTypes.DS_BARCODE,
                          "barcodes",
                          "BarcodeSet",
                          "Pac Bio Barcode Dataset XML")
    p.add_output_file_type(FileTypes.CHUNK,
                           "chunk_report_json",
                           "Chunk %sSet" % C.READ_TYPE,
                           "PacBio Chunked JSON %sSet" % C.READ_TYPE,
                           "%sset_chunked" % C.READ_TYPE_ABBREV)

    # max nchunks for this specific task
    p.add_int("pbcoretools.task_options.scatter_subread_max_nchunks",
              "max_nchunks", Constants.DEFAULT_NCHUNKS,
              "Max NChunks", "Maximum number of Chunks")

    p.add_str("pbcoretools.task_options.scatter_subreadset_chunk_key",
              "chunk_key", "$chunk:subreadset_id", "Chunk key",
              "Chunk key to use (format $chunk:{chunk-key}")
    return p

get_contract_parser = functools.partial(get_contract_parser_impl, Constants)

def run_main(chunk_output_json, dataset_xml, barcode_xml, max_nchunks,
             output_dir):
    return CU.write_subreadset_zmw_chunks_to_file(
        chunk_file=chunk_output_json,
        dataset_path=dataset_xml,
        max_total_chunks=max_nchunks,
        dir_name=output_dir,
        chunk_base_name="chunk_dataset",
        chunk_ext=FileTypes.DS_SUBREADS.ext,
        extra_chunk_keys={"$chunk.barcodeset_id":barcode_xml})


def _args_runner(args):
    return run_main(args.chunk_report_json, args.subreadset, args.barcodeset,
                    args.max_nchunks, os.path.dirname(args.chunk_report_json))


def _rtc_runner(rtc):
    output_dir = os.path.dirname(rtc.task.output_files[0])
    max_nchunks = rtc.task.max_nchunks
    return run_main(rtc.task.output_files[0], rtc.task.input_files[0],
                    rtc.task.input_files[1],
                    max_nchunks, output_dir)


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
