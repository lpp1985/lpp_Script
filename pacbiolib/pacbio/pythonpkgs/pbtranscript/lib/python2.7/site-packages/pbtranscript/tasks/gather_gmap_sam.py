"""
Gather GMAP sam output.
"""

import logging
import sys

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcoretools.chunking.gather import get_datum_from_chunks_by_chunk_key

from pbtranscript.__init__ import get_version
from pbtranscript.Utils import rmpath
from pbtranscript.collapsing import concatenate_sam, sort_sam

log = logging.getLogger(__name__)


class Constants(object):
    """Constants used in pbtranscript.tasks.gather_gmap_sam"""
    TOOL_ID = "pbtranscript.tasks.gather_gmap_sam"
    VERSION = get_version()
    DRIVER = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__
    CHUNK_KEY = "$chunk.sam_id"


def get_contract_parser():
    """Return contract parser.
    input idx 0: chunk json
    output idx 0: concatenated SAM file.
    """
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Gather GMAP SAM",
                            "Chunk Gather for GMAP SAM",
                            Constants.DRIVER,
                            is_distributed=True)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "Gather CHUNK Json",
                          "Gathered CHUNK Json with chunk key")
    p.add_output_file_type(FileTypes.SAM, "sam_out", "SAM",
                           "Gathered SAM", "Gathered SAM")
    return p


def run_main(chunk_json, sam_output, chunk_key):
    """run main"""
    chunks = load_pipeline_chunks_from_json(chunk_json)

    # Allow looseness
    if not chunk_key.startswith('$chunk.'):
        chunk_key = '$chunk.' + chunk_key
        log.warn("Prepending chunk key with '$chunk.' to '%s'", str(chunk_key))

    sam_files = get_datum_from_chunks_by_chunk_key(chunks, chunk_key)
    log.debug("Chunked SAM files are %s.", (', '.join(sam_files)))

    log.info("Concatenate chunked SAM files to %s.", sam_output)

    # concatenate sam files
    unsorted_sam_output = sam_output + ".unsorted.sam"
    concatenate_sam(sam_files, unsorted_sam_output)

    # then sort
    sort_sam(unsorted_sam_output, sam_output)

    # remove intermediate file
    rmpath(unsorted_sam_output)
    return 0


def args_runner(args):
    """Args runner"""
    raise NotImplementedError()


def rtc_runner(rtc):
    """
    :type rtc: pbcommand.models.ResolvedToolContract
    :return:
    """
    return run_main(chunk_json=rtc.task.input_files[0],
                    sam_output=rtc.task.output_files[0],
                    chunk_key=Constants.CHUNK_KEY)


def main():
    """main"""
    mp = get_contract_parser()
    return pbparser_runner(sys.argv[1:],
                           mp, args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
