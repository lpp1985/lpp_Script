
"""
Gather task for the pickle files output by IcePartial.
"""

import logging
import sys

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.models import get_gather_pbparser, FileTypes
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbcoretools.chunking.gather import get_datum_from_chunks_by_chunk_key

from pbtranscript.ice.IceUtils import combine_nfl_pickles

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbtranscript.tasks.gather_nfl_pickle"
    VERSION = "0.1.0"
    DRIVER = "python -m pbtranscript.tasks.gather_nfl_pickle  --resolved-tool-contract"
    CHUNK_KEY = "$chunk.pickle_id"
    OPT_CHUNK_KEY = "pbtranscript.task_options.pickle_gather_chunk_key"


def get_parser():
    p = get_gather_pbparser(Constants.TOOL_ID,
                            Constants.VERSION,
                            "Dev Pickle Gather",
                            "General Chunk Pickle Gather",
                            Constants.DRIVER,
                            is_distributed=True)
    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "GCHUNK Json",
                          "Gathered CHUNK Json with Pickle chunk key")
    p.add_output_file_type(FileTypes.PICKLE, "pickle_out", "Pickle",
                           "Gathered Pickle", "gathered")
    p.arg_parser.add_str(Constants.OPT_CHUNK_KEY,
                         "chunk_key",
                         default=Constants.CHUNK_KEY,
                         name="Chunk key",
                         description="Chunk key to use (format $chunk:{chunk-key}")
    return p


def run(chunk_input_json, output_file, chunk_key):
    chunks = load_pipeline_chunks_from_json(chunk_input_json)
    chunked_files = get_datum_from_chunks_by_chunk_key(chunks, chunk_key)
    _ = combine_nfl_pickles(chunked_files, output_file)
    return 0


def args_runner(args):
    raise NotImplementedError()


def rtc_runner(rtc):
    return run(
        chunk_input_json=rtc.task.input_files[0],
        output_file=rtc.task.output_files[0],
        chunk_key=Constants.CHUNK_KEY)


def main(argv=sys.argv):
    return pbparser_runner(argv[1:],
                           get_parser(),
                           args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
