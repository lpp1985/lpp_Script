import logging
import sys

from pbcommand.pb_io import load_pipeline_chunks_from_json

from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import (FileTypes, get_gather_pbparser)

from pbcoretools.chunking.gather import (gather_fasta,
                                     get_datum_from_chunks_by_chunk_key)

log = logging.getLogger(__name__)

TOOL_ID = "pbcoretools.tasks.gather_fasta"
CHUNK_KEY = "$chunk.fasta_id"


def get_contract_parser():
    driver = "python -m pbcoretools.tasks.gather_fasta --resolved-tool-contract "

    p = get_gather_pbparser(TOOL_ID, "0.1.3", "Gather Fasta",
                            "Gather Fasta", driver, is_distributed=True)

    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "Gather ChunkJson",
                          "Fasta Gather Chunk JSON")

    p.add_output_file_type(FileTypes.FASTA, "fasta_out", "Fasta Gathered",
                           "Fasta Gathered",
                           "gathered")

    # FIXME. There's an a bit of friction between the TC and the argrunning
    # layer here. For nproc and nchunks, chunk-key, the values are a bit
    # different between the two backends.
    p.arg_parser.add_str("pbcoretools.task_options.dev_scatter_chunk_key", "chunk_key",
                         "$chunk.fasta_id", "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def run_main(chunk_json, fasta_output, chunk_key):
    chunks = load_pipeline_chunks_from_json(chunk_json)

    # Allow looseness
    if not chunk_key.startswith('$chunk.'):
        chunk_key = '$chunk.' + chunk_key
        log.warn("Prepending chunk key with '$chunk.' to '{c}'".format(c=chunk_key))

    fastx_files = get_datum_from_chunks_by_chunk_key(chunks, chunk_key)
    _ = gather_fasta(fastx_files, fasta_output)

    return 0


def args_runner(args):
    return run_main(args.cjson_in, args.fasta_out, args.chunk_key)


def rtc_runner(rtc):
    """
    :type rtc: pbcommand.models.ResolvedToolContract
    :return:
    """
    # the input file is just a sentinel file
    return run_main(rtc.task.input_files[0], rtc.task.output_files[0], rtc.task.chunk_key)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp, args_runner,
                           rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
