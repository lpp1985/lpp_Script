from .. import tusks as pbfalcon
from pbcommand.pb_io import load_pipeline_chunks_from_json
from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import (FileTypes, get_gather_pbparser)
from pbcoretools.chunking.gather import (get_datum_from_chunks_by_chunk_key)
import logging
import os
import sys
cd = pbfalcon.cd

log = logging.getLogger(__name__)

TOOL_ID = "pbfalcon.tasks.gather0_run_daligner_jobs"
CHUNK_KEY = "$chunk.daligner_job_id"


def get_contract_parser():
    driver = "python -m pbfalcon.tasks.gather0_run_daligner_jobs --resolved-tool-contract "

    p = get_gather_pbparser(TOOL_ID, "0.1.3", "Gather Daligner",
                            "Gather Daligner Jobs", driver, is_distributed=False)

    p.add_input_file_type(FileTypes.CHUNK, "cjson_in", "Gather ChunkJson",
                          "Fasta Gather Chunk JSON")

    p.add_output_file_type(FileTypes.FOFN, "fofn_out", "FOFN of .las files, stage-0",
                           "FOFN Not Sure Yet",
                           "file_gathered")

    # FIXME. There's an a bit of friction between the TC and the argrunning
    # layer here. For nproc and nchunks, chunk-key, the values are a bit
    # different between the two backends.
    p.arg_parser.add_str("pbsmrtpipe.task_options.task_falcon_run_daligner_job_scatter_chunk_key", "chunk_key",
                         "$chunk.daligner_job_id", "Chunk key", "Chunk key to use (format $chunk:{chunk-key}")
    return p


def run_main(chunk_json, fofn_output, chunk_key):
  with cd(os.path.dirname(fofn_output)):
    chunks = load_pipeline_chunks_from_json(chunk_json)

    # Allow looseness
    if not chunk_key.startswith('$chunk.'):
        chunk_key = '$chunk.' + chunk_key
        log.warn("Prepending chunk key with '$chunk.' to '{c}'".format(c=chunk_key))

    fofn_files = get_datum_from_chunks_by_chunk_key(chunks, chunk_key)
    print("fofn_files:%s %s" %(repr(fofn_files), repr(fofn_output)))
    # Combine all into one.
    with open(fofn_output, 'w') as ofs:
        for fn in fofn_files:
            with open(fn) as ifs:
                ofs.write(ifs.read())


def args_runner(args):
    return run_main(args.cjson_in, args.fofn_out, args.chunk_key)


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
