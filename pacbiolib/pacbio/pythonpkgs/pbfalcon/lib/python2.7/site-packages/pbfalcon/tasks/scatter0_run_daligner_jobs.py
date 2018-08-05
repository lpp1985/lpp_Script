from pbcommand.utils import setup_log
from pbcommand.cli import pbparser_runner
from pbcommand.models import get_scatter_pbparser, FileTypes
from .. import chunk as CU
import logging
import os
import sys

log = logging.getLogger(__name__)

TOOL_ID = "pbfalcon.tasks.scatter0_run_daligner_jobs"


class Constants(object):
    DEFAULT_NCHUNKS = 24
    CHUNK_KEYS = ('$chunk.json_id', '$chunk.bash_id')

def get_contract_parser():
    driver = "python -m pbfalcon.tasks.scatter0_run_daligner_jobs --resolved-tool-contract "

    p = get_scatter_pbparser(TOOL_ID, "0.1.3", "Scatter Daligner",
                             "Scatter Daligner Jobs", driver, Constants.CHUNK_KEYS,
                             is_distributed=False)

    p.add_input_file_type(FileTypes.JSON, "config", "Config",
                          "Pac Bio ???")
    p.add_input_file_type(FileTypes.TXT, "bash", "Bash",
                          "Pac Bio ???")
    p.add_output_file_type(FileTypes.CHUNK, "cjson_out", "Chunk of .fasta for daligner, stage-0",
                           "Chunked JSON Filtered Fasta",
                           "fasta.chunked")
    return p

def run_main(json_file, bash_file, output_json, max_nchunks, chunk_keys):
    log.info("Running {f} into {n} chunks".format(f=bash_file, n=max_nchunks))
    output_dir = os.path.dirname(output_json)
    #CU.write_fasta_chunks_to_file(output_json, fasta_file, max_nchunks, output_dir, "scattered-daligner-jobs", "sh")

    CU.write_run_daligner_chunks_falcon(False, output_json, json_file,
                                        bash_file,
                                        max_nchunks, output_dir,
                                       "scattered-daligner-jobs", "sh", chunk_keys)
    return 0

def _args_run_to_random_fasta_file(args):
    raise NotImplementedError('hi') #return run_main(args.fasta_in, args.cjson_out, args.max_nchunks, args.chunk_key)

def _rtc_runner(rtc):
    # the chunk key isn't really something that can be tweaked here.
    return run_main(rtc.task.input_files[0], rtc.task.input_files[1], rtc.task.output_files[0], rtc.task.max_nchunks, Constants.CHUNK_KEYS)

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
