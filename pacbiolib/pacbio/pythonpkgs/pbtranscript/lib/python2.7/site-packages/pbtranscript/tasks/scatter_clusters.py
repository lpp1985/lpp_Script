
"""
Specialized scatter task for pickle file containing clusters
"""

import logging
import cPickle
import os
import sys

from pbcommand.models import get_scatter_pbparser, FileTypes, PipelineChunk
from pbcommand.pb_io import write_pipeline_chunks
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

log = logging.getLogger(__name__)


class Constants(object):
    TOOL_ID = "pbtranscript.tasks.scatter_clusters"
    DEFAULT_NCHUNKS = 24
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m pbtranscript.tasks.scatter_clusters --resolved-tool-contract "
    CHUNK_KEYS = ('$chunk.subreadset_id', '$chunk.contigset_id',
                  '$chunk.pickle_id', '$chunk.nfl_pickle_id')


def get_contract_parser():

    p = get_scatter_pbparser(Constants.TOOL_ID, Constants.VERSION,
                             "Scatter ContigSet",
                             __doc__, Constants.DRIVER_EXE,
                             chunk_keys=Constants.CHUNK_KEYS,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.DS_SUBREADS, "subreads_in", "SubreadSet In",
                          "PacBio ContigSet")
    p.add_input_file_type(FileTypes.DS_CONTIG, "fasta_in", "ContigSet",
                          "PacBio ContigSet")
    p.add_input_file_type(FileTypes.PICKLE, "pickle_in", "Pickle",
                          "Cluster pickle file")
    p.add_input_file_type(FileTypes.PICKLE, "nfl_pickle_in", "Pickle",
                          "Non-full-length pickle file")
    p.add_output_file_type(FileTypes.CHUNK, "cjson_out",
                           "Chunk JSON Filtered Fasta",
                           "Chunked JSON ContigSet",
                           "pickles.chunked")
    # max nchunks for this specific task
    p.add_int("pbsmrtpipe.task_options.dev_scatter_max_nchunks", "max_nchunks",
              default=Constants.DEFAULT_NCHUNKS,
              name="Max NChunks",
              description="Maximum number of Chunks")
    return p


def run_main(subreads_file, isoforms_file, cluster_pickle_file,
             nfl_pickle_file, output_json, max_nchunks):
    log.info("Running {f} into {n} chunks".format(f=cluster_pickle_file,
                                                  n=max_nchunks))
    uc = {}
    with open(cluster_pickle_file, 'rb') as f:
        a = cPickle.load(f)
        uc = a['uc']
    assert len(uc) > 0
    n_chunks = min(len(uc), max_nchunks)
    base_name = "cluster_chunk"
    dir_name = os.path.dirname(output_json)
    chunks = []
    for i in range(n_chunks):
        chunk_id = "_".join([base_name, str(i)])
        chunk_name = ".".join([chunk_id, "pickle"])
        chunk_pickle_file = os.path.join(dir_name, chunk_name)
        with open(chunk_pickle_file, 'wb') as f:
            cPickle.dump({
                '__chunk_i': i,
                '__chunk_n': n_chunks,
                'pickle_file': cluster_pickle_file,
            }, f)
        d = {
            '$chunk.subreadset_id': subreads_file,
            '$chunk.contigset_id': isoforms_file,
            '$chunk.nfl_pickle_id': nfl_pickle_file,
            '$chunk.pickle_id': chunk_pickle_file,
        }
        c = PipelineChunk(chunk_id, **d)
        chunks.append(c)
    write_pipeline_chunks(chunks, output_json,
        "created by pbtranscript.tasks.scatter_clusters")
    return 0


def _args_run(args):
    raise NotImplementedError()


def _rtc_runner(rtc):
    # the chunk key isn't really something that can be tweaked here.
    return run_main(rtc.task.input_files[0], rtc.task.input_files[1], rtc.task.input_files[2], rtc.task.input_files[3], rtc.task.output_files[0], rtc.task.max_nchunks)


def main(argv=sys.argv):
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           _args_run,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
