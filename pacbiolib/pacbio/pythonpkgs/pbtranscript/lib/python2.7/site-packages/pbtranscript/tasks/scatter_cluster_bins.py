"""
Specialized scatter task which spawns cluster_chunks.pickle.
"""

import logging
import sys
import os.path as op

from pbcommand.models import get_scatter_pbparser, FileTypes, PipelineChunk
from pbcommand.pb_io import write_pipeline_chunks
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbtranscript.tasks.TPickles import ChunkTasksPickle, ClusterChunkTask
log = logging.getLogger(__name__)


class Constants(object):
    """Constants used for scatter_cluster_bins."""
    TOOL_ID = "pbtranscript.tasks.scatter_cluster_bins"
    DEFAULT_NCHUNKS = 24
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID
    CHUNK_KEYS = ('$chunk.cluster_chunk_pickle_id', '$chunk.ccs_id')


def get_contract_parser():
    """Get scatter tool contract parser
    Input:
        idx 0 - cluster_chunks.pickle
        idx 1 - ccs
    Output:
        idx 0 - chunk.json
    """
    p = get_scatter_pbparser(Constants.TOOL_ID, Constants.VERSION,
                             "Scatter Cluster Bins",
                             __doc__, Constants.DRIVER_EXE,
                             chunk_keys=Constants.CHUNK_KEYS,
                             is_distributed=True)

    p.add_input_file_type(FileTypes.PICKLE, "cluster_chunks_pickle", "Pickle In",
                          "Cluster chunks pickle file") # input 0
    p.add_input_file_type(FileTypes.DS_CCS, "ccs_in", "ConsensusReadSet In",
                          "PacBio ConsensusReadSet") # input 1
    p.add_output_file_type(FileTypes.CHUNK, "cjson_out",
                           "Chunk JSON Cluster Bins",
                           "Chunked JSON Cluster Bins",
                           "ice_cluster.chunked") # output 0
    # max nchunks for this specific task
    p.add_int("pbsmrtpipe.task_options.dev_scatter_max_nchunks", "max_nchunks",
              Constants.DEFAULT_NCHUNKS, "Max NChunks", "Maximum number of Chunks")
    return p


def run_main(cluster_chunks_pickle_file, ccs_file, output_json_file, max_nchunks):
    """Scatter items in cluster_chunks_pickle
    Parameters:
      cluster_chunks_pickle_file -- ChunkTasksPickle of ClusterChunkTask objects.
      ccs_file -- ccs.consensusreadset.xml
      output_json_file -- chunk.json
      max_nchunks -- maximum # of chunks
    """
    p = ChunkTasksPickle.read(cluster_chunks_pickle_file)
    assert all([isinstance(r, ClusterChunkTask) for r in p])
    out_dir = op.dirname(output_json_file)

    # sort and group tasks
    groups = p.sort_and_group_tasks(max_nchunks=max_nchunks)

    # Writing chunk.json
    base_name = "spawned_cluster_chunk"
    chunks = []
    spawned_pickles = []
    for group_index in range(0, len(groups)):
        chunk_id = "_".join([base_name, 'group', str(group_index)])
        spawned_pickle_file = op.join(out_dir, chunk_id + ".pickle")
        d = {Constants.CHUNK_KEYS[0]: spawned_pickle_file,
             Constants.CHUNK_KEYS[1]: ccs_file}
        c = PipelineChunk(chunk_id, **d)
        chunks.append(c)
        spawned_pickles.append(spawned_pickle_file)

    log.info("Spawning %s into %d files", cluster_chunks_pickle_file, len(groups))
    p.spawn_pickles_by_groups(groups, spawned_pickles)
    log.debug("Spawned files: %s.", ", ".join(spawned_pickles))

#    n_chunks = len(p)
#    for i in range(0, n_chunks):
#        chunk_id = "_".join([base_name, str(i)])
#        spawned_pickle_file = op.join(out_dir, chunk_id + ".pickle")
#        d = {Constants.CHUNK_KEYS[0]: spawned_pickle_file,
#             Constants.CHUNK_KEYS[1]: ccs_file}
#        c = PipelineChunk(chunk_id, **d)
#        chunks.append(c)
#        spawned_pickles.append(spawned_pickle_file)
#
#    log.info("Spawning %s into %s files", cluster_chunks_pickle_file, str(n_chunks))
#    p.spawn_pickles(spawned_pickles)
#    log.debug("Spawned files: %s.", ", ".join(spawned_pickles))

    log.info("Writing chunk.json to %s", output_json_file)
    write_pipeline_chunks(chunks, output_json_file,
                          "created by %s" % Constants.TOOL_ID)
    return 0


def _args_run(args):
    raise NotImplementedError()


def _rtc_runner(rtc):
    return run_main(cluster_chunks_pickle_file=rtc.task.input_files[0],
                    ccs_file=rtc.task.input_files[1],
                    output_json_file=rtc.task.output_files[0],
                    max_nchunks=rtc.task.max_nchunks)


def main():
    """Main"""
    return pbparser_runner(sys.argv[1:],
                           get_contract_parser(),
                           _args_run,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
