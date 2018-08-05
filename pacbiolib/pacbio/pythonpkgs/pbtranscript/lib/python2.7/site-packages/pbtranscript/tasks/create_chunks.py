"""
Create chunk tasks for ICE, ice_partial, and ice_polish, after
binning (separate_flnc) is done.
  nfl.contigset will be scattered into chunks when chunk mode is true.
  All chunk tasks for ICE will be saved to ice_cluster.pickle as a list.
  All chunk tasks for ice_partial will be saved to ice_partial.pickle as a list.
  All chunk tasks for ice_polish will be saved to ice_polish.pickle as a list.
"""

import logging
import sys
import os.path as op

from pbcommand.models import FileTypes, get_pbparser, SymbolTypes
from pbcommand.cli import pbparser_runner

from pbcommand.pb_io.common import load_pipeline_chunks_from_json
from pbcommand.utils import setup_log

import pbcoretools.chunking.chunk_utils as CU
from pbcoretools.chunking.gather import get_datum_from_chunks_by_chunk_key

from pbtranscript.Utils import ln
from pbtranscript.tasks.TPickles import ClusterChunkTask, PartialChunkTask,\
        PolishChunkTask, ChunkTasksPickle, n_reads_in_contigsets
from pbtranscript.separate_flnc import SeparateFLNCBase


log = logging.getLogger(__name__)


class Constants(object):
    """Constants used for create_chunks."""
    TOOL_ID = "pbtranscript.tasks.create_chunks"
    DEFAULT_NCHUNKS = 24
    VERSION = "0.1.0"
    DRIVER_EXE = "python -m %s --resolved-tool-contract " % TOOL_ID


def get_contract_parser():
    """Get tool contract parser.
    Input:
        idx 0 - separate_flnc.pickle
        idx 1 - nfl.contigset

    Output:
        idx 0 - ice_cluster_chunks.pickle
        idx 1 - ice_partial_chunks.pickle
        idx 2 - ice_polish_chunks.pickle
    """
    p = get_pbparser(tool_id=Constants.TOOL_ID,
                     version=Constants.VERSION,
                     name=Constants.TOOL_ID,
                     description=__doc__,
                     driver_exe=Constants.DRIVER_EXE,
                     nproc=SymbolTypes.MAX_NPROC)

    p.add_input_file_type(FileTypes.PICKLE, "separate_flnc_pickle_in",
                          "Pickle In",
                          "Separate flnc pickle file") # input 0

    p.add_input_file_type(FileTypes.PICKLE, "nfl_contigset",
                          "ContigSet In",
                          "Nfl Contigset") # input 1

    p.add_output_file_type(FileTypes.PICKLE, "cluster_chunks",
                           "ICE CLUSTER CHUNK PICKLE",
                           "Pickle containing ICE chunk tasks",
                           "cluster_chunks") # output 0

    p.add_output_file_type(FileTypes.PICKLE, "partial_chunks",
                           "ICE PARTIAL CHUNK PICKLE",
                           "Pickle containing ice_partial chunk tasks",
                           "partial_chunks") # output 1

    p.add_output_file_type(FileTypes.PICKLE, "polish_chunks",
                           "ICE POLISH CHUNK PICKLE",
                           "Pickle containing ice_polish (quiver|arrow) chunk tasks",
                           "polish_chunks") # output 2
    return p


def _get_cluster_out_dir(flnc_file):
    """Return cluster out dir given flnc file."""
    return op.join(op.dirname(flnc_file), "cluster_out")


def create_cluster_pickle(flnc_files, out_pickle):
    """Create cluster chunk task pickle.
    Parameters:
      n_bins -- number of bins
      flnc_files -- full-length non-chimeric files in bins
      out_pickle -- output pickle for saving ClusterChunkTask objects
    """
    n_bins = len(flnc_files)
    log.info("Writing %s cluster chunk tasks to %s.", str(n_bins), out_pickle)
    p = ChunkTasksPickle()

    for i, flnc_file in enumerate(flnc_files):
        log.debug("Processing cluster bin index=%s.", i)
        cluster_out_dir = _get_cluster_out_dir(flnc_file)

        # Create Cluster chunk tasks.
        task_ = ClusterChunkTask(cluster_bin_index=i, flnc_file=flnc_file,
                                 cluster_out_dir=cluster_out_dir)
        p.append(task_)

    p.write(out_pickle)
    log.info("Saved %s cluster chunk tasks to %s.", str(n_bins), out_pickle)


def create_partial_pickle(flnc_files, chunked_nfl_files, out_pickle):
    """
    Parameters:
      flnc_files -- full-length non-chimeric files in bins
      chunked_nfl_files -- chunked non-chimeric files
      out_pickle -- output pickle for saving PolishChunkTask objects
    """
    n_bins = len(flnc_files)
    n_nfl_chunks = max(1, len(chunked_nfl_files))

    log.info("Writing %s ice_partial chunk tasks to %s.", str(n_bins * n_nfl_chunks), out_pickle)
    p = ChunkTasksPickle()

    for i, flnc_file in enumerate(flnc_files):
        log.debug("Processing cluster bin index=%s.", i)
        cluster_out_dir = _get_cluster_out_dir(flnc_file)

        for j, nfl_file in enumerate(chunked_nfl_files):
            # Create Partial chunk tasks.
            task_ = PartialChunkTask(cluster_bin_index=i, flnc_file=flnc_file,
                                     cluster_out_dir=cluster_out_dir,
                                     nfl_file=nfl_file,
                                     nfl_index=j, n_nfl_chunks=n_nfl_chunks)
            p.append(task_)

    p.write(out_pickle)
    log.info("Saved %s partial chunk tasks to %s.", str(n_bins * n_nfl_chunks), out_pickle)


def create_polish_pickle(n_polish_chunks_in_bins, flnc_files, out_pickle):
    """
    Parameters:
      n_polish_chunks_in_bins -- number of ice_polish chunks in each bin
      flnc_files -- full-length non-chimeric files in bins
      out_pickle -- output pickle for saving PolishChunkTask objects
    """
    n_bins = len(flnc_files)
    assert isinstance(n_polish_chunks_in_bins, list)
    assert len(n_polish_chunks_in_bins) == n_bins

    log.info("Writing %s ice_polish chunk tasks to %s.",
             str(sum(n_polish_chunks_in_bins)), out_pickle)
    p = ChunkTasksPickle()

    for i, flnc_file in enumerate(flnc_files):
        log.debug("Creating %s ice_polish chunks for bin index=%s.",
                  str(n_polish_chunks_in_bins[i]), str(i))
        cluster_out_dir = _get_cluster_out_dir(flnc_file)

        for j in range(0, n_polish_chunks_in_bins[i]):
            # Create Polish chunk tasks.
            task_ = PolishChunkTask(cluster_bin_index=i, flnc_file=flnc_file,
                                    cluster_out_dir=cluster_out_dir,
                                    polish_index=j,
                                    n_polish_chunks=n_polish_chunks_in_bins[i])
            p.append(task_)

    p.write(out_pickle)
    log.info("Saved %s polish chunk tasks to %s.", str(sum(n_polish_chunks_in_bins)), out_pickle)


def chunk_contigset(in_file, n_chunks, out_dir, out_chunk_json):
    """
    Chunk input contigset into n_chunks under out_dir, and
    write chunk info to out_chunk_json, return chunked files.
    """
    log.info("Splitting %s into at most %s chunks", in_file, str(n_chunks))
    CU.write_contigset_chunks_to_file(out_chunk_json, in_file, n_chunks,
                                      out_dir, "scattered-nfl", "contigset.xml")

    out_chunks = load_pipeline_chunks_from_json(out_chunk_json)
    chunked_files = get_datum_from_chunks_by_chunk_key(out_chunks, '$chunk.contigset_id')
    log.info("Splitted files are %s\n", ("\n".join(chunked_files)))

    # Return chunked files from out_chunk_json
    return chunked_files


def run_main(separate_flnc_pickle_file, nfl_contigset,
             cluster_chunk_pickle, partial_chunk_pickle, polish_chunk_pickle,
             max_nchunks):
    """
    Create chunk tasks for ICE, ice_partial and ice_polish, write each
    set of chunk tasks to output pickles.
    """
    log.info("Getting all binned flnc files from %s", separate_flnc_pickle_file)
    flnc_fns = SeparateFLNCBase.convert_pickle_to_sorted_flnc_files(separate_flnc_pickle_file)
    log.debug("Binned flnc files are: %s", ", ".join(flnc_fns))

    # Number of ICE chunk tasks equals to number of bins.
    n_bins = len(flnc_fns)
    assert n_bins > 0

    log.info("max_nchunks: %s", max_nchunks)
    n_nfl_chunks = max(1, int(max_nchunks))

    out_dir = op.dirname(cluster_chunk_pickle)
    nfl_chunk_json = op.join(out_dir, 'nfl_chunk.json')
    chunked_nfl_files = chunk_contigset(in_file=nfl_contigset, n_chunks=n_nfl_chunks,
                                        out_dir=out_dir, out_chunk_json=nfl_chunk_json)

    create_cluster_pickle(flnc_files=flnc_fns,
                          out_pickle=cluster_chunk_pickle)
    create_partial_pickle(flnc_files=flnc_fns,
                          chunked_nfl_files=chunked_nfl_files,
                          out_pickle=partial_chunk_pickle)

    # Total number of flnc reads in all bins
    n_reads_in_bins = n_reads_in_contigsets(flnc_fns)
    sum_n_flnc_reads = sum(n_reads_in_bins)
    n_polish_chunks_in_bins = [max(1, int(n * max_nchunks / (1.0 * sum_n_flnc_reads)))
                               for n in n_reads_in_bins]
    create_polish_pickle(n_polish_chunks_in_bins=n_polish_chunks_in_bins,
                         flnc_files=flnc_fns,
                         out_pickle=polish_chunk_pickle)

    # Make a soft link of nfl_contigset in the same directory as separate_flnc.pickle
    # for users' convenience
    dst_nfl_contigset = op.join(op.dirname(separate_flnc_pickle_file),
                                "isoseq_nfl.contigset.xml")
    log.info("Making a soft link of %s to %s.", nfl_contigset, dst_nfl_contigset)
    ln(nfl_contigset, dst_nfl_contigset)


def _args_run(args):
    """args run"""
    raise NotImplementedError()


def _rtc_runner(rtc):
    """resolved tool contract runner."""
    max_nchunks = Constants.DEFAULT_NCHUNKS
    if hasattr(rtc.task, 'max_nchunks'):
        max_nchunks = int(rtc.task.max_nchunks)
    return run_main(separate_flnc_pickle_file=rtc.task.input_files[0],
                    nfl_contigset=rtc.task.input_files[1],
                    cluster_chunk_pickle=rtc.task.output_files[0],
                    partial_chunk_pickle=rtc.task.output_files[1],
                    polish_chunk_pickle=rtc.task.output_files[2],
                    max_nchunks=max_nchunks)


def main():
    """Main"""
    argv = sys.argv
    mp = get_contract_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           _args_run,
                           _rtc_runner,
                           log,
                           setup_log)


if __name__ == '__main__':
    sys.exit(main())
