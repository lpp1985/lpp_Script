"""
Given a list of PolishChunkTask in ChunkTasksPickle, call
IceQuiver to polish consensus isoforms.
"""

import logging
import os.path as op
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.utils import setup_log
from pbcommand.models import FileTypes

from pbtranscript.Utils import mkdir
from pbtranscript.tasks.TPickles import ChunkTasksPickle, PolishChunkTask
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser)
from pbtranscript.ClusterOptions import SgeOptions
from pbtranscript.ice.IceQuiver import IceQuiver

log = logging.getLogger(__name__)

class Constants(BaseConstants):
    """Constants used."""
    TOOL_ID = "pbtranscript.tasks.ice_polish_cluster_bins"
    DRIVER_EXE = "python -m %s --resolved-tool-contract" % TOOL_ID
    PARSER_DESC = __doc__


def get_contract_parser():
    """
    input idx 0: polish_chunk_pickle_id
    input idx 1: sentinel.txt
    input idx 2: *.subreadset.xml
    output idx 0: chunk json
    """
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    p.add_input_file_type(FileTypes.PICKLE, "polish_chunk_pickle",
                          "PICKLE", "Polish Chunk Tasks Pickle") # input idx 0
    p.add_input_file_type(FileTypes.TXT, "sentinel_in", "Sentinel In",
                          "Setinel file") # input idx 1
    p.add_input_file_type(FileTypes.DS_SUBREADS, "subreads_in", "SubreadSet In",
                          "PacBio SubreadSet") # input idx 2

    # this file does nothing other than connect pbsmrtpipe tasks
    p.add_output_file_type(FileTypes.TXT, "polish done txt",
                           name="Polish Done Txt file",
                           description="Polish Done Txt file.",
                           default_name="polish_chunks_done")
    return p


class IceQuiverRTC(IceQuiver):

    """IceQuiver Resolved tool contract runner."""

    def __init__(self, root_dir, subread_set, nproc):
        tmp_dir = op.join(root_dir, "tmp")
        mkdir(tmp_dir)
        super(IceQuiverRTC, self).__init__(
            root_dir=root_dir,
            bas_fofn=subread_set,
            fasta_fofn=None,
            sge_opts=SgeOptions(
                unique_id=12345,
                use_sge=False,
                max_sge_jobs=0,
                blasr_nproc=nproc,
                quiver_nproc=nproc),
            prog_name="IceQuiver")

    def cluster_dir(self, cid):
        """"overwrite IceQuiver.cluster_dir"""
        dir_name = IceQuiver.cluster_dir(self, cid)
        mkdir(dir_name)
        return dir_name


def args_runner(args):
    """args runner"""
    raise NotImplementedError()


def resolved_tool_contract_runner(rtc):
    """resolved tool contract runner."""
    p = ChunkTasksPickle.read(rtc.task.input_files[0])
    assert all([isinstance(task, PolishChunkTask) for task in p])
    dummy_sentinel_file = rtc.task.input_files[1]

    subread_set = rtc.task.input_files[2]
    nproc = rtc.task.nproc
    tmp_dir = rtc.task.tmpdir_resources[0].path \
            if len(rtc.task.tmpdir_resources) > 0 else None

    with open(rtc.task.output_files[0], 'w') as writer:
        for task in p:
            log.info("Running ice_polish on cluster bin %s, polish chunk %s/%s",
                     str(task.cluster_bin_index),
                     str(task.polish_index), str(task.n_polish_chunks))
            log.debug("ice_quiver root_dir is %s", task.cluster_out_dir)
            log.debug("consensus_isoforms is %s", task.consensus_isoforms_file)

            task_runner(task=task, subread_set=subread_set, nproc=nproc, tmp_dir=tmp_dir)
            writer.write("ice_polish of cluster bin %s, polish chunk %s/%s in %s is DONE.\n" %
                         (task.cluster_bin_index, task.polish_index, task.n_polish_chunks,
                          task.cluster_out_dir))

def task_runner(task, subread_set, nproc, tmp_dir):
    """
    Given a PolishChunkTask object, run
    """
    flnc_pickle = task.flnc_pickle
    if not op.exists(flnc_pickle):
        raise IOError("flnc pickle %s does not exist." % flnc_pickle)

    nfl_pickle = task.nfl_pickle
    if not op.exists(nfl_pickle):
        raise IOError("nfl pickle %s does not exist." % nfl_pickle)

    iceq = IceQuiverRTC(root_dir=task.cluster_out_dir,
                        subread_set=subread_set,
                        nproc=nproc)
    iceq.validate_inputs()
    iceq.process_chunk_i(i=task.polish_index,
                         num_chunks=task.n_polish_chunks)


def main():
    """main"""
    mp = get_contract_parser()
    return pbparser_runner(
        argv=sys.argv[1:],
        parser=mp,
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)


if __name__ == "__main__":
    sys.exit(main())
