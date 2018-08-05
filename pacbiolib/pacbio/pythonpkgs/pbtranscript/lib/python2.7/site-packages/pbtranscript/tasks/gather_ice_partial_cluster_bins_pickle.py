"""
Given ChunkTasksPickle of PartialChunkTask objects,
this task merges ice_partial chunked output nfl pickles.
"""

import logging
import sys
from itertools import groupby

from pbcommand.models import FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log

from pbtranscript.ice.IceUtils import combine_nfl_pickles
from pbtranscript.ice.IceFiles import IceFiles
from pbtranscript.tasks.TPickles import ChunkTasksPickle, PartialChunkTask


log = logging.getLogger(__name__)


class Constants(object):
    """Constants used for TOOL_ID"""
    TOOL_ID = "pbtranscript.tasks.gather_ice_partial_cluster_bins_pickle"
    DRIVER = "python -m %s --resolved-tool-contract " % TOOL_ID
    PARSER_DESC = __doc__
    VERSION = "0.1.0"


def get_contract_parser():
    """Tool contract should have the following inputs and outputs.
    Input:
        idx 0 - ChunkTasksPickle of PartialChunkTask objects.
        idx 1 - senitel txt file
    Output:
        idx 0 - gather_nfl_pickles_done.txt, a sential file which does nothing
                but connect this task and the subsequence task ice_polish
    """
    p = get_pbparser(tool_id=Constants.TOOL_ID,
                     version=Constants.VERSION,
                     name=Constants.TOOL_ID,
                     description=__doc__,
                     driver_exe=Constants.DRIVER)
    p.add_input_file_type(FileTypes.PICKLE, "partial_chunks_pickle", "Pickle In",
                          "Partial chunks pickle file") # input 0
    p.add_input_file_type(FileTypes.TXT, "partial_sentinel_in", "Sentinel In",
                          "Setinel file") # input idx 1
    p.add_output_file_type(FileTypes.TXT, "gather nfl pickles done txt",
                           name="Gather nfl pickles Done Txt file",
                           description="Gather nfl pickles Done Txt file.",
                           default_name="gather_ice_partial_pickles_done")
    return p


def args_runner(args):
    """args runner"""
    raise NotImplementedError()


def resolved_tool_contract_runner(rtc):
    """Given resolved tool contract, run"""
    p = ChunkTasksPickle.read(rtc.task.input_files[0])
    p.sorted_by_attr(attr='cluster_bin_index')
    assert all([isinstance(task, PartialChunkTask) for task in p])

    with open(rtc.task.output_files[0], 'w') as writer:
        for i, group in groupby(p, lambda x: x.cluster_bin_index):
            gs = [g for g in group]
            nfl_pickles_of_bin_i = [g.nfl_pickle for g in gs]
            out_pickle = IceFiles(prog_name="", root_dir=gs[0].cluster_out_dir,
                                  no_log_f=True).nfl_all_pickle_fn
            log.info("Combining nfl pickles of cluster bin %s.", str(i))
            log.debug("nfl pickles are: %s.", (", ".join(nfl_pickles_of_bin_i)))
            log.debug("Output merged nfl pickle is %s.", out_pickle)
            combine_nfl_pickles(splitted_pickles=nfl_pickles_of_bin_i, out_pickle=out_pickle)
            writer.write("Merge nfl pickles of cluster bin %s DONE: %s\n" %
                         (i, out_pickle))


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
