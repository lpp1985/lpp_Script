
"""
Calls the ICE algorithm, which stands for 'Iteratively Clustering and Error
correction', to identify de novo consensus isoforms.
"""

import logging
import os
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptRunner import PBTranscript
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser,
                                              get_argument_parser,
                                              add_cluster_arguments)

log = logging.getLogger(__name__)

class Constants(BaseConstants):
    TOOL_ID = "pbtranscript.tasks.cluster"
    DRIVER_EXE = "python -m %s --resolved-tool-contract" % TOOL_ID
    PARSER_DESC = __doc__


def get_contract_parser():
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    add_cluster_arguments(p)
    p.tool_contract_parser.add_output_file_type(
        FileTypes.PICKLE, "pickle_fn",
        name="Clusters pickle file",
        description="Python pickle file of clusters",
        default_name="final_clusters")
    # rtc has 4 inputs:
    #    idx 0 - flnc.contigset,
    #    idx 1 - nfl.contigset,
    #    idx 2 - ccs,
    #    idx 3 - subreads
    #
    # rtc has 4 outputs:
    #    idx 0 - consensus_isoform.contigset,
    #    idx 1 - output.json,
    #    idx 2 - output.csv
    #    idx 3 - output.pickle

    return p


def args_runner(args):
    return PBTranscript(args, subCommand="cluster").start()


def resolved_tool_contract_to_args(rtc):
    args = [
        "--verbose",
        "cluster",
        "--nfl_fa", rtc.task.input_files[1],
        "--ccs_fofn", rtc.task.input_files[2],
        "--bas_fofn", rtc.task.input_files[3],
        "--blasr_nproc", str(rtc.task.nproc),
        "--quiver_nproc", str(rtc.task.nproc),
        "--summary", str(rtc.task.output_files[1]),   # output 1: JSON
        "--report", rtc.task.output_files[2],         # output 2: CSV
        "--pickle_fn", str(rtc.task.output_files[3]),  # output 3: pickle
    ]
    #if rtc.task.options.get(Constants.USE_FINER_QV_ID, False):
    #    args.append("--use_finer_qv")
    args.extend([
        "-d", os.path.dirname(rtc.task.output_files[1]),
        # NOTE old smrtpipe script has this:
        # --unique_id 91373 --cDNA_size under1k
        rtc.task.input_files[0],
        rtc.task.output_files[0],  # output 0: ContigSet
    ])
    log.info("ARGUMENTS: %s" % " ".join(args))
    return get_argument_parser().parse_args(args)


def resolved_tool_contract_runner(resolved_tool_contract):
    args = resolved_tool_contract_to_args(resolved_tool_contract)
    return args_runner(args)


def main(argv=sys.argv[1:]):
    mp = get_contract_parser()
    return pbparser_runner(
        argv=argv,
        parser=mp,
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
