
"""
Calls the ICE algorithm, which stands for 'Iteratively Clustering and Error
correction', to identify de novo consensus isoforms.
"""

import logging
import os.path as op
import os
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.models import FileTypes
from pbcommand.utils import setup_log

from pbtranscript.PBTranscriptRunner import PBTranscript
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser, add_fofn_arguments,
                                              add_cluster_summary_report_arguments,
                                              add_ice_post_quiver_hq_lq_arguments)
from pbtranscript.ClusterOptions import SgeOptions, IceQuiverHQLQOptions
from pbtranscript.ice.IceQuiver import IceQuiver
from pbtranscript.ice.IceQuiverPostprocess import IceQuiverPostprocess
from pbtranscript.ice.IceFiles import IceFiles

log = logging.getLogger(__name__)


class Constants(BaseConstants):
    TOOL_ID = "pbtranscript.tasks.ice_quiver_postprocess"
    DRIVER_EXE = "python -m pbtranscript.tasks.ice_quiver_postprocess --resolved-tool-contract"
    PARSER_DESC = __doc__


def args_runner(args):
    raise NotImplementedError()


def resolved_tool_contract_runner(rtc):
    opts = rtc.task.options
    final_pickle_fn = rtc.task.input_files[2]
    output_dir = op.dirname(final_pickle_fn)
    IceFiles.final_consensus_fa = property(
        lambda self: rtc.task.input_files[1])
    IceFiles.final_pickle_fn = property(lambda self: final_pickle_fn)
    IceFiles.nfl_all_pickle_fn = property(lambda self: rtc.task.input_files[3])
    ipq_opts = IceQuiverHQLQOptions(
        hq_isoforms_fa=rtc.task.output_files[2],
        hq_isoforms_fq=rtc.task.output_files[3],
        lq_isoforms_fa=rtc.task.output_files[4],
        lq_isoforms_fq=rtc.task.output_files[5],
        qv_trim_5=opts[Constants.QV_TRIM_FIVEPRIME_ID],
        qv_trim_3=opts[Constants.QV_TRIM_THREEPRIME_ID],
        hq_quiver_min_accuracy=opts[Constants.HQ_QUIVER_MIN_ACCURACY_ID])
    _jobs_log = op.join(output_dir, "submitted_quiver_jobs.txt")
    shell_scripts = []
    for file_name in os.listdir(op.join(output_dir, "quivered")):
        if file_name.endswith(".sh"):
            shell_scripts.append(file_name)
    with open(_jobs_log, 'w') as f:
        f.write("\n".join(["\t".join(["local", s]) for s in shell_scripts]))

    class _IceQuiverPostprocess(IceQuiverPostprocess):

        def validate_inputs(self): return True

        @property
        def submitted_quiver_jobs_log(self):
            return _jobs_log
    icep = _IceQuiverPostprocess(
        root_dir=output_dir,
        ipq_opts=ipq_opts,
        summary_fn=rtc.task.output_files[0],
        report_fn=rtc.task.output_files[1])
    icep.run()
    return 0


def get_contract_parser():
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    tcp = p.tool_contract_parser
    add_fofn_arguments(p.arg_parser.parser, bas_fofn=True,
                       tool_contract_parser=tcp)
    tcp.add_input_file_type(FileTypes.DS_CONTIG, "consensus_fa",
                            name="Consensus isoforms",
                            description="ContigSet of consensus isoforms")
    tcp.add_input_file_type(FileTypes.PICKLE, "cluster_pickle",
                            name="Clusters",
                            description="Cluster pickle file")
    tcp.add_input_file_type(FileTypes.PICKLE, "map_nofl_pickle",
                            name="Pickle file",
                            description="Pickle file for non-full-length read mapping")
    tcp.add_input_file_type(FileTypes.JSON, "json_in", name="JSON file",
                            description="Sentinel file from ice_quiver task")
    add_cluster_summary_report_arguments(p)
    add_ice_post_quiver_hq_lq_arguments(p)
    return p


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
