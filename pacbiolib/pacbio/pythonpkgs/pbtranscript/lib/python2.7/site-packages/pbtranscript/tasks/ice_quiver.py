
"""
Calls the ICE algorithm, which stands for 'Iteratively Clustering and Error
correction', to identify de novo consensus isoforms.
"""

import cPickle
import logging
import os.path as op
import json
import os
import sys

from pbcommand.cli.core import pbparser_runner
from pbcommand.utils import setup_log
from pbcommand.models.report import Report
from pbcommand.models import FileTypes

from pbtranscript.PBTranscriptRunner import PBTranscript
from pbtranscript.PBTranscriptOptions import (BaseConstants,
                                              get_base_contract_parser, add_fofn_arguments)
from pbtranscript.ClusterOptions import SgeOptions, IceQuiverHQLQOptions
from pbtranscript.ice.IceQuiver import IceQuiver
from pbtranscript.ice.IceQuiverAll import add_ice_quiver_all_arguments
from pbtranscript.ice.IceFiles import IceFiles

log = logging.getLogger(__name__)

class Constants(BaseConstants):
    TOOL_ID = "pbtranscript.tasks.ice_quiver"
    DRIVER_EXE = "python -m pbtranscript.tasks.ice_quiver --resolved-tool-contract"
    PARSER_DESC = __doc__


class IceQuiverRTC(IceQuiver):

    def __init__(self, root_dir, subread_set, nproc):
        tmp_dir = op.join(root_dir, "tmp")
        if not op.isdir(tmp_dir):
            os.makedirs(tmp_dir)
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
        dir_name = IceQuiver.cluster_dir(self, cid)
        if not op.isdir(dir_name):
            os.makedirs(dir_name)
        return dir_name


def args_runner(args):
    raise NotImplementedError()


def resolved_tool_contract_runner(rtc):
    opts = rtc.task.options
    # XXX to handle chunking I am simply re-using the old i/N arguments, but
    # embedded in the input pickle instead of passed on the command line
    final_pickle_fn = rtc.task.input_files[2]
    _tmp = cPickle.load(open(final_pickle_fn, 'rb'))
    i_chunk = 0
    n_chunks = 1
    if "__chunk_i" in _tmp:
        i_chunk = _tmp['__chunk_i']
        n_chunks = _tmp['__chunk_n']
        final_pickle_fn = _tmp['pickle_file']
    output_dir = os.path.dirname(final_pickle_fn)
    IceFiles.final_consensus_fa = property(
        lambda self: rtc.task.input_files[1])
    IceFiles.final_pickle_fn = property(lambda self: final_pickle_fn)
    IceFiles.nfl_all_pickle_fn = property(lambda self: rtc.task.input_files[3])
    iceq = IceQuiverRTC(
        root_dir=output_dir,
        subread_set=rtc.task.input_files[0],
        nproc=rtc.task.nproc)
    iceq.validate_inputs()
    iceq.process_chunk_i(i=i_chunk, num_chunks=n_chunks)
    with open(rtc.task.output_files[0], 'w') as f:
        report = Report.from_simple_dict(
            report_id="isoseq_ice_quiver",
            raw_d={'n_chunks': 1},
            namespace="ice_quiver")
        f.write(report.to_json())
    return 0


def get_contract_parser():
    p = get_base_contract_parser(Constants, default_level="DEBUG")
    arg_parser = p.arg_parser.parser
    tcp = p.tool_contract_parser
    add_fofn_arguments(arg_parser, bas_fofn=True,
                       tool_contract_parser=tcp)
    tcp.add_input_file_type(FileTypes.DS_CONTIG, "consensus_fa",
                            name="ContigSet", description="Consensus isoforms")
    tcp.add_input_file_type(FileTypes.PICKLE, "cluster_pickle",
                            name="Clusters",
                            description="Cluster pickle file")
    tcp.add_input_file_type(FileTypes.PICKLE, "map_nofl_pickle",
                            name="Pickle file",
                            description="Pickle file for non-full-length read mapping")
    # XXX this file does nothing other than connect this task to
    # ice_quiver_postprocess in pbsmrtpipe
    tcp.add_output_file_type(FileTypes.JSON, "json_out", "JSON file",
                             "JSON sentinel file", default_name="quiver_out")
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
