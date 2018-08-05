
"""
Separate wrapper for HGAP version of mapping stats report.
"""

import logging
import sys
import os
import os.path as op

from pbcommand.models import get_pbparser, FileTypes

from pbreports.report.mapping_stats import *
from pbreports.report.mapping_stats import Constants as BaseConstants
from pbreports.io.specs import *

__version__ = "0.1"
log = logging.getLogger(__name__)


class Constants(BaseConstants):
    TOOL_ID = "pbreports.tasks.mapping_stats_hgap"
    DRIVER_EXE = "python -m pbreports.report.mapping_stats_hgap --resolved-tool-contract"
    R_ID = "mapping_stats_hgap"

spec = load_spec(Constants.R_ID)


def to_report(alignment_file, output_dir, subreads_file):
    return spec.apply_view(MappingStatsCollector(alignment_file, subreads_file).to_report(output_dir, Constants.R_ID))


def _args_runner(args):
    return run_and_write_report(args.alignment_file,
                                args.report_json,
                                report_func=to_report,
                                subreads_file=args.subreads_file)


def _resolved_tool_contract_runner(resolved_contract):
    alignment_path = resolved_contract.task.input_files[0]
    subreads_path = resolved_contract.task.input_files[1]
    output_report = resolved_contract.task.output_files[0]
    return run_and_write_report(alignment_path, output_report,
                                report_func=to_report,
                                subreads_file=subreads_path)


def _get_parser():
    parser = get_pbparser(Constants.TOOL_ID, __version__,
                          spec.title, __doc__,
                          Constants.DRIVER_EXE)
    parser.add_input_file_type(FileTypes.DS_ALIGN_CCS, "alignment_file",
                               "ConsensusAlignment XML DataSet",
                               "BAM, SAM or ConsensusAlignment DataSet")
    parser.add_input_file_type(FileTypes.DS_SUBREADS, "subreads_file",
                               "Subreads XML DataSet",
                               "Unmapped BAM or Subreads DataSet")
    parser.add_output_file_type(FileTypes.REPORT, "report_json",
                                "Mapping Statistics Report",
                                "Summary of alignment results",
                                default_name=Constants.R_ID)
    return parser


if __name__ == '__main__':
    sys.exit(main(get_parser_func=_get_parser,
                  args_runner_func=_args_runner,
                  rtc_runner_func=_resolved_tool_contract_runner))
