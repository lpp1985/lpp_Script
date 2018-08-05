"""Summarize the Long Amplicon Analysis using the ZMW results"""

from collections import defaultdict
from pprint import pformat
import argparse
import logging
import csv
import os
import os.path as op
import sys

from pbcommand.models.report import Report, Table, Column
from pbcommand.models import FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.common_options import add_debug_option
from pbcommand.utils import setup_log

from pbreports.io.specs import *

log = logging.getLogger(__name__)

__version__ = '0.1.1'


class Constants(object):
    TOOL_ID = "pbreports.tasks.amplicon_analysis_input"
    DRIVER_EXE = "python -m pbreports.report.amplicon_analysis_input --resolved-tool-contract "
    R_ID = "amplicon_analysis_input"
    DATA_GOOD = "good"
    DATA_CHIMERA = "chimera"
    DATA_NOISE = "noise"
    T_R = "result_table"
    C_BC = "barcode_col"
    C_GOOD = "good"
    C_GOOD_PCT = "good_pct"
    C_CHIM = "chimera"
    C_CHIM_PCT = "chimera_pct"
    C_NOISE = "noise"
    C_NOISE_PCT = "noise_pct"

spec = load_spec(Constants.R_ID)


def create_table(summary_csv):
    """Long Amplicon Analysis results table"""

    columns = []
    columns.append(Column(Constants.C_BC))
    columns.append(Column(Constants.C_GOOD))
    columns.append(Column(Constants.C_GOOD_PCT))
    columns.append(Column(Constants.C_CHIM))
    columns.append(Column(Constants.C_CHIM_PCT))
    columns.append(Column(Constants.C_NOISE))
    columns.append(Column(Constants.C_NOISE_PCT))

    t = Table(Constants.T_R, columns=columns)

    COL_IDS = [Constants.C_GOOD, Constants.C_GOOD_PCT, Constants.C_CHIM,
               Constants.C_CHIM_PCT, Constants.C_NOISE, Constants.C_NOISE_PCT]

    def add_column(barcode_id, n_good, n_chimera, n_noise):
        pct_good = pct_chimera = pct_noise = 0
        total = n_good + n_chimera + n_noise
        if total > 0:
            pct_good = n_good / float(total)
            pct_chimera = n_chimera / float(total)
            pct_noise = n_noise / float(total)
        values = [n_good, pct_good, n_chimera, pct_chimera, n_noise, pct_noise]
        t.add_data_by_column_id(Constants.C_BC, bc_id)
        for column_id, value in zip(COL_IDS, values):
            t.add_data_by_column_id(column_id, value)

    with open(summary_csv) as csv_in:
        reader = csv.reader(csv_in, delimiter=',')
        reader.next()
        for rec in reader:
            assert len(rec) == 7, rec
            bc_id = rec[0]
            if bc_id == "All":
                continue
            add_column(bc_id, int(rec[1]), int(rec[3]), int(rec[5]))
    n_good = sum(t.get_column_by_id(Constants.C_GOOD).values)
    n_chimera = sum(t.get_column_by_id(Constants.C_CHIM).values)
    n_noise = sum(t.get_column_by_id(Constants.C_NOISE).values)
    add_column("All", n_good, n_chimera, n_noise)
    return t


def run_to_report(summary_csv):
    log.info("Generating PCR report v{v} from summary '{s}'".format(
             v=__version__,
             s=summary_csv))
    # Convert the data into a PBreports table
    table = create_table(summary_csv)
    # ids must be lowercase.
    r = Report(Constants.R_ID, tables=[table])
    return spec.apply_view(r)


def amplicon_analysis_input(summary_csv, report_json):
    log.info("Running {f} v{v}.".format(
        f=os.path.basename(__file__), v=__version__))
    report = run_to_report(summary_csv)
    log.info(pformat(report.to_dict()))
    report.write_json(report_json)
    return 0


def args_runner(args):
    amplicon_analysis_input(args.report_csv, args.report_json)
    return 0


def resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    amplicon_analysis_input(rtc.task.input_files[0], rtc.task.output_files[0])
    return 0


def _add_options_to_parser(p):
    p.add_input_file_type(
        FileTypes.CSV,
        file_id="report_csv",
        name="Consensus Summary CSV",
        description="Consensus summary CSV by barcode")
    p.add_output_file_type(
        FileTypes.REPORT,
        file_id="report_json",
        name="LAA Input Report",
        description="Summary of input amplicon quality",
        default_name="amplicon_input_report")


def add_options_to_parser(p):
    """
    API function for extending main pbreport arg parser (independently of
    tool contract interface).
    """
    p_wrap = _get_parser_core()
    p_wrap.arg_parser.parser = p
    p.description = __doc__
    add_debug_option(p)
    _add_options_to_parser(p_wrap)
    p.set_defaults(func=args_runner)
    return p


def _get_parser_core():
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        "Amplicon Analysis Input",
        __doc__,
        Constants.DRIVER_EXE)
    return p


def get_parser():
    p = _get_parser_core()
    _add_options_to_parser(p)
    return p


def main(argv=sys.argv):
    mp = get_parser()
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           setup_log)

# for 'python -m pbreports.report.amplicon_analysis_input ...'
if __name__ == "__main__":
    sys.exit(main())
