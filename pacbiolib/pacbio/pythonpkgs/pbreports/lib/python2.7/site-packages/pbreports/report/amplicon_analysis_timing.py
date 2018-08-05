
# TODO modernize this

"""Amplicon Analysis Timing Report"""

from collections import defaultdict
from pprint import pformat
import datetime
import logging
import os
import os.path as op
import re
import sys

from numpy import median

from pbcommand.models.report import Report, Table, Column
from pbcommand.models import FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.common_options import add_debug_option
from pbcommand.utils import setup_log

from pbreports.util import validate_nonempty_file
from pbreports.io.specs import *


log = logging.getLogger(__name__)

__version__ = '0.1.1'


class Constants(object):
    TOOL_ID = "pbreports.tasks.amplicon_analysis_timing"
    R_ID = "amplicon_analysis_timing"
    T_ID = "result_table"
    C_BC = "barcode_col"
    C_HOUR = "hour_col"
    C_MIN = "minute_col"
    C_SEC = "second_col"

LOG_LINE_REGEX = re.compile('^\d+-\d+-\d+\s+\d+:\d+:\d+')
LOG_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
spec = load_spec(Constants.R_ID)


def parse_log_file(log_file):

    # Parse the LAA log lines from the task's logfile
    time_stamps = defaultdict(list)
    with open(log_file) as handle:
        rollover_count = 0
        for line in handle:
            if not LOG_LINE_REGEX.match(line):
                continue

            # Parse the barcode and timestamp from each line
            datetime_part = line.split('.')[0]
            parts = line.strip().split()
            barcode = parts[2] if parts[1].endswith('Barcode') else None
            line_time = datetime.datetime.strptime(
                datetime_part, LOG_TIME_FORMAT)

            # Save the results
            if barcode is not None:
                time_stamps[barcode].append(line_time)
            time_stamps['All'].append(line_time)

    # Compute and return the difference between the min and max time for each
    # barcode
    time_extents = {}
    for k, v in time_stamps.iteritems():
        barcode = k if k == 'All' else k[1:-2]
        time_extents[barcode] = max(v) - min(v)
    return time_extents


def create_table(timings):
    """Long Amplicon Analysis Timing Result table"""

    columns = []
    columns.append(Column(Constants.C_BC))
    columns.append(Column(Constants.C_HOUR))
    columns.append(Column(Constants.C_MIN))
    columns.append(Column(Constants.C_SEC))

    t = Table(Constants.T_ID, columns=columns)

    seconds = []
    for barcode in sorted(timings):
        if barcode != 'All':
            data = timings[barcode]
            t.add_data_by_column_id(Constants.C_BC, barcode)
            t.add_data_by_column_id(Constants.C_HOUR,   data.seconds / 3600)
            t.add_data_by_column_id(Constants.C_MIN, data.seconds / 60)
            t.add_data_by_column_id(Constants.C_SEC, data.seconds)
            seconds.append(data.seconds)
    # Add the average time information
    seconds_sum = sum(seconds)
    avg_seconds = seconds_sum / len(timings)
    t.add_data_by_column_id(Constants.C_BC, 'Mean')
    t.add_data_by_column_id(Constants.C_HOUR,   avg_seconds / 3600)
    t.add_data_by_column_id(Constants.C_MIN, avg_seconds / 60)
    t.add_data_by_column_id(Constants.C_SEC, avg_seconds)
    # Add the median time information
    median_seconds = int(median(seconds))
    t.add_data_by_column_id(Constants.C_BC, 'Median')
    t.add_data_by_column_id(Constants.C_HOUR,   median_seconds / 3600)
    t.add_data_by_column_id(Constants.C_MIN, median_seconds / 60)
    t.add_data_by_column_id(Constants.C_SEC, median_seconds)
    # Add the total time information
    t.add_data_by_column_id(Constants.C_BC, 'Total')
    t.add_data_by_column_id(Constants.C_HOUR,   timings['All'].seconds / 3600)
    t.add_data_by_column_id(Constants.C_MIN, timings['All'].seconds / 60)
    t.add_data_by_column_id(Constants.C_SEC, timings['All'].seconds)

    log.debug(str(t))
    return t


def run_to_report(log_file):
    log.info("Generating Timing report v{v} from logfile '{l}'".format(
        v=__version__,
        l=log_file))

    # Parse the data into dictionaries
    timings = parse_log_file(log_file)

    # Convert the data into a PBreports table
    table = create_table(timings)

    # ids must be lowercase.
    r = Report(Constants.R_ID, tables=[table])

    return spec.apply_view(r)


def amplicon_analysis_timing(log_file, report_json):
    log.info("Running {f} v{v}.".format(
        f=os.path.basename(__file__), v=__version__))
    report = run_to_report(log_file)
    log.info(pformat(report.to_dict()))
    report.write_json(report_json)
    return 0


def args_runner(args):
    validate_nonempty_file(args.log_file)
    amplicon_analysis_timing(args.log_file, args.report_json)
    return 0


def resolved_tool_contract_runner(resolved_tool_contract):
    rtc = resolved_tool_contract
    amplicon_analysis_timing(rtc.task.input_files[0],
                             rtc.task.output_files[0])
    return 0


def _add_options_to_parser(p):
    p.add_input_file_type(
        FileTypes.TXT,
        file_id="log_file",
        name="AmpliconAnalysisLog",
        description="Amplicon Analysis Log File")
    p.add_output_file_type(
        FileTypes.JSON,
        file_id="report_json",
        name=spec.title,
        description="Timing Report JSON",
        default_name="timing_report")


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
    driver_exe = ("python -m "
                  "pbreports.report.amplicon_analysis_timing "
                  "--resolved-tool-contract ")
    p = get_pbparser(
        Constants.TOOL_ID,
        __version__,
        spec.title,
        __doc__,
        driver_exe)
    return p


def get_parser():
    p = _get_parser_core()
    _add_options_to_parser(p)
    return p


def main(argv=sys.argv):
    mp = get_parser()
    logging.basicConfig(level=logging.INFO)
    log.setLevel(logging.INFO)
    return pbparser_runner(argv[1:],
                           mp,
                           args_runner,
                           resolved_tool_contract_runner,
                           log,
                           # FIXME for some bizarre reason, calling setup_log
                           # here results in the input log file being written
                           # to when run in 'args' mode.  this has to be a bug
                           # somewhere but it might be in the core library...
                           lambda alog, **kw: alog)

# for 'python -m pbreports.report.amplicon_analysis_timing ...'
if __name__ == "__main__":
    sys.exit(main())
