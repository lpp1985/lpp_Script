
# TODO modernize or delete

"""
Counts the number of movies and cells in the input dataset.
"""

import os
import os.path as op
import sys
import argparse
import logging

from pbcommand.models.report import Report, Attribute
from pbcommand.models import FileTypes, get_pbparser
from pbcommand.cli import pbparser_runner
from pbcommand.utils import setup_log
from pbcore.io import openDataSet, BamReader

from pbreports.util import movie_to_cell, path_to_movie
from pbreports.io.specs import *

log = logging.getLogger(__name__)

__version__ = '2.0'


class Constants(object):
    TOOL_ID = "pbreports.tasks.overview"
    TOOL_NAME = "Overview report"
    DRIVER_EXE = "python -m pbreports.report.overview --resolved-tool-contract "
    R_ID = "overview"
    A_NCELLS = "ncells"
    A_NMOVIES = "nmovies"

spec = load_spec(Constants.R_ID)


def run(dataset_file):
    """Reads in the input.fofn and counts movies and cells. Outputs in XML."""

    with openDataSet(dataset_file) as ds:
        movies = None
        movies = set([])
        for file_name in ds.toExternalFiles():
            if type(ds).__name__ == "HdfSubreadSet":
                movies.add(path_to_movie(file_name))
            else:
                with BamReader(file_name) as bam:
                    for rg in bam.peer.header["RG"]:
                        movies.add(rg["PU"])
        cells = set([movie_to_cell(movie) for movie in movies])
        ncells_attr = Attribute(Constants.A_NCELLS, len(cells))
        nmovies_attr = Attribute(Constants.A_NMOVIES, len(movies))
        attrs = [ncells_attr, nmovies_attr]
        report = Report(Constants.R_ID, attributes=attrs)
        return spec.apply_view(report)


def make_report(input_ds, output_json):
    report = run(input_ds)
    report.write_json(output_json)
    return 0


def args_runner(args):
    return make_report(input_ds=args.input_ds, output_json=args.output_json)


def resolved_tool_contract_runner(rtc):
    return make_report(
        input_ds=rtc.task.input_files[0],
        output_json=rtc.task.output_files[0])


def get_parser():
    p = get_pbparser(
        tool_id=Constants.TOOL_ID,
        version=__version__,
        name=Constants.TOOL_NAME,
        description=__doc__,
        driver_exe=Constants.DRIVER_EXE)
    p.add_input_file_type(FileTypes.DS_SUBREADS, "input_ds",
                          name="SubreadSet",
                          description="Subread dataset")
    p.add_output_file_type(FileTypes.REPORT, "output_json",
                           name=spec.title,
                           description="Path to write report JSON output",
                           default_name="overview_report")
    return p


def main(argv=sys.argv):
    mp = get_parser()
    return pbparser_runner(
        argv=argv[1:],
        parser=mp,
        args_runner_func=args_runner,
        contract_runner_func=resolved_tool_contract_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == '__main__':
    sys.exit(main())
