
# TODO modernize this or delete

"""
Streaming Subreading Filtering Report

CSV file format is:

MovieName,HoleNumber,Start,End,Length,PassedFilter
m120404_091356_42139_c000325362550000001500000112311376_s1_p0,8,426,664,238,0
m120404_091356_42139_c000325362550000001500000112311376_s1_p0,8,705,1035,330,0
m120404_091356_42139_c000325362550000001500000112311376_s1_p0,8,1079,1401,322,0

The report generates:

1. 4 attributes
2. One plot in a plot group


Attributes

    total = readlengths.shape[0]
    nbases = np.sum(readlengths)
    mean = readlengths.mean()
    n50 = compute_n50(readlengths)

    attr_nbases = Attribute('filter_subread_nbases',
                            nbases, "Total Number of Bases")

    attr_total = Attribute('filter_subread_nreads', total,
                           "Number of Reads")

    attr_mean = Attribute('filter_subread_mean', int(mean),
                          name="Mean Subread length")

    attr_n50 = Attribute('filter_subread_n50', n50, name="N50")
"""

import os
import os.path as op
import sys
import argparse
import logging
import functools
import itertools
import math

import numpy as np

from pbcommand.models.report import Report, Attribute
from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbreports.io.validators import validate_dir, validate_file
from pbreports.plot.helper import get_green
from pbreports.util import compute_n50
from pbreports.report.streaming_utils import (PlotViewProperties,
                                              to_plot_groups,
                                              custom_subread_length_histogram)
from pbreports.model.aggregators import (BaseAggregator, SumAggregator,
                                         HistogramAggregator,
                                         MaxAggregator, MeanAggregator,
                                         CountAggregator)
from pbreports.io.specs import *

log = logging.getLogger(__name__)
__version__ = '1.2'


class Constants(object):
    R_ID = "filter_subread"

    # Table
    #T_MY_ID = ''

    # Plot Groups
    PG_SUBREAD_LENGTH = 'subread'

    # Plot
    P_POST_FILTER = 'post_filter'

    # Attributes
    A_NBASES = 'filter_subread_nbases'
    A_NREADS = 'filter_subread_nreads'
    A_MEAN = 'filter_subread_mean'
    A_N50 = 'filter_subread_n50'

    I_FILTER_SUBREADS_HIST = "filtered_subread_report.png"

spec = load_spec(Constants.R_ID)


class _BaseSubreadFilterException(Exception):
    pass


class NoSubreadsFound(Exception):
    pass


class NoSubreadsPassedFilter(_BaseSubreadFilterException):
    pass


class _BaseHistogram(BaseAggregator):

    def __init__(self, dx=100.0, nbins=1000, dtype=np.int32):
        """
        :param dx: float, int
        :param nbins: int
        :param dtype: type used to create the np array
        """
        self.dx = dx
        self.dtype = dtype
        self.bins = np.zeros(nbins, dtype=np.int32)

    @property
    def nbins(self):
        return len(self.bins)

    @property
    def bin_edges(self):
        """Used for plotting cdf

        ds = zip(self.bin_edges, self.bins)

        c = to_cdf(ds)

        plot(self.bin_edges, c)
        plot(self.bin_edges, self.bins)

        """
        return [self.dx * i for i in xrange(self.nbins)]

    def apply(self, record):
        """This will be readlengths"""
        raise NotImplemented

    def __repr__(self):
        x = self.dx * self.nbins
        _d = dict(k=self.__class__.__name__,
                  d=self.dx,
                  n=self.nbins,
                  x=x,
                  i=0)
        return "<{k} dx:{d} nbins:{n} min:{i} max:{x} >".format(**_d)


class SubreadLengthHistogram(_BaseHistogram):

    def apply(self, record):
        """This will be readlengths"""
        i = int(math.ceil(record.length / self.dx))
        # add
        try:
            self.bins[i] += 1
        except IndexError as e:
            log.error(e)
            x = self.dx * self.nbins
            _d = dict(v=record.length, i=i, d=self.dx, x=x, n=self.nbins)
            log.error(
                "Max value {v} dx:{d} nbins{n} max value {x}".format(**_d))


class MeanSubreadLengthAggregator(MeanAggregator):
    pass


class N50Aggregator(BaseAggregator):

    def __init__(self, record_field, values=()):
        self.record_field = record_field
        # this is not constant memory
        self.values = values if values else []

    def apply(self, record):
        v = getattr(record, self.record_field)
        self.values.append(v)

    @property
    def n50(self):
        return compute_n50(self.values)

    def __repr__(self):

        _d = dict(k=self.__class__.__name__,
                  v=len(self.values),
                  n=self.n50)
        return "<{k} n50:{n} nvalues:{v} >".format(**_d)

    @property
    def attribute(self):
        return self.n50


class Record(object):

    def __init__(self, movie_name, hole_number, start, end, length, passed_filter):
        self.movie_name = movie_name
        self.hole_number = hole_number
        self.start = start
        self.end = end
        self.length = length
        self.passed_filter = passed_filter

    def __repr__(self):
        _d = dict(m=self.movie_name,
                  h=self.hole_number,
                  s=self.start,
                  e=self.end,
                  l=self.length,
                  p=self.passed_filter)
        return "<{k} {m} {h} {s} {e} {l} {p}>".format(**_d)


def to_record(line):
    to_b = lambda x: True if x == '1' else False
    to_c = [str, int, int, int, int, to_b]

    try:
        a = [v(k) for k, v in zip(line.strip().split(','), to_c)]
        r = Record(*a)
        return r
    except (ValueError, TypeError) as e:
        msg = "Unable to process line '{l}'".format(l=line)
        sys.stderr.write(msg + "\n")
        raise


def _filter_record(filter_func, record):
    """Returns Bool"""
    return filter_func(record)


def null_filter(record):
    return True


def _multi_filter_record(filter_funcs, record):
    """Returns Bool"""
    for filter_func in filter_funcs:
        if not filter_func(record):
            # Bail out at the first chance
            return False

    return True


def _apply(filter_funcs, aggregators, record):
    """
    Run the filters and call apply method on the aggregator if
    the record passes filtering.

    The should be used with functools.partial


    my_func = functools.partial(_apply, [lambda r: True], [])

    This will be my_func(record)

    """
    if not isinstance(filter_funcs, (list, tuple)):
        filter_funcs = [filter_funcs]

    if not isinstance(aggregators, (list, tuple)):
        aggregators = [aggregators]

    if _multi_filter_record(filter_funcs, record):
        for aggregator in aggregators:
            aggregator.apply(record)


def applyer(row_to_record_func, iterable, funcs):
    """

    :params funcs: list of funcs that operate on Records

    """
    for it in iterable:
        record = row_to_record_func(it)
        for func in funcs:
            func(record)

        del record


def _to_attributes(nreads, nbases, mean_readlength, n50):
    """
    Returns a list of attributes
    """

    attr_nbases = Attribute(Constants.A_NBASES, nbases)

    attr_total = Attribute(Constants.A_NREADS, nreads)

    attr_mean = Attribute(Constants.A_MEAN, int(mean_readlength))

    attr_n50 = Attribute(Constants.A_N50, n50)

    attributes = [attr_mean, attr_n50, attr_nbases, attr_total]

    for attribute in attributes:
        log.debug(attribute)

    return attributes


def to_report(filtered_csv, output_dir, dpi=72, thumb_dpi=20):
    """
    Run Report
    """
    validate_file(filtered_csv)
    validate_dir(output_dir)

    aggregators = {'nbases': SumAggregator('length'),
                   'nreads': CountAggregator('length'),
                   'mean_subreadlength': MeanSubreadLengthAggregator('length'),
                   'max_readlength': MaxAggregator('length'),
                   'n50': N50Aggregator('length'),
                   'readlength_histogram': HistogramAggregator('length', 0, 100, nbins=10000),
                   'subread': SubreadLengthHistogram(dx=100)}

    passed_filter = lambda record: record.passed_filter is True

    passed_filter_func = functools.partial(
        _apply, [passed_filter], aggregators.values())

    all_subread_aggregators = {'raw_nreads': SumAggregator('length'),
                               'max_raw_readlength': MaxAggregator('length'),
                               'raw_readlength_histogram': HistogramAggregator('length', 0, 100, nbins=10000)}

    all_filter_func = functools.partial(
        _apply, [null_filter], all_subread_aggregators.values())

    funcs = [passed_filter_func, all_filter_func]

    with open(filtered_csv, 'r') as f:
        # read in header
        header = f.readline()
        # validate_header(header)
        applyer(to_record, f, funcs)

    for aggregator in itertools.chain(aggregators.values(), all_subread_aggregators.values()):
        log.info(aggregator)

    # Check if any reads are found
    if all_subread_aggregators['raw_nreads'].attribute == 0:
        raise NoSubreadsFound(
            "No subreads found in {f}".format(f=filtered_csv))

    # Now check
    if aggregators['nreads'].attribute == 0:
        msg = "No subreads passed the filter in {f}.".format(f=filtered_csv)
        raise NoSubreadsPassedFilter(msg)

    # this is where you change the plotting options
    plot_view = PlotViewProperties(Constants.P_POST_FILTER,
                                   Constants.PG_SUBREAD_LENGTH,
                                   custom_subread_length_histogram,
                                   Constants.I_FILTER_SUBREADS_HIST,
                                   xlabel=get_plot_xlabel(
                                       spec, Constants.PG_SUBREAD_LENGTH, Constants.P_POST_FILTER),
                                   ylabel="Subreads",
                                   rlabel="bp > Subread Length",
                                   thumb="filtered_subread_report_thmb.png",
                                   use_group_thumb=True,
                                   color=get_green(3),
                                   edgecolor=get_green(2))

    view_config_d = {'post_filter': plot_view}
    id_aggregators = {'post_filter': aggregators['subread']}

    plot_groups = to_plot_groups(view_config_d, output_dir, id_aggregators)

    to_a = lambda n: aggregators[n].attribute

    attributes = _to_attributes(to_a('nreads'),
                                to_a('nbases'),
                                to_a('mean_subreadlength'),
                                to_a('n50'))

    report = Report(Constants.R_ID,
                    plotgroups=plot_groups,
                    attributes=attributes)

    log.debug(str(report))

    return spec.apply_view(report)


def args_runner(args):
    log.info("Starting {f} version {v} report generation".format(
        f=__file__, v=__version__))
    report = to_report(args.filter_summary_csv, args.output, dpi=args.dpi)
    report.write_json(args.report)
    return 0


def get_parser():
    """Old Usage in pbpy:

    usage = filter_subread.py --debug filtered_subread_summary.csv --output /path/to/outputDir --report /path/to/outputDir/junk3.json

    filtered_subread_summary.csv
    """
    desc = ""
    parser = get_default_argparser_with_base_opts(
        version=__version__, description=__doc__)
    parser.add_argument("filter_summary_csv",
                        help="Path to Filter Subread Summary CSV file.",
                        type=validate_file)
    parser.add_argument('-o', '--output', dest='output', default=os.getcwd(),
                        type=validate_dir,
                        help='Output directory to write to Subread Hist plots to.')
    parser.add_argument('--dpi', type=int, dest='dpi', default=72,
                        help="dpi (dots/inch) for plots that were generated.")
    parser.add_argument('-r', '--report', dest='report', default=None,
                        help="Write the Json report to disk.")
    return parser


def main(argv=sys.argv[1:]):
    """Main point of Entry"""
    return pacbio_args_runner(
        argv=argv,
        parser=get_parser(),
        args_runner_func=args_runner,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main())
