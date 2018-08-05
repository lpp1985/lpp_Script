"""Grab bag utils for the new streaming reports that use aggregators.

This needs to be refactored and/or consolidated with the existing
plotting code!!!
"""
import os
import logging
import functools
import types

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pbcommand.models.report import Plot, PlotGroup

from pbreports.plot.helper import get_fig_axes_lpr, get_blue, get_green

log = logging.getLogger(__name__)


class PlotViewProperties(object):
    """
    This is essentially a container for *all* view related config.
    The idea is to isolate all the config to a single layer
    """

    def __init__(self, plot_id,
                 plot_group_id,
                 plot_func,
                 image_name,
                 title=None,
                 xlabel=None,
                 ylabel=None,
                 rlabel=None,
                 color=get_green(2),
                 edgecolor=get_green(3),
                 thumb=None,
                 use_group_thumb=False,
                 plot_group_title=None):
        """

        FIXME This should be easier to use... Mixing the Plot and PlotGroup
        """
        if not isinstance(plot_func, (types.FunctionType, functools.partial)):
            _d = dict(t=type(plot_func), f=plot_func)
            raise TypeError(
                "plot_func requies a function, Got type {t} for {f}".format(**_d))

        # This plotting function must have the signature (aggregator, plot_view, output_dir)
        # and must return a fig, ax tuple
        self.plot_func = plot_func

        self.plot_id = plot_id
        self.plot_group_id = plot_group_id
        self.image_name = image_name
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        # right x axis label
        self.rlabel = rlabel
        self._thumb = thumb
        self.color = color
        self.edgecolor = edgecolor
        # Within a PlotGroup is this plot that will be used (True/False)
        self.use_group_thumb = use_group_thumb

        # Used for the plot group
        self.plot_group_title = plot_group_title

    @property
    def thumb(self):
        if self._thumb is None:
            return _to_thumbnail_path(self.image_name)
        else:
            # custom thumb name
            return self._thumb

    def _plot_func_name(self):
        if isinstance(self.plot_func, functools.partial):
            return self.plot_func.func.__name__
        else:
            return self.plot_func.__name__

    def __repr__(self):
        _d = dict(k=self.__class__.__name__,
                  i=self.plot_id,
                  n=self._plot_func_name(),
                  g=self.plot_group_id)
        return "<{k} {i} group:{g} plot_func:{n} >".format(**_d)


def plot_aggregator_histogram(a, plot_view, output_dir):
    """
    h is the histogram (the h from h, bin_edges = np.histogram(data))
    bin_edges list of

    This does NOT save the image!

    :param a: Histogram Aggregator
    :param plot_view: Instance of PlotViewProperites


    :type plot_view: PlotViewProperties
    :type a: HistogramAggregator
    :type output_dir: str
    """

    h = a.bins
    bin_edges = [a.dx * i for i in xrange(len(a.bins) + 1)]

    # need to look up values in the histogram that aren't 0 to find
    # the max and min ranges to plot over. Use the indexes to look up
    # the actual values as a.dx * i
    first_index = 0
    for i, v in enumerate(a.bins):
        if v != 0:
            first_index = i
            break

    last_index = len(a.bins) - 1
    for i, v in enumerate(a.bins[::-1]):
        if v != 0:
            last_index = i
            break

    min_x = first_index * a.dx
    max_x = (a.nbins - last_index) * a.dx

    log.debug(a)
    log.debug((len(h), len(bin_edges)))
    log.debug(("Min, Max", min_x, max_x))

    #assert len(h) == (len(bin_edges) + 1)

    #dims = (12, 9)
    fig, ax = get_fig_axes_lpr()
    # need to inspect the data to find where to the plotting range
    # and the bin width. Should this just be dx?

    # this should be in PlotViewProperties
    ax.bar(bin_edges[:-1], h, width=a.dx, color=plot_view.color,
           edgecolor=plot_view.edgecolor)

    # Custom limits
    ax.set_xlim(min_x, max_x)

    ax.set_ylabel(plot_view.ylabel)
    ax.set_xlabel(plot_view.xlabel)

    if plot_view.title is not None:
        ax.set_title(plot_view.title)

    return fig, ax


def plot_aggregator_histogram_with_cdf(a, plot_view, output_dir, dpi=72):
    """

    :type a: HistogramAggregator
    :type plot_view: PlotViewProperties
    :type output_dir: str
    :type dpi: int

    """
    fig, ax = plot_aggregator_histogram(a, plot_view, output_dir)

    rax = ax.twinx()

    def to_cdf(points):
        _total = 0
        datum = []
        for x, y in points:
            _total += int(x * y)
            datum.append(_total)
        return datum

    cdf = to_cdf(zip(a.bin_edges, a.bins))
    max_cdf = max(cdf)
    sdf = [max_cdf - i for i in cdf]

    # Plot the data
    rax.plot(a.bin_edges, sdf, 'k')

    # get the min and max values from the histogram. This requires a lookup
    # to get the index of the histogram, then a lookup to the bin_edges
    # to get the actual value.
    first_index = 0
    for i, v in enumerate(a.bins):
        if v != 0:
            first_index = i
            break

    last_index = len(a.bins) - 1
    for i, v in enumerate(a.bins[::-1]):
        if v != 0:
            last_index = i
            break

    min_x = first_index * a.dx
    max_x = (a.nbins - last_index) * a.dx

    rax.set_xlim(min_x, max_x)
    if plot_view.rlabel is not None:
        rax.set_ylabel(plot_view.rlabel)

    return fig, ax


def _to_thumbnail_path(path):
    return path.replace(".png", "_thumb.png")


def _generate_plot(f, aggregator, plot_view, output_dir):
    fig, ax = f(aggregator, plot_view, output_dir)
    return fig, ax

generate_plot = functools.partial(_generate_plot, plot_aggregator_histogram)
generate_plot_with_cdf = functools.partial(
    _generate_plot, plot_aggregator_histogram_with_cdf)


def _custom_histogram_with_cdf(custom_rlabel, cdf_max_threshold, aggregator, plot_view, output_dir):
    """
    Custom histogram which will adjust the raxis to MB


    :param cdf_max_threshold: max threshold to use
    :param custom_rlabel: Label that will be used if the max threshold is exceed
    :param aggregator: HistogramAggregator instance
    :param plot_view: PlotView instance
    :param output_dir: Output directory to write images to

    :type aggregator: HistogramAggregator
    :type custom_rlabel: str
    :type cdf_max_threshold: int
    :type plot_view: PlotViewProperties

    """
    fig, ax = plot_aggregator_histogram(aggregator, plot_view, output_dir)

    rax = ax.twinx()

    def to_cdf(points):
        _total = 0
        datum = []
        for x, y in points:
            _total += int(x * y)
            datum.append(_total)
        return datum

    cdf = to_cdf(zip(aggregator.bin_edges, aggregator.bins))
    max_cdf = max(cdf)
    sdf = [max_cdf - i for i in cdf]

    use_custom_label = False
    if max_cdf >= cdf_max_threshold:
        use_custom_label = True
        tmp_cdf = [i / cdf_max_threshold for i in cdf]
        cdf = tmp_cdf
        max_cdf = max(cdf)
        sdf = [max_cdf - i for i in cdf]

    # Plot the data
    rax.plot(aggregator.bin_edges, sdf, 'k')

    # get the min and max values from the histogram. This requires a lookup
    # to get the index of the histogram, then a lookup to the bin_edges
    # to get the actual value.
    first_index = 0
    for i, v in enumerate(aggregator.bins):
        if v != 0:
            first_index = i
            break

    last_index = len(aggregator.bins) - 1
    for i, v in enumerate(aggregator.bins[::-1]):
        if v != 0:
            last_index = i
            break

    min_x = first_index * aggregator.dx
    max_x = (aggregator.nbins - last_index) * aggregator.dx

    rax.set_xlim(min_x, max_x)
    if plot_view.rlabel is not None:
        if use_custom_label:
            rax.set_ylabel(custom_rlabel)
        else:
            log.debug("CDF max threshold execeeded ({t}) Setting custom label to '{r}'".format(
                r=custom_rlabel, t=cdf_max_threshold))
            rax.set_ylabel(plot_view.rlabel)

    return fig, ax


custom_read_length_histogram = functools.partial(
    _custom_histogram_with_cdf, 'Mb > Read Length', 1000000)
custom_read_accuracy_histogram = functools.partial(
    _custom_histogram_with_cdf, 'Mb > Read Quality', 1000000)
custom_subread_length_histogram = functools.partial(
    _custom_histogram_with_cdf, 'Mb > Subread Length', 1000000)


def to_plot_groups(view_config_d, output_dir, id_to_aggregators):
    """
    How to handle custom rendering?

    For example, changing raxis to display MB instead of Bases.

    Make custom func on plot view?
    """

    plots_by_group_id = {}

    # Messy way to set the PlotGroup
    plot_group_id_to_title = {}
    plot_group_id_to_thumb = {}

    # this is the plot_id
    for id_, aggregator in id_to_aggregators.iteritems():
        plot_view = view_config_d[id_]
        # log.debug(plot_view)
        # log.debug(pformat(plot_view.__dict__))

        log.info("creating plot with func {f}".format(f=plot_view.plot_func))
        fig, ax = plot_view.plot_func(aggregator, plot_view, output_dir)

        log.debug("Saving image to {i}".format(i=plot_view.image_name))
        try:
            fig.tight_layout()
        except AttributeError as e:  # FIXME bug 25872
            log.warn("figure.tight_layout() not available")
            log.warn(str(e))
        except ValueError as e:
            log.error(str(e))
        fig.savefig(os.path.join(output_dir, plot_view.image_name))

        # Always write a thumb
        thumb_path = os.path.join(output_dir, plot_view.thumb)
        log.debug("Saving thumb to {t}".format(t=thumb_path))
        fig.savefig(thumb_path, dpi=20)
        plt.close(fig)

        # these are relative paths
        plot = Plot(id_, plot_view.image_name, title=plot_view.title,
                    caption=plot_view.title, thumbnail=plot_view.thumb)

        if plot_view.plot_group_id in plots_by_group_id:
            plots_by_group_id[plot_view.plot_group_id].append(plot)
        else:
            plots_by_group_id[plot_view.plot_group_id] = [plot]

        # Set Plot Group thumbnail in a goofy way
        if plot_view.use_group_thumb:
            # this must be relative
            plot_group_id_to_thumb[plot_view.plot_group_id] = plot_view.thumb

        if plot_view.plot_group_title is not None:
            plot_group_id_to_title[
                plot_view.plot_group_id] = plot_view.plot_group_title

    plot_groups = []
    for plot_group_id, plots in plots_by_group_id.iteritems():
        # This is awkward
        title = plot_group_id_to_title.get(plot_group_id, None)
        thumb = plot_group_id_to_thumb.get(plot_group_id, None)
        plot_group = PlotGroup(plot_group_id, plots=plots, thumbnail=thumb,
                               title=title)
        plot_groups.append(plot_group)

    return plot_groups


def get_percentile(h, bin_edges, percentile):
    """
    h and bin_edges are the same structures used from np.histogram

    :param percentile: int Percentile to compute
    :type percentile: int

    The bin width is assumed to be constant.
    """
    assert (percentile >= 0) and (percentile <= 100)

    total_integral = 0
    dx = bin_edges[1] - bin_edges[0]
    for x, y in zip(h, bin_edges[:-1]):
        v = x * dx
        total_integral += v

    max_integral = total_integral * (percentile / 100.0)
    assert max_integral <= total_integral

    total = 0
    for x, y in zip(h, bin_edges[:-1]):
        v = x * dx
        total += v
        if total >= max_integral:
            # return x, y, max_integral
            return y

    # should never get here
    raise ValueError("Unable to compute percentile {n}".format(n=percentile))
