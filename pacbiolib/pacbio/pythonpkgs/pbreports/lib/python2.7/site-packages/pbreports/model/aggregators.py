import math
import abc
import logging

log = logging.getLogger(__name__)


class BaseAggregator(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def apply(self, record):
        pass


class BaseAttribute(object):
    # This class is to be used for aggregators that
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def attribute(self):
        return ""


class BaseAggregatorAttribute(BaseAggregator, BaseAttribute):
    # this doesn't actually work as expected.
    pass


class CountAggregator(BaseAggregatorAttribute):

    def __init__(self, record_field, total=0):
        self.total = total
        self.record_field = record_field

    def apply(self, record):
        self.total += 1

    def __repr__(self):
        _d = dict(k=self.__class__.__name__, t=self.total, f=self.record_field)
        return "<{k} {f} total={t} >".format(**_d)

    @property
    def attribute(self):
        return self.total


class MinAggregator(BaseAggregatorAttribute):

    def __init__(self, record_field):
        self.record_field = record_field
        self.value = None

    def apply(self, record):
        v = getattr(record, self.record_field)
        if self.value is None:
            self.value = v

        if v < self.value:
            self.value = v

    def __repr__(self):
        _d = dict(k=self.__class__.__name__, t=self.value, f=self.record_field)
        return "<{k} {f} min={t} >".format(**_d)

    @property
    def attribute(self):
        return self.value


class MaxAggregator(BaseAggregatorAttribute):

    def __init__(self, record_field):
        self.record_field = record_field
        self.value = None

    def apply(self, record):
        v = getattr(record, self.record_field)
        if self.value is None:
            self.value = v

        if v > self.value:
            self.value = v

    def __repr__(self):
        _d = dict(k=self.__class__.__name__, t=self.value, f=self.record_field)
        return "<{k} {f} max={t}>".format(**_d)

    @property
    def attribute(self):
        return self.value


class SumAggregator(BaseAggregatorAttribute):

    def __init__(self, record_field, total=0):
        self.total = total
        self.record_field = record_field

    def apply(self, record):
        self.total += getattr(record, self.record_field)

    def __repr__(self):
        _d = dict(k=self.__class__.__name__,
                  t=self.total,
                  f=self.record_field)
        return "<{k} {f} total={t} >".format(**_d)

    @property
    def attribute(self):
        return self.total


class MeanAggregator(BaseAggregatorAttribute):

    def __init__(self, record_field):
        self.record_field = record_field
        self.total = 0
        self.nvalues = 0

    def apply(self, record):
        self.nvalues += 1
        v = getattr(record, self.record_field)
        self.total += v

    @property
    def mean(self):
        if self.nvalues == 0:
            return 0.0
        return self.total / self.nvalues

    def __repr__(self):
        _d = dict(n=self.nvalues, t=self.total,
                  k=self.__class__.__name__,
                  f=self.record_field,
                  m=self.mean)
        return "<{k} {f} mean={m} nvalue={n} total={t}>".format(**_d)

    @property
    def attribute(self):
        return self.mean


class HistogramAggregator(BaseAggregator):

    def __init__(self, record_field, min_value, dx, nbins=10):
        self.record_field = record_field
        self.min_value = min_value
        # bin width
        self.dx = dx
        self.bins = [0 for _ in xrange(nbins)]

    @property
    def nbins(self):
        return len(self.bins)

    @property
    def max_value(self):
        return self.nbins * self.dx

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
        """Adaptively compute the histogram. If there are not enough bins,
        more will be added."""
        v = getattr(record, self.record_field)
        # If value is larger than the current list of bins
        n = int(math.ceil(v / self.dx))

        max_v = (self.nbins - 1) * self.dx

        if v >= max_v:
            delta = v - max_v
            n_new_bins = delta / self.dx
            i = int(math.ceil(n_new_bins)) + 2
            # add more bins
            #log.info(("Adding more bins ", delta, n_new_bins, i))
            for _ in xrange(i):
                self.bins.append(0)

        #log.info((v, n, max_v ))
        #log.info("{k} {f} Adding value {v} index={n} to nbins {b} dx {x}".format(v=v, b=self.nbins, x=self.dx, n=n, f=self.record_field, k=self.__class__.__name__))

        self.bins[n] += 1

    def __repr__(self):
        _d = dict(k=self.__class__.__name__,
                  f=self.record_field,
                  n=self.min_value,
                  x=self.max_value,
                  dx=self.dx,
                  nbins=len(self.bins))
        return "<{k} {f} nbins={nbins} dx={dx} min={n} max={x} >".format(**_d)
