"""Implements an interval tree for looking of alignments that
overlap a reference region.

Mostly taken from Erik Garrison's MIT-Licensed C++ version, available
at https://github.com/ekg/intervaltree
"""

import collections

Interval = collections.namedtuple("Interval", ['start', 'stop'])


class IntervalTree(object):
    """Stores mapped reads as intervals that can be efficiently queried."""

    def __init__(self, ivals, depth=16, minbucket=64, leftextent=0,
                 rightextent=0, maxbucket=512):
        """Recursively create the interval tree from the Intervals in ivals."""

        self.intervals = []
        self.left = None
        self.right = None
        self.center = None

        depth -= 1
        if depth == 0 or (len(ivals) < minbucket and len(ivals) < maxbucket):
            ivals.sort(key=lambda x: x.start)
            self.intervals = ivals
        else:
            if leftextent == 0 and rightextent == 0:
                ivals.sort(key=lambda x: x.start)

            leftp, rightp, centerp = 0, 0, 0

            if leftextent or rightextent:
                leftp = leftextent
                rightp = rightextent
            else:
                leftp = ivals[0].start
                rightp = max(k.stop for k in ivals)

            centerp = ivals[len(ivals) // 2].start
            self.center = centerp

            lefts = []
            rights = []

            for ival in ivals:
                if ival.stop < self.center:
                    lefts.append(ival)
                elif ival.start > self.center:
                    rights.append(ival)
                else:
                    self.intervals.append(ival)

            if lefts:
                self.left = IntervalTree(lefts, depth, minbucket, leftp,
                                         centerp)
            if rights:
                self.right = IntervalTree(rights, depth, minbucket, centerp,
                                          rightp)

    def find_overlapping(self, start, stop, overlapping):
        """Get a list of intervals that overlap the interval defined
        by [start, stop].
        """
        if self.intervals and not stop < self.intervals[0].start:
            for interval in self.intervals:
                if interval.stop > start and interval.start < stop:
                    overlapping.append(interval)

        if self.left and start <= self.center:
            self.left.find_overlapping(start, stop, overlapping)

        if self.right and stop >= self.center:
            self.right.find_overlapping(start, stop, overlapping)
