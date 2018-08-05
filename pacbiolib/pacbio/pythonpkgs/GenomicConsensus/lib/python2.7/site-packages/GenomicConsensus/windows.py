#################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

# windows.py: logic for windows/intervals of the genome
#
#  NOTE that by convention:
#    (start, end) is an /interval/
#    (refId, start, end) is a /window/.
#  All windows/intervals use 0-based indexing and are half-open
#  (includes start, not end)
#
# Author: David Alexander

import numpy as np, math
from pbcore.io.rangeQueries import projectIntoRange
from ConsensusCore import CoveredIntervals

# TODO(lhepler): replace the above with the following:
# from ConsensusCore2 import CoveredIntervals

def intervalToPair(v):
    return (v.Begin, v.End)

def kCoveredIntervals(k, tStart, tEnd, winStart, winEnd):
    return map(intervalToPair, CoveredIntervals(k, tStart, tEnd, int(winStart), int(winEnd-winStart)))

def kSpannedIntervals(refWindow, k, start, end, minLength=0):
    """
    Find intervals in the window that are k-spanned by the reads.

    Given:
     `refWindow`: the window under consideration
     `k`: the number of reads that must span intervals to be returned
     `start`, `end`: numpy arrays of start and end coordinates for reads,
       where the extent of each read is [start, end).  Must be ordered
       so that `start` is sorted in ascending order.

    Find a maximal set of maximal disjoint intervals within
    refWindow such that each interval is spanned by at least k reads.
    Intervals are returned in sorted order, as a list of (start, end)
    tuples.

    Note that this is a greedy search procedure and may not always
    return the optimal solution, in some sense.  However it will
    always return the optimal solutions in the most common cases.
    """
    assert k >= 1
    winId, winStart_, winEnd_ = refWindow

    # Truncate to bounds implied by refWindow
    start = np.clip(start, winStart_, winEnd_)
    end = np.clip(end, winStart_, winEnd_)

    # Translate the start, end to coordinate system where
    # refWindow.start is 0.
    start = start - winStart_
    end   = end - winStart_
    winStart = 0
    winEnd   = winEnd_ - winStart_

    positions = np.arange(winEnd - winStart, dtype=int)
    coverage = projectIntoRange(start, end,
                                winStart, winEnd)
    x = -1
    y = 0
    intervalsFound = []

    while y < winEnd:
        # Step 1: let x be the first pos >= y that is k-covered
        eligible = np.flatnonzero((positions >= y) & (coverage >= k))
        if len(eligible) > 0:
            x = eligible[0]
        else:
            break

        # Step 2: extend the window [x, y) until [x, y) is no longer
        # k-spanned.  Do this by setting y to the k-th largest `end`
        # among reads covering x
        eligible = end[(start <= x)]
        eligible.sort()
        if len(eligible) >= k:
            y = eligible[-k]
        else:
            break

        intervalsFound.append((x, y))

    # Translate intervals back
    return [ (s + winStart_,
              e + winStart_)
             for (s, e) in intervalsFound
             if e - s >= minLength ]

def abut(intervals):
    """
    Abut adjacent intervals.  Useful for debugging...
    """
    output = []
    lastS = None
    lastE = None
    for (s, e) in intervals:
        if s == lastE:
            lastS, lastE = lastS, e
        else:
            if lastS is not None:
                output.append((lastS, lastE))
            lastS, lastE = s, e
    output.append((lastS, lastE))
    return output

def holes(refWindow, intervals):
    """
    Given a window and a set of disjoint subintervals, return the
    "holes", which are the intervals of the refWindow not covered by
    the given subintervals.
    """
    winId, winStart, winEnd = refWindow
    output = []
    intervals = sorted(intervals)
    lastE = winStart
    for (s, e) in intervals:
        if s > lastE:
            output.append((lastE, s))
        lastE = e
    if lastE < winEnd:
        output.append((lastE, winEnd))
    return output

def intersection(int1, int2):
    s1, e1 = int1
    s2, e2 = int2
    si, ei = max(s1, s2), min(e1, e2)
    if si < ei:
        return (si, ei)
    else:
        return None

def windowsIntersect(w1, w2):
    i1, s1, e1 = w1
    i2, s2, e2 = w2
    return (i1 == i2) and (e1 > s2) and (e2 > s1)

def subWindow(refWindow, subinterval):
    winId, winStart, winEnd = refWindow
    intS, intE = subinterval
    assert intS >= winStart
    assert intE <= winEnd
    return winId, intS, intE

def enumerateIntervals(bounds, stride):
    """
    Enumerate windows of size "stride", attempting to align window
    boundaries on multiple of stride.
    """
    def alignDown(chunk, x):
        return (x/chunk)*chunk
    def alignUp(chunk, x):
        return int(math.ceil(float(x)/chunk)*chunk)

    start, end = bounds
    roundStart = alignDown(stride, start)
    roundEnd   = alignUp  (stride, end)

    for s in xrange(roundStart, roundEnd, stride):
        roundWin = (s, s + stride)
        yield intersection(bounds, roundWin)
