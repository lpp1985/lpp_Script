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

# Author: David Alexander

import numpy as np

__all__ = [ "Consensus",
            "QuiverConsensus",
            "ArrowConsensus",
            "totalLength",
            "areContiguous",
            "join" ]

class Consensus(object):
    """
    A multiple sequence consensus corresponding to a
    (reference/scaffold) coordinate region
    """
    def __init__(self, refWindow, sequence, confidence):
        assert (len(sequence) ==
                len(confidence))
        self.refWindow  = refWindow
        self.sequence   = sequence
        self.confidence = confidence

    def __cmp__(self, other):
        return cmp(self.refWindow, other.refWindow)

    #
    # Functions for calling the consensus for regions of inadequate
    # coverage
    #

    @classmethod
    def nAsConsensus(cls, refWin, referenceSequence):
        length = len(referenceSequence)
        seq = np.empty(length, dtype="S1")
        seq.fill("N")
        conf = np.zeros(length, dtype=np.uint8)
        return cls(refWin, seq.tostring(), conf)

    @classmethod
    def referenceAsConsensus(cls, refWin, referenceSequence):
        conf = np.zeros(len(referenceSequence), dtype=np.uint8)
        return cls(refWin, referenceSequence, conf)

    @classmethod
    def lowercaseReferenceAsConsensus(cls, refWin, referenceSequence):
        conf = np.zeros(len(referenceSequence), dtype=np.uint8)
        return cls(refWin, referenceSequence.lower(), conf)

    @classmethod
    def noCallConsensus(cls, noCallStyle, refWin, refSequence):
        d = { "nocall"             : cls.nAsConsensus,
              "reference"          : cls.referenceAsConsensus,
              "lowercasereference" : cls.lowercaseReferenceAsConsensus}
        factory = d[noCallStyle]
        return factory(refWin, refSequence)


class QuiverConsensus(Consensus):
    """
    A QuiverConsensus object carries an additional field, `mms`, which
    is the ConsensusCore MultiReadMutationScorer object, which can be
    used to perform some post-hoc analyses (diploid, sample mixture, etc)
    """
    def __init__(self, refWindow, sequence, confidence, mms=None):
        super(QuiverConsensus, self).__init__(refWindow, sequence, confidence)
        self.mms = mms


class ArrowConsensus(Consensus):
    """
    An ArrowConsensus object carries an additional field, `ai`, which
    is the ConsensusCore2 abstract integrator object, which can be used
    to perform some post-hoc analyses (diploid, sample mixture, etc)
    """
    def __init__(self, refWindow, sequence, confidence, ai=None):
        super(ArrowConsensus, self).__init__(refWindow, sequence, confidence)
        self.ai = ai


def totalLength(consensi):
    """
    Total length of reference/scaffold coordinate windows
    """
    return sum(cssChunk.refWindow[2] - cssChunk.refWindow[1]
               for cssChunk in consensi)

def areContiguous(refWindows):
    """
    Predicate that determines whether the reference/scaffold windows
    are contiguous.
    """
    lastEnd = None
    lastId  = None
    for refWin in refWindows:
        id, start, end = refWin
        if ((lastId is not None and id != lastId) or
            (lastEnd is not None and start != lastEnd)):
            return False
        lastEnd = end
        lastId  = id
    return True

def join(consensi):
    """
    [Consensus] -> Consensus

    String together all the consensus objects into a single consensus.
    Will raise a ValueError if the reference windows are not
    contiguous.
    """
    assert len(consensi) >= 1
    sortedConsensi = sorted(consensi)
    if not areContiguous([cssChunk.refWindow for cssChunk in sortedConsensi]):
        raise ValueError, "Consensus chunks must be contiguous"

    joinedRefWindow  = (sortedConsensi[0].refWindow[0],
                        sortedConsensi[0].refWindow[1],
                        sortedConsensi[-1].refWindow[2])
    joinedSeq        = "".join([cssChunk.sequence for cssChunk in sortedConsensi])
    joinedConfidence = np.concatenate([cssChunk.confidence for cssChunk in sortedConsensi])

    return Consensus(joinedRefWindow,
                     joinedSeq,
                     joinedConfidence)


#
# Naming convention for consensus contigs
#
def consensusContigName(referenceName, algorithmName):
    return "%s|%s" % (referenceName, algorithmName)
