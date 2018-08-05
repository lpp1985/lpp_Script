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

from __future__ import absolute_import

import math, logging, numpy as np, random
from itertools import izip
from collections import Counter
from ..utils import *
from .. import reference
from ..options import options
from ..Worker import WorkerProcess, WorkerThread
from ..ResultCollector import ResultCollectorProcess, ResultCollectorThread
from ..consensus import *
from ..variants import *


#
# --------------- Configuration ----------------------
#

class PluralityConfig(object):
    """
    Plurality configuration options
    """
    def __init__(self,
                 minMapQV=10,
                 minCoverage=3,
                 minConfidence=40,
                 diploid=False,
                 noEvidenceConsensus="nocall"):
        self.minMapQV            = minMapQV
        self.minCoverage         = minCoverage
        self.minConfidence       = minConfidence
        self.noEvidenceConsensus = noEvidenceConsensus
        self.diploid             = diploid
        self.realignHomopolymers = False # not available yet


#
# -----------  The actual algorithm code -------------
#

def pluralityConsensusAndVariants(refWindow, referenceSequenceInWindow, alns,
                                  pluralityConfig):
    """
    Compute (Consensus, [Variant]) for this window, using the given
    `alns`, by applying a straightforward column-oriented consensus
    calling algorithm.

    If the consensus cannot be called for a base, "N" will be placed
    in the consensus sequence for that position.

    If `realignHomopolymers` is True, alignment gaps will be shuffled
    in homopolymer regions in an attempt to maximize variant detection
    sensitivity (not yet implemented, and may never be).
    """
    _, refStart, refEnd = refWindow
    windowSize = refEnd - refStart
    assert len(referenceSequenceInWindow) == windowSize

    #
    # Build up these arrays in reference coordinates.
    #
    consensusSequence_   = []
    consensusFrequency_  = []
    consensusConfidence_ = []
    effectiveCoverage_   = []
    alternateAllele_     = []        # DIPLOID ONLY
    alternateFrequency_  = []        # "
    heterozygousConfidence_ = []     # "

    noCallCss = Consensus.noCallConsensus(pluralityConfig.noEvidenceConsensus,
                                          refWindow, referenceSequenceInWindow)

    baseCallsMatrix = tabulateBaseCalls(refWindow, alns)

    for j in xrange(0, windowSize):
        counter = Counter(baseCallsMatrix[:, j])
        if "" in counter: counter.pop("")

        siteEffectiveCoverage = sum(counter.itervalues())
        if ((siteEffectiveCoverage == 0) or
            (siteEffectiveCoverage < pluralityConfig.minCoverage)):
            siteConsensusFrequency = siteEffectiveCoverage
            siteConsensusSequence  = noCallCss.sequence[j]
            top2                   = None
        else:
            # Not for production code:
            top2 = counter.most_common(2)
            siteConsensusSequence, siteConsensusFrequency = top2[0]

        # Replace explicit gaps with empty string
        if siteConsensusSequence == "-":
            siteConsensusSequence = ""

        consensusSequence_.append(siteConsensusSequence)
        consensusFrequency_.append(siteConsensusFrequency)
        effectiveCoverage_.append(siteEffectiveCoverage)

        if pluralityConfig.diploid:
            if top2 and len(top2) > 1:
                siteAlternateAllele, siteAlternateFrequency = top2[1]
            else:
                siteAlternateAllele    = "N"
                siteAlternateFrequency =  0
            if siteAlternateAllele == "-":
                siteAlternateAllele = ""
            alternateAllele_.append(siteAlternateAllele)
            alternateFrequency_.append(siteAlternateFrequency)
        else:
            siteAlternateAllele    = "N"
            siteAlternateFrequency = 0

        siteConsensusConfidence, siteHeterozygousConfidence = \
            posteriorConfidences(siteEffectiveCoverage,
                                 siteConsensusFrequency,
                                 siteAlternateFrequency,
                                 diploid=pluralityConfig.diploid)
        consensusConfidence_.append(siteConsensusConfidence)
        if pluralityConfig.diploid:
            heterozygousConfidence_.append(siteHeterozygousConfidence)

    #
    # Derive variants from reference-coordinates consensus
    #
    variants = _computeVariants(pluralityConfig,
                                refWindow,
                                referenceSequenceInWindow,
                                effectiveCoverage_,
                                consensusSequence_,
                                consensusFrequency_,
                                consensusConfidence_,
                                alternateAllele_,
                                alternateFrequency_,
                                heterozygousConfidence_)
    #
    # Now we need to put everything in consensus coordinates
    #
    consensusLens = map(len, consensusSequence_)
    consensusSequence = "".join(consensusSequence_)
    consensusConfidence = np.repeat(consensusConfidence_, consensusLens)
    css = Consensus(refWindow, consensusSequence, consensusConfidence)
    return (css, variants)


def varsFromRefAndRead(refId, refPos, refBase, readSeq, **kwargs):
    """
    Compute the haploid/heterozygous Variant[s] corresponding to a
    readSeq aligned against refSeq.

    Two variant scenario:
      REF:   G
      READ: AC
        => insertion(A), substitution(G->C)

    Required: refBase != readSeq
    Returned: List of Variant objects (length one or two)
    """
    assert refBase != readSeq
    vars = []
    readBefore, readAt = readSeq[:-1], readSeq[-1:]
    if readBefore:
        # Insertion
        vars.append(Variant(refId, refPos, refPos, "", readBefore, **kwargs))
    if readAt != refBase:
        vars.append(Variant(refId, refPos, refPos+1, refBase, readAt, **kwargs))
    return vars

def varsFromRefAndReads(refId, refPos, refBase,
                        readSeq1, readSeq2, **kwargs):
    """
    Heterozygous extension of the above
    """
    assert (refBase != readSeq1) or (refBase != readSeq2)
    vars = []
    readBefore1, readAt1 = readSeq1[:-1], readSeq1[-1:]
    readBefore2, readAt2 = readSeq2[:-1], readSeq2[-1:]
    if readBefore1 or readBefore2:
        vars.append(Variant(refId, refPos, refPos, "",
                            readBefore1, readBefore2, **kwargs))
    if readAt1 != refBase or readAt2 != refBase:
        vars.append(Variant(refId, refPos, refPos+1, refBase,
                            readAt1, readAt2, **kwargs))
    return vars


def _isMixedLengthVariant(v):
    return (v.isHeterozygous and
            len(v.readSeq1) != len(v.readSeq2))

def _isSameLengthVariant(v):
    return not _isMixedLengthVariant(v)

def _computeVariants(config,
                     refWindow,
                     refSequenceInWindow,
                     coverageArray,
                     consensusArray,
                     consensusFrequencyArray,
                     consensusConfidenceArray,
                     alternateAlleleArray=None,
                     alternateAlleleFrequency=None,
                     heterozygousConfidence=None):

    refId, refStart, refEnd = refWindow
    windowSize = refEnd - refStart
    assert len(refSequenceInWindow) == windowSize
    assert len(consensusArray) == windowSize
    if config.diploid:
        assert len(alternateAlleleArray) == windowSize
        assert len(alternateAlleleFrequency) == windowSize

    vars = []
    for j in xrange(windowSize):
        refPos = j + refStart
        refBase = refSequenceInWindow[j]
        cov  = coverageArray[j]
        cssBases = consensusArray[j]
        conf = consensusConfidenceArray[j]
        cssFreq = consensusFrequencyArray[j]
        if config.diploid:
            altBases = alternateAlleleArray[j]
            altFreq  = alternateAlleleFrequency[j]
            hetConf  = heterozygousConfidence[j]
        else:
            altBases = "N"
            altFreq  = 0

        if cov < config.minCoverage: continue

        if (config.diploid and hetConf > conf):
            #
            # Diploid variant[s]?
            #
            if (hetConf >= config.minConfidence) and (refBase != "N"):
                vs = varsFromRefAndReads(refId, refPos, refBase,
                                         cssBases, altBases,
                                         confidence=hetConf, coverage=cov,
                                         frequency1=cssFreq, frequency2=altFreq)
                vars = vars + vs

        else:
            #
            # Haploid variant[s]?
            #
            if (conf >= config.minConfidence) and \
               (refBase  != cssBases)         and \
               (refBase  != "N")              and \
               (cssBases != "N")              and \
               (cssBases == "" or cssBases.isupper()):

                vs = varsFromRefAndRead(refId, refPos, refBase, cssBases,
                                        confidence=conf, coverage=cov,
                                        frequency1=cssFreq)
                vars = vars + vs

    if config.diploid:
        vars = filter(_isSameLengthVariant, vars)
    return sorted(vars)

def tabulateBaseCalls(refWindow, alns, realignHomopolymers=False):
    """
    Go through the reads and build up the structured baseCallsMatrix
    table, which tabulates the read bases occurring at each reference
    coordinate in each read.  This code is somewhat tricky, read carefully.
    """
    _, refStart, refEnd = refWindow
    windowSize = refEnd - refStart

    baseCallsMatrix = np.zeros(shape=(len(alns), windowSize), dtype="S8")

    for i, aln in enumerate(alns):
        aln = aln.clippedTo(refStart, refEnd)
        alnRef    = aln.reference(orientation="genomic")
        alnRead   = aln.read(orientation="genomic")
        if realignHomopolymers:
            alnRef, alnRead =  normalizeHomopolymerGaps(alnRef, alnRead)

        # Idea: scan through the ref, read; for each non-gap character
        # in ref, record all non-gap characters seen in read since
        # last ref gap.
        readBases = []
        accum = []
        for (refBase, readBase) in izip(alnRef, alnRead):
            if readBase != "-":
                readBases.append(readBase)
            if refBase != "-":
                basesForRefPos = "".join(readBases) if readBases else "-"
                accum.append(basesForRefPos)
                readBases = []
        s, e = (aln.referenceStart - refStart,
                aln.referenceEnd   - refStart)
        baseCallsMatrix[i, s:e] = accum
    return baseCallsMatrix

#
# ------ HACKISH POSTERIOR PROBABILITY CALCULATION ----------
#

EPS = 0.05
LOGEPS = np.log(EPS)
LOG_O_M_EPS = np.log(1-EPS)
LOG_O_M_EPS_2 = np.log((1-EPS)/2)

def posteriorConfidences(depth, cssFreq, altFreq, diploid=False, cap=40):
    """
    Return crude approximations to the posterior probabilities of the
    genotypes s_1 and s_1/s_2, where s_1 and s_2 are the observed
    consensus and alternate allele.  The assumption here is that the
    probability of the genotype being anything other that s_1, s_2, or
    s_1/s_2 is vanishingly small.  Not really a very good assumption,
    but plurality is not our real algorithm anyway.
    """
    cssFreq = cssFreq+1
    altFreq = altFreq+1
    depth = depth + 2
    cssLL_ = cssFreq*LOG_O_M_EPS + (depth-cssFreq)*LOGEPS
    altLL_ = altFreq*LOG_O_M_EPS + (depth-altFreq)*LOGEPS
    cssL_ = np.exp(cssLL_)
    altL_ = np.exp(altLL_)
    if diploid:
        hetLL_ = (cssFreq+altFreq)*LOG_O_M_EPS_2 + (depth-cssFreq-altFreq)*LOGEPS
        hetL_ = np.exp(hetLL_)
        total =  cssL_ + altL_ + hetL_
        hetProb = hetL_/total
        hetConf = -10*np.log10(1.-hetProb) if (hetProb < 1) else cap
    else:
        total =  cssL_ + altL_
        hetConf = 0
    cssProb = cssL_/total
    cssConf = -10*np.log10(1.-cssProb) if (cssProb < 1) else cap
    return int(min(cap, cssConf)), int(min(cap, hetConf))

#
# --------------  Plurality Worker class --------------------
#

class PluralityWorker(object):

    @property
    def pluralityConfig(self):
        return self._algorithmConfig

    def onStart(self):
        random.seed(42)

    def onChunk(self, workChunk):
        referenceWindow = workChunk.window
        refSeqInWindow = reference.sequenceInWindow(referenceWindow)
        logging.info("Plurality operating on %s" %
                     reference.windowToString(referenceWindow))

        if not workChunk.hasCoverage:
            noCallCss = Consensus.noCallConsensus(options.noEvidenceConsensusCall,
                                                  referenceWindow, refSeqInWindow)
            return (referenceWindow, (noCallCss, []))

        alnHits = readsInWindow(self._inAlnFile, referenceWindow,
                                   depthLimit=options.coverage,
                                   minMapQV=options.minMapQV,
                                   strategy="long-and-strand-balanced",
                                   stratum=options.readStratum,
                                   barcode=options.barcode)
        return (referenceWindow,
                pluralityConsensusAndVariants(referenceWindow, refSeqInWindow,
                                              alnHits, self.pluralityConfig))

# define both process and thread-based plurality callers
class PluralityWorkerProcess(PluralityWorker, WorkerProcess): pass
class PluralityWorkerThread(PluralityWorker, WorkerThread): pass

#
#  --------------------- Plugin API --------------------------------
#

# Pluggable module API for algorithms:
#  - Algorithm lives in a package
#  - Package must never fail to import, even if some of
#    its dependencies are not installed.
#  - Package must provide a main module exposing these top level
#    variables/methods:
#    - name                            = str
#    - availability                    = (bool, str)
#    - configure                      -> options -> cmph5 -> algorithm specific config object;
#                                        (can raise IncompatibleDataException)
#    - slaveFactories                 -> bool -> (class, class)

__all__ = [ "name",
            "availability",
            "configure",
            "slaveFactories" ]

name = "plurality"
availability = (True, "OK")

def slaveFactories(threaded):
    if threaded:
        return (PluralityWorkerThread,  ResultCollectorThread)
    else:
        return (PluralityWorkerProcess, ResultCollectorProcess)

def configure(options, cmpH5):
    pluralityConfig = PluralityConfig(minMapQV=options.minMapQV,
                                      minCoverage=options.minCoverage,
                                      minConfidence=options.minConfidence,
                                      diploid=options.diploid,
                                      noEvidenceConsensus=options.noEvidenceConsensusCall)
    return pluralityConfig
