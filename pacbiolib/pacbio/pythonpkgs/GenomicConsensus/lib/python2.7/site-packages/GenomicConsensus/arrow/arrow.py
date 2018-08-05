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

# Authors: David Alexander, Lance Hepler

import logging, os.path
import ConsensusCore2 as cc, numpy as np

from .. import reference
from ..options import options
from ..Worker import WorkerProcess, WorkerThread
from ..ResultCollector import ResultCollectorProcess, ResultCollectorThread

from GenomicConsensus.consensus import Consensus, ArrowConsensus, join
from GenomicConsensus.windows import kSpannedIntervals, holes, subWindow
from GenomicConsensus.variants import filterVariants, annotateVariants
from GenomicConsensus.arrow.evidence import ArrowEvidence
from GenomicConsensus.arrow import diploid
from GenomicConsensus.utils import die

import GenomicConsensus.arrow.model as M
import GenomicConsensus.arrow.utils as U

def consensusAndVariantsForWindow(alnFile, refWindow, referenceContig,
                                  depthLimit, arrowConfig):
    """
    High-level routine for calling the consensus for a
    window of the genome given a cmp.h5.

    Identifies the coverage contours of the window in order to
    identify subintervals where a good consensus can be called.
    Creates the desired "no evidence consensus" where there is
    inadequate coverage.
    """
    winId, winStart, winEnd = refWindow
    logging.info("Arrow operating on %s" %
                 reference.windowToString(refWindow))

    if options.fancyChunking:
        # 1) identify the intervals with adequate coverage for arrow
        #    consensus; restrict to intervals of length > 10
        alnHits = U.readsInWindow(alnFile, refWindow,
                                  depthLimit=20000,
                                  minMapQV=arrowConfig.minMapQV,
                                  strategy="long-and-strand-balanced",
                                  stratum=options.readStratum,
                                  barcode=options.barcode)
        starts = np.fromiter((hit.tStart for hit in alnHits), np.int)
        ends   = np.fromiter((hit.tEnd   for hit in alnHits), np.int)
        intervals = kSpannedIntervals(refWindow, arrowConfig.minPoaCoverage,
                                      starts, ends, minLength=10)
        coverageGaps = holes(refWindow, intervals)
        allIntervals = sorted(intervals + coverageGaps)
        if len(allIntervals) > 1:
            logging.info("Usable coverage in %s: %r" %
                         (reference.windowToString(refWindow), intervals))

    else:
        allIntervals = [ (winStart, winEnd) ]

    # 2) pull out the reads we will use for each interval
    # 3) call consensusForAlignments on the interval
    subConsensi = []
    variants = []

    for interval in allIntervals:
        intStart, intEnd = interval
        intRefSeq = referenceContig[intStart:intEnd]
        subWin = subWindow(refWindow, interval)

        windowRefSeq = referenceContig[intStart:intEnd]
        alns = U.readsInWindow(alnFile, subWin,
                               depthLimit=depthLimit,
                               minMapQV=arrowConfig.minMapQV,
                               strategy="long-and-strand-balanced",
                               stratum=options.readStratum,
                               barcode=options.barcode)
        clippedAlns_ = [ aln.clippedTo(*interval) for aln in alns ]
        clippedAlns = U.filterAlns(subWin, clippedAlns_, arrowConfig)

        if len([ a for a in clippedAlns
                 if a.spansReferenceRange(*interval) ]) >= arrowConfig.minPoaCoverage:

            logging.debug("%s: Reads being used: %s" %
                          (reference.windowToString(subWin),
                           " ".join([str(hit.readName) for hit in alns])))

            alnsUsed = [] if options.reportEffectiveCoverage else None
            css = U.consensusForAlignments(subWin,
                                           intRefSeq,
                                           clippedAlns,
                                           arrowConfig,
                                           alnsUsed=alnsUsed)

            # Tabulate the coverage implied by these alignments, as
            # well as the post-filtering ("effective") coverage
            siteCoverage = U.coverageInWindow(subWin, alns)
            effectiveSiteCoverage = U.coverageInWindow(subWin, alnsUsed) if options.reportEffectiveCoverage else None

            variants_ = U.variantsFromConsensus(subWin, windowRefSeq,
                                                css.sequence, css.confidence, siteCoverage, effectiveSiteCoverage,
                                                options.aligner,
                                                ai=None)

            filteredVars =  filterVariants(options.minCoverage,
                                           options.minConfidence,
                                           variants_)
            # Annotate?
            if options.annotateGFF:
                annotateVariants(filteredVars, clippedAlns)

            variants += filteredVars

            # Dump?
            maybeDumpEvidence = \
                ((options.dumpEvidence == "all") or
                 (options.dumpEvidence == "outliers") or
                 (options.dumpEvidence == "variants") and (len(variants) > 0))
            if maybeDumpEvidence:
                refId, refStart, refEnd = subWin
                refName = reference.idToName(refId)
                windowDirectory = os.path.join(
                    options.evidenceDirectory,
                    refName,
                    "%d-%d" % (refStart, refEnd))
                ev = ArrowEvidence.fromConsensus(css)
                if options.dumpEvidence != "outliers":
                    ev.save(windowDirectory)
                elif (np.max(ev.delta) > 20):
                    # Mathematically I don't think we should be seeing
                    # deltas > 6 in magnitude, but let's just restrict
                    # attention to truly bonkers outliers.
                    ev.save(windowDirectory)

        else:
            css = ArrowConsensus.noCallConsensus(arrowConfig.noEvidenceConsensus,
                                                 subWin, intRefSeq)
        subConsensi.append(css)

    # 4) glue the subwindow consensus objects together to form the
    #    full window consensus
    css = join(subConsensi)

    # 5) Return
    return css, variants


class ArrowWorker(object):

    @property
    def arrowConfig(self):
        return self._algorithmConfig

    def onChunk(self, workChunk):
        referenceWindow  = workChunk.window
        refId, refStart, refEnd = referenceWindow

        refSeqInWindow = reference.sequenceInWindow(referenceWindow)

        # Quick cutout for no-coverage case
        if not workChunk.hasCoverage:
            noCallCss = ArrowConsensus.noCallConsensus(self.arrowConfig.noEvidenceConsensus,
                                                       referenceWindow, refSeqInWindow)
            return (referenceWindow, (noCallCss, []))

        # General case
        eWindow = reference.enlargedReferenceWindow(referenceWindow,
                                                    options.referenceChunkOverlap)
        _, eStart, eEnd = eWindow

        # We call consensus on the enlarged window and then map back
        # to the reference and clip the consensus at the implied
        # bounds.  This seems to be more reliable thank cutting the
        # consensus bluntly
        refContig = reference.byName[refId].sequence
        refSequenceInEnlargedWindow = refContig[eStart:eEnd]

        #
        # Get the consensus for the enlarged window.
        #
        css_, variants_ = \
            consensusAndVariantsForWindow(self._inAlnFile, eWindow,
                                          refContig, options.coverage, self.arrowConfig)

        #
        # Restrict the consensus and variants to the reference window.
        #
        ga = cc.Align(refSequenceInEnlargedWindow, css_.sequence)
        targetPositions = cc.TargetToQueryPositions(ga)
        cssStart = targetPositions[refStart-eStart]
        cssEnd   = targetPositions[refEnd-eStart]

        cssSequence    = css_.sequence[cssStart:cssEnd]
        cssQv          = css_.confidence[cssStart:cssEnd]
        variants       = [ v for v in variants_
                           if refStart <= v.refStart < refEnd ]

        consensusObj = Consensus(referenceWindow,
                                 cssSequence,
                                 cssQv)

        return (referenceWindow, (consensusObj, variants))



#
# Slave process/thread classes
#
class ArrowWorkerProcess(ArrowWorker, WorkerProcess): pass
class ArrowWorkerThread(ArrowWorker, WorkerThread): pass


#
# Plugin API
#
__all__ = [ "name",
            "availability",
            "configure",
            "slaveFactories" ]

name = "arrow"
availability = (True, "OK")

def configure(options, alnFile):
    if alnFile.readType != "standard":
        raise U.IncompatibleDataException(
            "The Arrow algorithm requires a BAM file containing standard (non-CCS) reads." )

    if options.diploid:
        logging.warn("Diploid analysis not yet supported under Arrow model.")

    # load parameters from file
    if options.parametersFile:
        logging.info("Loading model parameters from: ({0})".format(options.parametersFile))
        if not cc.LoadModels(options.parametersFile):
            die("Arrow: unable to load parameters from: ({0})".format(options.parametersFile))

    # test available chemistries
    supp = set(cc.SupportedChemistries())
    logging.info("Found consensus models for: ({0})".format(", ".join(sorted(supp))))

    used = set([])
    if options.parametersSpec != "auto":
        logging.info("Overriding model selection with: ({0})".format(options.parametersSpec))
        if not cc.OverrideModel(options.parametersSpec):
            die("Arrow: unable to override model with: ({0})".format(options.parametersSpec))
        used.add(options.parametersSpec)
    else:
        used.update(alnFile.sequencingChemistry)
        unsupp = used - supp
        if used - supp:
            die("Arrow: unsupported chemistries found: ({0})".format(", ".join(sorted(unsupp))))

    # All arrow models require PW except P6 and the first S/P1-C1
    for readGroup in alnFile.readGroupTable:
        if set([readGroup["SequencingChemistry"]]) - set(["P6-C4", "S/P1-C1/beta"]):
            if ("Ipd" not in readGroup["BaseFeatures"] or
                "PulseWidth" not in readGroup["BaseFeatures"]):
                die("Arrow model requires missing base feature: IPD or PulseWidth")

    logging.info("Using consensus models for: ({0})".format(", ".join(sorted(used))))

    return M.ArrowConfig(minMapQV=options.minMapQV,
                         noEvidenceConsensus=options.noEvidenceConsensusCall,
                         computeConfidence=(not options.fastMode),
                         minReadScore=options.minReadScore,
                         minHqRegionSnr=options.minHqRegionSnr,
                         minZScore=options.minZScore,
                         minAccuracy=options.minAccuracy)

def slaveFactories(threaded):
    # By default we use slave processes. The tuple ordering is important.
    if threaded:
        return (ArrowWorkerThread,  ResultCollectorThread)
    else:
        return (ArrowWorkerProcess, ResultCollectorProcess)
