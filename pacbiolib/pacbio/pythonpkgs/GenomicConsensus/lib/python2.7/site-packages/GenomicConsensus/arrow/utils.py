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

import numpy as np, itertools, logging, re, sys
from collections import Counter
from traceback import format_exception

from GenomicConsensus.variants import *
from GenomicConsensus.utils import *
from GenomicConsensus.consensus import ArrowConsensus
from pbcore.io.rangeQueries import projectIntoRange
import ConsensusCore2 as cc

def uniqueSingleBaseMutations(templateSequence, positions=None):
    """
    Return an iterator over all single-base mutations of a
    templateSequence that result in unique mutated sequences.
    """
    allBases = "ACGT"
    positions = positions or xrange(0, len(templateSequence))
    for tplStart in positions:
        tplBase     = templateSequence[tplStart]
        prevTplBase = templateSequence[tplStart-1] if (tplStart > 0) else None
        # snvs
        for subsBase in allBases:
            if subsBase != tplBase:
                yield cc.Mutation(cc.MutationType_SUBSTITUTION,
                                  tplStart,
                                  subsBase)
        # Insertions---only allowing insertions that are not cognate
        # with the previous base.
        for insBase in allBases:
            if insBase != prevTplBase:
                yield cc.Mutation(cc.MutationType_INSERTION,
                                  tplStart,
                                  insBase)
        # Deletion--only allowed if refBase does not match previous tpl base
        if tplBase != prevTplBase:
            yield cc.Mutation(cc.MutationType_DELETION, tplStart)

def allSingleBaseMutations(templateSequence, positions=None):
    """
    Same as ``uniqueSingleBaseMutations``, but no filtering as to
    whether the mutated sequences are unique.
    """
    allBases = "ACGT"
    positions = positions or xrange(0, len(templateSequence))
    for tplStart in positions:
        tplBase = templateSequence[tplStart]
        # snvs
        for subsBase in allBases:
            if subsBase != tplBase:
                yield cc.Mutation(cc.MutationType_SUBSTITUTION,
                                  tplStart,
                                  subsBase)
        # Insertions
        for insBase in allBases:
            yield cc.Mutation(cc.MutationType_INSERTION,
                              tplStart,
                              insBase)
        # Deletion
        yield cc.Mutation(cc.MutationType_DELETION, tplStart)

def nearbyMutations(mutations, tpl, neighborhoodSize):
    """
    Return mutations nearby the previously-tried mutations
    """
    mutationPositions = map(cc.Mutation.Start, mutations)
    nearbyPositions = set()
    for mp in mutationPositions:
        nearbyPositions.update(range(max(0, mp - neighborhoodSize),
                                     min(len(tpl), mp + neighborhoodSize)))
    return uniqueSingleBaseMutations(tpl, sorted(nearbyPositions))

def bestSubset(mutationsAndScores, separation):
    """
    Given a list of (mutation, score) tuples, this utility method
    greedily chooses the highest scoring well-separated elements.  We
    use this to avoid applying adjacent high scoring mutations, which
    are the rule, not the exception.  We only apply the best scoring one
    in each neighborhood, and then revisit the neighborhoods after
    applying the mutations.
    """
    input = mutationsAndScores[:]
    output = []

    while input:
        best = max(input, key=snd)
        output.append(best)
        nStart = best[0].Start() - separation
        nEnd   = best[0].Start() + separation
        for t in input[:]:
            if nStart <= t[0].Start() <= nEnd:
                input.remove(t)

    return output

def refineConsensus(ai, arrowConfig):
    """
    Given a MultiReadMutationScorer, identify and apply favorable
    template mutations.  Return (consensus, didConverge) :: (str, bool)
    """
    cfg = cc.PolishConfig(arrowConfig.maxIterations,
                          arrowConfig.mutationSeparation,
                          arrowConfig.mutationNeighborhood)
    polishResult = cc.Polish(ai, cfg)
    return str(ai), polishResult.hasConverged

def consensusConfidence(ai, positions=None):
    """
    Returns an array of QV values reflecting the consensus confidence
    at each position specified.  If the `positions` argument is
    omitted, confidence values are returned for all positions in the
    consensus (str(ai)).
    """
    return np.array(np.clip(cc.ConsensusQualities(ai), 0, 93), dtype=np.uint8)

def variantsFromAlignment(a, refWindow, cssQvInWindow=None, siteCoverage=None, effectiveSiteCoverage=None):
    """
    Extract the variants implied by a pairwise alignment to the
    reference.
    """
    variants = []
    refId, refStart, _ = refWindow
    refPos = refStart
    cssPos = 0
    tbl = zip(a.Transcript(),
              a.Target(),
              a.Query())

    # We don't call variants where either the reference or css is 'N'
    grouper = lambda row: "N" if (row[1]=="N" or row[2]=="N") else row[0]
    runs = itertools.groupby(tbl, grouper)

    for code, run in runs:
        assert code in "RIDMN"
        run = list(run)
        ref = "".join(map(snd, run))
        refLen = len(ref) - Counter(ref)["-"]
        css = "".join(map(third, run))
        cssLen = len(css) - Counter(css)["-"]
        variant = None

        if code == "M" or code == "N":
            pass
        elif code == "R":
            assert len(css)==len(ref)
            variant = Variant(refId, refPos, refPos+len(css), ref, css)
        elif code == "I":
            variant = Variant(refId, refPos, refPos, "", css)
        elif code == "D":
            variant = Variant(refId, refPos, refPos + len(ref), ref, "")

        if variant is not None:
            # HACK ALERT: variants at the first and last position
            # are not handled correctly
            if siteCoverage is not None and np.size(siteCoverage) > 0:
                refPos_ = min(refPos-refStart, len(siteCoverage)-1)
                variant.coverage = siteCoverage[refPos_]
            if effectiveSiteCoverage is not None and np.size(effectiveSiteCoverage) > 0:
                refPos_ = min(refPos-refStart, len(siteCoverage)-1)
                variant.annotate("effectiveCoverage", effectiveSiteCoverage[refPos_])
                #import ipdb
                #ipdb.set_trace()
            if cssQvInWindow is not None and np.size(cssQvInWindow) > 0:
                cssPos_ = min(cssPos, len(cssQvInWindow)-1)
                variant.confidence = cssQvInWindow[cssPos_]
            variants.append(variant)

        refPos += refLen
        cssPos += cssLen

    return variants

def referenceSpanWithinWindow(referenceWindow, aln):
    """
    Helper function for sorting reads by their reference span
    after restriction to a window.
    """
    _, winStart, winEnd = referenceWindow
    return min(winEnd, aln.referenceEnd) - \
           max(winStart, aln.referenceStart)

def lifted(queryPositions, mappedRead):
    """
    Lift a mappedRead into a new coordinate system by using the
    position translation table `queryPositions`
    """
    newStart = queryPositions[mappedRead.TemplateStart]
    newEnd   = queryPositions[mappedRead.TemplateEnd]
    copy = cc.MappedRead(mappedRead)
    copy.TemplateStart = newStart
    copy.TemplateEnd = newEnd
    return copy


_typeMap = { cc.MutationType_INSERTION    : "Ins",
             cc.MutationType_DELETION     : "Del",
             cc.MutationType_SUBSTITUTION : "Sub" }

def _shortMutationDescription(mut, tpl):
    """
    More compact and uniform mutation description strings
    Examples:

    201 Ins . > G
    201 Sub C > T
    201 Del C > .
    """
    _type = _typeMap[mut.Type]
    _pos = mut.Start()
    _oldBase = "." if mut.Type == cc.MutationType_INSERTION \
               else tpl[_pos]
    _newBase = "." if mut.Type == cc.MutationType_DELETION \
               else mut.Base
    return "%d %s %s > %s" % (_pos, _type, _oldBase, _newBase)

def scoreMatrix(ai):
    """
    Returns (rowNames, columnNames, S)

    where:
      - S is a matrix where S_{ij} represents the score delta
        of mutation j against read i
      - rowNames[i] is an identifier name for the the read i---presently
        we use the the row number within the cmp.h5, encoded as a string
      - columnNames[j] is an identifier for mutation j, encoding the
        position, type, and base change
    """
    css = str(ai)
    allMutations = sorted(allSingleBaseMutations(css))
    readNames = list(ai.ReadNames())
    numReads = len(readNames)
    shape = (numReads, len(allMutations))
    scoreMatrix = np.zeros(shape)
    for j, mut in enumerate(allMutations):
        mutScores = ai.LLs(mut)
        scoreMatrix[:, j] = mutScores
    baselineScores =  np.array(ai.LLs())
    columnNames = [ _shortMutationDescription(mut, css)
                    for mut in allMutations ]
    rowNames = readNames
    return (rowNames, columnNames, baselineScores, scoreMatrix)


def variantsFromConsensus(refWindow, refSequenceInWindow, cssSequenceInWindow,
                          cssQvInWindow=None, siteCoverage=None, effectiveSiteCoverage=None,
                          aligner="affine", ai=None):
    """
    Compare the consensus and the reference in this window, returning
    a list of variants.
    """
    refId, refStart, refEnd = refWindow

    if aligner == "affine":
        align = cc.AlignAffine
    else:
        align = cc.Align

    ga = align(refSequenceInWindow, cssSequenceInWindow)

    return variantsFromAlignment(ga, refWindow, cssQvInWindow, siteCoverage, effectiveSiteCoverage)


def filterAlns(refWindow, alns, arrowConfig):
    """
    Given alns (already clipped to the window bounds), filter out any
    that are incompatible with Arrow.

    By and large we avoid doing any filtering to avoid potential
    reference bias in variant calling.

    However at the moment the POA (and potentially other components)
    breaks when there is a read of zero length.  So we throw away
    reads that are "stumpy", where the aligner has inserted a large
    gap, such that while the alignment technically spans the window,
    it may not have any read content therein:

          Ref   ATGATCCAGTTACTCCGATAAA
          Read  ATG---------------TA-A
          Win.     [              )
    """
    return [ a for a in alns
             if a.readLength >= (arrowConfig.readStumpinessThreshold * a.referenceSpan) and
                min(a.hqRegionSnr) >= arrowConfig.minHqRegionSnr and
                a.readScore >= arrowConfig.minReadScore ]


def sufficientlyAccurate(mappedRead, poaCss, minAccuracy):
    if minAccuracy <= 0.0:
        return True
    s, e = mappedRead.TemplateStart, mappedRead.TemplateEnd
    tpl = poaCss[s:e]
    if mappedRead.Strand == cc.StrandType_FORWARD:
        pass
    elif mappedRead.Strand == cc.StrandType_REVERSE:
        tpl = reverseComplement(tpl)
    else:
        return False
    aln = cc.AlignLinear(tpl, mappedRead.Seq)
    nErrors = sum(1 for t in aln.Transcript() if t != 'M')
    tlen = len(tpl)
    acc = 1.0 - 1.0 * min(nErrors, tlen) / tlen
    return acc >= minAccuracy


def consensusForAlignments(refWindow, refSequence, alns, arrowConfig, draft=None, polish=True,
                           alnsUsed=None):
    """
    Call consensus on this interval---without subdividing the interval
    further.

    Returns an ArrowConsensus object.

    Requires that clipping has already been done.

    If `draft` is provided, it will serve as the starting
    point for polishing.  If not, the POA will be used to generate a
    draft starting point.

    If `polish` is False, the arrow polishing procedure will not be
    used, and the draft consensus will be returned.

    `alnsUsed` is an output parameter; if not None, it should be an
    empty list on entry; on return from this function, the list will
    contain the alns objects that were actually used to compute the
    consensus (those not filtered out).
    """
    _, refStart, refEnd = refWindow

    if alnsUsed is not None:
        assert alnsUsed == []

    if draft is None:
        # Compute the POA consensus, which is our initial guess, and
        # should typically be > 99.5% accurate
        fwdSequences = [ a.read(orientation="genomic", aligned=False)
                         for a in alns
                         if a.spansReferenceRange(refStart, refEnd) ]
        assert len(fwdSequences) >= arrowConfig.minPoaCoverage

        try:
            p = cc.PoaConsensus.FindConsensus(fwdSequences[:arrowConfig.maxPoaCoverage])
        except:
            logging.info("%s: POA could not be generated" % (refWindow,))
            return ArrowConsensus.noCallConsensus(arrowConfig.noEvidenceConsensus,
                                                  refWindow, refSequence)
        draft = p.Sequence

    ga = cc.Align(refSequence, draft)

    # Extract reads into ConsensusCore2-compatible objects, and map them into the
    # coordinates relative to the POA consensus
    mappedReads = [ arrowConfig.extractMappedRead(aln, refStart) for aln in alns ]
    queryPositions = cc.TargetToQueryPositions(ga)
    mappedReads = [ lifted(queryPositions, mr) for mr in mappedReads ]

    # Load the mapped reads into the mutation scorer, and iterate
    # until convergence.
    ai = cc.MultiMolecularIntegrator(draft, cc.IntegratorConfig(arrowConfig.minZScore))
    coverage = 0
    for i, mr in enumerate(mappedReads):
        if (mr.TemplateEnd <= mr.TemplateStart or
            mr.TemplateEnd - mr.TemplateStart < 2 or
            mr.Length() < 2):
            continue
        if not sufficientlyAccurate(mr, draft, arrowConfig.minAccuracy):
            tpl = draft[mr.TemplateStart:mr.TemplateEnd]
            if mr.Strand == cc.StrandType_FORWARD:
                pass
            elif mr.Strand == cc.StrandType_REVERSE:
                tpl = reverseComplement(tpl)
            else:
                tpl = "INACTIVE/UNMAPPED"
            logging.debug("%s: skipping read '%s' due to insufficient accuracy, (poa, read): ('%s', '%s')" % (refWindow, mr.Name, tpl, mr.Seq))
            continue
        if ai.AddRead(mr) == cc.State_VALID:
            coverage += 1
            if alnsUsed is not None:
                alnsUsed.append(alns[i])

    if coverage < arrowConfig.minPoaCoverage:
        logging.info("%s: Inadequate coverage to call consensus" % (refWindow,))
        return ArrowConsensus.noCallConsensus(arrowConfig.noEvidenceConsensus,
                                              refWindow, refSequence)

    if not polish:
        confidence = np.zeros(len(draft), dtype=int)
        return ArrowConsensus(refWindow, draft, confidence, ai)

    # Iterate until covergence
    _, converged = refineConsensus(ai, arrowConfig)
    if converged:
        arrowCss = str(ai)
        if arrowConfig.computeConfidence:
            confidence = consensusConfidence(ai)
        else:
            confidence = np.zeros(shape=len(arrowCss), dtype=int)
        return ArrowConsensus(refWindow,
                              arrowCss,
                              confidence,
                              ai)
    else:
        logging.info("%s: Arrow did not converge to MLE" % (refWindow,))
        return ArrowConsensus.noCallConsensus(arrowConfig.noEvidenceConsensus,
                                              refWindow, refSequence)


def coverageInWindow(refWin, hits):
    winId, winStart, winEnd = refWin
    a = np.array([(hit.referenceStart, hit.referenceEnd)
                  for hit in hits
                  if hit.referenceName == winId])
    if len(a) == 0:
        return np.zeros(winEnd - winStart, dtype=np.uint)
    else:
        tStart = a[:,0]
        tEnd   = a[:,1]
        cov = projectIntoRange(tStart, tEnd, winStart, winEnd)
        return cov
