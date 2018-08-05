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

from GenomicConsensus.arrow.utils import allSingleBaseMutations
from GenomicConsensus.variants import Variant

import numpy as np
import ConsensusCore2 as cc

# IUPAC reference:
#   http://www.bioinformatics.org/sms/iupac.html

_packIupac =   { ("A", "G") : "R" ,
                 ("G", "A") : "R" ,
                 ("C", "T") : "Y" ,
                 ("T", "C") : "Y" ,
                 ("G", "C") : "S" ,
                 ("C", "G") : "S" ,
                 ("A", "T") : "W" ,
                 ("T", "A") : "W" ,
                 ("G", "T") : "K" ,
                 ("T", "G") : "K" ,
                 ("A", "C") : "M" ,
                 ("C", "A") : "M" }

_unpackIupac = { "R" : ("A", "G") ,
                 "Y" : ("C", "T") ,
                 "S" : ("G", "C") ,
                 "W" : ("A", "T") ,
                 "K" : ("G", "T") ,
                 "M" : ("A", "C") }

def packIUPAC(bases):
    return _packIupac[bases]

def unpackIUPAC(iupacCode):
    return _unpackIupac[iupacCode]

def isHeterozygote(base):
    return (base in _unpackIupac)

def packMuts(cssBase, mut1, mut2):
    # Turn two muts (with same Start, End, LengthDiff) into a single mutation to
    # IUPAC.  The no-op mutation is coded as None.
    #
    # Example1: (_, Subs A, Subs T) -> Subs W
    # Example2: (_, Ins A, Ins T)   -> Ins W
    # Example3: (A, None, Subs T)   -> Subs W
    #
    nonNullMut = mut1 or mut2
    start   = nonNullMut.Start()
    mutType = nonNullMut.Type
    newBase1 = mut1.Base if mut1 else cssBase
    newBase2 = mut2.Base if mut2 else cssBase
    newBasePacked = packIUPAC((newBase1, newBase2))
    return cc.Mutation(mutType, start, newBasePacked)


def scoresForPosition(ai, pos):
    muts = allSingleBaseMutations(str(ai), positions=[pos])
    noMutScore = [0] * ai.NumReads()
    mutScores_ = [ ai.ReadLLs(mut)
                   for mut in muts ]
    mutScores = np.column_stack([noMutScore] + mutScores_).astype(np.float32)
    return mutScores


def variantsFromConsensus(refWindow, refSequenceInWindow, cssSequenceInWindow,
                          cssQvInWindow=None, siteCoverage=None, aligner="affine",
                          ai=None):
    """
    Compare the consensus and the reference in this window, returning
    a list of variants.

    Uses the integrator to identify heterozygous variants.
    """
    assert (cssQvInWindow is None) == (siteCoverage is None)  # Both or none

    refId, refStart, refEnd = refWindow

    if ai is not None:
        #
        # Hunting diploid variants:
        # 1. find confident heterozygous sites;
        # 2. build a "diploid consensus" using IUPAC encoding
        #    for het sites; mark cssQv accordingly
        # 3. align diploid consensus to reference
        # 4. extract and decorate variants
        #
        assert str(ai) == cssSequenceInWindow
        iupacMutations = []  # List of (Mutation, confidence)
        for pos in xrange(0, ai.Length()):
            ds = cc.IsSiteHeterozygous(scoresForPosition(ai, pos), 40)
            if ds:
                muts = [None] + list(allSingleBaseMutations(cssSequenceInWindow, positions=[pos]))
                mut0 = muts[ds.Allele0]
                mut1 = muts[ds.Allele1]
                cssBase = cssSequenceInWindow[pos]
                packedMut = packMuts(cssBase, mut0, mut1)
                iupacMutations.append((packedMut, 40))

        # Create diploidCss by applying mutations, meanwhile updating the
        # confidence vector accordingly.
        diploidCss = cc.ApplyMutations([pair[0] for pair in iupacMutations],
                                       cssSequenceInWindow)

        diploidQv  = list(cssQvInWindow) if cssQvInWindow is not None else None

        runningLengthDiff = 0
        for (mut, conf) in iupacMutations:
            start = mut.Start() + runningLengthDiff
            end   = mut.End() + runningLengthDiff
            diploidQv[start:end] = [conf]
        assert len(diploidCss) == len(diploidQv)

        cssSequenceInWindow = diploidCss
        cssQvInWindow = diploidQv

    vars = variantsFromAlignment(refWindow,
                                 refSequenceInWindow, cssSequenceInWindow,
                                 cssQvInWindow, siteCoverage)
    return vars


def variantsFromAlignment(refWindow, refSeq, cssSeq,
                          cssQV=None, refCoverage=None):
    """
    Extract the variants implied by a pairwise alignment of cssSeq to
    refSeq reference.  If cssQV, refCoverage are provided, they will
    be used to decorate the variants with those attributes.

    Arguments:
      - cssQV: QV array, same length as css
      - refCoverage: coverage array, sample length as reference window

    This is trickier than in the haploid case.  We have to break out
    diploid variants as single bases, in order to avoid implying
    phase.
    """
    variants = []
    refId, refStart, refEnd = refWindow

    aln = cc.AlignAffineIupac(refSeq, cssSeq);
    alnTarget = aln.Target()
    alnQuery = aln.Query()

    assert (cssQV is None) == (refCoverage is None)  # Both or none
    assert len(refSeq) == refEnd - refStart
    assert cssQV is None or len(cssSeq) == len(cssQV)
    assert refCoverage is None or len(refSeq) == len(refCoverage)

    transcript = [ X if (Q != "N" and T != "N") else "N"
                   for (X, T, Q) in zip(aln.Transcript(),
                                        alnTarget,
                                        alnQuery) ]
    variants = []
    runStart = -1
    runStartRefPos = None
    runX = None
    refPos = refStart
    for pos, (X, T, Q) in enumerate(zip(transcript,
                                        alnTarget,
                                        alnQuery)):
        if X != runX or isHeterozygote(Q):
            if runStart >= 0 and runX not in "MN":
                # Package up the run and dump a variant
                ref  = alnTarget[runStart:pos].replace("-", "")
                read = alnQuery [runStart:pos].replace("-", "")
                if isHeterozygote(read):
                    allele1, allele2 = unpackIUPAC(read)
                    var = Variant(refId, runStartRefPos, refPos, ref, allele1, allele2)
                else:
                    var = Variant(refId, runStartRefPos, refPos, ref, read)
                variants.append(var)
            runStart = pos
            runStartRefPos = refPos
            runX = X
        if T != "-": refPos += 1


    # This might be better handled within the loop above, just keeping
    # track of Qpos, Tpos
    if cssQV is not None:
        cssPosition = cc.TargetToQueryPositions(aln)
        for v in variants:
            # HACK ALERT: we are not really handling the confidence or
            # coverage for variants at last position of the window
            # correctly here.
            refPos_ = min(v.refStart-refStart, len(refCoverage)-1)
            cssPos_ = min(cssPosition[v.refStart-refStart], len(cssQV)-1)

            if refCoverage is not None: v.coverage   = refCoverage[refPos_]
            if cssQV       is not None: v.confidence = cssQV[cssPos_]

    return variants
