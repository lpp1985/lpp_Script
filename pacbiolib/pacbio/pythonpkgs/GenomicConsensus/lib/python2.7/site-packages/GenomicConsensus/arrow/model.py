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

import numpy as np, ConfigParser, collections, logging
from glob import glob
from os.path import join
from pkg_resources import resource_filename, Requirement

from GenomicConsensus.utils import die
from GenomicConsensus.arrow.utils import fst, snd
from pbcore.chemistry import ChemistryLookupError
from pbcore.io import CmpH5Alignment
import ConsensusCore2 as cc

__all__ = [ "ArrowConfig" ]


#
#  ArrowConfig: the kitchen sink class of arrow options
#

class ArrowConfig(object):
    """
    Arrow configuration options
    """
    def __init__(self,
                 minMapQV=10,
                 minPoaCoverage=3,
                 maxPoaCoverage=11,
                 mutationSeparation=10,
                 mutationNeighborhood=20,
                 maxIterations=40,
                 noEvidenceConsensus="nocall",
                 computeConfidence=True,
                 readStumpinessThreshold=0.1,
                 minReadScore=0.75,
                 minHqRegionSnr=3.75,
                 minZScore=-3.5,
                 minAccuracy=0.82):

        self.minMapQV                   = minMapQV
        self.minPoaCoverage             = minPoaCoverage
        self.maxPoaCoverage             = maxPoaCoverage
        self.mutationSeparation         = mutationSeparation
        self.mutationNeighborhood       = mutationNeighborhood
        self.maxIterations              = maxIterations
        self.noEvidenceConsensus        = noEvidenceConsensus
        self.computeConfidence          = computeConfidence
        self.readStumpinessThreshold    = readStumpinessThreshold
        self.minReadScore               = minReadScore
        self.minHqRegionSnr             = minHqRegionSnr
        self.minZScore                  = minZScore
        self.minAccuracy                = minAccuracy

    def extractMappedRead(self, aln, windowStart):
        """
        Given a clipped alignment, convert its coordinates into template
        space (starts with 0), bundle it up with its features as a
        MappedRead.
        """
        if isinstance(aln, CmpH5Alignment):
            die("Arrow does not support CmpH5 files!")

        assert aln.referenceSpan > 0

        def baseFeature(featureName):
            if aln.reader.hasBaseFeature(featureName):
                rawFeature = aln.baseFeature(featureName, aligned=False, orientation="native")
                return rawFeature.clip(0,255).astype(np.uint8)
            else:
                return np.zeros((aln.readLength,), dtype=np.uint8)

        name = aln.readName
        chemistry = aln.sequencingChemistry
        strand = cc.StrandType_REVERSE if aln.isReverseStrand else cc.StrandType_FORWARD
        read = cc.Read(name,
                       aln.read(aligned=False, orientation="native"),
                       cc.Uint8Vector(baseFeature("Ipd").tolist()),
                       cc.Uint8Vector(baseFeature("PulseWidth").tolist()),
                       cc.SNR(aln.hqRegionSnr),
                       chemistry)
        return cc.MappedRead(read,
                             strand,
                             int(aln.referenceStart - windowStart),
                             int(aln.referenceEnd   - windowStart))
