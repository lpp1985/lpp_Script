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
from .utils import CommonEqualityMixin

__all__ = [ "Variant" ]

class Variant(CommonEqualityMixin):
    """
    Variant objects represent homozygous/haploid OR heterozygous
    variants corresponding to a fixed window of a reference genome

    Internally we use Python-style half-open intervals zero-based
    [start, end) to delineate reference ranges.  An insertion has
     start==end, a SNP has start+1==end, etc.

    GFF files use 1-based indexing and open intervals [start, end).
    In a GFF both insertions and SNPs have start==end, which doesn't
    make too much sense to me, but so be it.

    VCF files use 1-based indexing as well, but do not record the
    "end"
    """
    def __init__(self, refId, refStart, refEnd, refSeq, readSeq1,
                 readSeq2=None, confidence=None, coverage=None,
                 frequency1=None, frequency2=None, annotations=None):
        self.refId       = refId
        self.refStart    = refStart
        self.refEnd      = refEnd
        self.refSeq      = refSeq
        self.readSeq1    = readSeq1
        self.readSeq2    = readSeq2
        self.confidence  = confidence
        self.coverage    = coverage
        self.frequency1  = frequency1
        self.frequency2  = frequency2
        self.annotations = annotations

    @property
    def isHeterozygous(self):
        return (self.readSeq2 != None)

    @property
    def variantType(self):
        lr = len(self.refSeq)
        l1 = len(self.readSeq1)
        l2 = len(self.readSeq2) if self.readSeq2 else None
        if lr == 0:
            return "Insertion"
        elif l1==0 or l2==0:
            return "Deletion"
        elif (l1==lr) and (l2==None or l2==lr):
            return "Substitution"
        else:
            return "Variant"

    def __str__(self):
        refSeq_ = self.refSeq or "."
        if self.isHeterozygous:
            readAlleles = "%s/%s" % (self.readSeq1 or ".",
                                     self.readSeq2 or ".")
        else:
            readAlleles = "%s" % (self.readSeq1 or ".")
        return "%s@%s:%d-%d %s -> %s" % \
            (self.variantType,
             self.refId,
             self.refStart,
             self.refEnd,
             refSeq_,
             readAlleles)

    def __repr__(self):
        return str(self)

    def __lt__(self, other):
        return ((self.refId, self.refStart, self.refEnd, self.readSeq1) <
                (other.refId, other.refStart, other.refEnd, other.readSeq1))

    def annotate(self, key, value):
        if self.annotations == None:
            self.annotations = []
        self.annotations.append((key, value))


def filterVariants(minCoverage, minConfidence, variants):
    return [ v for v in variants
             if ((v.coverage >= minCoverage) and
                 (v.confidence >= minConfidence)) ]

def annotateVariants(variants, alns):
    # Operates in place
    for v in variants:
        v.annotate("rows", ",".join(str(a.rowNumber) for a in alns))
