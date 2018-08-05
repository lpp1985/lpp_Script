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

import time
from pbcore.io import GffWriter, Gff3Record
from GenomicConsensus import __VERSION__, reference


def gffVariantSeq(var):
    if var.isHeterozygous:
        return "%s/%s" % (var.readSeq1 or ".",
                          var.readSeq2 or ".")
    else:
        return var.readSeq1 or "."

def gffVariantFrequency(var):
    if var.frequency1==None:
        return None
    elif var.isHeterozygous:
        return "%d/%d" % (var.frequency1, var.frequency2)
    else:
        return str(var.frequency1)

def toGffRecord(var):
    varType  = var.variantType
    gffType  = varType.lower()
    gffStart = (var.refStart + 1) if (var.refSeq != "") else var.refStart
    gffEnd   = var.refEnd         if (var.refSeq != "") else var.refStart
    gffFreq = gffVariantFrequency(var)

    record = Gff3Record(reference.idToFullName(var.refId), gffStart, gffEnd, gffType)
    record.reference  = var.refSeq or "."
    record.variantSeq = gffVariantSeq(var)
    if gffFreq:
        record.frequency  = gffFreq
    record.coverage   = var.coverage
    record.confidence = var.confidence
    if var.annotations:
        for (k, v) in var.annotations:
            record.put(k, v)
    return record

class VariantsGffWriter(object):

    ONTOLOGY_URL = \
        "http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12"

    def __init__(self, f, optionsDict, referenceEntries):
        self._gffWriter = GffWriter(f)
        self._gffWriter.writeHeader("##pacbio-variant-version 2.1")
        self._gffWriter.writeHeader("##date %s" % time.ctime())
        self._gffWriter.writeHeader("##feature-ontology %s" % self.ONTOLOGY_URL)
        self._gffWriter.writeHeader("##source GenomicConsensus %s" % __VERSION__)
        self._gffWriter.writeHeader("##source-commandline %s" % optionsDict["shellCommand"])
        self._gffWriter.writeHeader("##source-alignment-file %s" % optionsDict["inputFilename"])
        self._gffWriter.writeHeader("##source-reference-file %s" % optionsDict["referenceFilename"])
        # Reference groups.
        for entry in referenceEntries:
            self._gffWriter.writeHeader("##sequence-region %s 1 %d" \
                                            % (entry.name, entry.length))

    def writeVariants(self, variants):
        for var in variants:
            self._gffWriter.writeRecord(toGffRecord(var))

    def close(self):
        self._gffWriter.close()
