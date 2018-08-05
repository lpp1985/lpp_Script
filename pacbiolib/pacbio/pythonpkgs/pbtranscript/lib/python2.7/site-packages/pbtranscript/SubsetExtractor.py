"""Sub set input Fasta file."""
import logging
from pbcore.io import FastaReader, FastaWriter
from pbtranscript.io import ReadAnnotation
from collections import namedtuple

SubsetRules = namedtuple("SubsetRules", "FL nonChimeric")


class ReadsSubsetExtractor(object):

    """An object of this class can subset reads in an input Fasta file."""

    def __init__(self, inFN, outFN, rules, ignore_polyA, printReadLengthOnly=False):
        self.inFN = inFN
        self.outFN = outFN
        self.rules = rules
        self.printReadLengthOnly = printReadLengthOnly
        self.ignore_polyA = ignore_polyA

    def satisfy(self, annotation, rules):
        """Return whether or not this annotated read satisfy the
        subset rules.
        FL =
            0, reads must not be full length;
            1, reads must be full length;
            2, reads can either be full length or not full length.
        nonChimeric =
            0, reads must be chimeras;
            1, reads must be non-chimeric;
            2, reads can either be chimeric or non-chimeric.
        """
        if (rules.FL == 1 and not annotation.isFullLength) or \
           (rules.FL == 0 and annotation.isFullLength):
            return False
        if (rules.nonChimeric == 1 and annotation.chimera) or \
           (rules.nonChimeric == 0 and not annotation.chimera):
            return False
        return True

    def run(self):
        """Subset reads based on read annotation and subset rules."""
        infoMsg = "Extracting reads from {f} based on ".format(f=self.inFN)
        infoMsg += "rules(FullLength={fl}, nonChimeric={nc}).".format(
                   fl="true" if self.rules.FL != 0 else "false",
                   nc="true" if self.rules.nonChimeric != 0 else "false")
        logging.info(infoMsg)

        if not self.printReadLengthOnly:
            with FastaReader(self.inFN) as reader, \
                    FastaWriter(self.outFN) as writer:
                for r in reader:
                    annotation = ReadAnnotation.fromString(r.name,
                                                           self.ignore_polyA)
                    if self.satisfy(annotation, self.rules):
                        writer.writeRecord(r.name, r.sequence[:])
        else:  # print read length only, dont print read names and sequences
            with FastaReader(self.inFN) as reader, \
                    open(self.outFN, 'w') as writer:
                for r in reader:
                    annotation = ReadAnnotation.fromString(r.name,
                                                           self.ignore_polyA)
                    if self.satisfy(annotation, self.rules):
                        writer.write("{rl}\n".format(rl=len(r.sequence)))
