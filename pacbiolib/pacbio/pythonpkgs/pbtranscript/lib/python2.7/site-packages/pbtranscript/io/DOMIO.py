"""Streaming IO support for DOM files."""

__all__ = ["DOMRecord",
           "DOMReader"]


from pbcore.io import ReaderBase
from pbcore.io._utils import splitFileContents


class DOMRecord(object):

    """A DOMRecord models a record in HMMER DOM file. """

    def __init__(self, pid, sid, score, pStart, pEnd, pLen,
                 sStart, sEnd, sLen):
        self.pid = pid  # primer id, ex: F1
        self.sid = sid  # ex: m1302...._s1_p0/45/1687_1987_back
        self.score = float(score)  # phmmer hit score
        self.pStart = int(pStart)  # hit start pos in primer
        self.pEnd = int(pEnd)      # hit end pos in primer
        self.pLen = int(pLen)      # primer length
        self.sStart = int(sStart)  # hit end pos in read
        self.sEnd = int(sEnd)      # hit start pos in read
        self.sLen = int(sLen)      # read length

    def __str__(self):
        return "{pid} {sid} {score} {ps} {pe} {pl} {ss} {se} {sl}".\
            format(score=self.score, pid=self.pid, sid=self.sid,
                   ps=self.pStart, pe=self.pEnd, pl=self.pLen,
                   ss=self.sStart, se=self.sEnd, sl=self.sLen)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return self.pid == other.pid and self.sid == other.sid and \
            self.score == other.score and self.pStart == other.pStart and \
            self.pEnd == other.pEnd and self.pLen == other.pEnd and \
            self.sStart == other.sStart and self.sEnd == other.sEnd and \
            self.sLen == other.sLen

    @classmethod
    def fromString(cls, line):
        """Construct and return a DOMRecord object given a DOM line."""
        try:
            fields = line.strip().split()
            assert len(fields) == 23
            return DOMRecord(pid=fields[0], sid=fields[3],
                             score=float(fields[13]), pStart=int(fields[17]) - 1,
                             pEnd=int(fields[18]), pLen=int(fields[2]),
                             sStart=int(fields[15]) - 1, sEnd=int(fields[16]),
                             sLen=int(fields[5]))
        except AssertionError:
            errMsg = "String not recognized as a valid DOM record."
            raise ValueError(errMsg)


class DOMReader(ReaderBase):

    """
    Streaming reader for DOM files.

    Example:

    .. doctest::
        >>> from pbtranscript.io import DOMReader
        >>> filename = "../../../tests/data/test_DOMReader.dom"
        >>> r = DOMReader(filename)
        >>> for record in r:
        ...     print record
        F1 m131018_081703_42161_c100585152550000001823088404281404_s1_p0/45/ccs 23.7 0 31 31 2170 2201 3931
        R1 m131018_081703_42161_c100585152550000001823088404281404_s1_p0/45/ccs 16.2 0 25 25 3906 3931 3931
    """

    def __iter__(self):
        try:
            lines = splitFileContents(self.file, "\n")
            for line in lines:
                line = line.strip()
                if len(line) > 0 and line[0] != "#":
                    yield DOMRecord.fromString(line)
        except AssertionError:
            raise ValueError("Invalid DOM file.")
