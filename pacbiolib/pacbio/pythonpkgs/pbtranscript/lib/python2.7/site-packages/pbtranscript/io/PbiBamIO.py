# pylint: disable=abstract-method
"""
Defines classes for fetching data from a collection of original SMRTCell bam
files, and for writing bam files.
"""

import os.path as op
import numpy as np
import pysam

from pbcore.io import (SubreadSet, ConsensusReadSet,
                       ContigSet, FastaReader, openDataFile, openDataSet, IndexedBamReader)
from pbcore.io import BamAlignment

from pbtranscript.Utils import get_files_from_file_or_fofn


__all__ = ["BamCollection",
           "BamHeader",
           "BamWriter"]


class BamZmw(object):

    """Represent zmws in bam file, like Zmw in pbcore.io"""

    def __init__(self, bamRecords, isCCS):
        if isinstance(bamRecords, list):
            assert len(bamRecords) > 0
            assert isinstance(bamRecords[0], BamAlignment)
        else:
            assert isinstance(bamRecords, BamAlignment)
            bamRecords = [bamRecords]
        assert len(bamRecords) >= 1
        self.bamRecords = bamRecords
        self._isCCS = isCCS
        if self._isCCS:
            assert len(bamRecords) == 1
        if not all([b.holeNumber == self.holeNumber for b in bamRecords]):
            raise ValueError("%s's input reads must come from same zmw." %
                             self.__class__.__name__)

    @property
    def holeNumber(self):
        """Return hole numer of this zmw."""
        return self.bamRecords[0].holeNumber

    @property
    def movieName(self):
        """Return movie name of this Zmw."""
        return self.bamRecords[0].qName.split("/")[0]

    @property
    def isCCS(self):
        """True if this zmw has ccs reads; False otherwise."""
        return self._isCCS

    @property
    def zmwName(self):
        """Return zmw name in string."""
        return "/".join(self.bamRecords[0].qName.split("/")[0:2])

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__,
                             self.zmwName)

    @property
    def ccsRead(self):
        """Return ccs read of this zmw. None if this zmw only
        contains subreads."""
        if not self.isCCS:
            return None
        assert len(self.bamRecords) == 1
        return BamCCSZmwRead(self, self.bamRecords[0], self.holeNumber)

    @property
    def subreads(self):
        """Return a list of subreads of this zmw."""
        if self.isCCS:
            return []
        return [BamZmwRead(self, _record, self.holeNumber)
                for _record in self.bamRecords]

    def read(self, readStart=None, readEnd=None):
        """Return read of movie/zmw/readStart_readEnd.
           If either readStart or readEnd are None, return read
           of the zmw.
        """
        _record = None
        if readStart is None or readEnd is None:
            if len(self.bamRecords) != 1:
                # Reading full unrolled sequence in BamZmw is not supported
                raise ValueError("Must specify read start and end when" +
                                 "a BamZmw contains multiple bam records.")
            else:
                _record = self.bamRecords[0]
                readStart, readEnd = 0, _record.qLen
                if not self.isCCS:
                    readStart, readEnd = _record.qStart, _record.qEnd
        else:
            _record = _findBamRecord(self.bamRecords, readStart, readEnd)
            if _record is None:
                raise IndexError("Invalid slice of zmw")
        return BamZmwRead(self, _record, self.holeNumber, readStart, readEnd)


def _findBamRecord(bamRecords, readStart, readEnd):
    """
    Given a list of BamAlignment object each of which represents a read,
    return the first record which covers interval [readStart, readEnd).
    """
    for bamRecord in bamRecords:
        if bamRecord.qStart <= readStart and bamRecord.qEnd >= readEnd:
            return bamRecord
    return None


class BamZmwRead(object):

    """Represent zmw read in bam files, like ZmwRead in pbcore.io"""
    __slots__ = ["zmw", "bamRecord", "holeNumber", "readStart", "readEnd"]

    def __init__(self, zmw, bamRecord, holeNumber,
                 readStart=None, readEnd=None):
        assert isinstance(zmw, BamZmw)
        self.zmw = zmw
        self.bamRecord = bamRecord
        assert isinstance(self.bamRecord, BamAlignment)
        assert holeNumber == self.zmw.holeNumber
        if readStart is not None:
            self.readStart = readStart
        else:
            self.readStart = self.qStart
        if readEnd is not None:
            self.readEnd = readEnd
        else:
            self.readEnd = self.qEnd
        if self.readStart < self.qStart or \
           self.readEnd > self.qEnd:
            raise ValueError("[%d, %d) not in %s" % (self.readStart,
                                                     self.readEnd,
                                                     self.bamRecord))

    @property
    def qStart(self):
        """
        Subread start pos in coordinate of zmw polymerase read, inclusive.
        """
        if self.bamRecord.isCCS:
            return 0
        else:
            return self.bamRecord.qStart

    @property
    def qEnd(self):
        """
        Subread end pos in coordinate of zmw polymerase read, exclusive.
        """
        if self.bamRecord.isCCS:
            return self.bamRecord.qLen
        else:
            return self.bamRecord.qEnd

    @property
    def movieName(self):
        """Return movie name of this read."""
        return self.zmw.movieName

    @property
    def holeNumber(self):
        """Return hole number of this read."""
        return self.zmw.holeNumber

    @property
    def readName(self):
        """Return zmw name string."""
        return "%s/%d_%d" % (self.zmw.zmwName,
                             self.readStart,
                             self.readEnd)

    @property
    def alignedSegment(self):
        """Return a pysam.calignmentfile.AlignedSegment
        object which is associated with this record.
        """
        return self.bamRecord.peer

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__,
                             self.readName)

    def __len__(self):
        return self.readEnd - self.readStart

    @property
    def _offsetBegin(self):
        """Read start offset."""
        return self.readStart - self.qStart

    @property
    def _offsetEnd(self):
        """Read end offset."""
        return self.readEnd - self.qStart

    def basecalls(self):
        """Return base calls of this zmw read in string."""
        return self.bamRecord.read(False)[self._offsetBegin: self._offsetEnd]

    def qv(self, _qvName):
        """Return QVs of this zmw read."""
        qvName = _qvName.lower()
        qvs = []
        try:
            if qvName == "deletionqv":
                qvs = self.bamRecord.DeletionQV(False)
            elif qvName == "insertionqv":
                qvs = self.bamRecord.InsertionQV(False)
            elif qvName == "mergeqv":
                qvs = self.bamRecord.MergeQV(False)
            elif qvName == "substitutionqv":
                qvs = self.bamRecord.SubstitutionQV(False)
            elif qvName == "substitutiontag":
                qvs = self.bamRecord.SubstitutionTag(False)
            elif qvName == "deletiontag":
                qvs = self.bamRecord.DeletionTag(False)
            elif qvName == "qualityvalue":
                if self.bamRecord.peer.qual is None:
                    raise AttributeError("No Quality Value")
                qvs = np.fromstring(self.bamRecord.peer.qual,
                                    dtype=np.uint8) - 33
            elif qvName == "ipd": # although not used by isoseq
                qvs = self.bamRecord.IPD(False)
            elif qvName == "pulsewidth": # although not used by isoseq
                qvs = self.bamRecord.PulseWidth(False)
            return qvs[self._offsetBegin: self._offsetEnd]
        except AttributeError:
            raise AttributeError("%s does not contain QV %s." %
                                 (self.readName, _qvName))

    def Clip(self, readStart, readEnd):
        """
        Clip this read to [readStart:readEnd), and return a copy of
        pysam.calignmentfile.AlignedSegment object.
        Assume that read.bamRecord.peer is an unmapped AlignedSegment.
        """
        new_query_name = "%s/%s/%d_%d" % (self.movieName, self.holeNumber,
                                          readStart, readEnd)
        if not (readStart >= self.readStart and readStart <= readEnd and
                readEnd <= self.readEnd):
            raise ValueError("Unable to clip subread %s from read %s." %
                             (new_query_name, self.readName))

        s, e = readStart - self.readStart, readEnd - self.readStart
        QV_TAGS = ["iq", "dq", "dt", "st", "sq", "mq", "ip", "pw"]

        # Create an unaligned pysam.AlignedSegment object.
        ret = pysam.AlignedSegment()
        ret.query_name = new_query_name

        peer = self.bamRecord.peer
        ret.query_sequence = peer.query_sequence[s:e]
        ret.flag = peer.flag

        assert peer.reference_id == -1
        assert peer.reference_start == -1
        assert peer.cigartuples is None
        assert peer.next_reference_id == -1
        assert peer.next_reference_start == -1
        assert peer.template_length == 0

        ret.reference_id = peer.reference_id
        ret.reference_start = peer.reference_start
        ret.cigar = []
        ret.next_reference_id = peer.next_reference_id
        ret.next_reference_start = peer.next_reference_start
        ret.template_length = peer.template_length

        ret.mapping_quality = peer.mapping_quality
        if peer.query_qualities is None:
            ret.query_qualities = None
        else:
            ret.query_qualities = pysam.array_to_qualitystring(peer.query_qualities)[s:e]

        tags = peer.tags[::]
        for index, (tag_name, tag_val) in enumerate(tags):
            if tag_name in QV_TAGS:
                if self.__len__() != len(tag_val):
                    raise ValueError("%s's %s length %d ! = sequence length %d" %
                                     (peer.query_name, tag_name, len(tag_val), self.__len__()))
                tags[index] = (tag_name, tag_val[s:e])
            elif tag_name == 'qs':
                tags[index] = (tag_name, int(readStart))
            elif tag_name == 'qe':
                tags[index] = (tag_name, int(readEnd))

        ret.tags = tags
        return ret

    def __delitem__(self, dummy_name):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__delitem__"))

    def __setitem__(self, dummy_index, dummy_name):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__setitem__"))

    def __getitem__(self, key):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__getitem__"))


class BamCCSZmwRead(BamZmwRead):

    """Represent Zmw ccs read."""
    @property
    def readName(self):
        return "%s/ccs" % self.zmw.zmwName


class BamCollection(object):
    """
    Class wraps PacBio bam fofn or datasets.

    Can be initialized from a list of bam files, or an input.fofn
    file containing a list of bam files, or a subreads dataset or a
    ccs dataset.

    Note that ALL bam files MUST have pbi index generated.
    """

    def __init__(self, *args):
        if len(args) == 1:
            args = get_files_from_file_or_fofn(args[0])
        self._dataset = openDataFile(*args)
        # Implementation notes: find all the bam files, and group
        # them together by movieName
        self._header = BamHeader(ignore_pg=True)

        for bam in self._dataset.resourceReaders():
            if not isinstance(bam, IndexedBamReader):
                raise ValueError("%s in %s must have pbi index generated",
                                 bam.filename, str(self._dataset))
            self._header.add(bam.peer.header)
            for rg in bam.peer.header["RG"]: #readGroupTable:
                if rg['PL'] != "PACBIO":
                    raise IOError("Input BAM file %s for %s must be PacBio BAM.",
                                  bam.filename, self.__class__.__name__)
            for rg in bam.readGroupTable:
                assert rg.ReadType in ["CCS", "SUBREAD"]

    @property
    def movieNames(self):
        """Return movie names as a list of string."""
        return self._dataset.movieIds.keys()

    @property
    def header(self):
        """Return a BamHeader object which combines headers from
        all readers.
        """
        return self._header

    @property
    def isCCS(self):
        """Is CCS read set?"""
        return isinstance(self._dataset, ConsensusReadSet)

    @property
    def isSubread(self):
        """Is Subread read set?"""
        return isinstance(self._dataset, SubreadSet)

    def __getitem__(self, key):
        """
        Slice by movie name, zmw name, or zmw range name, using standard
        PacBio naming conventions.  Examples:

          - ["m110818_..._s1_p0"]             -> new dataset filtered on movie name
          - ["m110818_..._s1_p0/24480"]       -> BamZmw
          - ["m110818_..._s1_p0/24480/20_67"] -> ZmwRead
          - ["m110818_..._s1_p0/24480/ccs"]   -> CCSZmwRead
        """
        if not isinstance(key, str):
            raise KeyError("%s key %s is not a string." % (self.__class__.__name__, key))

        indices = key.rstrip("/").split("/")

        if len(indices) < 1:
            raise KeyError("%s invalid slice: %s" % (self.__class__.__name__, key))

        _movie = indices[0]
        if len(indices) == 1:
            ds2 = self._dataset.copy()
            ds2.filters.addRequirement(movie=[('==', _movie)])
            return ds2

        _hn = int(indices[1])
        # for efficiency: first select on hn, then movie name.
        _reads = self._dataset[self._dataset.index.holeNumber == _hn]
        _reads = [_read for _read in _reads if _read.movieName == _movie]
        if len(_reads) == 0:
            raise KeyError("Could not find %s in %s" % (key, str(self._dataset)))

        _zmw = BamZmw(bamRecords=_reads, isCCS=self.isCCS)
        if len(indices) == 2:
            return _zmw

        if len(indices) >= 3:
            if indices[2].lower() == "ccs":
                return _zmw.ccsRead
            else:
                start, end = (int(x) for x in indices[2].split("_"))
                return _zmw.read(start, end)

        raise KeyError("%s invalid slice: %s" % (self.__class__.__name__, key))

    def __iter__(self):
        """Iterators over Zmw, ZmwRead objects.
        """
        _qid_hns = np.unique(self._dataset.index[['qId', 'holeNumber']])
        for _movie in self.movieNames:
            ds2 = self._dataset.copy()
            ds2.filters.addRequirement(movie=[('==', _movie)])
            _movie_id = self._dataset.movieIds[_movie]
            for (__movie_id, _hn) in _qid_hns:
                if __movie_id == _movie_id:
                    _reads = ds2[ds2.index.holeNumber == _hn]
                    assert all([_read.movieName == _movie for _read in _reads])
                    if len(_reads) > 0:
                        yield BamZmw(bamRecords=_reads, isCCS=self.isCCS)

    def reads(self):
        """Iterate over all reads"""
        for r in self._dataset:
            yield r

    def subreads(self):
        """Iterate over all subreads."""
        if self.isSubread:
            for r in self._dataset:
                yield r

    def ccsReads(self):
        """Iterate over all ccs reads."""
        if self.isCCS:
            for r in self._dataset:
                yield r

    def __repr__(self):
        return "<%s of movies:\n%s>" % (self.__class__.__name__,
                                        "\n".join(self.movieNames))

    def __delitem__(self, dummy_name):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__delitem__"))

    def __setitem__(self, dummy_index, dummy_name):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__setitem__"))

    def __len__(self):
        """Return total number of zmws in all movies."""
        return len(np.unique(self._dataset.index[['qId', 'holeNumber']]))

    def close(self):
        """Close all readers."""
        self._dataset.close()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        return self.close()

HDHEADERTAG = 'HD'
RGHEADERTAG = 'RG'
SQHEADERTAG = 'SQ'
PGHEADERTAG = 'PG'
HEADERTAGS = [HDHEADERTAG, RGHEADERTAG, SQHEADERTAG, PGHEADERTAG]


class BamHeader(object):
    """
    Represent a bam header.
    """
    def __init__(self, headers=None, ignore_pg=False):
        """
        ignore_pg: whether or not to propogate @PG info.
        """
        self._dict = dict({tag:[] for tag in HEADERTAGS})
        self.ignore_pg = ignore_pg
        if headers is None:
            pass
        elif isinstance(headers, dict):
            self._dict = headers
        elif isinstance(headers, list):
            for _header in headers:
                self.add(_header)
        else:
            raise ValueError("<%s, __init__> does not support input type %s." %
                             (self.__class__.__name__, type(headers)))

        for tag in HEADERTAGS:
            if tag not in self._dict:
                self._dict[tag] = []

    @property
    def headerLine(self):
        """Return the header line.
        e,g,. bh.headerLine = {'SO': 'UNKNOWN', 'VN': '1.5', 'pb': '3.0b5'}
        """
        return self._dict[HDHEADERTAG]

    @property
    def readGroups(self):
        """Return all read groups.
        e.g., bh.readGroups = [
        {'DS': 'READTYPE=CCS;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;BINDINGKIT=100236500;SEQUENCINGKIT=001558034;BASECALLERVERSION=2.3', 'ID': '9e743229','PL': 'PACBIO', 'PU': 'm131018_081703_42161_c100585152550000001823088404281404_s1_p0'},
        {'DS': 'READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;BINDINGKIT=100256000;SEQUENCINGKIT=100254800;BASECALLERVERSION=2.1.0.0.127824', 'ID': '0829d831','PL': 'PACBIO', 'PU': 'm140802_070303_42161_c110036822550000001823106706241502_s1_p0'}
        ]
        """
        return self._dict[RGHEADERTAG]

    @property
    def programGroups(self):
        """Return all program groups."""
        return self._dict[PGHEADERTAG]

    @property
    def referenceSequences(self):
        """Return all reference sequences."""
        return self._dict[SQHEADERTAG]

    @property
    def header(self):
        """Return this BamHeader as a dict"""
        return self._dict

    @property
    def referenceLengthsDict(self):
        """Return {reference_name: reference_length} as a dict"""
        return dict({refseq['SN']:refseq['LN']
                     for refseq in self.referenceSequences})

    def containReadGroup(self, readGroupId):
        """Return whether or not contain this object contains
        a read group with ID=readGroupId.
        """
        return readGroupId in [rg['ID'] for rg in self.readGroups]

    def containReferenceSequence(self, refSeqId):
        """Return whether or not contain this object contains
        a reference sequence with ID=refSeqId.
        """
        return refSeqId in [refseq['SN'] for refseq in self.referenceSequences]

    def _addRG(self, rg):
        """Add a 'RG' entry to header, while rg['ID'] must be unique."""
        if not self.containReadGroup(rg['ID']):
            self.readGroups.append(rg)

    def _addPG(self, pg):
        """Add a 'PG' entry to header."""
        self.programGroups.append(pg)

    def _addSQ(self, sq):
        """Add a 'SQ' entry to header, while sq['SN'] must be unique."""
        if not self.containReferenceSequence(sq['SN']):
            self.referenceSequences.append(sq)

    def add(self, other):
        """Add/Merge another BamHeader or dict to this object."""
        _toadd = other
        if isinstance(other, BamHeader):
            _toadd = other.header
        else:
            if not isinstance(other, dict):
                raise TypeError("<%s, add> does not support type %s" %
                                (self.__class__.__name__, type(other)))
            _toadd = BamHeader(other).header

        if len(self.headerLine) == 0:
            self._dict[HDHEADERTAG] = _toadd[HDHEADERTAG]

        for rg in _toadd[RGHEADERTAG]:
            self._addRG(rg)

        if not self.ignore_pg:
            for pg in _toadd[PGHEADERTAG]:
                self._addPG(pg)
        for sq in _toadd[SQHEADERTAG]:
            self._addSQ(sq)

    def __repr__(self):
        maxn = 100
        return "<%s: %d RG, %d SN>\n" % (
            self.__class__.__name__, len(self.readGroups),
            len(self.referenceSequences)) + \
            "RGs [%s]\n" % (", ".join([rg['ID'] for rg in self.readGroups])[0:maxn]) + \
            "SNs [%s]\n" % (", ".join([sn['SN'] for sn in self.referenceSequences])[0:maxn])


class BamWriter(object):
    """
    Write unaligned bam records to a bam file.
    """

    def __init__(self, fname, header):
        self.filename = op.abspath(op.expanduser(fname))
        if isinstance(header, dict):
            self.peer = pysam.AlignmentFile(fname, "wb", header=header)
        elif isinstance(header, BamHeader):
            self.peer = pysam.AlignmentFile(fname, "wb", header=header.header)
        else:
            raise ValueError("<%s, __init__(fname, header)> " %
                             self.__class__.__name__ +
                             "header type must be either dict or BamHeader.")

    def close(self):
        """Close bam file."""
        self.peer.close()

    def write(self, record):
        """Write a record or a list of records.
        If record type is one of the followings,
            * pysam.calignmentfile.AlignedSegment
            * BamZmwRead
            * BamCCSZmwRead
        write this record.
        """
        # workaround for changes to the internal structure of pysam -
        # could probably be improved upon
        segment_class = getattr(pysam.calignmentfile, "AlignedSegment", None)
        if segment_class is None:
            segment_class = pysam.calignedsegment.AlignedSegment
        if isinstance(record, segment_class):
            self.peer.write(record)
        elif isinstance(record, BamCCSZmwRead) or \
             isinstance(record, BamZmwRead):
            self.peer.write(record.alignedSegment)
        else:
            raise ValueError("<%s, write()> does not support record type %s." %
                             (self.__class__.__name__, type(record)))

    def __repr__(self):
        return "<%s: %s>" % (self.__class__.__name__,
                             self.filename)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class CCSBamSequence(object):
    """
    Wrapper for a BAM record, used to substitute for a Fasta record.
    """
    def __init__(self, bam):
        self.name = bam.qname
        self.sequence = bam.seq


class CCSInput(object):
    """
    Wrapper class for handling multiple formats specifying CCS sequences.
    The old convention was to use .fasta, but we would like to be able to pass
    the classifier a ConsensusReadSet (i.e. .bam files) instead for use within
    pbsmrtpipe.
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self._is_fasta = False
        self.ext = op.splitext(file_name)[1].upper()
        if self.ext in [".FA", ".FASTA"]:
            self._dataset = FastaReader(file_name)
            self._is_fasta = True
        elif self.ext == ".BAM":
            self._dataset = openDataFile(file_name)
        else: # either contigset.xml or consensusreadset.xml
            assert self.ext == ".XML"
            self._dataset = openDataSet(file_name)
            if isinstance(self._dataset, ContigSet):
                self._is_fasta = True

    def __iter__(self):
        for rec in self._dataset:
            if not self._is_fasta:
                rec = CCSBamSequence(rec.peer)
            yield rec

    def close(self):
        """Close all datasets."""
        self._dataset.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __len__(self):
        if not self._is_fasta:
            return len(self._dataset)
        else:
            if self.ext in [".FA", ".FASTA"]:
                return len([r for r in FastaReader(self.file_name)])
            else: # contigset
                n = 0
                for rr in self._dataset.resourceReaders():
                    n += len([r for r in rr])
                return n

    def __delitem__(self, dummy_name):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__delitem__"))

    def __setitem__(self, dummy_index, dummy_name):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__setitem__"))

    def __getitem__(self, key):
        raise NotImplementedError("%s.%s" % (self.__class__.__name__,
                                             "__getitem__"))

