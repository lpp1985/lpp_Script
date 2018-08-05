"""
Streaming I/O support for FASTA files.

# Author: David Alexander
"""
from __future__ import absolute_import
from collections import namedtuple, OrderedDict, Sequence
from os.path import abspath, expanduser, isfile, getsize
from string import maketrans
import gzip
import mmap
import re
__all__ = [ "FastaRecord",
            "FastaReader",
            "FastaWriter",
            "FastaTable",
            "IndexedFastaReader",
            "splitFastaHeader",
            "splitFileContents",
            "complement",
            "reverseComplement",
            "ReaderBase", "WriterBase",
]

DNA_COMPLEMENT = maketrans('agcturyswkmbdhvnAGCTURYSWKMBDHV-N',
                           'tcgannnnnnnnnnnnTCGANNNNNNNNNNN-N')

def reverse(sequence):
    """Return the reverse of any sequence
    """
    return sequence[::-1]

def complement(sequence):
    """
    Return the complement of a sequence
    """
    if re.search('[^agcturyswkmbdhvnAGCTURYSWKMBDHVN-]', sequence):
        raise ValueError("Sequence contains invalid DNA characters - "
                         "only standard IUPAC nucleotide codes allowed")
    return sequence.translate(DNA_COMPLEMENT)

def reverseComplement(sequence):
    """
    Return the reverse-complement of a sequence
    NOTE: This only currently supports DNA
    """
    return complement(sequence)[::-1]

def splitFileContents(f, delimiter, BLOCKSIZE=8192):
    """
    Same semantics as f.read().split(delimiter), but with memory usage
    determined by largest chunk rather than entire file size
    """
    from cStringIO import StringIO
    remainder = StringIO()
    while True:
        block = f.read(BLOCKSIZE)
        if not block:
            break
        parts = block.split(delimiter)
        remainder.write(parts[0])
        for part in parts[1:]:
            yield remainder.getvalue()
            remainder = StringIO()
            remainder.write(part)
    yield remainder.getvalue()

def splitFastaHeader( name ):
    """
    Split a FASTA/FASTQ header into its id and comment components
    """
    nameParts = re.split('\s', name, maxsplit=1)
    id_ = nameParts[0]
    if len(nameParts) > 1:
        comment = nameParts[1].strip()
    else:
        comment = None
    return (id_, comment)

class FastaRecord(object):
    """
    A FastaRecord object models a named sequence in a FASTA file.
    """
    DELIMITER = ">"
    COLUMNS   = 60

    def __init__(self, header, sequence):
        try:
            assert "\n" not in header
            assert "\n" not in sequence
            assert self.DELIMITER not in sequence
            self._header = header
            self._sequence = sequence
            self._id, self._comment = splitFastaHeader(header)
        except AssertionError:
            raise ValueError("Invalid FASTA record data")

    @property
    def header(self):
        """
        The header of the sequence in the FASTA file, equal to the entire
        first line of the FASTA record following the '>' character.

        .. warning::

           You should almost certainly be using "id", not "header".
        """
        return self._header

    @property
    def name(self):
        """
        DEPRECATED: The name of the sequence in the FASTA file, equal to
        the entire FASTA header following the '>' character
        """
        return self._header

    @property
    def id(self):
        """
        The id of the sequence in the FASTA file, equal to the FASTA header
        up to the first whitespace.
        """
        return self._id

    @property
    def comment(self):
        """
        The comment associated with the sequence in the FASTA file, equal to
        the contents of the FASTA header following the first whitespace
        """
        return self._comment

    @property
    def sequence(self):
        """
        The sequence for the record as present in the FASTA file.
        (Newlines are removed but otherwise no sequence normalization
        is performed).
        """
        return self._sequence

    @classmethod
    def fromString(cls, s):
        """
        Interprets a string as a FASTA record.  Does not make any
        assumptions about wrapping of the sequence string.
        """
        try:
            lines = s.splitlines()
            assert len(lines) > 1
            assert lines[0][0] == cls.DELIMITER
            header = lines[0][1:]
            sequence = "".join(lines[1:])
            return FastaRecord(header, sequence)
        except AssertionError:
            raise ValueError("String not recognized as a valid FASTA record")

    def reverseComplement(self, preserveHeader=False):
        """
        Return a new FastaRecord with the reverse-complemented DNA sequence.
        Optionally, supply a name
        """
        rcSequence = reverseComplement(self.sequence)
        if preserveHeader:
            return FastaRecord(self.header, rcSequence)
        else:
            rcName = '{0} [revcomp]'.format(self.header.strip())
            return FastaRecord(rcName, rcSequence)

    def __len__(self):
        return len(self._sequence)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.header   == other.header and
                    self.sequence == other.sequence)
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "<FastaRecord: %s>" % self.header

    def __str__(self):
        """
        Output a string representation of this FASTA record, observing
        standard conventions about sequence wrapping.
        """
        return (">%s\n" % self.header) + \
            wrap(self.sequence, self.COLUMNS)

def isFileLikeObject(o):
    return hasattr(o, "read") and hasattr(o, "write")

def getFileHandle(filenameOrFile, mode="r"):
    """
    Given a filename not ending in ".gz", open the file with the
    appropriate mode.

    Given a filename ending in ".gz", return a filehandle to the
    unzipped stream.

    Given a file object, return it unless the mode is incorrect--in
    that case, raise an exception.
    """
    assert mode in ("r", "w")

    if isinstance(filenameOrFile, basestring):
        filename = abspath(expanduser(filenameOrFile))
        if filename.endswith(".gz"):
            return gzip.open(filename, mode)
        else:
            return open(filename, mode)
    elif isFileLikeObject(filenameOrFile):
        return filenameOrFile
    else:
        raise Exception("Invalid type to getFileHandle")

class ReaderBase(object):
    def __init__(self, f):
        """
        Prepare for iteration through the records in the file
        """
        self.file = getFileHandle(f, "r")
        if hasattr(self.file, "name"):
            self.filename = self.file.name
        else:
            self.filename = "(anonymous)"

    def close(self):
        """
        Close the underlying file
        """
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __repr__(self):
        return "<%s for %s>" % (type(self).__name__, self.filename)

class WriterBase(object):
    def __init__(self, f):
        """
        Prepare for output to the file
        """
        self.file = getFileHandle(f, "w")
        if hasattr(self.file, "name"):
            self.filename = self.file.name
        else:
            self.filename = "(anonymous)"

    def close(self):
        """
        Close the underlying file
        """
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __repr__(self):
        return "<%s for %s>" % (type(self).__name__, self.filename)

class FastaReader(ReaderBase):
    """
    Streaming reader for FASTA files, useable as a one-shot iterator
    over FastaRecord objects.  Agnostic about line wrapping.

    Example:

    .. doctest::

        >>> from pbcore.io import FastaReader
        >>> from pbcore import data
        >>> filename = data.getTinyFasta()
        >>> r = FastaReader(filename)
        >>> for record in r:
        ...     print record.header, len(record.sequence)
        ref000001|EGFR_Exon_2 183
        ref000002|EGFR_Exon_3 203
        ref000003|EGFR_Exon_4 215
        ref000004|EGFR_Exon_5 157
        >>> r.close()

    """
    DELIMITER = ">"

    def __iter__(self):
        try:
            parts = splitFileContents(self.file, ">")
            assert "" == next(parts)
            for part in parts:
                yield FastaRecord.fromString(">" + part)
        except AssertionError:
            raise ValueError("Invalid FASTA file {f}".format(f=self.filename))


class FastaWriter(WriterBase):
    """
    A FASTA file writer class

    Example:

    .. doctest::

        >>> from pbcore.io import FastaWriter
        >>> with FastaWriter("output.fasta.gz") as writer:
        ...     writer.writeRecord("dog", "GATTACA")
        ...     writer.writeRecord("cat", "CATTACA")

    (Notice that underlying file will be automatically closed after
    exit from the `with` block.)

    .. testcleanup::

        import os; os.unlink("output.fasta.gz")

    """
    def writeRecord(self, *args):
        """
        Write a FASTA record to the file.  If given one argument, it is
        interpreted as a ``FastaRecord``.  Given two arguments, they
        are interpreted as the name and the sequence.
        """
        if len(args) not in (1, 2):
            raise ValueError
        if len(args) == 1:
            record = args[0]
        else:
            header, sequence = args
            record = FastaRecord(header, str(sequence))
        self.file.write(str(record))
        self.file.write("\n")


##
## Utility functions for FastaReader
##
def wrap(s, columns):
    return "\n".join(s[start:start+columns]
                     for start in xrange(0, len(s), columns))



# ------------------------------------------------------------------------------
# IndexedFastaReader: random access Fasta class
#

FaiRecord = namedtuple("FaiRecord", ("id", "comment", "header", "length", "offset", "lineWidth", "stride"))

def faiFilename(fastaFilename):
    return fastaFilename + ".fai"

def loadFastaIndex(faidxFilename, fastaView):

    if not isfile(faidxFilename): # os.path.isfile
        raise IOError("Companion FASTA index (.fai) file not found or "
                      "malformatted! Use 'samtools faidx' to generate FASTA "
                      "index.")

    tbl = []
    # NB: We have to look back in the FASTA to find the full header;
    # only "id" makes it into the fai.
    offsetEnd = 0
    for line in open(faidxFilename):
        length, offset, lineWidth, blen = map(int, line.split()[-4:])
        newlineWidth = blen - lineWidth                                # 2 for DOS, 1 for UNIX
        header_    = fastaView[offsetEnd:offset]
        assert (header_[0] == ">" and header_[-1] == "\n")
        header     = header_[1:-newlineWidth]
        id, comment = splitFastaHeader(header)
        q, r = divmod(length, lineWidth)
        numNewlines = q + (r > 0)
        offsetEnd = offset + length + numNewlines*newlineWidth
        record = FaiRecord(id, comment, header, length, offset, lineWidth, blen)
        tbl.append(record)
    return tbl

def fileOffset(faiRecord, pos):
    """
    Find the in-file position (in bytes) corresponding to the position
    in the named contig, using the FASTA index.
    """
    q, r = divmod(pos, faiRecord.lineWidth)
    offset = faiRecord.offset + q*faiRecord.stride + r
    return offset

class MmappedFastaSequence(Sequence):
    """
    A string-like view of a contig sequence that is backed by a file
    using mmap.
    """
    def __init__(self, view, faiRecord):
        self.view = view
        self.faiRecord = faiRecord

    def __getitem__(self, spec):
        if isinstance(spec, slice):
            start, stop, stride = spec.indices(len(self))
            if stride != 1:
                raise ValueError, "Unsupported stride"
        elif spec < 0:
            start = self.faiRecord.length + spec
            stop = start + 1
            stride = 1
        else:
            start = spec
            stop = start + 1
            stride = 1
        if not (0 <= start <= stop <= self.faiRecord.length):
            raise IndexError, "Out of bounds"
        startOffset = fileOffset(self.faiRecord, start)
        endOffset   = fileOffset(self.faiRecord, stop)
        snip = self.view[startOffset:endOffset].translate(None, "\r\n")
        return snip

    def __len__(self):
        return self.faiRecord.length

    def __eq__(self, other):
        return (isinstance(other, MmappedFastaSequence) and
                self[:] == other[:])

    def __str__(self):
        return str(self[:])


class IndexedFastaRecord(object):

    COLUMNS   = 60

    def __init__(self, view, faiRecord):
        self.view = view
        self.faiRecord = faiRecord

    @property
    def name(self):
        return self.header

    @property
    def header(self):
        return self.faiRecord.header

    @property
    def id(self):
        return self.faiRecord.id

    @property
    def comment(self):
        return self.faiRecord.comment

    @property
    def sequence(self):
        return MmappedFastaSequence(self.view, self.faiRecord)

    def __len__(self):
        return self.faiRecord.length

    def __repr__(self):
        return "<IndexedFastaRecord: %s>" % self.header

    def __eq__(self, other):
        return (isinstance(other, IndexedFastaRecord) and
                self.header == other.header and
                self.sequence == other.sequence)

    def __str__(self):
        """
        Output a string representation of this FASTA record, observing
        standard conventions about sequence wrapping.
        """
        return (">%s\n" % self.header) + \
            wrap(self.sequence, self.COLUMNS)

class IndexedFastaReader(ReaderBase, Sequence):
    """
    Random-access FASTA file reader.

    Requires that the lines of the FASTA file be fixed-length and that
    there is a FASTA index file (generated by `samtools faidx`) with
    name `fastaFilename.fai` in the same directory.

    .. doctest::

        >>> from pbcore.io import FastaTable
        >>> from pbcore import data
        >>> filename = data.getFasta()
        >>> t = IndexedFastaReader(filename)
        >>> print t[:4] # doctest: +NORMALIZE_WHITESPACE
        [<IndexedFastaRecord: ref000001|EGFR_Exon_2>,
         <IndexedFastaRecord: ref000002|EGFR_Exon_3>,
         <IndexedFastaRecord: ref000003|EGFR_Exon_4>,
         <IndexedFastaRecord: ref000004|EGFR_Exon_5>]
        >>> t.close()

    """
    def __init__(self, filename):
        self.filename = abspath(expanduser(filename))
        self.file = open(self.filename, "r")
        self.faiFilename = faiFilename(self.filename)
        if getsize(self.filename) > 0:
            self.view = mmap.mmap(self.file.fileno(), 0,
                                  prot=mmap.PROT_READ)
            self.fai = loadFastaIndex(self.faiFilename, self.view)
        else:
            self.view = None
            self.fai = []
        self.contigLookup = self._loadContigLookup()

    def _loadContigLookup(self):
        contigLookup = dict()
        for (pos, faiRecord) in enumerate(self.fai):
            contigLookup[pos]              = faiRecord
            contigLookup[faiRecord.id]     = faiRecord
            contigLookup[faiRecord.header] = faiRecord
        return contigLookup

    def __getitem__(self, key):
        if key < 0:
            key = len(self) + key

        if isinstance(key, slice):
            indices = xrange(*key.indices(len(self)))
            return [ IndexedFastaRecord(self.view, self.contigLookup[i])
                     for i in indices ]
        elif key in self.contigLookup:
            return IndexedFastaRecord(self.view, self.contigLookup[key])
        else:
            raise IndexError, "Contig not in FastaTable"

    def __iter__(self):
        return (self[i] for i in xrange(len(self)))

    def __len__(self):
        return len(self.fai)

# old name for IndexedFastaReader was FastaTable
FastaTable = IndexedFastaReader
