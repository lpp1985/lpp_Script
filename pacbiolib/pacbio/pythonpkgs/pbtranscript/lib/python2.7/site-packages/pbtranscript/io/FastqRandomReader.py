"""
Define miscBio.FastqRandomReader.
It will be repaced by pbcore.io.FastqIO.FastqRandomReader.
"""
from pbcore.io.FastqIO import FastqRecord

__all__ = ["FastqRandomReader"]

class FastqRandomReader:

    """
        Like FastaRandomReader except works with fastq.

        This is meant to substitute for the Bio.SeqIO.to_dict method since some fastq files
        are too big to fit entirely to memory. The only requirement is that every id line
        begins with the symbol >. It is ok for the sequences to stretch multiple lines.
        The sequences, when read, are returned as Bio.SeqRecord objects.

        Example:
            r = FastqRandomReader('test.fastq')
            r['m131018_081703_42161_c100585152550000001823088404281404_s1_p0/61232/4045_63_CCS'] ==> this shows a FastqRecord
    """

    def __init__(self, fastq_filename):
        self.f = open(fastq_filename)
        self.d = {}
        self.locations = []
        self.locations_d_key_map = {}

        idx = -1
        while 1:
            line = self.f.readline()
            if len(line) == 0:
                break
            idx += 1
            if idx % 4 == 0:
                assert line.startswith('@')
                sid = line.strip()[1:].split(None, 1)[0]
                # the header MUST be just 1 line
                if sid in self.d:
                    raise ValueError("Bad fastq format: redundant record {sid}".
                                     format(sid=sid))
                self.d[sid] = self.f.tell()

    def __getitem__(self, k):
        if k not in self.d:
            errMsg = "key {k} not in {f}!".format(k=k, f=self.f.name)
            raise ValueError(errMsg)
        self.f.seek(self.d[k])
        name, seq, plus, qv = k, "", "", ""

        seq  = self.f.readline().strip()
        plus = self.f.readline().strip()
        qv = self.f.readline().strip()
        if plus != '+':
            raise ValueError("Bad fastq format: line 3 of {name}".
                             format(name=name) + " is not '+'")
        return FastqRecord(header=name, sequence=seq, qualityString=qv)

    def __len__(self):
        return len(self.d)

    def __delitem__(self, key):
        errMsg = "FastqRandomReader.__delitem__ not defined."
        raise NotImplementedError(errMsg)

    def __setitem__(self, key):
        errMsg = "FastqRandomReader.__setitem__ not defined."
        raise NotImplementedError(errMsg)

    def keys(self):
        """Return d.keys."""
        return self.d.keys()
