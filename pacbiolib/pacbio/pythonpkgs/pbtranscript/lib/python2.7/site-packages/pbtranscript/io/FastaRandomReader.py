"""
Define FastaRandomReader which accesses reads in FASTA
or ContigSet files randomly.

Note that a ContigSet can contain multiple FASTA files,
but can only contain FASTA files and filters in the
ContigSet will not be respected.
"""

from collections import namedtuple
import os.path as op
import logging

from pbcore.io.FastaIO import FastaRecord
from pbcore.io import ContigSet

__all__ = ["FastaRandomReader",
           "MetaSubreadFastaReader",
           "SubreadFastaReader"]

Interval = namedtuple('Interval', ['start', 'end'])


class FastaRandomReader(object):

    """
        A reader class to access sequences in a FASTA or ContigSet XML
        file randomely.

        Note that this class can take ContigSet file as input, however
        filters are not respected.

        This is meant to substitute for the Bio.SeqIO.to_dict method since some fasta files
        are too big to fit entirely to memory.

        The only requirement is that every id line begins with the symbol >.

        Example:
            r = FastaRandomReader('output/test.fna')
            r['6C_49273_NC_008578/2259031-2259297'] ==> this shows a FastaRecord

            r = FastaRandomReader('test.contigset.xml')
            r['m/zmw/s_e'] ==> this shows a FastaRecord
    """

    def __init__(self, *args):
        self.fasta_filenames = self.get_fasta_filenames(*args)
        self.fhandlers = self._open_files()
        self.d = {}
        self._init_index()

    def get_fasta_filenames(self, *args):
        """Return all FASTA file names as a list."""
        ret = []
        for fn in args:
            ext = op.splitext(fn)[1].upper()
            if ext in [".XML"]:
                try:
                    fns = ContigSet(fn).toExternalFiles()
                    ret.extend(fns)
                except IOError:
                    raise IOError("Could not open %s as ContigSet." % fn)
            else:
                ret.append(fn)
        for fn in ret:
            if op.splitext(fn)[1].upper() not in [".FA", ".FASTA"]:
                raise IOError("%s input must be FASTA or ContigSet."
                              % self.__class__.__name__)
        return ret

    def _open_files(self):
        return [open(fn) for fn in self.fasta_filenames]

    def _init_index(self):
        """Indexing reads in fasta_filenames"""
        for index, f in enumerate(self.fhandlers):
            while 1:
                line = f.readline()
                if len(line) == 0:
                    break
                if line.startswith('>'):
                    sid = line.strip()[1:].split(None, 1)[0]
                    # the header MUST be just 1 line
                    # if id in self.d:
                    #    print "duplicate id {0}!!".format(id)
                    self.d[sid] = (index, f.tell())

    def __getitem__(self, k):
        if k not in self.d:
            errMsg = "key {k} not in {f}!".format(k=k, f=",".join(self.fasta_filenames))
            logging.error(self.d.keys())
            raise ValueError(errMsg)
        return self._get_record(k)

    def _get_record(self, k):
        index, tell = self.d[k]
        f = self.fhandlers[index]
        f.seek(tell)
        content = ''
        for line in f:
            if line.startswith('>'):
                break
            content += line.strip()
        return FastaRecord(header=k, sequence=content)

    def __len__(self):
        return len(self.d)

    def __delitem__(self, key):
        errMsg = "%s.__delitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)

    def __setitem__(self, key, value):
        errMsg = "%s.__setitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)

    def keys(self):
        """Return d.keys."""
        return self.d.keys()


class MetaSubreadFastaReader(object):

    """Reader for reading PabBio subreads in a list of fasta files."""

    def __init__(self, fasta_filenames):
        self.meta_f = {}
        self.meta_zmw_d = {}
        # record just the zmw and let the subread reader handle it
        for fn in fasta_filenames:
            self.meta_f[fn] = SubreadFastaReader(fn)
            # combine all the keys
            for k in self.meta_f[fn].zmw_d.keys():
                self.meta_zmw_d[k] = fn

    def __getitem__(self, k):
        """
        k -- could be zmw or subread id
        """
        if k.count('/') == 2:
            zmw = k[:k.rfind('/')]
        else:
            zmw = k
        return self.meta_f[self.meta_zmw_d[zmw]][k]

    def __len__(self):
        return len(sum([len(d) for d in self.meta_zmw_d]))

    def __delitem__(self, key):
        errMsg = "%s.__delitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)

    def __setitem__(self, key, value):
        errMsg = "%s.__setitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)


class SubreadFastaReader(object):

    """Reader for reading PabBio subreads in a fasta file."""

    def __init__(self, fasta_filename):
        self.f = open(fasta_filename)
        self.d = {}
        self.zmw_d = {}
        while 1:
            line = self.f.readline()
            if len(line) == 0:
                break
            if line.startswith('>'):
                # ex: m140..._s1_p0/155/0_1673 RQ=0.845
                rid = line.strip()[1:].split(None, 1)[0]
                # the header MUST be just 1 line
                zmw = rid[:rid.rfind('/')]
                self.d[rid] = (rid, self.f.tell())
                if zmw not in self.zmw_d:
                    self.zmw_d[zmw] = []
                self.zmw_d[zmw].append((rid, self.f.tell()))

    def __getitem__(self, k):
        """
        k --- should be <movie>/<zmw> or <movie>/<zmw>/<start_end>
        If former, return a list of records associated with that ZMW
        If latter, return just that record but still in a list
        """
        if k.count('/') == 2:  # is a subread
            if k not in self.d:
                raise ValueError("key {0} not in dictionary!".format(k))
            locations = [self.d[k]]
        else:  # is a ZMW
            if k not in self.zmw_d:
                raise ValueError("key {0} not in dictionary!".format(k))
            locations = self.zmw_d[k]
        output = []
        for seqid, loc in locations:
            self.f.seek(loc)
            content = ''
            for line in self.f:
                if line.startswith('>'):
                    break
                content += line.strip()
            output.append(FastaRecord(header=seqid, sequence=content))
        return output

    def keys(self):
        """return keys (subreads)."""
        return self.d.keys()

    def __len__(self):
        return len(self.d)

    def __delitem__(self, key):
        errMsg = "%s.__delitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)

    def __setitem__(self, key, value):
        errMsg = "%s.__setitem__ not defined." % self.__class__.__name__
        raise NotImplementedError(errMsg)
