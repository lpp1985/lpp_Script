
"""
I/O handling for BAM alignment files.
"""

import warnings
import copy
import os.path as op
import os
import logging
import re

import numpy as np
import pysam

from pbcore.io import IndexedBamReader
import sys

log = logging.getLogger(__name__)


class MovieIdx(object):

    """ Provides a simple way to index into the numpy arrays by movie """

    def __init__(self, name, rOffs=0, rLen=0, sOffs=0, sLen=0):
        self.name = name
        # read indicies
        self.rOffs = rOffs
        self.rLen = rLen
        # sub-read indicies
        self.sOffs = sOffs
        self.sLen = sLen

    @property
    def shortname(self):
        # only matches Springfield movies
        match = re.search('m\d+_\d+_\w+\_(c\d+_s\d_p\d)', self.name)
        if match:
            # shorten the name (displayed to user)
            cellid = match.group(1)
            return cellid
        else:
            # astro movie names don't match, just display as is
            # m.shortname = name
            return self.name

    def __repr__(self):
        return "<MovieIdx name:{n}>".format(n=self.name)


def _getPct(percentile, vector):
    """Returns the specified percentile value of the numpy vector"""
    sorted_vector = np.sort(vector)
    index = np.ceil((percentile / 100.0) * len(sorted_vector))
    return sorted_vector[-1] if index >= len(sorted_vector) else sorted_vector[index]


class MovieAlignmentInfo(object):

    def __init__(self, bam_file_name, movie_name):
        self.bam_file_name = bam_file_name
        self.movie_name = movie_name
        self.datum = {}
        self.unrolled = {}
        self.max_subread = {}

    def as_tuple(self):
        return self.datum, self.unrolled, self.max_subread, set([self.movie_name])


def alignment_info_from_bam(bam_file_name):
    """
    Extract subread information from an indexed BAM file.  This should be
    relatively fast since it will not access the BAM records directly.
    """
    by_movie = {}
    last_zmw_id = None
    with IndexedBamReader(bam_file_name) as bam:
        if len(bam) > 0:
            identities = bam.identity
            subread_lengths = bam.aEnd - bam.aStart
            for i_aln, rgId in enumerate(bam.qId):
                movie_name = bam.readGroupInfo(rgId).MovieName
                if not movie_name in by_movie:
                    by_movie[movie_name] = MovieAlignmentInfo(bam_file_name,
                                                              movie_name)
                m = by_movie[movie_name]
                hole_number = bam.holeNumber[i_aln]
                qs, qe = bam.qStart[i_aln], bam.qEnd[i_aln]
                rstart, rend = bam.aStart[i_aln], bam.aEnd[i_aln]
                identity = None
                if (qs, qe) == (-1, -1):
                    qs = 0
                    # XXX This is only used to key subreads so the exact value is
                    # not important - still clumsy though
                    qe = rend - rstart

                # Compound ids
                zmw_id = (movie_name, hole_number)
                subread_id = (movie_name, hole_number, qs, qe)

                this_a = []
                this_a.append(subread_lengths[i_aln])

                this_a.append(identities[i_aln])
                this_a.append(bam.readQual[i_aln])

                this_a.append(1.0 if zmw_id != last_zmw_id else 0.0)  # isFirst

                # modStart, a value without a clear meaning, so just write some
                # garbage
                this_a.append(99999)

                last_zmw_id = zmw_id

                if subread_id in m.datum:
                    warnings.warn("Duplicate subread %s" % str(subread_id))

                # No Z-score
                m.datum[subread_id] = tuple(this_a)

                if zmw_id not in m.max_subread or subread_lengths[i_aln] > m.max_subread[zmw_id][1]:
                    m.max_subread[zmw_id] = (
                        subread_id, subread_lengths[i_aln])

                m.unrolled.setdefault(zmw_id, [99999, 0])
                m.unrolled[zmw_id][0] = min(m.unrolled[zmw_id][0], rstart)
                m.unrolled[zmw_id][1] = max(m.unrolled[zmw_id][1], rend)

    return by_movie  # datum, unrolled, max_subread


def from_alignment_file(aln_info):  # movie_name, alignment_file_name):
    columns = ["Length", "Concordance", "Read quality", "isFirst", "modStart"]
    datum, unrolled, max_subread, movie_names = aln_info.as_tuple()
    return movie_names, unrolled, datum, columns


class CrunchedAlignments(object):

    """
    Holds crunched data from the AlignmentSet that's ready to move on for
    display, provides some convenience data accessors.

    mnames -
        A set of movie names covered by this dataset
    reads -
        Alignment information for the reads
              key: (movie name, hole number)
            value: [length, accuracy, is first, mod start]

    subreads -
        Alignment information for the subreads
              key: (movie name, hole number, start, end)
            value: [start, end]

    cols -
        Column names represented in the subread numpy representation
    """

    def __init__(self, mnames, reads, subreads, cols):
        self._movieNames = mnames
        self._reads = reads
        self._subreads = subreads
        self._cols = cols
        # list of MovieIdx instances
        self._movies = []
        # np.array
        self._nReads = None
        # np.array
        self._nSubreads = None
        #  updates _movies, _nReads, _nSubreads
        self._numpify()

    @property
    def movies(self):
        return self._movies

    def _numpify(self):
        """ Create numpy representations of the data """
        # read offset, subread offset
        rOffs, sOffs = 0, 0
        # read list
        rl = []
        # subread list
        sl = []

        for name in self._movieNames:
            m = MovieIdx(name)

            # index the reads
            m.rOffs = rOffs
            l = [v[1] - v[0]
                 for k, v in self._reads.iteritems() if k[0] == name]
            m.rLen = len(l)
            rOffs = rOffs + m.rLen
            rl.extend(l)

            # index the sub-reads
            m.sOffs = sOffs
            l = [v for k, v in self._subreads.iteritems() if k[0] == name]
            m.sLen = len(l)
            sOffs = sOffs + m.sLen
            sl.extend(l)

            self._movies.append(m)

        self._nReads = np.array(rl)
        if len(sl[0]) != len(self._cols):
            msg = "Length of sl {n} len of columns {c} is incompatible.".format(
                n=len(sl[0]), c=len(self._cols))
            sys.stderr.write(msg + "\n")
            raise IndexError(msg)

        self._nSubreads = np.array(
            sl, dtype=[(col, np.float64) for col in self._cols])

    def reads(self, movie=None):
        """
        Numpy representation of reads, optionally by movie.
        """
        if movie:
            s = movie.rOffs
            e = movie.rOffs + movie.rLen
            return self._nReads[s:e]
        else:
            return self._nReads

    def subreads(self, movie=None):
        """
        Numpy representation of subreads, optionally by movie.  Array of alignment lengths.
        """
        if movie:
            s = movie.sOffs
            e = movie.sOffs + movie.sLen
            return self._nSubreads[s:e]
        else:
            return self._nSubreads
