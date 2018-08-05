
# TODO modernize or delete

"""Script to break down post-mapping statistics by movie.

OLD Pre-bax files

m130715_185638_SMRT1_c000000062559900001500000112311501_s1_p0.bas.h5

m{DATE_DATE}_{INSTRUMENT_NAME}_c{CHIP_STRIP_BARCODE}{CELL_NUMBER}_s{
SET_NUMBER}_{EXPIRED_OR_NOT}{STROBE_NUMBER}

Multi-part files (bax)

m130715_185638_SMRT1_c000000062559900001500000112311501_s1_p0.1.bax.h5

m{DATE_DATE}_{INSTRUMENT_NAME}_c{CHIP_STRIP_BARCODE}{CELL_NUMBER}_s{
SET_NUMBER}_{EXPIRED_OR_NOT}{STROBE_NUMBER}.{MOVIE_PART_ID}

Date: yymmdd_hhmmss of when the movie context was created
Instrument name: e.g. Richard, 42142, etc.
Chip strip barcode: 32 digits
Cell number: 0-7, indicates location of chip in strip
Set number: 1 or 2, with 150k instruments, now just 1
p or X: Indicates whether expired barcodes were used (X) or not (p)
Strobe number: Deprecated, always 0
Movie Part Id: 1-3

"""
import sys
import os
import re
import math
import time
import argparse
import functools
import logging
from pprint import pformat

import numpy as np

from pbcore.util.Process import backticks
from pbcore.io import CmpH5Reader, CmpH5Alignment

__version__ = '2.0'

log = logging.getLogger(__name__)


NAME_PARSER = re.compile(r'x([-0-9]+)_y([-0-9]+)_(\d+)-(\d\d\d\d)_([^|]+)')
MOVIE_PARSER = re.compile(r'm\d+_\d+_([^_]+)_([^_]+)_')

# Use to Validate the pls/bas/bax/rgn name.
MOVIE_NAME2 = re.compile(
    r'(m\d{6}_\d{6}_.+_c[0-9]{32}[0-7]_s(1|2)_(p|X)0).(pls|bas|rgn).h5')

MOVIE_NAME3 = re.compile(
    r'(m\d{6}_\d{6}_.+_c[0-9]{32}[0-7]_s(1|2)_(p|X)0).(1|2|3).(pls|bax|bas|rgn).h5')

# This will include the Movie part id.
MOVIE_PART_NAME = re.compile(
    r'(m\d{6}_\d{6}_.+_c[0-9]{32}[0-7]_s(1|2)_(p|X)0.(1|2|3)).(pls|bax|bas|rgn).h5')

# Old pre multi-part files
MOVIE_PARSER2 = re.compile(
    r'm(\d{6}_\d{6})_(.+)_c([0-9]{32})([0-7])_s(1|2)_(p|X)0.(pls|bas|rgn).h5')
# Multi-part bax era files
MOVIE_PARSER3 = re.compile(
    r'm(\d{6}_\d{6})_(.+)_c([0-9]{32})([0-7])_s(1|2)_(p|X)0.(1|2|3).(pls|bax|bas|rgn).h5')

# FROM MJ
MOVIE_ASTRO_NAME = re.compile(
    r'(m\d{6}_\d{6}_[A-Z,a-z]{3}_p(\d|\d{2})_b(\d{2}|\d)).(pls|bas|rgn).h5')

MOVIE_ASTRO_PARSER = re.compile(
    r'(m\d{6}_\d{6})_[A-Z,a-z]{3}_p(\d|\d{2})_b(\d{2}|\d).(pls|bas|rgn).h5')


PSEUDO_PARSER = re.compile(r'(_F[^|]+)|(\.f[^|]+)$')

# To Grab the internal LIMS RunCode
# /mnt/secondary-siv/testdata/LIMS/2311209/0001/Analysis_Results/m120201_042231_42129_c100275262550000001523007907041260_s1_p0.bas.h5
EXPT_RUN_PARSER = re.compile(r'/(\d\d\d\d\d\d\d)/(\d\d\d\d)/')

QV_THRESHOLDS = [6, 8, 10]
NULL_QV_STRING = ','.join(['0' for i in xrange(len(QV_THRESHOLDS))])
Z_THRESHOLD = 3.0

# length to reuse when converting Z to accuracy
FIXED_LENGTH = 500
# fixed parameters for Z (simple model)
Z_A = -77.27
Z_C = 0.08540
Z_S = 0.001210


def setup_log(alog, file_name=None, level=logging.DEBUG, str_formatter=None):
    """
    Util function for setting up logging.

    Due to how smrtpipe logs, the default behavior is that the stdout
    is where the logging is redirected. If a file name is given the log
    will be written to that file.

    :param log: (log instance) Log instance that handlers and filters will
    be added.
    :param file_name: (str, None), Path to file. If None, stdout will be used.
    :param level: (int) logging level
    """
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)

    if str_formatter is None:
        str_formatter = '[%(levelname)s] %(asctime)-15s [%(name)s %(funcName)s %(lineno)d] %(message)s'

    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    alog.addHandler(handler)
    alog.setLevel(level)


def _log_warn_once_per_movie():
    """If the zScore isn't found, only emit a single message per Movie"""
    movies = set()

    def _log_once(movie, message):
        if movie not in movies:
            log.warn(message)
            movies.add(movie)

    return _log_once

log_warn_once_per_movie = _log_warn_once_per_movie()


def _validate_resource(func, resource):
    """Validate the existence of a file/dir"""
    if func(resource):
        return os.path.abspath(resource)
    else:
        raise IOError("Unable to find {f}".format(f=resource))

validate_file = functools.partial(_validate_resource, os.path.isfile)
validate_dir = functools.partial(_validate_resource, os.path.isdir)
validate_output_dir = functools.partial(_validate_resource, os.path.isdir)


def _nfs_exists_check(ff):
    """Return whether a file or a dir ff exists or not.
    Call ls instead of python os.path.exists to eliminate NFS errors.
    """
    # this is taken from Yuan
    cmd = "ls %s" % ff
    _, rcode, _ = backticks(cmd)
    return rcode == 0


def fofn_to_files(fofn):
    """Util func to convert a bas/bax fofn file to a list of bas/bax files."""
    if os.path.exists(fofn):
        with open(fofn, 'r') as f:
            bas_files = {line.strip() for line in f.readlines()}

        for bas_file in bas_files:
            if not os.path.isfile(bas_file):
                # try one more time to find the file by
                # performing an NFS refresh
                found = _nfs_exists_check(bas_file)
                if not found:
                    raise IOError(
                        "Unable to find bas/bax file '{f}'".format(f=bas_file))

        return list(bas_files)
    else:
        raise IOError("Unable to find FOFN {f}".format(f=fofn))


def validate_fofn(fofn):
    """Validate existence of FOFN and files within the FOFN.

    :param fofn: (str) Path to File of file names.
    :raises: IOError if any file is not found.
    :return: (str) abspath of the input fofn

   :rtype: str
    """
    if os.path.isfile(fofn):
        file_names = fofn_to_files(os.path.abspath(fofn))
        log.debug("Found {n} files in FOFN {f}.".format(
            n=len(file_names), f=fofn))
        return os.path.abspath(fofn)
    else:
        raise IOError("Unable to find {f}".format(f=fofn))


def _is_cmp_h5(file_name):
    return file_name.lower().endswith('.cmp.h5')


def validate_file_or_none(path):
    if path is None:
        return path
    else:
        return validate_file(path)


def get_parser():
    desc = "Compare cmp.h5 to original bas/bax files."
    p = argparse.ArgumentParser(version=__version__, description=desc)
    p.add_argument('cmp_h5', type=validate_file,
                   help="Path to compare file CMP.H5 alignment file.")
    p.add_argument('output_csv', type=str, help="Output CSV file path.")
    p.add_argument('--fofn', help="FOFN of bas|bax files.",
                   required=True, type=validate_fofn)
    p.add_argument('--external', action='store_true',
                   help="Only output PacBio-external metrics")
    p.add_argument('--debug', action='store_true',
                   help="Send debug output to stdout")

    return p


class MovieParseException(Exception):
    pass


class ReadStats(object):

    def __init__(self):
        self.lengths = []
        self.zs = []
        self.accs = []
        self.n = 0
        self.reportAccuracy = False

    def __repr__(self):
        _d = dict(k=self.__class__.__name__,
                  n=self.n,
                  z=self.zscore,
                  a=self.accuracy,
                  l=self.length)
        return "<{k} n:{n} zscore:{z} accuracy:{a} length:{l} >".format(**_d)

    @property
    def zscore(self):
        if self.zs:
            return np.median(np.array(self.zs))
        else:
            return 0.0

    @property
    def length(self):
        if self.lengths:
            return np.median(np.array(self.lengths))
        else:
            return 0.0

    @property
    def accuracy(self):
        if self.accs:
            return np.median(np.array(self.accs))
        else:
            return 0.0

    def add(self, hit):
        l = abs(int(hit.attrib['end']) - int(hit.attrib['start']))
        self.lengths.append(l)
        z = float(hit.find('zScore').attrib['value'])
        self.zs.append(z)
        acc = float(hit.find('nCorrect').attrib['percent'])
        self.accs.append(acc)
        self.n += 1

    def addCmpAlnHit(self, hit):
        l = abs(hit.query_end - hit.query_start)
        self.lengths.append(l)
        z = hit.zScore
        if z is None:
            z = -1.0
        self.zs.append(z)
        nCorrect = l - hit.nMismatch - hit.nIns - hit.nDel
        acc = nCorrect / float(l) * 100.0
        self.accs.append(acc)
        self.n += 1

    def _quantile(self, v, quantile):
        if len(v) == 0:
            return 0.0
        n = len(v)
        nq = int(round(quantile * float(n)))
        if nq > n - 1:
            nq = n - 1
        v.sort()
        return v[nq]

    def _accFromZ(self, z):
        """derive accuracy from Z"""
        a = Z_S * z * math.sqrt(Z_A * Z_A + 1.0)
        y = math.exp(a - Z_C - Z_A / (FIXED_LENGTH + 20.0))
        acc = y / (1.0 + y)
        return acc

    def tostring(self, external=False):
        if not external:
            if not self.reportAccuracy:
                if self.n == 0:
                    return '0,0.0,0.0,0'
                return '%d,%.2f,%.2f,%.0f' % (self.n, self.zscore, self.accuray, self.length)

            if self.n == 0:
                return '0,0.0,0.0,0.0,0.0,0.0,0'
            zs = np.array(self.zs)
            z50 = np.median(zs)
            z95 = self._quantile(zs, 0.95)
            return '%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.0f' % (self.n, z50, 100.0 * self._accFromZ(z50), z95, 100.0 * self._accFromZ(z95), self.accuracy, self.length)
        else:
            if self.n == 0:
                return '0,0.00,0.0'
            meanAcc = np.mean(np.array(self.accs))
            meanRl = np.mean(np.array(self.lengths))
            return '%d,%.2f,%.1f' % (self.n, meanAcc, meanRl)


class MovieStats(object):

    def __init__(self, expt, chip, movie, inst, movieType='', setId='', partId='', cellId='', date='', time=''):
        """
        Container for all the movie statistics

        Note This is not quite consistent with the spec.

        :param expt: The 1234 of the internal LIMS run code 1234-123456
        :param chip: The 123456 of the internal LIMS run code 1234-123456
        :param movie: The movie name
        :param inst: The instrument
        :param movieType: Not sure what this is
        :param setId: 's1' (SET NUMBER)
        :param partId: 'p0' or x0 (Expired barcode? with the 0 (Strobe Number. always 0 now)
        :param date: the DATE part of m{DATE_TIME}_{INSTRUMENT}...
        :param time: the TIME part of m{DATE_TIME}_{INSTRUMENT}...

        """

        self.expt = expt
        if not self.expt:
            self.expt = '0000000'
        self.chip = chip
        if not self.chip:
            self.chip = '0000'

        self.movie = movie
        self.inst = inst
        self.movieType = movieType if movieType else 'NA'
        self.allStats = ReadStats()
        self.highZStats = ReadStats()
        self.highZStats.reportAccuracy = True
        self.nHighQVs = [0 for _ in xrange(len(QV_THRESHOLDS))]
        self.code = 'POSTALIGNMENT'
        self.setId = setId if setId else 'NA'
        self.partId = partId if partId else 'NA'
        self.cellId = cellId if cellId else 'NA'
        self.date = date if date else 'NA'

    def __repr__(self):
        _d = dict(k=self.__class__.__name__,
                  n=self.movie,
                  i=self.inst,
                  s=self.setId,
                  p=self.partId,
                  c=self.chip,
                  d=self.date,
                  t=self.movieType)
        return "<{k} {n} chip:{c} instrument:{i} set:{s} part:{p} type:{t} date:{d} >".format(**_d)

    def add(self, hit):
        if isinstance(hit, CmpH5Alignment):
            self._addCmpAlnHit(hit)
            return

        # should never get here.
        z = float(hit.find('zScore').attrib['value'])
        self.allStats.add(hit)
        if z > Z_THRESHOLD:
            self.highZStats.add(hit)

    def _addCmpAlnHit(self, hit):
        """

        :type hit: CmpH5Alignment
        """
        # pbcore supports zScore?
        try:
            z = hit.zScore
            self.allStats.addCmpAlnHit(hit)
            if z > Z_THRESHOLD:
                self.highZStats.addCmpAlnHit(hit)
        except ValueError:
            # probably want to disable this
            movie = hit.movieInfo[1]
            msg = "unable to get zScore from alignment ID {a} from movie".format(
                a=hit.AlnID, m=movie)
            log_warn_once_per_movie(movie, msg)

        # if not hit.hasPulseInfo:
        #    return

        for qv in hit.QualityValue():
            qv = qv / 10.0
            for i, qvThreshold in enumerate(QV_THRESHOLDS):
                if qv >= qvThreshold:
                    self.nHighQVs[i] += 1

    def tostring(self, external=False):
        """date,movie,runCode,expt,chip,inst,movieType,nReads,medianZ,
        medianAcc,medianLength,nReadsAboveZ%(z).0f,medianZaboveZ%(z).0f,
        AZ50,Z95,AZ95,medianAccAboveZ%(z).0f,medianLengthAboveZ%(z).0f,
        %(q)s,CODE' % { 'z':Z_THRESHOLD, 'q':qvHeader }
         or
        date,movie,inst,cellId,setId,strobeId,nReads,meanAcc,meanLength,
        %(q)s,CODE"""
        if not external:
            if not self.movie:
                return 'NA,NA,NA,NA,NA,NA,NA,0,0,%s,0' % NULL_QV_STRING
        else:
            if not self.movie:
                return 'NA,NA,NA,NA,NA,NA,0,0.0,0,%s,0' % NULL_QV_STRING
        qvString = ','.join([str(c) for c in self.nHighQVs])
        if not external:
            date = self.movie.split('_')[0][1:]
            runCode = '%s-%s' % (self.expt, self.chip)
            return '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' % \
                   (date, self.movie, runCode, self.expt,
                    self.chip, self.inst, self.movieType,
                    self.allStats.tostring(external),
                    self.highZStats.tostring(external),
                    qvString, self.code)
        else:
            return '%s,%s,%s,%s,%s,%s,%s,%s,%s' % \
                   (self.date, self.movie, self.inst, self.cellId,
                    self.setId, self.partId, self.allStats.tostring(external),
                    qvString, self.code)


def _parse_read_name(readName):
    search = NAME_PARSER.search(readName)
    if not search:
        msg = "Couldn't process read %s" % readName
        log.warn(msg)
        sys.stderr.write(msg + "\n")
        return 0, 0, None, None, None, None, None

    x = search.group(1)
    y = search.group(2)
    expt = search.group(3)
    chip = search.group(4)
    movie = search.group(5)

    # trim any pseudoread suffix from the movie
    movie = PSEUDO_PARSER.sub('', movie)
    search = MOVIE_PARSER.search(movie)

    if not search:
        msg = "Couldn't process read %s" % readName
        log.warn(msg)
        sys.stderr.write(msg + "\n")
        return 0, 0, None, None, None, None, None

    inst = search.group(1)
    movieType = search.group(2)

    return x, y, expt, chip, inst, movieType, movie


def _get_run_code_from_path(path):
    """Parse a PacBio internal path which gives expt and run IDs

    Returns a tuple of (x, y) or (None, None) if the value isn't found
    """
    search = EXPT_RUN_PARSER.search(path)
    if not search:
        return None, None

    return search.group(1), search.group(2)


def _get_post_mapping_from_movies(allMovies, cmp_h5):
    """
    Go through all movies post alignment.


    returns dict of {movie:MovieStats}
    """
    postMappingMovies = {}

    reader = CmpH5Reader(cmp_h5)

    for alignment in reader:
        # returns a tuple of
        # (2, 'm101210_151341_Jan_p1_b15', 100.0, 0.009999999776482582)
        movie_info = alignment.movieInfo

        movie = movie_info[1]
        if movie not in postMappingMovies:
            stats = allMovies[movie]
            postMappingMovies[movie] = MovieStats(stats.expt,
                                                  stats.chip, stats.movie,
                                                  stats.inst,
                                                  movieType=stats.movieType,
                                                  setId=stats.setId,
                                                  partId=stats.partId,
                                                  cellId=stats.cellId,
                                                  date=stats.date)
        postMappingMovies[movie].add(alignment)

    reader.close()

    return postMappingMovies


def _get_movie_stats_from_movie_files(movie_files):
    """
    :param movie_files: List of Movies

    Returns a dict of {movie_name: MovieStats}
    """

    allMovies = {}

    for movie_file in movie_files:
        log.info("Getting movie stats from {f}".format(f=movie_file))
        movie_stat = _movie_file_name_to_movie_stat(movie_file)
        log.info(movie_stat)
        allMovies[movie_stat.movie] = movie_stat

    return allMovies


def _get_internal_csv_header(z_threshold, qv_header):
    x = 'date,movie,runCode,expt,chip,inst,movieType,nReads,medianZ,medianAcc,medianLength,nReadsAboveZ%(z).0f,medianZaboveZ%(z).0f,AZ50,Z95,AZ95,medianAccAboveZ%(z).0f,medianLengthAboveZ%(''z).0f,%(q)s,CODE' % {
        'z': z_threshold, 'q': qv_header}
    return x


def _get_external_csv_header(qv_header):
    x = 'date,movie,inst,cellId,setId,strobeId,nReads,meanAcc,' \
        'meanLength,%(q)s,CODE' % {
            'q': qv_header}
    return x


def __movie_bax_name_to_movie_stat(file_name):
    """Parse a multi-part BAX file.

    Return a MovieStat instance

    :raises: MovieParseException
    """
    exp, runcode = _get_run_code_from_path(file_name)

    match = MOVIE_NAME3.search(file_name)
    if match:
        movie_name = match.groups()[0]
        search = MOVIE_PARSER3.search(file_name)

        if search:
            g = search.groups()
            date_time = g[0]
            date, dtime = date_time.split('_')
            inst = g[1]
            # Keeping backward compatibility. This is really the chip number
            cellId = 'c' + g[2]

            cellNumber = g[3]
            #
            setId = 's' + g[4]

            # p0 or X0
            partId = g[5] + '0'

            # bax.(1|2|3)
            movie_part_id = g[6]
            m = MovieStats(exp, runcode, movie_name, inst, cellId=cellId,
                           setId=setId, partId=partId, date=date, time=dtime)
            return m
        else:
            _d = dict(f=file_name, p=MOVIE_PARSER3.pattern)
            raise MovieParseException(
                "Unable to parse {f} with {p}".format(**_d))

    else:
        _d = dict(f=file_name, p=MOVIE_PARSER3.pattern)
        raise MovieParseException("Unable to parse {f} with {p}".format(**_d))


def __movie_bas_name_to_movie_stat(file_name):
    """Parse an old-style PRE multi-part movie file

    Return a MovieStat instance

    :raises: MovieParseException
    """
    exp, runcode = _get_run_code_from_path(file_name)
    match = MOVIE_NAME2.search(file_name)
    if match:
        movie_name = match.groups()[0]
        search = MOVIE_PARSER2.search(file_name)

        if search:
            g = search.groups()
            date_time = g[0]
            date, dtime = date_time.split('_')
            inst = g[1]
            # Keeping backward compatibility. This is really the chip number
            cellId = 'c' + g[2]

            cellNumber = g[3]
            #
            setId = 's' + g[4]

            # p0 or X0
            partId = g[5] + '0'

            m = MovieStats(exp, runcode, movie_name, inst, cellId=cellId,
                           setId=setId, partId=partId, date=date, time=dtime)
            return m
        else:
            _d = dict(f=file_name, p=MOVIE_PARSER2.pattern)
            raise MovieParseException(
                "Unable to parse {f} with {p}".format(**_d))

    else:
        _d = dict(f=file_name, p=MOVIE_NAME2.pattern)
        raise MovieParseException("Unable to parse {f} with {p}".format(**_d))


def __movie_astro_name_to_movie_stat(file_name):
    exp, runcode = _get_run_code_from_path(file_name)
    match = MOVIE_ASTRO_NAME.search(file_name)
    if match:
        movie_name = match.groups()[0]
        search = MOVIE_ASTRO_PARSER.search(file_name)
        if search:
            gs = search.groups()
            # Not quite sure how to handle this from a backward compatibility
            date_time = gs[0]
            # Not quite sure how to handle this from a backward compatibility
            cellId = 'c00000000000000000000000000000000'
            inst = "NA"
            cellNumber = gs[1]
            partId = "p1"
            date, dtime = date_time.split('_')

            m = MovieStats(exp, runcode, movie_name, inst,
                           cellId=cellId, setId="s1", date=date, time=dtime)
            return m

        else:
            _d = dict(f=file_name, p=MOVIE_ASTRO_NAME.pattern)
            raise MovieParseException(
                "Unable to parse {f} with {p}".format(**_d))
    else:
        _d = dict(f=file_name, p=MOVIE_ASTRO_NAME.pattern)
        raise MovieParseException("Unable to parse {f} with {p}".format(**_d))


def _movie_file_name_to_movie_stat(file_name):
    """
    Convert the a bas/bax/pls to a MovieStat instance.

    :type file_name: str

    :rtype: MovieStats
    """

    # Try Multi-part, fall back to old pre-bax style, then to Astro-style

    try:
        movie_stat = __movie_bax_name_to_movie_stat(file_name)
        return movie_stat
    except MovieParseException:
        try:
            movie_stat = __movie_bas_name_to_movie_stat(file_name)
            return movie_stat
        except MovieParseException:
            movie_stat = __movie_astro_name_to_movie_stat(file_name)
            return movie_stat


def _log_summary(movie_stats):
    for movie_stat in movie_stats:
        log.info(movie_stat)
        # log.info(movie_stat.allStats)
        # log.info(movie_stat.highZStats)


def run(cmp_h5, movie_files, output_csv, external_mode=False):
    """
    Run the analysis and create the output summary as a CSV

    :param cmp_h5: Path to .cmp.h5 alignment file
    :param movie_files: List of bax/bas files
    :param output_csv: Path to output CSV file
    :param external_mode: Used for internal or external metric display

    :type cmp_h5: str
    :type movie_files: list
    :type output_csv: str
    :type external_mode: bool

    :rtype: int
    """

    allMovies = _get_movie_stats_from_movie_files(movie_files)
    log.info(pformat(allMovies, indent=4))

    postMappingMovies = _get_post_mapping_from_movies(allMovies, cmp_h5)

    _log_summary(postMappingMovies.values())

    if len(allMovies) > 0:
        allMovieNames = allMovies.keys()
    else:
        allMovieNames = postMappingMovies.keys()

    allMovieNames.sort()

    qvHeader = ','.join(['nQVs>=%d' % q for q in QV_THRESHOLDS])

    with open(output_csv, 'w') as f:
        log.info("Writing output to CSV file {f}".format(f=output_csv))
        if not external_mode:
            f.write(_get_internal_csv_header(Z_THRESHOLD, qvHeader) + "\n")
        else:
            f.write(_get_external_csv_header(qvHeader) + "\n")

        for movie in allMovieNames:
            if movie in postMappingMovies:
                f.write(postMappingMovies[movie].tostring(
                    external_mode) + "\n")
            else:
                stats = allMovies[movie]
                stats.code = 'PREFILTER'
                f.write(stats.tostring(external_mode) + "\n")

    log.info("Completed writing CSV to {f}".format(f=output_csv))
    return 0


def main(argv=sys.argv):
    """Main point of entry"""
    p = get_parser()
    # the first arg with be the exe name
    args = p.parse_args(argv[1:])

    fofn = args.fofn
    cmp_h5 = args.cmp_h5
    movie_files = fofn_to_files(fofn)
    output_csv = args.output_csv
    external_mode = args.external

    if args.debug:
        setup_log(log, level=logging.DEBUG)
    else:
        log.addHandler(logging.NullHandler())

    log.debug(args)
    started_at = time.time()
    try:
        rcode = run(cmp_h5, movie_files, output_csv,
                    external_mode=external_mode)
    except Exception as e:
        rcode = -1
        log.error(e, exc_info=True)
        sys.stderr.write(str(e) + "\n")

    run_time = time.time() - started_at

    _d = dict(f=os.path.basename(__file__), x=rcode, s=run_time, v=__version__)
    log.info(
        "completed {f} v{v} with exit code {x} in {s:.2f} sec.".format(**_d))
    return rcode
