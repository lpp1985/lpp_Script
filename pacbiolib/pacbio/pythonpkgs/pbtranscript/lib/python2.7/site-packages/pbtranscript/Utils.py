"""Define util functions."""
import os
import os.path as op
import shutil
import logging
import sys
from time import sleep

from pbcore.io import openDataSet, ContigSet
from pbcore.util.Process import backticks


def revcmp(seq):
    """Given a sequence return its reverse complement sequence."""
    NTMAP = {'a': 't', 'c': 'g', 't': 'a', 'g': 'c',
             'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    return "".join([NTMAP[x] for x in seq])[::-1]


def realpath(f):
    """Return absolute, user expanded path."""
    if f is None:
        return None
    return op.abspath(op.expanduser(f))


def real_ppath(fn):
    """Return real 'python-style' path of a file.
    Consider files with white spaces in their paths, such as
    'res\ with\ space/out.sam' or 'res with space/out.sam',
    'res\ with\ space/out.sam' is unix-style file path.
    'res with space/out.sam' is python style file path.

    We need to convert all '\_' in path to ' ' so that python
    can handle files with space correctly, which means that
    'res\ with\ space/out.sam' will be converted to
    'res with space/out.sam'.
    """
    if fn is None:
        return None
    return op.abspath(op.expanduser(fn)).replace(r'\ ', ' ')


def real_upath(fn):
    """Return real 'unix-style' path of a file.
    Consider files with white spaces in their paths, such as
    'res\ with\ space/out.sam' or 'res with space/out.sam',
    'res\ with\ space/out.sam' is unix-style file path.
    'res with space/out.sam' is python style file path.

    We need to convert all ' ' to '\ ' so that unix can handle
    files with space correctly, which means that
    'res with space/out.sam' will be converted to
    'res\ with\ space/out.sam'.
    """
    if fn is None:
        return None
    return real_ppath(fn).replace(' ', r'\ ')

def nfs_exists(fn):
    """Detect whether a NFS file or a directory exists or not.
    In rare cases, a file f, which is created by a node X,
    may not be detected by node Y immediately due to NFS latency.
    In even more rare cases, node Y's shell script may be able
    to detect f's existence while Y's python can not, because of
    a python file system cache latency. So in order to eliminate
    this problem, first call 'ls' to trigger mount, then detect
    existence of fn using python.open(fn, 'r'), and try twice
    before finally give up.

    This script should return True, if fn is either
    an existing file or an existing directory created before
    nfs_exists is called. However, this script may return an
    arbitrary value, if fn is created or deleted at the same
    time or after nfs_exists is called.
    """
    # Call ls just to trigger mount, don't trust the return value.
    _o, _c, _m = backticks("ls {f}".format(f=fn))

    ERROR_NO_SUCH_FILE_OR_DIRECTORY = 2
    ERROR_IS_DIRECTORY = 21
    # Try to open fn with read-only mode
    try:
        with open(fn, 'r') as _reader:
            pass
        return True # fn is an existing file
    except IOError as err:
        if err.errno == ERROR_NO_SUCH_FILE_OR_DIRECTORY:
            # Wait 15 seconds and try it again
            sleep(15)
            try:
                with open(fn, 'r') as _reader:
                    pass
                return True # An newly detected existing file
            except IOError as err2:
                if err2.errno == ERROR_NO_SUCH_FILE_OR_DIRECTORY:
                    return False # fn does not exist
                elif err2.errno == ERROR_IS_DIRECTORY:
                    return True # fn is an existing directory
        elif err.errno == ERROR_IS_DIRECTORY:
            return True # fn is an existing directory
        else:
            return False # other IOErrors
    return False # other errors


def execute(cmd, errmsg="", errcls=RuntimeError):
    """Execute command and check exit code. If exit code is not
    0, raise RuntimeError with errmsg.
    """
    logging.debug("CMD: " + cmd)
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        msgs = [msg for msg in ["CMD failed: %s" % cmd, _msg, errmsg]
                    if len(msg) > 0]
        raise errcls("\n".join(msgs))


def mkdir(path):
    """Create a directory if it does not pre-exist,
    otherwise, pass."""
    if not op.exists(path) and not op.lexists(path):
        try:
            os.makedirs(path)
        except OSError as e:
            # "File exists error" can happen when python
            # fails to syncronize with NFS or multiple
            # processes are trying to make the same dir.
            if e.errno == 17:
                pass
            else:
                raise OSError(e)


def rmpath(path):
    """Remove a file or a directory"""
    if op.exists(path):
        execute("chmod +w -R %s && rm -rf %s" % (path, path))


def mknewdir(path):
    """Create a new directory if it does not pre-exist,
    otherwise, delete it and then re-create it."""
    rmpath(path)
    #shutil.rmtree(path)
    os.makedirs(path)


def touch(path):
    """touch a file."""
    if op.exists(path):
        os.utime(path, None)
    else:
        open(path, 'a').close()


def generateChunkedFN(out_dir, prefix, num_chunks):
    """Generate n chunked file names, e.g.
    outDir/$prefix.0, outDir/$prefix.1, ..., outDir/$prefix.num_chunks-1
    """
    return [op.join(out_dir, prefix + "." + str(i))
            for i in xrange(0, num_chunks)]


def get_files_from_fofn(fofn_filename):
    """Return a list of file names within a fofn file."""
    fns = []
    try:
        with open(fofn_filename, 'r') as fofn:
            for line in fofn:
                fns.append(realpath(line.strip()))
    except (IOError, OSError) as e:
        raise IOError("Failed to read from fofn file {fofn}.\n".
                      format(fofn=fofn_filename) + str(e))
    return fns


def write_files_to_fofn(file_names, fofn_filename):
    """Write files in list file_names to fofn_filename."""
    try:
        with open(fofn_filename, 'w') as fofn:
            for fn in file_names:
                fofn.write(str(fn) + "\n")
    except (IOError, OSError) as e:
        raise IOError("Failed to files to fofn file {fofn}.\n".
                      format(fofn=fofn_filename) + str(e))


def get_files_from_file_or_fofn(filename):
    """Return a list of file names in a file or a fofn."""
    if filename.endswith(".fofn"):
        return get_files_from_fofn(filename)
    else:
        return [filename]


def enum(**enums):
    """Simulate enum."""
    return type('Enum', (), enums)


FILE_FORMATS = enum(H5="H5", BAM="BAM", FASTA="FASTA", UNKNOWN="UNKNOWN")


def guess_file_format(filenames):
    """Given a list of files, return 'BAM' if all files ends with .bam;
    return 'H5' if all files ends with .h5; reutrn 'UNKNOWN' otherwise."""
    if isinstance(filenames, str):
        filenames = [filenames]
    elif not isinstance(filenames, list):
        raise TypeError("guess_file_format does not support %s" %
                        type(filenames))
    fns = []
    for fn in filenames:
        fns.extend(get_files_from_file_or_fofn(fn))

    if all([fn.endswith(".h5") for fn in fns]):
        return FILE_FORMATS.H5
    elif all([fn.endswith(".bam") for fn in fns]):
        return FILE_FORMATS.BAM
    elif all([fn.endswith(".xml") for fn in fns]):
        return FILE_FORMATS.BAM
    else:
        return FILE_FORMATS.UNKNOWN


def validate_fofn(fofn_filename):
    """Validate existence of FOFN and files within the FOFN.

    :param fofn: (str) Path to File of file names or None.
    :raises: IOError if any file is not found.
    :return: (str) input fofn or None
    """
    if fofn_filename is None:
        return None

    if nfs_exists(fofn_filename):
        if fofn_filename.endswith(".xml") or fofn_filename.endswith(".bam"):
            ds = openDataSet(fofn_filename, strict=True)
            return fofn_filename
        fns = get_files_from_fofn(fofn_filename)
        for fn in fns:
            if not nfs_exists(fn):
                raise IOError("Unable to find {f} in FOFN {fofn}.".
                              format(f=fn, fofn=fofn_filename))
        return fofn_filename
    else:
        raise IOError("Unable to find FOFN {fofn}.".
                      format(fofn=fofn_filename))


def setup_log(alog, file_name=None, level=logging.DEBUG, str_formatter=None):
    """
    Copied from mkocher's pbreports/utils.py.
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
        str_formatter = '[%(levelname)s] %(asctime)-15s ' + \
                        '[%(name)s %(funcName)s %(lineno)d] %(message)s'

    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    alog.addHandler(handler)
    alog.setLevel(level)


def now_str():
    """Return string of current time."""
    import datetime
    return str(datetime.datetime.now()).split(".")[0]


def phred_to_qv(phred):
    """Phred value to quality value."""
    return 10 ** -(phred / 10.0)


def cat_files(src, dst):
    """Concatenate files in src and save to dst.
       src --- source file names in a list
       dst --- destinate file name
    """
    if src is None or len(src) == 0:
        raise ValueError("src should contain at least one file.")
    if dst in src:
        raise IOError("Unable to cat a file and save to itself.")

    with open(real_ppath(dst), 'w') as writer:
        for src_f in src:
            with open(real_ppath(src_f), 'r') as reader:
                for line in reader:
                    writer.write(line.rstrip() + '\n')


def get_all_files_in_dir(dir_path, extension=None):
    """return all files in a directory."""
    fs = []
    for f in os.listdir(dir_path):
        if extension is None:
            fs.append(f)
        else:
            if f.endswith(extension):
                fs.append(f)
    return fs


class CIGAR(object):

    """Cigar string."""

    def __init__(self, cigar_str):
        self.cigar_str = cigar_str
        self.num_match = 0
        self.num_mismatch = 0
        self.num_insert = 0
        self.num_deletion = 0
        self.num_hardclip = 0
        self.num_softclip = 0
        self.num_padding = 0
        self.num_unknown = 0
        self.parse(self.cigar_str)

    def parse(self, cigar_str):
        """Parse cigar string."""
        s = cigar_str
        if s == "*" or s == "=":
            return
        while(len(s) > 0):
            i = 0
            while s[i].isdigit() and i < len(s):
                i += 1
            num = int(s[:i])
            action = s[i]
            s = s[i + 1:]
            if action == "M":
                self.num_match += num
            elif action == "X":
                self.num_mismatch += num
            elif action == "I":
                self.num_insert += num
            elif action == "D":
                self.num_deletion += num
            elif action == "H":
                self.num_hardclip += num
            elif action == "S":
                self.num_softclip += num
            elif action == "P":
                self.num_padding += num
            elif action == "N":
                self.num_unknown += num
            else:
                raise ValueError("Can not parse CIGAR string " +
                                 "{s}".format(s=cigar_str))
# using regular expression is 20% slower than naive way
#        pattern = r"^(\d+)(M|I|D|X|S|H|P|N)(.*)$"
#        while(len(s) > 0):
#            m = re.search(pattern, s)
#            if m:
#                num = int(m.groups()[0])
#                action = m.groups()[1]

    def __str__(self):
        return "Match = {m}, ".format(m=self.num_match) + \
               "Mismatch = {m}, ".format(m=self.num_mismatch) + \
               "Insert = {m}, ".format(m=self.num_insert) + \
               "Deletion = {m}, ".format(m=self.num_deletion) + \
               "HardClipping = {m}, ".format(m=self.num_hardclip) + \
               "SoftClipping = {m}, ".format(m=self.num_softclip) + \
               "Padding = {m}, ".format(m=self.num_padding)

    def match_seq(self, seq):
        """Return if this cigar string matches the sequence."""
        return self.cigar_str == "*" or self.cigar_str == "=" or \
            (self.num_match + self.num_insert + self.num_softclip == len(seq))


def cigar_match_seq(sam_str):
    """Return True if cigar length match sequence length, otherwise, False"""
    fields = sam_str.split('\t')
    cigar, seq = CIGAR(fields[5]), fields[9]
    if cigar.match_seq(seq):
        return True
    else:
        return False


def filter_sam(in_sam, out_sam):
    """Filter sam alignments with bad cigar string."""
    if not op.exists(in_sam):
        raise IOError("Unable to find input sam {f}".format(f=in_sam))
    if realpath(in_sam) == realpath(out_sam):
        raise IOError("in_sam and out_sam can not be identical.")
    with open(in_sam, 'r') as reader, \
            open(out_sam, 'w') as writer:
        for l, in_line in enumerate(reader):
            logging.debug("Processing {l}".format(l=in_line))
            if in_line.startswith("#") or \
               in_line.startswith("@") or \
               cigar_match_seq(in_line):
                writer.write(in_line)
            else:
                logging.warn("line {l}, cigar does not match sequence.".
                             format(l=l + 1))


def ln(src, dst):
    """if src and dst are identical, pass. Otherwise, create dst, a soft
    symbolic link pointing to src."""
    if realpath(src) != realpath(dst):
        if op.exists(dst) or op.lexists(dst):
            os.remove(dst)
        logging.debug("Creating a symbolic link {dst} pointing to {src}".
                      format(dst=dst, src=src))
        os.symlink(src, dst)


def mv(src, dst):
    """move src file to dst"""
    if realpath(src) != realpath(dst):
        execute("mv %s %s" % (src, dst))


def make_pbi(bam_fn, force_make=True):
    """Make *.bam.pbi from *.bam."""
    pbi_fn = "{bam}.pbi".format(bam=bam_fn)
    if not nfs_exists(pbi_fn) or force_make:
        logging.debug("create pbi index for {bam}.".format(bam=bam_fn))
        cmd = "pbindex {bam}".format(bam=bam_fn)
        _o, _c, _m = backticks(cmd)
        if _c != 0:
            raise RuntimeError("Failed to make pbi index for {bam}: {e}".
                               format(bam=bam_fn, e=_m))
    else:
        logging.debug("pbi index {pbi} already exists.".
                      format(pbi=pbi_fn))


def as_contigset(fasta_file, xml_file):
    if fasta_file == xml_file or xml_file is None:
        if not op.isfile(fasta_file) or op.getsize(fasta_file) == 0:
            return ContigSet()
        return ContigSet(fasta_file)
    file_size = op.getsize(fasta_file)

    fai_file = fasta_file + ".fai"
    if op.exists(fai_file):
        os.remove(fai_file)

    ds = ContigSet(fasta_file, generateIndices=True)
    ds.write(xml_file)
    if not file_size > 0:
        with open(fai_file, "w") as fai:
            fai.write("")
    return ds


def get_sample_name(input_sample_name):
    """
    Return sample name, only alphabet and digits are allowed. (no underscore, no space)
    If input_sample_name is None or empty string, return a random string
    starts with 'sample'
    """
    if input_sample_name is None or str(input_sample_name) == "None":
        input_sample_name = ""
    else:
        input_sample_name = str(input_sample_name)

    input_sample_name = "".join([x for x in input_sample_name
                                 if x.isalpha() or x.isdigit()])

    if len(input_sample_name) == 0:
        import binascii
        return str("sample"+binascii.b2a_hex(os.urandom(3)))
    else:
        return input_sample_name


def get_samtools_version():
    """
    Return samtools version string.
    If samtools supports '--version', return its version.
    Otherwise, return '0.1.19' which is the one used in SA2.x, SA3.0, SA3.1 and SA3.2.
    """
    default_sam_version = "0.1.19"
    cmd = 'samtools --version'
    _out, _code, _msg = backticks(cmd)
    try:
        if "samtools" in _out[0]:
            return str(_out[0][8:]).strip()
    except Exception:
        # samtools used in SA2.x and SA3.1, SA3.2 is
        # version 0.1.19, which does not support --version
        # samtools used in SA3.3 or up is 1.3.1
        return default_sam_version


def use_samtools_v_1_3_1():
    """Return True if samtools major version is >= 1.3.1"""
    try:
        versions = get_samtools_version().split('.')
        major, minor, last = 0, 0, 0
        major = int(versions[0])
        if len(versions) >= 2:
            minor = int(versions[1])
            if len(versions) >= 3:
                last = int(versions[2])
        return major >= 2 or (major == 1 and minor > 3) or (major == 1 and minor == 3 and last >= 1)
    except:
        return False
