
"""
Validate BAM alignment files against PacBio's internal spec.  It is assumed
that files can be read by the pbcore.io module, which will do some checking
implicitly.
"""

from __future__ import division, print_function
from collections import defaultdict
from functools import wraps
import argparse
import unittest
import hashlib
import logging
import os.path
import re
import sys

import pysam
from numpy import array

import pbcore.chemistry.chemistry
from pbcore.io.align._BamSupport import IncompatibleFile
import pbcore.io

from pbcoretools.pbvalidate.core import *

#log = logging.getLogger()


class Constants (object):
    PYSAM_VERSION = [int(x) for x in pysam.__version__.split(".")]
    # first some global definitions
    FLAG_SEGMENT_UNMAPPED = 4
    # XXX shouldn't these be in pbcore.io?
    ADAPTER_BEFORE = 1
    ADAPTER_AFTER = 2
    BARCODE_BEFORE = 4
    BARCODE_AFTER = 8
    FORWARD_PASS = 16
    REVERSE_PASS = 32
    DEFAULT_SN = [-1., -1., -1., -1.]
    # XXX supposedly the movie name should be lowercase, but this is not always
    # true - should it be?
    REGEX_QNAME = "([a-zA-Z0-9_]+)/([0-9]+)/((ccs)|([0-9]+_[0-9]+))$"
    REGEX_BASECALLS = "^[ACGTN]+$"
    REGEX_BAD_CIGAR = "M"
    # various tag keys and values
    CHEMISTRY_TAGS = ["BINDINGKIT", "SEQUENCINGKIT", "BASECALLERVERSION"]
    FRAMERATE_TAG = "FRAMERATEHZ"
    BASECALLER_TAG = "BASECALLERVERSION"
    READ_TYPE_TAG = "READTYPE"
    UNMAPPED_SORT = ["queryname", "unknown"]
    QSTART_QEND_TAGS = ["qs", "qe"]
    SC_TAG = "SC"
    SC_TAG_VALUES = ["A", "B", "L"]
    REQUIRED_ALIGNMENT_TAGS = ["zm", "np", "rq", "sn"]
    REQUIRED_BARCODE_TAGS = ["bc", "bq"]
    CODEC_NAMES = ["FRAMES", "CODECV1"]
    PULSE_FEATURE_KEYS = ["IPD", "PULSEWIDTH"]
    REQUIRES_CONTEXT_TAGS = ["SUBREAD", "CSS"]
    UNKNOWN_READ_TYPES = ["ZMW", "SCRAP"]
    EXPECTED_SUFFIX = {
        "standard": "subreads",  # really SUBREAD
        "CCS": "ccs",
        #"SUBREAD" : "subreads",
        "SCRAP": "scraps",
        "ZMW": "zmws",
    }
    # misc string constants
    READ_TYPE_SUBREAD = "SUBREAD"
    # XXX file-level readType.  the substitution of 'standard' for 'SUBREAD' in
    # pbcore.io seems really broken to me
    READ_TYPE_SUBREAD_FILE = "standard"  # XXX file-level readType
    READ_TYPE_CCS = "CCS"
    READ_TYPE_SCRAP = "SCRAP"
    READ_TYPE_MIXED = "mixed"
    SORT_ORDER_UNKNOWN = "unknown"


def _check_pysam_version():  # XXX unused, see below
    import pysam.version
    version = pysam.version.__version__.split(".")
    if (version < ["0", "8", "2"]):
        raise ImportError("pysam >= 0.8.2 required")

# workaround for lack of has_tag() in pysam < 0.8.2
# FIXME it would be nice if we could just require a newer version, but one of
# the other PacBio modules forces me to use 0.8.1


def _has_tag(aln, tag):
    if hasattr(aln, "has_tag"):
        return aln.has_tag(tag)
    try:
        value = aln.opt(tag)
    except KeyError as e:
        return False
    else:
        return True


class BAMError (ValidatorError):
    pass


class BAMReadError (RecordValidatorError):
    pass


class FileContentMismatchError (BAMError):
    MESSAGE_FORMAT = "File was expected to contain only %s reads, but this file contains %s reads"


class FileSuffixMissingError (BAMError):
    MESSAGE_FORMAT = "File contains %s reads, but lacks the expected suffix."


class MissingIndexError (BAMError):
    MESSAGE_FORMAT = "Missing corresponding .pbi index file"


class BadIndexError (BAMError):
    MESSAGE_FORMAT = "Couldn't access record %d directly (error: '%s'); " +\
        "this may indicate that the BAM file is truncated or the .pbi " +\
        "file is corrupt."


class FileNotAlignedError (BAMError):
    MESSAGE_FORMAT = "This file was expected to contain aligned reads, but " +\
        "does not specify any reference sequences."


class FileAlignedError (BAMError):
    MESSAGE_FORMAT = "This file was expected to contain unmapped reads, but " +\
        "a reference sequence is defined."


class UnsortedError (BAMError):
    MESSAGE_FORMAT = "This file has not been sorted by %s, or the " +\
        "header has not been updated."


class MixedReadUnsortedError (BAMError):  # XXX untested
    MESSAGE_FORMAT = "This file contains a mix of read types, but has " +\
        "not been sorted correctly."


class MixedReadSortingError (BAMError):  # XXX untested
    MESSAGE_FORMAT = "Improper sorting: subread %s follows read with type %s"


class BadSortingError (BAMError):  # XXX untested
    MESSAGE_FORMAT = "Improper sorting: either hole number or qStart for read %s " +\
        "is less than that of previous read %s"


class MissingPlatformError (BAMError):
    MESSAGE_FORMAT = "Missing platform (PL) for read group %s"


class WrongPlatformError (BAMError):
    MESSAGE_FORMAT = "Platform for read group %s is defined as '%s', not 'PACBIO'"

# XXX redundant with ReadGroupChemistryError - do we need both?


class BasecallerVersionError (BAMError):
    MESSAGE_FORMAT = "The basecaller version number '%s' in read group '%s' " +\
        "is not one of the allowed values; it is probably being misreported " +\
        "by an upstream program such as baz2bam or bax2bam."


# XXX basically redundant with FileContentMismatchError
class ReadTypeError (BAMError):
    MESSAGE_FORMAT = "Only '%s' reads are allowed, but read group %s has " +\
        "type '%s'"


class ReadGroupIdMismatchError (BAMError):
    MESSAGE_FORMAT = "Mismatch between specified and expected read group " +\
        "ID: %s in file, but computed as %s"


class ReadGroupChemistryError (BAMError):
    MESSAGE_FORMAT = "The chemistry information for read group %s is either " +\
        "missing or cannot be interpreted: %s"


class MissingCodecError (BAMError):
    MESSAGE_FORMAT = "The read group %s declares the tag '%s' for pulse " +\
        "features, but the encoding scheme is not specified or not recognized"


class AlignmentNotUniqueError (BAMReadError):
    MESSAGE_FORMAT = "QNAME %s occurs more than once (tStart = %s)"


class AlignmentUnmappedError (BAMReadError):
    MESSAGE_FORMAT = "The alignment of %s is not mapped to the reference " +\
        "sequence"


class ReadLengthError (BAMReadError):
    MESSAGE_FORMAT = "The length of the sequence of %s and the length " +\
        "indicated by qStart/qEnd disagree (%d versus %d)"


class MissingAlignmentTagError (BAMReadError):  # XXX untested
    MESSAGE_FORMAT = "Record for %s missing tag '%s' required for %s"

    def __init__(self, message, object_ref, tag_name):
        BAMReadError.__init__(self, message, object_ref)
        self.tag_name = tag_name

    def __hash__(self):
        """
        Here this incorporates the tag name, so in non-verbose mode every
        missing tag will be listed.
        """
        return hash(self.__class__.__name__ + "_" + self.tag_name)

    def test_name(self):
        return "test_" + self.__class__.__name__[:-5] + "_" + self.tag_name

    @classmethod
    def from_args(cls, object_ref, *args):
        tag_name = args[1]
        return cls(cls.MESSAGE_FORMAT % args, object_ref, tag_name=tag_name)


class MissingReadGroupError (BAMReadError):
    MESSAGE_FORMAT = "Can't find read group with ID '%s' for alignment of %s"


class TagValueError (MissingAlignmentTagError):
    MESSAGE_FORMAT = "Value '%s' of tag '%s' for qname %s is not one of the " +\
        "expected values or within the required range."


class BadQVsError (BAMReadError):
    MESSAGE_FORMAT = "The QV field '%s' for record %s contains entirely '!' " +\
        "characters; this will cause problems for mapping programs such as " +\
        "Blasr."


class QnameFormatError (BAMReadError):
    MESSAGE_FORMAT = "Query name '%s' is not in the expected format!"


class QnameMovieError (BAMReadError):
    MESSAGE_FORMAT = "Mismatch between query movie name and read group " +\
        "movie name: %s, %s"


class QnameHoleNumberError (BAMReadError):
    MESSAGE_FORMAT = "Mismatch between hole number in QNAME %s and ZM " +\
        "field (%s) in optional tags"


class QnameReadTypeError (BAMReadError):  # XXX untested
    MESSAGE_FORMAT = "Mismatch between read group type and query name format: %s, %s"


# XXX implicit in regex, remove?
class QnameRangeFormatError (BAMReadError):
    MESSAGE_FORMAT = "Can't convert '%s' in QNAME to integer range"


class QnameRangeError (BAMReadError):
    MESSAGE_FORMAT = "Range specified in QNAME %s conflicts with QS and QE tags"

# XXX redundant with pbcore.io and hence untested


class QnameRangeTagError (BAMReadError):
    MESSAGE_FORMAT = "QNAME %s specifies a range, but either QS or QE is undefined"


class NonNucleotideError (BAMReadError):  # XXX untested, requires .pbi
    MESSAGE_FORMAT = "SubstitutionTag or DeletionTag for alignment of %s has non-ACGTN characters"


class LocalContextFlagsError (BAMReadError):  # XXX untested
    MESSAGE_FORMAT = "The local context (CX tag) for alignment of %s contains " +\
        "at least one pair of mutually exclusive flags (integer value = %d)"


class PulseFeatureError (BAMReadError):
    MESSAGE_FORMAT = "Record %s is missing one or more pulse features " +\
        "declared in manifest: %s"


class UninitializedSNRError (BAMReadError):
    MESSAGE_FORMAT = "The HQRegionSNR field ('sn') for the alignment of %s " +\
        "is uninitialized"


class AlignmentCigarError (BAMReadError):
    pass


class AlignmentCigarMatchError (BAMReadError):
    MESSAGE_FORMAT = "The CIGAR string for the alignment of %s contains " +\
        "the ambiguous 'M' character to specify either a match or mismatch"


class UnmappedPropertiesError (BAMReadError):
    MESSAGE_FORMAT = "This file contains unmapped reads, but the alignment " +\
        "of %s has one or more conflicting properties"


class BadTranscriptError (BAMReadError):  # XXX untested
    MESSAGE_FORMAT = "Unexpected error attempting to extract transcript " +\
        "for %s: '%s'"

#
#
# Utility functions


def _get_key_value_pairs_dict(field):
    pairs = [item.split("=") for item in field.split(";")]
    return {k: v for (k, v) in pairs}


class ValidateFileName (ValidateFileObject):

    """
    Check whether the file name conforms to the convention expected given
    the contents.
    """

    def validate(self, file_obj):
        e = self._get_errors(file_obj)
        return e is None or (isinstance(e, list) and len(e) == 0)

    def _get_errors(self, file_obj):
        error = None
        file_name = os.path.basename(file_obj.filename)
        readType = file_obj.readType
        if readType == "unknown":
            # XXX for some reason pbcore.io treats all other read types as
            # "unknown", but the PacBio BAM spec does have rules for them too
            readTypes = file_obj.readGroupTable.ReadType
            for otherType in Constants.UNKNOWN_READ_TYPES:
                if all(readTypes == otherType):
                    readType = otherType
                    break
        for rtype, suffix in Constants.EXPECTED_SUFFIX.iteritems():
            if file_name.endswith(suffix + ".bam") and readType != rtype:
                error = FileContentMismatchError.from_args(file_obj,
                                                           rtype, readType)
            elif readType == rtype and not file_name.endswith(suffix + ".bam"):
                # XXX we may not actually want to enforce this...
                error = FileSuffixMissingError.from_args(file_obj, readType)
        if ".aligned_" in file_name and not file_obj.isMapped:
            error = FileNotAlignedError.from_args(file_obj)
        return error

    def to_errors(self, file_object):
        return self._get_errors(file_object)


class ValidateContents (ValidateFileName):

    """Enforce a particular content type, e.g. aligned subreads"""

    def __init__(self, aligned, content_type):
        self._aligned = aligned
        self._content_type = content_type

    def _get_errors(self, file_obj):
        errors = []
        if self._aligned and not file_obj.isMapped:
            errors.append(FileNotAlignedError.from_args(file_obj))
        # XXX explicitly checking for False, because None means no preference
        elif self._aligned == False and file_obj.isMapped:
            errors.append(FileAlignedError.from_args(file_obj))
        if self._content_type is not None:
            if ((file_obj.readType == Constants.READ_TYPE_SUBREAD_FILE and
                 self._content_type == Constants.READ_TYPE_SUBREAD) or
                    (file_obj.readType == self._content_type)):
                pass
            else:
                errors.append(FileContentMismatchError.from_args(file_obj,
                                                                 self._content_type, file_obj.readType))
        return errors


class ValidateIndex (ValidateFileName):
    """
    Check that an PacBio index file (.pbi) exists.
    """

    def _get_errors(self, file_obj):
        if os.path.exists(file_obj.filename + ".pbi"):
            return []
        else:
            return [MissingIndexError.from_args(file_obj)]


# TODO write test
class ValidateRandomAccess (ValidateRandomAccessBase):
    ERROR_CLASS = BadIndexError
    INDEX_ATTR = "pbi"


class ValidateSorting (ValidateFileName):

    """
    Enforce sorting of reads.  This first examines the @HD::SO tag, which
    should be explicitly set to either 'coordinate' or 'queryname'
    depending on whether reads are aligned or not, or 'unknown' if the
    file contains a mix of unaligned SUBREAD+CCS reads.  In the latter
    case, the sorting is also verified explicitly.
    """

    def __init__(self, aligned=False):
        self._expect_aligned = aligned

    def _get_errors(self, file_obj):
        # XXX should I also verify that mapped reads are sorted correctly?
        if file_obj.isMapped or self._expect_aligned:
            if not file_obj.isSorted:
                return [UnsortedError.from_args(file_obj, "position")]
        elif file_obj.readType == Constants.READ_TYPE_SUBREAD_FILE:
            if not file_obj.peer.header["HD"]["SO"] in Constants.UNMAPPED_SORT:
                return [UnsortedError.from_args(file_obj, "QNAME")]
        elif file_obj.readType == Constants.READ_TYPE_MIXED:
            if file_obj.peer.header["HD"]["SO"] != Constants.SORT_ORDER_UNKNOWN:
                return [MixedReadUnsortedError.from_args(file_obj)]
            last_aln = None
            for aln in file_obj:
                if last_aln is not None:
                    if (aln.readType == Constants.READ_TYPE_SUBREAD and
                            last_aln.readType != Constants.READ_TYPE_SUBREAD):
                        return [MixedReadSortingError.from_args(file_obj,
                                                                aln.readType, last_aln.readType)]
                    elif (aln.readType == Constants.READ_TYPE_CCS and
                          last_aln.readType == Constants.READ_TYPE_SUBREAD):
                        pass
                    elif last_aln.qName > aln.qName:
                        return [BadSortingError.from_args(file_obj, aln.qName,
                                                          last_aln.qName)]
                    last_aln = aln
        return []


class ValidateReadGroup (ValidateBase):

    """
    Base class for validation of @RG headers
    """

    def validate(self, rg):
        return len(self._get_errors(rg)) == 0

    def to_errors(self, rg):
        return self._get_errors(rg)


class ValidateReadGroupPlatform (ValidateReadGroup):

    """
    Check that the read group platform (PL tag) is set to 'PACBIO'
    """

    def _get_errors(self, rg):
        if (not "PL" in rg):
            return [MissingPlatformError.from_args(rg, rg['ID'])]
        elif (not "PACBIO" in rg["PL"]):  # XXX or == PACBIO?
            return [WrongPlatformError.from_args(rg, rg['ID'], rg['PL'])]
        return []


class ValidateReadGroupBasecaller (ValidateReadGroup):

    """
    Check that BASECALLERVERSION is one of the allowed values (as defined in
    pbcore.chemistry).
    """
    __VERSION_KEYS = pbcore.chemistry.chemistry._BARCODE_MAPPINGS.keys()
    VERSIONS = [k[-1] for k in __VERSION_KEYS]

    def _get_errors(self, rg):
        version = _get_basecaller_version(rg)
        if not version in self.VERSIONS:
            return [BasecallerVersionError.from_args(rg, version, rg['ID'])]
        return []


def _get_basecaller_version(rg):
    ds_dict = _get_key_value_pairs_dict(rg.get("DS", ""))
    version = ds_dict.get(Constants.BASECALLER_TAG, "None.None")
    return ".".join(version.split(".")[0:2])


def _rg_id_string(rg):
    ds_dict = _get_key_value_pairs_dict(rg.get("DS", ""))
    movieName = rg.get("PU")
    readType = ds_dict.get(Constants.READ_TYPE_TAG)
    return hashlib.md5(movieName + "//" + readType).hexdigest()[:8]


def _get_read_type(rg):
    ds_dict = _get_key_value_pairs_dict(rg.get("DS", ""))
    return ds_dict.get(Constants.READ_TYPE_TAG)


class ValidateReadGroupId (ValidateReadGroup):

    """
    Verify that re-computing the read group ID from the movieName and
    readType gives us the same ID as defined in the file.
    """

    def validate(self, rg):
        from pbcore.io.align._BamSupport import rgAsInt
        rgId = rg['ID']
        rgIdInt = rgAsInt(rgId)
        rgIdString = _rg_id_string(rg)
        rgIdInt_derived = rgAsInt(rgIdString)
        if (rgIdInt_derived != rgIdInt):
            return False
        return True

    def to_errors(self, rg):
        rgId = rg['ID']
        rgIdString = _rg_id_string(rg)
        return [ReadGroupIdMismatchError.from_args(rg, rgId, rgIdString)]


class ValidateReadGroupType (ValidateReadGroup):

    """Check that all read groups are of the requested type."""

    def __init__(self, read_type):
        self._read_type = read_type

    def _get_errors(self, rg):
        readType = _get_read_type(rg)
        if readType != self._read_type:
            return [ReadTypeError.from_args(rg, self._read_type, rg['ID'],
                                            readType)]
        return []


class ValidateReadGroupChemistry (ValidateReadGroup):

    """
    Check that the chemistry information in @RG is present and can be decoded.
    """

    def _get_errors(self, rg):
        ds_dict = _get_key_value_pairs_dict(rg.get("DS", ""))
        fields = []
        for tag in Constants.CHEMISTRY_TAGS:
            val = ".".join(ds_dict.get(tag, "None.None").split(".")[0:2])
            fields.append(val)
        if None in fields:
            return [ReadGroupChemistryError.from_args(rg, rg['ID'],
                                                      tuple(fields))]
        decoded = pbcore.chemistry.decodeTriple(*fields)
        if decoded == "unknown":
            return [ReadGroupChemistryError.from_args(rg, rg['ID'],
                                                      tuple(fields))]
        return []


class ValidateReadGroupPulseManifest (ValidateReadGroupChemistry):

    """Check whether the pulse features declared by the read group header
    specify the encoding used."""

    def _get_errors(self, rg):
        ds_dict = _get_key_value_pairs_dict(rg.get("DS", ""))
        for key in ds_dict.keys():
            key_ = key.upper()
            if key_ in Constants.PULSE_FEATURE_KEYS:
                return [MissingCodecError.from_args(rg, rg['ID'], key)]
            else:
                for pf in Constants.PULSE_FEATURE_KEYS:
                    if key_.startswith(pf):
                        ft, codec = key_.split(":")
                        if not codec in Constants.CODEC_NAMES:
                            return [MissingCodecError.from_args(rg, rg['ID'], key)]
        return []


class ValidateReadBase (ValidateRecord):

    """Base class for validation of individual alignment records"""

    def __init__(self, aligned=None, contents=None):
        self._aligned = aligned
        self._contents = contents

    def validate(self, aln):
        e = self._get_errors(aln)
        return e is None or (isinstance(e, list) and len(e) == 0)

    def to_errors(self, aln):
        return self._get_errors(aln)

    def to_error(self, aln):
        return self.to_errors()[0]


class ValidateReadUnique (ValidateReadBase):

    def __init__(self):
        self._tstarts = {}

    def _get_errors(self, aln):
        if aln.qName in self._tstarts:
            self._tstarts[aln.qName].append(aln.tStart)
            tstarts = ", ".join([str(x) for x in self._tstarts[aln.qName]])
            return [AlignmentNotUniqueError.from_args(aln, aln.qName, tstarts)]
        self._tstarts[aln.qName] = [aln.tStart]
        return []


class ValidateReadMapped (ValidateReadBase):

    """Confirm that an alignment is actually mapped to the reference (only
    suitable for some BAM files)"""

    def _get_errors(self, aln):
        # XXX only need to do validation if the file is mapped AND aligned
        # isn't explicitly set to False
        if aln.bam.isMapped and self._aligned != False:
            if ((aln.peer.flag & Constants.FLAG_SEGMENT_UNMAPPED != 0) or
                    # XXX this is a little confusing because the SAM and BAM
                    # conventions are different
                    aln.peer.pos < 0 or
                    aln.peer.rname == -1):
                return [AlignmentUnmappedError.from_args(aln, aln.qName)]
        return []


class ValidateReadUnmapped (ValidateReadBase):

    """
    For unmapped files, verify that an alignment has the expected properties -
    POS=0, RNAME=*, FLAG&4
    """

    def _get_errors(self, aln):
        # XXX like the mapped case, only need to check if the file and the
        # aligned parameter both say it's unmapped
        if not aln.bam.isMapped and not self._aligned:
            if (aln.peer.flag & Constants.FLAG_SEGMENT_UNMAPPED == 0 or
                    aln.peer.pos >= 0 or
                    aln.peer.rname != -1):  # '*'
                return [UnmappedPropertiesError.from_args(aln, aln.qName)]
        return []


class ValidateReadReadGroup (ValidateReadBase):

    """Make sure the read group for an alignment actually exists"""

    def _get_errors(self, aln):
        try:
            rg = aln.readGroupInfo
        except KeyError:
            return [MissingReadGroupError.from_args(aln, aln.peer.opt("RG"), aln.qName)]
        else:
            return []


def _requires_read_group_info(method):
    """
    Decorator function to check that the readGroupInfo property can be accessed
    before calling the decorated method.
    """
    @wraps(method)
    def f(self, aln, *args, **kwds):
        try:
            rg = aln.readGroupInfo
        except KeyError:
            return []
        else:
            return method(self, aln, *args, **kwds)
    return f


class ValidateReadQname (ValidateReadBase):

    """
    Check the format of the QNAME field, which should take the form::

        {movieName}/{holeNumber}/{qStart}_{qEnd}

    for unrolled reads and subreads, or::

        {movieName}/{holeNumber}/ccs

    for CCS reads.  The hole number, qStart, and qEnd must be integers.
    """

    def _get_errors(self, aln):
        # XXX the regex is convenient, but it may be checking too many things
        # at once - should we just use qName.split("/") and check each field
        # separately?
        fields = re.match(Constants.REGEX_QNAME, aln.qName)
        if not fields:  # don't even bother checking the rest
            return [QnameFormatError.from_args(aln, aln.qName)]

        @_requires_read_group_info
        def _get_errors(self_, aln_):
            rg = aln.readGroupInfo
            movieName = fields.group(1)
            holeNumber = fields.group(2)
            qrange = fields.group(3)  # either qStart_qEnd or 'ccs'
            err = []
            if aln.movieName != movieName:
                err.append(
                    QnameMovieError.from_args(aln, movieName, aln.movieName))
            # commented out since this is implicit in the regex
            # try :
            holeNumber = int(holeNumber)
            # except ValueError :
            #    err.append("Hole number '%s' is not an integer." % holeNumber)
            if holeNumber != aln.HoleNumber:
                err.append(
                    QnameHoleNumberError.from_args(aln, aln.qName, aln.HoleNumber))
            if ((rg.ReadType == Constants.READ_TYPE_CCS and qrange != "ccs") or
                    (rg.ReadType != Constants.READ_TYPE_CCS and qrange == "ccs")):
                err.append(
                    QnameReadTypeError.from_args(aln, rg.readType, qname))
            if qrange != "ccs":
                try:
                    qstart, qend = [int(x) for x in qrange.split("_")]
                except Exception:
                    err.append(QnameRangeFormatError.from_args(aln, aln.qName))
                else:
                    try:
                        if qstart != aln.qStart or qend != aln.qEnd:
                            err.append(
                                QnameRangeError.from_args(aln, aln.qName))
                    except Exception:
                        err.append(
                            QnameRangeTagError.from_args(aln, aln.qName))
            return err
        return _get_errors(self, aln)


class ValidateReadLength (ValidateReadBase):

    def _get_errors(self, aln):
        rg = aln.readGroupInfo
        if rg.ReadType != Constants.READ_TYPE_CCS:
            qlen = aln.qEnd - aln.qStart
            seq_len = len(aln.peer.seq)
            if seq_len != qlen:
                return [ReadLengthError.from_args(aln, aln.qName, seq_len, qlen)]
        return []


class ValidateReadTags (ValidateReadBase):

    """
    Check for mandatory tags and associated constraints in the 'optional'
    section.  The values may vary by read type, but the tags should be present
    in all cases.
    """

    @_requires_read_group_info
    def _get_errors(self, aln):
        rg = aln.readGroupInfo
        errors = []
        for tag in Constants.REQUIRED_ALIGNMENT_TAGS:
            if not _has_tag(aln.peer, tag):
                errors.append(
                    MissingAlignmentTagError.from_args(aln, aln.qName, tag, rg.ReadType))
            else:
                if tag == "np":
                    if rg.ReadType == Constants.READ_TYPE_SUBREAD and aln.numPasses != 1:
                        errors.append(
                            TagValueError.from_args(aln, aln.numPasses, "NP", aln.qName))
                elif tag == "rq":
                    # as of spec 3.0.1, this is a float between 0 and 1
                    if not 0 <= aln.peer.opt("rq") <= 1:
                        errors.append(
                            TagValueError.from_args(aln, aln.peer.opt("rq"), "rq",
                                                    aln.qName))
        return errors


class ValidateReadTagsMisc (ValidateReadTags):

    @_requires_read_group_info
    def _get_errors(self, aln):
        errors = []
        rg = aln.readGroupInfo
        if rg.readType == Constants.READ_TYPE_SUBREAD:
            if not _has_tag(aln.peer, "cx"):
                return [MissingAlignmentTagError.from_args(
                    aln, aln.qName, 'cx', Constants.READ_TYPE_SUBREAD)]
            localContext = aln.peer.opt("cx")
        if rg.ReadType != Constants.READ_TYPE_CCS:
            for tag in Constants.QSTART_QEND_TAGS:
                if not _has_tag(aln.peer, tag):
                    errors.append(MissingAlignmentTagError.from_args(aln,
                                                                     aln.qName, tag, rg.ReadType))
        if rg.ReadType == Constants.READ_TYPE_SCRAP:
            if not _has_tag(aln.peer, Constants.SC_TAG):
                errors.append(
                    MissingAlignmentTagError.from_args(aln, aln.qName,
                                                       Constants.SC_TAG, Constants.READ_TYPE_SCRAP))
            else:
                value = aln.peer.opt(Constants.SC_TAG)
                if not value in Constants.SC_TAG_VALUES:
                    errors.append(
                        TagValueError.from_args(aln, value, Constants.SC_TAG,
                                                aln.qName))
        return errors


class ValidateReadQVs (ValidateReadBase):

    def _get_errors(self, aln):
        e = []
        for tag in ['sq', 'dq', 'iq']:
            if not _has_tag(aln.peer, tag):
                continue
            field = aln.peer.opt(tag)
            if set(field) == set(['!']):
                e.append(BadQVsError.from_args(aln, tag, aln.qName))
        return e


class ValidateReadBaseInfo (ValidateReadBase):

    """
    Check that the values of specific tags are nucleotide sequences.
    """

    def validate(self, aln):
        # XXX should we do this without an index?
        if getattr(aln.bam, "pbi", None) is None:
            return None
        if not re.match(Constants.REGEX_BASECALLS, aln.DeletionTag()):
            return False
        elif not re.match(Constants.REGEX_BASECALLS, aln.SubstitutionTag()):
            return False
        return True

    def to_errors(self, aln):
        return [NonNucleotideError.from_args(aln, aln.qName)]


class ValidateReadSNR (ValidateReadBase):

    """
    Check for uninitialized values of the 'sn' field (signal-to-noise ratio
    over the high-quality region).
    """

    def validate(self, aln):
        if not _has_tag(aln.peer, 'sn'):
            return True
        if list(aln.peer.opt('sn')) == Constants.DEFAULT_SN:
            return False
        return True

    def to_errors(self, aln):
        return [UninitializedSNRError.from_args(aln, aln.qName)]


def _get_missing_pulse_tags(aln):
    from pbcore.io.align._BamSupport import PULSE_FEATURE_TAGS
    # XXX this may not be the behavior we want - it will only catch
    # pulse features that are present in ALL read groups
    pulseFeatures = aln.bam.pulseFeaturesAvailable()
    missing = []
    for featureKey in pulseFeatures:
        if not featureKey in PULSE_FEATURE_TAGS:
            continue
        featureTag = PULSE_FEATURE_TAGS[featureKey][0]
        if not _has_tag(aln.peer, featureTag):
            missing.append(featureTag)
    return missing


class ValidateReadPulseFeatures (ValidateReadBase):

    """
    Verify that any pulse feature declared in all @RG header manifests is
    present in each alignment.
    """

    def validate(self, aln):
        return len(_get_missing_pulse_tags(aln)) == 0

    def to_errors(self, aln):
        missing = _get_missing_pulse_tags(aln)
        return [PulseFeatureError.from_args(aln, aln.qName, ", ".join(sorted(missing)))]


class ValidateReadLocalContext (ValidateReadBase):

    """
    Check for mutually exclusive flags in the CX field - there is only
    one pair of these in the current spec.
    """

    def _get_errors(self, aln):
        if not aln.readType in Constants.REQUIRES_CONTEXT_TAGS:
            return []
        errors = []
        localContext = 0
        if aln.readType == Constants.READ_TYPE_SUBREAD:
            if not _has_tag(aln.peer, "cx"):
                return []  # XXX this is caught above
            localContext = aln.peer.opt("cx")
        else:
            # dummy value for CCS reads, which lack a CX tag
            localContext = 3
        if not (localContext & Constants.BARCODE_BEFORE or
                localContext & Constants.BARCODE_AFTER):
            for tag in Constants.REQUIRED_BARCODE_TAGS:
                if not _has_tag(aln.peer, tag):
                    errors.append(MissingAlignmentTagError.from_args(aln,
                                                                     aln.qName, tag, aln.readType))
        for x in [Constants.FORWARD_PASS]:
            if localContext & x and localContext & x * 2:
                return [LocalContextFlagsError.from_args(aln, aln.qName, localContext)]
        return errors

# XXX will this actually catch anything?


class ValidateReadCigar (ValidateReadBase):

    """Attempt to decode the CIGAR string"""

    def _get_errors(self, aln):
        if aln.isUnmapped:
            return []
        try:
            cigar = aln.unrolledCigar()
            for i, c in enumerate(cigar[:-1]):
                if ((c == 1 and cigar[i + 1] == 2) or
                        (c == 2 and cigar[i + 1] == 1)):
                    raise AlignmentCigarError("Adjacent insertion and " +
                                              "deletion in CIGAR string for QNAME %s" % aln.qName)
            refpos = aln.referencePositions()
        except IncompatibleFile as e:
            return []  # XXX already checked for 'M' earlier
        except ValidatorError as e:
            return [e]
        else:
            return []


class ValidateReadCigarMatches (ValidateReadBase):

    def _get_errors(self, aln):
        if aln.isUnmapped:
            return []
        cigar = aln.peer.cigarstring
        if re.search(Constants.REGEX_BAD_CIGAR, cigar):
            return [AlignmentCigarMatchError.from_args(aln, aln.qName)]
        return []


class ValidateReadTranscript (ValidateReadBase):

    """
    Attempt to extract the transcript (requires reference FASTA input).
    """

    def _get_errors(self, aln):
        if not aln.bam.isReferenceLoaded:
            return []
        try:
            t = aln.transcript()
        except Exception as e:
            return BadTranscriptError.from_args(aln, aln.qName, str(e))
        else:
            return []

#-----------------------------------------------------------------------
# runtime


def get_format_specific_args(parser):
    parser.add_argument("--unaligned", dest="aligned", action="store_false",
                        default=None,
                        help="Specify that the file should contain only " +
                             "unmapped alignments (DEFAULT: no requirement)")
    parser.add_argument("--unmapped", dest="aligned", action="store_false",
                        default=None, help="Alias for --unaligned")
    parser.add_argument("--aligned", dest="aligned", action="store_true",
                        default=None,
                        help="Specify that the file should contain only " +
                             "mapped alignments (DEFAULT: no requirement)")
    parser.add_argument("--mapped", dest="aligned", action="store_true",
                        default=None, help="Alias for --aligned")
    parser.add_argument("--contents", dest="contents", action="store",
                        choices=["SUBREAD", "CCS"], default=None,
                        help="Enforce read type")
    parser.add_argument("--reference", dest="reference", action="store",
                        default=None,
                        help="Path to optional reference FASTA file, used " +
                             "for additional validation of mapped BAM records")
    parser.add_argument("--permissive-headers",
                        dest="permissive_headers",
                        action="store_true",
                        help="Don't check chemistry/basecaller versions")
    return parser


def get_validators(aligned=None, contents=None,
                   include_file_validators=True,
                   validate_index=False,
                   permissive_headers=False):
    validators = []
    if include_file_validators:
        validators.extend([
            ValidateSorting(),
            ValidateContents(aligned, contents),
            # XXX not sure we want this...
            # ValidateFileName(),
            ValidateRandomAccess(),
        ])
        if validate_index:
            validators.append(ValidateIndex())
    validators.extend([
        ValidateReadGroupPlatform(),
        ValidateReadGroupId(),
    ])
    if not permissive_headers:
        validators.extend([
            ValidateReadGroupChemistry(),
            ValidateReadGroupBasecaller(),
        ])
    validators.extend([
        ValidateReadGroupPulseManifest(),
        ValidateReadUnique(),
        ValidateReadReadGroup(),
        ValidateReadQname(),
        ValidateReadLength(),
        # XXX disabling this, see bug 29224
        # ValidateReadQVs(),
        ValidateReadTags(),
        ValidateReadSNR(),
        ValidateReadCigar(),
        ValidateReadCigarMatches(),
        ValidateReadPulseFeatures(),
        # XXX is this actually done?
        # ValidateReadLocalContext(),
        ValidateReadMapped(aligned, contents),
        ValidateReadUnmapped(aligned, contents),
        ValidateReadTranscript(),
    ])
    # if contents is not None:
    #    validators.insert(2, ValidateReadType(contents))
    return validators


def validate_read_groups(ctx, validators, reader):
    """
    Separate loop over read groups (@RG) in BAM header.
    """
    for rg in reader.peer.header["RG"]:
        for v in validators:
            if isinstance(v, ValidateReadGroup):
                apply_validator_with_ctx(ctx, v, rg)


def _get_reader(file_name, reference=None):
    pbi_file = file_name + ".pbi"
    if os.path.isfile(pbi_file):
        reader_class = pbcore.io.IndexedBamReader
    else:
        reader_class = pbcore.io.BamReader

    def _reader(bam_file):
        return reader_class(bam_file, referenceFastaFname=reference)
    return _reader


def validate_bam(file_name,
                 reference=None,
                 aligned=None,
                 contents=None,
                 quick=False,
                 max_errors=None,
                 max_records=None,
                 validate_index=False,
                 permissive_headers=False):
    """
    Main API entry point for running BAM validation.

    Example:

    .. doctest::
        >>> from pbcoretools.pbvalidate.bam import validate_bam
        >>> from pbcore import data
        >>> bam_file = data.getBamAndCmpH5()[0]
        >>> errors, metrics = validate_bam(file_name=bam_file)
        >>> len(errors)
        231
        >>> print(errors[0])
        Mismatch between specified and expected read group ID: a9a22406c5 in file, but computed as b89a4406
        >>> unmapped_file = data.getUnalignedBam()
        >>> errors, metrics = validate_bam(file_name=unmapped_file)
        >>> len(errors)
        118
        >>> print(errors[0])
        This file has not been sorted by QNAME, or the header has not been updated.
        >>> errors, metrics = validate_bam(file_name=unmapped_file,
        ...     aligned=True, contents="CCS")
        >>> len(errors)
        120
    """
    validators = get_validators(aligned=aligned, contents=contents,
                                validate_index=validate_index,
                                permissive_headers=permissive_headers)
    e, m = run_validators(
        context_class=get_context_class(
            quick=quick,
            max_errors=max_errors,
            max_records=max_records),
        path=file_name,
        reader_class=_get_reader(file_name, reference),
        validators=validators,
        additional_validation_function=validate_read_groups)
    return e, m
