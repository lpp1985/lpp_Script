
"""
Validation of PacBio conventions for FASTA format.  Implicit in this is the
assumption that the input file will at least be processed by the parser in
pbcore.io without raising an exception.
"""

from __future__ import division, print_function
import argparse
import logging
import gzip
import re
import os.path as op
import sys

import pbcore.io

from pbcoretools.pbvalidate.core import *

log = logging.getLogger()


class Constants (object):
    OUTER_WHITESPACE = "^(\s)|(\s)$"
    # IUPAC nucleotide characters, minus '-' and '.'
    ILLEGAL_NUCLEOTIDES = "([^gatcuryswkmbdhvnGATCURYSWKMBDHVN]+)"
    ILLEGAL_NUC_STRICT = "([^gatcnGATCN]+)"
    ILLEGAL_IDENTIFIER = "([\,\:\"]+)"


class FastaError (ValidatorError):
    pass


class MissingIndexError (FastaError):
    MESSAGE_FORMAT = "Missing corresponding .fai index file"


class BadIndexError (FastaError):
    MESSAGE_FORMAT = "Couldn't access record %d directly (error: '%s'); " +\
        "this may indicate that the .fai file is corrupt."


class FastaRecordError (RecordValidatorError):
    pass


class EmptyLineError (FastaRecordError):
    MESSAGE_FORMAT = "Line %d is blank"


# XXX I would not have expected this to be a common problem, but multiple
# files in our reference library have it, and FastaReader crashes
class MissingSequenceError (FastaRecordError):
    MESSAGE_FORMAT = "No sequence found for record '%s'"


class WhitespaceError (FastaRecordError):
    MESSAGE_FORMAT = "Line %d ('%s') contains leading or trailing whitespace"


class NoWrappingError (FastaRecordError):
    MESSAGE_FORMAT = "Sequence '%s' does not have line wrapping"


class SeqWrappingError (FastaRecordError):
    MESSAGE_FORMAT = "Inconsistent line wrapping for sequence '%s'"


class GlobalWrappingError (FastaError):
    MESSAGE_FORMAT = "Inconsistent line wrapping conventions: please make " +\
        "sure all sequences are wrapped at the same line length (60-80 chars)"


class BadNucleotideError (FastaRecordError):
    MESSAGE_FORMAT = "Sequence '%s' contains illegal nucleotide " +\
        "character(s): %s"


class IdentifierAsteriskError (FastaRecordError):
    MESSAGE_FORMAT = "Identifier '%s' starts with an asterisk ('*')"


class BlankIdentifierError (FastaRecordError):
    MESSAGE_FORMAT = "Header '%s' has a blank identifier - you should not " +\
        "have a space after the '>' symbol"


# XXX see bug 24236
class ExtraGTError (FastaRecordError):
    MESSAGE_FORMAT = "Header '%s' has an additional '>' character - this " +\
        "may only appear at the start of the header line"


class BadIdentifierError (FastaRecordError):
    MESSAGE_FORMAT = "Identifier '%s' contains illegal character(s): %s"


class DuplicateIdError (FastaRecordError):
    MESSAGE_FORMAT = "The identifier '%s' appears more than once"


class ValidateFastaIndex(ValidateFileObject):

    """Check whether the input file was indexed by samtools"""

    def validate(self, file_obj):
        if isinstance(file_obj, pbcore.io.IndexedFastaReader):
            return True
        return False

    def to_error(self, file_obj):
        return MissingIndexError.from_args(file_obj)


# TODO write test
class ValidateRandomAccess (ValidateRandomAccessBase):
    ERROR_CLASS = BadIndexError
    INDEX_ATTR = "fai"


class ValidateFastaRecordHeader (ValidateRecord):

    """Check for non-blank FASTA identifier"""

    def validate(self, rec):
        if rec.id == "":
            return False
        return True

    def to_error(self, rec):
        return BlankIdentifierError.from_args(rec, rec.header)


class ValidateFastaIdentifier (ValidateRecord):

    """Check for illegal characters in the FASTA identifier"""

    def validate(self, rec):
        bad_header = re.findall(Constants.ILLEGAL_IDENTIFIER, rec.id)
        if bad_header:
            return False
        return True

    def to_error(self, rec):
        bad_header = re.findall(Constants.ILLEGAL_IDENTIFIER, rec.id)
        return BadIdentifierError.from_args(rec, rec.id, str(bad_header))


class ValidateFastaIdentifierStart (ValidateRecord):

    """Make sure the first character of the identifier isn't an asterisk"""

    def validate(self, rec):
        if len(rec.id) == 0:  # handled by ValidateFastaRecordHeader
            return True
        elif rec.id[0] == '*':
            return False
        return True

    def to_error(self, rec):
        return IdentifierAsteriskError.from_args(rec, rec.id)


class ValidateFastaIdentifierUnique (ValidateRecord):

    """Verify that a FASTA identifier is unique within the file"""

    def __init__(self):
        self._known_identifiers = set()

    def validate(self, rec):
        if rec.id in self._known_identifiers:
            return False
        self._known_identifiers.add(rec.id)
        return True

    def to_error(self, rec):
        return DuplicateIdError.from_args(rec, rec.id)


class ValidateFastaNucleotides (ValidateRecord):

    """
    Check for non-allowed nucleotide symbols.  By default all IUPAC codes
    will be accepted, but instantiating with strict=True will restrict the
    choices to acgtnACGTN.
    """

    def __init__(self, strict=False):
        self._use_regex = Constants.ILLEGAL_NUCLEOTIDES
        if strict:
            self._use_regex = Constants.ILLEGAL_NUC_STRICT

    def validate(self, rec):
        # already checked for trailing or leading whitespace
        bad_nuc = re.findall(self._use_regex, rec.sequence[:].strip())
        return len(bad_nuc) == 0

    def to_error(self, rec):
        bad_nuc = re.findall(self._use_regex, rec.sequence[:].strip())
        return BadNucleotideError.from_args(rec, rec.name, str(bad_nuc))


# FIXME this can only return a single error, but there may be multiple
# different problems with the raw format.
class ValidateFastaRaw (ValidateFile):

    """
    Since pbcore.io processes line wrapping automatically and doesn't care
    whether or not it's consistent, this validator will examine the raw file
    content line by line.
    """

    def __init__(self):
        self._errors = []

    def validate(self, path):
        all_seq_line_lengths = set([])
        single_line_lengths = set([])
        current_seq_line_lengths = []

        def check_current_sequence_lines(lines, label):
            if current_seq_line_lengths:
                if len(current_seq_line_lengths) == 1:  # all one line
                    pass #if current_seq_line_lengths[0] > 80:
                         #self._errors.append(
                         #   NoWrappingError.from_args(path, label))
                elif (current_seq_line_lengths[-1] >
                      current_seq_line_lengths[-2]):
                    # last line is longer than previous line
                    self._errors.append(
                        SeqWrappingError.from_args(path, label))
                elif len(set(current_seq_line_lengths[:-1])) > 1:
                    # multiple wrapping lengths in this sequence
                    self._errors.append(
                        SeqWrappingError.from_args(path, label))
                else:
                    # consistent wrapping for this sequence at least, save
                    # length for later
                    all_seq_line_lengths.add(current_seq_line_lengths[0])

        def _open(file_name):
            if file_name.endswith(".gz"):
                return gzip.open(file_name)
            else:
                # see https://www.python.org/dev/peps/pep-0278/
                return open(file_name, "rU")
        with _open(path) as f:
            prev_line_was_header = False
            prev_header = None
            for i, line in enumerate(f.readlines()):
                line = line.rstrip('\n')
                if line.strip() == "":
                    self._errors.append(EmptyLineError.from_args(path, i + 1))
                elif re.search(Constants.OUTER_WHITESPACE, line):
                    prev_line_was_header = False
                    self._errors.append(
                        WhitespaceError.from_args(path, i + 1, line))
                elif line.startswith(">"):
                    if prev_line_was_header:
                        self._errors.append(
                            MissingSequenceError.from_args(path, prev_header))
                    prev_line_was_header = True
                    check_current_sequence_lines(current_seq_line_lengths,
                                                 label=prev_header)
                    if ">" in line[1:]:
                        self._errors.append(ExtraGTError.from_args(path, line))
                    current_seq_line_lengths = []
                    prev_header = line[1:].strip()
                else:
                    prev_line_was_header = False
                    current_seq_line_lengths.append(len(line))
            if prev_line_was_header:
                self._errors.append(
                    MissingSequenceError.from_args(path, prev_header))
        check_current_sequence_lines(current_seq_line_lengths,
                                     label=prev_header)
        if len(all_seq_line_lengths) > 1:
            self._errors.append(GlobalWrappingError.from_args(path))
        return len(self._errors) == 0

    def to_error(self, path):
        return self._errors[0]


def get_format_specific_args(parser):
    pass
    # parser.add_argument("--strict", dest="strict", action="store_true",
    # help="Limit allowed nucleotide symbols to [actgnACTGN]")


def get_validators(strict=False, validate_raw_format=True,
                   validate_index=False):
    validators = [
        ValidateFastaRecordHeader(),
        ValidateFastaIdentifier(),
        ValidateFastaIdentifierStart(),
        ValidateFastaIdentifierUnique(),
        ValidateFastaNucleotides(strict=strict),
    ]
    if validate_raw_format:
        validators.insert(0, ValidateFastaRaw())
    if validate_index:
        validators.insert(0, ValidateFastaIndex())
    return validators


def _fasta_reader(file_name):
    idx_file = file_name + ".fai"
    if op.isfile(idx_file):
        return pbcore.io.IndexedFastaReader(file_name)
    else:
        return pbcore.io.FastaReader(file_name)


def validate_fasta(
        file_name,
        strict=False,
        quick=False,
        max_errors=None,
        validate_index=False,
        out=sys.stdout):
    """
    Main API entry point for validating Fasta files.

    Example:

    .. doctest::
        >>> from pbvalidate.Fasta import validate_fasta
        >>> from pbcore import data
        >>> filename = data.getTinyFasta()
        >>> errors, context = validate_fasta(file_name=filename)
        >>> print(len(errors))
        0
    """
    validators = get_validators(strict, validate_index=validate_index)
    errors, metrics = run_validators(
        context_class=get_context_class(quick, max_errors),
        path=file_name,
        reader_class=pbcore.io.FastaReader,
        validators=validators)
    return errors, metrics
