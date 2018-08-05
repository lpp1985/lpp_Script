
"""
Validation of PacBio dataset XML (and referenced files)
"""

import xml.etree.ElementTree as ET
from cStringIO import StringIO
from urlparse import urlparse
import xml.parsers.expat
import traceback
import itertools
import argparse
import logging
import os.path
import sys

try:
    from pyxb import exceptions_ as pyxbexceptions
except ImportError:
    class pyxbexceptions(object):

        class PyXBException(Exception):
            pass

        class ValidationError(Exception):
            pass

        class StructuralBadDocumentError(Exception):
            pass

from pbcore.io.dataset.DataSetReader import xmlRootType
from pbcore.io.dataset.DataSetIO import _dsIdToName
from pbcore.io.dataset import DataSet, DataSetValidator
from pbcore.io import IndexedBamReader, IndexedFastaReader
import pbcore.io

from pbcoretools.pbvalidate.core import (get_context_class, run_validators,
                                         ValidatorError, ValidateFile, ValidateFileObject)
from pbcoretools.pbvalidate import fasta
from pbcoretools.pbvalidate import bam

log = logging.getLogger(__name__)


class Constants(object):
    XML_NAMESPACE = "http://pacificbiosciences.com/PacBioBaseDataModel.xsd"


class DatasetTypes(object):
    BAM_DATASET = ["AlignmentSet", "ConsensusSet", "ConsensusAlignmentSet",
                   "SubreadSet"]
    FASTA_DATASET = ["BarcodeSet", "ContigSet", "ReferenceSet",
                     "GmapReferenceSet"]
    HDF5_DATASET = ["HdfSubreadSet"]
    ALL = BAM_DATASET + FASTA_DATASET + HDF5_DATASET


def _validate_read_groups(ctx, validators, reader):
    """
    Extra loop for validating just the read groups in .bam file headers.
    """
    if not DatasetReader.get_dataset_type(reader) in DatasetTypes.BAM_DATASET:
        return None
    try:
        bam_readers = reader.resourceReaders()
        for bam_reader in reader.resourceReaders():
            if bam_reader is None:
                log.warn("Skipping unopenable file")
                continue
            log.debug("Opened file: " + str(bam_reader))
            bam.validate_read_groups(ctx, validators, bam_reader)
    except IOError as e:
        # missing file, will be caught by ValidateResources
        return


class MissingIndexError(ValidatorError):
    MESSAGE_FORMAT = "Missing corresponding index files for the underlying " +\
        "raw data files."


class ReaderError(ValidatorError):
    MESSAGE_FORMAT = "Unexpected error reading dataset: %s.  This prevents " +\
        "any further validation functions from being run."


class MissingEncodingError(ValidatorError):
    MESSAGE_FORMAT = "This XML document is either missing the header " +\
        "or the encoding type is missing or wrong; all DataSet XMLs should " +\
        "explicitly specify UTF-8 encoding."


class XMLError(ValidatorError):
    MESSAGE_FORMAT = "XML schema error: %s"

# XXX currently redundant with underlying dataset API, untested


class MissingResourceIdError (ValidatorError):
    MESSAGE_FORMAT = "Found ExternalResource but no ResourceId is specified"


class MissingResourceError (ValidatorError):
    MESSAGE_FORMAT = "The external resource %s referenced by this dataset " +\
        "could not be located."


class ResourceOpenError (ValidatorError):
    MESSAGE_FORMAT = "The external resource %s referenced by this dataset " +\
        "is present, but could not be opened."


class DatasetTypeError(ValidatorError):
    MESSAGE_FORMAT = "The expected type was %s, but the loaded dataset has " +\
        "type %s.  (Note that this may render any additional validation " +\
        "errors irrelevant.)"


class FileNameError(ValidatorError):
    MESSAGE_FORMAT = "The dataset file %s is named incorrectly - datasets " +\
        "of type '%s' should have the extension '%s'."


class TimeStampedNameError(ValidatorError):
    MESSAGE_FORMAT = "This dataset does not contain the TimeStampedName " +\
        "attribute, which is a mandatory component of the current schema."


class NamespaceError(ValidatorError):
    MESSAGE_FORMAT = "The XML namespace '%s' for externalResources is " +\
        "different than the expected namespace '%s'; this may indicate " +\
        "that it is using obsolete schema."


class RootTagError(ValidatorError):
    MESSAGE_FORMAT = "The XML root tag '%s' does not match the declared " +\
        "dataset MetaType '%s'."


class NumRecordsError(ValidatorError):
    MESSAGE_FORMAT = "The number of records specified in the metadata (%s) "+\
        "is greater than the number of records in the data file(s) (%s)."


class ValidateXML(ValidateFile):

    def _get_errors(self, path):
        emsg = None
        try:
            DataSetValidator.validateFile(path, skipResources=True)
        except pyxbexceptions.StructuralBadDocumentError as e:
            emsg = "{t} ('<{n}>')".format(t=type(e).__name__,
                                          n=e.node.tagName)
        except pyxbexceptions.ValidationError as e:
            emsg = "{t}: {m}".format(t=type(e).__name__, m=e.details())
        except pyxbexceptions.PyXBException as e:
            emsg = "{t}: {m})".format(t=type(e).__name__, m=str(e.message))
        except Exception as e:
            emsg = str(e)
        if emsg is not None:
            return [XMLError.from_args(path, emsg)]
        return []

    def validate(self, path):
        return len(self._get_errors(path)) == 0

    def to_errors(self, path):
        return self._get_errors(path)


class ValidateRootTag(ValidateXML):

    def _get_errors(self, path):
        first = DataSet(path, strict=False)
        dsId = first.objMetadata.get('MetaType')
        xml_rt = xmlRootType(path)
        ds_name = _dsIdToName(dsId)
        if ds_name != xml_rt:
            if ds_name == "SubreadSet" and xml_rt == "HdfSubreadSet":
                return []
            return [RootTagError.from_args(path, xml_rt, _dsIdToName(dsId))]
        return []


class ValidateEncoding(ValidateXML):

    def __init__(self, *args, **kwds):
        self._has_xml_declaration = False
        super(ValidateEncoding, self).__init__(*args, **kwds)

    def _get_errors(self, path):
        self._has_xml_declaration = False
        e = []
        with open(path, 'r') as xmlfile:
            p = xml.parsers.expat.ParserCreate()

            def handle_xml_decl(version, encoding, standalone):
                if not self._has_xml_declaration:
                    if encoding is None or encoding.lower() != "utf-8":
                        e.append(MissingEncodingError.from_args(path))
                self._has_xml_declaration = True
            p.XmlDeclHandler = handle_xml_decl
            p.Parse(xmlfile.read())
        if not self._has_xml_declaration:
            e.append(MissingEncodingError.from_args(path))
        return e


class ValidateResources (ValidateFileObject):

    """
    Verify that the external resources specified in the XML file actually
    exist on the local filesystem
    """

    def _get_errors(self, file_obj):
        e = []
        for item in file_obj.externalResources:
            # XXX this is redundant
            if (not hasattr(item, "resourceId") or
                    item.resourceId is None):
                e.append(MissingResourceIdError.from_args(file_obj))
            else:
                continue
        return e

    def validate(self, file_obj):
        return len(self._get_errors(file_obj)) == 0

    def to_errors(self, file_obj):
        return self._get_errors(file_obj)


class ValidateResourcesOpen (ValidateResources):

    """
    Verify that the dataset object is capable of supplying open resource files.
    Note that since we assume ValidateResources is being run first, we can
    ignore any errors that result from the file(s) being absent entirely.
    """

    def _get_errors(self, file_obj):
        errors = []
        try:
            for r, f in itertools.izip(file_obj.externalResources,
                                       file_obj.resourceReaders()):
                if f is None:
                    errors.append(ResourceOpenError.from_args(file_obj,
                                                              urlparse(r.resourceId).path))
        except IOError as e:
            if e.filename is None or not os.path.exists(e.filename):
                log.info("File %s doesn't exist, skipping" % e.filename)
                return []
            log.warn("Encountered IOError opening %s" % e.filename)
            return [ResourceOpenError.from_args(file_obj, e.filename)]
        else:
            return errors


class ValidateIndex (ValidateResources):

    def _get_errors(self, file_obj):
        if not file_obj.isIndexed:
            return [MissingIndexError.from_args(file_obj)]
        return []


# TODO write test
class ValidateRandomAccess (ValidateResources):

    def _get_errors(self, file_obj):
        if len(file_obj.resourceReaders()) == 0 or not file_obj.isIndexed:
            return []
        errors = []
        for rr in file_obj.resourceReaders():
            if isinstance(rr, IndexedBamReader):
                errors.extend(bam.ValidateRandomAccess()._get_errors(rr))
            elif isinstance(rr, IndexedFastaReader):
                errors.extend(fasta.ValidateRandomAccess()._get_errors(rr))
            else:
                # logging.warn("Can't check indices for %s" % rr.filename)
                pass
        return errors


def _dataset_type(ds):
    return ds.objMetadata.get('MetaType').split(".")[-1]


class ValidateDatasetType (ValidateResources):

    """
    Verify that the opened dataset class name is the same as the user-supplied
    expected class name (if given).
    """

    def __init__(self, dataset_type=None):
        self.dataset_type = dataset_type

    def _get_errors(self, file_obj):
        dType = _dataset_type(file_obj)
        if self.dataset_type is None or self.dataset_type == "any":
            return []
        elif self.dataset_type != dType:
            # XXX see pbcore.io.dataset.DataSetIO:HdfSubreadSet - not sure
            # I understand what's going on here but I think it is a patch for
            # bug 27976
            if self.dataset_type == "HdfSubreadSet" and dType == "SubreadSet":
                return []
            return [DatasetTypeError.from_args(
                DatasetReader.get_dataset_object(file_obj),
                self.dataset_type,
                dType)]
        return []


class ValidateFileName (ValidateResources):

    """
    Check for consistency with file name conventions enforced when writing
    dataset files from the pbcore API in strict mode.
    """

    def __init__(self, file_name=None):
        self.file_name = file_name

    def _get_errors(self, file_obj):
        if self.file_name is not None:
            dataset_type = DatasetReader.get_dataset_type(file_obj)
            extension = ".%s.xml" % dataset_type.lower()
            if not self.file_name.endswith(extension):
                return [FileNameError.from_args(
                    DatasetReader.get_dataset_object(file_obj),
                    os.path.basename(self.file_name),
                    dataset_type,
                    extension)]
        return []


class ValidateMetadata(ValidateResources):

    """
    Check that the metadata in the XML file contains tags expected by the
    current schema.
    """

    def _get_errors(self, file_obj):
        ds = DatasetReader.get_dataset_object(file_obj)
        if not "TimeStampedName" in ds.objMetadata:
            return [TimeStampedNameError.from_args(ds)]
        return []


class ValidateNamespace(ValidateResources):

    def _get_errors(self, file_obj):
        ds = DatasetReader.get_dataset_object(file_obj)
        ns = ds.externalResources.namespace
        if ns != Constants.XML_NAMESPACE:
            return [NamespaceError.from_args(ds, ns, Constants.XML_NAMESPACE)]
        return []


class ValidateFileProxy (ValidateFileObject):

    """
    Wrapper for calling a file validator repeatedly on all files in the
    dataset.  Since it is assumed that ValidateResourcesOpen will be run first,
    failure to open a resource will be ignored with a log warning.
    """
    validator_class = None

    def __init__(self, **kwds):
        self._validator = self.validator_class(**kwds)
        self._errors = set([])

    def validate(self, file_obj):
        try:
            for _reader in file_obj.resourceReaders():
                if _reader is None:
                    # XXX if this happens, the file simply isn't present -
                    # which is handled separately by ValidateResourcesOpen
                    return True
                log.debug("Opened file: " + str(_reader))
                if not self._validator.validate(_reader):
                    errors_ = self._validator.to_errors(_reader)
                    self._errors.update(set(errors_))
        except IOError as e:
            #log.warn("Can't open file %s" % e.filename)
            return True
        else:
            return len(self._errors) == 0

    def to_errors(self, file_obj):
        return list(self._errors)


class ValidateContents (ValidateFileProxy):

    """Wrapper for pbvalidate.bam.ValidateContents"""
    validator_class = bam.ValidateContents


class ValidateSorting (ValidateFileProxy):

    """Wrapper for pbvalidate.bam.ValidateSorting"""
    validator_class = bam.ValidateSorting


class ValidateFastaRaw (ValidateFileObject):

    """Wrapper for pbvalidate.fasta.ValidateFastaRaw"""

    def __init__(self, **kwds):
        self._validator = fasta.ValidateFastaRaw(**kwds)
        self._errors = set([])

    def validate(self, file_obj):
        for extRes in file_obj.externalResources:
            path = urlparse(extRes.resourceId).path
            log.debug("Validating %s" % path)
            if not self._validator.validate(path):
                errors_ = self._validator.to_errors(path)
                self._errors.update(set(errors_))
        return len(self._errors) == 0

    def to_errors(self, file_obj):
        return list(self._errors)


class ValidateNumRecords(ValidateResources):

    def _get_errors(self, file_obj):
        ds = DatasetReader.get_dataset_object(file_obj)
        nr_metadata = ds.numRecords
        nr_actual = 0
        if ds.isIndexed:
            for rr in ds.resourceReaders():
                nr_actual += len(rr)
        else:
            for rr in ds.resourceReaders():
                nr_actual += len([rec for rec in rr])
        if nr_metadata > nr_actual:
            return [NumRecordsError.from_args(ds, nr_metadata, nr_actual)]
        return []



class DatasetReader (object):

    """
    Proxy for opening a dataset and iterating over records while avoiding an
    IOError if the resources can't be opened (since we already validate this
    separately).
    """

    def __init__(self, reader_class, file_name):
        self.reader_class = reader_class
        self.file_name = file_name
        self._reader = self.reader_class(self.file_name)

    @property
    def reader_name(self):
        return type(self._reader).__name__

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._reader.close()
        return False

    def __getattr__(self, name):
        return getattr(self._reader, name)

    def __iter__(self):
        try:
            for rec in self._reader:
                yield rec
        except IOError as e:
            pass

    def __str__(self):
        return str(self._reader)

    def __repr__(self):
        return repr(self._reader)

    @staticmethod
    def get_dataset_type(file_obj):
        """
        Return the class name of the actual DataSet-derived object, which
        could be either file_obj or its _reader attribute.  This allows the
        relevant validator classes to be called with either a DatasetReader or
        the underlying DataSet.
        """
        if isinstance(file_obj, DatasetReader):
            return file_obj.reader_name
        return type(file_obj).__name__

    @staticmethod
    def get_dataset_object(file_obj):
        if isinstance(file_obj, DatasetReader):
            return file_obj._reader
        return file_obj


def validate_dataset(
        file_name,
        dataset_type=None,
        reference=None,
        quick=False,
        max_errors=None,
        max_records=None,
        contents=None,
        aligned=None,
        validate_index=False,
        strict=False,
        permissive_headers=False):
    assert os.path.isfile(os.path.realpath(file_name))
    ds = None
    ReaderClass = getattr(pbcore.io, str(dataset_type), pbcore.io.openDataSet)
    log.debug("ReaderClass: %s" % ReaderClass.__name__)
    try:
        # XXX suppressing logging errors temporarily
        #logging.disable(logging.CRITICAL)
        try:
            ds = ReaderClass(file_name, strict=True)
        finally:
            pass #logging.disable(logging.NOTSET)
    except Exception as e:
        # XXX in strict mode the reader will cough up an IOError if the
        # requested dataset type does not agree with the XML.  if this happens
        # there's no point doing any additional validation.
        if False: #True:
            # XXX actually, it can cough up other errors too if there is
            # something wrong with the underlying files and it tries to read
            # them immediately.  Still treating this as a validation error, but
            # it may indicate bugs.
            _, _, ex_traceback = sys.exc_info()
            tb_lines = traceback.format_exception(e.__class__, e, ex_traceback)
            log.error("\n".join(tb_lines))
        errors = [ReaderError.from_args(file_name, str(e))]
        return errors, {}
    log.debug("Dataset type: %s" % ds.__class__.__name__)
    actual_dataset_type = _dataset_type(ds)
    log.debug("Actual type:  %s" % actual_dataset_type)
    if isinstance(ds, pbcore.io.SubreadSet) and contents is None:
        contents = "SUBREAD"
    elif isinstance(ds, pbcore.io.ConsensusReadSet) and contents is None:
        contents = "CCS"
    elif isinstance(ds, pbcore.io.AlignmentSet):
        pass
    validators = [
        ValidateEncoding(),
        ValidateRootTag(),
        ValidateResources(),
        ValidateDatasetType(dataset_type),
        ValidateMetadata(),
        ValidateNamespace(),
        ValidateRandomAccess(),
    ]
    if not actual_dataset_type in DatasetTypes.HDF5_DATASET:
        validators.extend([
            ValidateResourcesOpen(),
            ValidateNumRecords(),
        ])
        if validate_index:
            validators.append(ValidateIndex())
    if strict:
        validators.extend([
            ValidateXML(),
            ValidateFileName(file_name),
        ])
    additional_validation_function = None
    opened_class_name = ds.__class__.__name__
    # XXX not sure this is ideal - what if it opens as a ReferenceSet but we
    # asked for an AlignmentSet?  This is caught by ValidateDatasetType, but
    # we'd still check for Fasta file errors.
    if opened_class_name in DatasetTypes.FASTA_DATASET:
        validators_ = fasta.get_validators(validate_raw_format=False)
        validators_.insert(0, ValidateFastaRaw())
        validators.extend(validators_)
    elif opened_class_name in DatasetTypes.BAM_DATASET:
        validators_ = bam.get_validators(aligned=aligned,
                                         contents=contents,
                                         include_file_validators=False,
                                         permissive_headers=permissive_headers)
        validators_.insert(0, ValidateSorting())
        validators_.insert(0, ValidateContents(aligned=aligned,
                                               content_type=contents))
        validators.extend(validators_)
        additional_validation_function = _validate_read_groups

    def ReaderClass_wrapper(*args, **kwds):
        logging.disable(logging.CRITICAL)
        try:
            return DatasetReader(ReaderClass, *args, **kwds)
        finally:
            logging.disable(logging.NOTSET)
    context_class = get_context_class(
        quick=quick,
        max_errors=max_errors,
        max_records=max_records)
    errors, metrics = run_validators(
        context_class=context_class,
        path=file_name,
        reader_class=ReaderClass_wrapper,
        validators=validators,
        additional_validation_function=additional_validation_function)
    return errors, metrics


def get_parser():
    parser = argparse.ArgumentParser()
    return parser


def get_format_specific_args(parser):
    pass
