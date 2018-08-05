
from __future__ import division, print_function
import sys
import traceback
from collections import namedtuple
import logging

_log = logging.getLogger(__name__)


def _log_traceback(ex, ex_traceback):
    tb_lines = traceback.format_exception(ex.__class__, ex, ex_traceback)
    tb_text = ''.join(tb_lines)
    _log.error(tb_text)


class ValidatorError(AssertionError):

    """
    Base class for violations of a PacBio file format spec.  The errors are
    hashed by the class name, which facilitates hiding many similar errors
    in the output.  The 'message' class attribute should be a format string;
    each sublcass should then be instantiated with appropriate arguments
    to be interpolated (if applicable).
    """
    MESSAGE_FORMAT = None

    def __init__(self, message, object_ref=None):
        super(ValidatorError, self).__init__(message)
        self.object_ref = object_ref

    def __repr__(self):
        _d = dict(k=self.__class__.__name__, n=self.message)
        return "<{k} {n} >".format(**_d)

    @classmethod
    def from_args(cls, object_ref, *args):
        return cls(cls.MESSAGE_FORMAT % args, object_ref)

    def test_name(self):
        """Generate a name for pseudo-unit test, which I am abusing to
        produce JUnit XML output"""
        return "test_" + self.__class__.__name__[:-5]

    # the next two methods are trickery for saving in a dict or set
    def __hash__(self):
        return hash(self.__class__.__name__)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


class RecordValidatorError (ValidatorError):

    """
    Sub-class used for validation of individual records in a file - no extra
    functionality, just easily identifiable by the context manager.
    """
    pass


class ValidateBase(object):

    """Base class for all validation classes"""

    def validate(self, item):
        raise NotImplementedError

    def __repr__(self):
        return "<{k} >".format(k=self.__class__.__name__)

    def to_error(self, item):
        # override if you want a custom error
        m = "{r} failed validation".format(f=item)
        return ValidatorError(m)

    def to_errors(self, item):
        """
        Alternate API: return a list instead (for validators that check
        related properties with multiple possible errors)
        """
        return [self.to_error(item)]


class ValidateRecord(ValidateBase):

    """Base class for validating an individual record in a file"""

    def validate(self, record):
        return True

    def to_error(self, record):
        # override if you want a custom error
        m = "Record {r} failed".format(r=record)
        return ValidatorError(m)


class ValidateFile(ValidateBase):

    """Base class for validating external file properties (name, etc.)"""

    def validate(self, path):
        raise NotImplementedError

    def to_error(self, path):
        m = "File validator error {p}".format(p=path)
        return ValidatorError(m)


class ValidateFileObject (ValidateBase):

    """Base class for validating global internal file properties (such as
    header information"""

    def validate(self, file_obj):
        raise NotImplementedError

    def to_error(self, file_obj):
        m = "File properties validator error {p}".format(p=file_obj)
        return ValidatorError(m)


# TODO write test
class ValidateRandomAccessBase (ValidateFileObject):
    """
    Base class for validating indexed file types (FASTA, BAM, and equivalent
    datasets) where we expect a __getitem__ method.  This just checking that
    the last expected item is addressible.
    """
    ERROR_CLASS = None
    INDEX_ATTR = None

    def _get_errors(self, file_obj):
        if hasattr(file_obj, self.INDEX_ATTR):
            n_records = len(file_obj)
            if n_records > 0:
                try:
                    last_record = file_obj[n_records - 1]
                except Exception as e:
                    emsg = "{t}: {m}".format(t=type(e).__name__, m=str(e))
                    return [self.ERROR_CLASS.from_args(file_obj, n_records - 1,
                                                       emsg)]
        return []

    def validate(self, file_obj):
        return len(self._get_errors(file_obj)) == 0

    def to_errors(self, file_obj):
        return self._get_errors(file_obj)


class ValidationStopException(BaseException):
    # Use this to exit
    pass

# Store validate
ValidMetric = namedtuple("ValidMetric", ["metric_type", "object_ref", ])


class ValidatorErrorContext(object):

    def __init__(self, errors, metrics):
        self._errors = errors
        self._other_errors = []
        self.metrics = metrics

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self._errors  # and not self._other_errors

    def __repr__(self):
        _d = dict(k=self.__class__.__name__, n=len(self._errors))
        return "<{k} {n} errors >".format(**_d)

    def add_validation_error(self, error):
        # This is for subclasses to override and raise a
        # ValidationStopException if
        self._errors.append(error)

    def add_other_error(self, error):
        self._other_errors.append(error)

    def add_validate_metric(self, metric):
        key = metric.metric_type.__name__
        if key in self.metrics:
            self.metrics[key] += 1
        else:
            self.metrics[key] = 1


class ValidatorContextFirstError(ValidatorErrorContext):

    def add_validation_error(self, error):
        self._errors.append(error)
        raise ValidationStopException(
            "Found first error. Failed to validate {e}".format(e=self._errors))


class ValidatorContextFirstBadRecord(ValidatorErrorContext):

    def __init__(self, *args, **kwds):
        super(ValidatorContextFirstBadRecord, self).__init__(*args, **kwds)
        self._records = set([])

    def add_validation_error(self, error):
        if isinstance(error, RecordValidatorError):
            self._records.add(error.object_ref)  # assume it's hashable...
            if len(self._records) > 1:
                raise ValidationStopException("At least two bad records, " +
                                              "stopping now.")
        self._errors.append(error)


class ValidatorContextMaxErrors(ValidatorErrorContext):

    """
    Context manager that will raise ValidationStopException after the specified
    number of errors has been recorded.
    """

    def __init__(self, errors, metrics, max_errors):
        super(ValidatorContextMaxErrors, self).__init__(errors, metrics)
        self.max_errors = max_errors

    def add_validation_error(self, error):
        self._errors.append(error)
        if len(self._errors) >= self.max_errors:
            _d = dict(m=self.max_errors, e=error)
            raise ValidationStopException(
                "Exceeded max errors threshold {m}. Failed to validate. Last Validation {e}".format(**_d))


class ValidatorContextMaxRecords(ValidatorErrorContext):

    """
    Context manager that will raise ValidationStopException after the specified
    number of records have been examined (regardless of whether they fail
    validation or not).
    """

    def __init__(self, errors, metrics, max_records):
        super(ValidatorContextMaxRecords, self).__init__(errors, metrics)
        self.max_records = max_records
        self._records = set([])

    def _save_record(self, object_ref):
        self._records.add(object_ref)
        if len(self._records) >= self.max_records:
            raise ValidationStopException("Reached maximum record count, " +
                                          "stopping now.")

    def add_validation_error(self, error):
        if isinstance(error, RecordValidatorError):
            self._save_record(error.object_ref)
        self._errors.append(error)

    def add_validate_metric(self, metric):
        if issubclass(metric.metric_type, ValidateRecord):  # yuck!
            self._save_record(metric.object_ref)
        super(ValidatorContextMaxRecords, self).add_validate_metric(metric)


class ValidatorContextUniqueErrors (ValidatorErrorContext):

    """Context for collecting only unique errors and ignoring repeated
    occurrences."""

    def __init__(self, errors, metrics):
        ValidatorErrorContext.__init__(self, errors, metrics)
        self._known_errors = set()
        self.n_skipped = 0

    def add_validation_error(self, error):
        if not error in self._known_errors:
            self._errors.append(error)
            self._known_errors.add(error)
        else:
            self.n_skipped += 1


class ValidatorContextFailFirst(ValidatorContextMaxErrors):

    """Immediately exit after the first failure"""

    def __init__(self, errors, metrics):
        super(ValidatorContextFailFirst, self).__init__(
            errors, metrics, max_errors=1)


def apply_validator(v, item):
    """
    Applies validation to an item and logs the error

    :param v: validator instance
    :param item: record to validate
    :return:
    """
    try:
        #_log.debug("Validating {v} {i}".format(v=v, i=item))
        # XXX using str(item) can be prohibitively slow for .bam files when
        # the reference is loaded!
        #_log.debug("Validating {v} {i}".format(v=v, i=repr(item)))
        if v.validate(item):
            return ValidMetric(v.__class__, object_ref=item)
        else:
            return v.to_errors(item)
    except Exception as e:
        msg = "Unexpected {t} error validating {i}: {e}".format(i=item, e=e,
                                                                t=type(e).__name__)
        _log.error(msg)
        _, _, ex_traceback = sys.exc_info()
        _log_traceback(e, ex_traceback)
        return ValidatorError(msg)


def apply_validator_with_ctx(ctx, v, item):
    """
    Applies the validator to the item and pushes the results to the context.

    :param ctx: Validation Context
    :param v: Validator
    :param item: Record to validate
    :return:
    """

    result = apply_validator(v, item)
    if isinstance(result, ValidMetric):
        ctx.add_validate_metric(result)
    else:
        assert result is not None
        if isinstance(result, list):
            for result_ in result:
                ctx.add_validation_error(result_)
        else:
            ctx.add_validation_error(result)


def get_context_class(quick=False, max_errors=None, max_records=None):
    if quick:
        # return ValidatorContextFirstBadRecord
        def _context_class(*args, **kwds):
            kwds['max_records'] = 100
            return ValidatorContextMaxRecords(*args, **kwds)
        return _context_class
    elif max_errors is not None:
        def _context_class(*args, **kwds):
            kwds['max_errors'] = max_errors
            return ValidatorContextMaxErrors(*args, **kwds)
        return _context_class
    elif max_records is not None:
        def _context_class(*args, **kwds):
            kwds['max_records'] = max_records
            return ValidatorContextMaxRecords(*args, **kwds)
        return _context_class
    return ValidatorErrorContext


def run_validators(context_class, path, reader_class, validators,
                   additional_validation_function=None):
    """Runs a validations in a context and returns a tuple of errors, metrics

    Errors are a list of Error instances
    Metrics are a dict of {ValidMetric: Int} number of a specific validation metric

    :param context_class: Context to run in
    :type context_class: ValidatorErrorContext
    :param path: path to file
    :type path: str
    :param reader_class: ReaderBase
    :param validators: list of validators
    :param additional_validation_function:
    :return:
    """

    errors = []
    # Metric -> Int
    metrics = {}
    try:
        with context_class(errors, metrics) as ctx:
            for v in validators:
                if isinstance(v, ValidateFile):
                    apply_validator_with_ctx(ctx, v, path)

            # FIXME this should probably be refactored to catch reader errors
            # separately from logic errors elsewhere
            with reader_class(path) as reader:
                for v in validators:
                    if isinstance(v, ValidateFileObject):
                        apply_validator_with_ctx(ctx, v, reader)
                if additional_validation_function is not None:
                    additional_validation_function(ctx, validators, reader)
                for record in reader:
                    for v in validators:
                        if isinstance(v, ValidateRecord):
                            apply_validator_with_ctx(ctx, v, record)
    except ValidationStopException as e:
        # This doesn't necessarily mean the validation failed. This determined
        # from errors list
        _log.debug("Stopping validation due to {e}".format(e=e))
        return errors, metrics
    except Exception as e:
        # this means it failed
        msg = "Unexpected context validator error for file '{p}' Error {e}".format(
            p=path, e=e)
        _log.error(msg)
        _, _, ex_traceback = sys.exc_info()
        _log_traceback(e, ex_traceback)
        errors.append(ValidatorError(msg))

    return errors, metrics


def run_validators_fail_quick(path, reader_class, validators,
                              additional_validation_function=None):
    """Run validators and fail after the first error

    :param path: Path to string
    :param validators:
    """

    # this could also be done using functools. Doing this explicitly for the docstring
    # func = functools.partial(run_validators, ValidatorContextFailFirst)
    return run_validators(ValidatorContextFailFirst, path, reader_class, validators,
                          additional_validation_function=additional_validation_function)


def run_validators_expect_errors(path, reader_class, validators,
                                 expected_error_types,
                                 additional_validation_function=None,
                                 quick=False):

    errors, metrics = run_validators(ValidatorErrorContext, path, reader_class, validators,
                                     additional_validation_function=additional_validation_function, quick=quick)

    # check that expected errors are in errors returned
    # if the errors are found, then return the proper metrics and errors
    # FIXME Nat to implement
    return errors, metrics
