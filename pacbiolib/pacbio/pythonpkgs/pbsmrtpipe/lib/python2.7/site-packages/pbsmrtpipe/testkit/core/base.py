import logging
import unittest

log = logging.getLogger(__name__)

# Env var to define
# PYSIV_JOB_DIRNAME = 'PYSIV_JOB_DIRNAME'
# PYSIV_INPUT_XML = "PYSIV_INPUT_XML"
# PYSIV_WORKFLOW_XML = "PYSIV_WORKFLOW_XML"


class _TestBase(unittest.TestCase):

    """

    Abstract Base class for all unittest classes.

    The Job directory name *must* be defined in as env variable, otherwise an
    Exception will be raise in the class setup method.

    Call class setup to assign job_path
    1. Check for core directories to exist
    2. Check for core files to exist

    :raises ValueError, IOError:
    """

    # job_path = None
    # job_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'tests', 'data', 'job_01'

    # This is used a flag to show the decorator was applied and the appropriate
    # methods where added to the class.
    #_was_monkey_patch = False

    # ONLY Classes that are 'patched' are consider test classes.
    # this is the cls._is_test flag

    # Call class setup to assign job_path
    # 1. Check for core directories to exist
    # 2. Check for core files to exist
    # 3. Check for

#    @classmethod
#    def setUpClass(cls):
#        job_path = os.environ.get(PYSIV_JOB_DIRNAME, None)
#        # these are arbitrarily named
#        INPUT_XML = os.environ.get(PYSIV_INPUT_XML, None)
#        WORKFLOW_XML = os.environ.get(PYSIV_WORKFLOW_XML, None)
#
#        to_msg = lambda x: "Job is NOT defined. You need to define ENV Var '{x}'.".format(x=x)
#
#        if job_path is None:
#            msg = to_msg(PYSIV_JOB_DIRNAME)
#            log.error(msg)
#            raise ValueError(msg)
#        else:
#            cls.job_path = job_path
#            log.info("Using {x} as {f}.".format(x=PYSIV_JOB_DIRNAME, f=cls.job_path))
#
#        if os.path.exists(job_path):
#            log.info("{c} :: Using Job directory {x}".format(x=cls.job_path, c=cls.__name__))
#        else:
#            msg = "Job directory does NOT exist. {x}".format(x=job_path)
#            log.error(msg)
#            raise IOError(msg)
#
#        if INPUT_XML is None:
#            msg = to_msg(PYSIV_INPUT_XML)
#            log.error(msg)
#            raise ValueError(msg)
#        else:
#            cls.INPUT_XML = INPUT_XML
#            log.info("Using {x} as {f}".format(x=PYSIV_INPUT_XML, f=cls.INPUT_XML))
#
#        if WORKFLOW_XML is None:
#            msg = to_msg(PYSIV_WORKFLOW_XML)
#            log.error(msg)
#            raise ValueError(msg)
#        else:
#            cls.WORKFLOW_XML = WORKFLOW_XML
#            log.info("Using {x} as {f}".format(x=PYSIV_WORKFLOW_XML, f=cls.WORKFLOW_XML))
#
#        log.info("Completed setUpClass of {c}".format(c=cls.__name__))
#        return cls


# These are the only non-abstract Base classes that should be used.
# They should *always* be 'patched' via a class decorator and the
# _is_test cls variable should be set to True
class TestBase(_TestBase):
    _is_test = True
    job_dir = None
    service_access_layer = None
    job_id = None


class TestTaskBase(TestBase):
    pass
