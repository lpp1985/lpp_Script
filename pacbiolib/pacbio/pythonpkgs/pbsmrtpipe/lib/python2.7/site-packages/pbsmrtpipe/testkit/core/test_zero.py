import os
import logging

from .base import TestBase

log = logging.getLogger(__name__)


class TestSmrtPipeReturnCode(TestBase):
    # this is necessary for the test to be registered
    _is_test = True

    def setUp(self):
        self.pbsmrtpipe_log = os.path.join(self.job_dir, 'logs', 'pbsmrtpipe.log')

    def test_a_workflow_was_successful(self):
        """Test if SmrtPipe has a exited with returncode of 0."""
        # the test_a_X is done so the tests are run first

        # this should be fixed. The exit code should be clear.
        successful_str = 'Workflow was Successful'

        state = False
        error_msg = "Workflow did NOT successfully complete in {j}".format(j=self.job_dir)

        if os.path.exists(self.pbsmrtpipe_log):
            with open(self.pbsmrtpipe_log, 'r') as f:
                for line in f:
                    if successful_str in line:
                        state = True
                        break
        else:
            error_msg += " Unable to find log file '{f}'.".format(f=self.pbsmrtpipe_log)
            log.warn(error_msg)

        self.assertTrue(state, error_msg)

    def test_no_errors_in_log(self):
        """Check if ERROR or CRITICAL are in the smrtpipe.log"""

        # True means the job has no errors
        state = True
        error_msg = "ERROR or CRITICAL was found in pbsmrtpipe.log"

        if os.path.exists(self.pbsmrtpipe_log):
            with open(self.pbsmrtpipe_log, 'r') as f:
                for i, line in enumerate(f):
                    if '[ERROR]' in line or 'CRITICAL' in line:
                        state = False
                        error_msg = "(Line {i})".format(i=i) + line
                        break

        else:
            error_msg = "Unable to find {f}".format(f=self.pbsmrtpipe_log)

        self.assertTrue(state, error_msg)
