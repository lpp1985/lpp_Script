import os
import logging

from .base import TestBase
from pbsmrtpipe.testkit.base import monkey_patch

log = logging.getLogger(__name__)


@monkey_patch
class TestCoreResources(TestBase):

    """
    Test to see of the core job directory structure and files exist
    """
    FILES = ('job.stdout', 'job.stderr')
    DIRS = ('html', 'workflow', 'tasks', 'logs')

    HTML_DIRS = ('css', 'images', 'js')
    HTML_FILES = ('index.html', 'settings.html', 'workflow.html', 'datastore.html')

    WORKFLOW_FILES = ('datastore.json', 'entry-points.json', 'report-tasks.json', 'options-task.json',
                      'options-workflow.json', 'workflow.dot', 'workflow.png', 'workflow.svg')

    def test_job_path_exists(self):
        self.assertTrue(os.path.exists(self.job_dir))
