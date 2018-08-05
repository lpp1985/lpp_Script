"""Module providing support for sending progress events across HTTP"""
import os
import logging
import urllib
import urllib2
import json
import base64

log = logging.getLogger(__name__)


class DefaultProgressErrorHandler(urllib2.HTTPDefaultErrorHandler):

    def http_error_default(self, req, fp, code, msg, headers):
        result = urllib2.HTTPError(req.get_full_url(), code, msg, headers, fp)
        result.status = code
        log.error(result)
        return result


class ProgressPublisher(object):

    """
    Service implementation for relaying progress
    messages to a JSON-HTTP service
    """

    # Allowed States
    # FIXME Convert to Enum
    PROGRESS_STARTED = "Started"
    PROGRESS_UPDATE = "Updated"
    PROGRESS_FAILED = "Failed"
    PROGRESS_WAITING = "Waiting"
    PROGRESS_SUCCESS = "Completed"

    # Default constants/values
    DEFAULT_USER_AGENT = 'SMRTpipe/1.0 +http://www.pacificbiosciences.com'

    # this needs to be removed
    #DEFAULT_MARTIN_HOST = "smrtpipe:8081"

    #DEFAULT_PROGRESS_URL = "%(host)s/jobs/%%(job_id)s/status"

    def __init__(self, progress_url, user, password, debug=False, user_agent=DEFAULT_USER_AGENT):
        """
        host (str): Host to send updates to
        debug (bool): use debug mode
        user_agent (str): user Agent to add to HTTP headers
        """
        # this should have the form "%(host)s/jobs/%%(job_id)s/status"
        self.progress_url = progress_url

        # Authentication
        self.user = user
        self.password = password

        self.user_agent = user_agent

        self._debug = debug

    def started(self, message='Started', value=0, jobStage=None, moduleName=None):
        self.post_progress(code=ProgressPublisher.PROGRESS_STARTED, message=message,
                           value=value, jobStage=jobStage, moduleName=moduleName)

    def update(self, message='Update', value=0, jobStage=None, moduleName=None):
        self.post_progress(code=ProgressPublisher.PROGRESS_UPDATE, message=message,
                           value=value, jobStage=jobStage, moduleName=moduleName)

    def waiting(self, message='Waiting', value=100, jobStage=None, moduleName=None):
        self.post_progress(code=ProgressPublisher.PROGRESS_WAITING, message=message,
                           value=value, jobStage=jobStage, moduleName=moduleName)

    def failed(self, message='Failed', value=100, jobStage=None, moduleName=None):
        self.post_progress(code=ProgressPublisher.PROGRESS_FAILED, message=message,
                           value=value, jobStage=jobStage, moduleName=moduleName)

    def success(self, message='Success', value=100, jobStage=None, moduleName=None):
        self.post_progress(code=ProgressPublisher.PROGRESS_SUCCESS, message=message,
                           value=value, jobStage=jobStage, moduleName=moduleName)

    def post_progress_event(self, event):
        self._post_progress(event)

    def post_progress(self, code=None, message=None, value=0,
                      jobStage=None, moduleName=None):
        """
        Kwargs:
            event (ProgressEvent): event to send to host
            code ():
            message (str, None):
            value (int):
            jobStage (str, None):
            moduleName (str, None):
        """
        # don't try to post progress if it's obvious we're running outside
        # of an automation system

        _event = ProgressEvent(code=code, message=message, value=value,
                               jobStage=jobStage, moduleName=moduleName)
        self._post_progress(_event)

    def _post_progress(self, event):
        #url = self._url % {'job_id': self._jobid}
        url = self.progress_url
        data = urllib.urlencode({'progress': event.to_json()})

        request = urllib2.Request(url, data=data)
        request.add_header('User-Agent', self.user_agent)
        key = 'Basic %s' % (base64.b64encode("{u}:{p}".format(u=self.user, p=self.password)))
        request.add_header('Authorization', key)

        opener = urllib2.build_opener(DefaultProgressErrorHandler())

        try:
            response = opener.open(request)
            response.close()
            log.info("Job Progress '%s' event POSTED to %s" % (event._code, str(url)))
            log.debug("The data string is %s" % (data))
        except Exception as e:
            log.info("Job Progress '%s' event POSTED to %s fail: %s" % (event._code, str(url), str(e)))
            log.debug("Request is \n%s" % request_to_string(request))
            pass


class ProgressEvent(object):

    def __init__(self, code=None, message=None, value=0, jobStage=None, moduleName=None):
        """
        Kwargs:
            code (int):
            message (str):
            value (int):
            jobStage (str):
            moduleName (str):
        """
        self._code = code
        self._message = message
        self._value = value
        self._stage = jobStage
        self._module = moduleName

    def to_dict(self):
        return {'code': self._code, 'message': self._message,
                'value': self._value, 'jobStage': self._stage,
                'moduleName': self._module}

    def to_json(self):
        return json.dumps(self.to_dict())


def request_to_string(request):
    """for debugging"""
    buffer = []
    buffer.append('Method: %s' % request.get_method())
    buffer.append('Host: %s' % request.get_host())
    buffer.append('Selector: %s' % request.get_selector())
    buffer.append('Data: %s' % request.get_data())
    return os.linesep.join(buffer)
