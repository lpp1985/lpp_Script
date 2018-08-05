import ConfigParser
import logging
import functools
import json
import abc
import os

log = logging.getLogger(__name__)

EXE = "pbsmrtpipe"

__all__ = ['ButlerTask',
           'ButlerWorkflow',
           'config_parser_to_butler']

__author__ = "Michael Kocher"


class TestkitCfgParserError(ValueError):
    pass


class Constants(object):

    """Allowed values in cfg file."""
    CFG_TASK = 'pbsmrtpipe:task'
    CFG_WORKFLOW = 'pbsmrtpipe:pipeline'

    CFG_JOB_ID = "id"

    CFG_ENTRY_POINTS = 'entry_points'

    CFG_PRESET_XML = 'preset_xml'
    CFG_PRESET_JSON = 'preset_json'
    CFG_WORKFLOW_XML = 'pipeline_xml'
    CFG_TASK_ID = 'task_id'
    CFG_BASE_EXE = 'base_exe'

    CFG_OUTPUT_DIR = 'output_dir'

    CFG_DEBUG = 'debug'
    CFG_MOCK = 'mock'


class Butler(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, job_id, output_dir, entry_points, preset_json,
                 preset_xml, debug,
                 force_distribute=None, force_chunk=None, base_exe=EXE):
        self.output_dir = output_dir
        self.entry_points = entry_points
        self.preset_json = preset_json
        self.preset_xml = preset_xml
        self.debug_mode = debug

        # None means no override, True|False means override
        self.force_distribute = force_distribute
        self.force_chunk = force_chunk

        # this needs to be set in the Butler.cfg file.
        self.job_id = job_id

        self.base_exe = base_exe

    def __repr__(self):
        _d = dict(k=self.__class__.__name__, p=self.prefix)
        return "<{k} {p} >".format(**_d)

    @abc.abstractproperty
    def prefix(self):
        # Used in the repr
        return ""

    def to_cmd(self):
        return _to_pbsmrtpipe_cmd(self.prefix, self.output_dir,
                                  self.entry_points, self.preset_json,
                                  self.preset_xml,
                                  self.debug_mode, self.force_distribute,
                                  self.force_chunk, self.base_exe)



class ButlerWorkflow(Butler):

    def __init__(self, job_id, output_dir, pipeline_id, workflow_xml, entry_points, preset_json_path, preset_xml_path, debug, force_distribute=None, force_chunk=None, base_exe=EXE):
        super(ButlerWorkflow, self).__init__(job_id, output_dir, entry_points, preset_json_path, preset_xml_path, debug, force_distribute=force_distribute, force_chunk=force_chunk, base_exe=base_exe)
        assert [workflow_xml, pipeline_id].count(None) == 1
        self.workflow_xml = workflow_xml
        self.pipeline_id = pipeline_id

    @property
    def prefix(self):
        if self.pipeline_id is not None:
            return "pipeline-id {i}".format(i=self.pipeline_id)
        else:
            return "pipeline {i}".format(i=self.workflow_xml)

    @staticmethod
    def from_json(file_name, force_distribute=None, force_chunk=None):
        with open(file_name) as json_f:
            d = json.load(json_f)
            assert d.get('jobType', "pbsmrtpipe") == "pbsmrtpipe"
            return ButlerWorkflow(
                job_id=d['testId'],
                output_dir=d.get('outputDir', "job_output"),
                pipeline_id=d.get("pipelineId", None),
                workflow_xml=d.get("workflowXml", None),
                entry_points={e['entryId']:e['path'] for e in d['entryPoints']},
                preset_xml_path=d.get('presetXml', None),
                preset_json_path=d.get("presetJson", None),
                debug=d.get("debug", False),
                force_distribute=force_distribute,
                force_chunk=force_chunk)


class ButlerTask(Butler):

    def __init__(self, job_id, output_dir, task_id, entry_points, preset_json, preset_xml, debug, force_distribute=None, force_chunk=None):
        super(ButlerTask, self).__init__(job_id, output_dir, entry_points, preset_json, preset_xml, debug, force_distribute=force_distribute, force_chunk=force_chunk)
        self.task_id = task_id

    @property
    def prefix(self):
        return "task {i}".format(i=self.task_id)


def _to_pbsmrtpipe_cmd(prefix_mode, output_dir, entry_points_d, preset_json, preset_xml, debug, force_distribute, force_chunk, base_exe=EXE):
    ep_str = " ".join([" -e '" + ":".join([k, v]) + "'" for k, v in entry_points_d.iteritems()])
    d_str = '--debug' if debug else " "
    p_str = " " if preset_xml is None else "--preset-xml={p}".format(p=preset_xml)
    j_str = " " if preset_json is None else "--preset-json={j}".format(j=preset_json)
    m_str = ' '

    force_distribute_str = ''
    if isinstance(force_distribute, bool):
        m = {True: '--force-distributed', False: '--local-only'}
        force_distribute_str = m[force_distribute]

    force_chunk_str = ''
    if isinstance(force_chunk, bool):
        m = {True: '--force-chunk-mode', False: '--disable-chunk-mode'}
        force_chunk_str = m[force_chunk]

    _d = dict(x=base_exe, e=ep_str, d=d_str, p=p_str, j=j_str, m=prefix_mode, o=output_dir, k=m_str,
              f=force_distribute_str, c=force_chunk_str)
    cmd = "{x} {m} {c} {d} {e} {p} {j} {k} {f} --output-dir={o}"
    return cmd.format(**_d)


to_task_cmd = functools.partial(_to_pbsmrtpipe_cmd, 'task')
to_workflow_cmd = functools.partial(_to_pbsmrtpipe_cmd, 'pipeline')


def _parse_or_default(section, key, p, default):
    if p.has_option(section, key):
        return p.get(section, key)
    return default


def _parse_preset_xml(section_name, p, base_dir):
    v = _parse_or_default(section_name, Constants.CFG_PRESET_XML, p, None)
    if v is None:
        return None
    else:
        p = v if os.path.isabs(v) else os.path.join(base_dir, v)
        if os.path.exists(p):
            return p
        else:
            raise IOError("Unable to find preset XML '{p}'".format(p=p))


def _parse_preset_json(section_name, p, base_dir):
    v = _parse_or_default(section_name, Constants.CFG_PRESET_JSON, p, None)
    if v is None:
        return None
    else:
        p = v if os.path.isabs(v) else os.path.join(base_dir, v)
        if os.path.exists(p):
            return p
        else:
            raise IOError("Unable to find preset XML '{p}'".format(p=p))


def _parse_debug_mode(section_name, p):
    return bool(_parse_or_default(section_name, Constants.CFG_DEBUG, p, True))


def _parse_entry_points(p, root_dir_name):
    """

    Files may be defined relative to the butler.cfg file or absolute paths

    """
    ep_d = {}
    ep_keys = p.options(Constants.CFG_ENTRY_POINTS)

    for ep_key in ep_keys:
        v = p.get(Constants.CFG_ENTRY_POINTS, ep_key)
        if not os.path.isabs(v):
            v = os.path.join(root_dir_name, v)

        ep_d[ep_key] = os.path.abspath(v)

    return ep_d


def _parse_entry_points_and_presets(section_name, p, root_dir):
    return _parse_entry_points(p, root_dir), _parse_preset_xml(section_name, p, root_dir), _parse_preset_json(section_name, p, root_dir)


def _to_parse_workflow_config(job_output_dir, base_dir):
    """

    :param job_output_dir: Job output directory
    :param base_dir:  base directory of the butler.cfg file
    :return:
    """
    def _parse_workflow_config(p):
        ep_d, preset_xml, preset_json = _parse_entry_points_and_presets(Constants.CFG_WORKFLOW, p, base_dir)
        x = p.get(Constants.CFG_WORKFLOW, Constants.CFG_WORKFLOW_XML)

        if not os.path.isabs(x):
            x = os.path.join(base_dir, x)

        if not os.path.exists(x):
            raise IOError("Unable to find pipeline XML '{x}'".format(x=x))

        d = _parse_debug_mode(Constants.CFG_WORKFLOW, p)
        workflow_xml = os.path.abspath(x)

        # FIXME. This should be defined in cfg file.
        default_job_id = os.path.basename(base_dir)
        job_id = _parse_or_default(Constants.CFG_WORKFLOW, Constants.CFG_JOB_ID, p, default_job_id)
        base_exe = _parse_or_default(Constants.CFG_WORKFLOW, Constants.CFG_BASE_EXE, p, EXE)

        return ButlerWorkflow(job_id, job_output_dir, None, workflow_xml, ep_d, preset_json, preset_xml, d, base_exe=base_exe)

    return _parse_workflow_config


def _to_parse_task_config(output_dir, base_dir):
    def _parse_task_config(p):
        # FIXME. This should be defined in cfg file.
        default_job_id = os.path.basename(base_dir)
        ep_d, preset_xml, preset_json = _parse_entry_points_and_presets(Constants.CFG_TASK, p, base_dir)
        job_id = _parse_or_default(Constants.CFG_TASK, Constants.CFG_JOB_ID, p, default_job_id)
        task_id = p.get(Constants.CFG_TASK, Constants.CFG_TASK_ID)
        d = _parse_debug_mode(Constants.CFG_TASK, p)
        b = ButlerTask(job_id, output_dir, task_id, ep_d, preset_json, preset_xml, d, force_distribute=False)

        return b

    return _parse_task_config


def config_parser_to_butler(file_path):
    """

    :param file_path: path to butler config file
    :return: Butler instance

    :rtype: Butler
    """
    if file_path.endswith(".json"):
        return ButlerWorkflow.from_json(file_path)
    # this is weak. Needs error handling.
    p = ConfigParser.ConfigParser()
    _ = p.read(file_path)

    # paths within the config file can be relative the butler.cfg file

    base_dir = os.path.dirname(file_path)

    # pbsmrtpipe will make the directory if it doesn't exist
    default_output_dir = os.path.join(base_dir, 'job_output')
    output_dir = _parse_or_default(Constants.CFG_WORKFLOW, Constants.CFG_OUTPUT_DIR, p, default_output_dir)
    # output_dir must be defined relative to testkit.cfg or an absolute path
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(base_dir, output_dir)

    if p.has_section(Constants.CFG_WORKFLOW):
        func = _to_parse_workflow_config(output_dir, base_dir)
    elif p.has_section(Constants.CFG_TASK):
        func = _to_parse_task_config(output_dir, base_dir)
    else:
        _d = dict(x=Constants.CFG_WORKFLOW, y=Constants.CFG_TASK, f=file_path)
        raise TestkitCfgParserError("Expected section {x} or {y} in {f}".format(**_d))

    butler = func(p)

    return butler
