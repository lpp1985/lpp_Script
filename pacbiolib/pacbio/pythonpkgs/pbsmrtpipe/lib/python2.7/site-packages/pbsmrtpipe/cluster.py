import importlib
import os
import logging
from string import Template
import types
import warnings

import pbsmrtpipe

# not sure this is the greatest idea. This is public API entry point
# however this will require package data modification to change templates
# For example
# .virtualenvs/dev_pbsmrtpipe/lib/python2.7/site-packages/pbsmrtpipe/cluster_templates/my_sge
from pbsmrtpipe.cluster_templates import CLUSTER_TEMPLATE_DIR

log = logging.getLogger(__name__)

# Built-in Cluster Templates
_PACKAGE_DIR = os.path.dirname(os.path.dirname(pbsmrtpipe.__file__))

__all__ = ['load_installed_cluster_templates',
           'load_installed_cluster_templates_by_name',
           'ClusterTemplate',
           'ClusterTemplateRender',
           'CLUSTER_TEMPLATE_DIR']


class Constants(object):
    # Each template must be implemented.
    START = "start"
    STOP = "stop"

    @classmethod
    def all(cls):
        return cls.START, cls.STOP

    @classmethod
    def is_valid(cls, name):
        if name in cls.all():
            return True

        return False


def _get_cluster_files(template_dir, suffix='.tmpl'):
    """Returns a dict of template name type to file path

    Example {start:/path/to/start.tmpl}

    TODO: add more error checking and informative error messages.
    """
    if not os.path.exists(template_dir):
        raise IOError("Unable to find cluster templates root directory '{t}'.".format(t=template_dir))

    paths = [os.path.join(template_dir, x + suffix) for x in Constants.all()]

    for p in paths:
        if not os.path.isfile(p):
            msg = "cluster template loading error. Unable to find template file '{f}'".format(f=p)
            log.exception(msg)
            raise IOError(msg)

    return dict(zip(Constants.all(), paths))


def _template_file_to_str(template_file):
    with open(template_file) as f:
        s = f.read()
    return s


def load_cluster_templates_from_dir(cluster_model_dir):
    """
    Load tmpl files from dir and return list of ClusterTemplate instances.

    The directory should contain 'start.tmpl' and 'stop.tmpl' Cluster templates

    For example, /path/to/cluster_templates/my_sge

    :returns: A list of ClusterTemplate
    :rtype: ClusterTemplateRender

    FIXME this shouldn't be a public method
    """

    t_name_to_file_name_dct = _get_cluster_files(cluster_model_dir)
    cluster_tmpls = [ClusterTemplate(name, _template_file_to_str(file_name)) for name, file_name in t_name_to_file_name_dct.iteritems()]

    return ClusterTemplateRender(cluster_tmpls)


def load_installed_cluster_templates_by_name(name):
    path = os.path.join(_PACKAGE_DIR, 'pbsmrtpipe', 'cluster_templates', name)
    return ClusterTemplateRender.from_dir(path)


def load_installed_cluster_templates_by_module_name(module_name):
    """Load cluster templates from a python module root. Example "pbsmrtpipe.cluster_templates.sge")

    :rtype: ClusterTemplateRender
    """
    m = importlib.import_module(module_name)
    name = os.path.basename(os.path.dirname(m.__file__))
    return load_installed_cluster_templates_by_name(name)


def load_cluster_templates(path_or_module_name):
    """
    Load cluster templates from an abspath /path/to/sge or
    python module name "pbsmrtpipe.cluster_templates.sge"
    """
    if os.path.isdir(path_or_module_name):
        return load_cluster_templates_from_dir(path_or_module_name)
    elif isinstance(path_or_module_name, str):
        return load_installed_cluster_templates_by_module_name(path_or_module_name)
    else:
        raise ValueError("Unable to load cluster manager from '{p}'".format(p=path_or_module_name))


def load_installed_cluster_templates():
    """

    Returns a list of {root_tmpl_name: ClusterTemplateRender}

    """
    id_to_cluster_render = {}
    for x in os.listdir(CLUSTER_TEMPLATE_DIR):
        d = os.path.join(CLUSTER_TEMPLATE_DIR, x)
        if os.path.isdir(d):
            cluster_render = ClusterTemplateRender.from_dir(d)
            id_to_cluster_render[x] = cluster_render

    return id_to_cluster_render


class ClusterTemplate(object):

    def __init__(self, name, template_str):
        if not Constants.is_valid(name):
            raise ValueError("Unsupported value '{x}'. Supported types {t}".format(x=name, t=Constants.all()))

        self.name = name
        # this is the actual template string, *NOT* the file
        self.template_str = template_str

    def __str__(self):
        return "{s}".format(s=self.template_str)

    def __repr__(self):
        return "" .join(["<", str(self), ">"])

    def __hash__(self):
        return hash(self.name)


class ClusterTemplateRender(object):
    # TODO FIXME This should be refactored to just use functions.

    def __init__(self, templates):
        assert all([isinstance(t, ClusterTemplate) for t in templates])
        self._templates = {t.name: t for t in templates}

    def __repr__(self):
        t = self._templates[Constants.START]
        i = t.template_str
        _d = dict(k=self.__class__.__name__, i=i)
        return "<{k} start:'{i}' >".format(**_d)

    @property
    def cluster_templates(self):
        return self._templates

    def _template_names(self):
        return self._templates.keys()

    def _get_template_by_name(self, template_name):
        return self._templates[template_name]

    @staticmethod
    def from_dir(cluster_manager_dir):
        return load_cluster_templates_from_dir(cluster_manager_dir)

    def render(self, template_name, shell_script, job_id, stdout=None, stderr=None, nproc=1, extras=None):
        """
        :param template_name: (str) name of template type (e.g., start, stop)
        :param shell_script: (str) path to shell script
        :param job_id: (int, str) For SGE, the job id will be correct format
        :param stdout: (str) path to stdout
        :param stderr: (str) path to stderr
        :param nproc: (int) number of processors to use
        :param extras: (None, str) extra options (e.g., '-l h_rt=24:0:0')
        :return: (str) qsub command
        """
        extras_types = (types.NoneType, basestring)
        if not isinstance(extras, extras_types):
            raise TypeError("'extras' expected types {x} got {t}.".format(t=type(extras), x=extras_types))

        if template_name not in self._template_names():
            raise ValueError("Unable to find template name '{t}'".format(t=template_name))

        # Job id's for SGE must begin with a digit
        if isinstance(job_id, int):
            job_id = "j" + str(job_id)
            warnings.warn("SGE requires the job id to begin with a character. Changing job to {i}".format(i=job_id))

        cluster_template = self._get_template_by_name(template_name)

        t = Template(cluster_template.template_str)
        # log.info(t)
        # log.info(t.template)

        extras_str = "" if extras is None else extras
        d = dict(CMD=shell_script, JOB_ID=job_id,
                 STDOUT_FILE=stdout,
                 STDERR_FILE=stderr,
                 EXTRAS=extras_str,
                 NPROC=str(nproc))

        # log.info(d)
        s = t.substitute(**d)

        return s


def validate_cluster_manager(cluster_manager):
    """
    Returns the cluster manager directory to templates.

    If a relative path is given 'mytemplate', the code will look for
    pbsmrtpipe/cluster_templates/mytemplate

    The required template files (e.g., start, stop) are also validated.

    This should be removed.
    """
    e_msg = "Unable to load cluster templates from {d}".format(d=cluster_manager)
    if os.path.isabs(cluster_manager):
        if os.path.isdir(cluster_manager):
            _ = load_cluster_templates_from_dir(cluster_manager)
            return cluster_manager

    resolved_cluster_manager = os.path.join(CLUSTER_TEMPLATE_DIR, cluster_manager)

    if os.path.exists(resolved_cluster_manager):
        if os.path.isdir(resolved_cluster_manager):
            _ = load_cluster_templates_from_dir(resolved_cluster_manager)
            return resolved_cluster_manager
    e_msg += " or from {r}".format(r=resolved_cluster_manager)
    raise IOError(e_msg)
