"""Load Registered tasks, file types, pipelines or Chunk Operators.

Not Particularly thrilled with this model, howver, it is centralized.
"""

import os
import importlib
import logging
import functools
import warnings

from pbcommand.models.common import REGISTERED_FILE_TYPES
import pbsmrtpipe.constants as GlobalConstants
from pbcommand.pb_io import load_tool_contract_from

log = logging.getLogger(__name__)

# Loading Caches
# ToolContracts that are converted to MetaTasks (old model)
_REGISTERED_TOOL_CONTRACTS = None
_REGISTERED_PIPELINES = None
_REGISTERED_OPERATORS = None


def _load_all_tool_contracts_from(dir_name):
    # the old MetaTask model is making this
    # a bit convoluted
    import pbsmrtpipe.pb_io as IO
    mtasks = {}
    for file_name in os.listdir(dir_name):
        if file_name.endswith('.json'):
            f = os.path.join(dir_name, file_name)
            # sanity check using pbcommand
            _ = load_tool_contract_from(f)
            # Old layer to use MetaTask
            mtask = IO.tool_contract_to_meta_task_from_file(f)
            mtasks[mtask.task_id] = mtask

    return mtasks


def _load_all_tool_contracts(module_name, registered_tasks_d, filter_filename_func, processing_func):

    m = importlib.import_module(module_name)

    d = os.path.dirname(m.__file__)
    log.debug("Loading static meta tasks from {m}".format(m=d))

    for x in os.listdir(d):
        if filter_filename_func(x):
            json_file = os.path.join(d, x)
            try:
                meta_task = processing_func(json_file)
                log.debug(meta_task)
                if meta_task.task_id in registered_tasks_d:
                     log.warn("MetaTask {i} already loaded".format(i=meta_task.task_id))
                registered_tasks_d[meta_task.task_id] = meta_task
            except Exception as e:
                log.error("Failed loading Static Task from '{x}'".format(x=json_file))
                log.error(e.message)
                raise

    return registered_tasks_d


def _get_env_path_if_defined(env_var):
    """Get a Config env variable directory, or return None"""
    path = os.getenv(env_var)
    if path is not None:
        if os.path.isdir(path):
            return os.path.abspath(path)
        else:
            warnings.warn("Skipping loading contracts from {e} Enable to find {p}".format(e=env_var, p=path))

    return None


def load_all_tool_contracts():
    """
    This name is a bit of misnomer. This loads the TCs, then converts to MetaTask

    Loads all Tool Contracts.

    1. from pbsmrtpipe.registered_tool_contracts
    2. from pbsmrtpipe.regiesteried_tool_contracts_sa3
    3. from all json files in dir defined by by the env var PB_TC_DIR
    """

    import pbsmrtpipe.pb_io as IO

    # this is gross.
    global _REGISTERED_TOOL_CONTRACTS

    if _REGISTERED_TOOL_CONTRACTS is None:
        _REGISTERED_TOOL_CONTRACTS = {}

    def filter_by(name, path):
        return path.endswith(".json") and name in path

    filter_contracts = functools.partial(filter_by, "tool_contract")

    rtasks = _load_all_tool_contracts("pbsmrtpipe.registered_tool_contracts_sa3", _REGISTERED_TOOL_CONTRACTS, filter_contracts, IO.tool_contract_to_meta_task_from_file)
    rtasks = _load_all_tool_contracts("pbsmrtpipe.registered_tool_contracts", rtasks, filter_contracts, IO.tool_contract_to_meta_task_from_file)

    tc_path = _get_env_path_if_defined(GlobalConstants.ENV_TC_DIR)
    if tc_path is not None:
        tcs_mtasks = _load_all_tool_contracts_from(tc_path)
        rtasks.update(tcs_mtasks)

    return rtasks


def __load_chunk_operators_from_dir(path):
    import pbsmrtpipe.pb_io as IO

    operators = []

    for x in os.listdir(path):
        if x.endswith(".xml"):
            p = os.path.join(path, x)
            operator = IO.parse_operator_xml(p)
            operators.append(operator)

    return {op.idx: op for op in operators}


def _load_chunk_operators_from_env(env_name=GlobalConstants.ENV_CHK_OPT_DIR):
    dir_name = _get_env_path_if_defined(env_name)
    operators_d = {}
    if dir_name is not None:
        operators_d = __load_chunk_operators_from_dir(dir_name)
        log.debug("Loaded {o} operators from {d}".format(o=len(operators_d), d=dir_name))
    return operators_d


def _load_xml_chunk_operators_from_python_module_name(name):
    m = importlib.import_module(name)
    path = os.path.dirname(m.__file__)
    return __load_chunk_operators_from_dir(path)


def load_all_installed_chunk_operators():
    """:rtype dict[str, ChunkOperator]"""

    module_name = "pbsmrtpipe.chunk_operators"
    chunk_operator_env_name = GlobalConstants.ENV_CHK_OPT_DIR

    global _REGISTERED_OPERATORS
    if _REGISTERED_OPERATORS is None:
        d1 = _load_xml_chunk_operators_from_python_module_name(module_name)
        d2 = _load_chunk_operators_from_env(env_name=chunk_operator_env_name)
        # Values defined in the ENV will overwrite operators defined in python module loading
        d1.update(d2)
        _REGISTERED_OPERATORS = d1
        log.debug("Loaded {o} total chunk operators from module {m} ENV {e} (if defined)".format(o=len(d1), m=module_name, e=chunk_operator_env_name))

    return _REGISTERED_OPERATORS


def _load_pipelines_from_python_module_name(name):
    # FIXME This is terrible form. Need to update the loading to be configurable to
    # dynamically load pipelines
    # it's a dict and sometimes is a list of registered resources
    # m = importlib.import_module(name)

    global _REGISTERED_PIPELINES

    if _REGISTERED_PIPELINES is None:
        import pbsmrtpipe.pb_pipelines.pb_pipelines_dev
        import pbsmrtpipe.pb_pipelines.pb_pipelines_sa3
        import pbsmrtpipe.pb_pipelines.pb_pipelines_falcon

        from pbsmrtpipe.models import REGISTERED_PIPELINES
        _REGISTERED_PIPELINES = REGISTERED_PIPELINES

    return _REGISTERED_PIPELINES


def load_resolved_pipeline_template_jsons_from_dir(dir_name):
    """
    :rtype: list[Pipeline]
    """
    import pbsmrtpipe.pb_io as IO

    pipelines = []
    if os.path.exists(dir_name):
        for file_name in os.listdir(dir_name):
            if file_name.endswith(".json"):
                try:
                    p = IO.load_pipeline_template_from(os.path.join(dir_name, file_name))
                    pipelines.append(p)
                except Exception as e:
                    log.warn("Unable to load Resolved Pipeline Template from {}. {}".format(dir_name, str(e)))
    else:
        log.warn("Unable to load Resolved Pipeline Template from {}. Path does not exist.".format(dir_name))

    return pipelines


def _env_load_resolved_pipeline_template_json_from_env(env=GlobalConstants.ENV_PT_DIR):

    global _REGISTERED_PIPELINES

    dir_value = os.getenv(env)
    # print "ENV {}".format(dir_value)

    if dir_value:
        dir_names = dir_value.split(":")
        for dir_name in dir_names:
            # print "Loading from ", dir_name
            # this needs to be able to reference existing pipeline templates
            resolved_pipeline_templates = load_resolved_pipeline_template_jsons_from_dir(dir_name)
            for resolved_pipeline_template in resolved_pipeline_templates:
                _REGISTERED_PIPELINES[resolved_pipeline_template.idx] = resolved_pipeline_template

    return _REGISTERED_PIPELINES


def load_all_installed_pipelines():
    _ = _load_pipelines_from_python_module_name("")
    _ = _env_load_resolved_pipeline_template_json_from_env(env=GlobalConstants.ENV_PT_DIR)
    return _REGISTERED_PIPELINES


def load_all_registered_file_types():
    return REGISTERED_FILE_TYPES


def load_all():
    """
    Load all resources and return a tuple of (MetaTasks, FileTypes, ChunkOperators, Pipelines)

    :note: This will only be loaded once and cached.
    """
    meta_tasks = load_all_tool_contracts()
    operators = load_all_installed_chunk_operators()
    pipelines = load_all_installed_pipelines()

    from pbsmrtpipe.core import REGISTERED_FILE_TYPES
    return meta_tasks, REGISTERED_FILE_TYPES, operators, pipelines


def load_and_validate_chunk_operators():
    from .models import validate_operator

    rtasks = load_all_tool_contracts()
    chunk_operators = load_all_installed_chunk_operators()
    for operator_id, chunk_operator in chunk_operators.iteritems():
        # this will raise if invalid
        validate_operator(chunk_operator, rtasks)
