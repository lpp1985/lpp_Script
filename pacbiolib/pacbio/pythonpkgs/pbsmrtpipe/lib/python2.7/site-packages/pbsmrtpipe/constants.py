import functools
# This is where all the rdf imports should be. Other modules should import them
# from here because of the plugin registry.
import os
import re
import socket


ENV_PRESET = 'PB_SMRTPIPE_XML_PRESET'
# Extra Directory for JSON Tool Contracts
ENV_TC_DIR = "PB_TOOL_CONTRACT_DIR"

# Extra Directory for JSON/Avro Pipeline Templates
ENV_PT_DIR = "PB_PIPELINE_TEMPLATE_DIR"

# Chunk Operators
ENV_CHK_OPT_DIR = "PB_CHUNK_OPERATOR_DIR"


PBSMRTPIPE_PID_KILL_FILE_SCRIPT = ".pbsmrtpipe-terminate.sh"

# Map of exception types to exit codes.
EXCEPTION_TO_EXIT_CODE = {KeyboardInterrupt: 7, IOError: 2, socket.error: 7}
TERM_FILE = ".TERMINATED"
EXIT_SUCCESS = 0
EXIT_FAILURE = 1
EXIT_TERMINATED = 7
# For Unknown error
DEFAULT_EXIT_CODE = EXIT_FAILURE

DEEP_DEBUG = False

# Global Env vars that are necessary to do *anything*. This is stupid.
SEYMOUR_HOME = 'SEYMOUR_HOME'

SLOG_PREFIX = 'status.'

DATASTORE_VERSION = "0.2.1"
CHUNK_API_VERSION = "0.1.0"

ENTRY_PREFIX = "$entry:"

# Generic Id
RX_TASK_ID = re.compile(r'^([A-z0-9_]*)\.tasks\.([A-z0-9_]*)$')
RX_PIPELINE_ID = re.compile(r'^([A-z0-9_]*)\.pipelines\.([A-z0-9_\.]*)$')
RX_FILE_TYPES_ID = re.compile(r'^([A-z0-9_]*)\.files\.([A-z0-9_\.]*)$')
RX_TASK_OPTION_ID = re.compile(r'^([A-z0-9_]*)\.task_options\.([A-z0-9_\.]*)')

# Bindings format
# Only Task includes the :0
RX_BINDING_TASK = re.compile(r'^([A-z0-9_]*)\.tasks\.([A-z0-9_]*):(\d*)$')
# Advanced bindings referencing an instance of a task
# {namespace}:tasks:{instance-id}:{in_out_index}
RX_BINDING_TASK_ADVANCED = re.compile(r'^([A-z0-9_]*)\.tasks\.([A-z0-9_]*):(\d*):(\d*)$')

# {pipeline_id}:{task_id}:{instance_id}
RX_BINDING_PIPELINE_TASK = re.compile(r'^([A-z0-9_]*).pipelines.([A-z0-9_]*):([A-z0-9_]*).tasks.(\w*):([0-9]*)$')
# {pipeline_id}:$entry:{entry_label}
RX_BINDING_PIPELINE_ENTRY = re.compile(r'^([A-z0-9_]*).pipelines.([A-z0-9_]*):\$entry:([A-z0-9_]*)$')
# Only Entry points
RX_ENTRY = re.compile(r'^\$entry:([A-z0-9_]*)$')
# to be consistent with the new naming scheme
RX_BINDING_ENTRY = re.compile(r'^\$entry:([A-z0-9_]*)$')

RX_VALID_BINDINGS = (RX_BINDING_PIPELINE_TASK,
                     RX_BINDING_PIPELINE_ENTRY,
                     RX_BINDING_TASK,
                     RX_BINDING_TASK_ADVANCED,
                     RX_BINDING_ENTRY)

# This should really use a semantic version lib
RX_VERSION = re.compile(r'(\d*).(\d*).(\d*)')

# Chunk Key $chunk.my_label_id
RX_CHUNK_KEY = re.compile(r'^\$chunk\.([A-z0-9_]*)')
RX_CHUNK_ID = re.compile(r'(^[A-z0-9_]*)')


TASK_MANIFEST_JSON = 'task-manifest.json'
RUNNABLE_TASK_JSON = "runnable-task.json"
TASK_MANIFEST_VERSION = '0.3.0'

RESOLVED_TOOL_CONTRACT_JSON = "resolved-tool-contract.json"
RESOLVED_TOOL_CONTRACT_AVRO = 'resolved-tool-contract.avro'
TOOL_CONTRACT_JSON = "tool-contract.json"

# ***** DEFAULT PIPELINE LEVEL OPTIONS ******
# Global hard limit on the maximum number of chunks per task are created
MAX_NCHUNKS = 128
MAX_NPROC = 16
MAX_TOTAL_NPROC = None
MAX_NWORKERS = 100
CHUNKED_MODE = False
# Only if the CLUSTER_MANAGER_DIR is defined
DISTRIBUTED_MODE = True
CLUSTER_MANAGER_DIR = None
TMP_DIR = os.getenv('TMP_DIR', '/tmp')
EXIT_ON_FAILIURE = False
DEBUG_MODE = False


class PacBioNamespaces(object):
    # File Types
    #PBSMRTPIPE_FILE_PREFIX = 'pbsmrtpipe.files'
    # NEW File Type Identifier style Prefix
    NEW_PBSMRTPIPE_FILE_PREFIX = "PacBio.FileTypes"
    # New DataSet Identifier Prefix
    DATASET_FILE_PREFIX = "PacBio.DataSet"
    # Task Ids
    PBSMRTPIPE_TASK_PREFIX = 'pbsmrtpipe.tasks'
    # Task Options
    PBSMRTPIPE_TASK_OPTS_PREFIX = 'pbsmrtpipe.task_options'
    # Workflow Level Options
    PBSMRTPIPE_OPTS_PREFIX = 'pbsmrtpipe.options'
    # Constants
    PBSMRTPIPE_CONSTANTS_PREFIX = 'pbsmrtpipe.constants'
    # Pipelines
    PBSMRTPIPE_PIPELINES = "pbsmrtpipe.pipelines"
    # OptionTypes (this should really be in pbcommand
    OPTION_TYPE = "pbsmrtpipe.option_types"


def __to_type(prefix, name):
    return ".".join([prefix, name])

to_constant_ns = functools.partial(__to_type, PacBioNamespaces.PBSMRTPIPE_CONSTANTS_PREFIX)
to_file_ns = functools.partial(__to_type, PacBioNamespaces.NEW_PBSMRTPIPE_FILE_PREFIX)
to_ds_ns = functools.partial(__to_type, PacBioNamespaces.DATASET_FILE_PREFIX)
to_task_option_ns = functools.partial(__to_type, PacBioNamespaces.PBSMRTPIPE_TASK_OPTS_PREFIX)
to_task_ns = functools.partial(__to_type, PacBioNamespaces.PBSMRTPIPE_TASK_PREFIX)
to_workflow_option_ns = functools.partial(__to_type, PacBioNamespaces.PBSMRTPIPE_OPTS_PREFIX)
to_pipeline_ns = functools.partial(__to_type, PacBioNamespaces.PBSMRTPIPE_PIPELINES)
to_opt_type_ns = functools.partial(__to_type, PacBioNamespaces.OPTION_TYPE)
