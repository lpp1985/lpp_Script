
class RequiredExeNotFoundError(Exception):

    """External required subprocess is not found"""
    pass


class WorkflowBaseException(Exception):
    pass


class TaskIdNotFound(KeyError):

    """Task id not found in Registered Tasks"""
    pass


class MalformedBindingGraphError(ValueError):

    """Binding graph is invalid form"""
    pass


class InvalidEntryPointError(ValueError):

    """RequiredEntry point has not be resolved"""
    pass


class InvalidDependencyInjectError(ValueError):

    """General exception for handling DI errors"""
    pass


class MalformedDependencyInjectionFileMetadataError(ValueError):

    """Unsupported or malformed file metadata"""


class MalformedBindingError(ValueError):

    """Binding to a task is incompatible"""


class MalformedBindingStrError(ValueError):

    """Binding str is not well-formed"""


class MalformedEntryStrError(ValueError):

    """Entry point format is malformed"""


class MalformedMetaTaskError(ValueError):

    """MetaTask was not correctly formed or is missing class vars"""


class BindingFileTypeIncompatiblyError(TypeError):

    """
    FileType descriptions in Task are bound to FileTypes that
    are incompatible"""
    pass


class PipelineTemplateIdNotFoundError(KeyError):

    """Unable to find pipeline template id in registry"""
    pass


class MalformedPipelineError(ValueError):

    """Pipeline definition is not valid"""
    pass


class MalformedChunkOperatorError(ValueError):

    """Invalid Chunk Operator Definition"""
    pass


class MalformedChunkKeyError(ValueError):

    """Chunk Key does NOT adhere to the spec"""
    pass


class PipelineRuntimeError(RuntimeError):

    """A generic error occurred when running the pipeline"""


class PipelineRuntimeKeyboardInterrupt(KeyboardInterrupt):

    """User killed the workflow"""


class WorkflowError(WorkflowBaseException):
    pass


class TaskExecutionError(WorkflowBaseException):

    """Task failed during pipeline running"""
    pass


class BaseChunkError(WorkflowBaseException):
    pass


class ChunkScatteringError(BaseChunkError):

    """Unable to create or process output of scattered Task"""
    pass


class ChunkGatheringError(BaseException):

    """Failed to Gather Chunked tasks outputs"""
    pass


class TaskChunkingError(BaseChunkError):

    """General error for creating chunked task instances"""
    pass
