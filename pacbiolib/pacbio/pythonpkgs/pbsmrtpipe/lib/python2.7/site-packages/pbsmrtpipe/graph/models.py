import logging
import json
import datetime

from pbcommand.models import FileType

from pbsmrtpipe.models import (TaskStates, ToolContractMetaTask,
                               MetaTask, MetaScatterTask,
                               ScatterToolContractMetaTask, MetaGatherTask,
                               GatherToolContractMetaTask)
from pbsmrtpipe.pb_io import strip_entry_prefix
from pbsmrtpipe.utils import validate_type_or_raise

log = logging.getLogger(__name__)


class DotShapeConstants(object):

    """Commonly used shapes so I don't have to keep looking this up"""
    ELLIPSE = "ellipse"
    RECTANGLE = "rectangle"
    OCTAGON = "octagon"
    TRIPLE_OCTAGON = 'tripleoctagon'
    # Use for scatter task
    TRIANGLE_UP = 'triangle'
    # Use for gather
    TRIANGLE_DOWN = 'invtriangle'

    DIAMOND = 'diamond'
    PARALLELOGRAM = 'parallelogram'


class DotColorConstants(object):

    """Colors so I don't have to keep looking this up

    http://www.graphviz.org/doc/info/colors.html
    """
    AQUA = "aquamarine"
    AQUA_DARK = 'aquamarine3'
    CYAN = 'cyan'
    RED = "red"
    BLUE = 'blue'
    ORANGE = 'orange'
    WHITE = 'azure'
    PURPLE = 'mediumpurple'
    PURPLE_DARK = 'mediumpurple4'
    GREY = 'grey'


class DotStyleConstants(object):
    DOTTED = 'dotted'
    FILLED = 'filled'


class Constants(object):
    TASK_FAILED_COLOR = DotColorConstants.RED


class ConstantsNodes(object):
    FILE_ATTR_IS_RESOLVED = 'is_resolved'
    FILE_ATTR_PATH = 'path'
    FILE_ATTR_RESOLVED_AT = 'resolved_at'

    TASK_ATTR_STATE = 'state'
    TASK_ATTR_NPROC = 'nproc'
    TASK_ATTR_ROPTS = 'ropts'
    TASK_ATTR_CMDS = 'cmds'
    TASK_ATTR_EMESSAGE = "error_message"
    TASK_ATTR_RUN_TIME = 'run_time'
    TASK_ATTR_UPDATED_AT = 'updated_at'
    TASK_ATTR_CREATED_AT = 'created_at'

    # label for identifiying tasks that have a companion chunk task
    # defined by the chunk operator
    TASK_ATTR_IS_CHUNKABLE = 'is_chunkable'

    # Was the chunked applied to the scattered chunk and the new
    # chunked N-tasks were created
    TASK_ATTR_WAS_CHUNKED = 'was_chunked'

    # Was the Gather'ing process applied to the
    TASK_ATTR_WAS_GATHERED = "was_gathered"

    # Label scatterable TaskBindingNode
    TASK_ATTR_IS_SCATTERABLE = "is_scatterable"

    TASK_ATTR_COMPANION_CHUNK_TASK_TYPE_ID = "companion_chunk_task_id"
    # Chunk operator
    TASK_ATTR_OPERATOR_ID = 'operator_id'
    # Chunk Group, applies to
    TASK_ATTR_CHUNK_GROUP_ID = "chunk_group_id"

    # Chunk keys to store on TaskScatterBindingNode
    TASK_ATTR_CHUNK_KEYS = "chunk_keys"


class _NodeLike(object):

    """Base Graph Node type"""
    NODE_ATTRS = {}


class _NodeEqualityMixin(object):

    def __repr__(self):
        return ''.join(['<', str(self), '>'])

    def __str__(self):
        return "{k}_{i}".format(k=self.__class__.__name__, i=self.idx)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.idx == other.idx:
                return True
        return False

    def __ne__(self, other):
        return not self == other


class _DotAbleMixin(object):
    DOT_SHAPE = DotShapeConstants.ELLIPSE
    DOT_COLOR = DotColorConstants.WHITE


class _FileLike(_NodeLike):
    # Attributes initialized at the graph level
    NODE_ATTRS = {ConstantsNodes.FILE_ATTR_IS_RESOLVED: False,
                  ConstantsNodes.FILE_ATTR_PATH: None,
                  ConstantsNodes.FILE_ATTR_RESOLVED_AT: None}


class _TaskLike(_NodeLike):

    """Attributes initialized at the graph level


    IS_CHUNKABLE => the chunk operator has found a companion scatter task for the original task
    IS_CHUNK_RUNNING => the companion scatter task is running
    WAS_CHUNKED => the scatter/chunking process has been applied to the task
    scatter/chunking are sloppy equivalents
    All of the chunked
    """
    NODE_ATTRS = {ConstantsNodes.TASK_ATTR_STATE: TaskStates.CREATED,
                  ConstantsNodes.TASK_ATTR_ROPTS: {},
                  ConstantsNodes.TASK_ATTR_NPROC: 1,
                  ConstantsNodes.TASK_ATTR_CMDS: [],
                  ConstantsNodes.TASK_ATTR_RUN_TIME: None,
                  ConstantsNodes.TASK_ATTR_CREATED_AT: lambda: datetime.datetime.now(),
                  ConstantsNodes.TASK_ATTR_UPDATED_AT: lambda: datetime.datetime.now(),
                  ConstantsNodes.TASK_ATTR_EMESSAGE: None,
                  ConstantsNodes.TASK_ATTR_IS_CHUNKABLE: False,
                  ConstantsNodes.TASK_ATTR_WAS_CHUNKED: False,
                  ConstantsNodes.TASK_ATTR_COMPANION_CHUNK_TASK_TYPE_ID: None,
                  ConstantsNodes.TASK_ATTR_OPERATOR_ID: None,
                  ConstantsNodes.TASK_ATTR_CHUNK_GROUP_ID: None,
                  ConstantsNodes.TASK_ATTR_WAS_GATHERED: False,
                  ConstantsNodes.TASK_ATTR_CHUNK_KEYS: []
                  }


class _ChunkLike(object):
    # Must have self.chunk_id
    pass


class EntryPointNode(_NodeEqualityMixin, _DotAbleMixin, _TaskLike):
    # this is like a Task
    # This abstraction needs to be deleted
    NODE_ATTRS = {ConstantsNodes.FILE_ATTR_IS_RESOLVED: False,
                  ConstantsNodes.FILE_ATTR_PATH: None,
                  ConstantsNodes.TASK_ATTR_RUN_TIME: 1}

    DOT_COLOR = DotColorConstants.PURPLE
    DOT_SHAPE = DotShapeConstants.DIAMOND

    def __init__(self, idx, file_klass):
        """

        :type idx: str
        :type file_klass: FileType

        :param idx:
        :param file_klass:
        :return:
        """
        self.idx = idx
        # FileType instance
        self.file_klass = file_klass
        # adding this to consistent with the Task data model
        self.instance_id = 0

    def __str__(self):
        _d = dict(k=self.__class__.__name__, i=self.idx, f=self.file_klass.file_type_id)
        return "{k} {i} {f}".format(**_d)


class TaskBindingNode(_NodeEqualityMixin, _DotAbleMixin, _TaskLike):

    """ Standard base Task Node """
    DOT_COLOR = DotColorConstants.AQUA
    DOT_SHAPE = DotShapeConstants.OCTAGON

    def __init__(self, meta_task, instance_id):
        """

        :type meta_task: MetaTask
        :type instance_id: int

        :param meta_task:
        :param instance_id:
        """

        self.meta_task = validate_type_or_raise(meta_task, (MetaTask, ToolContractMetaTask))
        self.instance_id = validate_type_or_raise(instance_id, int)

    @property
    def url(self):
        # for backwards compatible with the driver
        return str(self)

    @property
    def idx(self):
        return self.meta_task.task_id

    @property
    def nid(self):
        return repr(self)

    def __str__(self):
        _d = dict(k=self.__class__.__name__, i=self.idx, n=self.instance_id)
        return "{k} {i}-{n}".format(**_d)


class TaskChunkedBindingNode(TaskBindingNode):

    """Chunked "instances" of a Task node, must have chunk_id"""

    DOT_SHAPE = DotShapeConstants.TRIPLE_OCTAGON
    DOT_COLOR = DotColorConstants.AQUA_DARK

    def __init__(self, meta_task, instance_id, chunk_id, chunk_group_id, operator_id):
        super(TaskChunkedBindingNode, self).__init__(meta_task, instance_id)
        # Chunk Operator Id
        self.operator_id = operator_id
        self.chunk_id = chunk_id
        # str(uuid)
        self.chunk_group_id = chunk_group_id

    def __str__(self):
        _d = dict(k=self.__class__.__name__,
                  i=self.idx, n=self.instance_id, c=self.chunk_id, u=self.chunk_group_id)
        return "{k}-{c} {i}-{n} group {u}".format(**_d)


class TaskScatterBindingNode(TaskBindingNode):

    """Scattered Task that produces a chunk.json output file type

    This will have a 'companion' task that shares the same input file type signature
    """

    DOT_SHAPE = DotShapeConstants.OCTAGON
    DOT_COLOR = DotColorConstants.ORANGE

    def __init__(self, scatter_meta_task, original_nid, original_task_type_id, instance_id, chunk_group_id):
        validate_type_or_raise(scatter_meta_task, (MetaScatterTask, ScatterToolContractMetaTask))
        super(TaskScatterBindingNode, self).__init__(scatter_meta_task, instance_id)
        # Keep track of the original task type that was chunked
        self.original_task_id = original_task_type_id
        self.original_nid = original_nid
        self.chunk_group_id = chunk_group_id


class TaskGatherBindingNode(TaskBindingNode):

    """Gathered Task node. Consumes a gathered chunk.json and emits a single
    file type
    """
    DOT_SHAPE = DotShapeConstants.OCTAGON
    DOT_COLOR = DotColorConstants.GREY

    def __init__(self, meta_task, instance_id, chunk_key):
        validate_type_or_raise(meta_task, (MetaGatherTask, GatherToolContractMetaTask))
        super(TaskGatherBindingNode, self).__init__(meta_task, instance_id)
        # Keep track of the chunk_key that was passed to the exe.
        # Perhaps this should be in the meta task instance?
        self.chunk_key = chunk_key


class _BindingFileNode(_NodeEqualityMixin, _DotAbleMixin, _FileLike):
    # Grab from meta task
    ATTR_NAME = "input_types"
    # Used as a label in dot
    DIRECTION = "IN"

    DOT_COLOR = DotColorConstants.WHITE
    DOT_SHAPE = DotShapeConstants.ELLIPSE

    def __init__(self, meta_task, instance_id, index, file_type_instance):
        """

        :type meta_task: MetaTask
        :type instance_id: int
        :type index: int
        :type file_type_instance: FileType
        """

        # FileTypes are unique by (In/Out, file-type, instance-id)
        # Trying to use int which are more friendly and make graph easier
        # to view, but have to manually assigned in a centralized place

        # This was not a great design choice
        self.meta_task = validate_type_or_raise(meta_task, MetaTask)

        # task index (int)
        self.instance_id = validate_type_or_raise(instance_id, int)

        # positional index of input/output
        self.index = validate_type_or_raise(index, int)

        # this is a little odd. The input/output type are not necessarily identical
        self.file_klass = validate_type_or_raise(file_type_instance, FileType)

    @property
    def task_instance_id(self):
        # this is the {file klass}-{Instance id}
        types_ = getattr(self.meta_task, self.__class__.ATTR_NAME)
        file_type = types_[self.index]
        return "-".join([file_type.file_type_id, str(self.instance_id)])

    @property
    def idx(self):
        # the fundamental id used in the graph
        return "{n}.{i}".format(n=self.task_instance_id, i=self.index)

    def __str__(self):
        _d = dict(k=self.__class__.__name__,
                  d=self.__class__.DIRECTION,
                  i=self.idx,
                  f=self.file_klass.file_type_id,
                  n=self.meta_task.task_id,
                  m=self.instance_id)
        return "{k} {n}-{m} {i}".format(**_d)


class BindingInFileNode(_BindingFileNode):
    DIRECTION = "in"
    # Grab from meta task
    ATTR_NAME = "input_types"

    DOT_SHAPE = DotShapeConstants.ELLIPSE


class BindingChunkInFileNode(BindingInFileNode):

    """Chunked file type that belongs to a Chunked Group

    This should always be generated from a Chunk.json
    """

    def __init__(self, meta_task, instance_id, index, file_type_instance, chunk_id, chunk_group_id):
        super(BindingChunkInFileNode, self).__init__(meta_task, instance_id, index, file_type_instance)
        self.chunk_id = chunk_id
        self.chunk_group_id = chunk_group_id


class BindingOutFileNode(_BindingFileNode):
    DIRECTION = "out"
    ATTR_NAME = "output_types"

    DOT_SHAPE = DotShapeConstants.OCTAGON


class BindingChunkOutFileNode(BindingOutFileNode):

    def __init__(self, meta_task, instance_id, index, file_type_instance, chunk_id, chunk_group_id):
        super(BindingChunkOutFileNode, self).__init__(meta_task, instance_id, index, file_type_instance)
        self.chunk_id = chunk_id
        self.chunk_group_id = chunk_group_id


class EntryOutBindingFileNode(_NodeEqualityMixin, _DotAbleMixin, _FileLike):
    DOT_SHAPE = DotShapeConstants.RECTANGLE
    DOT_COLOR = DotColorConstants.WHITE

    def __init__(self, entry_id, file_klass):
        self.entry_id = strip_entry_prefix(entry_id)
        # FileType instance
        self.file_klass = file_klass
        self.instance_id = 0
        self.index = 0
        self.direction = 'out'

    def __str__(self):
        _d = dict(k=self.__class__.__name__,
                  i=self.entry_id,
                  f=self.file_klass.file_type_id)
        return "{k} {f} {i}".format(**_d)

    @property
    def idx(self):
        _d = dict(n=self.entry_id, i=self.index)
        return "{n}.{i}".format(**_d)


VALID_FILE_NODE_CLASSES = (BindingInFileNode, BindingOutFileNode, EntryOutBindingFileNode)
VALID_TASK_NODE_CLASSES = (TaskBindingNode, EntryPointNode)
# FIXME
VALID_ALL_TASK_NODE_CLASSES = (TaskBindingNode, EntryPointNode, TaskChunkedBindingNode, TaskGatherBindingNode, TaskScatterBindingNode)
