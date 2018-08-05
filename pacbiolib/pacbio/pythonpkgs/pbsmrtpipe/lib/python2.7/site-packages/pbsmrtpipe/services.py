"""Utils for Updating state/progress and results to WebServices"""

# keeping this for backward compatibility
from pbcommand.services import ServiceAccessLayer
# These are hidden methods for now
from pbcommand.services.service_access_layer import (log_pbsmrtpipe_progress,
                                                     add_datastore_file, LogLevels)
