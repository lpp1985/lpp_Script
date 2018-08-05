#################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

# Author: David Alexander, Jim Drake

import cProfile, logging, os.path
from multiprocessing import Process
from threading import Thread
from .options import options
from .reference import windowToString
from .io.utils import loadCmpH5, loadBam

class Worker(object):
    """
    Base class for compute worker that read reference coordinates
    from the task queue, perform variant calling, then push results
    back to another queue, to be written to a GFF file by a collector.

    All tasks that are O(genome length * coverage depth) should be
    distributed to compute workers, leaving the collector
    worker only O(genome length) work to do.
    """
    def __init__(self, workQueue, resultsQueue, algorithmConfig):
        self._workQueue = workQueue
        self._resultsQueue = resultsQueue
        self._algorithmConfig = algorithmConfig

    def _run(self):
        if options.usingBam:
            self._inAlnFile = loadBam(options.inputFilename, options.referenceFilename)
        else:
            self._inAlnFile = loadCmpH5(options.inputFilename, options.referenceFilename,
                                        disableChunkCache=options.disableHdf5ChunkCache)
        self.onStart()

        while True:
            datum = self._workQueue.get()
            if datum is None:
                # Sentinel indicating end of input.  Place a sentinel
                # on the results queue and end this worker process.
                self._resultsQueue.put(None)
                break
            else:
                if datum.hasCoverage:
                    msg = "%s received work unit, coords=%s"
                else:
                    msg = "%s received work unit, coords=%s (inadequate coverage)"
                logging.debug(msg % (self.name, windowToString(datum.window)))

                result = self.onChunk(datum)
                self._resultsQueue.put(result)

        self.onFinish()


    def run(self):
        if options.pdb:
            import ipdb
            with ipdb.launch_ipdb_on_exception():
                self._run()

        elif options.doProfiling:
            cProfile.runctx("self._run()",
                            globals=globals(),
                            locals=locals(),
                            filename=os.path.join(options.temporaryDirectory,
                                                  "profile-%s.out" % (self.name)))
        else:
            self._run()

    #==
    # Begin overridable interface
    #==

    def onStart(self):
        pass

    def onChunk(self, workChunk):
        """
        This function is the heart of the matter.

        workChunk -> result
        """
        pass

    def onFinish(self):
        pass

class WorkerProcess(Worker, Process):
    """Worker that executes as a process."""
    def __init__(self, *args):
        Process.__init__(self)
        super(WorkerProcess,self).__init__(*args)
        self.daemon = True

class WorkerThread(Worker, Thread):
    """Worker that executes as a thread (for debugging purposes only)."""
    def __init__(self, *args):
        Thread.__init__(self)
        super(WorkerThread,self).__init__(*args)
        self.daemon = True
        self.exitcode = 0
