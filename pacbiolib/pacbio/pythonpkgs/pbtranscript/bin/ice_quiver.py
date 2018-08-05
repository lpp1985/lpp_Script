#! python
###############################################################################
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
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

"""
Overview:
    pbtranscript cluster contains two main components:
    * (1) ICE (iterative clustering and error correction) to predict
      unpolished consensus isoforms.
    * (2) Polish, to use nfl reads and quiver to polish those predicted
      unpolished isoforms. Polish contains three steps:
      + (2.1) IceAllPartials (ice_partial.py all)
              Align and assign nfl reads to unpolished isoforms, and
              save results to a pickle file.
      + (2.2) IceQuiver (ice_quiver.py i and ice_quiver.py merge)
              Call quiver to polish each isoform based on alignments
              created by mapping its associated fl and nfl reads to
              this isoform.
      + (2.3) IceQuiverPostprocess (ice_quiver.py postprocess)
              Collect and post process quiver results, and classify
              HQ/LQ isoforms.

    In order to handle subtasks by SMRTPipe instead of pbtranscript
    itself, we will refactor the polish phase including
    (2.1) (2.2) and (2.3). The refactor of (2.1) is described in
    ice_partial.py.

    (2.2) IceQuiver will be refactored to
       + (2.2.1) IceQuiverI (ice_quiver.py i)
                 Split all unpolished isoforms into N chunks and
                 call Quiver to polish isoforms of the i-th chunk
                 at a time
       + (2.2.2) IceQuiverMerge (ice_quiver.py merge)
                 When all splitted quiver jobs are done,
                 collect all submitted jobs and save to
                 root_dir/log/submitted_quiver_jobs.txt
    (2.3) IceQuiverPostprocess will be renamed from ice_post_quiver.py to:
       + (2.3.1) ice_quiver.py postprocess

    e.g., IceQuiver (2.2) =
               ice_quiver_i.py root_dir 0 N   --> process the 0 th chunk
            +  ...
            +  ice_quiver_i.py root_dir N-1 N --> process the N-1 th chunk
            +  ice_quiver_merge.py root_dir N --> merge all polished isoforms
                                                  from N chunks

    Thus, the task of (2.2) + (2.3) can either be done by running:
        ice_quiver.py all ...
    or by
        running quiver on each chunk of consensus isoforms independently:
            ice_quiver.py i root_dir {i} N \
                          --bas_fofn=bas_fofn \
                          --use_sge? \
                          --max_sge_jobs=max_sge_jobs \
                          --unique_id=unique_id \
                          --quiver_nproc=quiver_nproc \
                          --blasr_nproc=blasr_nproc
            , for i = 0, ..., N-1
        and then collecting all polisehd consensus isoforms:
            ice_quiver.py merge root_dir N
        and finally post-processing polished isoforms to find LQ/HQ isoforms.
            ice_quiver.py postprocess root_dir

    Hierarchy:
        pbtranscript = iceiterative

        pbtranscript --quiver = iceiterative + \
                                ice_polish.py

        ice_polish.py =  ice_make_fasta_fofn.py + \
                         ice_partial.py all + \
                         ice_quiver.py all

        ice_partial.py all = ice_partial.py split + \
                             ice_partial.py i + \
                             ice_partial.py merge

        (ice_partial.py one --> only apply ice_partial on a given input fasta)

        ice_quiver.py all = ice_quiver.py i + \
                            ice_quiver.py merge + \
                            ice_quiver.py postprocess

Alternative way to call this script:
    python -m pbtranscript.ice_quiver
"""
import logging
import sys
import os.path as op
from pbcore.util.ToolRunner import PBMultiToolRunner
from pbtranscript.__init__ import get_version
from pbtranscript.PBTranscriptOptions import _wrap_parser
from pbtranscript.ClusterOptions import SgeOptions, \
    IceQuiverHQLQOptions
from pbtranscript.ice.IceQuiverAll import IceQuiverAll, \
    add_ice_quiver_all_arguments
from pbtranscript.ice.IceQuiverI import IceQuiverI, \
    add_ice_quiver_i_arguments
from pbtranscript.ice.IceQuiverMerge import IceQuiverMerge, \
    add_ice_quiver_merge_arguments
from pbtranscript.ice.IceQuiverPostprocess import \
    IceQuiverPostprocess, add_ice_quiver_postprocess_arguments


class IceQuiverRunner(PBMultiToolRunner):

    """ice_quiver runner, subcommands include 'all', 'i', 'merge' and
    'postprocess' """

    def __init__(self):
        desc = "Toolkit for assigning non-full-length reads to isoforms."
        super(IceQuiverRunner, self).__init__(desc)
        subparsers = self.subParsers

        parser = subparsers.add_parser('all',
                                       description=IceQuiverAll.desc)
        add_ice_quiver_all_arguments(_wrap_parser(parser))

        parser = subparsers.add_parser('i',
                                       description=IceQuiverI.desc)
        add_ice_quiver_i_arguments(parser)

        parser = subparsers.add_parser('merge',
                                       description=IceQuiverMerge.desc)
        add_ice_quiver_merge_arguments(parser)

        parser = subparsers.add_parser('postprocess',
                                       description=IceQuiverPostprocess.desc)
        add_ice_quiver_postprocess_arguments(parser)

    def getVersion(self):
        """Return version string."""
        return get_version()

    def run(self):
        """Execute ice_quiver.py all|i|merge|postprocess."""
        cmd = self.args.subCommand
        logging.info("Running {f} {cmd} v{v}.".format(f=op.basename(__file__),
                                                      cmd=cmd, v=get_version()))
        cmd_str = ""
        try:
            args = self.args
            obj = None
            if cmd == "all":
                sge_opts = SgeOptions(unique_id=args.unique_id,
                                      use_sge=args.use_sge,
                                      max_sge_jobs=args.max_sge_jobs,
                                      blasr_nproc=args.blasr_nproc,
                                      quiver_nproc=args.quiver_nproc)
                ipq_opts = IceQuiverHQLQOptions(
                    hq_isoforms_fa=args.hq_isoforms_fa,
                    hq_isoforms_fq=args.hq_isoforms_fq,
                    lq_isoforms_fa=args.lq_isoforms_fa,
                    lq_isoforms_fq=args.lq_isoforms_fq,
                    qv_trim_5=args.qv_trim_5,
                    qv_trim_3=args.qv_trim_3,
                    hq_quiver_min_accuracy=args.hq_quiver_min_accuracy)
                obj = IceQuiverAll(root_dir=args.root_dir,
                                   bas_fofn=args.bas_fofn,
                                   fasta_fofn=None,
                                   sge_opts=sge_opts,
                                   ipq_opts=ipq_opts,
                                   tmp_dir=args.tmp_dir)
            elif cmd == "i":
                sge_opts = SgeOptions(unique_id=args.unique_id,
                                      use_sge=args.use_sge,
                                      max_sge_jobs=args.max_sge_jobs,
                                      blasr_nproc=args.blasr_nproc,
                                      quiver_nproc=args.quiver_nproc)
                obj = IceQuiverI(root_dir=args.root_dir, i=args.i, N=args.N,
                                 bas_fofn=args.bas_fofn,
                                 fasta_fofn=None,
                                 sge_opts=sge_opts,
                                 tmp_dir=args.tmp_dir)
            elif cmd == "merge":
                obj = IceQuiverMerge(root_dir=args.root_dir, N=args.N)
            elif cmd == "postprocess":
                ipq_opts = IceQuiverHQLQOptions(
                    hq_isoforms_fa=args.hq_isoforms_fa,
                    hq_isoforms_fq=args.hq_isoforms_fq,
                    lq_isoforms_fa=args.lq_isoforms_fa,
                    lq_isoforms_fq=args.lq_isoforms_fq,
                    qv_trim_5=args.qv_trim_5,
                    qv_trim_3=args.qv_trim_3,
                    hq_quiver_min_accuracy=args.hq_quiver_min_accuracy)
                obj = IceQuiverPostprocess(root_dir=args.root_dir,
                                           ipq_opts=ipq_opts,
                                           use_sge=args.use_sge,
                                           quit_if_not_done=args.quit_if_not_done,
                                           summary_fn=args.summary_fn,
                                           report_fn=args.report_fn)
            else:
                raise ValueError("Unknown command passed to {f}: {cmd}.".
                                 format(f=op.basename(__file__), cmd=cmd))
            cmd_str = obj.cmd_str()
            logging.info("Running CMD: {cmd_str}".format(cmd_str=cmd_str))
            obj.run()
        except:
            logging.exception("Exiting {cmd_str} with return code 1.".
                              format(cmd_str=cmd_str))
            return 1
        return 0


def main():
    """Main function."""
    runner = IceQuiverRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
