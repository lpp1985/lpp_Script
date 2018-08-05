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

"""Convert an ice CCS FASTA file to a FASTQ file."""

import logging
import sys
import os.path as op
from pbcore.util.ToolRunner import PBToolRunner
from pbtranscript.__init__ import get_version
from pbtranscript.ice.IceUtils import ice_fa2fq


def add_ice_fa2fq_arguments(parser):
    """Add arguments for ice_fa2fq.py"""
    parser.add_argument("in_fa", type=str,
                        help="An input FASTA file containing ice ccs reads " +
                             "(e.g., isoseq_flnc.fasta).")
    parser.add_argument("ccs_fofn", type=str,
                        help="An input ccs.h5 or ccs FOFN file.")
    parser.add_argument("out_fq", type=str, default=None,
                        help="An output FASTQ file containing ice ccs reads " +
                             "with Quality Values.")
    return parser


class IceFAToFQRunner(PBToolRunner):

    """ice_fa2fq runner."""

    def __init__(self):
        desc = "Convert ice ccs reads in a FASTA file to a FASTQ file, " + \
               "QVs will be read from ccs.h5 files in a FOFN."
        super(IceFAToFQRunner, self).__init__(desc)
        self.parser = add_ice_fa2fq_arguments(self.parser)

    def getVersion(self):
        """Return version string."""
        return get_version()

    def validate_inputs(self, in_fa, ccs_fofn):
        """Validate inputs."""
        if not op.exists(in_fa):
            raise IOError("Input fasta {f} does not exist.".
                          format(f=in_fa))
        if not op.exists(ccs_fofn):
            raise IOError("Input ccs_fofn {f} does not exist.".
                          format(f=ccs_fofn))

    def cmd_str(self, in_fa, ccs_fofn, out_fq):
        """Return a cmd string."""
        cmd = "ice_fa2fq.py {in_fa} {ccs_fofn} {out_fq} ".\
              format(in_fa=in_fa, ccs_fofn=ccs_fofn, out_fq=out_fq)
        return cmd

    def run(self):
        """Execute ice_fa2fq.py."""
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=get_version()))
        cmd_str = ""
        try:
            args = self.args
            in_fa, ccs_fofn, out_fq = args.in_fa, args.ccs_fofn, \
                                      args.out_fq

            self.validate_inputs(in_fa=in_fa,
                                 ccs_fofn=ccs_fofn)

            cmd_str = self.cmd_str(in_fa=in_fa, ccs_fofn=ccs_fofn,
                                   out_fq=out_fq)

            ice_fa2fq(in_fa=in_fa, ccs_fofn=ccs_fofn, out_fq=out_fq)

        except:
            logging.exception("Exiting {cmd} with return code 1.".
                              format(cmd=cmd_str))
            return 1
        return 0


def main():
    """Main function."""
    runner = IceFAToFQRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
