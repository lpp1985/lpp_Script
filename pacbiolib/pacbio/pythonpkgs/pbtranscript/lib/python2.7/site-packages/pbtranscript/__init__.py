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

"""Define version of pbtranscript."""
from __future__ import absolute_import


# $Date: 2015/03/21 $
# $Revision: Update version number from 0.2.0 to 1.0.0 $


_changelist = "$Change: 177900 $"


def _get_changelist(perforce_str):
    """Extract change list from p4 str"""
    import re
    rx = re.compile(r'Change: (\d+)')
    match = rx.search(perforce_str)
    if match is None:
        v = 'UnknownChangelist'
    else:
        try:
            v = int(match.group(1))
        except (TypeError, IndexError):
            v = "UnknownChangelist"
    return v


def get_changelist():
    """Return changelist"""
    return _get_changelist(_changelist)


VERSION = (1, 0, 0, get_changelist())


def get_version():
    """Return the version as a string. "O.7"

    This uses a major.minor

    Each python module of the system (e.g, butler, detective, siv_butler.py)
    will use this version +  individual changelist. This allows top level
    versioning, and sub-component to be versioned based on a p4 changelist.

    .. note:: This should be improved to be compliant with PEP 386.
    """
    return ".".join([str(i) for i in VERSION])
