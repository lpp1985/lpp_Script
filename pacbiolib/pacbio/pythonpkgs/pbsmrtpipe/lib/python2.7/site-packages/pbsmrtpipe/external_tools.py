"""Externally/subprocess calls.

The dot to image should be done within pygraphviz if the dependency isn't too
difficult to consistent to install.
"""
import os
import logging
import functools

from pbsmrtpipe.exceptions import RequiredExeNotFoundError
from pbsmrtpipe.engine import backticks
from pbcommand.utils import which

log = logging.getLogger(__name__)

_SUPPORTED_IMAGE_TYPES = 'png svg eps'.split()
DOT_EXE = 'dot'


def _dot_to_image(image_type, dot_file, image_file):
    assert image_type.lower() in _SUPPORTED_IMAGE_TYPES

    if not os.path.exists(dot_file):
        raise IOError("Unable to find {f}".format(f=dot_file))

    if which(DOT_EXE) is None:
        raise RequiredExeNotFoundError("Unable to find required external exe '{x}'".format(x=DOT_EXE))

    cmd_str = "{e} -T{t} {i} -o {o}"
    d = dict(e=DOT_EXE, t=image_type, i=dot_file, o=image_file)
    cmd = cmd_str.format(**d)
    rcode, stdout, stderr, run_time = backticks(cmd)

    state = True if rcode == 0 else False
    return state

# For backward compatibility
dot_to_image = _dot_to_image

dot_file_to_png = functools.partial(_dot_to_image, 'png')
dot_file_to_svg = functools.partial(_dot_to_image, 'svg')
dot_file_to_eps = functools.partial(_dot_to_image, 'eps')
