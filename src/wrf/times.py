from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .util import extract_times


def get_times(wrfnc, timeidx=0):
    return extract_times(wrfnc, timeidx)
