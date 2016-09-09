from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .util import extract_times


def get_times(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
              meta=True, _key=None):
    return extract_times(wrfnc, timeidx, method, squeeze, cache, 
                         meta=meta, do_xtime=False)


def get_xtimes(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
               meta=True, _key=None):
    return extract_times(wrfnc, timeidx, method, squeeze, cache, 
                         meta=meta, do_xtime=True)

