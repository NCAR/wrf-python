from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

#from .extension import computedbz,computetk
from .extension import _dbz, _tk
from .constants import Constants
from .util import extract_vars
from .metadecorators import copy_and_set_metadata


@copy_and_set_metadata(copy_varname="T", name="dbz", 
                       description="radar reflectivity",
                       units="dBz")
def get_dbz(wrfnc, timeidx=0, method="cat", 
            squeeze=True, cache=None, meta=True, _key=None,
            use_varint=False, use_liqskin=False):
    """ Return the dbz
    
    use_varint - do variable intercept (if False, constants are used.  Otherwise, 
    intercepts are calculated using a technique from Thompson, Rasmussen, 
    and Manning (2004, Monthly Weather Review, Vol. 132, No. 2, pp. 519-542.)
    
    use_liqskin - do liquid skin for snow (frozen particles above freezing scatter
    as liquid)
    
    """
    varnames = ("T", "P", "PB", "QVAPOR", "QRAIN")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    qr = ncvars["QRAIN"]
    
    try:
        snowvars = extract_vars(wrfnc, timeidx, "QSNOW", 
                                method, squeeze, cache, meta=False,
                                _key=_key)
    except KeyError:
        qs = np.zeros(qv.shape, qv.dtype)
    else:
        qs = snowvars["QSNOW"]
    
    try:
        graupvars = extract_vars(wrfnc, timeidx, "QGRAUP", 
                                 method, squeeze, cache, meta=False,
                                 _key=_key)
    except KeyError:
        qg = np.zeros(qv.shape, qv.dtype)
    else:
        qg = graupvars["QGRAUP"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    # If qsnow is not all 0, set sn0 to 1
    sn0 = 1 if qs.any() else 0
    ivarint = 1 if use_varint else 0
    iliqskin = 1 if use_liqskin else 0
    
    return _dbz(full_p, tk, qv, qr, qs, qg, sn0, ivarint, iliqskin)


@copy_and_set_metadata(copy_varname="T", name="max_dbz", 
                       remove_dims=("bottom_top",),
                       description="maximum radar reflectivity",
                       units="dBz",
                       MemoryOrder="XY")
def get_max_dbz(wrfnc, timeidx=0, method="cat", 
                squeeze=True, cache=None, meta=True, _key=None,
                do_varint=False, do_liqskin=False):
    return np.amax(get_dbz(wrfnc, timeidx, method, squeeze, cache, meta,
                           _key, do_varint, do_liqskin), 
                  axis=-3)

