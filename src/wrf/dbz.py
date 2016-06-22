from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as n

#from .extension import computedbz,computetk
from .extension import computedbz, _tk
from .constants import Constants
from .util import extract_vars
from .metadecorators import copy_and_set_metadata


@copy_and_set_metadata(copy_varname="T", name="dbz", 
                       description="radar reflectivity",
                       units="dBz")
def get_dbz(wrfnc, timeidx=0, method="cat", 
            squeeze=True, cache=None, meta=True,
            do_varint=False, do_liqskin=False):
    """ Return the dbz
    
    do_varint - do variable intercept (if False, constants are used.  Otherwise, 
    intercepts are calculated using a technique from Thompson, Rasmussen, 
    and Manning (2004, Monthly Weather Review, Vol. 132, No. 2, pp. 519-542.)
    
    do_liqskin - do liquid skin for snow (frozen particles above freezing scatter
    as liquid)
    
    """
    varnames = ("T", "P", "PB", "QVAPOR", "QRAIN")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    qr = ncvars["QRAIN"]
    
    try:
        snowvars = extract_vars(wrfnc, timeidx, "QSNOW", 
                                method, squeeze, cache, meta=False)
    except KeyError:
        qs = n.zeros(qv.shape, "float")
    else:
        qs = snowvars["QSNOW"]
    
    try:
        graupvars = extract_vars(wrfnc, timeidx, "QGRAUP", 
                                 method, squeeze, cache, meta=False)
    except KeyError:
        qg = n.zeros(qv.shape, "float")
    else:
        qg = graupvars["QGRAUP"]
    
    # If qsnow is all 0, set sn0 to 1
    sn0 = 0
    if (n.any(qs != 0)):
        sn0 = 1
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    ivarint = 0
    if do_varint:
        ivarint = 1
        
    iliqskin = 0
    if do_liqskin:
        iliqskin = 1
    
    return computedbz(full_p,tk,qv,qr,qs,qg,sn0,ivarint,iliqskin)


@copy_and_set_metadata(copy_varname="T", name="max_dbz", 
                       remove_dims=("bottom_top",),
                       description="maximum radar reflectivity",
                       units="dBz",
                       MemoryOrder="XY")
def get_max_dbz(wrfnc, timeidx=0, method="cat", 
                squeeze=True, cache=None, meta=True,
                do_varint=False, do_liqskin=False):
    return n.amax(get_dbz(wrfnc, timeidx, method, 
                          squeeze, cache, meta,
                          do_varint, do_liqskin), 
                  axis=-3)

