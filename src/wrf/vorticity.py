from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .extension import _avo, _pvo
from .util import extract_vars, extract_global_attrs
from .metadecorators import copy_and_set_metadata


@copy_and_set_metadata(copy_varname="T", name="avo", 
                       description="absolute vorticity",
                       units="10-5 s-1")
def get_avo(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
            meta=True):
    ncvars = extract_vars(wrfnc, timeidx, ("U", "V", "MAPFAC_U",
                                           "MAPFAC_V", "MAPFAC_M",
                                           "F"),
                          method, squeeze, cache, meta=False)
    
    attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    u = ncvars["U"]
    v = ncvars["V"]
    msfu = ncvars["MAPFAC_U"]
    msfv = ncvars["MAPFAC_V"]
    msfm = ncvars["MAPFAC_M"]
    cor = ncvars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    return _avo(u, v, msfu, msfv, msfm, cor, dx, dy)


@copy_and_set_metadata(copy_varname="T", name="pvo", 
                       description="potential vorticity",
                       units="PVU")
def get_pvo(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
            meta=True):
    ncvars = extract_vars(wrfnc, timeidx, ("U", "V", "T", "P",
                                           "PB", "MAPFAC_U",
                                           "MAPFAC_V", "MAPFAC_M",
                                           "F"),
                          method, squeeze, cache, meta=False)
    attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    
    u = ncvars["U"]
    v = ncvars["V"]
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    msfu = ncvars["MAPFAC_U"]
    msfv = ncvars["MAPFAC_V"]
    msfm = ncvars["MAPFAC_M"]
    cor = ncvars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    full_t = t + 300
    full_p = p + pb
    
    return _pvo(u, v, full_t, full_p, msfu, msfv, msfm, cor, dx, dy)
    