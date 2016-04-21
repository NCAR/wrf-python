from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .extension import computepw,computetv,computetk
from .constants import Constants
from .util import extract_vars
from .metadecorators import copy_and_set_metadata

__all__ = ["get_pw"]

@copy_and_set_metadata(copy_varname="T", name="pw", 
                       remove_dims=("bottom_top",),
                       description="precipitable water",
                       MemoryOrder="XY",
                       units="kg m-2")
def get_pw(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
           meta=True):
    varnames=("T", "P", "PB", "PH", "PHB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)
    
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    qv = ncvars["QVAPOR"]
    
    # Change this to use real virtual temperature!
    full_p   =  p + pb
    ht = (ph + phb)/Constants.G
    full_t  =  t + Constants.T_BASE
    
    tk = computetk(full_p, full_t)
    tv = computetv(tk, qv)
    
    return computepw(full_p, tv, qv, ht)
    
    
    
    