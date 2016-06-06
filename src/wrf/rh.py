from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants   
#from .extension import computerh, computetk
from .extension import _rh, _tk
from .util import extract_vars
from .metadecorators import copy_and_set_metadata

__all__ = ["get_rh", "get_rh_2m"]

@copy_and_set_metadata(copy_varname="T", name="rh", 
                       description="relative humidity",
                       delete_attrs=("units",))
def get_rh(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
           meta=True):
    varnames=("T", "P", "PB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qvapor = ncvars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    qvapor[qvapor < 0] = 0
    tk = _tk(full_p, full_t)
    rh = _rh(qvapor, full_p, tk)
    
    return rh

@copy_and_set_metadata(copy_varname="T2", name="rh2", 
                       description="2m relative humidity",
                       delete_attrs=("units",))
def get_rh_2m(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None,
              meta=True):
    varnames=("T2", "PSFC", "Q2")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)
    t2 = ncvars["T2"]
    psfc = ncvars["PSFC"]
    q2 = ncvars["Q2"]
    
    q2[q2 < 0] = 0
    rh = _rh(q2, psfc, t2)
    
    return rh

