from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .extension import computetd
from .decorators import convert_units
from .metadecorators import copy_and_set_metadata
from .util import extract_vars

__all__ = ["get_dp", "get_dp_2m"]

@copy_and_set_metadata(copy_varname="QVAPOR", name="td", 
                       description="dew point temperature")
@convert_units("temp", "c")
def get_dp(wrfnc, timeidx=0, units="c",
           method="cat", squeeze=True, cache=None):
    
    varnames=("P", "PB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          nometa=True)
    
    p = ncvars["P"]
    pb = ncvars["PB"]
    qvapor = ncvars["QVAPOR"]
    
    # Algorithm requires hPa
    full_p = .01*(p + pb)
    qvapor[qvapor < 0] = 0
    
    td = computetd(full_p, qvapor)
    return td

@copy_and_set_metadata(copy_varname="Q2", name="td2", 
                       description="2m dew point temperature")
@convert_units("temp", "c")
def get_dp_2m(wrfnc, timeidx=0, units="c",
              method="cat", squeeze=True, cache=None):
    varnames=("PSFC", "Q2")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          nometa=True)

    # Algorithm requires hPa
    psfc = .01*(ncvars["PSFC"])
    q2 = ncvars["Q2"]
    q2[q2 < 0] = 0
    
    td = computetd(psfc, q2)
    
    return td

