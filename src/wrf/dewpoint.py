from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

#from .extension import computetd
from .extension import _td
from .decorators import convert_units
from .metadecorators import copy_and_set_metadata
from .util import extract_vars


@copy_and_set_metadata(copy_varname="QVAPOR", name="td", 
                       description="dew point temperature")
@convert_units("temp", "c")
def get_dp(wrfnc, timeidx=0, method="cat", squeeze=True, 
           cache=None, meta=True, units="c"):
    
    varnames=("P", "PB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)
    
    p = ncvars["P"]
    pb = ncvars["PB"]
    qvapor = ncvars["QVAPOR"]
    
    # Algorithm requires hPa
    full_p = .01*(p + pb)
    qvapor[qvapor < 0] = 0
    
    td = _td(full_p, qvapor)
    return td

@copy_and_set_metadata(copy_varname="Q2", name="td2", 
                       description="2m dew point temperature")
@convert_units("temp", "c")
def get_dp_2m(wrfnc, timeidx=0, method="cat", squeeze=True, 
              cache=None, meta=True, 
              units="c"):
    varnames=("PSFC", "Q2")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)

    # Algorithm requires hPa
    psfc = .01*(ncvars["PSFC"])
    q2 = ncvars["Q2"]
    q2[q2 < 0] = 0
    
    td = _td(psfc, q2)
    
    return td

