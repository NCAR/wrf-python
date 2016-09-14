from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

#from .extension import computeslp, computetk
from .extension import _slp, _tk
from .constants import Constants
from .destag import destagger
from .decorators import convert_units
from .metadecorators import copy_and_set_metadata
from .util import extract_vars


@copy_and_set_metadata(copy_varname="T", name="slp",
                       remove_dims=("bottom_top",), 
                       description="sea level pressure",
                       MemoryOrder="XY")
@convert_units("pressure", "hpa")
def get_slp(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True, _key=None,
            units="hPa"):
    varnames=("T", "P", "PB", "QVAPOR", "PH", "PHB")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)

    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qvapor = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    qvapor[qvapor < 0] = 0.
    
    full_ph = (ph + phb) / Constants.G
    
    destag_ph = destagger(full_ph, -3)
    
    tk = _tk(full_p, full_t)
    slp = _slp(destag_ph, tk, full_p, qvapor)
    
    return slp

