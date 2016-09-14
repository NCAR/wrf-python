from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np
import numpy.ma as ma

from .extension import _tk, _cape
from .destag import destagger
from .constants import Constants, ConversionFactors
from .util import extract_vars
from .metadecorators import set_cape_metadata

@set_cape_metadata(is2d=True)
def get_2dcape(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
               meta=True, _key=None, missing=Constants.DEFAULT_FILL):
    """Return the 2d fields of cape, cin, lcl, and lfc"""

    varnames = ("T", "P", "PB", "QVAPOR", "PH","PHB", "HGT", "PSFC")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    ter = ncvars["HGT"]
    psfc = ncvars["PSFC"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    z = geopt_unstag/Constants.G
    
    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc 
    
    i3dflag = 0
    ter_follow = 1
    
    cape_cin = _cape(p_hpa, tk, qv, z, ter, psfc_hpa, missing, i3dflag, 
                     ter_follow)
    
    left_dims = cape_cin.shape[1:-3]
    right_dims = cape_cin.shape[-2:]
    
    resdim = (4,) + left_dims + right_dims
    
    # Make a new output array for the result
    result = np.zeros(resdim, cape_cin.dtype)
    
    # Cape 2D output is not flipped in the vertical, so index from the 
    # end
    result[0,...,:,:] = cape_cin[0,...,-1,:,:]
    result[1,...,:,:] = cape_cin[1,...,-1,:,:]
    result[2,...,:,:] = cape_cin[1,...,-2,:,:]
    result[3,...,:,:] = cape_cin[1,...,-3,:,:]
    
    return ma.masked_values(result, missing)


@set_cape_metadata(is2d=False)
def get_3dcape(wrfnc, timeidx=0, method="cat", 
               squeeze=True, cache=None, meta=True,
               _key=None, missing=Constants.DEFAULT_FILL):
    """Return the 3d fields of cape and cin"""
    varnames = ("T", "P", "PB", "QVAPOR", "PH", "PHB", "HGT", "PSFC")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    ter = ncvars["HGT"]
    psfc = ncvars["PSFC"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    z = geopt_unstag/Constants.G
    
    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc 
    
    i3dflag = 1
    ter_follow = 1
    
    cape_cin = _cape(p_hpa, tk, qv, z, ter, psfc_hpa, missing, i3dflag, 
                     ter_follow)
    
    
    return ma.masked_values(cape_cin, missing)
    
    
    