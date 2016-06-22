from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as n

#from .extension import computectt, computetk
from .extension import computectt, _tk
from .constants import Constants, ConversionFactors
from .destag import destagger
from .decorators import convert_units 
from .metadecorators import copy_and_set_metadata
from .util import extract_vars


@copy_and_set_metadata(copy_varname="T", name="ctt",
                       remove_dims=("bottom_top",),
                       description="cloud top temperature",
                       MemoryOrder="XY")
@convert_units("temp", "c")
def get_ctt(wrfnc, timeidx=0, method="cat", 
            squeeze=True, cache=None, meta=True,
            units="c"):
    """Return the cloud top temperature.
    
    """
    varnames = ("T", "P", "PB", "PH", "PHB", "HGT", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    ter = ncvars["HGT"]
    qv = ncvars["QVAPOR"] * 1000.0 # g/kg
    
    haveqci = 1
    try:
        icevars = extract_vars(wrfnc, timeidx, "QICE", 
                               method, squeeze, cache, meta=False)
    except KeyError:
        qice = n.zeros(qv.shape, qv.dtype)
        haveqci = 0
    else:
        qice = icevars["QICE"] * 1000.0 #g/kg
    
    try:
        cldvars = extract_vars(wrfnc, timeidx, "QCLOUD", 
                               method, squeeze, cache, meta=False)
    except KeyError:
        raise RuntimeError("'QCLOUD' not found in NetCDF file")
    else:
        qcld = cldvars["QCLOUD"] * 1000.0 #g/kg
    
    full_p = p + pb
    p_hpa = full_p * ConversionFactors.PA_TO_HPA
    full_t = t + Constants.T_BASE
    tk = _tk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    ght = geopt_unstag / Constants.G
    
    ctt = computectt(p_hpa, tk, qice, qcld, qv, ght, ter, haveqci)
    
    return ctt
