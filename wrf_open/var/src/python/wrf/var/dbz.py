import numpy as n

from wrf.var.extension import computedbz,computetk
from wrf.var.constants import Constants
from wrf.var.util import extract_vars

__all__ = ["get_dbz", "get_max_dbz"]

def get_dbz(wrfnc, timeidx=0, do_varint=False, do_liqskin=False):
    """ Return the dbz
    
    do_varint - do variable intercept (if False, constants are used.  Otherwise, 
    intercepts are calculated using a technique from Thompson, Rasmussen, 
    and Manning (2004, Monthly Weather Review, Vol. 132, No. 2, pp. 519-542.)
    
    do_liqskin - do liquid skin for snow (frozen particles above freezing scatter
    as liquid)
    
    """
    ncvars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "QVAPOR", 
                                              "QRAIN"))
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    qr = ncvars["QRAIN"]
    
    try:
        snowvars = extract_vars(wrfnc, timeidx, vars="QSNOW")
    except KeyError:
        qs = n.zeros(qv.shape, "float")
    else:
        qs = snowvars["QSNOW"]
    
    try:
        graupvars = extract_vars(wrfnc, timeidx, vars="QGRAUP")
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
    tk = computetk(full_p, full_t)
    
    ivarint = 0
    if do_varint:
        ivarint = 1
        
    iliqskin = 0
    if do_liqskin:
        iliqskin = 1
    
    return computedbz(full_p,tk,qv,qr,qs,qg,sn0,ivarint,iliqskin)

def get_max_dbz(wrfnc, timeidx=0, do_varint=False, do_liqskin=False):
    return n.amax(get_dbz(wrfnc, do_varint, do_liqskin, timeidx), 
                  axis=-3)

