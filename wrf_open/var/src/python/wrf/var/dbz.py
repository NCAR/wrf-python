import numpy as n

from wrf.var.extension import computedbz,computetk
from wrf.var.constants import Constants

__all__ = ["get_dbz", "get_max_dbz"]

def get_dbz(wrfnc, do_varint=False, do_liqskin=False, timeidx=0):
    """ Return the dbz
    
    do_varint - do variable intercept (if False, constants are used.  Otherwise, 
    intercepts are calculated using a technique from Thompson, Rasmussen, 
    and Manning (2004, Monthly Weather Review, Vol. 132, No. 2, pp. 519-542.)
    
    do_liqskin - do liquid skin for snow (frozen particles above freezing scatter
    as liquid)
    
    """
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    
    qv = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    qr = wrfnc.variables["QRAIN"][timeidx,:,:,:]
    
    if "QSNOW" in wrfnc.variables:
        qs = wrfnc.variables["QSNOW"][timeidx,:,:,:]
    else:
        qs = n.zeros((qv.shape[0], qv.shape[1], qv.shape[2]), "float")
        
    if "QGRAUP" in wrfnc.variables:
        qg = wrfnc.variables["QGRAUP"][timeidx,:,:,:]
    else:
        qg = n.zeros((qv.shape[0], qv.shape[1], qv.shape[2]), "float")
    
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

def get_max_dbz(wrfnc, do_varint=False, do_liqskin=False, timeidx=0):
    return n.amax(get_dbz(wrfnc, do_varint, do_liqskin, timeidx), 
                  axis=0)

