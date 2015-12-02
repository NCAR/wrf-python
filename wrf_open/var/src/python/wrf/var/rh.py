
from wrf.var.constants import Constants   
from wrf.var.extension import computerh, computetk

__all__ = ["get_rh", "get_rh_2m"]

def get_rh(wrfnc, timeidx=0):
    t = wrfnc.variables["T"][timeidx,:,:,:]
    #t00 = wrfnc.variables["T00"][timeidx]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    qvapor = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    qvapor[qvapor < 0] = 0
    tk = computetk(full_p, full_t)
    rh = computerh(qvapor, full_p, tk)
    
    return rh

def get_rh_2m(wrfnc, timeidx=0):
    t2 = wrfnc.variables["T2"][timeidx,:,:]
    psfc = wrfnc.variables["PSFC"][timeidx,:,:]
    q2 = wrfnc.variables["Q2"][timeidx,:,:]
    
    q2[q2 < 0] = 0
    rh = computerh(q2, psfc, t2)
    
    return rh

