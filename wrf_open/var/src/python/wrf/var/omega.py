
from wrf.var.constants import Constants
from wrf.var.destagger import destagger
from wrf.var.extension import computeomega,computetk

__all__ = ["get_omega"]

def get_omega(wrfnc, timeidx=0):
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    w = wrfnc.variables["W"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    qv = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    
    wa = destagger(w, 0)
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    omega = computeomega(qv,tk,wa,full_p)
    
    return omega
    