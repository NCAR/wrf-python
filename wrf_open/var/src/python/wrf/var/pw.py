
from wrf.var.extension import computepw,computetv,computetk
from wrf.var.constants import Constants

__all__ = ["get_pw"]

def get_pw(wrfnc, timeidx=0):
    
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    ph = wrfnc.variables["PH"][timeidx,:,:,:]
    phb = wrfnc.variables["PHB"][timeidx,:,:,:]
    qv = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    
    # Change this to use real virtual temperature!
    full_p   =  p + pb
    ht = (ph + phb)/Constants.G
    full_t  =  t + Constants.T_BASE
    
    tk = computetk(full_p, full_t)
    tv = computetv(tk,qv)
    
    return computepw(full_p,tv,qv,ht)
    
    
    
    