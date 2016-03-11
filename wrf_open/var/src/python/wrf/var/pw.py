
from wrf.var.extension import computepw,computetv,computetk
from wrf.var.constants import Constants
from wrf.var.util import extract_vars

__all__ = ["get_pw"]

def get_pw(wrfnc, timeidx=0):
    ncvars = extract_vars(wrfnc, timeidx, varnames=("T", "P", "PB", "PH", 
                                                    "PHB", "QVAPOR"))
    
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    qv = ncvars["QVAPOR"]
    
    # Change this to use real virtual temperature!
    full_p   =  p + pb
    ht = (ph + phb)/Constants.G
    full_t  =  t + Constants.T_BASE
    
    tk = computetk(full_p, full_t)
    tv = computetv(tk,qv)
    
    return computepw(full_p,tv,qv,ht)
    
    
    
    