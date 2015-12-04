
from wrf.var.extension import computepw,computetv,computetk
from wrf.var.constants import Constants
from wrf.var.util import extract_vars

__all__ = ["get_pw"]

def get_pw(wrfnc, timeidx=0):
    vars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "PH", "PHB", 
                                              "QVAPOR"))
    
    t = vars["T"]
    p = vars["P"]
    pb = vars["PB"]
    ph = vars["PH"]
    phb = vars["PHB"]
    qv = vars["QVAPOR"]
    
    # Change this to use real virtual temperature!
    full_p   =  p + pb
    ht = (ph + phb)/Constants.G
    full_t  =  t + Constants.T_BASE
    
    tk = computetk(full_p, full_t)
    tv = computetv(tk,qv)
    
    return computepw(full_p,tv,qv,ht)
    
    
    
    