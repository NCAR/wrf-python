
from wrf.var.constants import Constants
from wrf.var.destagger import destagger
from wrf.var.extension import computeomega,computetk
from wrf.var.util import extract_vars

__all__ = ["get_omega"]

def get_omega(wrfnc, timeidx=0):
    ncvars = extract_vars(wrfnc, timeidx, vars=("T", "P", "W", "PB", "QVAPOR"))
    t = ncvars["T"]
    p = ncvars["P"]
    w = ncvars["W"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    wa = destagger(w, 0)
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    omega = computeomega(qv,tk,wa,full_p)
    
    return omega
    