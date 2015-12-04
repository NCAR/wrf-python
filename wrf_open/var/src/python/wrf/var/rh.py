
from wrf.var.constants import Constants   
from wrf.var.extension import computerh, computetk
from wrf.var.util import extract_vars

__all__ = ["get_rh", "get_rh_2m"]

def get_rh(wrfnc, timeidx=0):
    vars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "QVAPOR"))
    t = vars["T"]
    p = vars["P"]
    pb = vars["PB"]
    qvapor = vars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    qvapor[qvapor < 0] = 0
    tk = computetk(full_p, full_t)
    rh = computerh(qvapor, full_p, tk)
    
    return rh

def get_rh_2m(wrfnc, timeidx=0):
    vars = extract_vars(wrfnc, timeidx, vars=("T2", "PSFC", "Q2"))
    t2 = vars["T2"]
    psfc = vars["PSFC"]
    q2 = vars["Q2"]
    
    q2[q2 < 0] = 0
    rh = computerh(q2, psfc, t2)
    
    return rh

