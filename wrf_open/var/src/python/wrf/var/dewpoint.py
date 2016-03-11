from wrf.var.extension import computetd
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars

__all__ = ["get_dp", "get_dp_2m"]

@convert_units("temp", "c")
def get_dp(wrfnc, timeidx=0, units="c"):
    
    ncvars = extract_vars(wrfnc, timeidx, varnames=("P", "PB", "QVAPOR"))
    
    p = ncvars["P"]
    pb = ncvars["PB"]
    qvapor = ncvars["QVAPOR"]
    
    # Algorithm requires hPa
    full_p = .01*(p + pb)
    qvapor[qvapor < 0] = 0
    
    td = computetd(full_p, qvapor)
    return td
    
@convert_units("temp", "c")
def get_dp_2m(wrfnc, timeidx=0, units="c"):
    ncvars = extract_vars(wrfnc, timeidx, varnames=("PSFC", "Q2"))

    # Algorithm requires hPa
    psfc = .01*(ncvars["PSFC"])
    q2 = ncvars["Q2"]
    q2[q2 < 0] = 0
    
    td = computetd(psfc, q2)
    
    return td

