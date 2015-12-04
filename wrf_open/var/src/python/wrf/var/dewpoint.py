from wrf.var.extension import computetd
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars

__all__ = ["get_dp", "get_dp_2m"]

@convert_units("temp", "c")
def get_dp(wrfnc, units="c", timeidx=0):
    
    vars = extract_vars(wrfnc, timeidx, vars=("P", "PB", "QVAPOR"))
    
    p = vars["P"]
    pb = vars["PB"]
    qvapor = vars["QVAPOR"]
    
    # Algorithm requires hPa
    full_p = .01*(p + pb)
    qvapor[qvapor < 0] = 0
    
    td = computetd(full_p, qvapor)
    return td
    
@convert_units("temp", "c")
def get_dp_2m(wrfnc, units="c", timeidx=0):
    vars = extract_vars(wrfnc, timeidx, vars=("PSFC", "Q2"))

    # Algorithm requires hPa
    psfc = .01*(vars["PSFC"])
    q2 = vars["Q2"]
    q2[q2 < 0] = 0
    
    td = computetd(psfc, q2)
    
    return td

