from wrf.var.extension import computeslp, computetk
from wrf.var.constants import Constants
from wrf.var.destag import destagger
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars

__all__ = ["get_slp"]

@convert_units("pressure", "hpa")
def get_slp(wrfnc, timeidx=0, units="hpa"):
    ncvars = extract_vars(wrfnc, timeidx, varnames=("T", "P", "PB", "QVAPOR",
                                              "PH", "PHB"))

    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qvapor = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    qvapor[qvapor < 0] = 0.
    full_ph = (ph + phb) / Constants.G
    
    destag_ph = destagger(full_ph, -3)
    
    tk = computetk(full_p, full_t)
    slp = computeslp(destag_ph, tk, full_p, qvapor)
    
    return slp

