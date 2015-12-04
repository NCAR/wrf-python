
import numpy as n

from wrf.var.extension import computectt, computetk
from wrf.var.constants import Constants, ConversionFactors
from wrf.var.destagger import destagger
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars

__all__ = ["get_ctt"]

@convert_units("temp", "c")
def get_ctt(wrfnc, units="c", timeidx=0):
    """Return the cloud top temperature.
    
    """
    vars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "PH" ,"PHB",
                                              "HGT", "QVAPOR"))
    t = vars["T"]
    p = vars["P"]
    pb = vars["PB"]
    ph = vars["PH"]
    phb = vars["PHB"]
    ter = vars["HGT"]
    qv = vars["QVAPOR"] * 1000.0 # g/kg
    
    haveqci = 1
    try:
        icevars = extract_vars(wrfnc, timeidx, vars="QICE")
    except KeyError:
        qice = n.zeros(qv.shape, qv.dtype)
        haveqci = 0
    else:
        qice = icevars["QICE"] * 1000.0 #g/kg
    
    try:
        cldvars = extract_vars(wrfnc, timeidx, vars="QCLOUD")
    except KeyError:
        raise RuntimeError("'QCLOUD' not found in NetCDF file")
    else:
        qcld = cldvars["QCLOUD"] * 1000.0 #g/kg
    
    full_p = p + pb
    p_hpa = full_p * ConversionFactors.PA_TO_HPA
    full_t = t + Constants.T_BASE
    tk = computetk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, 0)
    ght = geopt_unstag / Constants.G
    
    ctt = computectt(p_hpa,tk,qice,qcld,qv,ght,ter,haveqci)
    
    return ctt
