
import numpy as n

from wrf.var.extension import computectt, computetk
from wrf.var.constants import Constants, ConversionFactors
from wrf.var.destagger import destagger
from wrf.var.decorators import convert_units

__all__ = ["get_ctt"]

@convert_units("temp", "c")
def get_ctt(wrfnc, units="c", timeidx=0):
    """Return the cloud top temperature.
    
    """
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    ph = wrfnc.variables["PH"][timeidx,:,:,:]
    phb = wrfnc.variables["PHB"][timeidx,:,:,:]
    ter = wrfnc.variables["HGT"][timeidx,:,:]
    qv = wrfnc.variables["QVAPOR"][timeidx,:,:,:] * 1000.0 # g/kg
    
    haveqci = 1
    if "QICE" in wrfnc.variables:
        qice = wrfnc.variables["QICE"][timeidx,:,:,:] * 1000.0 #g/kg
    else:
        qice = n.zeros(qv.shape, qv.dtype)
        haveqci = 0
        
    if "QCLOUD" in wrfnc.variables:
        qcld = wrfnc.variables["QCLOUD"][timeidx,:,:,:] * 1000.0 #g/kg
    else:
        raise RuntimeError("'QCLOUD' not found in NetCDF file")
    
    full_p = p + pb
    p_hpa = full_p * ConversionFactors.PA_TO_HPA
    full_t = t + Constants.T_BASE
    tk = computetk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, 0)
    ght = geopt_unstag / Constants.G
    
    ctt = computectt(p_hpa,tk,qice,qcld,qv,ght,ter,haveqci)
    
    return ctt