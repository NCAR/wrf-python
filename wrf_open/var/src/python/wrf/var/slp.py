from wrf.var.extension import computeslp, computetk
from wrf.var.constants import Constants
from wrf.var.destagger import destagger
from wrf.var.decorators import convert_units

__all__ = ["get_slp"]

@convert_units("pressure", "hpa")
def get_slp(wrfnc, units="hpa", timeidx=0):

    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    qvapor = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    ph = wrfnc.variables["PH"][timeidx,:,:,:]
    phb = wrfnc.variables["PHB"][timeidx,:,:,:]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    qvapor[qvapor < 0] = 0.
    full_ph = (ph + phb) / Constants.G
    
    destag_ph = destagger(full_ph, 0)
    
    tk = computetk(full_p, full_t)
    slp = computeslp(destag_ph, tk, full_p, qvapor)
    
    return slp

