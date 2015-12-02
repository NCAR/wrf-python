
from wrf.var.constants import Constants
from wrf.var.decorators import convert_units

__all__ = ["get_pressure"]

@convert_units("pressure", "pa")
def get_pressure(wrfnc, units="hpa", timeidx=0):

    if "P" in wrfnc.variables:
        p = wrfnc.variables["P"][timeidx,:,:,:]
        pb = wrfnc.variables["PB"][timeidx,:,:,:]
        pres = p + pb
    elif "PRES" in wrfnc.variables:
        pres = wrfnc.variables["PRES"][timeidx,:,:,:]
    
    return pres

    
    