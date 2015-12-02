
from wrf.var.decorators import convert_units

__all__ = ["get_terrain"]

@convert_units("height", "m")
def get_terrain(wrfnc, units="m", timeidx=0):

    if "HGT" in wrfnc.variables:
        hgt = wrfnc.variables["HGT"][timeidx,:,:]
    elif "HGT_M":
        hgt = wrfnc.variables["HGT_M"][timeidx,:,:] 
    
    return hgt
    
        