
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars

__all__ = ["get_terrain"]

@convert_units("height", "m")
def get_terrain(wrfnc, timeidx=0, units="m"):
    
    try:
        hgt_vars = extract_vars(wrfnc, timeidx, vars="HGT")
    except KeyError:
        try:
            hgt_m_vars = extract_vars(wrfnc, timeidx, vars="HGT_M")
        except KeyError:
            raise RuntimeError("height variable not found in NetCDF file")
        else:
            hgt = hgt_m_vars["HGT_M"]
    else:
        hgt = hgt_vars["HGT"]
    
    return hgt
    
        