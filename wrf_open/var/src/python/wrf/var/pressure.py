
from wrf.var.constants import Constants
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars

__all__ = ["get_pressure"]

@convert_units("pressure", "pa")
def get_pressure(wrfnc, units="hpa", timeidx=0):

    try:
        p_vars = extract_vars(wrfnc, timeidx, vars=("P", "PB"))
    except KeyError:
        try:
            pres_vars = extract_vars(wrfnc, timeidx, vars="PRES")
        except:
            raise RuntimeError("pressure variable not found in NetCDF file")
        else:
            pres = pres_vars["PRES"]
    else:
        p = p_vars["P"]
        pb = p_vars["PB"]
        pres = p + pb
    
    return pres

    
    