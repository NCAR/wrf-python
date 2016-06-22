from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .decorators import convert_units
from .metadecorators import copy_and_set_metadata
from .util import extract_vars, either


@copy_and_set_metadata(copy_varname=either("P", "PRES"), name="pressure", 
                       description="pressure")
@convert_units("pressure", "pa")
def get_pressure(wrfnc, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True,
                 units="pa"):
    
    varname = either("P", "PRES")(wrfnc)
    if varname == "P":
        p_vars = extract_vars(wrfnc, timeidx, ("P", "PB"), 
                              method, squeeze, cache, meta=False)
        p = p_vars["P"]
        pb = p_vars["PB"]
        pres = p + pb
    else:
        pres = extract_vars(wrfnc, timeidx, "PRES", 
                            method, squeeze, cache, meta=False)["PRES"]
    
    return pres

def get_pressure_hpa(wrfnc, timeidx=0, method="cat", squeeze=True, 
                     cache=None, meta=True,
                     units="hpa"):
    return get_pressure(wrfnc, timeidx, method, squeeze, cache, meta, units)


    
    