
from .decorators import convert_units, copy_and_set_metadata
from .util import extract_vars, either

__all__ = ["get_terrain"]

# Need to handle either
@copy_and_set_metadata(copy_varname=either("HGT", "HGT_M"), name="terrain", 
                       description="terrain height")
@convert_units("height", "m")
def get_terrain(wrfnc, timeidx=0, units="m", method="cat", squeeze=True, 
              cache=None):
    varname = either("HGT", "HGT_M")(wrfnc)
    return extract_vars(wrfnc, timeidx, varname, 
                        method, squeeze, cache)[varname]

    
        