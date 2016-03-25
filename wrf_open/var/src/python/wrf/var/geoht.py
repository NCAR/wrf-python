from .constants import Constants
from .destag import destagger
from .decorators import convert_units, set_height_metadata
from .util import extract_vars, either

__all__ = ["get_geopt", "get_height"]

def _get_geoht(wrfnc, timeidx, height=True, msl=True,
               method="cat", squeeze=True, cache=None):
    """Return the geopotential in units of m2 s-2 if height is False,
    otherwise return the geopotential height in meters.  If height is True,
    then if msl is True the result will be in MSL, otherwise AGL (the terrain
    height is subtracted).
    
    """
    
    varname = either("PH", "GHT")(wrfnc)
    if varname == "PH":
        ph_vars = extract_vars(wrfnc, timeidx, varnames=("PH", "PHB", "HGT"))
        ph = ph_vars["PH"]
        phb = ph_vars["PHB"]
        hgt = ph_vars["HGT"]
        geopt = ph + phb
        geopt_unstag = destagger(geopt, -3)
    else:
        ght_vars = extract_vars(wrfnc, timeidx, varnames=("GHT", "HGT_M"))
        geopt_unstag = ght_vars["GHT"] * Constants.G
        hgt = ght_vars["HGT_M"]
    
    if height:
        if msl:
            return geopt_unstag / Constants.G
        else:
            # Due to broadcasting with multifile/multitime, the 2D terrain
            # array needs to be reshaped to a 3D array so the right dims
            # line up
            new_dims = list(hgt.shape)
            new_dims.insert(-2,1)
            hgt = hgt.reshape(new_dims)
            
            return (geopt_unstag / Constants.G) - hgt
    else:
        return geopt_unstag

@set_height_metadata(geopt=True)
def get_geopt(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None):
    return _get_geoht(wrfnc, timeidx, False, True, method, squeeze, cache)

@set_height_metadata(geopt=False)
@convert_units("height", "m")
def get_height(wrfnc, timeidx=0, msl=True, units="m",
               method="cat", squeeze=True, cache=None):
    
    return _get_geoht(wrfnc, timeidx, True, msl, method, squeeze, cache)

