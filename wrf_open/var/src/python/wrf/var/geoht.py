from wrf.var.constants import Constants
from wrf.var.destagger import destagger
from wrf.var.decorators import convert_units

__all__ = ["get_geopt", "get_height"]

def _get_geoht(wrfnc, height=True, msl=True, timeidx=0):
    """Return the geopotential in units of m2 s-2 if height is False,
    otherwise return the geopotential height in meters.  If height is True,
    then if msl is True the result will be in MSL, otherwise AGL (the terrain
    height is subtracted).
    
    """

    if "PH" in wrfnc.variables:
        ph = wrfnc.variables["PH"][timeidx,:,:,:]
        phb = wrfnc.variables["PHB"][timeidx,:,:,:]
        hgt = wrfnc.variables["HGT"][timeidx,:,:]
        geopt = ph + phb
        geopt_unstag = destagger(geopt, 0)
    elif "GHT" in wrfnc.variables: # met_em files
        geopt_unstag = wrfnc.variables["GHT"][timeidx,:,:,:] * Constants.G
        hgt = destagger(wrfnc.variables["HGT_U"][timidx,:,:], 1)
    
    if height:
        if msl:
            return geopt_unstag / Constants.G
        else:
            return (geopt_unstag / Constants.G) - hgt
    else:
        return geopt_unstag

def get_geopt(wrfnc, timeidx=0):
    return _get_geoht(wrfnc, False, timeidx=timeidx)

@convert_units("height", "m")
def get_height(wrfnc, msl=True, units="m", timeidx=0):
    z = _get_geoht(wrfnc, True, msl, timeidx)
    return z

