from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants
from .destag import destagger
from .decorators import convert_units
from .metadecorators import set_height_metadata
from .util import extract_vars, either

def _get_geoht(wrfnc, timeidx, method="cat", squeeze=True, 
               cache=None, meta=True, _key=None,
               height=True, msl=True):
    """Return the geopotential in units of m2 s-2 if height is False,
    otherwise return the geopotential height in meters.  If height is True,
    then if msl is True the result will be in MSL, otherwise AGL (the terrain
    height is subtracted).
    
    """
    
    varname = either("PH", "GHT")(wrfnc)
    if varname == "PH":
        ph_vars = extract_vars(wrfnc, timeidx, ("PH", "PHB", "HGT"),
                               method, squeeze, cache, meta=False,
                               _key=_key)
        ph = ph_vars["PH"]
        phb = ph_vars["PHB"]
        hgt = ph_vars["HGT"]
        geopt = ph + phb
        geopt_unstag = destagger(geopt, -3)
    else:
        ght_vars = extract_vars(wrfnc, timeidx, ("GHT", "HGT_M"),
                                method, squeeze, cache, meta=False,
                                _key=_key)
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
def get_geopt(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
              meta=True, _key=None):
    return _get_geoht(wrfnc, timeidx, method, squeeze, cache, meta, _key,
                      False, True)


@set_height_metadata(geopt=False)
@convert_units("height", "m")
def get_height(wrfnc, timeidx=0, method="cat", squeeze=True, 
               cache=None, meta=True, _key=None,
               msl=True, units="m"):
    
    return _get_geoht(wrfnc, timeidx, method, squeeze, cache, meta, _key,
                      True, msl)

