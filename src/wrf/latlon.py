from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .util import extract_vars
from .latlonutils import (_lat_varname, _lon_varname, _ll_to_xy, _xy_to_ll)
from .metadecorators import set_latlon_metadata


def get_lat(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True,
            stagger=None):
    
    varname = _lat_varname(wrfnc, stagger)
    lat_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           meta)
    
    return lat_var[varname]

        
def get_lon(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True,
            stagger=None):
    
    varname = _lon_varname(wrfnc, stagger)
    lon_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           meta)
    
    return lon_var[varname]


# Can either use wrfnc as a single file or sequence, or provide 
# projection parameters (which don't allow for moving domains)
@set_latlon_metadata(xy=True) 
def ll_to_xy(wrfnc, latitude, longitude, timeidx=0, stagger=None, method="cat", 
             squeeze=True, cache=None, meta=True):
    return _ll_to_xy(latitude, longitude, wrfnc, timeidx, stagger, method, 
                     squeeze, cache, **{})


@set_latlon_metadata(xy=True) 
def ll_to_xy_proj(latitude, longitude, meta=True, squeeze=True, **projparams):
    return _ll_to_xy(latitude, longitude, None, 0, squeeze, "cat", True, None, 
                     **projparams)
    
    
@set_latlon_metadata(xy=False) 
def xy_to_ll(wrfnc, x, y, timeidx=0, stagger=None, method="cat", squeeze=True, 
             cache=None, meta=True):
    return _xy_to_ll(x, y, wrfnc, timeidx, stagger, method, squeeze, cache, 
                     **{})
 
    
@set_latlon_metadata(xy=False) 
def xy_to_ll_proj(x, y, meta=True, squeeze=True, **projparams):
    return _xy_to_ll(x, y, None, 0, None, "cat", squeeze, None, **projparams)
    
    
    