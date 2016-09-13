from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .util import extract_vars, get_id
from .latlonutils import (_lat_varname, _lon_varname, _ll_to_xy, _xy_to_ll)
from .metadecorators import set_latlon_metadata


def get_lat(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True, _key=None,
            stagger=None):
    
    varname = _lat_varname(wrfnc, stagger)
    lat_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           meta, _key)
    
    return lat_var[varname]

        
def get_lon(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True, _key=None,
            stagger=None):
    
    varname = _lon_varname(wrfnc, stagger)
    lon_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           meta, _key)
    
    return lon_var[varname]

# TODO:  Do we need the user to know about method, squeeze, cache for this?

# Can either use wrfnc as a single file or sequence, or provide 
# projection parameters (which don't allow for moving domains)
@set_latlon_metadata(xy=True) 
def ll_to_xy(wrfnc, latitude, longitude, timeidx=0, stagger=None, method="cat", 
             squeeze=True, cache=None, meta=True, as_int=True):
    _key = get_id(wrfnc)
    return _ll_to_xy(latitude, longitude, wrfnc, timeidx, stagger, method, 
                     squeeze, cache, _key, as_int, **{})


@set_latlon_metadata(xy=True) 
def ll_to_xy_proj(latitude, longitude, meta=True, squeeze=True, as_int=True,
                  map_proj=None, truelat1=None, truelat2=None, stand_lon=None, 
                  ref_lat=None, ref_lon=None, pole_lat=None, pole_lon=None, 
                  known_x=None, known_y=None, dx=None, dy=None, 
                  latinc=None, loninc=None):

    loc = locals()
    projparams = {name : loc[name] for name in ("map_proj", "truelat1", 
                                            "truelat2", "stand_lon", "ref_lat",
                                            "ref_lon", "pole_lat", "pole_lon",
                                            "known_x", "known_y", "dx", "dy",
                                            "latinc", "loninc")}

    return _ll_to_xy(latitude, longitude, None, 0, squeeze, "cat", True, None,
                     None, as_int, **projparams)
    
    
@set_latlon_metadata(xy=False) 
def xy_to_ll(wrfnc, x, y, timeidx=0, stagger=None, method="cat", squeeze=True, 
             cache=None, meta=True):
    _key = get_id(wrfnc)
    return _xy_to_ll(x, y, wrfnc, timeidx, stagger, method, squeeze, cache, 
                     _key, **{})
 
    
@set_latlon_metadata(xy=False) 
def xy_to_ll_proj(x, y, meta=True, squeeze=True, map_proj=None, truelat1=None, 
                  truelat2=None, stand_lon=None, ref_lat=None, ref_lon=None, 
                  pole_lat=None, pole_lon=None, known_x=None, known_y=None, 
                  dx=None, dy=None, latinc=None, loninc=None):
    loc = locals()
    projparams = {name : loc[name] for name in ("map_proj", "truelat1", 
                                            "truelat2", "stand_lon", "ref_lat",
                                            "ref_lon", "pole_lat", "pole_lon",
                                            "known_x", "known_y", "dx", "dy",
                                            "latinc", "loninc")}
    return _xy_to_ll(x, y, None, 0, None, "cat", squeeze, None, None,
                     **projparams)
    
    
    