from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from collections import Iterable

import numpy as np

from .constants import Constants, ProjectionTypes
from .extension import _lltoxy, _xytoll
from .util import (extract_vars, extract_global_attrs, 
                   either, is_moving_domain, is_multi_time_req,
                   iter_left_indexes, is_mapping, is_multi_file)
from .py3compat import viewkeys, viewitems

def _lat_varname(wrfnc, stagger):
    if stagger is None or stagger.lower() == "m":
        varname = either("XLAT", "XLAT_M")(wrfnc)
    elif stagger.lower() == "u" or stagger.lower() == "v":
        varname = "XLAT_{}".format(stagger.upper())
    else:
        raise ValueError("invalid 'stagger' value")
    
    return varname
    
def _lon_varname(wrfnc, stagger):
    if stagger is None or stagger.lower() == "m":
        varname = either("XLONG", "XLONG_M")(wrfnc)
    elif stagger.lower() == "u" or stagger.lower() == "v":
        varname = "XLONG_{}".format(stagger.upper())
    else:
        raise ValueError("invalid 'stagger' value")
    
    return varname

def _get_proj_params(wrfnc, timeidx, stagger, method, squeeze, cache, _key):
    if timeidx < 0:
        raise ValueError("'timeidx' must be greater than 0")
    
    attrs = extract_global_attrs(wrfnc, attrs=("MAP_PROJ", "TRUELAT1",
                                               "TRUELAT2", "STAND_LON",
                                               "DX", "DY"))
    map_proj = attrs["MAP_PROJ"]
    truelat1 = attrs["TRUELAT1"]
    truelat2 = attrs["TRUELAT2"]
    stdlon = attrs["STAND_LON"]
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    if map_proj == ProjectionTypes.LAT_LON:
        pole_attrs = extract_global_attrs(wrfnc, attrs=("POLE_LAT", 
                                                        "POLE_LON"))
        pole_lat = pole_attrs["POLE_LAT"]
        pole_lon = pole_attrs["POLE_LON"]
        latinc = (dy*360.0)/2.0 / Constants.PI/Constants.WRF_EARTH_RADIUS
        loninc = (dx*360.0)/2.0 / Constants.PI/Constants.WRF_EARTH_RADIUS
    else:
        pole_lat = 90.0
        pole_lon = 0.0
        latinc = 0.0
        loninc = 0.0
    
    latvar = _lat_varname(wrfnc, stagger)
    lonvar = _lon_varname(wrfnc, stagger)
    
    lat_timeidx = timeidx
    
    is_moving = is_moving_domain(wrfnc, latvar=latvar, lonvar=lonvar)
    
    # Only need one file and one time if the domain is not moving
    if not is_moving:
        if is_multi_time_req(timeidx):
            lat_timeidx = 0
            
        if is_multi_file(wrfnc):
            if not is_mapping(wrfnc):
                wrfnc = next(iter(wrfnc))  # only need one file
            else:
                first_entry = next(iter(viewkeys(wrfnc)))
                wrfnc = wrfnc[first_entry]
                key = _key[first_entry]
                return _get_proj_params(wrfnc, timeidx, stagger, 
                                        method, squeeze, cache, key)
            
    xlat = extract_vars(wrfnc, lat_timeidx, (latvar,), method, squeeze, cache,
                           meta=False, _key=_key)[latvar]
    xlon = extract_vars(wrfnc, lat_timeidx, (lonvar,), method, squeeze, cache,
                           meta=False, _key=_key)[lonvar]
    
    ref_lat = np.ravel(xlat[..., 0, 0])
    ref_lon = np.ravel(xlon[..., 0, 0])
    
    # Note: fortran index
    known_x = 1.0
    known_y = 1.0
    
    return (map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon,
            pole_lat, pole_lon, known_x, known_y, dx, dy, latinc, loninc)


def _dict_keys_to_upper(d):
    return {key.upper(): val for key, val in viewitems(d)}        

# known_x and known_y are 0-based
def _kwarg_proj_params(projparams):
    projparams = _dict_keys_to_upper(projparams)
    
    map_proj = projparams.get("MAP_PROJ"),
    truelat1 = projparams.get("TRUELAT1"),
    truelat2 = projparams.get("TRUELAT2"),
    stdlon = projparams.get("STDLON"),
    ref_lat = projparams.get("REF_LAT"),
    ref_lon = projparams.get("REF_LON"),
    pole_lat = projparams.get("POLE_LAT", 90.0),
    pole_lon = projparams.get("POLE_LON", 0.0),
    known_x = projparams.get("KNOWN_X"), # Use 0-based
    known_y = projparams.get("KNOWN_Y"), # Use 0-based
    dx = projparams.get("DX"),
    dy = projparams.get("DY"),
    latinc = projparams.get("LATINC"),
    loninc = projparams.get("LONINC")
    
    # Sanity checks
    # Required args for all projections
    for name, var in viewitems({"MAP_PROJ" : map_proj, 
                                "STDLON" : stdlon, 
                                "REF_LAT" : ref_lat, 
                                "REF_LON" : ref_lon, 
                                "KNOWN_X" : known_x, 
                                "KNOWN_Y" : known_y,
                                "DX" : dx}):
        if var is None:
            raise ValueError("'{}' argument required".format(name))
    
    # Fortran wants 1-based indexing
    known_x = known_x + 1
    known_y = known_y + 1
    
    if map_proj in (ProjectionTypes.LAMBERT_CONFORMAL, 
                    ProjectionTypes.POLAR_STEREOGRAPHIC, 
                    ProjectionTypes.MERCATOR):
        if truelat1 is None:
            raise ValueError("'TRUELAT1' argument required")
    else:
        if truelat1 is None:
            truelat1 = 0.0
    
    # Map projection 6 (lat lon) required latinc, loninc, and dy
    if map_proj == ProjectionTypes.LAT_LON:
        if latinc is None:
            raise ValueError("'LATINC' argument required")
        
        if loninc is None:
            raise ValueError("'LONINC' argument required")
        
        if dy is None:
            raise ValueError("'DY' argument required")
    else:
        latinc = 0.0
        loninc = 0.0
        dy = 0.0
    
    return (map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon, pole_lat,
            pole_lon, known_x, known_y, dx, dy, latinc, loninc)


# Will return 0-based indexes
def _ll_to_xy(latitude, longitude, wrfnc=None, timeidx=0,
           stagger=None, method="cat", squeeze=True, cache=None,
           _key=None, as_int=True, **projparms):
    
    if wrfnc is not None:
        (map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon,
        pole_lat, pole_lon, known_x, known_y, dx, dy, latinc,
        loninc) = _get_proj_params(wrfnc, timeidx, stagger, method, squeeze, 
                                    cache, _key)
    else:
        (map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon,
        pole_lat, pole_lon, known_x, known_y, dx, dy, latinc,
        loninc) = _kwarg_proj_params(**projparams)
    
    if isinstance(latitude, Iterable):
        lats = np.asarray(latitude)
        lons = np.asarray(longitude)
        
        if lats.ndim > 1:
            lats = lats.ravel()
        
        if lons.ndim > 1:
            lons = lons.ravel()
        
        if (lats.size != lons.size):
            raise ValueError("'latitude' and 'longitude' "
                             "must be the same length")
        
        if ref_lat.size == 1:
            outdim = [2, lats.size]
            #outdim = [lats.size, 2]
            extra_dims = [outdim[1]]
        else:
            # Moving domain will have moving ref_lats/ref_lons
            outdim = [2, ref_lat.size, lats.size]
            #outdim = [lats.size, ref_lat.size, 2]
            extra_dims = outdim[1:]
            
        result = np.empty(outdim, np.float64)
        
        for left_idxs in iter_left_indexes(extra_dims):
            #left_and_slice_idxs = left_idxs + (slice(None), )
            # Left indexes is a misnomer, since these will be on the right
            x_idxs = (0,) + left_idxs
            y_idxs = (1,) + left_idxs
            if ref_lat.size == 1:
                ref_lat_val = ref_lat[0]
                ref_lon_val = ref_lon[0]
            else:
                ref_lat_val = ref_lat[left_idxs[-1]]
                ref_lon_val = ref_lon[left_idxs[-1]]
            
            lat = lats[left_idxs[0]]
            lon = lons[left_idxs[0]]
            
            xy = _lltoxy(map_proj, truelat1, truelat2, stdlon,
               ref_lat_val, ref_lon_val, pole_lat, pole_lon,
               known_x, known_y, dx, dy, latinc, loninc,
               lat, lon)
            
            # Note:  comes back from fortran as y,x
            result[x_idxs] = xy[1]
            result[y_idxs] = xy[0]
            
    else:
        result = np.empty((2,), np.float64)
        
        fort_out = _lltoxy(map_proj, truelat1, truelat2, stdlon,
               ref_lat, ref_lon, pole_lat, pole_lon,
               known_x, known_y, dx, dy, latinc, loninc,
               latitude, longitude)
        
        # Note, comes back from fortran as y,x.  So, need to swap them.
        result[0] = fort_out[1]
        result[1] = fort_out[0]
        
    
    # Make indexes 0-based
    result = result - 1
    
    if as_int:
        result = np.rint(result).astype(int)
    
    return result

# X and Y should be 0-based
def _xy_to_ll(x, y, wrfnc=None, timeidx=0, stagger=None, 
           method="cat", squeeze=True, cache=None, _key=None,
           **projparams):
        
    if wrfnc is not None:
        (map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon,
        pole_lat, pole_lon, known_x, known_y, dx, dy, latinc,
        loninc) = _get_proj_params(wrfnc, timeidx, stagger, method, squeeze, 
                                    cache, _key)
    else:
        (map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon,
        pole_lat, pole_lon, known_x, known_y, dx, dy, latinc,
        loninc) = _kwarg_proj_params(**projparams)
        
    
    if isinstance(x, Iterable):
        x_arr = np.asarray(x)
        y_arr = np.asarray(y)
        
        # Convert 0-based x, y to 1-based
        x_arr = x_arr + 1
        y_arr = y_arr + 1
        
        if x_arr.ndim > 1:
            x_arr = x_arr.ravel()
        
        if y_arr.ndim > 1:
            y_arr = y_arr.ravel()
        
        if (x_arr.size != y_arr.size):
            raise ValueError("'x' and 'y' must be the same length")
            
        if ref_lat.size == 1:
            #outdim = [x_arr.size, 2]
            #extra_dims = [outdim[0]]
            outdim = [2, x_arr.size]
            extra_dims = [outdim[1]]
        else:
            # Moving domain will have moving ref_lats/ref_lons
            #outdim = [x_arr.size, ref_lat.size, 2]
            #extra_dims = outdim[0:2]
            outdim = [2, ref_lat.size, x_arr.size]
            extra_dims = outdim[1:]
            
        result = np.empty(outdim, np.float64)
        
        for left_idxs in iter_left_indexes(extra_dims):
            #left_and_slice_idxs = left_idxs + (slice(None), )
            lat_idxs = (0,) + left_idxs
            lon_idxs = (1,) + left_idxs
            
            if ref_lat.size == 1:
                ref_lat_val = ref_lat[0]
                ref_lon_val = ref_lon[0]
            else:
                ref_lat_val = ref_lat[left_idxs[-1]]
                ref_lon_val = ref_lon[left_idxs[-1]]
            
            x_val = x_arr[left_idxs[0]]
            y_val = y_arr[left_idxs[0]]
            
            ll = _xytoll(map_proj, truelat1, truelat2, stdlon, ref_lat_val, 
                         ref_lon_val, pole_lat, pole_lon, known_x, known_y,
                         dx, dy, latinc, loninc, x_val, y_val)
            
            #result[left_and_slice_idxs] = ll[:]
            result[lat_idxs] = ll[0]
            result[lon_idxs] = ll[1]
            
    else:
        # Convert 0-based to 1-based for Fortran
        x_val = x + 1
        y_val = y + 1
        
        result = _xytoll(map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon,
                      pole_lat, pole_lon, known_x, known_y, dx, dy, latinc, 
                      loninc, x_val, y_val)
        
    return result
