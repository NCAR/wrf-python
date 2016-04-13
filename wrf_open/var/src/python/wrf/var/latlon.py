from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from collections import Iterable

import numpy as np

from .config import xarray_enabled
from .constants import Constants
from .extension import computeij, computell
from .util import (extract_vars, extract_global_attrs, 
                   either, _is_moving_domain, _is_multi_time_req,
                   iter_left_indexes)
from .metadecorators import set_latlon_metadata

if xarray_enabled():
    from xarray import DataArray

__all__ = ["get_lat", "get_lon", "get_ij", "get_ll"]

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

def get_lat(wrfnc, timeidx=0, stagger=None,
            method="cat", squeeze=True, cache=None):
    
    varname = _lat_varname(wrfnc, stagger)
    lat_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           nometa=False)
    
    return lat_var[varname]
        
def get_lon(wrfnc, timeidx=0, stagger=None,
            method="cat", squeeze=True, cache=None):
    
    varname = _lon_varname(wrfnc, stagger)
    lon_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           nometa=False)
    
    return lon_var[varname]

def _get_proj_params(wrfnc, timeidx, stagger, method, squeeze, cache):
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
    
    if map_proj == 6:
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
    
    latvar = _lat_varname(stagger)
    lonvar = _lon_varname(stagger)
    
    lat_timeidx = timeidx
    # Only need all the lats/lons if it's a moving domain file/files
    if _is_multi_time_req(timeidx):
        if not _is_moving_domain(wrfnc, latvar=latvar, lonvar=lonvar):
            lat_timeidx = 0
        
    xlat = get_lat(wrfnc, lat_timeidx, stagger, method, squeeze, cache)
    xlon = get_lon(wrfnc, lat_timeidx, stagger, method, squeeze, cache)
    
    ref_lat = np.ravel(xlat[...,0,0])
    ref_lon = np.ravel(xlon[...,0,0])
    
    # Note: fortran index
    known_i = 1.0
    known_j = 1.0
    
    return (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
             pole_lat,pole_lon,known_i,known_j,dx,latinc,
             loninc)
    

@set_latlon_metadata(ij=True) 
def get_ij(wrfnc, latitude, longitude, timeidx=0,
           stagger=None, method="cat", squeeze=True, cache=None):
    
    (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
    pole_lat,pole_lon,known_i,known_j,dx,latinc,
    loninc) = _get_proj_params(wrfnc, timeidx)
    
    if isinstance(wrfnc, Iterable) and not isinstance(wrfnc, str):
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
            outdim = [lats.size, 2]
            extra_dims = outdim[0]
        else:
            # Moving domain will have moving ref_lats/ref_lons
            outdim = [lats.size, ref_lat.size, 2]
            extra_dims = outdim[0:2]

        res = np.empty(outdim, np.float64)
        
        for left_idxs in iter_left_indexes(extra_dims):
            left_and_slice_idxs = left_idxs + [slice(None, None, None)]
            
            if ref_lat.size == 1:
                ref_lat_val = ref_lat[0]
                ref_lon_val = ref_lon[0]
            else:
                ref_lat_val = ref_lat[left_idxs[-1]]
                ref_lon_val = ref_lon[left_idxs[-1]]
            
            lat = lats[left_idxs[0]]
            lon = lons[left_idxs[0]]
            
            ij = computeij(map_proj, truelat1, truelat2, stdlon,
               ref_lat_val, ref_lon_val, pole_lat, pole_lon,
               known_i, known_j, dx, latinc, loninc,
               lat, lon)
            
            res[left_and_slice_idxs] = ij[:]
            
    else:
        
        res = computeij(map_proj, truelat1, truelat2, stdlon,
               ref_lat, ref_lon, pole_lat, pole_lon,
               known_i, known_j, dx, latinc, loninc,
               latitude, longitude)
    
    return res


@set_latlon_metadata(ij=False)
def get_ll(wrfnc, i, j, timeidx=0):
    if isinstance(wrfnc, Iterable) and not isinstance(wrfnc, str):
        raise TypeError("'get_ll' is only applicabe for single files")
    
    (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
    pole_lat,pole_lon,known_i,known_j,dx,latinc,
    loninc) = _get_proj_params(wrfnc, timeidx)
    
    if isinstance(wrfnc, Iterable) and not isinstance(wrfnc, str):
        i_arr = np.asarray(i)
        j_arr = np.asarray(j)
        
        if i_arr.ndim > 1:
            i_arr = i_arr.ravel()
        
        if j_arr.ndim > 1:
            j_arr = j_arr.ravel()
        
        if (i_arr.size != j_arr.size):
            raise ValueError("'i' and 'j' "
                             "must be the same length")
            
        if ref_lat.size == 1:
            outdim = [i_arr.size, 2]
            extra_dims = outdim[0]
        else:
            # Moving domain will have moving ref_lats/ref_lons
            outdim = [i_arr.size, ref_lat.size, 2]
            extra_dims = outdim[0:2]

        res = np.empty(outdim, np.float64)
        
        for left_idxs in iter_left_indexes(extra_dims):
            left_and_slice_idxs = left_idxs + [slice(None, None, None)]
            
            if ref_lat.size == 1:
                ref_lat_val = ref_lat[0]
                ref_lon_val = ref_lon[0]
            else:
                ref_lat_val = ref_lat[left_idxs[-1]]
                ref_lon_val = ref_lon[left_idxs[-1]]
            
            i_val = i_arr[left_idxs[0]]
            j_val = j_arr[left_idxs[0]]
            
            ll = computell(map_proj, truelat1, truelat2,
                           stdlon, ref_lat_val, ref_lon_val,
                           pole_lat, pole_lon, known_i, known_j,
                           dx, latinc, loninc,
                           i_val, j_val)
            
            res[left_and_slice_idxs] = ll[:]
            
    else:
        res = computell(map_proj, truelat1, truelat2,
                        stdlon, ref_lat, ref_lon,
                        pole_lat, pole_lon, known_i, known_j,
                        dx, latinc, loninc,
                        i_val, j_val)
        
    return res
    
    
    