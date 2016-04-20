from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from collections import Iterable

import numpy as np

from .constants import Constants
from .extension import computeij, computell
from .util import (extract_vars, extract_global_attrs, 
                   either, _is_moving_domain, _is_multi_time_req,
                   iter_left_indexes, _is_mapping, _is_multi_file,
                   viewkeys)

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
    
    latvar = _lat_varname(wrfnc, stagger)
    lonvar = _lon_varname(wrfnc, stagger)
    
    lat_timeidx = timeidx
    
    is_moving = _is_moving_domain(wrfnc, latvar=latvar, lonvar=lonvar)
    
    # Only need one file and one time if the domain is not moving
    if not is_moving:
        if _is_multi_time_req(timeidx):
            lat_timeidx = 0
            
        if _is_multi_file(wrfnc):
            if not _is_mapping(wrfnc):
                wrfnc = next(iter(wrfnc))  # only need one file
            else:
                wrfnc = wrfnc[next(iter(viewkeys(wrfnc)))]
                return _get_proj_params(wrfnc, timeidx, stagger, 
                                        method, squeeze, cache)
            
    xlat = extract_vars(wrfnc, lat_timeidx, (latvar,), method, squeeze, cache,
                           nometa=True)[latvar]
    xlon = extract_vars(wrfnc, lat_timeidx, (lonvar,), method, squeeze, cache,
                           nometa=True)[lonvar]
    
    ref_lat = np.ravel(xlat[...,0,0])
    ref_lon = np.ravel(xlon[...,0,0])
    
    # Note: fortran index
    known_i = 1.0
    known_j = 1.0
    
    return (map_proj, truelat1, truelat2, stdlon, ref_lat, ref_lon,
            pole_lat, pole_lon, known_i, known_j, dx, latinc, loninc)
    

def ll_to_ij(wrfnc, latitude, longitude, timeidx=0,
           stagger=None, method="cat", squeeze=True, cache=None):
    
    (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
    pole_lat,pole_lon,known_i,known_j,dx,latinc,
    loninc) = _get_proj_params(wrfnc, timeidx, stagger, method, squeeze, cache)
    
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
            outdim = [lats.size, 2]
            extra_dims = [outdim[0]]
        else:
            # Moving domain will have moving ref_lats/ref_lons
            outdim = [lats.size, ref_lat.size, 2]
            extra_dims = outdim[0:2]
            
        res = np.empty(outdim, np.float64)
        
        for left_idxs in iter_left_indexes(extra_dims):
            left_and_slice_idxs = left_idxs + (slice(None), )
            
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
    
    if squeeze:
        res = res.squeeze()
        
    return res

def ij_to_ll(wrfnc, i, j, timeidx=0,
           stagger=None, method="cat", squeeze=True, cache=None):
    
    (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
    pole_lat,pole_lon,known_i,known_j,dx,latinc,
    loninc) = _get_proj_params(wrfnc, timeidx, stagger, method, squeeze, cache)
    
    if isinstance(i, Iterable):
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
            extra_dims = [outdim[0]]
        else:
            # Moving domain will have moving ref_lats/ref_lons
            outdim = [i_arr.size, ref_lat.size, 2]
            extra_dims = outdim[0:2]
            
        res = np.empty(outdim, np.float64)
        
        for left_idxs in iter_left_indexes(extra_dims):
            left_and_slice_idxs = left_idxs + (slice(None), )
            
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
        i_val = i
        j_val = j
        
        res = computell(map_proj, truelat1, truelat2,
                        stdlon, ref_lat, ref_lon,
                        pole_lat, pole_lon, known_i, known_j,
                        dx, latinc, loninc,
                        i_val, j_val)
    
    if squeeze:
        res = res.squeeze()
        
    return res