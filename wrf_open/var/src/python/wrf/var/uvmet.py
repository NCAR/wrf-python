from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from math import fabs, log, tan, sin, cos

import numpy as np

from .extension import computeuvmet
from .destag import destagger
from .constants import Constants
from .wind import _calc_wspd_wdir
from .decorators import convert_units
from .metadecorators import set_wind_metadata
from .util import extract_vars, extract_global_attrs, either

__all__=["get_uvmet", "get_uvmet10", "get_uvmet_wspd_wdir", 
         "get_uvmet10_wspd_wdir"]


@convert_units("wind", "mps")
def _get_uvmet(wrfnc, timeidx=0, method="cat", squeeze=True, 
               cache=None, meta=True,
               ten_m=False, units ="mps"):
    """ Return a tuple of u,v with the winds rotated in to earth space"""
    
    if not ten_m:
        varname = either("U", "UU")(wrfnc)
        u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                              meta=False)
        
        u = destagger(u_vars[varname], -1)
        
        varname = either("V", "VV")(wrfnc)
        v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                              meta=False)
        v = destagger(v_vars[varname], -2)
    else:
        varname = either("U10", "UU")(wrfnc)
        u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                              meta=False)
        u = (u_vars[varname] if varname == "U10" else 
             destagger(u_vars[varname][...,0,:,:], -1)) 
        
        varname = either("V10", "VV")(wrfnc)
        v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                              meta=False)
        v = (v_vars[varname] if varname == "V10" else 
             destagger(v_vars[varname][...,0,:,:], -2))
    
    map_proj_attrs = extract_global_attrs(wrfnc, attrs="MAP_PROJ")
    map_proj = map_proj_attrs["MAP_PROJ"]
    
    # 1 - Lambert
    # 2 - Polar Stereographic
    # 3 - Mercator
    # 6 - Lat/Lon
    # Note:  NCL has no code to handle other projections (0,4,5) as they 
    # don't appear to be supported any longer
    
    if map_proj in (0,3,6):
        # No rotation needed for Mercator and Lat/Lon, but still need
        # u,v aggregated in to one array
        
        end_idx = -3 if not ten_m else -2
        resdim = list(u.shape[0:end_idx]) + [2] + list(u.shape[end_idx:])
    
        # Make a new output array for the result
        res = np.empty(resdim, u.dtype)
        
        # For 2D array, this makes (...,0,:,:) and (...,1,:,:)
        # For 3D array, this makes (...,0,:,:,:) and (...,1,:,:,:)
        idx0 = tuple([Ellipsis] + [0] + [slice(None)]*(-end_idx))
        idx1 = tuple([Ellipsis] + [1] + [slice(None)]*(-end_idx))
    
        res[idx0] = u[:]
        res[idx1] = v[:]
        
        return res
    elif map_proj in (1,2):
        lat_attrs = extract_global_attrs(wrfnc, attrs=("TRUELAT1",
                                                       "TRUELAT2"))
        radians_per_degree = Constants.PI/180.0
        # Rotation needed for Lambert and Polar Stereographic
        true_lat1 = lat_attrs["TRUELAT1"]
        true_lat2 = lat_attrs["TRUELAT2"]
        
        try:
            lon_attrs = extract_global_attrs(wrfnc, attrs="STAND_LON")
        except AttributeError:
            try:
                cen_lon_attrs = extract_global_attrs(wrfnc, attrs="CEN_LON")
            except AttributeError:
                raise RuntimeError("longitude attributes not found in NetCDF")
            else:
                cen_lon = cen_lon_attrs["CEN_LON"]
        else:
            cen_lon = lon_attrs["STAND_LON"]
        
        
        varname = either("XLAT_M", "XLAT")(wrfnc)
        xlat_var = extract_vars(wrfnc, timeidx, varname, 
                                method, squeeze, cache, meta=False)
        lat = xlat_var[varname]
        
        varname = either("XLONG_M", "XLONG")(wrfnc)
        xlon_var = extract_vars(wrfnc, timeidx, varname, 
                                method, squeeze, cache, meta=False)
        lon = xlon_var[varname]
        
        if map_proj == 1:
            if((fabs(true_lat1 - true_lat2) > 0.1) and
                    (fabs(true_lat2 - 90.) > 0.1)): 
                cone = (log(cos(true_lat1*radians_per_degree)) 
                    - log(cos(true_lat2*radians_per_degree)))
                cone = (cone / 
                        (log(tan((45.-fabs(true_lat1/2.))*radians_per_degree)) 
                    - log(tan((45.-fabs(true_lat2/2.))*radians_per_degree)))) 
            else:
                cone = sin(fabs(true_lat1)*radians_per_degree)
        else:
            cone = 1
        
        res = computeuvmet(u, v ,lat, lon, cen_lon, cone)
        
        return res

    
@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="uvmet",
                   description="earth rotated u,v",
                   two_d=False, 
                   wspd_wdir=False)      
def get_uvmet(wrfnc, timeidx=0, method="cat", squeeze=True, 
              cache=None, meta=True,
              units="mps"):
    
    return _get_uvmet(wrfnc, timeidx, method, squeeze, cache, meta,
                      False, units)


@set_wind_metadata(copy_varname=either("PSFC", "F"), 
                   name="uvmet10",
                   description="10m earth rotated u,v",
                   two_d=True, 
                   wspd_wdir=False)
def get_uvmet10(wrfnc, timeidx=0, method="cat", squeeze=True, 
                cache=None, meta=True,
                units="mps"):
    
    return _get_uvmet(wrfnc, timeidx, method, squeeze, cache, meta, 
                      True, units)


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="uvmet_wspd_wdir",
                   description="earth rotated wspd,wdir",
                   two_d=False, 
                   wspd_wdir=True)
def get_uvmet_wspd_wdir(wrfnc, timeidx=0, method="cat", squeeze=True, 
                        cache=None, meta=True,
                        units="mps"):
    
    uvmet = _get_uvmet(wrfnc, timeidx, False, units, method, squeeze, cache,
                       meta)
    
    return _calc_wspd_wdir(uvmet[...,0,:,:,:], uvmet[...,1,:,:,:], False, units)
    

@set_wind_metadata(copy_varname=either("PSFC", "F"), 
                   name="uvmet10_wspd_wdir",
                   description="10m earth rotated wspd,wdir",
                   two_d=True, 
                   wspd_wdir=True)
def get_uvmet10_wspd_wdir(wrfnc, timeidx=0, method="cat", squeeze=True, 
                          cache=None, meta=True,
                          units="mps"):
    
    uvmet10 = _get_uvmet(wrfnc, timeidx, method, squeeze, cache, meta, 
                         True, units)
    
    return _calc_wspd_wdir(uvmet10[...,0,:,:], uvmet10[...,1,:,:], True, units)
    

            
        
            
            