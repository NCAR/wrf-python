from math import fabs, log, tan, sin, cos

from .extension import computeuvmet
from .destag import destagger
from .constants import Constants
from .wind import _calc_wspd_wdir
from .decorators import convert_units, set_wind_metadata
from .util import extract_vars, extract_global_attrs, either

__all__=["get_uvmet", "get_uvmet10", "get_uvmet_wspd_wdir", 
         "get_uvmet10_wspd_wdir"]

@set_wind_metadata(wspd_wdir=False)
@convert_units("wind", "mps")
def get_uvmet(wrfnc, timeidx=0, ten_m=False, units ="mps", 
              method="cat", squeeze=True, cache=None):
    """ Return a tuple of u,v with the winds rotated in to earth space"""
    
    if not ten_m:
        varname = either("U", "UU")(wrfnc)
        u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache)
        u = destagger(u_vars[varname], -1)
        
        varname = either("V", "VV")(wrfnc)
        v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache)
        v = destagger(v_vars[varname], -2)
    else:
        varname = either("U10", "UU")
        u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache)
        u = (u_vars[varname] if varname == "U10" else 
             destagger(u_vars[varname][...,0,:,:], -1)) 
        
        varname = either("V10", "VV")
        v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache)
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
        # No rotation needed for Mercator and Lat/Lon
        return u,v
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
                                method, squeeze, cache)
        lat = xlat_var[varname]
        
        varname = either("XLONG_M", "XLONG")
        xlon_var = extract_vars(wrfnc, timeidx, varname, 
                                method, squeeze, cache)
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
        
        res = computeuvmet(u,v,lat,lon,cen_lon,cone)
        
        return res
            

def get_uvmet10(wrfnc, timeidx=0, units="mps",
                method="cat", squeeze=True, cache=None):
    return get_uvmet(wrfnc, timeidx, True, units)

@set_wind_metadata(wspd_wdir=True)
def get_uvmet_wspd_wdir(wrfnc, timeidx=0, units="mps",
                        method="cat", squeeze=True, cache=None):
    
    uvmet = get_uvmet(wrfnc, timeidx, False, units, method, squeeze, cache)
    
    return _calc_wspd_wdir(uvmet[0,:], uvmet[1,:], False, units)
    

@set_wind_metadata(wspd_wdir=True)
def get_uvmet10_wspd_wdir(wrfnc, timeidx=0, units="mps",
                          method="cat", squeeze=True, cache=None):
    
    uvmet10 = get_uvmet10(wrfnc, timeidx, units="mps", method, squeeze, cache)
    
    return _calc_wspd_wdir(uvmet10[0,:], uvmet10[1,:], True, units)
    

            
        
            
            