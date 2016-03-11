from collections import Iterable

from wrf.var.constants import Constants
from wrf.var.extension import computeij, computell
from wrf.var.util import extract_vars, extract_global_attrs

__all__ = ["get_lat", "get_lon", "get_ij", "get_ll"]

def get_lat(wrfnc, timeidx=0):
    
    try:
        lat_vars = extract_vars(wrfnc, timeidx, varnames="XLAT")
    except KeyError:
        try:
            latm_vars = extract_vars(wrfnc, timeidx, varnames="XLAT_M")
        except:
            raise RuntimeError("Latitude variable not found in NetCDF file")
        else:
            xlat = latm_vars["XLAT_M"]
    else:
        xlat = lat_vars["XLAT"]
    
    return xlat
        
def get_lon(wrfnc, timeidx=0):
    try:
        lon_vars = extract_vars(wrfnc, timeidx, varnames="XLONG")
    except KeyError:
        try:
            lonm_vars = extract_vars(wrfnc, timeidx, varnames="XLONG_M")
        except:
            raise RuntimeError("Latitude variable not found in NetCDF file")
        else:
            xlon = lonm_vars["XLONG_M"]
    else:
        xlon = lon_vars["XLONG"]
        
    return xlon

def _get_proj_params(wrfnc, timeidx):
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
    
    xlat = get_lat(wrfnc, timeidx)
    xlon = get_lon(wrfnc, timeidx)
    
    ref_lat = xlat[0,0]
    ref_lon = xlon[0,0]
    
    known_i = 1.0
    known_j = 1.0
    
    return (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
             pole_lat,pole_lon,known_i,known_j,dx,latinc,
             loninc)
    

# TODO:  longitude and latitude can also be lists
def get_ij(wrfnc, longitude, latitude, timeidx=0):
    if isinstance(wrfnc, Iterable) and not isinstance(wrfnc, str):
        raise TypeError("'get_ij' is only applicabe for single files")
        
    (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
    pole_lat,pole_lon,known_i,known_j,dx,latinc,
    loninc) = _get_proj_params(wrfnc, timeidx)
    
    return computeij(map_proj,truelat1,truelat2,stdlon,
               ref_lat,ref_lon,pole_lat,pole_lon,
               known_i,known_j,dx,latinc,loninc,latitude,longitude)
    
# TODO:  i and j can also be lists
def get_ll(wrfnc, i, j, timeidx=0):
    if isinstance(wrfnc, Iterable) and not isinstance(wrfnc, str):
        raise TypeError("'get_ll' is only applicabe for single files")
    
    (map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
    pole_lat,pole_lon,known_i,known_j,dx,latinc,
    loninc) = _get_proj_params(wrfnc, timeidx)
    
    return computell(map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
             pole_lat,pole_lon,known_i,known_j,dx,latinc,
             loninc,i,j)
    
    
    