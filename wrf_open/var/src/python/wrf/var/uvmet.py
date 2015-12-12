from math import fabs, log, tan, sin, cos

from wrf.var.extension import computeuvmet
from wrf.var.destagger import destagger
from wrf.var.constants import Constants
from wrf.var.wind import _calc_wspd_wdir
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars, extract_global_attrs

__all__=["get_uvmet", "get_uvmet10", "get_uvmet_wspd_wdir", 
         "get_uvmet10_wspd_wdir"]

@convert_units("wind", "mps")
def get_uvmet(wrfnc, ten_m=False, units ="mps", timeidx=0):
    """ Return a tuple of u,v with the winds rotated in to earth space"""
    
    if not ten_m:
        try:
            u_vars = extract_vars(wrfnc, timeidx, vars="U")
        except KeyError:
            try:
                uu_vars = extract_vars(wrfnc, timeidx, vars="UU")
            except KeyError:
                raise RuntimeError("No valid wind data found in NetCDF file")
            else:
                u = destagger(uu_vars["UU"], -1) # support met_em files
        else:
            u = destagger(u_vars["U"], -1)   
        
        try:
            v_vars = extract_vars(wrfnc, timeidx, vars="V")
        except KeyError:
            try:
                vv_vars = extract_vars(wrfnc, timeidx, vars="VV")
            except KeyError:
                raise RuntimeError("No valid wind data found in NetCDF file")
            else:
                v = destagger(vv_vars["VV"], -2) # support met_em files
        else:
            v = destagger(v_vars["V"], -2) 
            
    else:
        try:
            u_vars = extract_vars(wrfnc, timeidx, vars="U10")
        except KeyError:
            try:
                uu_vars = extract_vars(wrfnc, timeidx, vars="UU")
            except KeyError:
                raise RuntimeError("No valid wind data found in NetCDF file")
            else:
                # For met_em files, this just takes the lowest level winds
                # (3rd dimension from right is bottom_top)
                u = destagger(uu_vars["UU"][...,0,:,:], -1) # support met_em files
        else:
            u = u_vars["U10"] 
        
        try:
            v_vars = extract_vars(wrfnc, timeidx, vars="V10")
        except KeyError:
            try:
                vv_vars = extract_vars(wrfnc, timeidx, vars="VV")
            except KeyError:
                raise RuntimeError("No valid wind data found in NetCDF file")
            else:
                # For met_em files, this just takes the lowest level winds
                # (3rd dimension from right is bottom_top)
                v = destagger(vv_vars["VV"][...,0,:,:], -2) # support met_em files
        else:
            v = v_vars["V10"]
    
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
        lat_attrs = extract_global_attrs(wrfnc, attrs=("CEN_LAT",
                                                        "TRUELAT1",
                                                        "TRUELAT2"))
        radians_per_degree = Constants.PI/180.0
        # Rotation needed for Lambert and Polar Stereographic
        cen_lat  = lat_attrs["CEN_LAT"]
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
              
        try:
            xlat_m_vars = extract_vars(wrfnc, timeidx, vars="XLAT_M")
        except KeyError:
            try:
                xlat_vars = extract_vars(wrfnc, timeidx, vars="XLAT")
            except KeyError:
                raise RuntimeError("xlat not found in NetCDF file")
            else:
                lat = xlat_vars["XLAT"]
        else:
            lat = xlat_m_vars["XLAT_M"]
            
        try:
            xlon_m_vars = extract_vars(wrfnc, timeidx, vars="XLONG_M")
        except KeyError:
            try:
                xlon_vars = extract_vars(wrfnc, timeidx, vars="XLONG")
            except KeyError:
                raise RuntimeError("xlong not found in NetCDF file")
            else:
                lon = xlon_vars["XLONG"]
        else:
            lon = xlon_m_vars["XLONG_M"]
            
        if map_proj == 1:
            if((fabs(true_lat1 - true_lat2) > 0.1) and
                    (fabs(true_lat2 - 90.) > 0.1)): 
                cone = (log(cos(true_lat1*radians_per_degree)) 
                        - log(cos(true_lat2*radians_per_degree)))
                cone = cone / (log(tan((45.-fabs(true_lat1/2.))*radians_per_degree)) 
                        - log(tan((45.-fabs(true_lat2/2.))*radians_per_degree))) 
            else:
                cone = sin(fabs(true_lat1)*radians_per_degree)
        else:
            cone = 1
        
        res = computeuvmet(u,v,lat,lon,cen_lon,cone)
        
        return res
            
    
def get_uvmet10(wrfnc, units="mps", timeidx=0):
    return get_uvmet(wrfnc, True, units, timeidx)

def get_uvmet_wspd_wdir(wrfnc, units="mps", timeidx=0):
    u,v = get_uvmet(wrfnc, False, units, timeidx)
    return _calc_wspd_wdir(u, v, units)

def get_uvmet10_wspd_wdir(wrfnc, units="mps", timeidx=0):
    u,v = get_uvmet10(wrfnc, units="mps", timeidx=0)
    return _calc_wspd_wdir(u, v, units)

            
        
            
            