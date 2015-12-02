from math import fabs, log, tan, sin, cos

from wrf.var.extension import computeuvmet
from wrf.var.destagger import destagger
from wrf.var.constants import Constants
from wrf.var.wind import _calc_wspd_wdir
from wrf.var.decorators import convert_units

__all__=["get_uvmet", "get_uvmet10", "get_uvmet_wspd_wdir", 
         "get_uvmet10_wspd_wdir"]

@convert_units("wind", "mps")
def get_uvmet(wrfnc, ten_m=False, units ="mps", timeidx=0):
    """ Return a tuple of u,v with the winds rotated in to earth space"""
    
    if not ten_m:
        if "U" in wrfnc.variables:
            u = destagger(wrfnc.variables["U"][timeidx,:,:,:],  2)
        elif "UU" in wrfnc.variables:
            u = destagger(wrfnc.variables["UU"][timeidx,:,:,:], 2) # support met_em files
            
        if "V" in wrfnc.variables:
            v = destagger(wrfnc.variables["V"][timeidx,:,:,:], 1)
        elif "VV" in wrfnc.variables:
            v = destagger(wrfnc.variables["VV"][timeidx,:,:,:], 1) 
    else:
        if "U10" in wrfnc.variables:
            u = wrfnc.variables["U10"][timeidx,:,:]
        elif "UU" in wrfnc.variables:
            u = destagger(wrfnc.variables["UU"][timeidx,0,:,:], 1) # support met_em files
        
        if "V10" in wrfnc.variables:
            v = wrfnc.variables["V10"][timeidx,:,:]
        elif "VV" in wrfnc.variables:
            v = destagger(wrfnc.variables["VV"][timeidx,0,:,:], 0) # support met_em files
        
    map_proj = wrfnc.getncattr("MAP_PROJ")
    
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
        radians_per_degree = Constants.PI/180.0
        # Rotation needed for Lambert and Polar Stereographic
        cen_lat  = wrfnc.getncattr("CEN_LAT")
        if "STAND_LON" in wrfnc.ncattrs():
            cen_lon = wrfnc.getncattr("STAND_LON")
        else:
            cen_lon = wrfnc.getncattr("CEN_LON")
            
        true_lat1 = wrfnc.getncattr("TRUELAT1")
        true_lat2 = wrfnc.getncattr("TRUELAT2")
        
        if "XLAT_M" in wrfnc.variables:
            lat = wrfnc.variables["XLAT_M"][timeidx,:,:]
        else:
            lat = wrfnc.variables["XLAT"][timeidx,:,:]
            
        if "XLONG_M" in wrfnc.variables:
            lon = wrfnc.variables["XLONG_M"][timeidx,:,:]
        else:
            lon = wrfnc.variables["XLONG"][timeidx,:,:]
            
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
        
        if u.ndim == 3:
            return res
        else:
            return res[:,0,:,:]
            
    
def get_uvmet10(wrfnc, units="mps", timeidx=0):
    return get_uvmet(wrfnc, True, units, timeidx)

def get_uvmet_wspd_wdir(wrfnc, units="mps", timeidx=0):
    u,v = get_uvmet(wrfnc, False, units, timeidx)
    return _calc_wspd_wdir(u, v, units)

def get_uvmet10_wspd_wdir(wrfnc, units="mps", timeidx=0):
    u,v = get_uvmet10(wrfnc, units="mps", timeidx=0)
    return _calc_wspd_wdir(u, v, units)

            
        
            
            