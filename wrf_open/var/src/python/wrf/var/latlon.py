from wrf.var.extension import computeij, computell

__all__ = ["get_lat", "get_lon", "get_ij", "get_ll"]

def get_lat(wrfnc, timeidx=0):
    if "XLAT" in wrfnc.variables:
        xlat = wrfnc.variables["XLAT"][timeidx,:,:]
    elif "XLAT_M" in wrfnc.variables:
        xlat = wrfnc.variables["XLAT_M"][timeidx,:,:]
    
    return xlat
        
def get_lon(wrfnc, timeidx=0):
    if "XLONG" in wrfnc.variables:
        xlon = wrfnc.variables["XLONG"][timeidx,:,:]
    elif "XLONG_M" in wrfnc.variables:
        xlon = wrfnc.variables["XLONG_M"][timeidx,:,:]
        
    return xlon

def get_ij(wrfnc, longitude, latitude, timeidx=0):
    map_proj = wrfnc.getncattr("MAP_PROJ")
    truelat1 = wrfnc.getncattr("TRUELAT1")
    truelat2 = wrfnc.getncattr("TRUELAT2")
    stdlon = wrfnc.getncattr("STAND_LON")
    dx = wrfnc.getncattr("DX")
    dy = wrfnc.getncattr("DY")
    stdlon = wrfnc.getncattr("STAND_LON")
    
    if map_proj == 6:
        pole_lat = wrfnc.getncattr("POLE_LAT")
        pole_lon = wrfnc.getncattr("POLE_LON")
        latinc = (dy*360.0)/2.0/3.141592653589793/6370000.
        loninc = (dx*360.0)/2.0/3.141592653589793/6370000.
    else:
        pole_lat = 90.0
        pole_lon = 0.0
        latinc = 0.0
        loninc = 0.0
    
    if "XLAT" in wrfnc.variables:
        xlat = wrfnc.variables["XLAT"][timeidx,:,:]
    elif "XLAT_M" in wrfnc.variables:
        xlat = wrfnc.variables["XLAT_M"][timeidx,:,:]
        
    if "XLONG" in wrfnc.variables:
        xlon = wrfnc.variables["XLONG"][timeidx,:,:]
    elif "XLONG_M" in wrfnc.variables:
        xlon = wrfnc.variables["XLONG_M"][timeidx,:,:]
    
    ref_lat = xlat[0,0]
    ref_lon = xlon[0,0]
    
    known_i = 1.0
    known_j = 1.0
    
    return computeij(map_proj,truelat1,truelat2,stdlon,
               ref_lat,ref_lon,pole_lat,pole_lon,
               known_i,known_j,dx,latinc,loninc,latitude,longitude)
    

def get_ll(wrfnc, i, j, timeidx=0):
    
    map_proj = wrfnc.getncattr("MAP_PROJ")
    truelat1 = wrfnc.getncattr("TRUELAT1")
    truelat2 = wrfnc.getncattr("TRUELAT2")
    stdlon = wrfnc.getncattr("STAND_LON")
    dx = wrfnc.getncattr("DX")
    dy = wrfnc.getncattr("DY")
    stdlon = wrfnc.getncattr("STAND_LON")
    
    if map_proj == 6:
        pole_lat = wrfnc.getncattr("POLE_LAT")
        pole_lon = wrfnc.getncattr("POLE_LON")
        latinc = (dy*360.0)/2.0/3.141592653589793/6370000.
        loninc = (dx*360.0)/2.0/3.141592653589793/6370000.
    else:
        pole_lat = 90.0
        pole_lon = 0.0
        latinc = 0.0
        loninc = 0.0
    
    if "XLAT" in wrfnc.variables:
        xlat = wrfnc.variables["XLAT"][timeidx,:,:]
    elif "XLAT_M" in wrfnc.variables:
        xlat = wrfnc.variables["XLAT_M"][timeidx,:,:]
        
    if "XLONG" in wrfnc.variables:
        xlon = wrfnc.variables["XLONG"][timeidx,:,:]
    elif "XLONG_M" in wrfnc.variables:
        xlon = wrfnc.variables["XLONG_M"][timeidx,:,:]
    
    ref_lat = xlat[0,0]
    ref_lon = xlon[0,0]
    
    known_i = 1.0
    known_j = 1.0
    
    return computell(map_proj,truelat1,truelat2,stdlon,ref_lat,ref_lon,
             pole_lat,pole_lon,known_i,known_j,dx,latinc,
             loninc,i,j)
    
    
    