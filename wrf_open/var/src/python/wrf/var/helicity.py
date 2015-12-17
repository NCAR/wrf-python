from wrf.var.constants import Constants

from wrf.var.extension import computesrh, computeuh
from wrf.var.destag import destagger
from wrf.var.util import extract_vars, extract_global_attrs

__all__ = ["get_srh", "get_uh"]

def get_srh(wrfnc, timeidx=0, top=3000.0):
    # Top can either be 3000 or 1000 (for 0-1 srh or 0-3 srh)
    
    ncvars = extract_vars(wrfnc, timeidx, vars=("HGT", "PH", "PHB"))
    
    ter = ncvars["HGT"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    
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

    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    
    z = geopt_unstag / Constants.G
    
    # Re-ordering from high to low
    u1 = u[...,::-1,:,:] 
    v1 = v[...,::-1,:,:]
    z1 = z[...,::-1,:,:]
    
    srh = computesrh(u1, v1, z1, ter, top)
    
    return srh

def get_uh(wrfnc, timeidx=0, bottom=2000.0, top=5000.0):
    
    ncvars = extract_vars(wrfnc, timeidx, vars=("W", "PH", "PHB", "MAPFAC_M"))
    
    wstag = ncvars["W"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    mapfct = ncvars["MAPFAC_M"]
    
    attrs  = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    dx = attrs["DX"]
    dy = attrs["DY"]
    
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
    
    zp = ph + phb
    
    uh = computeuh(zp, mapfct, u, v, wstag, dx, dy, bottom, top)
    
    return uh

    
    
    
    