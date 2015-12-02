from wrf.var.constants import Constants

from wrf.var.extension import computesrh, computeuh
from wrf.var.destagger import destagger

__all__ = ["get_srh", "get_uh"]

def get_srh(wrfnc, top=3000.0, timeidx=0):
    # Top can either be 3000 or 1000 (for 0-1 srh or 0-3 srh)
    
    if "U" in wrfnc.variables:
        u = destagger(wrfnc.variables["U"][timeidx,:,:,:], 2)
    elif "UU" in wrfnc.variables:
        u = destagger(wrfnc.variables["UU"][timeidx,:,:,:], 2) # support met_em files
        
    if "V" in wrfnc.variables:
        v = destagger(wrfnc.variables["V"][timeidx,:,:,:], 1)
    elif "VV" in wrfnc.variables:
        v = destagger(wrfnc.variables["VV"][timeidx,:,:,:], 1) 
    
    ter = wrfnc.variables["HGT"][timeidx,:,:]
    ph = wrfnc.variables["PH"][timeidx,:,:,:]
    phb = wrfnc.variables["PHB"][timeidx,:,:,:]
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, 0)
    
    z = geopt_unstag / Constants.G
    
    # Re-ordering from high to low
    u1 = u[::-1,:,:]
    v1 = v[::-1,:,:]
    z1 = z[::-1,:,:]
    
    srh = computesrh(u1, v1, z1, ter, top)
    
    return srh

def get_uh(wrfnc, bottom=2000.0, top=5000.0, timeidx=0):
    
    if "U" in wrfnc.variables:
            u = destagger(wrfnc.variables["U"][timeidx,:,:,:], 2)
    elif "UU" in wrfnc.variables:
        u = destagger(wrfnc.variables["UU"][timeidx,:,:,:], 2) # support met_em files
        
    if "V" in wrfnc.variables:
        v = destagger(wrfnc.variables["V"][timeidx,:,:,:], 1)
    elif "VV" in wrfnc.variables:
        v = destagger(wrfnc.variables["VV"][timeidx,:,:,:], 1) 
    
    wstag = wrfnc.variables["W"][timeidx,:,:,:]
    ph = wrfnc.variables["PH"][timeidx,:,:,:]
    phb = wrfnc.variables["PHB"][timeidx,:,:,:]
    zp = ph + phb
    
    mapfct = wrfnc.variables["MAPFAC_M"][timeidx,:,:]
    dx = wrfnc.getncattr("DX")
    dy = wrfnc.getncattr("DY")
       
       
    uh = computeuh(zp, mapfct, u, v, wstag, dx, dy, bottom, top)
    
    return uh

    
    
    
    