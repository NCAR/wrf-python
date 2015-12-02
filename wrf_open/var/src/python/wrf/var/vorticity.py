from wrf.var.extension import computeavo, computepvo

__all__ = ["get_avo", "get_pvo"]

def get_avo(wrfnc, timeidx=0):
    u = wrfnc.variables["U"][timeidx,:,:,:]
    v = wrfnc.variables["V"][timeidx,:,:,:]
    msfu = wrfnc.variables["MAPFAC_U"][timeidx,:,:]
    msfv = wrfnc.variables["MAPFAC_V"][timeidx,:,:]
    msfm = wrfnc.variables["MAPFAC_M"][timeidx,:,:]
    cor = wrfnc.variables["F"][timeidx,:,:]
    dx = wrfnc.getncattr("DX")
    dy = wrfnc.getncattr("DY")
    
    return computeavo(u,v,msfu,msfv,msfm,cor,dx,dy)


def get_pvo(wrfnc, timeidx=0):
    u = wrfnc.variables["U"][timeidx,:,:,:]
    v = wrfnc.variables["V"][timeidx,:,:,:]
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    msfu = wrfnc.variables["MAPFAC_U"][timeidx,:,:]
    msfv = wrfnc.variables["MAPFAC_V"][timeidx,:,:]
    msfm = wrfnc.variables["MAPFAC_M"][timeidx,:,:]
    cor = wrfnc.variables["F"][timeidx,:,:]
    dx = wrfnc.getncattr("DX")
    dy = wrfnc.getncattr("DY")
    
    full_t = t + 300
    full_p = p + pb
    
    return computepvo(u,v,full_t,full_p,msfu,msfv,msfm,cor,dx,dy)
    