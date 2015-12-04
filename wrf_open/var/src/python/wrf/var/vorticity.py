from wrf.var.extension import computeavo, computepvo
from wrf.var.util import extract_vars, extract_global_attrs

__all__ = ["get_avo", "get_pvo"]

def get_avo(wrfnc, timeidx=0):
    ncvars = extract_vars(wrfnc, timeidx, vars=("U", "V", "MAPFAC_U",
                                              "MAPFAC_V", "MAPFAC_M",
                                              "F"))
    
    attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    u = ncvars["U"]
    v = ncvars["V"]
    msfu = ncvars["MAPFAC_U"]
    msfv = ncvars["MAPFAC_V"]
    msfm = ncvars["MAPFAC_M"]
    cor = ncvars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    return computeavo(u,v,msfu,msfv,msfm,cor,dx,dy)


def get_pvo(wrfnc, timeidx=0):
    ncvars = extract_vars(wrfnc, timeidx, vars=("U", "V", "T", "P",
                                              "PB", "MAPFAC_U",
                                              "MAPFAC_V", "MAPFAC_M",
                                              "F"))
    attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    
    u = ncvars["U"]
    v = ncvars["V"]
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    msfu = ncvars["MAPFAC_U"]
    msfv = ncvars["MAPFAC_V"]
    msfm = ncvars["MAPFAC_M"]
    cor = ncvars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    full_t = t + 300
    full_p = p + pb
    
    return computepvo(u,v,full_t,full_p,msfu,msfv,msfm,cor,dx,dy)
    