from wrf.var.extension import computeavo, computepvo
from wrf.var.util import extract_vars, extract_global_attrs

__all__ = ["get_avo", "get_pvo"]

def get_avo(wrfnc, timeidx=0):
    vars = extract_vars(wrfnc, timeidx, vars=("U", "V", "MAPFAC_U",
                                              "MAPFAC_V", "MAPFAC_M",
                                              "F"))
    
    attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    u = vars["U"]
    v = vars["V"]
    msfu = vars["MAPFAC_U"]
    msfv = vars["MAPFAC_V"]
    msfm = vars["MAPFAC_M"]
    cor = vars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    return computeavo(u,v,msfu,msfv,msfm,cor,dx,dy)


def get_pvo(wrfnc, timeidx=0):
    vars = extract_vars(wrfnc, timeidx, vars=("U", "V", "T", "P",
                                              "PB", "MAPFAC_U",
                                              "MAPFAC_V", "MAPFAC_M",
                                              "F"))
    attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    
    u = vars["U"]
    v = vars["V"]
    t = vars["T"]
    p = vars["P"]
    pb = vars["PB"]
    msfu = vars["MAPFAC_U"]
    msfv = vars["MAPFAC_V"]
    msfm = vars["MAPFAC_M"]
    cor = vars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    full_t = t + 300
    full_p = p + pb
    
    return computepvo(u,v,full_t,full_p,msfu,msfv,msfm,cor,dx,dy)
    