from .extension import computeavo, computepvo
from .util import extract_vars, extract_global_attrs
from .decorators import copy_and_set_metadata

__all__ = ["get_avo", "get_pvo"]

@copy_and_set_metadata(copy_varname="F", name="avo", 
                       description="absolute vorticity",
                       units="10-5 s-1")
def get_avo(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None):
    ncvars = extract_vars(wrfnc, timeidx, varnames=("U", "V", "MAPFAC_U",
                                              "MAPFAC_V", "MAPFAC_M",
                                              "F"),
                          method, squeeze, cache)
    
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


@copy_and_set_metadata(copy_varname="T", name="pvo", 
                       description="potential vorticity",
                       units="PVU")
def get_pvo(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None):
    ncvars = extract_vars(wrfnc, timeidx, varnames=("U", "V", "T", "P",
                                              "PB", "MAPFAC_U",
                                              "MAPFAC_V", "MAPFAC_M",
                                              "F"),
                          method, squeeze, cache)
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
    