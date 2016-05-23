from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants

from .extension import computesrh, computeuh
from .destag import destagger
from .util import extract_vars, extract_global_attrs, either
from .metadecorators import copy_and_set_metadata

__all__ = ["get_srh", "get_uh"]

@copy_and_set_metadata(copy_varname="HGT", name="srh", 
                       description="storm relative helicity",
                       units="m-2/s-2")
def get_srh(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True,
            top=3000.0):
    # Top can either be 3000 or 1000 (for 0-1 srh or 0-3 srh)
    
    ncvars = extract_vars(wrfnc, timeidx, ("HGT", "PH", "PHB"),
                          method, squeeze, cache, meta=False)
    
    ter = ncvars["HGT"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    
    # As coded in NCL, but not sure this is possible
    varname = either("U", "UU")(wrfnc)
    u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False)
    u = destagger(u_vars[varname], -1) 
    
    varname = either("V", "VV")(wrfnc)
    v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False)
    v = destagger(v_vars[varname], -2)

    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    
    z = geopt_unstag / Constants.G
    
    # Re-ordering from high to low
    u1 = u[...,::-1,:,:] 
    v1 = v[...,::-1,:,:]
    z1 = z[...,::-1,:,:]
    
    srh = computesrh(u1, v1, z1, ter, top)
    
    return srh

@copy_and_set_metadata(copy_varname="MAPFAC_M", name="updraft_helicity", 
                       description="updraft helicity",
                       units="m-2/s-2")
def get_uh(wrfnc, timeidx=0, method="cat", squeeze=True, 
           cache=None, meta=True,
           bottom=2000.0, top=5000.0):
    
    ncvars = extract_vars(wrfnc, timeidx, ("W", "PH", "PHB", "MAPFAC_M"),
                          method, squeeze, cache)
    
    wstag = ncvars["W"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    mapfct = ncvars["MAPFAC_M"]
    
    attrs  = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    # As coded in NCL, but not sure this is possible
    varname = either("U", "UU")(wrfnc)
    u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False)
    u = destagger(u_vars[varname], -1) 
    
    varname = either("V", "VV")(wrfnc)
    v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False)
    v = destagger(v_vars[varname], -2) 
    
    zp = ph + phb
    
    uh = computeuh(zp, mapfct, u, v, wstag, dx, dy, bottom, top)
    
    return uh

    
    
    
    