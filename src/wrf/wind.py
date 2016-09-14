from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

from .constants import Constants
from .destag import destagger
from .util import extract_vars, either
from .decorators import convert_units
from .metadecorators import set_wind_metadata


@convert_units("wind", "m s-1")
def _calc_wspd(u, v, units="m s-1"):
    return np.sqrt(u**2 + v**2)


def _calc_wdir(u, v):
    wdir = 270.0 - np.arctan2(v,u) * (180.0/Constants.PI)
    return np.remainder(wdir, 360.0)


def _calc_wspd_wdir(u, v, two_d, units):
    wspd = _calc_wspd(u, v, units)
    wdir = _calc_wdir(u, v)

    idx_end = -2 if two_d else -3
    
    outdims = list(wspd.shape[0:idx_end]) + [2] + list(wspd.shape[idx_end:])

    result = np.zeros(outdims, wspd.dtype)
    
    idxs0 = ((0,Ellipsis, slice(None), slice(None), slice(None)) 
            if not two_d else 
            (1, Ellipsis, slice(None), slice(None)))
    
    idxs1 = ((1, Ellipsis, slice(None), slice(None), slice(None)) 
            if not two_d else 
            (0, Ellipsis, slice(None), slice(None)))
    
    result[idxs0] = wspd[:]
    result[idxs1] = wdir[:]
    
    return result


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="ua",
                   description="destaggered u-wind component",
                   wind_ncvar=True, 
                   two_d=False, 
                   wspd_wdir=False)
@convert_units("wind", "m s-1")
def get_u_destag(wrfnc, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None,
                 units="m s-1"):
    varname = either("U", "UU")(wrfnc)
    u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = destagger(u_vars[varname], -1)
    
    return u


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="va",
                   description="destaggered v-wind component",
                   two_d=False,
                   wind_ncvar=True, 
                   wspd_wdir=False)
@convert_units("wind", "m s-1")
def get_v_destag(wrfnc, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None,
                 units="m s-1"):
    varname = either("V", "VV")(wrfnc)
    v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = destagger(v_vars[varname], -2)
    
    return v


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="wa",
                   description="destaggered w-wind component",
                   two_d=False,
                   wind_ncvar=True, 
                   wspd_wdir=False)
@convert_units("wind", "m s-1")
def get_w_destag(wrfnc, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None,
                 units="m s-1"):
    w_vars = extract_vars(wrfnc, timeidx, "W", method, squeeze, cache, 
                          meta=False, _key=_key)
    w = destagger(w_vars["W"], -3)
    return w


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="wspd_wdir",
                   description="wspd,wdir in projection space",
                   two_d=False, 
                   wspd_wdir=True)
def get_destag_wspd_wdir(wrfnc, timeidx=0, method="cat", 
                         squeeze=True, cache=None, meta=True, _key=None,
                         units="m s-1"):
    varname = either("U", "UU")(wrfnc)
    u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = destagger(u_vars[varname], -1)
    
    varname = either("V", "VV")(wrfnc)
    v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = destagger(v_vars[varname], -2)
    
    return _calc_wspd_wdir(u, v, False, units)


@set_wind_metadata(copy_varname=either("PSFC", "F"), 
                   name="wspd_wdir10",
                   description="10m wspd,wdir in projection space",
                   two_d=False, 
                   wspd_wdir=True)
def get_destag_wspd_wdir10(wrfnc, timeidx=0, method="cat", 
                           squeeze=True, cache=None, meta=True, _key=None, 
                           units="m s-1"):
    
    varname = either("U10", "UU")(wrfnc)
    u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = (u_vars[varname] if varname == "U10" else 
         destagger(u_vars[varname][...,0,:,:], -1)) 
    
    varname = either("V10", "VV")(wrfnc)
    v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = (v_vars[varname] if varname == "V10" else 
         destagger(v_vars[varname][...,0,:,:], -2))
    
    return _calc_wspd_wdir(u, v, True, units)

