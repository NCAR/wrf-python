from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants
from .extension import _tk, _eth, _tv, _wetbulb
from .decorators import convert_units
from .metadecorators import copy_and_set_metadata
from .util import extract_vars


@copy_and_set_metadata(copy_varname="T", name="theta", 
                       description="potential temperature")
@convert_units("temp", "k")
def get_theta(wrfnc, timeidx=0, method="cat", squeeze=True, 
              cache=None, meta=True, _key=None,
              units="K"):
    varnames = ("T",)
    
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    full_t = t + Constants.T_BASE

    return full_t


@copy_and_set_metadata(copy_varname="T", name="temp", 
                       description="temperature")
@convert_units("temp", "k")
def get_temp(wrfnc, timeidx=0, method="cat", squeeze=True, 
             cache=None, meta=True, _key=None,
             units="K"):
    """Return the temperature in Kelvin or Celsius"""
    
    varnames=("T", "P", "PB")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    return tk


@copy_and_set_metadata(copy_varname="T", name="theta_e", 
                       description="equivalent potential temperature")
@convert_units("temp", "K")
def get_eth(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True, _key=None,
            units="K"):
    "Return equivalent potential temperature (Theta-e) in Kelvin"
    
    varnames=("T", "P", "PB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    eth = _eth(qv, tk, full_p)
    
    return eth


@copy_and_set_metadata(copy_varname="T", name="tv", 
                       description="virtual temperature")
@convert_units("temp", "k")
def get_tv(wrfnc, timeidx=0, method="cat", squeeze=True, 
           cache=None, meta=True, _key=None,
           units="K"):
    "Return the virtual temperature (tv) in Kelvin or Celsius"
    
    varnames=("T", "P", "PB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    tv = _tv(tk, qv)
    
    return tv
    
    
@copy_and_set_metadata(copy_varname="T", name="twb", 
                       description="wetbulb temperature")
@convert_units("temp", "k")
def get_tw(wrfnc, timeidx=0, method="cat", squeeze=True, 
           cache=None, meta=True, _key=None,
           units="K"):
    "Return the wetbulb temperature (tw)"
    
    varnames=("T", "P", "PB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    
    tk = _tk(full_p, full_t)
    tw = _wetbulb(full_p, tk, qv)
    
    return tw


def get_tk(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None, 
           meta=True, _key=None):
    return get_temp(wrfnc, timeidx, method, squeeze, cache, meta, _key, 
                    units="K")


def get_tc(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None,
           meta=True, _key=None):
    return get_temp(wrfnc, timeidx, method, squeeze, cache, meta, _key,
                    units="degC")
    
    

