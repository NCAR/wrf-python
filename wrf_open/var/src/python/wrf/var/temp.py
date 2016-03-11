
from wrf.var.constants import Constants
from wrf.var.extension import computetk, computeeth, computetv, computewetbulb
from wrf.var.decorators import convert_units, set_metadata
from wrf.var.util import extract_vars

__all__ = ["get_theta", "get_temp", "get_eth", "get_tv", "get_tw",
           "get_tk", "get_tc"]

@set_metadata()
@convert_units("temp", "k")
def get_theta(wrfnc, timeidx=0, units="k"):
    ncvars = extract_vars(wrfnc, timeidx, varnames="T")
    t = ncvars["T"]
    full_t = t + Constants.T_BASE

    return full_t

@convert_units("temp", "k")
def get_temp(wrfnc, timeidx=0, units="k"):
    """Return the temperature in Kelvin or Celsius"""
    
    ncvars = extract_vars(wrfnc, timeidx, varnames=("T", "P", "PB"))
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    return tk

@convert_units("temp", "k")
def get_eth(wrfnc, timeidx=0, units="k"):
    "Return equivalent potential temperature (Theta-e) in Kelvin"
    
    ncvars = extract_vars(wrfnc, timeidx, varnames=("T", "P", "PB", "QVAPOR"))
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    eth = computeeth(qv, tk, full_p)
    
    return eth
    
@convert_units("temp", "k")
def get_tv(wrfnc, timeidx=0, units="k"):
    "Return the virtual temperature (tv) in Kelvin or Celsius"
    
    ncvars = extract_vars(wrfnc, timeidx, varnames=("T", "P", "PB", "QVAPOR"))
    
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    tv = computetv(tk,qv)
    
    return tv
    

@convert_units("temp", "k")
def get_tw(wrfnc, timeidx=0, units="k"):
    "Return the wetbulb temperature (tw)"
    
    ncvars = extract_vars(wrfnc, timeidx, varnames=("T", "P", "PB", "QVAPOR"))
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    
    tk = computetk(full_p, full_t)
    tw = computewetbulb(full_p,tk,qv)
    
    return tw

def get_tk(wrfnc, timeidx=0):
    return get_temp(wrfnc, timeidx, units="k")

def get_tc(wrfnc, timeidx=0):
    return get_temp(wrfnc, timeidx, units="c")
    
    

