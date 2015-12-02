
from wrf.var.constants import Constants
from wrf.var.extension import computetk, computeeth, computetv, computewetbulb
from wrf.var.decorators import convert_units

__all__ = ["get_theta", "get_temp", "get_eth", "get_tv", "get_tw"]

@convert_units("temp", "k")
def get_theta(wrfnc, units="k", timeidx=0):
    t = wrfnc.variables["T"][timeidx,:,:,:]
    full_t = t + Constants.T_BASE

    return full_t

@convert_units("temp", "k")
def get_temp(wrfnc, units="k", timeidx=0):
    """Return the temperature in Kelvin or Celsius"""
    
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    return tk

@convert_units("temp", "k")
def get_eth(wrfnc, units="k", timeidx=0):
    "Return equivalent potential temperature (Theta-e) in Kelvin"
    
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    qv = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    eth = computeeth ( qv, tk, full_p )
    
    return eth
    
@convert_units("temp", "k")
def get_tv(wrfnc, units="k", timeidx=0):
    "Return the virtual temperature (tv) in Kelvin or Celsius"
    
    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    qv = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    tv = computetv(tk,qv)
    
    return tv
    

@convert_units("temp", "k")
def get_tw(wrfnc, units="k", timeidx=0):
    "Return the wetbulb temperature (tw)"

    t = wrfnc.variables["T"][timeidx,:,:,:]
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    qv = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    
    tk = computetk(full_p, full_t)
    tw = computewetbulb(full_p,tk,qv)
    
    return tw
    
    

