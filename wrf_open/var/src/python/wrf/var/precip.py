
import numpy as n

__all__ = ["get_accum_precip", "get_precip_diff"]

def get_accum_precip(wrfnc, timeidx=0):
    rainc = wrfnc.variables["RAINC"][timeidx,:,:]
    rainnc = wrfnc.variables["RAINNC"][timeidx,:,:]
    
    rainsum = rainc + rainnc
    
    return rainsum

def get_precip_diff(wrfnc1, wrfnc2, timeidx=0):
    rainc1 = wrfnc1.variables["RAINC"][timeidx,:,:]
    rainnc1 = wrfnc1.variables["RAINNC"][timeidx,:,:]
    
    rainc2 = wrfnc2.variables["RAINC"][timeidx,:,:]
    rainnc2 = wrfnc2.variables["RAINNC"][timeidx,:,:]
    
    rainsum1 = rainc1 + rainnc1
    rainsum2 = rainc2 + rainnc2
    
    return (rainsum1 - rainsum2)

# TODO:  Handle bucket flipping
