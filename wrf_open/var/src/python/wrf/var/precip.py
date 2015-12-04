

from wrf.var.util import extract_vars

__all__ = ["get_accum_precip", "get_precip_diff"]

def get_accum_precip(wrfnc, timeidx=0):
    ncvars = extract_vars(wrfnc, timeidx, vars=("RAINC", "RAINNC"))
    rainc = ncvars["RAINC"]
    rainnc = ncvars["RAINNC"]
    
    rainsum = rainc + rainnc
    
    return rainsum

def get_precip_diff(wrfnc1, wrfnc2, timeidx=0):
    vars1 = extract_vars(wrfnc1, timeidx, vars=("RAINC", "RAINNC"))
    vars2 = extract_vars(wrfnc2, timeidx, vars=("RAINC", "RAINNC"))
    rainc1 = vars1["RAINC"]
    rainnc1 = vars1["RAINNC"]
    
    rainc2 = vars2["RAINC"]
    rainnc2 = vars2["RAINNC"]
    
    rainsum1 = rainc1 + rainnc1
    rainsum2 = rainc2 + rainnc2
    
    return (rainsum1 - rainsum2)

# TODO:  Handle bucket flipping
