from collections import Iterable

import numpy as n

__all__ = ["extract_vars", "extract_global_attrs", "extract_dim",
           "combine_files", "is_standard_wrf_var"]

def _is_multi_time(timeidx):
    if timeidx == -1:
        return True
    return False

def _is_multi_file(wrfnc):
    if isinstance(wrfnc, Iterable) and not isinstance(wrfnc, str):
        return True
    return False

def _get_attr(wrfnc, attr):
    val = getattr(wrfnc, attr)
    
    # PyNIO puts single values in to an array
    if isinstance(val,n.ndarray):
        if len(val) == 1:
            return val[0] 
    return val
        
def extract_global_attrs(wrfnc, attrs):
    if isinstance(attrs, str):
        attrlist = [attrs]
    else:
        attrlist = attrs
        
    multifile = _is_multi_file(wrfnc)
    
    if multifile:
        wrfnc = wrfnc[0]
    
    return {attr:_get_attr(wrfnc, attr) for attr in attrlist}

def extract_dim(wrfnc, dim):
    d = wrfnc.dimensions[dim]
    if not isinstance(d, int):
        return len(d) #netCDF4
    return d # PyNIO
        
def combine_files(wrflist, var, timeidx):
    
    multitime = _is_multi_time(timeidx)
    numfiles = len(wrflist)  
    
    if not multitime:
        time_idx_or_slice = timeidx
    else:
        time_idx_or_slice = slice(None, None, None)
        
    first_var = wrflist[0].variables[var][time_idx_or_slice, :]
    outdims = [numfiles]
    outdims += first_var.shape
    outarr = n.zeros(outdims, first_var.dtype)
    outarr[0,:] = first_var[:]
    
    for idx, wrfnc in enumerate(wrflist[1:], start=1):
        outarr[idx,:] = wrfnc.variables[var][time_idx_or_slice, :]
            
    return outarr

def extract_vars(wrfnc, timeidx, vars):
    if isinstance(vars, str):
        varlist = [vars]
    else:
        varlist = vars
        
    multitime = _is_multi_time(timeidx)
    multifile = _is_multi_file(wrfnc)
    
    # Single file, single time
    if not multitime and not multifile:
        return {var:wrfnc.variables[var][timeidx,:] for var in varlist}
    elif multitime and not multifile:
        return {var:wrfnc.variables[var][:] for var in varlist}
    else:
        return {var:combine_files(wrfnc, var, timeidx) for var in varlist}
    
def is_standard_wrf_var(wrfnc, var):
    multifile = _is_multi_file(wrfnc)
    if multifile:
        wrfnc = wrfnc[0]
    return var in wrfnc.variables
        
    
    
    