from collections import Iterable
from itertools import product
import datetime as dt

import numpy as n

__all__ = ["extract_vars", "extract_global_attrs", "extract_dim",
           "combine_files", "is_standard_wrf_var", "extract_times",
           "iter_left_indexes", "get_left_indexes", "get_right_slices",
           "is_staggered"]

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
    if _is_multi_file(wrfnc):
        wrfnc = wrfnc[0]
    
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
        
    first_var = n.squeeze(wrflist[0].variables[var][time_idx_or_slice, :])
    outdims = [numfiles]
    outdims += first_var.shape
    outarr = n.zeros(outdims, first_var.dtype)
    outarr[0,:] = first_var[:]
    
    for idx, wrfnc in enumerate(wrflist[1:], start=1):
        outarr[idx,:] = n.squeeze(wrfnc.variables[var][time_idx_or_slice, :])
            
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
        return {var:n.squeeze(wrfnc.variables[var][:]) for var in varlist}
    else:
        return {var:combine_files(wrfnc, var, timeidx) for var in varlist}

def _make_time(timearr):
    return dt.datetime.strptime("".join(timearr[:]), "%Y-%m-%d_%H:%M:%S")

def _file_times(wrfnc,timeidx):
    multitime = _is_multi_time(timeidx)
    if multitime:
        times = wrfnc.variables["Times"][:,:]
        for i in xrange(times.shape[0]):
            yield _make_time(times[i,:])
    else:
        times = wrfnc.variables["Times"][timeidx,:]
        yield _make_time(times)
        

def extract_times(wrfnc,timeidx):
    multi_file = _is_multi_file(wrfnc)
    if not multi_file:
        wrf_list = [wrfnc]
    else:
        wrf_list = wrfnc
    
    return [file_time 
            for wrf_file in wrf_list 
            for file_time in _file_times(wrf_file,timeidx) ]        
        
    
def is_standard_wrf_var(wrfnc, var):
    multifile = _is_multi_file(wrfnc)
    if multifile:
        wrfnc = wrfnc[0]
    return var in wrfnc.variables

def is_staggered(var, wrfnc):
    we = extract_dim(wrfnc, "west_east")
    sn = extract_dim(wrfnc, "south_north")
    bt = extract_dim(wrfnc, "bottom_top")
    
    if (var.shape[-1] != we or var.shape[-2] != sn or var.shape[-3] != bt):
        return True
    
    return False

def get_left_indexes(ref_var, expected_dims):
    """Returns the extra left side dimensions for a variable with an expected
    shape.
    
    For example, if a 2D variable contains an additional left side dimension
    for time, this will return the time dimension size.
    
    """
    extra_dim_num = ref_var.ndim - expected_dims
    
    if (extra_dim_num == 0):
        return []
    
    return [ref_var.shape[x] for x in xrange(extra_dim_num)] 

def iter_left_indexes(dims):
    """A generator which yields the iteration tuples for a sequence of 
    dimensions sizes.
    
    For example, if an array shape is (3,3), then this will yield:
    
    (0,0), (0,1), (1,0), (1,1)
    
    Arguments:
    
        - dims - a sequence of dimensions sizes (e.g. ndarry.shape)
    
    """
    arg = [xrange(dim) for dim in dims]
    for idxs in product(*arg):
        yield idxs
        
def get_right_slices(var, right_ndims, fixed_val=0):
    extra_dim_num = var.ndim - right_ndims
    if extra_dim_num == 0:
        return [slice(None,None,None)] * right_ndims
    
    return [fixed_val]*extra_dim_num + [slice(None,None,None)]*right_ndims
    
    
    
    




        
    
    
    