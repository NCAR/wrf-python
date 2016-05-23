from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from collections import Iterable

import wrapt 
import numpy as np
import numpy.ma as ma

from .units import do_conversion, check_units
from .util import (iter_left_indexes, viewitems, from_args, npvalues, range2)
from .config import xarray_enabled

if xarray_enabled():
    from xarray import DataArray

__all__ = ["convert_units", "handle_left_iter", "uvmet_left_iter", 
           "handle_casting", "handle_extract_transpose"]
  
def convert_units(unit_type, alg_unit):
    """A decorator that applies unit conversion to a function's output array.
    
    Arguments:
    
        - unit_type - the unit type category (wind, pressure, etc)
        - alg_unit - the units that the function returns by default
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        desired_units = from_args(wrapped, "units", *args, **kwargs)["units"]
        check_units(desired_units, unit_type)
        
        # Unit conversion done here
        return do_conversion(wrapped(*args, **kwargs), unit_type, 
                             alg_unit, desired_units)
        
    return func_wrapper


def _calc_out_dims(outvar, left_dims):
    left_dims = [x for x in left_dims]
    right_dims = [x for x in outvar.shape]
    return left_dims + right_dims 

def handle_left_iter(ref_var_expected_dims, ref_var_idx=-1,
                   ref_var_name=None,
                   ignore_args=None, ignore_kargs=None):
    """Decorator to handle iterating over leftmost dimensions when using 
    multiple files and/or multiple times.
    
    Arguments:
        - ref_var_expected_dims - the number of dimensions that the Fortran 
        algorithm is expecting for the reference variable
        - ref_var_idx - the index in args used as the reference variable for 
        calculating leftmost dimensions
        - ref_var_name - the keyword argument name for kwargs used as the 
        reference varible for calculating leftmost dimensions
        - alg_out_fixed_dims - additional fixed dimension sizes for the 
        numerical algorithm (e.g. uvmet has a fixed left dimsize of 2)
        - ignore_args - indexes of any arguments which should be passed 
        directly without any slicing
        - ignore_kargs - keys of any keyword arguments which should be passed
        directly without slicing
        - ref_stag_dim - in some cases the reference variable dimensions
        have a staggered dimension that needs to be corrected when calculating
        the output dimensions (e.g. avo)
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        _ignore_args = ignore_args if ignore_args is not None else ()
        _ignore_kargs = ignore_kargs if ignore_kargs is not None else ()
        
        if ref_var_idx >= 0:
            ref_var = args[ref_var_idx]
        else:
            ref_var = kwargs[ref_var_name]
        
        ref_var_shape = ref_var.shape
        extra_dim_num = ref_var.ndim - ref_var_expected_dims
        
        # No special left side iteration, return the function result
        if (extra_dim_num == 0):
            return wrapped(*args, **kwargs)
        
        # Start by getting the left-most 'extra' dims
        extra_dims = [ref_var_shape[x] for x in range2(extra_dim_num)]            
        
        out_inited = False
        for left_idxs in iter_left_indexes(extra_dims):
            # Make the left indexes plus a single slice object
            # The single slice will handle all the dimensions to
            # the right (e.g. [1,1,:])
            left_and_slice_idxs = left_idxs + (slice(None), )
            
            # Slice the args if applicable
            new_args = [arg[left_and_slice_idxs] 
                        if i not in _ignore_args else arg 
                        for i,arg in enumerate(args)]
                
            # Slice the kwargs if applicable
            new_kargs = {key:(val[left_and_slice_idxs] 
                         if key not in _ignore_kargs else val)
                         for key,val in viewitems(kwargs)}
            
            # Call the numerical routine
            res = wrapped(*new_args, **new_kargs)
            
            if isinstance(res, np.ndarray):
                # Output array
                if not out_inited:
                    outdims = _calc_out_dims(res, extra_dims)
                    if not isinstance(res, ma.MaskedArray):
                        output = np.empty(outdims, ref_var.dtype)
                        masked = False
                    else:
                        output = ma.MaskedArray(
                                        np.zeros(outdims, ref_var.dtype),
                                        mask=np.zeros(outdims, np.bool_),
                                        fill_value=res.fill_value)
                        masked = True
                    
                    out_inited = True 
                
                if not masked:
                    output[left_and_slice_idxs] = res[:]
                else:
                    output.data[left_and_slice_idxs] = res.data[:]
                    output.mask[left_and_slice_idxs] = res.mask[:]
                
            else:   # This should be a list or a tuple (cape)
                if not out_inited:
                    outdims = _calc_out_dims(res[0], extra_dims)
                    if not isinstance(res[0], ma.MaskedArray):
                        output = [np.empty(outdims, ref_var.dtype) 
                                  for i in range2(len(res))]
                        masked = False
                    else:
                        output = [ma.MaskedArray(
                                    np.zeros(outdims, ref_var.dtype),
                                    mask=np.zeros(outdims, np.bool_),
                                    fill_value=res[0].fill_value) 
                                  for i in range2(len(res))]
                        masked = True
                    
                    out_inited = True
                
                for i,outarr in enumerate(res):
                    if not masked:
                        output[i][left_and_slice_idxs] = outarr[:]
                    else:
                        output[i].data[left_and_slice_idxs] = outarr.data[:]
                        output[i].mask[left_and_slice_idxs] = outarr.mask[:]
            
        
        return output
    
    return func_wrapper


def handle_casting(ref_idx=0, arg_idxs=None, karg_names=None, dtype=np.float64):
    """Decorator to handle casting to/from required dtype used in 
    numerical routine.
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        _arg_idxs = arg_idxs if arg_idxs is not None else ()
        _karg_names = karg_names if karg_names is not None else ()
        
        orig_type = args[ref_idx].dtype
        
        new_args = [arg.astype(dtype)
                    if i in _arg_idxs else arg 
                    for i,arg in enumerate(args)]
        
        new_kargs = {key:(val.astype(dtype)
                    if key in _karg_names else val)
                    for key,val in viewitems(kwargs)}
        
        res = wrapped(*new_args, **new_kargs) 
        
        if isinstance(res, np.ndarray):
            if res.dtype == orig_type:
                return res
            return res.astype(orig_type)
        else:   # got back a sequence of arrays
            return tuple(arr.astype(orig_type) 
                         if arr.dtype != orig_type else arr
                         for arr in res)
    
    return func_wrapper


def _extract_and_transpose(arg, do_transpose):
    if xarray_enabled():
        if isinstance(arg, DataArray):
            arg = npvalues(arg)
    
    if do_transpose:
        if isinstance(arg, np.ndarray):
            if not arg.flags.f_contiguous and arg.ndim > 1:
                return arg.T
        
    return arg

    
def handle_extract_transpose(do_transpose=True):
    """Decorator to extract the data array from a DataArray and also
    transposes the view of the data if the data is not fortran contiguous.
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        new_args = [_extract_and_transpose(arg, do_transpose) for arg in args]
        
        new_kargs = {key:_extract_and_transpose(val) 
                     for key,val in viewitems(kwargs)}
        
        res = wrapped(*new_args, **new_kargs) 
        
        if isinstance(res, np.ndarray):
            if res.flags.f_contiguous and res.ndim > 1:
                return res.T
        elif isinstance(res, Iterable):
            return tuple(x.T if x.flags.f_contiguous and x.ndim > 1 else x 
                         for x in res)
        
        return res
    
    return func_wrapper
    


    
    


