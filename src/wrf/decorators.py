from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from collections import Iterable, OrderedDict

import wrapt 
import numpy as np
import numpy.ma as ma

from .units import do_conversion, check_units
from .util import iter_left_indexes, from_args, npvalues, combine_dims
from .py3compat import viewitems, viewvalues, isstr
from .config import xarray_enabled
from .constants import Constants

if xarray_enabled():
    from xarray import DataArray
  
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
        

def left_iteration(ref_var_expected_dims,
                   ref_var_right_ndims,
                   insert_dims=None,
                   ref_var_idx=None,
                   ref_var_name=None,
                   ignore_args=None, 
                   ignore_kargs=None,
                   outviews="outview",
                   alg_dtype=np.float64,
                   cast_output=True):
    """Decorator to handle iterating over leftmost dimensions when using 
    multiple files and/or multiple times.
    
    Arguments:
        - ref_var_expected_dims - the number of dimensions that the Fortran 
        algorithm is expecting for the reference variable
        - ref_var_right_ndims - the number of dims from the right to keep for
        the reference variable when making the output.  Can also be a 
        combine_dims instance if the sizes are determined from multiple 
        variables.
        - insert_dims - a sequence of dimensions to insert between the left 
        dimenions (e.g. time) and the kept right dimensions
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
        - outviews - a single key or sequence of keys indicating the keyword
        argument used as the output variable(s)
        - algtype - the data type used in the numerical routine
        - cast_output - cast the final output to the ref_var data type
        
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        _ignore_args = ignore_args if ignore_args is not None else ()
        _ignore_kargs = ignore_kargs if ignore_kargs is not None else ()
        _outkeys = [outviews] if isstr(outviews) else outviews
        
        if ref_var_idx is not None:
            ref_var = args[ref_var_idx]
        else:
            ref_var = kwargs[ref_var_name]
        
        ref_var_dtype = ref_var.dtype
        ref_var_shape = ref_var.shape
        extra_dim_num = ref_var.ndim - ref_var_expected_dims
        
        # No special left side iteration, return the function result
        if (extra_dim_num == 0):
            return wrapped(*args, **kwargs)
        
        # Start by getting the left-most 'extra' dims
        extra_dims = ref_var_shape[0:extra_dim_num]         
        
        mid_dims = () if insert_dims is None else tuple(insert_dims)
        
        if not isinstance(ref_var_right_ndims, combine_dims):
            right_dims = ref_var_shape[-ref_var_right_ndims:]
        else:
            right_dims = ref_var_right_ndims(*args)
        
        left_dims = extra_dims

        outdims = left_dims + mid_dims + right_dims
        
        if "outview" not in kwargs:
            outd = OrderedDict((outkey, np.empty(outdims, alg_dtype))
                                   for outkey in _outkeys)
        
        mask_output = False
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
            
            # Skip the possible empty/missing arrays for the join method
            skip_missing = False
            for arg in new_args:
                if isinstance(arg, DataArray):
                    arr = npvalues(arg)
                elif isinstance(arg, np.ndarray):
                    arr = arg
                else:
                    continue
                
                if isinstance(arr, np.ma.MaskedArray):
                    if arr.mask.all():
                        
                        for output in viewvalues(outd):
                            output[left_and_slice_idxs] = (
                                            Constants.DEFAULT_FILL)
                        skip_missing = True
                        mask_output = True
            
            if skip_missing:
                continue
            
            
            # Insert the output views if one hasn't been provided
            if "outview" not in new_kargs:
                for outkey,output in viewitems(outd):
                    outview = output[left_and_slice_idxs]
                    new_kargs[outkey] = outview
                    
            result = wrapped(*new_args, **new_kargs)
            
            # Make sure the result is the same data as what got passed in 
            # Can delete this once everything works
            if (result.__array_interface__["data"][0] != 
                outview.__array_interface__["data"][0]):
                raise RuntimeError("output array was copied")
        
        
        if len(outd) == 1:
            output = next(iter(viewvalues(outd)))
        else:
            output = tuple(arr for arr in viewvalues(outd))
        
        
        
        if cast_output:
            if isinstance(output, np.ndarray):
                output = output.astype(ref_var_dtype)
            else:
                output = tuple(arr.astype(ref_var_dtype) for arr in output)
        
        # Mostly when used with join
        if mask_output:
            if isinstance(output, np.ndarray):
                output = ma.masked_values(output, Constants.DEFAULT_FILL)
            else:
                output = tuple(ma.masked_values(arr, Constants.DEFAULT_FILL) 
                               for arr in output)
        
        return output
    
    return func_wrapper


def cast_type(ref_idx=0, arg_idxs=None, karg_names=None, 
                   alg_dtype=np.float64,
                   outviews="outview"):
    """Decorator to handle casting to/from required dtype used in 
    numerical routine.
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        _arg_idxs = arg_idxs if arg_idxs is not None else ()
        _karg_names = karg_names if karg_names is not None else ()
        
        # Handle output views if applicable
        _outkeys = [outviews] if isstr(outviews) else outviews
        _outviews = from_args(wrapped, _outkeys, *args, **kwargs)
        
        has_outview = False
        for outkey in _outkeys:
            _outview = _outviews[outkey]
            if _outview is not None:
                has_outview = True
            
        
        orig_type = args[ref_idx].dtype
        
        new_args = [arg.astype(alg_dtype)
                    if i in _arg_idxs else arg 
                    for i,arg in enumerate(args)]
        
        new_kargs = {key:(val.astype(alg_dtype)
                    if key in _karg_names else val)
                    for key,val in viewitems(kwargs)}
        
        result = wrapped(*new_args, **new_kargs) 
        
        # Do nothing for supplied output views
        if not has_outview:
            if isinstance(result, np.ndarray):
                if result.dtype == orig_type:
                    return result
                return result.astype(orig_type)
            elif isinstance(result, Iterable):   # got back a sequence of arrays
                return tuple(arr.astype(orig_type) 
                             if arr.dtype != orig_type else arr
                             for arr in result)
        
        return result
    
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

    
def extract_and_transpose(do_transpose=True, outviews="outview"):
    """Decorator to extract the data array from a DataArray and also
    transposes the view of the data if the data is not fortran contiguous.
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        # Handle output views if applicable
        _outkeys = [outviews] if isstr(outviews) else outviews
        _outviews = from_args(wrapped, _outkeys, *args, **kwargs)
        
        has_outview = False
        for outkey in _outkeys:
            _outview = _outviews[outkey]
            if _outview is not None:
                has_outview = True
        
        new_args = [_extract_and_transpose(arg, do_transpose) for arg in args]
        
        new_kargs = {key:_extract_and_transpose(val, do_transpose) 
                     for key,val in viewitems(kwargs)}
        
        result = wrapped(*new_args, **new_kargs) 
        
        # Do nothing for supplied output views
        if has_outview:
            return result
        
        if isinstance(result, np.ndarray):
            if result.flags.f_contiguous and result.ndim > 1:
                return result.T
        elif isinstance(result, Iterable):
            return tuple(x.T if x.flags.f_contiguous and x.ndim > 1 else x 
                         for x in result)
        
        return result
    
    return func_wrapper


def check_args(refvaridx, refvarndim, rightdims, stagger=None,
               refstagdim=None):
    """
    
    refvaridx - reference variable to check for presence of left dimensions
    refvarndim - the number of dimensions for the references variable 
                 expected by the fortran routine
    rightdims - the expected number of right dimensions for each argument
    stagger - the dimension which is staggered for each argument.  Use
              None to indicate no staggering.
    refstagdim - If the reference variable is staggered, indicate the dim that
                 is staggered
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        refvar = args[refvaridx]
        try:
            _ndim = refvar.ndim
        except AttributeError:
            raise ValueError("argument {} is not an arraylike "
                             "object".format(refvaridx))
        else:
            extra_dims = refvar.ndim - refvarndim
        
        # Always use unstaggered as the basis of comparison
        if refstagdim is not None:
            _refshape = list(refvar.shape)
            _refshape[refstagdim] -= 1
            _refshape = tuple(_refshape)
        else:
            _refshape = refvar.shape
        
        if stagger is None:
            _stagger = [None]*len(rightdims)
        else:
            _stagger = stagger
        
        for i,ndim in enumerate(rightdims):
            if ndim is None:
                continue
            
            var = args[i]
            
            try:
                _ = var.ndim
            except AttributeError:
                raise ValueError("argument {} is not an arraylike "
                                 "object".format(i))
            
            right_var_ndims = rightdims[i]
            
            # Check that the number of dims is correct
            if (var.ndim - extra_dims != right_var_ndims):
                raise ValueError("invalid number of dimensions for argument "
                                 "{} (got {}, expected {}).".format(i, 
                                                var.ndim,
                                                right_var_ndims + extra_dims))
            
            # Add 1 to the reference staggered dim index before doing the check
            if _stagger[i] is not None:
                ref_shape = list(_refshape)
                ref_shape[_stagger[i]] += 1
                ref_shape = tuple(ref_shape)
            else:
                ref_shape = _refshape
                
            ref_right_sizes = ref_shape[extra_dims:]
            
            # Check that right dimensions are lined up
            if (var.shape[-right_var_ndims:] != 
                ref_right_sizes[-right_var_ndims:]):
                
                raise ValueError("invalid shape for argument "
                                 "{} (got {}, expected {})".format(i, 
                                        var.shape[-right_var_ndims:],
                                        ref_right_sizes[-right_var_ndims:]))
        
        return wrapped(*args, **kwargs)
    
    return func_wrapper


    


    
    


