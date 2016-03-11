from functools import wraps
from inspect import getargspec

import numpy as n
import numpy.ma as ma

from wrf.var.units import do_conversion, check_units
from wrf.var.destag import destagger
from wrf.var.util import iter_left_indexes
from wrf.var.config import xarray_enabled

if xarray_enabled():
    from xarray import DataArray



__all__ = ["convert_units", "handle_left_iter", "uvmet_left_iter", 
           "handle_casting" "set_metadata"]

# TODO:  In python 3.x, getargspec is deprecated
class from_args(object):
    def __init__(self, argname):
        self.argname
        
    def __call__(self, func, *args, **kargs):
        """Parses the function args and kargs looking for the desired argument 
        value. Otherwise, the value is taken from the default keyword argument 
        using the arg spec.
        
        """
        # If units are provided to the method call, use them.  
        # Otherwise, need to parse the argspec to find what the default 
        # value is since wraps does not preserve this.
        argspec = getargspec(func)
        
        if self.argname not in argspec.args and self.argname not in kargs:
            return None
            
        try:
            result_idx = argspec.args.index(self.argname)
        except ValueError:
            result_idx = None
            
        if (self.argname in kargs):
            result = kargs[self.argname]
        elif (len(args) > result_idx and result_idx is not None):
            result = args[result_idx]
        else:    
            arg_idx_from_right = (len(argspec.args) - 
                                  argspec.args.index(self.argname))
            result = argspec.defaults[-arg_idx_from_right]
            
        return result
    
def convert_units(unit_type, alg_unit):
    """A decorator that applies unit conversion to a function's output array.
    
    Arguments:
    
        - unit_type - the unit type category (wind, pressure, etc)
        - alg_unit - the units that the function returns by default
    
    """
    def convert_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            
            desired_units = from_args("units")(func, *args, **kargs)
            check_units(desired_units, unit_type)
            
            # Unit conversion done here
            return do_conversion(func(*args, **kargs), unit_type, 
                                 alg_unit, desired_units)
        return func_wrapper
    
    return convert_decorator


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
        - ref_var_name - the keyword argument name for kargs used as the 
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
    def indexing_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            ignore_args = ignore_args if ignore_args is not None else ()
            ignore_kargs = ignore_kargs if ignore_kargs is not None else ()
            
            if ref_var_idx >= 0:
                ref_var = args[ref_var_idx]
            else:
                ref_var = kargs[ref_var_name]
            
            ref_var_shape = ref_var.shape
            extra_dim_num = ref_var.ndim - ref_var_expected_dims
            
            # No special left side iteration, return the function result
            if (extra_dim_num == 0):
                return func(*args, **kargs)
            
            # Start by getting the left-most 'extra' dims
            extra_dims = [ref_var_shape[x] for x in xrange(extra_dim_num)]            
            
            out_inited = False
            for left_idxs in iter_left_indexes(extra_dims):
                # Make the left indexes plus a single slice object
                # The single slice will handle all the dimensions to
                # the right (e.g. [1,1,:])
                left_and_slice_idxs = ([x for x in left_idxs] + 
                                       [slice(None, None, None)])
                # Slice the args if applicable
                new_args = [arg[left_and_slice_idxs] 
                            if i not in ignore_args else arg 
                            for i,arg in enumerate(args)]
                    
                # Slice the kargs if applicable
                new_kargs = {key:(val[left_and_slice_idxs] 
                             if key not in ignore_kargs else val)
                             for key,val in kargs.iteritems()}
                
                # Call the numerical routine
                res = func(*new_args, **new_kargs)
                
                if isinstance(res, n.ndarray):
                    # Output array
                    if not out_inited:
                        outdims = _calc_out_dims(res, extra_dims)
                        if not isinstance(res, ma.MaskedArray):
                            output = n.zeros(outdims, ref_var.dtype)
                            masked = False
                        else:
                            output = ma.MaskedArray(
                                            n.zeros(outdims, ref_var.dtype),
                                            mask=n.zeros(outdims, n.bool_),
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
                            output = [n.zeros(outdims, ref_var.dtype) 
                                      for i in xrange(len(res))]
                            masked = False
                        else:
                            output = [ma.MaskedArray(
                                        n.zeros(outdims, ref_var.dtype),
                                        mask=n.zeros(outdims, n.bool_),
                                        fill_value=res[0].fill_value) 
                                      for i in xrange(len(res))]
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
    
    return indexing_decorator

def uvmet_left_iter():
    """Decorator to handle iterating over leftmost dimensions when using 
    multiple files and/or multiple times with the uvmet product.
    
    """
    def indexing_decorator(func):
        @wraps(func)
        def func_wrapper(*args):
            u = args[0]
            v = args[1]
            lat = args[2]
            lon = args[3]
            cen_long  = args[4]
            cone = args[5]
            
            if u.ndim == lat.ndim:
                num_right_dims = 2
                is_3d = False
            else:
                num_right_dims = 3
                is_3d = True
            
            is_stag = False
            if ((u.shape[-1] != lat.shape[-1]) or 
                (u.shape[-2] != lat.shape[-2])):
                is_stag = True
            
            if is_3d:
                extra_dim_num = u.ndim - 3
            else:
                extra_dim_num = u.ndim - 2
                
            if is_stag:
                u = destagger(u,-1)
                v = destagger(v,-2)
            
            # No special left side iteration, return the function result
            if (extra_dim_num == 0):
                return func(u,v,lat,lon,cen_long,cone)
            
            # Start by getting the left-most 'extra' dims
            outdims = [u.shape[x] for x in xrange(extra_dim_num)]
            extra_dims = list(outdims) # Copy the left-most dims for iteration
            
            # Append the right-most dimensions
            outdims += [2] # For u/v components
            
            outdims += [u.shape[x] for x in xrange(-num_right_dims,0,1)]
            
            output = n.zeros(outdims, u.dtype)
            
            for left_idxs in iter_left_indexes(extra_dims):
                # Make the left indexes plus a single slice object
                # The single slice will handle all the dimensions to
                # the right (e.g. [1,1,:])
                left_and_slice_idxs = ([x for x in left_idxs] + 
                                       [slice(None, None, None)])
                        
                
                new_u = u[left_and_slice_idxs]
                new_v = v[left_and_slice_idxs]
                new_lat = lat[left_and_slice_idxs]
                new_lon = lon[left_and_slice_idxs]
                
                # Call the numerical routine
                res = func(new_u, new_v, new_lat, new_lon, cen_long, cone)
                
                # Note:  The 2D version will return a 3D array with a 1 length
                # dimension.  Numpy is unable to broadcast this without 
                # sqeezing first.
                res = n.squeeze(res) 
                
                output[left_and_slice_idxs] = res[:]
                
            
            return output
        
        return func_wrapper
    
    return indexing_decorator

def handle_casting(ref_idx=0, arg_idxs=None, karg_names=None, dtype=n.float64):
    """Decorator to handle casting to/from required dtype used in 
    numerical routine.
    
    """
    def cast_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            arg_idxs = arg_idxs if arg_idxs is not None else ()
            karg_names = karg_names if karg_names is not None else ()
            
            orig_type = args[ref_idx].dtype
            
            new_args = [arg.astype(dtype)
                        if i in arg_idxs else arg 
                        for i,arg in enumerate(args)]
            
            new_kargs = {key:(val.astype(dtype)
                        if key in karg_names else val)
                        for key,val in kargs.iteritems()}
            
            res = func(*new_args, **new_kargs) 
            
            if isinstance(res, n.ndarray):
                if res.dtype == orig_type:
                    return res
                return res.astype(orig_type)
            else:   # got back a sequence of arrays
                return tuple(arr.astype(orig_type) 
                             if arr.dtype != orig_type else arr
                             for arr in res)
        
        return func_wrapper
    
    return cast_decorator

def _extract_and_transpose(arg):
    if xarray_enabled():
        if isinstance(arg, DataArray):
            arg = arg.values
    
    if isinstance(arg, n.ndarray):
        if not arg.flags.F_CONTIGUOUS:
            return arg.T
        
    return arg
    
def handle_extract_transpose():
    """Decorator to extract the data array from a DataArray and also
    transposes if the data is not fortran contiguous.
    
    """
    def extract_transpose_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            
            new_args = [_extract_and_transpose(arg) for arg in args]
            
            new_kargs = {key:_extract_and_transpose(val) 
                         for key,val in kargs.iteritems()}
                
            res = func(*new_args, **new_kargs) 
            
            if isinstance(res, n.ndarray):
                if res.flags.F_CONTIGUOUS:
                    return res.T
            else:
                return tuple(x.T if x.flags.F_CONTIGUOUS else x for x in res)
            
            return res
        
        return func_wrapper
    
    return extract_transpose_decorator

def set_metadata(copy_from=None, ignore=None, **fixed_kargs):
    """Decorator to set the attributes for a WRF method.
    
    """
    def attr_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            if not xarray_enabled():
                return func(*args, **kargs)
            
            result = func(*args, **kargs)
            
            units = from_args("units")(func, *args, **kargs)
            if units is not None:
                result.attrs["units"] = units
            
            for argname, val in fixed_kargs.iteritems():
                result[argname] = val
            
            return result
        
        return func_wrapper
    
    return attr_decorator
    
    


