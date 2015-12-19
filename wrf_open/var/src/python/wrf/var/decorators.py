from functools import wraps
from inspect import getargspec
from itertools import product

import numpy as n

from wrf.var.units import do_conversion, check_units
from wrf.var.destag import destagger

__all__ = ["convert_units", "handle_left_iter", "uvmet_left_iter", 
           "handle_casting"]

def convert_units(unit_type, alg_unit):
    """A decorator that applies unit conversion to a function's output array.
    
    Arguments:
    
        - unit_type - the unit type category (wind, pressure, etc)
        - alg_unit - the units that the function returns by default
    
    """
    def convert_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            # If units are provided to the method call, use them.  
            # Otherwise, need to parse the argspec to find what the default 
            # value is since wraps does not preserve this.
            if ("units" in kargs):
                desired_units = kargs["units"]
            else:
                argspec = getargspec(func)
                arg_idx_from_right = (len(argspec.args) 
                                      - argspec.args.index("units"))
                desired_units = argspec.defaults[-arg_idx_from_right]
                
            check_units(desired_units, unit_type)
            
            # Unit conversion done here
            return do_conversion(func(*args, **kargs), unit_type, 
                                 alg_unit, desired_units)
        return func_wrapper
    
    return convert_decorator

def _left_indexes(dims):
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

def _calc_out_dims(outvar, left_dims):
    left_dims = [x for x in left_dims]
    right_dims = [x for x in outvar.shape]
    return left_dims + right_dims 

def handle_left_iter(ref_var_expected_dims, ref_var_idx=-1,
                   ref_var_name=None,
                   ignore_args=(), ignore_kargs={}):
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
            for left_idxs in _left_indexes(extra_dims):
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
                        output = n.zeros(outdims, ref_var.dtype)
                        out_inited = True
                    
                    output[left_and_slice_idxs] = res[:]
                    
                else:   # This should be a list or a tuple (cape)
                    if not out_inited:
                        outdims = _calc_out_dims(res[0], extra_dims)
                        output = [n.zeros(outdims, ref_var.dtype) 
                                  for i in xrange(len(res))]
                        out_inited = True
                    
                    for i,outarr in enumerate(res):
                        (output[i])[left_and_slice_idxs] = outarr[:]
                
            
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
            
            for left_idxs in _left_indexes(extra_dims):
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

def handle_casting(ref_idx=0, arg_idxs=(), karg_names=None,dtype=n.float64):
    """Decorator to handle casting to/from required dtype used in 
    numerical routine.
    
    """
    def cast_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
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


