from functools import wraps
from inspect import getargspec
from itertools import product

import numpy as n

from wrf.var.units import do_conversion, check_units

__all__ = ["convert_units", "handle_left_iter"]

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

def handle_left_iter(alg_out_num_dims, ref_var_expected_dims, ref_var_idx=-1,
                   ref_var_name=None, alg_out_fixed_dims=(), 
                   ignore_args=(), ignore_kargs={}):
    """Decorator to handle iterating over leftmost dimensions when using 
    multiple files and/or multiple times.
    
    Arguments:
        - alg_out_num_dims - the return dimension count from the numerical 
        algorithm
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
    
    """
    def indexing_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            if ref_var_idx >= 0:
                ref_var = args[ref_var_idx]
            else:
                ref_var = kargs[ref_var_name]
            
            ref_var_shape = ref_var.shape
            extra_dim_num = len(ref_var_shape) - ref_var_expected_dims
            
            # No special left side iteration, return the function result
            if (extra_dim_num == 0):
                return func(*args, **kargs)
            
            # Start by getting the left-most 'extra' dims
            outdims = [ref_var_shape[x] for x in xrange(extra_dim_num)]
            extra_dims = list(outdims) # Copy the left-most dims for iteration
            
            # Append the right-most dimensions
            # Note: numerical algorithm output dim count might be less than 
            # or greater than reference variable
            if alg_out_num_dims <= ref_var_expected_dims:
                outdims += [ref_var_shape[x] 
                        for x in xrange(-alg_out_num_dims,0,1)]
            else: 
                # numerical algorithm has greater dim count than reference, 
                # Note: user must provide this information to decorator
                right_dims = [dim for dim in alg_out_fixed_dims]
                right_dims += [ref_var_shape[x] 
                        for x in xrange(-alg_out_num_dims,0,1)]
                outdims += right_dims
                
            
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
                        output = n.zeros(outdims, ref_var.dtype)
                        out_inited = True
                    
                    output[left_and_slice_idxs] = res[:]
                    
                else:   # This should be a list or a tuple (cape)
                    if not out_inited:
                        output = [n.zeros(outdims, ref_var.dtype)]
                        out_inited = True
                    
                    for outarr in output:
                        outarr[left_and_slice_idxs] = res[:]
                
            
            return output
        
        return func_wrapper
    
    return indexing_decorator


