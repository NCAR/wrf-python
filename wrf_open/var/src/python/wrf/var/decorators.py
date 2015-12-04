from functools import wraps
from inspect import getargspec
from collections import Iterable

from wrf.var.units import do_conversion, check_units

__all__ = ["convert_units"]

def convert_units(unit_type, alg_unit):
    def convert_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            # If units are provided to the method call, use them.  
            # Otherwise, need to parse the argspec to find what the default 
            # value is.
            if ("units" in kargs):
                desired_units = kargs["units"]
            else:
                argspec = getargspec(func)
                print argspec
                arg_idx_from_right = len(argspec.args) - argspec.args.index("units")
                desired_units = argspec.defaults[-arg_idx_from_right]
                
                #print desired_idx
                #desired_units = argspec.defaults[desired_idx]
                print desired_units
                
            check_units(desired_units, unit_type)
            
            # Unit conversion done here
            return do_conversion(func(*args, **kargs), unit_type, 
                                 alg_unit, desired_units)
        return func_wrapper
    
    return convert_decorator


