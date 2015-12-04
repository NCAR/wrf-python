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

def combine_list_and_times(alg_out_dim):
    def combine_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kargs):
            # Multiple times?
            multitime = False
            if "timeidx" in kargs:
                if kargs["timeidx"] == 0:
                    multitime = True
                    
            # Multiple files?
            multifile = False
            if isinstance(args[0], Iterable) and not isinstance(args[0], str):
                multifile = True
            
            # Single file, single time
            if not multitime and not multifile:
                return func(*args, **kargs)
            
            # Get the dimensions
            if multifile:
                wrffile = args[0][0]
            else:
                wrffile = args[0]
            
            # TODO:  Add PyNIO support
            dims = wrffile.dimensions
            we_size = len(dims["west_east"])
            sn_size = len(dims["south_north"])
            bt_size = len(dims["bottom_top"])
            time_size = len(dims["Time"])
            
            
            if alg_out_dim == 2:
                pass
            elif algout_dim == 3:
                pass
            elif algout_dim == 4:
                pass
            else:
                raise RuntimeError("invalid algorithm output dimsize")
            
            
                
             
            return func(*args, **kargs)
            
        return func_wrapper
    
    return combine_decorator
