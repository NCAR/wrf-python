from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

import wrapt 

from .destag import destagger
from .util import iter_left_indexes, range2

__all__ = ["uvmet_left_iter"]

# Placed in separate module to resolve a circular dependency with destagger 
# module

def uvmet_left_iter():
    """Decorator to handle iterating over leftmost dimensions when using 
    multiple files and/or multiple times with the uvmet product.
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
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
            return wrapped(u, v, lat, lon, cen_long, cone)
        
        # Start by getting the left-most 'extra' dims
        outdims = [u.shape[x] for x in range2(extra_dim_num)]
        extra_dims = list(outdims) # Copy the left-most dims for iteration
        
        # Append the right-most dimensions
        outdims += [2] # For u/v components
        
        outdims += [u.shape[x] for x in range2(-num_right_dims,0,1)]
        
        output = np.empty(outdims, u.dtype)
        
        for left_idxs in iter_left_indexes(extra_dims):
            # Make the left indexes plus a single slice object
            # The single slice will handle all the dimensions to
            # the right (e.g. [1,1,:])
            left_and_slice_idxs = tuple([x for x in left_idxs] + [slice(None)])
                    
            new_u = u[left_and_slice_idxs]
            new_v = v[left_and_slice_idxs]
            new_lat = lat[left_and_slice_idxs]
            new_lon = lon[left_and_slice_idxs]
            
            # Call the numerical routine
            res = wrapped(new_u, new_v, new_lat, new_lon, cen_long, cone)
            
            # Note:  The 2D version will return a 3D array with a 1 length
            # dimension.  Numpy is unable to broadcast this without 
            # sqeezing first.
            res = np.squeeze(res) 
            
            output[left_and_slice_idxs] = res[:]
            
        return output
    
    return func_wrapper

