from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

import wrapt 

from .destag import destagger
from .util import iter_left_indexes, py3range, npvalues
from .config import xarray_enabled
from .constants import Constants

if xarray_enabled():
    from xarray import DataArray

__all__ = ["uvmet_left_iter"]


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
        outdims = u.shape[0:extra_dim_num]
        extra_dims = list(outdims) # Copy the left-most dims for iteration
        
        # Append the right-most dimensions
        outdims += [2] # For u/v components
        
        #outdims += [u.shape[x] for x in py3range(-num_right_dims,0,1)]
        outdims += list(u.shape[-num_right_dims:])
        
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
            result = wrapped(new_u, new_v, new_lat, new_lon, cen_long, cone)
            
            # Note:  The 2D version will return a 3D array with a 1 length
            # dimension.  Numpy is unable to broadcast this without 
            # sqeezing first.
            result = np.squeeze(result) 
            
            output[left_and_slice_idxs] = result[:]
            
        return output
    
    return func_wrapper

def uvmet_left_iter_nocopy(alg_dtype=np.float64):
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
        
        orig_dtype = u.dtype
        
        if u.ndim == lat.ndim:
            num_right_dims = 2
            is_3d = False
        else:
            num_right_dims = 3
            is_3d = True
        
        has_missing = False
        u_arr = u
        if isinstance(u, DataArray):
            u_arr = npvalues(u)
            
        v_arr = v
        if isinstance(v, DataArray):
            v_arr = npvalues(v)
          
        umissing = Constants.DEFAULT_FILL  
        if isinstance(u_arr, np.ma.MaskedArray):
            has_missing = True
            umissing = u_arr.fill_value
        
        vmissing = Constants.DEFAULT_FILL 
        if isinstance(v_arr, np.ma.MaskedArray):
            has_missing = True
            vmissing = v_arr.fill_value
            
        uvmetmissing = umissing
           
        is_stag = False
        if (u.shape[-1] != lat.shape[-1] or u.shape[-2] != lat.shape[-2]):
            is_stag = True
            # Sanity check
            if (v.shape[-1] == lat.shape[-1] or v.shape[-2] == lat.shape[-2]):
                raise ValueError("u is staggered but v is not")
        
        if (v.shape[-1] != lat.shape[-1] or v.shape[-2] != lat.shape[-2]):
            is_stag = True
            # Sanity check
            if (u.shape[-1] == lat.shape[-1] or u.shape[-2] == lat.shape[-2]):
                raise ValueError("v is staggered but u is not")
        
        if is_3d:
            extra_dim_num = u.ndim - 3
        else:
            extra_dim_num = u.ndim - 2
            
        
        # No special left side iteration, return the function result
        if (extra_dim_num == 0):
            return wrapped(u, v, lat, lon, cen_long, cone, isstag=is_stag,
                           has_missing=has_missing, umissing=umissing,
                           vmissing=vmissing, uvmetmissing=uvmetmissing)
        
        # Start by getting the left-most 'extra' dims
        outdims = u.shape[0:extra_dim_num]
        extra_dims = list(outdims) # Copy the left-most dims for iteration
        
        # Append the right-most dimensions
        if not is_3d:
            outdims += (2,1) # Fortran routine needs 3 dimensions
        else:
            outdims += (2,)
        
        #outdims += [u.shape[x] for x in py3range(-num_right_dims,0,1)]
        outdims += u.shape[-num_right_dims:]
        
        output = np.empty(outdims, alg_dtype)
        
        for left_idxs in iter_left_indexes(extra_dims):

            left_and_slice_idxs = left_idxs + (slice(None),)
                    
            new_u = u[left_and_slice_idxs]
            new_v = v[left_and_slice_idxs]
            new_lat = lat[left_and_slice_idxs]
            new_lon = lon[left_and_slice_idxs]
            outview = output[left_and_slice_idxs]
            
            # Call the numerical routine
            result = wrapped(new_u, new_v, new_lat, new_lon, cen_long, cone,
                             isstag=is_stag, has_missing=has_missing, 
                             umissing=umissing, vmissing=vmissing, 
                             uvmetmissing=uvmetmissing, outview=outview)
            
            # Make sure the result is the same data as what got passed in 
            # Can delete this once everything works
            if (result.__array_interface__["data"][0] != 
                outview.__array_interface__["data"][0]):
                raise RuntimeError("output array was copied")
            
        
        output = output.astype(orig_dtype)
        
        if has_missing:
            output = np.ma.masked_values(output, uvmetmissing)
        
        return output.squeeze()
    
    return func_wrapper

