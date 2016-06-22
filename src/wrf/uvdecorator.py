from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

import wrapt 

#from .destag import destagger
from .util import iter_left_indexes, py3range, npvalues
from .config import xarray_enabled
from .constants import Constants

if xarray_enabled():
    from xarray import DataArray


# def uvmet_left_iter():
#     """Decorator to handle iterating over leftmost dimensions when using 
#     multiple files and/or multiple times with the uvmet product.
#     
#     """
#     @wrapt.decorator
#     def func_wrapper(wrapped, instance, args, kwargs):
#         u = args[0]
#         v = args[1]
#         lat = args[2]
#         lon = args[3]
#         cen_long  = args[4]
#         cone = args[5]
#         
#         if u.ndim == lat.ndim:
#             num_right_dims = 2
#             is_3d = False
#         else:
#             num_right_dims = 3
#             is_3d = True
#         
#         is_stag = False
#         if ((u.shape[-1] != lat.shape[-1]) or 
#             (u.shape[-2] != lat.shape[-2])):
#             is_stag = True
#         
#         if is_3d:
#             extra_dim_num = u.ndim - 3
#         else:
#             extra_dim_num = u.ndim - 2
#             
#         if is_stag:
#             u = destagger(u,-1)
#             v = destagger(v,-2)
#         
#         # No special left side iteration, return the function result
#         if (extra_dim_num == 0):
#             return wrapped(u, v, lat, lon, cen_long, cone)
#         
#         # Start by getting the left-most 'extra' dims
#         outdims = u.shape[0:extra_dim_num]
#         extra_dims = list(outdims) # Copy the left-most dims for iteration
#         
#         # Append the right-most dimensions
#         outdims += [2] # For u/v components
#         
#         #outdims += [u.shape[x] for x in py3range(-num_right_dims,0,1)]
#         outdims += list(u.shape[-num_right_dims:])
#         
#         output = np.empty(outdims, u.dtype)
#         
#         for left_idxs in iter_left_indexes(extra_dims):
#             # Make the left indexes plus a single slice object
#             # The single slice will handle all the dimensions to
#             # the right (e.g. [1,1,:])
#             left_and_slice_idxs = tuple([x for x in left_idxs] + [slice(None)])
#                     
#             new_u = u[left_and_slice_idxs]
#             new_v = v[left_and_slice_idxs]
#             new_lat = lat[left_and_slice_idxs]
#             new_lon = lon[left_and_slice_idxs]
#             
#             # Call the numerical routine
#             result = wrapped(new_u, new_v, new_lat, new_lon, cen_long, cone)
#             
#             # Note:  The 2D version will return a 3D array with a 1 length
#             # dimension.  Numpy is unable to broadcast this without 
#             # sqeezing first.
#             result = np.squeeze(result) 
#             
#             output[left_and_slice_idxs] = result[:]
#             
#         return output
#     
#     return func_wrapper

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
        
        lat_lon_fixed = False
        if lat.ndim == 2:
            lat_lon_fixed = True
            
        if lon.ndim == 2 and not lat_lon_fixed:
            raise ValueError("'lat' and 'lon' shape mismatch")
        
        num_left_dims_u = u.ndim - 2
        num_left_dims_lat = lat.ndim - 2
        
        if (num_left_dims_lat > num_left_dims_u):
            raise ValueError("number of 'lat' dimensions is greater than 'u'")
        
        if lat_lon_fixed:
            mode = 0 # fixed lat/lon
        else:
            if num_left_dims_u == num_left_dims_lat:
                mode = 1 # lat/lon same as u
            else:
                mode = 2 # probably 3D with 2D lat/lon plus Time
        
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
        
        is_stag = 0
        if (u.shape[-1] != lat.shape[-1] or u.shape[-2] != lat.shape[-2]):
            is_stag = 1
            # Sanity check
            if (v.shape[-1] == lat.shape[-1] or v.shape[-2] == lat.shape[-2]):
                raise ValueError("u is staggered but v is not")
        
        if (v.shape[-1] != lat.shape[-1] or v.shape[-2] != lat.shape[-2]):
            is_stag = 1
            # Sanity check
            if (u.shape[-1] == lat.shape[-1] or u.shape[-2] == lat.shape[-2]):
                raise ValueError("v is staggered but u is not")
        
        
        
        # No special left side iteration, return the function result
        if (num_left_dims_u == 0):
            return wrapped(u, v, lat, lon, cen_long, cone, isstag=is_stag,
                           has_missing=has_missing, umissing=umissing,
                           vmissing=vmissing, uvmetmissing=uvmetmissing)

        # Initial output is time,nz,2,ny,nx to create contiguous views
        outdims = u.shape[0:num_left_dims_u]
        extra_dims = tuple(outdims) # Copy the left-most dims for iteration
        
        outdims += (2,)
        
        outdims += lat.shape[-2:]
        
        outview_array = np.empty(outdims, alg_dtype)
        
        for left_idxs in iter_left_indexes(extra_dims):
            left_and_slice_idxs = left_idxs + (slice(None),)
            
            if mode == 0:
                lat_left_and_slice = (slice(None),)
            elif mode == 1:
                lat_left_and_slice = left_and_slice_idxs
            elif mode == 2:
                # Only need the left-most
                lat_left_and_slice = tuple(left_idx 
                            for left_idx in left_idxs[0:num_left_dims_lat])
            
            
            new_u = u[left_and_slice_idxs]
            new_v = v[left_and_slice_idxs]
            new_lat = lat[lat_left_and_slice]
            new_lon = lon[lat_left_and_slice]
            outview = outview_array[left_and_slice_idxs]
            
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
            
        # Need to reshape this so that u_v is left dim, then time (or others), 
        # then nz, ny, nz
        output_dims = (2,)
        output_dims += extra_dims
        output_dims += lat.shape[-2:]
        output = np.empty(output_dims, orig_dtype)
        
        output[0,:] = outview_array[...,0,:,:].astype(orig_dtype)
        output[1,:] = outview_array[...,1,:,:].astype(orig_dtype)
        
        if has_missing:
            output = np.ma.masked_values(output, uvmetmissing)
        
        return output
    
    return func_wrapper

