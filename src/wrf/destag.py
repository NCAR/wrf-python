from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .decorators import extract_and_transpose
from .metadecorators import set_destag_metadata


@set_destag_metadata()
@extract_and_transpose(do_transpose=False)
def destagger(var, stagger_dim, meta=False):
    """Return the variable on the unstaggered grid.
    
    This function destaggers the variable by taking the average of the 
    values located on either side of the grid box. 
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A variable 
            on a staggered grid.
        
        stagger_dim (:obj:`int`): The dimension index to destagger.
            Negative values can be used to choose dimensions referenced 
            from the right hand side (-1 is the rightmost dimension).
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is False.
    
    Returns:
         
        :class:`xarray.DataArray` or :class:`numpy.ndarray`:
        The destaggered variable.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
    
    """
    var_shape = var.shape
    num_dims = var.ndim
    stagger_dim_size = var_shape[stagger_dim]
    
    # Dynamically building the range slices to create the appropriate 
    # number of ':'s in the array accessor lists.
    # For example, for a 3D array, the calculation would be 
    # result = .5 * (var[:,:,0:stagger_dim_size-2] 
    #                    + var[:,:,1:stagger_dim_size-1])
    # for stagger_dim=2.  So, full slices would be used for dims 0 and 1, but 
    # dim 2 needs the special slice.  
    full_slice = slice(None)
    slice1 = slice(0, stagger_dim_size - 1, 1)
    slice2 = slice(1, stagger_dim_size, 1)
    
    # default to full slices
    dim_ranges_1 = [full_slice] * num_dims
    dim_ranges_2 = [full_slice] * num_dims
    
    # for the stagger dim, insert the appropriate slice range
    dim_ranges_1[stagger_dim] = slice1
    dim_ranges_2[stagger_dim] = slice2
    
    result = .5*(var[dim_ranges_1] + var[dim_ranges_2])
    
    return result

    