from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .decorators import extract_and_transpose
from .metadecorators import set_destag_metadata


@set_destag_metadata()
@extract_and_transpose(do_transpose=False)
def destagger(var, stagger_dim, meta=False):
    """ De-stagger the variable.  
    
    Arguments:
        - var is a numpy array for the variable 
        - stagger_dim is the dimension of the numpy array to de-stagger 
        (e.g. 0, 1, 2).  Note:  negative values are acceptable to choose
        a dimensions from the right hand side (e.g. -1, -2, -3)
        - meta - set to True to include 'var' metadata
    
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

    