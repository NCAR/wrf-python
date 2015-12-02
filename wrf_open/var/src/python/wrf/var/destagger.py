
import numpy as n

__all__ = ["destagger", "destagger_windcomp", "destagger_winds"]

def destagger(var, stagger_dim):
    """ De-stagger the variable.  
    
    Arguments:
        - var is a numpy array for the variable 
        - stagger_dim is the dimension of the numpy array to de-stagger 
        (e.g. 0, 1, 2)
    
    """
    var_shape = var.shape
    num_dims = len(var_shape)
    stagger_dim_size = var_shape[stagger_dim]
    
    # Dynamically building the range slices to create the appropriate 
    # number of ':'s in the array accessor lists.
    # For example, for a 3D array, the calculation would be 
    # result = .5 * (var[:,:,0:stagger_dim_size-2] + var[:,:,1:stagger_dim_size-1])
    # for stagger_dim=2.  So, full slices would be used for dims 0 and 1, but 
    # dim 2 needs the special slice.  
    full_slice = slice(None, None, None)
    slice1 = slice(0, stagger_dim_size - 1, 1)
    slice2 = slice(1, stagger_dim_size, 1)
    
    # default to full slices
    dim_ranges_1 = [full_slice for x in xrange(num_dims)]
    dim_ranges_2 = [full_slice for x in xrange(num_dims)]
    
    # for the stagger dim, insert the appropriate slice range
    dim_ranges_1[stagger_dim] = slice1
    dim_ranges_2[stagger_dim] = slice2
    
    result = .5*(var[dim_ranges_1] + var[dim_ranges_2])
    
    return result

def destagger_windcomp(wrfnc, comp, timeidx=0):
    if comp.lower() == "u":
        wrfvar = "U"
        stagdim = 2
    elif comp.lower() == "v":
        wrfvar = "V"
        stagdim = 1
    elif comp.lower() == "w":
        wrfvar = "W"
        stagdim = 0
        
    wind_data = wrfnc.variables[wrfvar][timeidx,:,:,:]
    return destagger(wind_data, stagdim)

def destagger_winds(wrfnc, timeidx=0):
    return (destagger_windcomp(wrfnc, "u", timeidx),
            destagger_windcomp(wrfnc, "v", timeidx),
            destagger_windcomp(wrfnc, "w", timeidx))
    