from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

#from .extension import computedbz,computetk
from .extension import _dbz, _tk
from .constants import Constants
from .util import extract_vars
from .metadecorators import copy_and_set_metadata


@copy_and_set_metadata(copy_varname="T", name="dbz", 
                       description="radar reflectivity",
                       units="dBZ")
def get_dbz(wrfin, timeidx=0, method="cat", 
            squeeze=True, cache=None, meta=True, _key=None,
            use_varint=False, use_liqskin=False):
    """Return the simulated radar reflectivity.
    
    This functions extracts the necessary variables from the NetCDF file 
    object in order to perform the calculation.
    
    Args:
    
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types.
        
        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The 
            desired time index. This value can be a positive integer, 
            negative integer, or 
            :data:`wrf.ALL_TIMES` (an alias for None) to return 
            all times in the file or sequence. The default is 0.
        
        method (:obj:`str`, optional): The aggregation method to use for 
            sequences.  Must be either 'cat' or 'join'.  
            'cat' combines the data along the Time dimension.  
            'join' creates a new dimension for the file index.  
            The default is 'cat'.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        cache (:obj:`dict`, optional): A dictionary of (varname, ndarray) 
            that can be used to supply pre-extracted NetCDF variables to the 
            computational routines.  It is primarily used for internal 
            purposes, but can also be used to improve performance by 
            eliminating the need to repeatedly extract the same variables 
            used in multiple diagnostics calculations, particularly when using 
            large sequences of files. 
            Default is None.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        _key (:obj:`int`, optional): A caching key. This is used for internal 
            purposes only.  Default is None.
            
        use_varint (:obj:`bool`, optional): When set to False, 
            the intercept parameters are assumed constant 
            (as in MM5's Reisner-2 bulk microphysical scheme). 
            When set to True, the variable intercept 
            parameters are used as in the more recent version of Reisner-2 
            (based on Thompson, Rasmussen, and Manning, 2004, Monthly weather 
            Review, Vol. 132, No. 2, pp. 519-542.).  
            
        use_liqskin (:obj:`bool`, optional): When set to True, frozen particles 
            that are at a temperature above freezing are assumed to scatter 
            as a liquid particle.  Set to False to disable.
   
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The simulated 
        radar reflectivity. 
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    varnames = ("T", "P", "PB", "QVAPOR", "QRAIN")
    ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    qr = ncvars["QRAIN"]
    
    try:
        snowvars = extract_vars(wrfin, timeidx, "QSNOW", 
                                method, squeeze, cache, meta=False,
                                _key=_key)
    except KeyError:
        qs = np.zeros(qv.shape, qv.dtype)
    else:
        qs = snowvars["QSNOW"]
    
    try:
        graupvars = extract_vars(wrfin, timeidx, "QGRAUP", 
                                 method, squeeze, cache, meta=False,
                                 _key=_key)
    except KeyError:
        qg = np.zeros(qv.shape, qv.dtype)
    else:
        qg = graupvars["QGRAUP"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    # If qsnow is not all 0, set sn0 to 1
    sn0 = 1 if qs.any() else 0
    ivarint = 1 if use_varint else 0
    iliqskin = 1 if use_liqskin else 0
    
    return _dbz(full_p, tk, qv, qr, qs, qg, sn0, ivarint, iliqskin)


@copy_and_set_metadata(copy_varname="T", name="max_dbz", 
                       remove_dims=("bottom_top",),
                       description="maximum radar reflectivity",
                       units="dBZ",
                       MemoryOrder="XY")
def get_max_dbz(wrfin, timeidx=0, method="cat", 
                squeeze=True, cache=None, meta=True, _key=None,
                use_varint=False, use_liqskin=False):
    """Return the maximum simulated radar reflectivity.
    
    This function returns the maximum reflectivity found in the column for 
    each grid point.
    
    This functions extracts the necessary variables from the NetCDF file 
    object in order to perform the calculation.
    
    Args:
    
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types.
        
        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The 
            desired time index. This value can be a positive integer, 
            negative integer, or 
            :data:`wrf.ALL_TIMES` (an alias for None) to return 
            all times in the file or sequence. The default is 0.
        
        method (:obj:`str`, optional): The aggregation method to use for 
            sequences.  Must be either 'cat' or 'join'.  
            'cat' combines the data along the Time dimension.  
            'join' creates a new dimension for the file index.  
            The default is 'cat'.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        cache (:obj:`dict`, optional): A dictionary of (varname, ndarray) 
            that can be used to supply pre-extracted NetCDF variables to the 
            computational routines.  It is primarily used for internal 
            purposes, but can also be used to improve performance by 
            eliminating the need to repeatedly extract the same variables 
            used in multiple diagnostics calculations, particularly when using 
            large sequences of files. 
            Default is None.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        _key (:obj:`int`, optional): A caching key. This is used for internal 
            purposes only.  Default is None.
            
        use_varint (:obj:`bool`, optional): When set to False, 
            the intercept parameters are assumed constant 
            (as in MM5's Reisner-2 bulk microphysical scheme). 
            When set to True, the variable intercept 
            parameters are used as in the more recent version of Reisner-2 
            (based on Thompson, Rasmussen, and Manning, 2004, Monthly weather 
            Review, Vol. 132, No. 2, pp. 519-542.).  
            
        use_liqskin (:obj:`bool`, optional): When set to True, frozen particles 
            that are at a temperature above freezing are assumed to scatter 
            as a liquid particle.  Set to False to disable.
   
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The maximum 
        simulated radar reflectivity. 
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    return np.amax(get_dbz(wrfin, timeidx, method, squeeze, cache, meta,
                           _key, use_varint, use_liqskin), 
                  axis=-3)

