from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants
from .extension import _tk, _rh, _cloudfrac
from .metadecorators import set_cloudfrac_metadata
from .util import extract_vars

@set_cloudfrac_metadata()
def get_cloudfrac(wrfin, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None):
    """Return the cloud fraction.
    
    The leftmost dimension of the returned array represents three different 
    quantities:
        
        - return_val[0,...] will contain LOW level cloud fraction
        - return_val[1,...] will contain MID level cloud fraction
        - return_val[2,...] will contain HIGH level cloud fraction
    
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
   
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        cloud fraction array whose leftmost dimension is 3 (LOW=0, MID=1, 
        HIGH=2). 
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    
    ncvars = extract_vars(wrfin, timeidx, ("P", "PB", "QVAPOR", "T"), 
                          method, squeeze, cache, meta=False,
                          _key=_key)
    
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    t = ncvars["T"]
    
    full_p = p + pb
    full_t = t + Constants.T_BASE
    
    tk = _tk(full_p, full_t)
    rh = _rh(qv, full_p, tk)
    
    return _cloudfrac(full_p, rh)
