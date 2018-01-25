from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants
from .destag import destagger
from .extension import _omega, _tk
from .util import extract_vars
from .metadecorators import copy_and_set_metadata


@copy_and_set_metadata(copy_varname="T", name="omega", 
                       description="omega",
                       units="Pa s-1")
def get_omega(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
              meta=True, _key=None):
    """Return Omega.
    
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
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: Omega. 
        If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
    
    """
    varnames=("T", "P", "W", "PB", "QVAPOR")
    ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    w = ncvars["W"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    wa = destagger(w, -3)
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    omega = _omega(qv, tk, wa, full_p)
    
    return omega
    