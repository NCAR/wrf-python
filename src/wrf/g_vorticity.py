from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .extension import _avo, _pvo
from .util import extract_vars, extract_global_attrs
from .metadecorators import copy_and_set_metadata


@copy_and_set_metadata(copy_varname="T", name="avo", 
                       description="absolute vorticity",
                       units="10-5 s-1")
def get_avo(wrfin, timeidx=0, method="cat", squeeze=True, cache=None, 
            meta=True, _key=None):
    """Return the absolute vorticity.
    
    This function extracts the necessary variables from the NetCDF file 
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
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The absolute 
        vorticity. If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
    
    """
    ncvars = extract_vars(wrfin, timeidx, ("U", "V", "MAPFAC_U",
                                           "MAPFAC_V", "MAPFAC_M",
                                           "F"),
                          method, squeeze, cache, meta=False, _key=_key)
    
    attrs = extract_global_attrs(wrfin, attrs=("DX", "DY"))
    u = ncvars["U"]
    v = ncvars["V"]
    msfu = ncvars["MAPFAC_U"]
    msfv = ncvars["MAPFAC_V"]
    msfm = ncvars["MAPFAC_M"]
    cor = ncvars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    return _avo(u, v, msfu, msfv, msfm, cor, dx, dy)


@copy_and_set_metadata(copy_varname="T", name="pvo", 
                       description="potential vorticity",
                       units="PVU")
def get_pvo(wrfin, timeidx=0, method="cat", squeeze=True, cache=None, 
            meta=True, _key=None):
    """Return the potential vorticity.
    
    This function extracts the necessary variables from the NetCDF file 
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
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The potential 
        vorticity. If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
    
    """
    ncvars = extract_vars(wrfin, timeidx, ("U", "V", "T", "P",
                                           "PB", "MAPFAC_U",
                                           "MAPFAC_V", "MAPFAC_M",
                                           "F"),
                          method, squeeze, cache, meta=False, _key=_key)
    attrs = extract_global_attrs(wrfin, attrs=("DX", "DY"))
    
    u = ncvars["U"]
    v = ncvars["V"]
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    msfu = ncvars["MAPFAC_U"]
    msfv = ncvars["MAPFAC_V"]
    msfm = ncvars["MAPFAC_M"]
    cor = ncvars["F"]
    
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    full_t = t + 300
    full_p = p + pb
    
    return _pvo(u, v, full_t, full_p, msfu, msfv, msfm, cor, dx, dy)
    