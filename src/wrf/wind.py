from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

from .extension import _wspd, _wdir
from .destag import destagger
from .util import extract_vars, either
from .decorators import convert_units
from .metadecorators import set_wind_metadata


@convert_units("wind", "m s-1")
def _calc_wspd(u, v, units="m s-1"):
    """Return the wind speed.
    
    Args:
    
        u (:class:`numpy.ndarray`): The u component of the wind.
        
        v (:class:`numpy.ndarray`): The v component of the wind.
        
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'uvmet_wspd_wdir'.  
            Default is 'm s-1'.
            
    Returns:
    
        :class:`numpy.ndarray`: The wind speed.
    
    """
    return _wspd(u, v)


def _calc_wdir(u, v):
    """Return the wind direction.
    
    Args:
    
        u (:class:`numpy.ndarray`): The u component of the wind.
        
        v (:class:`numpy.ndarray`): The v component of the wind.
            
    Returns:
    
        :class:`numpy.ndarray`: The wind direction.
    
    """
    return _wdir(u, v)


def _calc_wspd_wdir(u, v, two_d, units):
    """Return the wind speed and wind direction.
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain WSPD
        - return_val[1,...] will contain WDIR
    
    Args:
    
        u (:class:`numpy.ndarray`): The u component of the wind.
        
        v (:class:`numpy.ndarray`): The v component of the wind.
        
        two_d (:obj:`bool`): Set to True if the u,v wind components are 
            for a two-dimensional array (no height dimension).  Otherwise,
            set to False.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'uvmet_wspd_wdir'.  
            Default is 'm s-1'.
            
    Returns:
    
        :class:`numpy.ndarray`: The wind speed and wind direction, whose 
        leftmost dimension is 2 (0=WSPD, 1=WDIR).
    
    """

    wspd = _calc_wspd(u, v, units)
    wdir = _calc_wdir(u, v)
    
    try:
        fill = wspd.fill_value
    except AttributeError:
        fill = None

    idx_end = -2 if two_d else -3
    
    outdims = [2] + list(wspd.shape[0:idx_end]) + list(wspd.shape[idx_end:])
    
    result = np.zeros(outdims, wspd.dtype)
    
    idxs0 = ((0,Ellipsis, slice(None), slice(None), slice(None)) 
            if not two_d else 
            (0, Ellipsis, slice(None), slice(None)))
    
    idxs1 = ((1, Ellipsis, slice(None), slice(None), slice(None)) 
            if not two_d else 
            (1, Ellipsis, slice(None), slice(None)))
    
    result[idxs0] = wspd[:]
    result[idxs1] = wdir[:]
    
    if fill is not None:
        result = np.ma.masked_equal(result, fill)
        
    return result


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="ua",
                   description="destaggered u-wind component",
                   wind_ncvar=True, 
                   two_d=False, 
                   wspd_wdir=False)
@convert_units("wind", "m s-1")
def get_u_destag(wrfin, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None,
                 units="m s-1"):
    """Return the u-component of the wind on mass points.
    
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
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'wspd_wdir'.  
            Default is 'm s-1'.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The u-component 
        of the wind.  
        If  xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    varname = either("U", "UU")(wrfin)
    u_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = destagger(u_vars[varname], -1)
    
    return u


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="va",
                   description="destaggered v-wind component",
                   two_d=False,
                   wind_ncvar=True, 
                   wspd_wdir=False)
@convert_units("wind", "m s-1")
def get_v_destag(wrfin, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None,
                 units="m s-1"):
    """Return the v-component of the wind on mass points.
    
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
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'wspd_wdir'.  
            Default is 'm s-1'.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The v-component 
        of the wind.  
        If  xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    varname = either("V", "VV")(wrfin)
    v_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = destagger(v_vars[varname], -2)
    
    return v


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="wa",
                   description="destaggered w-wind component",
                   two_d=False,
                   wind_ncvar=True, 
                   wspd_wdir=False)
@convert_units("wind", "m s-1")
def get_w_destag(wrfin, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None,
                 units="m s-1"):
    """Return the w-component of the wind on mass points.
    
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
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'wspd_wdir'.  
            Default is 'm s-1'.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The w-component 
        of the wind.  
        If  xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    w_vars = extract_vars(wrfin, timeidx, "W", method, squeeze, cache, 
                          meta=False, _key=_key)
    w = destagger(w_vars["W"], -3)
    return w


@set_wind_metadata(copy_varname=either("P", "PRES"), 
                   name="wspd_wdir",
                   description="wspd,wdir in projection space",
                   two_d=False, 
                   wspd_wdir=True)
def get_destag_wspd_wdir(wrfin, timeidx=0, method="cat", 
                         squeeze=True, cache=None, meta=True, _key=None,
                         units="m s-1"):
    """Return the wind speed and wind direction for the wind in the projected
    coordinate space.  
    
    The wind speed and direction from this function will be relative to the 
    grid.  This function should not be used to compare with observations, 
    instead use :meth:`wrf.uvmet10_wspd_wdir` to get the earth relative wind 
    speed and direction.
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain WSPD
        - return_val[1,...] will contain WDIR
        
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
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'wspd_wdir'.  
            Default is 'm s-1'.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind speed and 
        wind direction for the wind in projected space, whose 
        leftmost dimensions is 2 (0=WSPD, 1=WDIR).  If 
        xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    varname = either("U", "UU")(wrfin)
    u_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = destagger(u_vars[varname], -1)
    
    varname = either("V", "VV")(wrfin)
    v_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = destagger(v_vars[varname], -2)
    
    return _calc_wspd_wdir(u, v, False, units)


@set_wind_metadata(copy_varname=either("PSFC", "F"), 
                   name="wspd_wdir10",
                   description="10m wspd,wdir in projection space",
                   two_d=False, 
                   wspd_wdir=True)
def get_destag_wspd_wdir10(wrfin, timeidx=0, method="cat", 
                           squeeze=True, cache=None, meta=True, _key=None, 
                           units="m s-1"):
    """Return the wind speed and wind direction for the 10m wind in 
    projected coordinate space.
    
    The wind speed and direction from this function will be relative to the 
    grid.  This function should not be used to compare with observations, 
    instead use :meth:`wrf.uvmet10_wspd_wdir` to get the earth relative wind 
    speed and direction. 
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain WSPD10
        - return_val[1,...] will contain WDIR10
        
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
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 
            'wspd_wdir10'. Default is 'm s-1'.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind speed and 
        wind direction for the 10m wind in projected space, whose 
        leftmost dimensions is 2 (0=WSPD10, 1=WDIR10).  If 
        xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    
    varname = either("U10", "UU")(wrfin)
    u_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = (u_vars[varname] if varname == "U10" else 
         destagger(u_vars[varname][...,0,:,:], -1)) 
    
    varname = either("V10", "VV")(wrfin)
    v_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = (v_vars[varname] if varname == "V10" else 
         destagger(v_vars[varname][...,0,:,:], -2))
    
    return _calc_wspd_wdir(u, v, True, units)

