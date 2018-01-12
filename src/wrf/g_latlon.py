from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .util import extract_vars, get_id
from .latlonutils import (_lat_varname, _lon_varname, _ll_to_xy, _xy_to_ll)
from .metadecorators import set_latlon_metadata


def get_lat(wrfin, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True, _key=None,
            stagger=None):
    """Return the two dimensional latitude coordinate variable.
    
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
            
        stagger (:obj:`str`): By default, the latitude is returned on the mass
            grid, but a staggered grid can be chosen with the following 
            options:
            
                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component, 
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component, 
                  which has a staggered south_north (y) dimension.
   
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        two dimensional latitude coordinate variable.
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    
    varname = _lat_varname(wrfin, stagger)
    lat_var = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                           meta, _key)
    
    return lat_var[varname]

        
def get_lon(wrfin, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True, _key=None,
            stagger=None):
    """Return the two dimensional longitude coordinate variable.
    
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
            
        stagger (:obj:`str`): By default, the longitude is returned on the mass
            grid, but a staggered grid can be chosen with the following 
            options:
            
                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component, 
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component, 
                  which has a staggered south_north (y) dimension.
   
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        two dimensional longitude coordinate variable.
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    
    varname = _lon_varname(wrfin, stagger)
    lon_var = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                           meta, _key)
    
    return lon_var[varname]

# TODO:  Do we need the user to know about method, squeeze, cache for this?

@set_latlon_metadata(xy=True) 
def ll_to_xy(wrfin, latitude, longitude, timeidx=0,  
             squeeze=True, meta=True, stagger=None, as_int=True):
    """Return the x,y coordinates for a specified latitude and longitude.
    
    The *latitude* and *longitude* arguments can be a single value or a 
    sequence of values.
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain the X (west_east) values.
        - return_val[1,...] will contain the Y (south_north) values.
    
    Args:
    
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types.
            
        latitude (:obj:`float` or sequence): A single latitude or a sequence 
            of latitude values to be converted.
            
        longitude (:obj:`float` or sequence): A single longitude or a sequence 
            of latitude values to be converted.
        
        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The 
            desired time index. This value can be a positive integer, 
            negative integer, or 
            :data:`wrf.ALL_TIMES` (an alias for None) to return 
            all times in the file or sequence. The default is 0.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        stagger (:obj:`str`): By default, the latitude is returned on the mass
            grid, but a staggered grid can be chosen with the following 
            options:
            
                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component, 
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component, 
                  which has a staggered south_north (y) dimension.
        
        as_int (:obj:`bool`): Set to True to return the x,y values as 
            :obj:`int`, otherwise they will be returned as :obj:`float`.
        
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        x,y coordinate value(s) whose leftmost dimension is 2 (0=X, 1=Y).
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    _key = get_id(wrfin)
    return _ll_to_xy(latitude, longitude, wrfin, timeidx, stagger, "cat", 
                     squeeze, None, _key, as_int, **{})


@set_latlon_metadata(xy=True) 
def ll_to_xy_proj(latitude, longitude, meta=True, squeeze=True, as_int=True,
                  map_proj=None, truelat1=None, truelat2=None, stand_lon=None, 
                  ref_lat=None, ref_lon=None, pole_lat=None, pole_lon=None, 
                  known_x=None, known_y=None, dx=None, dy=None, 
                  latinc=None, loninc=None):
    """Return the x, y coordinates for a specified latitude and longitude.
    
    The *latitude* and *longitude* arguments can be a single value or a 
    sequence of values.  This version of the ll_to_xy routine allows users 
    to manually specify projection parameters.  
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain the X (west_east) values.
        - return_val[1,...] will contain the Y (south_north) values.
    
    Args:
            
        latitude (:obj:`float` or sequence): A single latitude or a sequence 
            of latitude values to be converted.
            
        longitude (:obj:`float` or sequence): A single longitude or a sequence 
            of latitude values to be converted.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
        
        as_int (:obj:`bool`): Set to True to return the x,y values as 
            :obj:`int`, otherwise they will be returned as :obj:`float`.
            
        map_proj (:obj:`int`): Model projection [1=Lambert Conformal, 
            2=Polar Stereographic, 3=Mercator, 6=Lat-Lon].  Required.
        
        truelat1 (:obj:`float`): True latitude 1.  Required for 
            map_proj = 1, 2, 3 (defaults to 0 otherwise).
        
        truelat2 (:obj:`float`): True latitude 2.  Optional for 
            map_proj = 1 (defaults to 0 otherwise).
        
        stand_lon (:obj:`float`): Standard longitude. Required.
        
        ref_lat (:obj:`float`): A reference latitude.  Required. 
        
        ref_lon (:obj:`float`): A reference longitude.  Required.
        
        known_x (:obj:`float`): The known x-coordinate associated with 
            *ref_lon*. Required.
        
        known_y (:obj:`float`): The known y-coordinate associated with 
            *ref_lat*.  Required.
        
        pole_lat (:obj:`float`): Pole latitude. Optional for 
            *map_proj* = 6 (defaults to 90 otherwise).
        
        pole_lon (:obj:`float`): Pole longitude. Optional for 
            *map_proj* = 6 (defaults to 0 otherwise).
        
        dx (:obj:`float`): The x spacing in meters at the true latitude.  
            Required for *map_proj* = 1, 2, 3 (defaults to 0 otherwise).
        
        dy (:obj:`float`) - The y spacing in meters at the true latitude.  
            Required for *map_proj* = 1, 2, 3 (defaults to 0 otherwise).
        
        latinc (:obj:`float`): Required for *map_proj* = 6. Defined as:
            
            .. code-block:: python
            
                latinc = (dy*360.0)/2.0/Constants.PI/Constants.WRF_EARTH_RADIUS
        
        loninc (:obj:`float`): Required for *map_proj* = 6. Defined as:
            
            .. code-block:: python
            
                loninc = (dx*360.0)/2.0/Constants.PI/Constants.WRF_EARTH_RADIUS

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        x,y coordinate value(s) whose leftmost dimension is 2 (0=X, 1=Y).
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    loc = locals()
    projparams = {name : loc[name] for name in ("map_proj", "truelat1", 
                                            "truelat2", "stand_lon", "ref_lat",
                                            "ref_lon", "pole_lat", "pole_lon",
                                            "known_x", "known_y", "dx", "dy",
                                            "latinc", "loninc")}

    return _ll_to_xy(latitude, longitude, None, 0, True, "cat", squeeze, None,
                     None, as_int, **projparams)
    
    
@set_latlon_metadata(xy=False) 
def xy_to_ll(wrfin, x, y, timeidx=0, stagger=None, squeeze=True, meta=True):
    """Return the latitude and longitude for specified x,y coordinates.
    
    The *x* and *y* arguments can be a single value or a sequence of values.
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain the latitude values.
        - return_val[1,...] will contain the longitude values.
    
    Args:
    
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types.
            
        x (:obj:`float` or sequence): A single x-coordinate or a sequence 
            of x-coordinate values to be converted.
            
        y (:obj:`float` or sequence): A single y-coordinate or a sequence 
            of y-coordinate values to be converted.
        
        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The 
            desired time index. This value can be a positive integer, 
            negative integer, or 
            :data:`wrf.ALL_TIMES` (an alias for None) to return 
            all times in the file or sequence. The default is 0.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        stagger (:obj:`str`): By default, the latitude is returned on the mass
            grid, but a staggered grid can be chosen with the following 
            options:
            
                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component, 
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component, 
                  which has a staggered south_north (y) dimension.
        
        as_int (:obj:`bool`): Set to True to return the x,y values as 
            :obj:`int`, otherwise they will be returned as :obj:`float`.
        
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        latitude and longitude values whose leftmost dimension is 2 
        (0=latitude, 1=longitude).
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    _key = get_id(wrfin)
    return _xy_to_ll(x, y, wrfin, timeidx, stagger, "cat", True, None, 
                     _key, **{})
 
    
@set_latlon_metadata(xy=False) 
def xy_to_ll_proj(x, y, meta=True, squeeze=True, map_proj=None, truelat1=None, 
                  truelat2=None, stand_lon=None, ref_lat=None, ref_lon=None, 
                  pole_lat=None, pole_lon=None, known_x=None, known_y=None, 
                  dx=None, dy=None, latinc=None, loninc=None):
    """Return the latitude and longitude for the specified x,y coordinates.
    
    The *x* and *y* arguments can be a single value or a 
    sequence of values.  This version of the xy_to_ll routine allows users 
    to manually specify map projection parameters.  
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain the latitude values.
        - return_val[1,...] will contain the longitude values.
    
    Args:
            
        x (:obj:`float` or sequence): A single x-coordinate or a sequence 
            of x-coordinate values to be converted.
            
        y (:obj:`float` or sequence): A single y-coordinate or a sequence 
            of y-coordinate values to be converted.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
        
        as_int (:obj:`bool`): Set to True to return the x,y values as 
            :obj:`int`, otherwise they will be returned as :obj:`float`.
            
        map_proj (:obj:`int`): Model projection [1=Lambert Conformal, 
            2=Polar Stereographic, 3=Mercator, 6=Lat-Lon].  Required.
        
        truelat1 (:obj:`float`): True latitude 1.  Required for 
            map_proj = 1, 2, 3 (defaults to 0 otherwise).
        
        truelat2 (:obj:`float`): True latitude 2.  Optional for 
            map_proj = 1 (defaults to 0 otherwise).
        
        stand_lon (:obj:`float`): Standard longitude. Required.
        
        ref_lat (:obj:`float`): A reference latitude.  Required. 
        
        ref_lon (:obj:`float`): A reference longitude.  Required.
        
        known_x (:obj:`float`): The known x-coordinate associated with 
            *ref_lon*. Required.
        
        known_y (:obj:`float`): The known y-coordinate associated with 
            *ref_lat*.  Required.
        
        pole_lat (:obj:`float`): Pole latitude. Optional for 
            *map_proj* = 6 (defaults to 90 otherwise).
        
        pole_lon (:obj:`float`): Pole longitude. Optional for 
            *map_proj* = 6 (defaults to 0 otherwise).
        
        dx (:obj:`float`): The x spacing in meters at the true latitude.  
            Required for *map_proj* = 1, 2, 3 (defaults to 0 otherwise).
        
        dy (:obj:`float`) - The y spacing in meters at the true latitude.  
            Required for *map_proj* = 1, 2, 3 (defaults to 0 otherwise).
        
        latinc (:obj:`float`): Required for *map_proj* = 6. Defined as:
            
            .. code-block:: python
            
                latinc = (dy*360.0)/2.0/Constants.PI/Constants.WRF_EARTH_RADIUS
        
        loninc (:obj:`float`): Required for *map_proj* = 6. Defined as:
            
            .. code-block:: python
            
                loninc = (dx*360.0)/2.0/Constants.PI/Constants.WRF_EARTH_RADIUS

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        latitude and longitude values whose leftmost dimension is 2 
        (0=latitude, 1=longitude).
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    """
    loc = locals()
    projparams = {name : loc[name] for name in ("map_proj", "truelat1", 
                                            "truelat2", "stand_lon", "ref_lat",
                                            "ref_lon", "pole_lat", "pole_lon",
                                            "known_x", "known_y", "dx", "dy",
                                            "latinc", "loninc")}
    return _xy_to_ll(x, y, None, 0, None, "cat", squeeze, None, None,
                     **projparams)

    
    
    
    
    
    
    