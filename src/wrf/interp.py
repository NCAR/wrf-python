from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np
import numpy.ma as ma

from .extension import (_interpz3d, _vertcross, _interpline, _smooth2d, 
                        _monotonic, _vintrp)

from .metadecorators import set_interp_metadata
from .util import extract_vars, is_staggered, get_id, to_np
from .py3compat import py3range
from .interputils import get_xy, get_xy_z_params, to_xy_coords
from .constants import Constants, ConversionFactors
from .terrain import get_terrain
from .geoht import get_height
from .temp import get_theta, get_temp, get_eth
from .pressure import get_pressure


#  Note:  Extension decorator is good enough to handle left dims
@set_interp_metadata("horiz")
def interplevel(field3d, vert, desiredlev, missing=Constants.DEFAULT_FILL, 
                meta=True):
    """Return the three-dimensional field interpolated to a horizontal plane 
    at the specified vertical level.

    Args:
    
        field3d (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            three-dimensional field to interpolate, with the rightmost 
            dimensions of nz x ny x nx.
    
        vert (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            three-dimensional array for the vertical coordinate, typically 
            pressure or height. This array must have the same dimensionality
            as *field3d*.
        
        desiredlev (:obj:`float`): The desired vertical level.  
            Must be in the same units as the *vert* parameter.
        
        missing (:obj:`float`): The fill value to use for the output.  
            Default is :data:`wrf.Constants.DEFAULT_FILL`.
        
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        interpolated variable.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    Example:
    
        Interpolate Geopotential Height to 500 hPa
    
        .. code-block:: python
        
            from netCDF4 import Dataset
            from wrf import getvar, interplevel
        
            wrfin = Dataset("wrfout_d02_2010-06-13_21:00:00")
            
            p = getvar(wrfin, "pressure")
            ht = getvar(wrfin, "z", units="dm")
            
            ht_500 = interplevel(ht, p, 500.0)
        
    
    """
    # Some fields like uvmet have an extra left dimension for the product
    # type, we'll handle that iteration here.
    multi = True if field3d.ndim - vert.ndim == 1 else False
    
    if not multi:
        result = _interpz3d(field3d, vert, desiredlev, missing)
    else:
        outshape = field3d.shape[0:-3] + field3d.shape[-2:]
        result = np.empty(outshape, dtype=field3d.dtype)
            
        for i in py3range(field3d.shape[0]):
            result[i,:] = (
                _interpz3d(field3d[i,:], vert, desiredlev, missing)[:])
            
    return ma.masked_values (result, missing)


@set_interp_metadata("cross")
def vertcross(field3d, vert, levels=None, missing=Constants.DEFAULT_FILL,
              wrfin=None, timeidx=0, stagger=None, projection=None, 
              pivot_point=None, angle=None,
              start_point=None, end_point=None,
              latlon=False, cache=None, meta=True):
    """Return the vertical cross section for a three-dimensional field.
    
    The cross section is defined by a horizontal line through the domain.  
    This horizontal line is defined by either including the 
    *pivot_point* and *angle* parameters, or the *start_point* and 
    *end_point* parameters. The *pivot_point*, *start_point*, and *end_point* 
    coordinates can be defined in either x,y or latitude,longitude space.  
    If latitude,longitude coordinates are used, then a WRF input file or 
    map projection must also be specified.
    
    The vertical levels for the cross section are fixed if *levels* is not 
    specified, and are determined by dividing the vertical coordinate in to 
    grid boxes of roughly 1% of the maximum vertical distance from top to 
    bottom.  Otherwise, the *levels* argument can be used to specify specific 
    vertical levels.  If all vertical levels are desired, use the raw 
    :meth:`wrf.interp2dxy` function.
    
    See Also:
    
        :meth:`wrf.interp2dxy`
    
    Args:
    
        field3d (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            three-dimensional field to interpolate, whose 
            rightmost dimensions are nz x ny x nx.
            
        vert (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            three-dimensional variable for the vertical coordinate, typically 
            pressure or height. This array must have the same dimensionality
            as *field3d*
            
        levels (sequence, optional): A sequence of :obj:`float` for the desired
            vertical levels in the output array.  Must be in the same units 
            as *vert*.  If None, a fixed set of vertical levels is provided.  
            Default is None.
            
        missing (:obj:`float`): The fill value to use for the output.  
            Default is :data:`wrf.Constants.DEFAULT_FILL`.
            
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types. This is used
            to obtain the map projection when using latitude,longitude 
            coordinates. Default is None.
            
        timeidx (:obj:`int`, optional): The 
            desired time index when obtaining map boundary information 
            from moving nests. This value can be a positive or negative integer.
            Only required when *wrfin* is specified and the nest is moving.
            Currently, :data:`wrf.ALL_TIMES` is not supported.
            Default is 0.
            
        stagger (:obj:`str`): If using latitude, longitude coordinate pairs 
            for *start_point*, *end_point*, or *pivot_point*, 
            set the appropriate grid staggering type for *field3d*. By default,
            the mass grid is used.  The options are:
            
                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component, 
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component, 
                  which has a staggered south_north (y) dimension.
            
        projection (:class:`wrf.WrfProj` subclass, optional): The map 
            projection object to use when working with latitude, longitude 
            coordinates, and must be specified if *wrfin* is None. Default 
            is None.
            
        pivot_point (:class:`wrf.CoordPair`, optional): A coordinate pair for 
            the pivot point, which indicates the location through which 
            the plane will pass. Must also specify *angle*.  The coordinate 
            pair can be in x,y grid coordinates or latitude, longitude 
            coordinates.  If using latitude, longitude coordinates, then 
            either *wrfin* or *projection* must be specified to obtain the 
            map projection.  Default is None.
        
        angle (:obj:`float`, optional): Only valid for cross sections where 
            a plane will be plotted through 
            a given point on the model domain. 0.0 represents a S-N cross 
            section.  90.0 is a W-E cross section. 
            
        start_point (:class:`wrf.CoordPair`, optional): A coordinate pair 
            which indicates the start location through which the plane will 
            pass. Must also specify *end_point*. The coordinate 
            pair can be in x,y grid coordinates or latitude, longitude 
            coordinates.  If using latitude, longitude coordinates, then 
            either *wrfin* or *projection* must be specified to obtain the 
            map projection.  Default is None.
            
        end_point (:class:`wrf.CoordPair`, optional): A coordinate pair 
            which indicates the end location through which the plane will 
            pass. Must also specify *end_point*. The coordinate 
            pair can be in x,y grid coordinates or latitude, longitude 
            coordinates.  If using latitude, longitude coordinates, then 
            either *wrfin* or *projection* must be specified to obtain the 
            map projection.  Default is None.
            
        latlon (:obj:`bool`, optional): Set to True to also interpolate the 
            two-dimensional latitude and longitude coordinates along the same 
            horizontal line and include this information in the metadata 
            (if enabled).  This can be helpful for plotting.  Default is False.
             
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
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`:
        The interpolated variable.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    """
    if timeidx is None:
        raise ValueError("'timeidx' must be a positive or negative integer")
    
    
    # Some fields like uvmet have an extra left dimension for the product
    # type, we'll handle that iteration here.
    multi = True if field3d.ndim - vert.ndim == 1 else False
    
    try:
        xy = cache["xy"]
        var2dz = cache["var2dz"]
        z_var2d = cache["z_var2d"]
    except (KeyError, TypeError):
        # Convert Lat/Lon points to grid points
        start_point_xy = None
        end_point_xy = None
        pivot_point_xy = None
        
        if pivot_point is not None:
            if pivot_point.lat is not None and pivot_point.lon is not None:
                xy_coords = to_xy_coords(pivot_point, wrfin, timeidx, 
                                         stagger, projection)
                pivot_point_xy = (xy_coords.x, xy_coords.y)
            else:
                pivot_point_xy = (pivot_point.x, pivot_point.y)
                
        if start_point is not None and end_point is not None:
            if start_point.lat is not None and start_point.lon is not None:
                xy_coords = to_xy_coords(start_point, wrfin, timeidx, 
                                         stagger, projection)
                start_point_xy = (xy_coords.x, xy_coords.y)
            else:
                start_point_xy = (start_point.x, start_point.y)
                
            if end_point.lat is not None and end_point.lon is not None:
                xy_coords = to_xy_coords(end_point, wrfin, timeidx, 
                                         stagger, projection)
                end_point_xy = (xy_coords.x, xy_coords.y)
            else:
                end_point_xy = (end_point.x, end_point.y)
                
        xy, var2dz, z_var2d = get_xy_z_params(to_np(vert), pivot_point_xy, 
                                              angle, start_point_xy, 
                                              end_point_xy, levels)
    
    if not multi:
        result = _vertcross(field3d, xy, var2dz, z_var2d, missing)
    else:
        outshape = field3d.shape[0:-3] + (z_var2d.shape[0], xy.shape[0])
        result = np.empty(outshape, dtype=field3d.dtype)
            
        for i in py3range(field3d.shape[0]):
            result[i,:] = _vertcross(field3d[i,:], xy, var2dz, z_var2d, 
                                     missing)[:]
    
    return ma.masked_values(result, missing)


@set_interp_metadata("line")
def interpline(field2d, pivot_point=None, 
               wrfin=None, timeidx=0, stagger=None, projection=None,
               angle=None, start_point=None,
               end_point=None, latlon=False, 
               cache=None, meta=True):
    """Return the two-dimensional field interpolated along a line.
      
    This line is defined by either including the 
    *pivot_point* and *angle* parameters, or the *start_point* and 
    *end_point* parameters. The *pivot_point*, *start_point*, and *end_point* 
    coordinates can be defined in either x,y or latitude,longitude space.  
    If latitude,longitude coordinates are used, then a WRF input file or 
    map projection must also be specified.
    
    Args:
    
        field2d (:class:`xarray.DataArray` or :class:`numpy.ndarray`):
            A two-dimensional field.
    
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types. This is used
            to obtain the map projection when using latitude,longitude 
            coordinates. Should not be used when working with x,y 
            coordinates. Default is None.
            
        timeidx (:obj:`int`, optional): The 
            desired time index when obtaining map boundary information 
            from moving nests. This value can be a positive or negative integer.
            Only required when *wrfin* is specified and the nest is moving.
            Currently, :data:`wrf.ALL_TIMES` is not supported.
            Default is 0. 
            
        stagger (:obj:`str`): If using latitude, longitude coordinate pairs 
            for *start_point*, *end_point*, or *pivot_point*, 
            set the appropriate grid staggering type for *field2d*. By default,
            the mass grid is used.  The options are:
            
                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component, 
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component, 
                  which has a staggered south_north (y) dimension.
            
        projection (:class:`wrf.WrfProj`, optional): The map 
            projection object to use when working with latitude, longitude 
            coordinates, and must be specified if *wrfin* is None. Should 
            not be used when working with x,y coordinates.  Default 
            is None.
            
        pivot_point (:class:`wrf.CoordPair`, optional): A coordinate pair for 
            the pivot point, which indicates the location through which 
            the plane will pass. Must also specify *angle*.  The coordinate 
            pair can be in x,y grid coordinates or latitude, longitude 
            coordinates.  If using latitude, longitude coordinates, then 
            either *wrfin* or *projection* must be specified to obtain the 
            map projection.  Default is None.
        
        angle (:obj:`float`, optional): Only valid for cross sections where 
            a plane will be plotted through 
            a given point on the model domain. 0.0 represents a S-N cross 
            section.  90.0 is a W-E cross section. 
            
        start_point (:class:`wrf.CoordPair`, optional): A coordinate pair 
            which indicates the start location through which the plane will 
            pass. Must also specify *end_point*. The coordinate 
            pair can be in x,y grid coordinates or latitude, longitude 
            coordinates.  If using latitude, longitude coordinates, then 
            either *wrfin* or *projection* must be specified to obtain the 
            map projection.  Default is None.
            
        end_point (:class:`wrf.CoordPair`, optional): A coordinate pair 
            which indicates the end location through which the plane will 
            pass. Must also specify *end_point*. The coordinate 
            pair can be in x,y grid coordinates or latitude, longitude 
            coordinates.  If using latitude, longitude coordinates, then 
            either *wrfin* or *projection* must be specified to obtain the 
            map projection.  Default is None.
            
        latlon (:obj:`bool`, optional): Set to True to also interpolate the 
            two-dimensional latitude and longitude coordinates along the same 
            horizontal line and include this information in the metadata 
            (if enabled).  This can be helpful for plotting.  Default is False.
            
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
    
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`:
        The interpolated variable.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
        
    """
    if timeidx is None:
        raise ValueError("'timeidx' must be a positive or negative integer")
    
    try:
        xy = cache["xy"]
    except (KeyError, TypeError):
        start_point_xy = None
        end_point_xy = None
        pivot_point_xy = None
        
        if pivot_point is not None:
            if pivot_point.lat is not None and pivot_point.lon is not None:
                xy_coords = to_xy_coords(pivot_point, wrfin, timeidx, 
                                         stagger, projection)
                pivot_point_xy = (xy_coords.x, xy_coords.y)
            else:
                pivot_point_xy = (pivot_point.x, pivot_point.y)
                
        if start_point is not None and end_point is not None:
            if start_point.lat is not None and start_point.lon is not None:
                xy_coords = to_xy_coords(start_point, wrfin, timeidx, 
                                         stagger, projection)
                start_point_xy = (xy_coords.x, xy_coords.y)
            else:
                start_point_xy = (start_point.x, start_point.y)
                
            if end_point.lat is not None and end_point.lon is not None:
                xy_coords = to_xy_coords(end_point, wrfin, timeidx, 
                                         stagger, projection)
                end_point_xy = (xy_coords.x, xy_coords.y)
            else:
                end_point_xy = (end_point.x, end_point.y)
                
        xy = get_xy(field2d, pivot_point_xy, angle, start_point_xy, 
                    end_point_xy)
        
    return _interpline(field2d, xy)


@set_interp_metadata("vinterp")
def vinterp(wrfin, field, vert_coord, interp_levels, extrapolate=False, 
            field_type=None, log_p=False, timeidx=0, method="cat", 
            squeeze=True, cache=None, meta=True):
    """Return the field vertically interpolated to the given the type of 
    surface and a set of new levels.
    
    Args:
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`,\
        or an iterable):
            WRF-ARW NetCDF data as a :class:`netCDF4.Dataset`, 
            :class:`Nio.NioFile` or an iterable sequence of the 
            aforementioned types.
            
        field (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            three-dimensional field.
            
        vert_coord (:obj:`str`): A string indicating the vertical coordinate 
            type to interpolate to.
            
            Valid strings are: 
                * 'pressure', 'pres', 'p': pressure [hPa]
                * 'ght_msl': grid point height msl [km]
                * 'ght_agl': grid point height agl [km]
                * 'theta', 'th': potential temperature [K]
                * 'theta-e', 'thetae', 'eth': equivalent potential temperature \
                [K]
                
        interp_levels (sequence): A 1D sequence of vertical levels to 
            interpolate to.
            
        extrapolate (:obj:`bool`, optional): Set to True to extrapolate 
            values below ground.  Default is False.
            
        field_type (:obj:`str`, optional): 
            The type of field.  Default is None.
            
            Valid strings are:
                * 'none': None
                * 'pressure', 'pres', 'p': pressure
                * 'z', 'ght': geopotential height
                * 'tc': temperature [degC]
                * 'tk': temperature [K]
                * 'theta', 'th': potential temperature [K]
                * 'theta-e', 'thetae', 'eth': equivalent potential temperature
            
        log_p (:obj:`bool`, optional): Use the log of the pressure for 
            interpolation instead of pressure. Default is False.
            
        timeidx (:obj:`int`, optional):
            The time index to use when extracting auxiallary variables used in 
            the interpolation.  This value must be set to match the same value
            used when the `field` variable was extracted.  Default is 0.
            
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
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`:
        The interpolated variable.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    """
    _key = get_id(wrfin)
    
    # Remove case sensitivity
    field_type = field_type.lower() if field_type is not None else "none"
    vert_coord = vert_coord.lower() if vert_coord is not None else "none"
        
    valid_coords = ("pressure", "pres", "p", "ght_msl", 
                    "ght_agl", "theta", "th", "theta-e", "thetae", "eth")
    
    valid_field_types = ("none", "pressure", "pres", "p", "z",
                         "tc", "tk", "theta", "th", "theta-e", "thetae", 
                         "eth", "ght")
    
    icase_lookup = {"none" : 0,
                    "p" : 1,
                    "pres" : 1,
                    "pressure" : 1,
                    "z" : 2,
                    "ght" : 2,
                    "tc" : 3, 
                    "tk" : 4,
                    "theta" : 5,
                    "th" : 5,
                    "theta-e" : 6,
                    "thetae" : 6,
                    "eth" : 6}
    
    # These constants match what's in the fortran code.  
    rgas = Constants.RD
    ussalr = Constants.USSALR
    sclht = Constants.SCLHT
    
    # interp_levels might be a list or tuple, make a numpy array
    if not isinstance(interp_levels, np.ndarray):
        interp_levels = np.asarray(interp_levels, np.float64)
        
    # TODO: Check if field is staggered
    if is_staggered(wrfin, field):
        raise RuntimeError("Please unstagger field in the vertical")
    
    # Check for valid coord
    if vert_coord not in valid_coords:
        raise RuntimeError("'%s' is not a valid vertical "
                           "coordinate type" %  vert_coord)
    
    # Check for valid field type
    if field_type not in valid_field_types:
        raise RuntimeError("'%s' is not a valid field type" % field_type)
    
    log_p_int = 1 if log_p else 0
    
    icase = 0
    extrap = 0
    
    if extrapolate:
        extrap = 1
        icase = icase_lookup[field_type]
    
    # Extract vriables
    #timeidx = -1 # Should this be an argument?
    ncvars = extract_vars(wrfin, timeidx, ("PSFC", "QVAPOR", "F"), 
                          method, squeeze, cache, meta=False, _key=_key)
    
    sfp = ncvars["PSFC"] * ConversionFactors.PA_TO_HPA
    qv = ncvars["QVAPOR"]
    coriolis = ncvars["F"]
    
    terht = get_terrain(wrfin, timeidx, units="m", 
                        method=method, squeeze=squeeze, cache=cache,
                        meta=False, _key=_key)
    t = get_theta(wrfin, timeidx,  units="k", 
                  method=method, squeeze=squeeze, cache=cache, 
                  meta=False, _key=_key)
    tk = get_temp(wrfin, timeidx, units="k",  
                  method=method, squeeze=squeeze, cache=cache, 
                  meta=False, _key=_key)
    p = get_pressure(wrfin, timeidx, units="pa",  
                     method=method, squeeze=squeeze, cache=cache,
                     meta=False, _key=_key)
    ght = get_height(wrfin, timeidx, msl=True, units="m", 
                     method=method, squeeze=squeeze, cache=cache,
                     meta=False, _key=_key)
    ht_agl = get_height(wrfin, timeidx, msl=False, units="m",
                        method=method, squeeze=squeeze, cache=cache,
                        meta=False, _key=_key)
    
    smsfp = _smooth2d(sfp, 3)        

    vcor = 0
    
    if vert_coord in ("pressure", "pres", "p"):
        vcor = 1
        vcord_array = p * ConversionFactors.PA_TO_HPA
        
    elif vert_coord == "ght_msl":
        vcor = 2
        vcord_array = np.exp(-ght/sclht)
        
    elif vert_coord == "ght_agl":
        vcor = 3
        vcord_array = np.exp(-ht_agl/sclht)
    
    elif vert_coord in ("theta", "th"):
        vcor = 4
        idir = 1
        icorsw = 0
        delta = 0.01
        
        p_hpa = p * ConversionFactors.PA_TO_HPA
        
        vcord_array = _monotonic(t, p_hpa, coriolis, idir, delta, icorsw)
        
        # We only extrapolate temperature fields below ground 
        # if we are interpolating to pressure or height vertical surfaces.
        
        icase = 0 
        
    elif vert_coord in ("theta-e", "thetae", "eth"):
        vcor = 5
        icorsw = 0
        idir = 1
        delta = 0.01
        
        eth = get_eth(wrfin, timeidx, method=method, squeeze=squeeze, 
            cache=cache, meta=False, _key=_key)
        
        p_hpa = p * ConversionFactors.PA_TO_HPA
        
        vcord_array = _monotonic(eth, p_hpa, coriolis, idir, delta, icorsw)
        # We only extrapolate temperature fields below ground if we are
        # interpolating to pressure or height vertical surfaces
        icase = 0
    
    # Set the missing value
    if isinstance(field, ma.MaskedArray):
        missing = field.fill_value
    else:
        missing = Constants.DEFAULT_FILL
    
    if (field.shape != p.shape):
        raise ValueError("'field' shape does not match other variable shapes. "
                         "Verify that the 'timeidx' parameter matches the "
                         "same value used when extracting the 'field' "
                         "variable.")
            
    res = _vintrp(field, p, tk, qv, ght, terht, sfp, smsfp,
                  vcord_array, interp_levels,
                  icase, extrap, vcor, log_p_int, missing)
    
    return ma.masked_values(res, missing)


    
    
    
    
