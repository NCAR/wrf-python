from __future__ import (absolute_import, division, print_function)

import numpy as np
import numpy.ma as ma

from .extension import (_interpz3d, _vertcross, _interpline, _smooth2d,
                        _monotonic, _vintrp, _interpz3d_lev2d)

from .metadecorators import set_interp_metadata
from .util import (extract_vars, is_staggered, get_id, to_np, get_iterable,
                   is_moving_domain, is_latlon_pair)
from .py3compat import py3range
from .interputils import get_xy, get_xy_z_params, to_xy_coords
from .constants import Constants, default_fill, ConversionFactors
from wrf.g_terrain import get_terrain
from wrf.g_geoht import get_height
from wrf.g_temp import get_theta, get_temp, get_eth
from wrf.g_pressure import get_pressure


#  Note:  Extension decorator is good enough to handle left dims
@set_interp_metadata("horiz")
def interplevel(field3d, vert, desiredlev, missing=default_fill(np.float64),
                squeeze=True, meta=True):
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

        desiredlev (:obj:`float`, 1D sequence, or :class:`numpy.ndarray`): The
            desired vertical level(s). This can be a single value (e.g. 500),
            a sequence of values (e.g. [1000, 850, 700, 500, 250]), or a
            multidimensional array where the right two dimensions (ny x nx)
            must match *field3d*, and any leftmost dimensions match
            field3d.shape[:-3] (e.g. planetary boundary layer).
            Must be in the same units as the *vert* parameter.

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(numpy.float64)`.

        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

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


        Interpolate Relative Humidity to Boundary Layer Heights

        .. code-block:: python

            from netCDF4 import Dataset
            from wrf import getvar, interplevel

            wrfin = Dataset("wrfout_d02_2010-06-13_21:00:00")

            rh = getvar(wrfin, "rh")
            height = getvar(wrfin, "height_agl")
            pblh = getvar(wrfin, "PBLH")

            rh_pblh = interplevel(rh, height, pblh)


    """

    _desiredlev = np.asarray(desiredlev)
    if _desiredlev.ndim == 0:
        _desiredlev = np.array([desiredlev], np.float64)
        levsare2d = False
    else:
        levsare2d = _desiredlev.ndim >= 2

    if not levsare2d:
        result = _interpz3d(field3d, vert, _desiredlev, missing)
    else:
        result = _interpz3d_lev2d(field3d, vert, _desiredlev, missing)

    masked = ma.masked_values(result, missing)

    if not meta:
        if squeeze:
            return masked.squeeze()

    return masked


@set_interp_metadata("cross")
def vertcross(field3d, vert, levels=None, missing=default_fill(np.float64),
              wrfin=None, timeidx=0, stagger=None, projection=None,
              ll_point=None,
              pivot_point=None, angle=None,
              start_point=None, end_point=None,
              latlon=False, autolevels=100, cache=None, meta=True):
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
            Default is :data:`wrf.default_fill(numpy.float64)`.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. This is used
            to obtain the map projection when using latitude,longitude
            coordinates. Default is None.

        timeidx (:obj:`int`, optional): The
            desired time index when obtaining map boundary information
            from moving nests. This value can be a positive or negative
            integer. Only required when *wrfin* is specified and the nest is
            moving. Currently, :data:`wrf.ALL_TIMES` is not supported.
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

        ll_point (:class:`wrf.CoordPair`, sequence of :class:`wrf.CoordPair`, \
        optional): The lower left latitude, longitude point for your domain,
            and must be specified
            if *wrfin* is None. If the domain is a moving nest, this should be
            a sequence of :class:`wrf.CoordPair`. Default is None.

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
            (if enabled). This can be helpful for plotting.  Default is False.

            Note:

                Currently, *field3d* must be of type :class:`xarray.DataArray`
                and contain coordinate information in order to generate the
                latitude and longitude coordinates along the line if
                *latlon* is set to True. Otherwise, a warning will be issued,
                and the latitude and longitude information will not be
                present.

        autolevels(:obj:`int`, optional): The number of evenly spaced
            automatically chosen vertical levels to use when *levels*
            is None. Default is 100.

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

        if (latlon is True or is_latlon_pair(start_point) or
                is_latlon_pair(pivot_point)):

            if wrfin is not None:
                is_moving = is_moving_domain(wrfin)
            else:
                is_moving = False

            if timeidx is None:
                if wrfin is not None:
                    # Moving nests aren't supported with ALL_TIMES because the
                    # domain could move outside of the line, which causes
                    # crashes or different line lengths.
                    if is_moving:
                        raise ValueError("Requesting all times with a moving "
                                         "nest is not supported when using "
                                         "lat/lon cross sections because the "
                                         "domain could move outside of the "
                                         "cross section. You must request "
                                         "each time individually.")
                    else:
                        # Domain not moving, just use 0
                        _timeidx = 0

                # If using grid coordinates, then don't care about lat/lon
                # coordinates. Just use 0.
                else:
                    _timeidx = 0
            else:
                if is_moving:
                    _timeidx = timeidx
                else:
                    # When using non-moving nests, set the time to 0
                    # to avoid problems downstream
                    _timeidx = 0

        if pivot_point is not None:
            if pivot_point.lat is not None and pivot_point.lon is not None:
                xy_coords = to_xy_coords(pivot_point, wrfin, _timeidx,
                                         stagger, projection, ll_point)
                pivot_point_xy = (xy_coords.x, xy_coords.y)
            else:
                pivot_point_xy = (pivot_point.x, pivot_point.y)

        if start_point is not None and end_point is not None:
            if start_point.lat is not None and start_point.lon is not None:
                xy_coords = to_xy_coords(start_point, wrfin, _timeidx,
                                         stagger, projection, ll_point)
                start_point_xy = (xy_coords.x, xy_coords.y)
            else:
                start_point_xy = (start_point.x, start_point.y)

            if end_point.lat is not None and end_point.lon is not None:
                xy_coords = to_xy_coords(end_point, wrfin, _timeidx,
                                         stagger, projection, ll_point)
                end_point_xy = (xy_coords.x, xy_coords.y)
            else:
                end_point_xy = (end_point.x, end_point.y)

        xy, var2dz, z_var2d = get_xy_z_params(to_np(vert), pivot_point_xy,
                                              angle, start_point_xy,
                                              end_point_xy, levels,
                                              autolevels)

    if not multi:
        result = _vertcross(field3d, xy, var2dz, z_var2d, missing)
    else:
        outshape = field3d.shape[0:-3] + (z_var2d.shape[0], xy.shape[0])
        result = np.empty(outshape, dtype=field3d.dtype)

        for i in py3range(field3d.shape[0]):
            result[i, :] = _vertcross(field3d[i, :], xy, var2dz, z_var2d,
                                      missing)[:]

    return ma.masked_values(result, missing)


@set_interp_metadata("line")
def interpline(field2d, wrfin=None, timeidx=0, stagger=None, projection=None,
               ll_point=None,
               pivot_point=None, angle=None, start_point=None,
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
            from moving nests. This value can be a positive or negative
            integer. Only required when *wrfin* is specified and the nest is
            moving. Currently, :data:`wrf.ALL_TIMES` is not supported.
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

        ll_point (:class:`wrf.CoordPair`, sequence of :class:`wrf.CoordPair`, \
        optional): The lower left latitude, longitude point for your domain,
            and must be specified
            if *wrfin* is None. If the domain is a moving nest, this should be
            a sequence of :class:`wrf.CoordPair`. Default is None.

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

            Note:

                Currently, *field2d* must be of type :class:`xarray.DataArray`
                and contain coordinate information in order to generate the
                latitude and longitude coordinates along the line if
                *latlon* is set to True. Otherwise, a warning will be issued,
                and the latitude and longitude information will not be
                present.

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

    try:
        xy = cache["xy"]
    except (KeyError, TypeError):
        start_point_xy = None
        end_point_xy = None
        pivot_point_xy = None

        if (latlon is True or is_latlon_pair(start_point) or
                is_latlon_pair(pivot_point)):

            if wrfin is not None:
                is_moving = is_moving_domain(wrfin)
            else:
                is_moving = False

            if timeidx is None:
                if wrfin is not None:
                    # Moving nests aren't supported with ALL_TIMES because the
                    # domain could move outside of the line, which causes
                    # crashes or different line lengths.
                    if is_moving:
                        raise ValueError("Requesting all times with a moving "
                                         "nest is not supported when using a "
                                         "lat/lon line because the domain "
                                         "could move outside of line. "
                                         "You must request each time "
                                         "individually.")
                    else:
                        # Domain not moving, just use 0
                        _timeidx = 0

                # If using grid coordinates, then don't care about lat/lon
                # coordinates. Just use 0.
                else:
                    _timeidx = 0
            else:
                if is_moving:
                    _timeidx = timeidx
                else:
                    # When using non-moving nests, set the time to 0
                    # to avoid problems downstream
                    _timeidx = 0

        if pivot_point is not None:
            if pivot_point.lat is not None and pivot_point.lon is not None:
                xy_coords = to_xy_coords(pivot_point, wrfin, _timeidx,
                                         stagger, projection, ll_point)
                pivot_point_xy = (xy_coords.x, xy_coords.y)
            else:
                pivot_point_xy = (pivot_point.x, pivot_point.y)

        if start_point is not None and end_point is not None:
            if start_point.lat is not None and start_point.lon is not None:
                xy_coords = to_xy_coords(start_point, wrfin, _timeidx,
                                         stagger, projection, ll_point)
                start_point_xy = (xy_coords.x, xy_coords.y)
            else:
                start_point_xy = (start_point.x, start_point.y)

            if end_point.lat is not None and end_point.lon is not None:
                xy_coords = to_xy_coords(end_point, wrfin, _timeidx,
                                         stagger, projection, ll_point)
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
                * 'theta-e', 'thetae', 'eth': equivalent potential \
                temperature [K]

        interp_levels (sequence): A 1D sequence of vertical levels to
            interpolate to. Values must be in the same units as specified
            above for the *vert_coord* parameter.

        extrapolate (:obj:`bool`, optional): Set to True to extrapolate
            values below ground.  This is only performed when *vert_coord* is
            a pressure or height type, and the *field_type* is a pressure type
            (with height vertical coordinate), a height type (with pressure as
            the vertical coordinate), or a temperature type (with either height
            or pressure as the vertical coordinate). If those conditions are
            not met, or *field_type* is None, then the lowest model level
            will be used. Extrapolation is performed using standard atmosphere.
            Default is False.

        field_type (:obj:`str`, optional):
            The type of field.  Default is None.

            Valid strings are:
                * 'none': None
                * 'pressure', 'pres', 'p': pressure [Pa]
                * 'pressure_hpa', 'pres_hpa', 'p_hpa': pressure [hPa]
                * 'z', 'ght': geopotential height [m]
                * 'z_km', 'ght_km': geopotential height [km]
                * 'tc': temperature [degC]
                * 'tk': temperature [K]
                * 'theta', 'th': potential temperature [K]
                * 'theta-e', 'thetae', 'eth': equivalent potential temperature

        log_p (:obj:`bool`, optional): Set to True to use the log of the
            vertical coordinate for interpolation. This is mainly intended
            for pressure vertical coordinate types, but note that the log
            will still be taken for any vertical coordinate type when
            this is set to True. Default is False.

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

    _wrfin = get_iterable(wrfin)

    # Remove case sensitivity
    field_type = field_type.lower() if field_type is not None else "none"
    vert_coord = vert_coord.lower() if vert_coord is not None else "none"

    valid_coords = ("pressure", "pres", "p", "ght_msl",
                    "ght_agl", "theta", "th", "theta-e", "thetae", "eth")

    valid_field_types = ("none", "pressure", "pres", "p",
                         'pressure_hpa', 'pres_hpa', 'p_hpa', "z",
                         "tc", "tk", "theta", "th", "theta-e", "thetae",
                         "eth", "ght", 'z_km', 'ght_km')

    icase_lookup = {"none": 0,
                    "p": 1,
                    "pres": 1,
                    "pressure": 1,
                    "p_hpa": 1,
                    "pres_hpa": 1,
                    "pressure_hpa": 1,
                    "z": 2,
                    "ght": 2,
                    "z_km": 2,
                    "ght_km": 2,
                    "tc": 3,
                    "tk": 4,
                    "theta": 5,
                    "th": 5,
                    "theta-e": 6,
                    "thetae": 6,
                    "eth": 6}

    in_unitmap = {"p_hpa": 1.0/ConversionFactors.PA_TO_HPA,
                  "pres_hpa": 1.0/ConversionFactors.PA_TO_HPA,
                  "pressure_hpa": 1.0/ConversionFactors.PA_TO_HPA,
                  "z_km": 1.0/ConversionFactors.M_TO_KM,
                  "ght_km": 1.0/ConversionFactors.M_TO_KM,
                  }

    out_unitmap = {"p_hpa": ConversionFactors.PA_TO_HPA,
                   "pres_hpa": ConversionFactors.PA_TO_HPA,
                   "pressure_hpa": ConversionFactors.PA_TO_HPA,
                   "z_km": ConversionFactors.M_TO_KM,
                   "ght_km": ConversionFactors.M_TO_KM,
                   }

    # These constants match what's in the fortran code.
    rgas = Constants.RD
    ussalr = Constants.USSALR
    sclht = Constants.SCLHT

    # interp_levels might be a list or tuple, make a numpy array
    if not isinstance(interp_levels, np.ndarray):
        interp_levels = np.asarray(interp_levels, np.float64)

    if len(interp_levels) == 0:
        raise ValueError("'interp_levels' contains no values")

    # Check if field is staggered
    if is_staggered(_wrfin, field):
        raise ValueError("Please unstagger field in the vertical")

    # Check for valid coord
    if vert_coord not in valid_coords:
        raise ValueError("'{}' is not a valid vertical "
                         "coordinate type".format(vert_coord))

    # Check for valid field type
    if field_type not in valid_field_types:
        raise ValueError("'{}' is not a valid field type".format(field_type))

    log_p_int = 1 if log_p else 0

    icase = 0
    extrap = 0

    if extrapolate:
        extrap = 1
        icase = icase_lookup[field_type]

    # Extract variables
    ncvars = extract_vars(_wrfin, timeidx, ("PSFC", "QVAPOR"),
                          method, squeeze, cache, meta=False, _key=_key)

    sfp = ncvars["PSFC"] * ConversionFactors.PA_TO_HPA
    qv = ncvars["QVAPOR"]

    terht = get_terrain(_wrfin, timeidx, units="m",
                        method=method, squeeze=squeeze, cache=cache,
                        meta=False, _key=_key)
    tk = get_temp(_wrfin, timeidx, units="k",
                  method=method, squeeze=squeeze, cache=cache,
                  meta=False, _key=_key)
    p = get_pressure(_wrfin, timeidx, units="pa",
                     method=method, squeeze=squeeze, cache=cache,
                     meta=False, _key=_key)
    ght = get_height(_wrfin, timeidx, msl=True, units="m",
                     method=method, squeeze=squeeze, cache=cache,
                     meta=False, _key=_key)

    smsfp = _smooth2d(sfp, 3, 2.0)

    vcor = 0

    if vert_coord in ("pressure", "pres", "p"):
        vcor = 1
        vcord_array = p * ConversionFactors.PA_TO_HPA

    elif vert_coord == "ght_msl":
        vcor = 2
        vcord_array = np.exp(-ght/sclht)

    elif vert_coord == "ght_agl":
        ht_agl = get_height(_wrfin, timeidx, msl=False, units="m",
                            method=method, squeeze=squeeze, cache=cache,
                            meta=False, _key=_key)

        vcor = 3
        vcord_array = np.exp(-ht_agl/sclht)

    elif vert_coord in ("theta", "th"):
        t = get_theta(_wrfin, timeidx,  units="k",
                      method=method, squeeze=squeeze, cache=cache,
                      meta=False, _key=_key)

        coriolis = extract_vars(_wrfin, timeidx, "F",
                                method, squeeze, cache, meta=False,
                                _key=_key)["F"]

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

        eth = get_eth(_wrfin, timeidx, method=method, squeeze=squeeze,
                      cache=cache, meta=False, _key=_key)

        coriolis = extract_vars(_wrfin, timeidx, "F",
                                method, squeeze, cache, meta=False,
                                _key=_key)["F"]

        p_hpa = p * ConversionFactors.PA_TO_HPA

        vcord_array = _monotonic(eth, p_hpa, coriolis, idir, delta, icorsw)
        # We only extrapolate temperature fields below ground if we are
        # interpolating to pressure or height vertical surfaces
        icase = 0

    # Set the missing value
    if isinstance(field, ma.MaskedArray):
        missing = field.fill_value
    else:
        missing = default_fill(np.float64)

    if (field.shape != p.shape):
        raise ValueError("'field' shape does not match other variable shapes. "
                         "Verify that the 'timeidx' parameter matches the "
                         "same value used when extracting the 'field' "
                         "variable.")

    # Some field types are in different units than the Fortran routine
    # expects

    conv_factor = in_unitmap.get(field_type)

    if conv_factor is not None:
        field_ = field * conv_factor
    else:
        field_ = field

    res = _vintrp(field_, p, tk, qv, ght, terht, sfp, smsfp,
                  vcord_array, interp_levels,
                  icase, extrap, vcor, log_p_int, missing)

    conv_factor = out_unitmap.get(field_type)

    if conv_factor is not None:
        res_ = res * conv_factor
    else:
        res_ = res

    return ma.masked_values(res_, missing)
