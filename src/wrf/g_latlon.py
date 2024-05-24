from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict

import numpy as np

from .util import extract_vars, get_id, get_iterable, is_mapping, to_np
from .py3compat import viewkeys
from .latlonutils import _lat_varname, _lon_varname, _ll_to_xy, _xy_to_ll
from .metadecorators import set_latlon_metadata
from .constants import Constants, ProjectionTypes
from .config import xarray_enabled

if xarray_enabled():
    from xarray import DataArray


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


def _llxy_mapping(wrfin, x_or_lat, y_or_lon, func, timeidx, stagger,
                  squeeze, meta, as_int=None):

    """Return the x,y/lat,lon coordinates for a dictionary input.

    The leftmost dimension(s) for the result is:

        - return_val[key,...,0,...] will contain the x/lat values.
        - return_val[key,...,1,...] will contain the y/lon values.

    Nested dictionaries are allowed.

    Args:

        wrfin (:obj:`dict`): A mapping of key name to a WRF NetCDF file object
            or sequence of WRF NetCDF file objects.

        x_or_lat (:obj:`float` or sequence): A single latitude/x value or a
            sequence of latitude/x values to be converted.

        y_or_lon (:obj:`float` or sequence): A single longitude/y value or a
            sequence of longitude/y values to be converted.

        func (function): Either the xy_to_ll or ll_to_xy function.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. The default is 0.

        stagger (:obj:`str`): By default, the values are returned on the mass
            grid, but a staggered grid can be chosen with the following
            options:

                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component,
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component,
                  which has a staggered south_north (y) dimension.

        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        as_int (:obj:`bool`, optional): Set to True to return the x,y values as
            :obj:`int`, otherwise they will be returned as :obj:`float`. This
            is only used when *func* is ll_to_xy.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        lat,lon/x,y coordinate value(s) whose leftmost dimensions are the
        dictionary keys, followed by a dimension of size
        2 (0=X, 1=Y)/(0=lat, 1=lon).
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    keynames = []
    # This might not work once mapping iterators are implemented
    numkeys = len(wrfin)

    key_iter = iter(viewkeys(wrfin))
    first_key = next(key_iter)
    keynames.append(first_key)

    first_args = [wrfin[first_key], x_or_lat, y_or_lon, timeidx, squeeze,
                  meta, stagger]
    if as_int is not None:
        first_args.append(as_int)

    first_array = func(*first_args)

    # Create the output data numpy array based on the first array
    outdims = [numkeys]
    outdims += first_array.shape
    outdata = np.empty(outdims, first_array.dtype)
    outdata[0, :] = first_array[:]

    idx = 1
    while True:
        try:
            key = next(key_iter)
        except StopIteration:
            break
        else:
            keynames.append(key)

            args = [wrfin[first_key], x_or_lat, y_or_lon, timeidx, squeeze,
                    meta, stagger]
            if as_int is not None:
                args.append(as_int)

            result_array = func(*args)

            if outdata.shape[1:] != result_array.shape:
                raise ValueError("data sequences must have the "
                                 "same size for all dictionary keys")
            outdata[idx, :] = to_np(result_array)[:]
            idx += 1

    if xarray_enabled() and meta:
        outname = str(first_array.name)
        # Note: assumes that all entries in dict have same coords
        outcoords = OrderedDict(first_array.coords)

        # First find and store all the existing key coord names/values
        # This is applicable only if there are nested dictionaries.
        key_coordnames = []
        coord_vals = []
        existing_cnt = 0
        while True:
            key_coord_name = "key_{}".format(existing_cnt)

            if key_coord_name not in first_array.dims:
                break

            key_coordnames.append(key_coord_name)
            coord_vals.append(to_np(first_array.coords[key_coord_name]))

            existing_cnt += 1

        # Now add the key coord name and values for THIS dictionary.
        # Put the new key_n name at the bottom, but the new values will
        # be at the top to be associated with key_0 (left most).  This
        # effectively shifts the existing 'key_n' coordinate values to the
        # right one dimension so *this* dicionary's key coordinate values
        # are at 'key_0'.
        key_coordnames.append(key_coord_name)
        coord_vals.insert(0, keynames)

        # make it so that key_0 is leftmost
        outdims = key_coordnames + list(first_array.dims[existing_cnt:])

        # Create the new 'key_n', value pairs
        for coordname, coordval in zip(key_coordnames, coord_vals):
            outcoords[coordname] = coordval

        outattrs = OrderedDict(first_array.attrs)

        outarr = DataArray(outdata, name=outname, coords=outcoords,
                           dims=outdims, attrs=outattrs)
    else:
        outarr = outdata

    return outarr


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

        as_int (:obj:`bool`): Set to False to return the x,y values as
            :obj:`float`, otherwise they will be returned as :obj:`int`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        x,y coordinate value(s) whose leftmost dimension is 2 (0=X, 1=Y).
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    if is_mapping(wrfin):
        return _llxy_mapping(wrfin, latitude, longitude, ll_to_xy,
                             timeidx, stagger, squeeze, meta, as_int)
    _key = get_id(wrfin)
    _wrfin = get_iterable(wrfin)

    return _ll_to_xy(latitude, longitude, _wrfin, timeidx, stagger, "cat",
                     squeeze, None, _key, as_int, **{})


def _set_defaults(projparams):
    """Check projection parameters and set defaults.

    Throws an exception if projection parameters required by WPS are not
    provided, along with any other parameters required by the map projection
    transformation routines.

    For parameters not used by the projection type, defaults are used so that
    the None values don't pass through to fortran downstream.

    Args:

        projparams (:obj:`dict`): Projection parameters dictionary.

    Returns:
        :obj:`dict`: The projection parameters with default values inserted
        where applicable.

    """
    _params = dict(projparams)

    map_proj = _params.get("map_proj")
    # All projections require these arguments
    if map_proj is None:
        raise ValueError("'map_proj' cannot be None")

    if _params.get("ref_lat") is None:
        raise ValueError("'ref_lat' cannot be None")

    if _params.get("ref_lon") is None:
        raise ValueError("'ref_lon' cannot be None")

    if _params.get("known_x") is None:
        raise ValueError("'known_x' cannot be None")

    if _params.get("known_y") is None:
        raise ValueError("'known_y' cannot be None")

    if _params.get("dx") is None:
        raise ValueError("'dx' cannot be None")

    # Requires truelat1,stand_lon, truelat2, dx, dy
    if map_proj == ProjectionTypes.LAMBERT_CONFORMAL:
        if _params.get("truelat1") is None:
            raise ValueError("'truelat1' cannot be None")

        if _params.get("stand_lon") is None:
            raise ValueError("'stand_lon' cannot be None")

        if _params.get("truelat2") is None:
            _params["truelat2"] = _params["truelat1"]

    # Requires truelat1, stand_lon
    elif map_proj == ProjectionTypes.POLAR_STEREOGRAPHIC:
        if _params.get("truelat1") is None:
            raise ValueError("'truelat1' cannot be None")

        if _params.get("stand_lon") is None:
            raise ValueError("'stand_lon' cannot be None")

    # Requires truelat1
    elif map_proj == ProjectionTypes.MERCATOR:
        if _params.get("truelat1") is None:
            raise ValueError("'truelat1' cannot be None")

        if _params.get("stand_lon") is None:
            _params["stand_lon"] = 0.0

    # Requires pole_lat, pole_lon, stand_lon
    elif map_proj == ProjectionTypes.LAT_LON:
        if _params.get("stand_lon") is None:
            raise ValueError("'stand_lon' cannot be None")

        if _params.get("dy") is None:
            raise ValueError("'dy' cannot be None")

        if _params.get("pole_lat") is None:
            raise ValueError("'pole_lat' cannot be None")

        if _params.get("pole_lon") is None:
            raise ValueError("'pole_lon' cannot be None")

        if _params.get("latinc") is None:
            _params["latinc"] = ((_params["dy"]*360.0)/2.0 /
                                 Constants.PI/Constants.WRF_EARTH_RADIUS)

        if _params.get("loninc") is None:
            _params["loninc"] = ((_params["dx"]*360.0)/2.0 /
                                 Constants.PI/Constants.WRF_EARTH_RADIUS)

    else:
        raise ValueError("invalid 'map_proj' value of {}".format(map_proj))

    # Set these to defaults if not used so that the Fortran routines
    # don't crash
    if _params.get("truelat1") is None:
        _params["truelat1"] = 0.

    if _params.get("truelat2") is None:
        _params["truelat2"] = 0.

    if _params.get("pole_lat") is None:
        _params["pole_lat"] = 90.0

    if _params.get("pole_lon") is None:
        _params["pole_lon"] = 0.0

    if _params.get("dx") is None:
        _params["dx"] = 0.0

    if _params.get("dy") is None:
        _params["dy"] = 0.0

    if _params.get("latinc") is None:
        _params["latinc"] = 0.

    if _params.get("loninc") is None:
        _params["loninc"] = 0.

    return _params


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

        as_int (:obj:`bool`): Set to False to return the x,y values as
            :obj:`float`, otherwise they will be returned as :obj:`int`.
            Default is True.

        map_proj (:obj:`int`): Model projection [1=Lambert Conformal,
            2=Polar Stereographic, 3=Mercator, 6=Lat-Lon].  Required.

        truelat1 (:obj:`float`): Latitude of true scale 1.  Required for
            map_proj = 1, 2, 3 (defaults to 0 otherwise).

        truelat2 (:obj:`float`): Latitude of true scale 2.  Optional for
            map_proj = 1 (defaults to 0 otherwise).

        stand_lon (:obj:`float`): Standard longitude. Required for *map_proj* =
            1, 2, 6 (defaults to 0 otherwise).

        ref_lat (:obj:`float`): A reference latitude.  Required.

        ref_lon (:obj:`float`): A reference longitude.  Required.

        known_x (:obj:`float`): The known x-coordinate associated with
            *ref_lon*. Required.

        known_y (:obj:`float`): The known y-coordinate associated with
            *ref_lat*.  Required.

        pole_lat (:obj:`float`): Pole latitude. Required for
            *map_proj* = 6 (use 90 for no rotation).

        pole_lon (:obj:`float`): Pole longitude. Required for
            *map_proj* = 6 (use 0 for no rotation).

        dx (:obj:`float`): The x spacing in meters at the true latitude.
            Required for all map projections.

        dy (:obj:`float`) - The y spacing in meters at the true latitude.
            Required for *map_proj* = 6 (defaults to 0 otherwise).

        latinc (:obj:`float`): Optional for *map_proj* = 6. Default is:

            .. code-block:: python

                latinc = (dy*360.0)/2.0/Constants.PI/Constants.WRF_EARTH_RADIUS

        loninc (:obj:`float`): Optional for *map_proj* = 6. Default is:

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
    _projparams = {name: loc[name] for name in ("map_proj", "truelat1",
                                                "truelat2", "stand_lon",
                                                "ref_lat", "ref_lon",
                                                "pole_lat", "pole_lon",
                                                "known_x", "known_y", "dx",
                                                "dy", "latinc", "loninc")}

    projparams = _set_defaults(_projparams)

    return _ll_to_xy(latitude, longitude, None, 0, True, "cat", squeeze, None,
                     None, as_int, **projparams)


@set_latlon_metadata(xy=False)
def xy_to_ll(wrfin, x, y, timeidx=0, squeeze=True, meta=True, stagger=None):
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

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        latitude and longitude values whose leftmost dimension is 2
        (0=latitude, 1=longitude).
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    if is_mapping(wrfin):
        return _llxy_mapping(wrfin, x, y, xy_to_ll,
                             timeidx, stagger, squeeze, meta)

    _key = get_id(wrfin)
    _wrfin = get_iterable(wrfin)
    return _xy_to_ll(x, y, _wrfin, timeidx, stagger, "cat", True, None,
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

        map_proj (:obj:`int`): Model projection [1=Lambert Conformal,
            2=Polar Stereographic, 3=Mercator, 6=Lat-Lon].  Required.

        truelat1 (:obj:`float`): Latitude of true scale 1.  Required for
            map_proj = 1, 2, 3 (defaults to 0 otherwise).

        truelat2 (:obj:`float`): Latitude of true scale 2.  Optional for
            map_proj = 1 (defaults to 0 otherwise).

        stand_lon (:obj:`float`): Standard longitude. Required for *map_proj* =
            1, 2, 6 (defaults to 0 otherwise).

        ref_lat (:obj:`float`): A reference latitude.  Required.

        ref_lon (:obj:`float`): A reference longitude.  Required.

        known_x (:obj:`float`): The known x-coordinate associated with
            *ref_lon*. Required.

        known_y (:obj:`float`): The known y-coordinate associated with
            *ref_lat*.  Required.

        pole_lat (:obj:`float`): Pole latitude. Required for
            *map_proj* = 6 (use 90 for no rotation).

        pole_lon (:obj:`float`): Pole longitude. Required for
            *map_proj* = 6 (use 0 for no rotation).

        dx (:obj:`float`): The x spacing in meters at the true latitude.
            Required for all map projections.

        dy (:obj:`float`) - The y spacing in meters at the true latitude.
            Required for *map_proj* = 6 (defaults to 0 otherwise).

        latinc (:obj:`float`): Optional for *map_proj* = 6. Default is:

            .. code-block:: python

                latinc = (dy*360.0)/2.0/Constants.PI/Constants.WRF_EARTH_RADIUS

        loninc (:obj:`float`): Optional for *map_proj* = 6. Default is:

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
    _projparams = {name: loc[name] for name in ("map_proj", "truelat1",
                                                "truelat2", "stand_lon",
                                                "ref_lat", "ref_lon",
                                                "pole_lat", "pole_lon",
                                                "known_x", "known_y", "dx",
                                                "dy", "latinc", "loninc")}

    projparams = _set_defaults(_projparams)
    return _xy_to_ll(x, y, None, 0, None, "cat", squeeze, None, None,
                     **projparams)
