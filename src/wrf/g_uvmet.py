from __future__ import (absolute_import, division, print_function)

from math import fabs, log, tan, sin, cos

import numpy as np

from .extension import _uvmet
from .destag import destagger
from .constants import Constants
from .g_wind import _calc_wspd_wdir
from .decorators import convert_units
from .metadecorators import set_wind_metadata
from .util import extract_vars, extract_global_attrs, either


@convert_units("wind", "m s-1")
def _get_uvmet(wrfin, timeidx=0, method="cat", squeeze=True,
               cache=None, meta=True, _key=None,
               ten_m=False, units="m s-1"):
    """Return the u,v wind components rotated to earth coordinates.

    The leftmost dimension of the returned array represents two different
    quantities:

        - return_val[0,...] will contain U_EARTH
        - return_val[1,...] will contain V_EARTH

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

        ten_m (:obj:`bool`, optional): Set to True to use the 10m surface
            winds, rather than the three-dimensional wind field.  Default is
            False.

        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar`
            product table for a list of available units for 'uvmet'.  Default
            is 'm s-1'.

    Returns:

        :class:`numpy.ndarray`: The u,v components of the wind rotated to
        earth coordinates, whose leftmost dimensions is 2 (0=U, 1=V).  The
        *meta* parameter is ignored for this function, so only a
        :class:`numpy.ndarray` is returned.

    """

    if not ten_m:
        varname = either("U", "UU")(wrfin)
        u_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                              meta=False, _key=_key)

        u = destagger(u_vars[varname], -1)

        varname = either("V", "VV")(wrfin)
        v_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                              meta=False, _key=_key)
        v = destagger(v_vars[varname], -2)
    else:
        varname = either("U10", "UU")(wrfin)
        u_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                              meta=False, _key=_key)
        u = (u_vars[varname] if varname == "U10" else
             destagger(u_vars[varname][..., 0, :, :], -1))

        varname = either("V10", "VV")(wrfin)
        v_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                              meta=False, _key=_key)
        v = (v_vars[varname] if varname == "V10" else
             destagger(v_vars[varname][..., 0, :, :], -2))

    map_proj_attrs = extract_global_attrs(wrfin, attrs="MAP_PROJ")
    map_proj = map_proj_attrs["MAP_PROJ"]

    # 1 - Lambert
    # 2 - Polar Stereographic
    # 3 - Mercator
    # 6 - Lat/Lon
    # Note:  NCL has no code to handle other projections (0,4,5) as they
    # don't appear to be supported any longer

    if map_proj in (0, 3, 6):
        # No rotation needed for Mercator and Lat/Lon, but still need
        # u,v aggregated in to one array

        end_idx = -3 if not ten_m else -2
        resdim = (2,) + u.shape[0:end_idx] + u.shape[end_idx:]

        # Make a new output array for the result
        result = np.empty(resdim, u.dtype)

        # For 2D array, this makes (0,...,:,:) and (1,...,:,:)
        # For 3D array, this makes (0,...,:,:,:) and (1,...,:,:,:)
        idx0 = (0,) + (Ellipsis,) + (slice(None),)*(-end_idx)
        idx1 = (1,) + (Ellipsis,) + (slice(None),)*(-end_idx)

        try:
            fill = u.fill_value
        except AttributeError:
            result[idx0] = u[:]
            result[idx1] = v[:]
        else:
            result[idx0] = np.ma.filled(u[:], fill)
            result[idx1] = np.ma.filled(v[:], fill)
            result = np.ma.masked_values(result, fill)

        return result
    elif map_proj in (1, 2):
        lat_attrs = extract_global_attrs(wrfin, attrs=("TRUELAT1",
                                                       "TRUELAT2"))
        radians_per_degree = Constants.PI/180.0
        # Rotation needed for Lambert and Polar Stereographic
        true_lat1 = lat_attrs["TRUELAT1"]
        true_lat2 = lat_attrs["TRUELAT2"]

        try:
            lon_attrs = extract_global_attrs(wrfin, attrs="STAND_LON")
        except AttributeError:
            try:
                cen_lon_attrs = extract_global_attrs(wrfin, attrs="CEN_LON")
            except AttributeError:
                raise RuntimeError("longitude attributes not found in NetCDF")
            else:
                cen_lon = cen_lon_attrs["CEN_LON"]
        else:
            cen_lon = lon_attrs["STAND_LON"]

        varname = either("XLAT_M", "XLAT")(wrfin)
        xlat_var = extract_vars(wrfin, timeidx, varname,
                                method, squeeze, cache, meta=False,
                                _key=_key)
        lat = xlat_var[varname]

        varname = either("XLONG_M", "XLONG")(wrfin)
        xlon_var = extract_vars(wrfin, timeidx, varname,
                                method, squeeze, cache, meta=False,
                                _key=_key)
        lon = xlon_var[varname]

        if map_proj == 1:
            if((fabs(true_lat1 - true_lat2) > 0.1) and
                    (fabs(true_lat2 - 90.) > 0.1)):
                cone = (log(cos(true_lat1*radians_per_degree)) -
                        log(cos(true_lat2*radians_per_degree)))
                cone = (cone /
                        (log(tan((45.-fabs(true_lat1/2.))*radians_per_degree))
                         - log(tan((45.-fabs(true_lat2/2.)) *
                                   radians_per_degree))))
            else:
                cone = sin(fabs(true_lat1)*radians_per_degree)
        else:
            cone = 1

        result = _uvmet(u, v, lat, lon, cen_lon, cone)

        if squeeze:
            result = result.squeeze()

        return result


@set_wind_metadata(copy_varname=either("P", "PRES"),
                   name="uvmet",
                   description="earth rotated u,v",
                   two_d=False,
                   wspd_wdir=False)
def get_uvmet(wrfin, timeidx=0, method="cat", squeeze=True,
              cache=None, meta=True, _key=None,
              units="m s-1"):
    """Return the u,v wind components rotated to earth coordinates.

    The leftmost dimension of the returned array represents two different
    quantities:

        - return_val[0,...] will contain U_EARTH
        - return_val[1,...] will contain V_EARTH

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
            product table for a list of available units for 'uvmet'.  Default
            is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The u,v components
        of the wind rotated to earth coordinates, whose leftmost dimensions is
        2 (0=U_EARTH, 1=V_EARTH).
        If xarray is enabled and the *meta* parameter is True,
        then the result will be a :class:`xarray.DataArray` object.  Otherwise,
        the result will be a :class:`numpy.ndarray` object with no metadata.

    """

    return _get_uvmet(wrfin, timeidx, method, squeeze, cache, meta, _key,
                      False, units)


@set_wind_metadata(copy_varname=either("PSFC", "F"),
                   name="uvmet10",
                   description="10m earth rotated u,v",
                   two_d=True,
                   wspd_wdir=False)
def get_uvmet10(wrfin, timeidx=0, method="cat", squeeze=True,
                cache=None, meta=True, _key=None,
                units="m s-1"):
    """Return the u,v components for the 10m winds rotated to earth
    coordinates.

    The leftmost dimension of the returned array represents two different
    quantities:

        - return_val[0,...] will contain U10_EARTH
        - return_val[1,...] will contain V10_EARTH

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
            product table for a list of available units for 'uvmet10'.
            Default is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The u,v components
        of the 10m wind rotated to earth coordinates, whose leftmost dimensions
        is 2 (0=U10_EARTH, 1=V10_EARTH).
        If xarray is enabled and the *meta* parameter is
        True, then the result will be a :class:`xarray.DataArray` object.
        Otherwise, the result will be a :class:`numpy.ndarray` object with no
        metadata.

    """

    return _get_uvmet(wrfin, timeidx, method, squeeze, cache, meta, _key,
                      True, units)


@set_wind_metadata(copy_varname=either("P", "PRES"),
                   name="uvmet_wspd_wdir",
                   description="earth rotated wspd,wdir",
                   two_d=False,
                   wspd_wdir=True)
def get_uvmet_wspd_wdir(wrfin, timeidx=0, method="cat", squeeze=True,
                        cache=None, meta=True, _key=None,
                        units="m s-1"):
    """Return the wind speed and wind direction for the wind rotated to
    earth coordinates.

    This function should be used when comparing with observations.

    The leftmost dimension of the returned array represents two different
    quantities:

        - return_val[0,...] will contain WSPD_EARTH
        - return_val[1,...] will contain WDIR_EARTH

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
            product table for a list of available units for 'uvmet_wspd_wdir'.
            Default is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind speed and
        wind direction for the wind rotated to earth coordinates, whose
        leftmost dimensions is 2 (0=WSPD_EARTH, 1=WDIR_EARTH).  If
        xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    uvmet = _get_uvmet(wrfin, timeidx, method, squeeze,
                       cache, meta, _key, False, units="m s-1")

    return _calc_wspd_wdir(uvmet[0, ..., :, :, :], uvmet[1, ..., :, :, :],
                           False, units)


@set_wind_metadata(copy_varname=either("PSFC", "F"),
                   name="uvmet10_wspd_wdir",
                   description="10m earth rotated wspd,wdir",
                   two_d=True,
                   wspd_wdir=True)
def get_uvmet10_wspd_wdir(wrfin, timeidx=0, method="cat", squeeze=True,
                          cache=None, meta=True, _key=None,
                          units="m s-1"):
    """Return the wind speed and wind direction for the 10m wind rotated to
    earth coordinates.

    This function should be used when comparing with observations.

    The leftmost dimension of the returned array represents two different
    quantities:

        - return_val[0,...] will contain WSPD10_EARTH
        - return_val[1,...] will contain WDIR10_EARTH

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
            'uvmet10_wspd_wdir'. Default is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind speed and
        wind direction for the 10m wind rotated to earth coordinates, whose
        leftmost dimensions is 2 (0=WSPD10_EARTH, 1=WDIR10_EARTH).  If
        xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    uvmet10 = _get_uvmet(wrfin, timeidx, method, squeeze, cache, meta, _key,
                         True, units="m s-1")

    return _calc_wspd_wdir(uvmet10[0, ..., :, :], uvmet10[1, ..., :, :],
                           True, units)


def get_uvmet_wspd(wrfin, timeidx=0, method="cat", squeeze=True,
                   cache=None, meta=True, _key=None,
                   units="m s-1"):
    """Return the wind speed for the wind rotated to earth coordinates.

    This function should be used when comparing with observations.

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
            product table for a list of available units for 'uvmet_wspd_wdir'.
            Default is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind speed
        for the wind rotated to earth coordinates.  If
        xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_uvmet_wspd_wdir(wrfin, timeidx, method, squeeze,
                                 cache, meta, _key, units)[0, :]

    if meta:
        result.attrs["description"] = "earth rotated wspd"

    return result


def get_uvmet_wdir(wrfin, timeidx=0, method="cat", squeeze=True,
                   cache=None, meta=True, _key=None,
                   units="m s-1"):
    """Return the wind direction for the wind rotated to earth coordinates.

    This function should be used when comparing with observations.

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
            product table for a list of available units for 'uvmet_wspd_wdir'.
            Default is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind direction
        for the wind rotated to earth coordinates.  If
        xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_uvmet_wspd_wdir(wrfin, timeidx, method, squeeze,
                                 cache, meta, _key, units)[1, :]

    if meta:
        result.attrs["description"] = "earth rotated wdir"

    return result


def get_uvmet10_wspd(wrfin, timeidx=0, method="cat", squeeze=True,
                     cache=None, meta=True, _key=None,
                     units="m s-1"):
    """Return the wind speed for the 10m wind rotated to earth coordinates.

    This function should be used when comparing with observations.

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
            product table for a list of available units for 'uvmet_wspd_wdir'.
            Default is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind speed
        for the 10m wind rotated to earth coordinates.  If
        xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_uvmet10_wspd_wdir(wrfin, timeidx, method, squeeze,
                                   cache, meta, _key, units)[0, :]
    if meta:
        result.attrs["description"] = "10m earth rotated wspd"

    return result


def get_uvmet10_wdir(wrfin, timeidx=0, method="cat", squeeze=True,
                     cache=None, meta=True, _key=None,
                     units="m s-1"):
    """Return the wind direction for the 10m wind rotated to earth coordinates.

    This function should be used when comparing with observations.

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
            product table for a list of available units for 'uvmet_wspd_wdir'.
            Default is 'm s-1'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wind direction
        for the 10m wind rotated to earth coordinates.  If
        xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_uvmet10_wspd_wdir(wrfin, timeidx, method, squeeze,
                                   cache, meta, _key, units)[1, :]

    if meta:
        result.attrs["description"] = "10m earth rotated wdir"

    return result
