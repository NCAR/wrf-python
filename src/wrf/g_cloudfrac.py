from __future__ import (absolute_import, division, print_function)

import numpy as np
import numpy.ma as ma

from .constants import Constants, default_fill
from .extension import _tk, _rh, _cloudfrac
from .metadecorators import set_cloudfrac_metadata
from .util import extract_vars
from .g_geoht import _get_geoht


@set_cloudfrac_metadata()
def get_cloudfrac(wrfin, timeidx=0, method="cat", squeeze=True,
                  cache=None, meta=True, _key=None,
                  vert_type="height_agl", low_thresh=None, mid_thresh=None,
                  high_thresh=None, missing=default_fill(np.float64)):
    """Return the cloud fraction for low, mid, and high level clouds.

    The leftmost dimension of the returned array represents three different
    quantities:

        - return_val[0,...] will contain LOW level cloud fraction
        - return_val[1,...] will contain MID level cloud fraction
        - return_val[2,...] will contain HIGH level cloud fraction

    If the vertical coordinate type is 'height_agl' or 'height_msl', the
    default cloud levels are defined as:

        300 m <= low_cloud < 2000 m
        2000 m <= mid_cloud < 6000 m
        6000 m <= high_cloud

    For 'pressure', the default cloud levels are defined as:

        97000 Pa <= low_cloud < 80000 Pa
        80000 Pa <= mid_cloud < 45000 Pa
        45000 Pa <= high_cloud

    Note that the default low cloud levels are chosen to
    exclude clouds near the surface (fog). If you want fog included, set
    *low_thresh* to ~99500 Pa if *vert_type* is set to 'pressure', or 15 m if
    using 'height_msl' or 'height_agl'. Keep in mind that the lowest mass grid
    points are slightly above the ground, and in order to find clouds, the
    *low_thresh* needs to be set to values that are slightly greater than
    (less than) the lowest height (pressure) values.

    When using 'pressure' or 'height_agl' for *vert_type*, there is a
    possibility that the lowest WRF level will be higher than the low_cloud or
    mid_cloud threshold, particularly for mountainous regions.  When this
    happens, a fill value will be used in the output.

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

        vert_type (:obj:`str`, optional):  The type of vertical coordinate used
            to determine cloud type thresholds.  Must be 'height_agl',
            'height_msl', or 'pres'. The default is 'height_agl'.

        low_thresh (:obj:`float`, optional): The lower bound for what is
            considered a low cloud.  If *vert_type* is 'pres', the default is
            97000 Pa.  If *vert_type* is 'height_agl' or 'height_msl', then the
            default is 300 m.

        mid_thresh (:obj:`float`, optional): The lower bound for what is
            considered a mid level cloud.  If *vert_type* is 'pres', the
            default is 80000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 2000 m.

        high_thresh (:obj:`float`, optional): The lower bound for what is
            considered a high level cloud.  If *vert_type* is 'pres', the
            default is 45000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 6000 m.

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

    if vert_type.lower() == "pres" or vert_type.lower() == "pressure":
        v_coord = full_p
        _low_thresh = 97000. if low_thresh is None else low_thresh
        _mid_thresh = 80000. if mid_thresh is None else mid_thresh
        _high_thresh = 45000. if high_thresh is None else high_thresh
        vert_inc_w_height = 0
    elif (vert_type.lower() == "height_msl"
          or vert_type.lower() == "height_agl"):
        is_msl = vert_type.lower() == "height_msl"
        v_coord = _get_geoht(wrfin, timeidx, method, squeeze,
                             cache, meta=False, _key=_key, height=True,
                             msl=is_msl)
        _low_thresh = 300. if low_thresh is None else low_thresh
        _mid_thresh = 2000. if mid_thresh is None else mid_thresh
        _high_thresh = 6000. if high_thresh is None else high_thresh
        vert_inc_w_height = 1
    else:
        raise ValueError("'vert_type' must be 'pres', 'height_msl', "
                         "or 'height_agl'")

    cfrac = _cloudfrac(v_coord, rh, vert_inc_w_height,
                       _low_thresh, _mid_thresh, _high_thresh, missing)

    return ma.masked_values(cfrac, missing)


def get_low_cloudfrac(wrfin, timeidx=0, method="cat", squeeze=True,
                      cache=None, meta=True, _key=None,
                      vert_type="height_agl", low_thresh=None,
                      mid_thresh=None, high_thresh=None,
                      missing=default_fill(np.float64)):
    """Return the cloud fraction for the low level clouds.

    If the vertical coordinate type is 'height_agl' or 'height_msl', the
    default cloud levels are defined as:

        300 m <= low_cloud < 2000 m
        2000 m <= mid_cloud < 6000 m
        6000 m <= high_cloud

    For 'pressure', the default cloud levels are defined as:

        97000 Pa <= low_cloud < 80000 Pa
        80000 Pa <= mid_cloud < 45000 Pa
        45000 Pa <= high_cloud

    Note that the default low cloud levels are chosen to
    exclude clouds near the surface (fog). If you want fog included, set
    *low_thresh* to ~99500 Pa if *vert_type* is set to 'pressure', or 15 m if
    using 'height_msl' or 'height_agl'. Keep in mind that the lowest mass grid
    points are slightly above the ground, and in order to find clouds, the
    *low_thresh* needs to be set to values that are slightly greater than
    (less than) the lowest height (pressure) values.

    When using 'pressure' or 'height_agl' for *vert_type*, there is a
    possibility that the lowest WRF level will be higher than the low_cloud or
    mid_cloud threshold, particularly for mountainous regions.  When this
    happens, a fill value will be used in the output.

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

        vert_type (:obj:`str`, optional):  The type of vertical coordinate used
            to determine cloud type thresholds.  Must be 'height_agl',
            'height_msl', or 'pres'. The default is 'height_agl'.

        low_thresh (:obj:`float`, optional): The lower bound for what is
            considered a low cloud.  If *vert_type* is 'pres', the default is
            97000 Pa.  If *vert_type* is 'height_agl' or 'height_msl', then the
            default is 300 m.

        mid_thresh (:obj:`float`, optional): The lower bound for what is
            considered a mid level cloud.  If *vert_type* is 'pres', the
            default is 80000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 2000 m.

        high_thresh (:obj:`float`, optional): The lower bound for what is
            considered a high level cloud.  If *vert_type* is 'pres', the
            default is 45000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 6000 m.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        cloud fraction array for low level clouds.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_cloudfrac(wrfin, timeidx, method, squeeze,
                           cache, meta, _key, vert_type, low_thresh,
                           mid_thresh, high_thresh, missing)[0, :]

    if meta:
        result.attrs["description"] = "low clouds"

    return result


def get_mid_cloudfrac(wrfin, timeidx=0, method="cat", squeeze=True,
                      cache=None, meta=True, _key=None,
                      vert_type="height_agl", low_thresh=None,
                      mid_thresh=None, high_thresh=None,
                      missing=default_fill(np.float64)):
    """Return the cloud fraction for the mid level clouds.

    If the vertical coordinate type is 'height_agl' or 'height_msl', the
    default cloud levels are defined as:

        300 m <= low_cloud < 2000 m
        2000 m <= mid_cloud < 6000 m
        6000 m <= high_cloud

    For 'pressure', the default cloud levels are defined as:

        97000 Pa <= low_cloud < 80000 Pa
        80000 Pa <= mid_cloud < 45000 Pa
        45000 Pa <= high_cloud

    Note that the default low cloud levels are chosen to
    exclude clouds near the surface (fog). If you want fog included, set
    *low_thresh* to ~99500 Pa if *vert_type* is set to 'pressure', or 15 m if
    using 'height_msl' or 'height_agl'. Keep in mind that the lowest mass grid
    points are slightly above the ground, and in order to find clouds, the
    *low_thresh* needs to be set to values that are slightly greater than
    (less than) the lowest height (pressure) values.

    When using 'pressure' or 'height_agl' for *vert_type*, there is a
    possibility that the lowest WRF level will be higher than the low_cloud or
    mid_cloud threshold, particularly for mountainous regions.  When this
    happens, a fill value will be used in the output.

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

        vert_type (:obj:`str`, optional):  The type of vertical coordinate used
            to determine cloud type thresholds.  Must be 'height_agl',
            'height_msl', or 'pres'. The default is 'height_agl'.

        low_thresh (:obj:`float`, optional): The lower bound for what is
            considered a low cloud.  If *vert_type* is 'pres', the default is
            97000 Pa.  If *vert_type* is 'height_agl' or 'height_msl', then the
            default is 300 m.

        mid_thresh (:obj:`float`, optional): The lower bound for what is
            considered a mid level cloud.  If *vert_type* is 'pres', the
            default is 80000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 2000 m.

        high_thresh (:obj:`float`, optional): The lower bound for what is
            considered a high level cloud.  If *vert_type* is 'pres', the
            default is 45000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 6000 m.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        cloud fraction array for mid level clouds.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_cloudfrac(wrfin, timeidx, method, squeeze,
                           cache, meta, _key,
                           vert_type, low_thresh, mid_thresh,
                           high_thresh, missing)[1, :]

    if meta:
        result.attrs["description"] = "mid clouds"

    return result


def get_high_cloudfrac(wrfin, timeidx=0, method="cat", squeeze=True,
                       cache=None, meta=True, _key=None,
                       vert_type="height_agl", low_thresh=None,
                       mid_thresh=None, high_thresh=None,
                       missing=default_fill(np.float64)):
    """Return the cloud fraction for the high level clouds.

    If the vertical coordinate type is 'height_agl' or 'height_msl', the
    default cloud levels are defined as:

        300 m <= low_cloud < 2000 m
        2000 m <= mid_cloud < 6000 m
        6000 m <= high_cloud

    For 'pressure', the default cloud levels are defined as:

        97000 Pa <= low_cloud < 80000 Pa
        80000 Pa <= mid_cloud < 45000 Pa
        45000 Pa <= high_cloud

    Note that the default low cloud levels are chosen to
    exclude clouds near the surface (fog). If you want fog included, set
    *low_thresh* to ~99500 Pa if *vert_type* is set to 'pressure', or 15 m if
    using 'height_msl' or 'height_agl'. Keep in mind that the lowest mass grid
    points are slightly above the ground, and in order to find clouds, the
    *low_thresh* needs to be set to values that are slightly greater than
    (less than) the lowest height (pressure) values.

    When using 'pressure' or 'height_agl' for *vert_type*, there is a
    possibility that the lowest WRF level will be higher than the low_cloud or
    mid_cloud threshold, particularly for mountainous regions.  When this
    happens, a fill value will be used in the output.

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

        vert_type (:obj:`str`, optional):  The type of vertical coordinate used
            to determine cloud type thresholds.  Must be 'height_agl',
            'height_msl', or 'pres'. The default is 'height_agl'.

        low_thresh (:obj:`float`, optional): The lower bound for what is
            considered a low cloud.  If *vert_type* is 'pres', the default is
            97000 Pa.  If *vert_type* is 'height_agl' or 'height_msl', then the
            default is 300 m.

        mid_thresh (:obj:`float`, optional): The lower bound for what is
            considered a mid level cloud.  If *vert_type* is 'pres', the
            default is 80000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 2000 m.

        high_thresh (:obj:`float`, optional): The lower bound for what is
            considered a high level cloud.  If *vert_type* is 'pres', the
            default is 45000 Pa.  If *vert_type* is 'height_agl' or
            'height_msl', then the default is 6000 m.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        cloud fraction array for high level clouds.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_cloudfrac(wrfin, timeidx, method, squeeze,
                           cache, meta, _key,
                           vert_type, low_thresh, mid_thresh,
                           high_thresh, missing)[2, :]

    if meta:
        result.attrs["description"] = "high clouds"

    return result
