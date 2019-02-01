from __future__ import (absolute_import, division, print_function)

import warnings

from .constants import Constants
from .destag import destagger
from .decorators import convert_units
from .metadecorators import set_height_metadata
from .util import extract_vars, either


def _get_geoht(wrfin, timeidx, method="cat", squeeze=True,
               cache=None, meta=True, _key=None,
               height=True, msl=True, stag=False):
    """Return the geopotential or geopotential height.

    If *height* is False, then geopotential is returned in units of
    [m2 s-2].  If *height* is True, then geopotential height is
    returned in units of [m].  If *msl* is True, then geopotential height
    is return as Mean Sea Level (MSL).  If *msl* is False, then geopotential
    height is returned as Above Ground Level (AGL).

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

        height (:obj:`bool`, optional): Set to True to return geopotential
            height instead of geopotential.  Default is True.

        msl (:obj:`bool`, optional): Set to True to return geopotential height
            as Mean Sea Level (MSL).  Set to False to return the
            geopotential height as Above Ground Level (AGL) by subtracting
            the terrain height.  Default is True.

        stag (:obj:`bool`, optional): Set to True to use the vertical
            staggered grid, rather than the mass grid. Default is False.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        geopotential or geopotential height.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    varname = either("PH", "GHT")(wrfin)
    if varname == "PH":
        ph_vars = extract_vars(wrfin, timeidx, ("PH", "PHB", "HGT"),
                               method, squeeze, cache, meta=False,
                               _key=_key)
        ph = ph_vars["PH"]
        phb = ph_vars["PHB"]
        hgt = ph_vars["HGT"]
        geopt = ph + phb
        if not stag:
            geopt_unstag = destagger(geopt, -3)
        else:
            geopt_unstag = geopt
    else:
        ght_vars = extract_vars(wrfin, timeidx, ("GHT", "HGT_M"),
                                method, squeeze, cache, meta=False,
                                _key=_key)
        geopt_unstag = ght_vars["GHT"] * Constants.G
        hgt = ght_vars["HGT_M"]

        if stag:
            warnings.warn("file contains no vertically staggered geopotential "
                          "height variable, returning unstaggered result "
                          "instead")
    if height:
        if msl:
            return geopt_unstag / Constants.G
        else:
            # Due to broadcasting with multifile/multitime, the 2D terrain
            # array needs to be reshaped to a 3D array so the right dims
            # line up
            new_dims = list(hgt.shape)
            new_dims.insert(-2, 1)
            hgt = hgt.reshape(new_dims)

            return (geopt_unstag / Constants.G) - hgt
    else:
        return geopt_unstag


@set_height_metadata(geopt=True, stag=False)
def get_geopt(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
              meta=True, _key=None):
    """Return the geopotential.

    The geopotential is returned in units of [m2 s-2].

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
        geopotential.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    return _get_geoht(wrfin, timeidx, method, squeeze, cache, meta, _key,
                      False, True)


@set_height_metadata(geopt=False, stag=False)
@convert_units("height", "m")
def get_height(wrfin, timeidx=0, method="cat", squeeze=True,
               cache=None, meta=True, _key=None,
               msl=True, units="m"):
    """Return the geopotential height.

    If *msl* is True, then geopotential height is returned as Mean Sea Level
    (MSL).  If *msl* is False, then geopotential height is returned as
    Above Ground Level (AGL) by subtracting the terrain height.

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

        msl (:obj:`bool`, optional): Set to True to return geopotential height
            as Mean Sea Level (MSL).  Set to False to return the
            geopotential height as Above Ground Level (AGL) by subtracting
            the terrain height.  Default is True.

        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar`
            product table for a list of available units for 'z'.  Default
            is 'm'.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        geopotential height.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    return _get_geoht(wrfin, timeidx, method, squeeze, cache, meta, _key,
                      True, msl)


@set_height_metadata(geopt=True, stag=True)
def get_stag_geopt(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
                   meta=True, _key=None):
    """Return the geopotential for the vertically staggered grid.

    The geopotential is returned in units of [m2 s-2].

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
        geopotential.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    return _get_geoht(wrfin, timeidx, method, squeeze, cache, meta, _key,
                      False, True, stag=True)


@set_height_metadata(geopt=False, stag=True)
@convert_units("height", "m")
def get_stag_height(wrfin, timeidx=0, method="cat", squeeze=True,
                    cache=None, meta=True, _key=None, msl=True, units="m"):
    """Return the geopotential height for the vertically staggered grid.

    If *msl* is True, then geopotential height is returned as Mean Sea Level
    (MSL).  If *msl* is False, then geopotential height is returned as
    Above Ground Level (AGL) by subtracting the terrain height.

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

        msl (:obj:`bool`, optional): Set to True to return geopotential height
            as Mean Sea Level (MSL).  Set to False to return the
            geopotential height as Above Ground Level (AGL) by subtracting
            the terrain height.  Default is True.

        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar`
            product table for a list of available units for 'z'.  Default
            is 'm'.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        geopotential height.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    return _get_geoht(wrfin, timeidx, method, squeeze, cache, meta, _key,
                      True, msl, stag=True)


@set_height_metadata(geopt=False, stag=False)
@convert_units("height", "m")
def get_height_agl(wrfin, timeidx=0, method="cat", squeeze=True,
                   cache=None, meta=True, _key=None, units="m"):
    """Return the geopotential height (AGL).

    The geopotential height is returned as Above Ground Level (AGL) by
    subtracting the terrain height.

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

        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar`
            product table for a list of available units for 'height_agl'.
            Default is 'm'.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        geopotential height.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    return _get_geoht(wrfin, timeidx, method, squeeze, cache, meta, _key,
                      True, False)
