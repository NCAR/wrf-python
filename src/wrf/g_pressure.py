from __future__ import (absolute_import, division, print_function)

from .decorators import convert_units
from .metadecorators import copy_and_set_metadata
from .util import extract_vars, either


@copy_and_set_metadata(copy_varname=either("P", "PRES"), name="pressure",
                       description="pressure")
@convert_units("pressure", "pa")
def get_pressure(wrfin, timeidx=0, method="cat", squeeze=True,
                 cache=None, meta=True, _key=None,
                 units="Pa"):

    """Return the pressure in the specified units.

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
            product table for a list of available units for 'pres'.  Default
            is 'Pa'.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The pressure in
        the specified units.
        If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """
    varname = either("P", "PRES")(wrfin)
    if varname == "P":
        p_vars = extract_vars(wrfin, timeidx, ("P", "PB"),
                              method, squeeze, cache, meta=False,
                              _key=_key)
        p = p_vars["P"]
        pb = p_vars["PB"]
        pres = p + pb
    else:
        pres = extract_vars(wrfin, timeidx, "PRES",
                            method, squeeze, cache, meta=False,
                            _key=_key)["PRES"]

    return pres


def get_pressure_hpa(wrfin, timeidx=0, method="cat", squeeze=True,
                     cache=None, meta=True, _key=None):
    """Return the pressure in [hPa].

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
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The pressure in
        [hPa].
        If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """
    return get_pressure(wrfin, timeidx, method, squeeze, cache, meta, _key,
                        units="hPa")
