from __future__ import (absolute_import, division, print_function)

import numpy as np
import numpy.ma as ma

from .extension import _tk, _cape
from .destag import destagger
from .constants import default_fill, Constants, ConversionFactors
from .util import extract_vars
from .metadecorators import set_cape_metadata


@set_cape_metadata(is2d=True)
def get_2dcape(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
               meta=True, _key=None, missing=default_fill(np.float64)):
    """Return the two-dimensional fields of MCAPE, MCIN, LCL, and LFC.

    The leftmost dimension of the returned array represents four different
    quantities:

        - return_val[0,...] will contain MCAPE [J kg-1]
        - return_val[1,...] will contain MCIN [J kg-1]
        - return_val[2,...] will contain LCL [m]
        - return_val[3,...] will contain LFC [m]

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        cape, cin, lcl, and lfc values as an array whose
        leftmost dimension is 4 (0=CAPE, 1=CIN, 2=LCL, 3=LFC).
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """

    varnames = ("T", "P", "PB", "QVAPOR", "PH", "PHB", "HGT", "PSFC")
    ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)

    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    ter = ncvars["HGT"]
    psfc = ncvars["PSFC"]

    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)

    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    z = geopt_unstag/Constants.G

    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc

    i3dflag = 0
    ter_follow = 1

    cape_cin = _cape(p_hpa, tk, qv, z, ter, psfc_hpa, missing, i3dflag,
                     ter_follow)

    left_dims = cape_cin.shape[1:-3]
    right_dims = cape_cin.shape[-2:]

    resdim = (4,) + left_dims + right_dims

    # Make a new output array for the result
    result = np.zeros(resdim, cape_cin.dtype)

    # Cape 2D output is not flipped in the vertical, so index from the
    # end
    result[0, ..., :, :] = cape_cin[0, ..., -1, :, :]
    result[1, ..., :, :] = cape_cin[1, ..., -1, :, :]
    result[2, ..., :, :] = cape_cin[1, ..., -2, :, :]
    result[3, ..., :, :] = cape_cin[1, ..., -3, :, :]

    return ma.masked_values(result, missing)


@set_cape_metadata(is2d=False)
def get_3dcape(wrfin, timeidx=0, method="cat",
               squeeze=True, cache=None, meta=True,
               _key=None, missing=default_fill(np.float64)):
    """Return the three-dimensional CAPE and CIN.

    The leftmost dimension of the returned array represents two different
    quantities:

        - return_val[0,...] will contain CAPE [J kg-1]
        - return_val[1,...] will contain CIN [J kg-1]

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        CAPE and CIN as an array whose leftmost dimension is 2 (0=CAPE, 1=CIN).
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    varnames = ("T", "P", "PB", "QVAPOR", "PH", "PHB", "HGT", "PSFC")
    ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                          meta=False, _key=_key)
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    ter = ncvars["HGT"]
    psfc = ncvars["PSFC"]

    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)

    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    z = geopt_unstag/Constants.G

    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc

    i3dflag = 1
    ter_follow = 1

    cape_cin = _cape(p_hpa, tk, qv, z, ter, psfc_hpa, missing, i3dflag,
                     ter_follow)

    return ma.masked_values(cape_cin, missing)


def get_cape2d_only(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
                    meta=True, _key=None, missing=default_fill(np.float64)):
    """Return the two-dimensional field of MCAPE (Max Convective Available
    Potential Energy).

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        2D MCAPE field.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_2dcape(wrfin, timeidx, method, squeeze, cache,
                        meta, _key, missing)[0, :]

    if meta:
        result.attrs["description"] = "mcape"
        result.attrs["units"] = "J kg-1"

    return result


def get_cin2d_only(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
                   meta=True, _key=None, missing=default_fill(np.float64)):
    """Return the two-dimensional field of MCIN (Max Convective Inhibition).

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        2D MCIN field.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_2dcape(wrfin, timeidx, method, squeeze, cache,
                        meta, _key, missing)[1, :]

    if meta:
        result.attrs["description"] = "mcin"
        result.attrs["units"] = "J kg-1"

    return result


def get_lcl(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
            meta=True, _key=None, missing=default_fill(np.float64)):
    """Return the two-dimensional field of LCL (Lifted Condensation Level).

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        2D LCL field.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_2dcape(wrfin, timeidx, method, squeeze, cache,
                        meta, _key, missing)[2, :]

    if meta:
        result.attrs["description"] = "lcl"
        result.attrs["units"] = "m"

    return result


def get_lfc(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
            meta=True, _key=None, missing=default_fill(np.float64)):
    """Return the two-dimensional field of LFC (Level of Free Convection).

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        2D LFC field.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_2dcape(wrfin, timeidx, method, squeeze, cache,
                        meta, _key, missing)[3, :]

    if meta:
        result.attrs["description"] = "lfc"
        result.attrs["units"] = "m"

    return result


def get_3dcape_only(wrfin, timeidx=0, method="cat",
                    squeeze=True, cache=None, meta=True,
                    _key=None, missing=default_fill(np.float64)):
    """Return the three-dimensional field of CAPE (Convective Available
    Potential Energy).

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        3D CAPE field.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_3dcape(wrfin, timeidx, method, squeeze, cache, meta,
                        _key, missing)[0, :]

    if meta:
        result.attrs["description"] = "cape"
        result.attrs["units"] = "J kg-1"

    return result


def get_3dcin_only(wrfin, timeidx=0, method="cat",
                   squeeze=True, cache=None, meta=True,
                   _key=None, missing=default_fill(np.float64)):
    """Return the three-dimensional field of CIN (Convective Inhibition).

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

        missing (:obj:`float`): The fill value to use for the output.
            Default is :data:`wrf.default_fill(np.float64)`.

    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The
        3D CIN field.
        If xarray is enabled and the *meta* parameter is True, then the result
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will
        be a :class:`numpy.ndarray` object with no metadata.

    """
    result = get_3dcape(wrfin, timeidx, method, squeeze, cache, meta,
                        _key, missing)[1, :]

    if meta:
        result.attrs["description"] = "cin"
        result.attrs["units"] = "J kg-1"

    return result
