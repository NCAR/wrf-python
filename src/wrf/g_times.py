from __future__ import (absolute_import, division, print_function)

from .util import extract_times


def get_times(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
              meta=True, _key=None):
    """Return a sequence of time objects.

    This function reads the 'Times' variable and creates a sequence of
    :class:`datetime.datetime` objects.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence.

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

        meta (:obj:`bool`, optional): Set to False to disable metadata.

        _key (:obj:`int`, optional): A caching key. This is used for internal
            purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: A sequence of
        :class:`datetime.datetime` objects.  If *meta* is True, the sequence
        will be of type :class:`xarray.DataArray`, otherwise the sequence is
        :class:`numpy.ndarray`.

    """
    return extract_times(wrfin, timeidx, method, squeeze, cache,
                         meta=meta, do_xtime=False)


def get_xtimes(wrfin, timeidx=0, method="cat", squeeze=True, cache=None,
               meta=True, _key=None):
    """Return a sequence of time objects.

    This function reads the 'XTIME' variable and creates a sequence of
    :obj:`float` objects.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence.

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

        meta (:obj:`bool`, optional): Set to False to disable metadata.

        _key (:obj:`int`, optional): A caching key. This is used for internal
            purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: A sequence of
        :obj:`float` objects.  If *meta* is True, the sequence will be of type
        :class:`xarray.DataArray`, otherwise the sequence is
        :class:`numpy.ndarray`.

    Raises:

        :class:`KeyError`: Raised if the 'XTIME' variable is not present
            in the NetCDF file.

    """
    return extract_times(wrfin, timeidx, method, squeeze, cache,
                         meta=meta, do_xtime=True)
