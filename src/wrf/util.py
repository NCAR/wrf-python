from __future__ import (absolute_import, division, print_function)

import os
from sys import version_info
from copy import copy
from collections import OrderedDict
from collections.abc import Iterable, Mapping
from itertools import product, tee
from types import GeneratorType
import datetime as dt
from inspect import getmodule

try:
    from inspect import signature
except ImportError:
    from inspect import getargspec

try:
    from inspect import getargvalues
except ImportError:
    from inspect import getgeneratorlocals

# Needed for Python 3.4 to use apply_defaults
try:
    from inspect import _empty, _VAR_POSITIONAL, _VAR_KEYWORD
except ImportError:
    pass

import numpy as np
import numpy.ma as ma

from .config import xarray_enabled
from .constants import default_fill, ALL_TIMES
from .py3compat import (viewitems, viewkeys, isstr, py3range)
from .cache import cache_item, get_cached_item
from .geobnds import GeoBounds, NullGeoBounds
from .coordpair import CoordPair
from .projection import getproj


if xarray_enabled():
    from xarray import DataArray


_COORD_PAIR_MAP = {"XLAT": ("XLAT", "XLONG"),
                   "XLONG": ("XLAT", "XLONG"),
                   "XLAT_M": ("XLAT_M", "XLONG_M"),
                   "XLONG_M": ("XLAT_M", "XLONG_M"),
                   "XLAT_U": ("XLAT_U", "XLONG_U"),
                   "XLONG_U": ("XLAT_U", "XLONG_U"),
                   "XLAT_V": ("XLAT_V", "XLONG_V"),
                   "XLONG_V": ("XLAT_V", "XLONG_V"),
                   "CLAT": ("CLAT", "CLONG"),
                   "CLONG": ("CLAT", "CLONG")}


_COORD_VARS = ("XLAT", "XLONG", "XLAT_M", "XLONG_M", "XLAT_U", "XLONG_U",
               "XLAT_V", "XLONG_V", "CLAT", "CLONG")

_LAT_COORDS = ("XLAT", "XLAT_M", "XLAT_U", "XLAT_V", "CLAT")

_LON_COORDS = ("XLONG", "XLONG_M", "XLONG_U", "XLONG_V", "CLONG")

_TIME_COORD_VARS = ("XTIME",)


def is_time_coord_var(varname):
    """Return True if the input variable name is a time coordinate.

    Args:

        varname (:obj:`str`): The input variable name.

    Returns:

        :obj:`bool`: True if the input variable is a time coordinate,
        otherwise False.

    """
    return varname in _TIME_COORD_VARS


def get_coord_pairs(coord_varname):
    """Return a :obj:`tuple` for the variable names of the coordinate pair used
    for the 2D curvilinear coordinate variable.

    For example, the 'XLAT' variable will have coordinate variables of
    ('XLAT', 'XLONG') since the 'XLAT' variable itself is two-dimensional.

    Args:

        coord_varname (:obj:`str`): The coordinate variable name.

    Returns:

        :obj:`bool`: True if the time index is :data:`wrf.ALL_TIMES` or
        :obj:`None`, otherwise False.

    """
    return _COORD_PAIR_MAP[coord_varname]


def is_multi_time_req(timeidx):
    """Return True if the requested time index is for :data:`wrf.ALL_TIMES` or
    :obj:`None`.

    Args:

        timeidx (:obj:`int`, :data:`wrf.ALL_TIMES`, or :obj:`None`): The
            requested time index.

    Returns:

        :obj:`bool`: True if the time index is :data:`wrf.ALL_TIMES` or
        :obj:`None`, otherwise False.

    """
    return timeidx is None


def is_multi_file(wrfin):
    """Return True if the input argument is an iterable.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

    Returns:

        :obj:`bool`: True if the input is an iterable.  False if the input
        is a single NetCDF file object.

    """
    return (isinstance(wrfin, Iterable) and not isstr(wrfin))


def has_time_coord(wrfnc):
    """Return True if the input file or sequence contains the time
    coordinate variable.

    The time coordinate is named 'XTIME'.

    Args:

        wrfnc (:class:`netCDF4.Dataset` or :class:`Nio.NioFile`): A single
            NetCDF file object.

    Returns:

        :obj:`bool`: True if the netcdf file contains the time coordinate
        variable, False otherwise.

    """
    return "XTIME" in wrfnc.variables


def is_mapping(wrfin):
    """Return True if the input file or sequence is a mapping type.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

    Returns:

        :obj:`bool`: True if the input is a mapping type, False otherwise.

    """
    return isinstance(wrfin, Mapping)


def _generator_copy(gen):
    """Return a copy of a generator.

    This function instantiates a new generator object using the arguments
    passed to the original.

    Args:

        gen (:class:`types.GeneratorType`): A generator object.

    Note:

        In order for this to work correctly, the generator cannot modify
        the original construction arguments.

    Returns:

        :class:`types.GeneratorType`: A copy of the generator object.

    """
    funcname = gen.__name__

    # This is for generator expressions.  Only solution is to tee the
    # original generator.
    if funcname == "<genexpr>":
        return tee(gen)

    try:
        argvals = getargvalues(gen.gi_frame)
    except NameError:
        argvals = getgeneratorlocals(gen)
    module = getmodule(gen.gi_frame)

    if module is not None:
        try:
            try:
                argd = {key: argvals.locals[key] for key in argvals.args}
                res = module.get(funcname)(**argd)
            except AttributeError:
                res = getattr(module, funcname)(**argd)
        except Exception:
            # This is the old way it used to work, but it looks like this was
            # fixed by Python.
            try:
                res = module.get(funcname)(**argvals.locals)
            except AttributeError:
                res = getattr(module, funcname)(**argvals.locals)
    else:
        # Created in jupyter or the python interpreter
        import __main__

        try:
            argd = {key: argvals.locals[key] for key in argvals.args}
            res = getattr(__main__, funcname)(**argd)
        except Exception:
            # This was the old way it used to work, but appears to have
            # been fixed by Python.
            res = getattr(__main__, funcname)(**argvals.locals)

    return res


def test():
    q = [1, 2, 3]
    for i in q:
        yield i


class TestGen(object):
    def __init__(self, count=3):
        self._total = count
        self._i = 0

    def __iter__(self):
        return self

    def next(self):
        if self._i >= self._total:
            raise StopIteration
        else:
            val = self._i
            self._i += 1
            return val

    # Python 3
    def __next__(self):
        return self.next()


def latlon_coordvars(ncvars):
    """Return the first found latitude and longitude coordinate names from a
    NetCDF variable dictionary.

    This function searches the dictionary structure for NetCDF variables
    and returns the first found latitude and longitude coordinate
    names (typically 'XLAT' and 'XLONG').

    Args:

        ncvars (:obj:`dict`): A NetCDF variable dictionary object.

    Returns:

        :obj:`tuple`: The latitude and longitude coordinate name pairs as
        (lat_coord_name, lon_coord_name).

    """
    lat_coord = None
    lon_coord = None

    for name in _LAT_COORDS:
        if name in viewkeys(ncvars):
            lat_coord = name
            break

    for name in _LON_COORDS:
        if name in viewkeys(ncvars):
            lon_coord = name
            break

    return lat_coord, lon_coord


def is_coordvar(varname):
    """Returns True if the variable is a coordinate variable.

    Args:

        varname (:obj:`str`):  The variable name.

    Returns:

        :obj:`bool`: True if the variable is a coordinate variable,
        otherwise False.

    """
    return varname in _COORD_VARS


class IterWrapper(Iterable):
    """A wrapper class for generators and custom iterable classes that returns
    a new iterator to the start of the sequence when
    :meth:`IterWrapper.__iter__` is called.

    If the wrapped object is a generator, a copy of the generator is
    constructed and returned when :meth:`IterWrapper.__iter__` is called.
    If the wrapped object is a custom type, then the :meth:`copy.copy` is
    called and a new instance is returned.  In both cases, the original
    iterable object is unchanged.

    Note:

        Do not increment the wrapped iterable outside of this wrapper.

    """
    def __init__(self, wrapped):
        """Initialize an :class:`wrf.IterWrapper` object.

        Args:

            wrapped (an iterable object): Any iterable object that contains the
                *__iter__* method.

        """
        self._wrapped = wrapped

    def __iter__(self):
        """Return an iterator to the start of the sequence.

        Returns:

        An iterator to the start of the sequence.

        """
        if isinstance(self._wrapped, GeneratorType):

            gen_copy = _generator_copy(self._wrapped)
            # If a tuple comes back, then this is a generator expression,
            # so store the first tee'd item, then return the other
            if isinstance(gen_copy, tuple):
                self._wrapped = gen_copy[0]
                return gen_copy[1]

            return gen_copy
        else:
            obj_copy = copy(self._wrapped)
            return obj_copy.__iter__()


def get_iterable(wrfseq):
    """Returns a resettable iterable object.

    In this context, resettable means that when :meth:`object.__iter__`
    is called, the iterable returned always points to the first element
    in the sequence, similar to how the list and tuple behave.

    Args:

        wrfseq (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

    Returns:

        iterable: A resettable iterator object.

    """
    if not is_multi_file(wrfseq):
        return wrfseq
    else:
        if not is_mapping(wrfseq):

            if isinstance(wrfseq, (list, tuple, IterWrapper)):
                return wrfseq
            else:
                return IterWrapper(wrfseq)  # generator/custom iterable class

        else:
            if isinstance(wrfseq, dict):
                return wrfseq
            else:
                return dict(wrfseq)  # generator/custom iterable dict class


# Helper to extract masked arrays from DataArrays that convert to NaN
def to_np(array):
    """Return the :class:`numpy.ndarray` contained in an
    :class:`xarray.DataArray` instance.

    If the :class:`xarray.DataArray` instance does not contain a *_FillValue*
    or *missing_value* attribute, then this routine simply returns the
    :attr:`xarray.DataArray.values` attribute.  If the
    :class:`xarray.DataArray` object contains a *_FillValue* or *missing_value*
    attribute, then this routine returns a :class:`numpy.ma.MaskedArray`
    instance, where the NaN values (used by xarray to represent missing data)
    are replaced with the fill value.

    If the object passed in to this routine is not an
    :class:`xarray.DataArray` instance, then this routine simply returns the
    passed in object.  This is useful in situations where you do not know
    if you have an :class:`xarray.DataArray` or a :class:`numpy.ndarray` and
    simply want a :class:`numpy.ndarray` returned.

    Args:

        array (:class:`xarray.DataArray`, :class:`numpy.ndarray`, or any \
        object): Can be any object type, but is generally
            used with :class:`xarray.DataArray` or :class:`numpy.ndarray`.

    Returns:

        :class:`numpy.ndarray` or :class:`numpy.ma.MaskedArray`: The
        extracted array or the *array* object if *array* is not a
        :class:`xarray.DataArray` object..

    """

    try:
        fill_value = array.attrs["_FillValue"]
    except AttributeError:
        result = array  # Not a DataArray
    except KeyError:
        result = array.values  # Does not have missing values
    else:
        result = ma.masked_invalid(array.values, copy=False)
        result.set_fill_value(fill_value)

    return result


# Helper utilities for metadata
class either(object):
    """A callable class that determines which variable is present in the
    file.

    This is used in situations where the same variable type has different names
    depending on the type of file used.  For example, in a WRF output file,
    'P' is used for pressure, whereas in a met_em file, pressure is named
    'PRES'.

    Methods:

        __call__(wrfin): Return the variable that is present in the file.

            Args:

                wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
                iterable): WRF-ARW NetCDF
                    data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
                    or an iterable sequence of the aforementioned types.

            Returns:

                :obj:`str`: The variable name that is present in the file.

    Attributes:

        varnames (sequence): A sequence of possible variable names.

    """
    def __init__(self, *varnames):
        """Initialize an :class:`either` object.

        Args:

            *varnames (sequence): A sequence of possible variable names.

        """
        self.varnames = varnames

    def __call__(self, wrfin):
        if is_multi_file(wrfin):
            if not is_mapping(wrfin):
                wrfin = next(iter(wrfin))
            else:
                entry = wrfin[next(iter(viewkeys(wrfin)))]
                return self(entry)

        for varname in self.varnames:
            if varname in wrfin.variables:
                return varname

        raise ValueError("{} are not valid variable names".format(
            self.varnames))


# This should look like:
# [(0, (-3,-2)), # variable 1
# (1, -1)] # variable 2
class combine_dims(object):
    """A callable class that mixes dimension sizes from different function
    arguments.

    This callable object is used for determining the output size for the
    extension module functions.  The class is initialized with a sequence of
    pairs, where the first value is the function argument index.  The second
    value is the dimension index(es) to use.  The second value can be a
    single integer or a sequence if multiple dimensions are used.

    Methods:

        __call__(*args): Return a tuple with the combined dimension sizes.

            Args:

                *args: The function arguments from which to extract the
                    dimensions sizes.

            Returns:

                :obj:`tuple`:  The shape for the combined dimensions.

    Attributes:

        pairs (sequence): A sequence representing how to combine the
            dimensions.

    Example:

        .. code-block:: python

            # Take the -3, and -2 dimension sizes from argument 0
            # Then take the -1 dimension size from argument 1
            pairs = [(0, (-3, -2), (1, -1)]

            combo = combine_dims(pairs)

    """

    def __init__(self, pairs):
        """Initialize a :class:`combine_dims` object.

        Args:

            pairs (sequence): A sequence where each element is a pair
                (:obj:`tuple`), with the first element being the function
                argument index and the second value being either an integer
                or sequence for the dimension size indexes to use.

        """
        self.pairs = pairs

    def __call__(self, *args):
        result = []
        for pair in self.pairs:
            var = args[pair[0]]
            dimidxs = pair[1]
            if isinstance(dimidxs, Iterable):
                for dimidx in dimidxs:
                    result.append(var.shape[dimidx])
            else:
                result.append(var.shape[dimidxs])

        return tuple(result)


class from_var(object):
    """A callable class that retrieves attributes from the function argument.

    If the function argument is not of type :class:`xarray.DataArray`, then
    None will be returned.

    It is assumed that the function has been wrapped using the :mod:`wrapt`
    module.

    Methods:

        __call__(wrapped, *args, **kwargs): Return the attribute found in \
        the function arguments.

            Args:
                wrapped: The wrapped function, as used by the :mod:`wrapt`
                    module.

                *args: The function arguments.

                **kwargs: The function keyword arguments.

            Returns:

                :obj:`object`: The requested attribute.

    Attributes:

        varname (:obj:`str`): The variable name.

        attribute (:obj:`str`): The attribute name.

    """
    def __init__(self, varname, attribute):
        """Initialize a :class:`from_var` object.

        Args:

            varname (:obj:`str`): The variable name.

            attribute (:obj:`str`): The attribute name.

        """
        self.varname = varname
        self.attribute = attribute

    def __call__(self, wrapped, *args, **kwargs):
        vard = from_args(wrapped, (self.varname,), *args, **kwargs)

        var = None
        if vard is not None:
            var = vard[self.varname]

        try:
            return var.attrs.get(self.attribute, None)
        except AttributeError:
            return None


def _corners_moved(wrfnc, first_ll_corner, first_ur_corner, latvar, lonvar):
    """Return True if the corner points have moved.

    This function is used to test for a moving domain, since the WRF output
    does not set any flags in the file for this.  The test will be performed
    for all time steps in the NetCDF file.

    Args:

        wrfnc (:class:`netCDF4.Dataset` or :class:`Nio.NioFile`): A single
            NetCDF file object.

        first_ll_corner (:obj:`tuple`): A (latitude, longitude) pair for the
            lower left corner found in the initial file.

        first_ur_corner (:obj:`tuple`): A (latitude, longitude) pair for the
            upper right corner found in the initial file.

        latvar (:obj:`str`): The latitude variable name to use.

        lonvar (:obj:`str`: The longitude variable name to use.


    Returns:

        :obj:`bool`:  True if the corner points have moved, False otherwise.

    """
    lats = wrfnc.variables[latvar][:]
    lons = wrfnc.variables[lonvar][:]

    # Need to check all times
    for i in py3range(lats.shape[-3]):
        start_idxs = [0] * len(lats.shape)  # PyNIO does not support ndim
        start_idxs[-3] = i
        start_idxs = tuple(start_idxs)

        end_idxs = [-1] * len(lats.shape)
        end_idxs[-3] = i
        end_idxs = tuple(end_idxs)

        if (first_ll_corner[0] != lats[start_idxs]
                or first_ll_corner[1] != lons[start_idxs]
                or first_ur_corner[0] != lats[end_idxs]
                or first_ur_corner[1] != lons[end_idxs]):
            return True

    return False


def is_moving_domain(wrfin, varname=None, latvar=either("XLAT", "XLAT_M"),
                     lonvar=either("XLONG", "XLONG_M"), _key=None):

    """Return True if the domain is a moving nest.

    This function is used to test for a moving domain, since the WRF output
    does not set any flags in the file for this.  The test will be performed
    for all files in any sequences and across all times in each file.

    This result is cached internally, so this potentially lengthy check is
    only done one time for any given *wrfin* parameter.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        varname (:obj:`str`, optional): A specific NetCDF variable to test,
            which must contain a 'coordinates' attribute.  If unspecified,
            The *latvar* and *lonvar* parameters are used.  Default is None.

        first_ur_corner (:obj:`tuple`): A (latitude, longitude) pair for the
            upper right corner found in the initial file.

        latvar (:obj:`str` or :class:`either`, optional): The latitude variable
            name.  Default is :class:`either('XLAT', 'XLAT_M')`.

        lonvar (:obj:`str` or :class:`either`, optional): The latitude variable
            name.  Default is :class:`either('XLAT', 'XLAT_M')`.

        _key (:obj:`int`, optional): A caching key. This is used for internal
            purposes only.  Default is None.

    Returns:

        :obj:`bool`:  True if the domain is a moving nest, False otherwise.

    """

    if isinstance(latvar, either):
        latvar = latvar(wrfin)

    if isinstance(lonvar, either):
        lonvar = lonvar(wrfin)

    # In case it's just a single file
    if not is_multi_file(wrfin):
        wrfin = [wrfin]

    # Slow, but safe. Compare the corner points to the first item and see
    # any move.  User iterator protocol in case wrfin is not a list/tuple.
    if not is_mapping(wrfin):
        wrf_iter = iter(wrfin)
        first_wrfnc = next(wrf_iter)
    else:
        # Currently only checking the first dict entry.
        dict_key = next(iter(viewkeys(wrfin)))
        entry = wrfin[dict_key]
        key = _key[dict_key] if _key is not None else None
        return is_moving_domain(entry, varname, latvar, lonvar, key)

    # The long way of checking all lat/lon corner points.  Doesn't appear
    # to be a shortcut in the netcdf files.
    if varname is not None:
        try:
            coord_str = getattr(first_wrfnc.variables[varname], "coordinates")
            # scipy.io.netcdf stores attributes as bytes rather than str
            if isinstance(coord_str, str):
                coord_names = coord_str.split()
            else:
                coord_names = coord_str.decode().split()
        except AttributeError:
            # Variable doesn't have a coordinates attribute, use the
            # arguments
            lon_coord = lonvar
            lat_coord = latvar
        else:
            for name in coord_names:
                if name in _LAT_COORDS:
                    lat_coord = name
                    continue

                if name in _LON_COORDS:
                    lon_coord = name
                    continue
    else:
        lon_coord = lonvar
        lat_coord = latvar

    # See if there is a cached value
    product = "is_moving_{}_{}".format(lat_coord, lon_coord)
    moving = get_cached_item(_key, product)
    if moving is not None:
        return moving

    # Need to search all the files
    lats = first_wrfnc.variables[lat_coord][:]
    lons = first_wrfnc.variables[lon_coord][:]

    zero_idxs = [0]*len(lats.shape)  # PyNIO doesn't have ndim
    last_idxs = list(zero_idxs)
    last_idxs[-2:] = [-1]*2

    zero_idxs = tuple(zero_idxs)
    last_idxs = tuple(last_idxs)

    lat0 = lats[zero_idxs]
    lat1 = lats[last_idxs]
    lon0 = lons[zero_idxs]
    lon1 = lons[last_idxs]

    ll_corner = (lat0, lon0)
    ur_corner = (lat1, lon1)

    # Need to check if the first file is moving, might be a single
    # file with multiple times
    if _corners_moved(first_wrfnc, ll_corner, ur_corner, lat_coord, lon_coord):
        cache_item(_key, product, True)
        return True

    # Now check any additional files
    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            if _corners_moved(wrfnc, ll_corner, ur_corner,
                              lat_coord, lon_coord):

                cache_item(_key, product, True)
                return True

    cache_item(_key, product, False)

    return False


def _get_global_attr(wrfnc, attr):
    """Return the global attribute.

    Args:

        wrfnc (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`): A single
            WRF NetCDF file object.

        attr (:obj:`str`): The attribute name.

    Returns:

        :obj:`object`:  The global attribute.

    """
    val = getattr(wrfnc, attr, None)

    # PyNIO puts single values in to an array
    if isinstance(val, np.ndarray):
        if len(val) == 1:
            return val[0]
    return val


def extract_global_attrs(wrfin, attrs):
    """Return the global attribute(s).

    If the *wrfin* parameter is a sequence, then only the first element is
    used, so the entire sequence must have the same global attributes.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        attrs (:obj:`str`, sequence): The attribute name
            or a sequence of attribute names.

    Returns:

        :obj:`dict`:  A mapping of attribute_name to value.

    """
    if isstr(attrs):
        attrlist = [attrs]
    else:
        attrlist = attrs

    multifile = is_multi_file(wrfin)

    if multifile:
        if not is_mapping(wrfin):
            wrfin = next(iter(wrfin))
        else:
            entry = wrfin[next(iter(viewkeys(wrfin)))]
            return extract_global_attrs(entry, attrs)

    return {attr: _get_global_attr(wrfin, attr) for attr in attrlist}


def extract_dim(wrfin, dim):
    """Return the dimension size for the specified dimension name.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        dim (:obj:`str`): The dimension name.

    Returns:

        :obj:`int`:  The dimension size.

    """
    if is_multi_file(wrfin):
        if not is_mapping(wrfin):
            wrfin = next(iter(wrfin))
        else:
            entry = wrfin[next(iter(viewkeys(wrfin)))]
            return extract_dim(entry, dim)

    d = wrfin.dimensions[dim]
    if not isinstance(d, int):
        try:
            return len(d)  # netCDF4
        except TypeError:  # scipy.io.netcdf
            # Scipy can't handled unlimited dimensions, so now we have to
            # figure it out
            try:
                s = wrfin.variables["P"].shape
            except Exception:
                raise ValueError("unsupported NetCDF reader")
            else:
                return s[-4]

    return d  # PyNIO


def _combine_dict(wrfdict, varname, timeidx, method, meta, _key):
    """Return an array object from a mapping input.

    The resulting array is the combination of all sequences associated with
    each key in the dictionary. The leftmost dimension will be the keys. The
    rightmost dimensions are the dimensions for the aggregated sequences of
    arrays, using either the 'cat' or 'join' *method* parameter value.
    If dictionaries are nested, then the outermost dictionary keys will be the
    leftmost dimension, followed by each subsequent dictionary's keys.

    If the order of the keys (and leftmost dimension ordering) matters, it is
    recommended that an :class:`OrderedDict` be used instead of a normal
    dictionary.  Otherwise, the leftmost dimension will be ordered by the
    iteration order.

    Args:

        wrfdict (mapping): A mapping of key to an array or
            key to a sequence of arrays.

        varname (:obj:`str`) : The variable name.

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

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """
    keynames = []
    numkeys = len(wrfdict)

    key_iter = iter(viewkeys(wrfdict))
    first_key = next(key_iter)
    keynames.append(first_key)

    is_moving = is_moving_domain(wrfdict, varname, _key=_key)

    _cache_key = _key[first_key] if _key is not None else None

    first_array = _extract_var(wrfdict[first_key], varname,
                               timeidx, is_moving=is_moving, method=method,
                               squeeze=False, cache=None, meta=meta,
                               _key=_cache_key)

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
            _cache_key = _key[key] if _key is not None else None
            vardata = _extract_var(wrfdict[key], varname, timeidx,
                                   is_moving=is_moving, method=method,
                                   squeeze=False, cache=None, meta=meta,
                                   _key=_cache_key)

            if outdata.shape[1:] != vardata.shape:
                raise ValueError("data sequences must have the "
                                 "same size for all dictionary keys")
            outdata[idx, :] = to_np(vardata)[:]
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


def _find_coord_names(coords):
    """Return the coordinate variables names found in a
    :attr:`xarray.DataArray.coords` mapping.

    Args:

        coords (mapping): A :attr:`xarray.DataArray.coords` mapping object.

    Returns:

        :obj:`tuple`: The latitude, longitude, and xtime variable names used
        in the coordinate mapping.

    """
    try:
        lat_coord = [name for name in _COORD_VARS[0::2] if name in coords][0]
    except IndexError:
        lat_coord = None

    try:
        lon_coord = [name for name in _COORD_VARS[1::2] if name in coords][0]
    except IndexError:
        lon_coord = None

    try:
        xtime_coord = [name for name in _TIME_COORD_VARS if name in coords][0]
    except IndexError:
        xtime_coord = None

    return lat_coord, lon_coord, xtime_coord


def _find_max_time_size(wrfseq):
    """Return the maximum number of times found in a sequence of
    WRF files.

    Args:

        wrfseq (sequence): A sequence of WRF NetCDF file objects.

    Returns:

        :obj:`int`: The maximum number of times found in a file.

    """
    wrf_iter = iter(wrfseq)

    max_times = 0
    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            t = extract_dim(wrfnc, "Time")
            max_times = t if t >= max_times else max_times

    return max_times


def _get_coord_names(wrfin, varname):

    # Need only the first real file
    if is_multi_file(wrfin):
        if not is_mapping(wrfin):
            wrfnc = next(iter(wrfin))
        else:
            entry = wrfin[next(iter(viewkeys(wrfin)))]
            return _get_coord_names(entry, varname)
    else:
        wrfnc = wrfin

    lat_coord = None
    lon_coord = None
    time_coord = None

    var = wrfnc.variables[varname]

    # WRF variables will have a coordinates attribute.  MET_EM files have
    # a stagger attribute which indicates the coordinate variable.
    try:
        # WRF files
        coord_attr = getattr(var, "coordinates")
    except AttributeError:

        if is_coordvar(varname):
            # Coordinate variable (most likely XLAT or XLONG)
            lat_coord, lon_coord = get_coord_pairs(varname)
            time_coord = None

            if has_time_coord(wrfnc):
                time_coord = "XTIME"

        elif is_time_coord_var(varname):
            lon_coord = None
            lat_coord = None
            time_coord = None
        else:
            try:
                # met_em files or old WRF files
                stag_attr = getattr(var, "stagger")
            except AttributeError:
                lon_coord = None
                lat_coord = None

                # Let's just check for xlat and xlong in this case
                if "XLAT" in wrfnc.variables:
                    lat_coord = "XLAT"
                    lon_coord = "XLONG"
            else:
                # For met_em files, use the stagger name to get the lat/lon var
                lat_coord = "XLAT_{}".format(stag_attr)
                lon_coord = "XLONG_{}".format(stag_attr)

                # If this coord name is missing, it might be an old WRF file
                if lat_coord not in wrfnc.variables:
                    lat_coord = None
                    lon_coord = None

                    if "XLAT" in wrfnc.variables:
                        lat_coord = "XLAT"
                        lon_coord = "XLONG"
    else:
        if isinstance(coord_attr, str):
            coord_names = coord_attr.split()
        else:
            coord_names = coord_attr.decode().split()

        for name in coord_names:
            if name in _LAT_COORDS:
                lat_coord = name
                continue

            if name in _LON_COORDS:
                lon_coord = name
                continue

            if name in _TIME_COORD_VARS:
                time_coord = name
                continue

        if time_coord is not None:
            # Make sure the time variable wasn't removed
            try:
                _ = wrfnc.variables[time_coord]
            except KeyError:
                time_coord = None

    return lat_coord, lon_coord, time_coord


def _build_data_array(wrfnc, varname, timeidx, is_moving_domain, is_multifile,
                      _key):
    """Return a :class:`xarray.DataArray` object for the desired variable in
    a single NetCDF file object.

    Args:

        wrfnc (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`): A single
            WRF NetCDF file object.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. The default is 0.

        is_moving_domain (:obj:`bool`): A boolean type that indicates if the
            NetCDF file object came from a moving nest.

        is_multifile (:obj:`bool`): A boolean type that indicates if the NetCDF
            file object came from a sequence.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray`:  An array object that contains metadata.

    """

    # Note:  wrfnc is always a single netcdf file object
    # is_moving_domain and is_multifile are arguments indicating if the
    # single file came from a sequence, and if that sequence is has a moving
    # domain.  Both arguments are used mainly for coordinate extraction and
    # caching.
    multitime = is_multi_time_req(timeidx)
    time_idx_or_slice = timeidx if not multitime else slice(None)
    var = wrfnc.variables[varname]
    if len(var.shape) > 1:
        data = var[time_idx_or_slice, :]
    else:
        data = var[time_idx_or_slice]

    # Want to preserve the time dimension
    if not multitime:
        if len(var.shape) > 1:
            data = data[np.newaxis, :]
        else:
            data = data[np.newaxis]

    attrs = OrderedDict()
    for dkey, val in viewitems(var.__dict__):
        # scipy.io adds these but don't want them
        if dkey in ("data", "_shape", "_size", "_typecode", "_attributes",
                    "maskandscale", "dimensions"):
            continue

        _dkey = dkey if isinstance(dkey, str) else dkey.decode()
        if isstr(val):
            _val = val
        else:
            if isinstance(val, bytes):
                _val = val.decode()  # scipy.io.netcdf
            else:
                _val = val

        attrs[_dkey] = _val

    dimnames = var.dimensions[-data.ndim:]

    lat_coord = lon_coord = time_coord = None

    try:
        if dimnames[-2] == "south_north" and dimnames[-1] == "west_east":
            lat_coord, lon_coord, time_coord = _get_coord_names(wrfnc, varname)
    except IndexError:
        pass

    coords = OrderedDict()

    # Handle lat/lon coordinates and projection information if available
    if lon_coord is not None and lat_coord is not None:
        # Using a cache for coordinate variables so the extraction only happens
        # once.
        lon_coord_dimkey = lon_coord + "_dim"
        lon_coord_valkey = lon_coord + "_val"
        lat_coord_dimkey = lat_coord + "_dim"
        lat_coord_valkey = lat_coord + "_val"

        lon_coord_dims = get_cached_item(_key, lon_coord_dimkey)
        lon_coord_vals = get_cached_item(_key, lon_coord_valkey)
        if lon_coord_dims is None or lon_coord_vals is None:
            lon_var = wrfnc.variables[lon_coord]
            lon_coord_dims = lon_var.dimensions
            lon_coord_vals = lon_var[:]

            # Only cache here if the domain is not moving, otherwise
            # caching is handled by cat/join
            if not is_moving_domain:
                cache_item(_key, lon_coord_dimkey, lon_coord_dims)
                cache_item(_key, lon_coord_valkey, lon_coord_vals)

        lat_coord_dims = get_cached_item(_key, lat_coord_dimkey)
        lat_coord_vals = get_cached_item(_key, lat_coord_valkey)
        if lat_coord_dims is None or lat_coord_vals is None:
            lat_var = wrfnc.variables[lat_coord]
            lat_coord_dims = lat_var.dimensions
            lat_coord_vals = lat_var[:]

            # Only cache here if the domain is not moving, otherwise
            # caching is done in cat/join
            if not is_moving_domain:
                cache_item(_key, lat_coord_dimkey, lat_coord_dims)
                cache_item(_key, lat_coord_valkey, lat_coord_vals)

        time_coord_vals = None
        if time_coord is not None:
            # If not from a multifile sequence, then cache the time
            # coordinate.  Otherwise, handled in cat/join/
            if not is_multifile:
                time_coord_vals = get_cached_item(_key, time_coord)

                if time_coord_vals is None:
                    time_coord_vals = wrfnc.variables[time_coord][:]

                    if not is_multifile:
                        cache_item(_key, time_coord, time_coord_vals)
            else:
                time_coord_vals = wrfnc.variables[time_coord][:]

        if multitime:
            if is_moving_domain:
                # Special case with a moving domain in a multi-time file,
                # otherwise the projection parameters don't change
                coords[lon_coord] = lon_coord_dims, lon_coord_vals
                coords[lat_coord] = lat_coord_dims, lat_coord_vals

            else:
                coords[lon_coord] = (lon_coord_dims[1:],
                                     lon_coord_vals[0, :])
                coords[lat_coord] = (lat_coord_dims[1:],
                                     lat_coord_vals[0, :])

            if time_coord is not None:
                coords[time_coord] = (lon_coord_dims[0], time_coord_vals)
        else:
            coords[lon_coord] = (lon_coord_dims[1:],
                                 lon_coord_vals[timeidx, :])
            coords[lat_coord] = (lat_coord_dims[1:],
                                 lat_coord_vals[timeidx, :])

            if time_coord is not None:
                coords[time_coord] = (lon_coord_dims[0],
                                      [time_coord_vals[timeidx]])

    proj_params = get_proj_params(wrfnc)
    proj = getproj(**proj_params)
    attrs["projection"] = proj

    if dimnames[0] == "Time":
        t = extract_times(wrfnc, timeidx, meta=False, do_xtime=False)
        if not multitime:
            t = [t]
        coords[dimnames[0]] = t

    data_array = DataArray(data, name=varname, dims=dimnames, coords=coords,
                           attrs=attrs)

    return data_array


def _find_forward(wrfseq, varname, timeidx, is_moving, meta, _key):
    """Find and return the array object within a sequence for a specific time
    index.

    Args:

        wrfseq (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int`): The desired time index. Must be positive.

        is_moving (:obj:`bool`): A boolean type that indicates if the
            sequence is a moving nest.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """

    wrf_iter = iter(wrfseq)
    comboidx = 0

    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            numtimes = extract_dim(wrfnc, "Time")

            if timeidx < comboidx + numtimes:
                filetimeidx = timeidx - comboidx

                if meta:
                    return _build_data_array(wrfnc, varname, filetimeidx,
                                             is_moving, True, _key)
                else:
                    var = wrfnc.variables[varname]
                    if len(var.shape) > 1:
                        result = var[filetimeidx, :]
                        return result[np.newaxis, :]  # So that nosqueeze works
                    else:
                        result = var[filetimeidx]
                        return result[np.newaxis]  # So that nosqueeze works

            else:
                comboidx += numtimes

    raise IndexError("timeidx {} is out of bounds".format(timeidx))


def _find_reverse(wrfseq, varname, timeidx, is_moving, meta, _key):
    """Find and return the array object within a sequence for a specific time
    index.

    The sequence is searched in reverse.

    Args:

        wrfseq (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int`): The desired time index. Must be negative.

        is_moving (:obj:`bool`): A boolean type that indicates if the
            sequence is a moving nest.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """

    try:
        revwrfseq = reversed(wrfseq)
    except TypeError:
        revwrfseq = reversed(list(wrfseq))

    wrf_iter = iter(revwrfseq)
    revtimeidx = -timeidx - 1

    comboidx = 0

    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            numtimes = extract_dim(wrfnc, "Time")

            if revtimeidx < comboidx + numtimes:
                # Finds the "forward" sequence index, then counts that
                # number back from the back of the ncfile times,
                # since the ncfile  needs to be iterated backwards as well.
                filetimeidx = numtimes - (revtimeidx - comboidx) - 1

                if meta:
                    return _build_data_array(wrfnc, varname, filetimeidx,
                                             is_moving, True, _key)
                else:
                    result = wrfnc.variables[varname][filetimeidx, :]
                    return result[np.newaxis, :]  # So that nosqueeze works
            else:
                comboidx += numtimes

    raise IndexError("timeidx {} is out of bounds".format(timeidx))


def _find_arr_for_time(wrfseq, varname, timeidx, is_moving, meta, _key):
    """Find and return the array object within a sequence for a specific time
    index.

    The sequence is searched in forward or reverse based on the time index
    chosen.

    Args:

        wrfseq (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int`): The desired time index.

        is_moving (:obj:`bool`): A boolean type that indicates if the
            sequence is a moving nest.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """
    if timeidx >= 0:
        return _find_forward(wrfseq, varname, timeidx, is_moving, meta, _key)
    else:
        return _find_reverse(wrfseq, varname, timeidx, is_moving, meta, _key)


def _cat_files(wrfseq, varname, timeidx, is_moving, squeeze, meta, _key):
    """Return an array object from a sequence of files using the concatenate
    method.

    The concatenate method aggregates all files in the sequence along the
    'Time' dimension, which will be the leftmost dimension.  No sorting is
    performed, so all files in the sequence must be sorted prior to calling
    this method.


    Args:

        wrfseq (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. The default is 0.

        is_moving (:obj:`bool`): A boolean type that indicates if the
            sequence is a moving nest.

        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """
    if is_moving is None:
        is_moving = is_moving_domain(wrfseq, varname, _key=_key)

    file_times = extract_times(wrfseq, ALL_TIMES, meta=False, do_xtime=False)

    multitime = is_multi_time_req(timeidx)

    # For single times, just need to find the ncfile and appropriate
    # time index, and return that array
    if not multitime:
        return _find_arr_for_time(wrfseq, varname, timeidx, is_moving, meta,
                                  _key)

    # If all times are requested, need to build a new array and cat together
    # all of the arrays in the sequence
    wrf_iter = iter(wrfseq)

    if xarray_enabled() and meta:
        first_var = _build_data_array(next(wrf_iter), varname,
                                      ALL_TIMES, is_moving, True, _key)
    else:
        first_var = (next(wrf_iter)).variables[varname][:]

    outdims = [len(file_times)]

    # Making a new time dim, so ignore this one
    outdims += first_var.shape[1:]

    outdata = np.empty(outdims, first_var.dtype)

    numtimes = first_var.shape[0]
    startidx = 0
    endidx = numtimes

    if first_var.ndim > 1:
        outdata[startidx:endidx, :] = first_var[:]
    else:
        outdata[startidx:endidx] = first_var[:]

    if xarray_enabled() and meta:
        latname, lonname, timename = _find_coord_names(first_var.coords)

        timecached = False
        latcached = False
        loncached = False

        outxtimes = None
        outlats = None
        outlons = None

        timekey = timename+"_cat" if timename is not None else None
        latkey = latname + "_cat" if latname is not None else None
        lonkey = lonname + "_cat" if lonname is not None else None

        if timename is not None:
            outxtimes = get_cached_item(_key, timekey)
            if outxtimes is None:
                outxtimes = np.empty(outdims[0])
                outxtimes[startidx:endidx] = to_np(
                    first_var.coords[timename][:])
            else:
                timecached = True

        if is_moving:
            outcoorddims = outdims[0:1] + outdims[-2:]

            if latname is not None:
                # Try to pull from the coord cache
                outlats = get_cached_item(_key, latkey)
                if outlats is None:
                    outlats = np.empty(outcoorddims, first_var.dtype)
                    outlats[startidx:endidx, :] = to_np(
                        first_var.coords[latname][:])
                else:
                    latcached = True

            if lonname is not None:
                outlons = get_cached_item(_key, lonkey)
                if outlons is None:
                    outlons = np.empty(outcoorddims, first_var.dtype)
                    outlons[startidx:endidx, :] = to_np(
                        first_var.coords[lonname][:])
                else:
                    loncached = True

    startidx = endidx
    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            vardata = wrfnc.variables[varname][:]

            numtimes = vardata.shape[0]

            endidx = startidx + numtimes

            if vardata.ndim > 1:
                outdata[startidx:endidx, :] = vardata[:]
            else:
                outdata[startidx:endidx] = vardata[:]

            if xarray_enabled() and meta:
                if timename is not None and not timecached:
                    xtimedata = wrfnc.variables[timename][:]
                    outxtimes[startidx:endidx] = xtimedata[:]

                if is_moving:
                    if latname is not None and not latcached:
                        latdata = wrfnc.variables[latname][:]
                        outlats[startidx:endidx, :] = latdata[:]

                    if lonname is not None and not loncached:
                        londata = wrfnc.variables[lonname][:]
                        outlons[startidx:endidx, :] = londata[:]

            startidx = endidx

    if xarray_enabled() and meta:

        # Cache the coords if applicable
        if not latcached and outlats is not None:
            cache_item(_key, latkey, outlats)
        if not loncached and outlons is not None:
            cache_item(_key, lonkey, outlons)
        if not timecached and outxtimes is not None:
            cache_item(_key, timekey, outxtimes)

        outname = first_var.name
        outattrs = OrderedDict(first_var.attrs)
        outcoords = OrderedDict(first_var.coords)
        outdimnames = list(first_var.dims)

        if "Time" not in outdimnames:
            outdimnames.insert(0, "Time")

        if not multitime:
            file_times = [file_times[timeidx]]

        outcoords[outdimnames[0]] = file_times

        outcoords["datetime"] = outdimnames[0], file_times

        if timename is not None:
            outxtimes = outxtimes[:]
            outcoords[timename] = outdimnames[0], outxtimes

        # If the domain is moving, need to create the lat/lon coords
        # since they can't be copied
        if is_moving:
            outlatdims = [outdimnames[0]] + outdimnames[-2:]

            if latname is not None:
                outlats = outlats[:]
                outcoords[latname] = outlatdims, outlats
            if lonname is not None:
                outlons = outlons[:]
                outcoords[lonname] = outlatdims, outlons

        outdata = outdata[:]

        outarr = DataArray(outdata, name=outname, coords=outcoords,
                           dims=outdimnames, attrs=outattrs)

    else:
        outdata = outdata[:]
        outarr = outdata

    return outarr


def _get_numfiles(wrfseq):
    """Return the number of files in the sequence.

    This function will first try to call the builtin :meth:`len` function, but
    if that fails, the entire squence will be iterated over and counted.

    Args:

        wrfseq (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

    Returns:

        :obj:`int`: The number of files in the sequence.

    """
    try:
        return len(wrfseq)
    except TypeError:
        wrf_iter = iter(wrfseq)
        return sum(1 for _ in wrf_iter)


def _join_files(wrfseq, varname, timeidx, is_moving, meta, _key):
    """Return an array object from a sequence of files using the join
    method.

    The join method creates a new leftmost dimension for the file/sequence
    index.  In situations where there are multiple files with multiple times,
    and the last file contains less times than the previous files, the
    remaining arrays will be arrays filled with missing values.  There are
    checks in place within the wrf-python algorithms to look for these missing
    arrays, but be careful when calling compiled routines outside of
    wrf-python.

    In general, join is rarely used, so the concatenate method should be used
    for most cases.

    Args:

        wrfseq (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. The default is 0.

        is_moving (:obj:`bool`): A boolean type that indicates if the
            sequence is a moving nest.

        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """
    if is_moving is None:
        is_moving = is_moving_domain(wrfseq, varname, _key=_key)
    multitime = is_multi_time_req(timeidx)
    numfiles = _get_numfiles(wrfseq)
    maxtimes = _find_max_time_size(wrfseq)

    time_idx_or_slice = timeidx if not multitime else slice(None)
    file_times_less_than_max = False
    file_idx = 0

    # wrfseq might be a generator
    wrf_iter = iter(wrfseq)
    wrfnc = next(wrf_iter)
    numtimes = extract_dim(wrfnc, "Time")

    if xarray_enabled() and meta:
        first_var = _build_data_array(wrfnc, varname, ALL_TIMES, is_moving,
                                      True, _key)
        time_coord = np.full((numfiles, maxtimes), np.datetime64("NaT"),
                             "datetime64[ns]")
        time_coord[file_idx, 0:numtimes] = first_var.coords["Time"][:]
    else:
        first_var = wrfnc.variables[varname][:]

    if numtimes < maxtimes:
        file_times_less_than_max = True

    # Out dimensions will be the number of files, maxtimes, then the
    # non-time shapes from the first variable
    outdims = [numfiles]
    outdims += [maxtimes]
    outdims += first_var.shape[1:]

    # For join, always need to start with full masked values
    outdata = np.full(outdims, default_fill(first_var.dtype), first_var.dtype)
    if first_var.ndim > 1:
        outdata[file_idx, 0:numtimes, :] = first_var[:]
    else:
        outdata[file_idx, 0:numtimes] = first_var[:]

    # Create the secondary coordinate arrays
    if xarray_enabled() and meta:
        latname, lonname, timename = _find_coord_names(first_var.coords)
        outcoorddims = outdims[0:2] + outdims[-2:]

        timecached = False
        latcached = False
        loncached = False

        outxtimes = None
        outlats = None
        outlons = None

        timekey = timename+"_join" if timename is not None else None
        latkey = latname + "_join" if latname is not None else None
        lonkey = lonname + "_join" if lonname is not None else None

        if timename is not None:
            outxtimes = get_cached_item(_key, timekey)
            if outxtimes is None:
                outxtimes = np.full(outdims[0:2],
                                    default_fill(first_var.dtype),
                                    first_var.dtype)
                outxtimes[file_idx, 0:numtimes] = first_var.coords[timename][:]
            else:
                timecached = True

        if is_moving:
            if latname is not None:
                outlats = get_cached_item(_key, latkey)
                if outlats is None:
                    outlats = np.full(outcoorddims,
                                      default_fill(first_var.dtype),
                                      first_var.dtype)
                    outlats[file_idx, 0:numtimes, :] = (
                        first_var.coords[latname][:])
                else:
                    latcached = True

            if lonname is not None:
                outlons = get_cached_item(_key, lonkey)
                if outlons is None:
                    outlons = np.full(outcoorddims,
                                      default_fill(first_var.dtype),
                                      first_var.dtype)
                    outlons[file_idx, 0:numtimes, :] = (
                        first_var.coords[lonname][:])
                else:
                    loncached = True

    file_idx = 1
    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            numtimes = extract_dim(wrfnc, "Time")
            if numtimes < maxtimes:
                file_times_less_than_max = True
            outvar = wrfnc.variables[varname][:]

            if not multitime:
                outvar = outvar[np.newaxis, :]

            if outvar.ndim > 1:
                outdata[file_idx, 0:numtimes, :] = outvar[:]
            else:
                outdata[file_idx, 0:numtimes] = outvar[:]

            if xarray_enabled() and meta:
                # For join, the times are a function of fileidx
                file_times = extract_times(wrfnc, ALL_TIMES, meta=False,
                                           do_xtime=False)
                time_coord[file_idx, 0:numtimes] = np.asarray(
                    file_times, "datetime64[ns]")[:]

                if timename is not None and not timecached:
                    xtimedata = wrfnc.variables[timename][:]
                    outxtimes[file_idx, 0:numtimes] = xtimedata[:]

                if is_moving:
                    if latname is not None and not latcached:
                        latdata = wrfnc.variables[latname][:]
                        outlats[file_idx, 0:numtimes, :] = latdata[:]

                    if lonname is not None and not loncached:
                        londata = wrfnc.variables[lonname][:]
                        outlons[file_idx, 0:numtimes, :] = londata[:]

            # Need to update coords here
            file_idx += 1

    # If any of the output files contain less than the max number of times,
    # then a mask array is needed to flag all the missing arrays with
    # missing values
    if file_times_less_than_max:
        outdata = np.ma.masked_values(outdata, default_fill(outdata.dtype))

    if xarray_enabled() and meta:
        # Cache the coords if applicable
        if not latcached and outlats is not None:
            cache_item(_key, latkey, outlats)
        if not loncached and outlons is not None:
            cache_item(_key, lonkey, outlons)
        if not timecached and outxtimes is not None:
            cache_item(_key, timekey, outxtimes)

        outname = first_var.name
        outcoords = OrderedDict(first_var.coords)
        outattrs = OrderedDict(first_var.attrs)
        # New dimensions
        outdimnames = ["file"] + list(first_var.dims)
        outcoords["file"] = [i for i in py3range(numfiles)]

        # Time needs to be multi dimensional, so use the default dimension
        del outcoords["Time"]

        time_coord = time_coord[:, time_idx_or_slice]
        if not multitime:
            time_coord = time_coord[:, np.newaxis]
        outcoords["datetime"] = outdimnames[0:2], time_coord

        if isinstance(outdata, np.ma.MaskedArray):
            outattrs["_FillValue"] = default_fill(outdata.dtype)
            outattrs["missing_value"] = default_fill(outdata.dtype)

        if timename is not None:
            outxtimes = outxtimes[:, time_idx_or_slice]
            if not multitime:
                outxtimes = outxtimes[:, np.newaxis]
            outcoords[timename] = outdimnames[0:2], outxtimes[:]

        # If the domain is moving, need to create the lat/lon coords
        # since they can't be copied
        if is_moving:
            outlatdims = outdimnames[0:2] + outdimnames[-2:]

            if latname is not None:
                outlats = outlats[:, time_idx_or_slice, :]
                if not multitime:
                    outlats = outlats[:, np.newaxis, :]
                outcoords[latname] = outlatdims, outlats
            if lonname is not None:
                outlons = outlons[:, time_idx_or_slice, :]
                if not multitime:
                    outlons = outlons[:, np.newaxis, :]
                outcoords[lonname] = outlatdims, outlons

        if not multitime:
            outdata = outdata[:, timeidx, :]
            outdata = outdata[:, np.newaxis, :]

        outarr = DataArray(outdata, name=outname, coords=outcoords,
                           dims=outdimnames, attrs=outattrs)

    else:
        if not multitime:
            outdata = outdata[:, timeidx, :]
            outdata = outdata[:, np.newaxis, :]

        outarr = outdata

    return outarr


def combine_files(wrfin, varname, timeidx, is_moving=None,
                  method="cat", squeeze=True, meta=True,
                  _key=None):
    """Combine and return an array object for the sequence of WRF output
    files.

    Two aggregation methodologies are available to combine the sequence:

        - 'cat': Concatenate the files along the 'Time' dimension.  The Time
          dimension will be the leftmost dimension.  No sorting is performed,
          so files must be properly ordered in the sequence prior to calling
          this function.

        - 'join': Join the files by creating a new leftmost dimension for the
          file index. In situations where there are multiple files with
          multiple times, and the last file contains less times than the
          previous files, the remaining arrays will be arrays filled with
          missing values. There are checks in place within the wrf-python
          algorithms to look for these missing arrays, but be careful when
          calling compiled routines outside of wrf-python.


    Args:

        wrfin (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. The default is 0.

        is_moving (:obj:`bool`): A boolean type that indicates if the
            sequence is a moving nest.

        method (:obj:`str`, optional): The aggregation method to use for
            sequences.  Must be either 'cat' or 'join'.
            'cat' combines the data along the Time dimension.
            'join' creates a new dimension for the file index.
            The default is 'cat'.

        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """

    # Handles generators, single files, lists, tuples, custom classes
    wrfseq = get_iterable(wrfin)

    # Dictionary is unique
    if is_mapping(wrfseq):
        outarr = _combine_dict(wrfseq, varname, timeidx, method, meta, _key)
    elif method.lower() == "cat":
        outarr = _cat_files(wrfseq, varname, timeidx, is_moving,
                            squeeze, meta, _key)
    elif method.lower() == "join":
        outarr = _join_files(wrfseq, varname, timeidx, is_moving, meta, _key)
    else:
        raise ValueError("method must be 'cat' or 'join'")

    return outarr.squeeze() if squeeze else outarr


def _extract_var(wrfin, varname, timeidx, is_moving,
                 method, squeeze, cache, meta, _key):
    """Extract a variable from a NetCDF file object or a sequence of NetCDF
    file objects.

    Args:

        wrfin (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. The default is 0.

        is_moving (:obj:`bool`): A boolean type that indicates if the
            sequence is a moving nest.

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

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.

    """

    if cache is not None:
        try:
            cache_var = cache[varname]
        except KeyError:
            pass
        else:
            if not meta:
                return to_np(cache_var)

            return cache_var

    multitime = is_multi_time_req(timeidx)
    multifile = is_multi_file(wrfin)

    if is_time_coord_var(varname):
        return extract_times(wrfin, timeidx, method, squeeze, cache,
                             meta, do_xtime=True)

    if not multifile:
        if xarray_enabled() and meta:
            if is_moving is None:
                is_moving = is_moving_domain(wrfin, varname, _key=_key)
            result = _build_data_array(wrfin, varname, timeidx, is_moving,
                                       multifile, _key)
        else:
            if not multitime:
                result = wrfin.variables[varname][timeidx, :]
                result = result[np.newaxis, :]  # So that no squeeze works
            else:
                result = wrfin.variables[varname][:]
    else:
        # Squeeze handled in this routine, so just return it
        return combine_files(wrfin, varname, timeidx, is_moving,
                             method, squeeze, meta, _key)

    return result.squeeze() if squeeze else result


def extract_vars(wrfin, timeidx, varnames, method="cat", squeeze=True,
                 cache=None, meta=True, _key=None):
    """Extract variables from a NetCDF file object or a sequence of NetCDF
    file objects.

    Args:

        wrfin (iterable): An iterable type, which includes lists, tuples,
            dictionaries, generators, and user-defined classes.

        varnames (sequence of :obj:`str`) : A sequence of variable names.

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

        _key (:obj:`int`, optional): Cache key for the coordinate variables.
            This is used for internal purposes only.  Default is None.

    Returns:

        :obj:`dict`: A mapping of variable name to an array object. If xarray
        is enabled and the *meta* parameter is True, then the array object will
        be a :class:`xarray.DataArray` object.  Otherwise, the array object
        will be a :class:`numpy.ndarray` object with no metadata.

    """
    if isstr(varnames):
        varlist = [varnames]
    else:
        varlist = varnames

    return {var: _extract_var(wrfin, var, timeidx, None,
                              method, squeeze, cache, meta, _key)
            for var in varlist}


def npbytes_to_str(var):
    """Return a :obj:`bytes` object for the raw character array.

    Args:

        var (:class:`numpy.ndarray`): An array of characters.

    Returns:

        :obj:`bytes`: A string of bytes.

    """
    return (bytes(c).decode("utf-8") for c in var[:])


def _make_time(timearr):
    """Return a :class:`datetime.datetime` object for the array of characters.

    Args:

        timearr (:class:`numpy.ndarray`): An array of characters.

    Returns:

        :class:`datetime.datetime`: A datetime object.

    """
    try:
        return dt.datetime.strptime("".join(npbytes_to_str(timearr)),
                                    "%Y-%m-%d_%H:%M:%S")
    except ValueError:
        return np.datetime64("NaT")


def _file_times(wrfin, do_xtime):
    """Yield a time object for the times found in a sequence of files.

    If *do_xtime* to True, a :class:`datetime.datetime` object is yielded.
    Otherwise, a :obj:`float` object is yielded.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        do_xtime (:obj:`bool`): Set to True to parse the 'XTIME' variable
            instead of the 'Times' variable.

    Yields:

        :class:`datetime.datetime` or :obj:`float`: A
        :class:`datetime.datetime` object if *do_xtime* is False,
        otherwise a :obj:`float`.

    """
    if not do_xtime:
        times = wrfin.variables["Times"][:, :]
        for i in py3range(times.shape[0]):
            yield _make_time(times[i, :])
    else:
        xtimes = wrfin.variables["XTIME"][:]
        for i in py3range(xtimes.shape[0]):
            yield xtimes[i]


def _extract_time_map(wrfin, timeidx, do_xtime, meta=False):
    """Return a mapping of key to a sequence of time objects.

    This function is used when *wrfin* is a mapping.

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

        do_xtime (:obj:`bool`): Set to True to parse the 'XTIME' variable
            instead of the 'Times' variable.

        meta (:obj:`bool`, optional): Set to False to disable metadata.

    Returns:

        :obj:`dict`: A mapping of key to a sequence of time objects.  If
        *meta* is True, the sequence will be of type :class:`xarray.DataArray`,
        otherwise the sequence is :class:`numpy.ndarray`.

    """
    return {key: extract_times(wrfseq, timeidx, do_xtime, meta)
            for key, wrfseq in viewitems(wrfin)}


def extract_times(wrfin, timeidx, method="cat", squeeze=True, cache=None,
                  meta=False, do_xtime=False):

    """Return a sequence of time objects.

    If *do_xtime*  is False, the 'XTIME' variable is used and each time object
    is a :obj:`float`.  Otherwise, the 'Times' variable is used, and each
    time object is a :class:`datetime.datetime` object.

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

        do_xtime (:obj:`bool`): Set to True to parse the 'XTIME' variable
            instead of the 'Times' variable.  Default is False.

    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: A sequence of time
        objects.  If *meta* is True, the sequence will be of type
        :class:`xarray.DataArray`, otherwise the sequence is
        :class:`numpy.ndarray`.

    """
    if is_mapping(wrfin):
        return _extract_time_map(wrfin, timeidx, do_xtime)

    multitime = is_multi_time_req(timeidx)
    multi_file = is_multi_file(wrfin)
    if not multi_file:
        wrf_list = [wrfin]
    else:
        wrf_list = wrfin

    dt = "datetime64[ns]" if not do_xtime else np.float64
    fill_value = (np.datetime64('NaT') if not do_xtime else
                  default_fill(np.float64))

    try:
        if method.lower() == "cat":
            time_list = [file_time
                         for wrf_file in wrf_list
                         for file_time in _file_times(wrf_file, do_xtime)]
            time_arr = np.asarray(time_list, dtype=dt)

        elif method.lower() == "join":
            time_list = [[file_time
                          for file_time in _file_times(wrf_file, do_xtime)]
                         for wrf_file in wrf_list]

            num_rows = len(time_list)
            num_cols = len(time_list[0])

            time_arr = np.full((num_rows, num_cols), fill_value, dtype=dt)
            for i, row in enumerate(time_list):
                if len(row) == num_cols:
                    time_arr[i, :] = row[:]
                else:
                    for j, val in enumerate(row):
                        time_arr[i, j] = val

            time_arr = ma.masked_values(time_arr, fill_value)

        else:
            raise ValueError("invalid method argument '{}'".format(method))
    except KeyError:
        return None  # Thrown for pre-3.7 XTIME not existing

    if xarray_enabled() and meta:
        outattrs = OrderedDict()
        outcoords = None

        if method.lower() == "cat":
            outdimnames = ["Time"]
        else:
            outdimnames = ["fileidx", "Time"]
            outattrs["missing_value"] = fill_value
            outattrs["_FillValue"] = fill_value

        if not do_xtime:
            outname = "times"
            outattrs["description"] = "model times [np.datetime64]"

        else:
            ncfile = next(iter(wrf_list))
            var = ncfile.variables["XTIME"]
            outattrs.update(var.__dict__)

            outname = "XTIME"

        outarr = DataArray(time_arr, name=outname, coords=outcoords,
                           dims=outdimnames, attrs=outattrs)
    else:
        outarr = time_arr

    if not multitime:
        return outarr[timeidx]

    return outarr


def is_standard_wrf_var(wrfin, varname):
    """Return True if the variable is a standard WRF variable and not a
    diagnostic.

    If *wrfin* is a sequence, only the first file is used.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        varname (:obj:`str`): The variable name.

    Returns:

        :obj:`bool`: True if the variable is a standard WRF variable,
        otherwise False.

    """
    multifile = is_multi_file(wrfin)
    if multifile:
        if not is_mapping(wrfin):
            wrfin = next(iter(wrfin))
        else:
            entry = wrfin[next(iter(viewkeys(wrfin)))]
            return is_standard_wrf_var(entry, varname)

    return varname in wrfin.variables


def is_staggered(wrfin, var):
    """Return True if the variable is on a staggered grid.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        var (array): An array object which contains a :attr:`shape`
            attribute.

    Returns:

        :obj:`bool`: True if the variable is on a staggered grid, otherwise
        False.

    """

    we = extract_dim(wrfin, "west_east")
    sn = extract_dim(wrfin, "south_north")
    bt = extract_dim(wrfin, "bottom_top")

    if (var.shape[-1] != we or var.shape[-2] != sn or var.shape[-3] != bt):
        return True

    return False


def get_left_indexes(var, expected_dims):
    """Returns a tuple for the extra leftmost dimension sizes.

    For example, if an algorithm expects a 3 dimensional variable, but the
    variable includes an additional left dimension for Time, and
    this Time dimension has 3 values, then this function will return (3,).

    Args:

        var (array): An array object that contains the :attr:`ndim`
            and :attr:`shape` attributes.

        expected_dims (:obj:`int`): The expected number of dimensions (usually
            for a computational algorithm).

    Returns:

        :obj:`tuple`: The shape for the extra leftmost dimensions.

    """
    extra_dim_num = var.ndim - expected_dims

    if (extra_dim_num == 0):
        return []

    return var.shape[0:extra_dim_num]


def iter_left_indexes(dims):
    """Yield the iteration tuples for a sequence of dimensions sizes.

    For example, if *dims* is (2,2), then this will yield:

    (0,0), (0,1), (1,0), (1,1)

    This is primarily used to iterate over the leftmost index values.

    Args:

        dims (indexable sequence): A sequence of dimension sizes.

    Yields:

        :obj:`tuple`: The leftmost indexing iteration sizes.

    """
    arg = [py3range(dim) for dim in dims]
    for idxs in product(*arg):
        yield idxs


def get_right_slices(var, right_ndims, fixed_val=0):
    """Return an indexing tuple where the left dimensions are held to a
    fixed value and the right dimensions are set to slice objects.

    For example, if *var* is a 5D variable, and the desired indexing sequence
    for a numpy array is (0,0,0,:,:), then *right_ndims* should be set to 2
    and *fixed_val* set to 0.

    Args:

        var (:class:`numpy.ndarray`): A numpy array.

        right_ndims (:obj:`int`): The number of right dimensions to be sliced.

        fixed_val (:obj:`int`): The value to hold the left dimensions to.

    Returns:

        :obj:`tuple`: An indexing tuple that can be used to index a
        :class:`numpy.ndarray`.

    """
    extra_dim_num = var.ndim - right_ndims
    if extra_dim_num == 0:
        return [slice(None)] * right_ndims

    return tuple([fixed_val]*extra_dim_num +
                 [slice(None)]*right_ndims)


def get_proj_params(wrfin):
    """Return a tuple of latitude, longitude, and projection parameters from
    a WRF output file object or a sequence of WRF output file objects.

    Args:

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence.  Default is 0.

        varname (:obj:`str`, optional): The variable name to extract the
            coordinate variable names from.  Default is None, which will
            use the default coordinate variable names ('XLAT', 'XLONG').

    Returns:

        :obj:`tuple`: A tuple of the latitude coordinate variable,
        longitude coordinate, and global projection attributes.

    """
    proj_params = extract_global_attrs(wrfin,
                                       attrs=("MAP_PROJ",
                                              "CEN_LAT", "CEN_LON",
                                              "TRUELAT1", "TRUELAT2",
                                              "MOAD_CEN_LAT", "STAND_LON",
                                              "POLE_LAT", "POLE_LON",
                                              "DX", "DY"))

    return proj_params


def from_args(func, argnames, *args, **kwargs):
    """Return a mapping of argument name to value for the called function.

    This function parses the function args and kwargs to obtain the
    desired argument value. If the argument has not been passed in, the value
    is taken from the default keyword argument value.

    This func is usually called from within a decorator.

    Note:

        This function currently does not work with functions that contain
        variable length args or kwargs arguments.

    Args:

        func (function): The function to examine (usually the function that is
            wrapped).

        argnames (iterable): An iterable sequence of argument names.

        *args: The positional arguments.

        **kwargs: The keyword arguments.

    Returns:

        :obj:`dict`: A mapping of argument name to argument value.

    """
    if isstr(argnames):
        arglist = [argnames]
    else:
        arglist = argnames

    result = OrderedDict()
    for argname in arglist:
        arg_loc = arg_location(func, argname, args, kwargs)

        if arg_loc is not None:
            result[argname] = arg_loc[0][arg_loc[1]]
        else:
            result[argname] = None

    return result


def _args_to_list2(func, args, kwargs):
    """Return all of the function arguments, including defaults, as a list.

    The result can then be passed to the function via *result.  This version
    uses :meth:`inspect.argspec`, so is only applicable for Python 2.7.

    Note:

        This function currently does not work with functions that contain
        variable length args or kwargs arguments.

    Args:

        func (function): The function to examine (usually the function
            that is wrapped).

        args (:obj:`tuple`): The positional arguments.

        kwargs (:obj:`dict`):  The keyword arguments.

    Returns:

        :obj:`list`: A list of all argument values, including defaults.

    """
    argspec = getargspec(func)

    # Build the full tuple with defaults filled in
    outargs = [None]*len(argspec.args)
    if argspec.defaults is not None:
        for i, default in enumerate(argspec.defaults[::-1], 1):
            outargs[-i] = default

    # Add the supplied args
    for i, arg in enumerate(args):
        outargs[i] = arg

    # Fill in the supplied kargs
    for argname, val in viewitems(kwargs):
        argidx = argspec.args.index(argname)
        outargs[argidx] = val

    return outargs


# For Python 3.4, this will run the BoundArguments.apply_defaults method that
# was introduced in Python 3.5.
def _apply_defaults(bound):
    """Set default values for missing arguments.

    For variable-positional arguments (*args) the default is an
    empty tuple.

    For variable-keyword arguments (**kwargs) the default is an
    empty dict.

    Args:

        bound (:class:`inspect.BoundArguments`): A BoundArguments object.

    Returns:

        None

    """
    arguments = bound.arguments

    new_arguments = []
    for name, param in bound._signature.parameters.items():
        try:
            new_arguments.append((name, arguments[name]))
        except KeyError:
            if param.default is not _empty:
                val = param.default
            elif param.kind is _VAR_POSITIONAL:
                val = ()
            elif param.kind is _VAR_KEYWORD:
                val = {}
            else:
                # This BoundArguments was likely produced by
                # Signature.bind_partial().
                continue
            new_arguments.append((name, val))

    bound.arguments = OrderedDict(new_arguments)


def _args_to_list3(func, args, kwargs):
    """Return all of the function arguments, including defaults, as a list.

    The result can then be passed to the function via *result.  This version
    uses :meth:`inspect.signature`, so is only applicable for Python 3.4+.

    Note:

        This function currently does not work with functions that contain
        variable length args or kwargs arguments.

    Args:

        func (function): The function to examine (usually the function
            that is wrapped).

        args (:obj:`tuple`): The positional arguments.

        kwargs (:obj:`dict`):  The keyword arguments.

    Returns:

        :obj:`list`: A list of all argument values, including defaults.

    """
    sig = signature(func)
    bound = sig.bind(*args, **kwargs)
    try:
        bound.apply_defaults()
    except AttributeError:
        _apply_defaults(bound)

    return [x for x in bound.arguments.values()]


# Note:  Doesn't allow for **kwargs or *args
def args_to_list(func, args, kwargs):
    """Return all of the function arguments, including defaults, as a list.

    The result can then be passed to the function via *result*.

    Note:

        This function currently does not work with functions that contain
        variable length args or kwargs arguments.

    Args:

        func (function): The function to examine (usually the function
            that is wrapped).

        args (:obj:`tuple`): The positional arguments.

        kwargs (:obj:`dict`):  The keyword arguments.

    Returns:

        :obj:`list`: A list of all argument values, including defaults.

    """
    if version_info > (3,):
        _args_to_list = _args_to_list3
    else:
        _args_to_list = _args_to_list2

    return _args_to_list(func, args, kwargs)


def _arg_location2(func, argname, args, kwargs):
    """Return the function arguments as a single list along with the
    index within that list for a specified argument name.

    This function parses the args, kargs and signature looking for the
    location of *argname*, and returns a list containing all arguments, along
    with the argument location in that list.

    This function requires :meth:`inspect.getargspec`, so it is only
    applicable for Python 2.7.

    Args:

        func (function): The function to examine (usually the function
            that is wrapped).

        argname (:obj:`str`): The argument name to locate.

        args (:obj:`tuple`): The positional arguments.

        kwargs (:obj:`dict`):  The keyword arguments.

    Returns:

        :obj:`tuple`: A tuple containing the list of all argument values along
        with the index for location of *argname*.

    """
    argspec = getargspec(func)

    list_args = _args_to_list2(func, args, kwargs)

    # Return the new sequence and location
    if argname not in argspec.args and argname not in kwargs:
        return None

    result_idx = argspec.args.index(argname)

    return list_args, result_idx


def _arg_location3(func, argname, args, kwargs):
    """Return the function arguments as a single list along with the
    index within that list for a specified argument name.

    This function parses the args, kargs and signature looking for the
    location of *argname*, and returns a list containing all arguments, along
    with the argument location in that list.

    This function requires :meth:`inspect.signature`, so it is only
    applicable for Python 3.4 and higher.

    Args:

        func (function): The function to examine (usually the function
            that is wrapped).

        argname (:obj:`str`): The argument name to locate.

        args (:obj:`tuple`): The positional arguments.

        kwargs (:obj:`dict`):  The keyword arguments.

    Returns:

        :obj:`tuple`: A tuple containing the list of all argument values along
        with the index for location of *argname*.

    """
    sig = signature(func)
    params = list(sig.parameters.keys())

    list_args = _args_to_list3(func, args, kwargs)

    try:
        result_idx = params.index(argname)
    except ValueError:
        return None

    return list_args, result_idx


def arg_location(func, argname, args, kwargs):
    """Return the function arguments as a single list along with the
    index within that list for a specified argument name.

    This function parses the args, kargs and signature looking for the
    location of *argname*, and returns a list containing all arguments, along
    with the argument location in that list.

    Args:

        func (function): The function to examine (usually the function
            that is wrapped).

        argname (:obj:`str`): The argument name to locate.

        args (:obj:`tuple`): The positional arguments.

        kwargs (:obj:`dict`):  The keyword arguments.

    Returns:

        :obj:`tuple`: A tuple containing the list of all argument values along
        with the index for location of *argname*.

    """
    if version_info > (3,):
        _arg_location = _arg_location3
    else:
        _arg_location = _arg_location2

    return _arg_location(func, argname, args, kwargs)


def psafilepath():
    """Return the full path to the 'psadilookup.dat' file.

    The 'psadilookup.dat' file contains the lookup table for the cape
    routines.

    Returns:

        :obj:`str`:  The full path to the 'psadilookup.dat' file.

    """
    return os.path.join(os.path.dirname(__file__), "data", "psadilookup.dat")


def get_filepath(obj):
    """Return the file path for the specified object.

    This is used to return the file path for a netcdf object. If the
    particular object does not have the appropriate file path information,
    then one is created based on the timestep in the file.

    Args:

        obj: An object.

    Returns:

        :obj:`str`: A string for a file path.

    """
    try:
        path = obj.filepath()
    except AttributeError:
        try:
            path = obj.file.path
        except AttributeError:
            # Let's make up a filename from the first file time
            found = False
            times = extract_times(obj, None, meta=False, do_xtime=False)
            for t in times:
                path = "wrfout_{}".format(str(t))
                found = True
                break

            if not found:
                raise ValueError("file contains no path information")

    return path


def get_id(obj, prefix=''):
    """Return the cache id.

    The cache id is used as a caching key for various routines. If the
    object type is a mapping, then the result will also be a
    mapping of each key to the object id for the value.

    Args:

        obj (:obj:`object`): Any object type.

        prefix (:obj:`str`): A string to help with recursive calls.

    Returns:

        :obj:`int` or :obj:`dict`: If the *obj* parameter is not a mapping,
        then the object id is returned.  Otherwise, a mapping of each
        key to the object id for the value is returned.

    """
    if not is_multi_file(obj):
        return hash(prefix + get_filepath(obj))

    # For sequences, the hashing string will be the list ID and the
    # path for the first file in the sequence
    if not is_mapping(obj):
        _obj = get_iterable(obj)
        _next = next(iter(_obj))
        return get_id(_next, prefix + str(id(obj)))

    # For each key in the mapping, recursively call get_id until
    # until a non-mapping is found
    return {key: get_id(val, prefix) for key, val in viewitems(obj)}


def geo_bounds(var=None, wrfin=None, varname=None, timeidx=0, method="cat",
               squeeze=True, cache=None):
    """Return the geographic boundaries for the variable or file(s).

    When using a :class:`xarray.DataArray` as the *var* parameter, the variable
    must contain latitude and longitude coordinates.  If these coordinate
    dimensions are greater than two dimensions, then an array of
    :class:`wrf.GeoBounds` objects will be returned with the same shape as the
    leftmost dimensions of the coordinate arrays.

    When using a WRF file, or sequence of WRF files, by supplying the
    *wrfin* parameter, an array of :class:`wrf.GeoBounds` objects will be
    returned if the domain is moving and :data:`wrf.ALL_TIMES` is selected as
    the *timeidx* parameter when using *wrfin*. Otherwise, a single
    :class:`wrf.GeoBounds` object is returned.

    Args:

        var (:class:`xarray.DataArray`, optional): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`, optional): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index when *wrfin* is not None. This value can be a
            positive integer, negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0. This value is
            ignored when *var* is used.

        method (:obj:`str`, optional): The aggregation method to use for
            sequences when *wrfin* is not None.  Must be either 'cat' or
            'join'.  'cat' combines the data along the Time dimension.
            'join' creates a new dimension for the file index.
            The default is 'cat'.

        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Only used when *wrfin* is used. Default is True.

        cache (:obj:`dict`, optional): A dictionary of (varname, ndarray)
            that can be used to supply pre-extracted NetCDF variables to the
            computational routines.  It is primarily used for internal
            purposes, but can also be used to improve performance by
            eliminating the need to repeatedly extract the same variables
            used in multiple diagnostics calculations, particularly when using
            large sequences of files. Only used when *wrfin* is used.
            Default is None.

    Returns:

        :class:`wrf.GeoBounds`:  The domain geographic bounds.

    """

    if var is None and wrfin is None:
        raise ValueError("'var' or 'wrfin' parameter is required")

    # Getting lat/lon from xarray coordinates
    if var is not None:
        if not xarray_enabled():
            raise ValueError("xarray is not installed or is disabled")

        is_moving = None
        try:
            var_coords = var.coords
        except AttributeError:
            raise ValueError("'var' object does not contain coordinate "
                             "attributes")

        latname, lonname, _ = _find_coord_names(var_coords)
        try:
            lats = to_np(var_coords[latname])
        except KeyError:
            raise ValueError("'var' object does not contain a latitude "
                             "coordinate")
        try:
            lons = to_np(var_coords[lonname])
        except KeyError:
            raise ValueError("'var' object does not contain a longitude "
                             "coordinate")

    # Getting lat/lon from the file
    elif wrfin is not None:
        _key = get_id(wrfin)
        is_moving = is_moving_domain(wrfin, varname=varname,
                                     latvar=either("XLAT", "XLAT_M"),
                                     lonvar=either("XLONG", "XLONG_M"),
                                     _key=_key)
        if varname is not None:
            if xarray_enabled():
                var = extract_vars(wrfin, timeidx, varname, method, squeeze,
                                   cache, meta=True, _key=_key)[varname]
                return geo_bounds(var)
            else:
                lat_coord, lon_coord, _ = _get_coord_names(wrfin, varname)
        else:
            lat_coord = either("XLAT", "XLAT_M")(wrfin)
            lon_coord = either("XLONG", "XLONG_M")(wrfin)

        # If requesting all times but the domain isn't moving, just
        # extract one time
        _timeidx = timeidx
        if timeidx is None and not is_moving:
            _timeidx = 0

        coord_data = extract_vars(wrfin, _timeidx, (lat_coord, lon_coord),
                                  method, squeeze, cache, meta=False,
                                  _key=_key)
        lats = coord_data[lat_coord]
        lons = coord_data[lon_coord]

    # Moving domains
    if lats.ndim > 2:
        # Requesting all times, works for 'cat' and 'join' data
        # and always for xarray data
        extra_dims = lats.shape[0:-2]
        out_geo = np.full(extra_dims, NullGeoBounds(), np.object)

        for left_idxs in iter_left_indexes(extra_dims):
            latlon_idx = left_idxs + (slice(None),)
            out_geo[left_idxs] = GeoBounds(lats=lats[latlon_idx],
                                           lons=lons[latlon_idx])
        return out_geo

    # Non-moving domains
    return GeoBounds(lats=lats, lons=lons)


def _get_wrf_proj_geobnds(var, wrfin, varname, timeidx, method, squeeze,
                          cache):
    """Return the :class:`wrf.WrfProj` subclass and :class:`wrf.GeoBounds`.

    Args:

        var (:class:`xarray.DataArray`): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0.

        method (:obj:`str`): The aggregation method to use for
            sequences.  Must be either 'cat' or 'join'.
            'cat' combines the data along the Time dimension.
            'join' creates a new dimension for the file index.
            The default is 'cat'.

        squeeze (:obj:`bool`): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

        cache (:obj:`dict`): A dictionary of (varname, ndarray)
            that can be used to supply pre-extracted NetCDF variables to the
            computational routines.  It is primarily used for internal
            purposes, but can also be used to improve performance by
            eliminating the need to repeatedly extract the same variables
            used in multiple diagnostics calculations, particularly when using
            large sequences of files.
            Default is None.

    Returns:

        :obj:`tuple`: A tuple of :class:`wrf.WrfProj`
        and :class:`wrf.GeoBounds`

    """
    # Using a variable
    if var is not None:
        if not xarray_enabled():
            raise ValueError("xarray is not installed or is disabled")

        geobnds = geo_bounds(var)
        try:
            wrf_proj = var.attrs["projection"]
        except AttributeError:
            raise ValueError("variable does not contain projection "
                             "information")
    else:
        geobnds = geo_bounds(wrfin=wrfin, varname=varname, timeidx=timeidx,
                             method=method, cache=cache)
        proj_params = get_proj_params(wrfin)
        wrf_proj = getproj(**proj_params)

    return wrf_proj, geobnds


def _get_proj_obj(ob_type, var, wrfin, varname, timeidx, method, squeeze,
                  cache, **kwargs):
    """Return the desired mapping object for the plotting type.

    Args:

        ob_type (:obj:`str`): Must be 'cartopy', 'basemap', or 'pyngl'.

        var (:class:`xarray.DataArray`): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0.

        method (:obj:`str`): The aggregation method to use for
            sequences.  Must be either 'cat' or 'join'.
            'cat' combines the data along the Time dimension.
            'join' creates a new dimension for the file index.
            The default is 'cat'.

        squeeze (:obj:`bool`): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

        cache (:obj:`dict`): A dictionary of (varname, ndarray)
            that can be used to supply pre-extracted NetCDF variables to the
            computational routines.  It is primarily used for internal
            purposes, but can also be used to improve performance by
            eliminating the need to repeatedly extract the same variables
            used in multiple diagnostics calculations, particularly when using
            large sequences of files.
            Default is None.

        **kwargs: Additional keyword arguments for the plotting type.

    Returns:

        mapping object:  Will be either :class:`cartopy.crs.Projection`,
        :class:`matplotlib.mpl_toolkits.basemap.Basemap` or
        :class:`Ngl.Resources`.

    """

    wrf_proj, geobnds = _get_wrf_proj_geobnds(var, wrfin, varname, timeidx,
                                              method, squeeze, cache)

    if ob_type == "cartopy":
        proj_obj = wrf_proj.cartopy()
    elif ob_type == "basemap":
        try:
            _ = len(geobnds)
        except TypeError:  # Only a single object
            proj_obj = wrf_proj.basemap(geobnds, **kwargs)
        else:
            proj_obj = np.empty(geobnds.shape, np.object)

            for idxs, geobnd_val in np.ndenumerate(geobnds):
                proj_obj[idxs] = wrf_proj.basemap(geobnd_val, **kwargs)
    elif ob_type == "pyngl":
        try:
            _ = len(geobnds)
        except TypeError:  # Only a single object
            proj_obj = wrf_proj.pyngl(geobnds, **kwargs)
        else:
            proj_obj = np.empty(geobnds.shape, np.object)

            for idxs, geobnd_val in np.ndenumerate(geobnds):
                proj_obj[idxs] = wrf_proj.pyngl(geobnd_val, **kwargs)

    return proj_obj


def latlon_coords(var, as_np=False):
    """Return the latitude and longitude coordinates from a
    :class:`xarray.DataArray` object.

    Args:

        var (:class:`xarray.DataArray`):  A variable.

        as_np (:obj:`bool`): Set to True to return the coordinates as
            :class:`numpy.ndarray` objects instead of
            :class:`xarray.DataArray` objects.

    Returns:

        :obj:`tuple`: The latitude and longitude coordinate variables.

    """

    if not xarray_enabled():
        raise ValueError("xarray is not installed or is disabled")

    try:
        var_coords = var.coords
    except AttributeError:
        raise ValueError("'var' object does not contain coordinate "
                         "attributes")

    latname, lonname, _ = _find_coord_names(var_coords)
    try:
        lats = var_coords[latname]
    except KeyError:
        raise ValueError("'var' object does not contain a latitude "
                         "coordinate")
    try:
        lons = var_coords[lonname]
    except KeyError:
        raise ValueError("'var' object does not contain a longitude "
                         "coordinate")

    if as_np:
        return to_np(lats), to_np(lons)

    return lats, lons


def get_cartopy(var=None, wrfin=None, varname=None, timeidx=0, method="cat",
                squeeze=True, cache=None):
    """Return a :class:`cartopy.crs.Projection` subclass for the
    map projection.

    Args:

        var (:class:`xarray.DataArray`, optional): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        geobounds (:class:`wrf.GeoBounds`, optional):  The geobounds to
            get the extents.  If set to None and using the *var* parameter,
            the geobounds will be taken from the variable.  If using a
            file, then the geobounds will be taken from the native grid.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`, optional): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0.

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

    Returns:

        :class:`cartopy.crs.Projection`: A Projection subclass for the
        map projection.

    See Also:

        :class:`cartopy.crs.Projection`

    """
    return _get_proj_obj("cartopy", var, wrfin, varname, timeidx, method,
                         squeeze, cache)


def get_basemap(var=None, wrfin=None, varname=None, timeidx=0, method="cat",
                squeeze=True, cache=None, **kwargs):
    """Return a :class:`matplotlib.mpl_toolkits.basemap.Basemap` object
        for the map projection.

    Args:

        var (:class:`xarray.DataArray`, optional): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        geobounds (:class:`wrf.GeoBounds`, optional):  The geobounds to
            get the extents.  If set to None and using the *var* parameter,
            the geobounds will be taken from the variable.  If using a
            file, then the geobounds will be taken from the native grid.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`, optional): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0.

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

        **kwargs: Keyword arguments for creating a
            :class:`matplotlib.mpl_toolkits.basemap.Basemap`.  By default,
            the domain bounds will be set to the native projection, the
            resolution will be set to 'l', and the other projection
            parameters will be set by the information in the file.

    Returns:

        :class:`cartopy.crs.Projection`: A Projection subclass for the
        map projection.

    Returns:

            :class:`matplotlib.mpl_toolkits.basemap.Basemap`: A Basemap
            object for the projection.

        See Also:

            :class:`matplotlib.mpl_toolkits.basemap.Basemap`

    """
    return _get_proj_obj("basemap", var, wrfin, varname, timeidx, method,
                         squeeze, cache, **kwargs)


def get_pyngl(var=None, wrfin=None, varname=None, timeidx=0, method="cat",
              squeeze=True, cache=None, **kwargs):
    """Return a :class:`Ngl.Resources` object for the map projection.

    Args:

        var (:class:`xarray.DataArray`, optional): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        geobounds (:class:`wrf.GeoBounds`, optional):  The geobounds to
            get the extents.  If set to None and using the *var* parameter,
            the geobounds will be taken from the variable.  If using a
            file, then the geobounds will be taken from the native grid.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`, optional): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0.

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

        **kwargs: Additional PyNGL resources to set while creating the
            :class:`Ngl.Resources` object.

    Returns:

        :class:`Ngl.Resources`: A dict-like object that contains the
        PyNGL resources for the map projection.

    See Also:

        `PyNGL <https://www.pyngl.ucar.edu/>`_
    """
    return _get_proj_obj("pyngl", var, wrfin, varname, timeidx, method,
                         squeeze, cache)


def cartopy_xlim(var=None, geobounds=None, wrfin=None, varname=None, timeidx=0,
                 method="cat", squeeze=True, cache=None):
    """Return the x-axis limits in the projected coordinates.

    For some map projections, like :class`wrf.RotatedLatLon`, the
    :meth:`cartopy.GeoAxes.set_extent` method does not work correctly.  This
    method is equivalent to:

    .. code-block:: python

        pc = crs.PlateCarree()
        xs, ys, _  = self._cartopy().transform_points(pc,
                             np.array([geobounds.bottom_left.lon,
                                       geobounds.top_right.lon]),
                             np.array([geobounds.bottom_left.lat,
                                       geobounds.top_right.lat])).T


        _xlimits = xs.tolist()
        _ylimits = ys.tolist()

        return (_xlimits, _ylimits)[0]

    Args:

        var (:class:`xarray.DataArray`, optional): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        geobounds (:class:`wrf.GeoBounds`, optional):  The geobounds to
            get the extents.  If set to None and using the *var* parameter,
            the geobounds will be taken from the variable.  If using a
            file, then the geobounds will be taken from the native grid.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`, optional): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0.

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

    Returns:

        :obj:`list`:  A list of [start_x, end_x] in the projected coordinate
        system.

    """
    wrf_proj, native_geobnds = _get_wrf_proj_geobnds(var, wrfin, varname,
                                                     timeidx, method, squeeze,
                                                     cache)
    if geobounds is not None:
        return wrf_proj.cartopy_xlim(geobounds)

    return wrf_proj.cartopy_xlim(native_geobnds)


def cartopy_ylim(var=None, geobounds=None, wrfin=None, varname=None, timeidx=0,
                 method="cat", squeeze=True, cache=None):
    """Return the y-axis limits in the projected coordinates.

    For some map projections, like :class`wrf.RotatedLatLon`, the
    :meth:`cartopy.GeoAxes.set_extent` method does not work correctly.  This
    method is equivalent to:

    .. code-block:: python

        pc = crs.PlateCarree()
        xs, ys, _  = self._cartopy().transform_points(pc,
                             np.array([geobounds.bottom_left.lon,
                                       geobounds.top_right.lon]),
                             np.array([geobounds.bottom_left.lat,
                                       geobounds.top_right.lat])).T


        _xlimits = xs.tolist()
        _ylimits = ys.tolist()

        return (_xlimits, _ylimits)[1]

    Args:

        var (:class:`xarray.DataArray`, optional): A :class:`xarray.DataArray`
            variable that includes latitude,longitude coordinate information.
            If not used, then *wrfin* must be provided.

        geobounds (:class:`wrf.GeoBounds`, optional):  The geobounds to
            get the extents.  If set to None and using the *var* parameter,
            the geobounds will be taken from the variable.  If using a
            file, then the geobounds will be taken from the native grid.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. If not used,
            then *var* must be provided.

        varname (:obj:`str`, optional): If using *wrfin*, then this will be the
            variable name to use to determine the geobounds.  The variable
            can be a coordinate variable, or a regular variable that contains
            coordinate attributes.  If None,
            then the 'XLAT', 'XLAT_M', 'XLONG', 'XLONG_M' variables
            will be used.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Default is 0.

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

    Returns:

        :obj:`list`:  A list of [start_y, end_y] in the projected coordinate
        system.

    """
    wrf_proj, native_geobnds = _get_wrf_proj_geobnds(var, wrfin, varname,
                                                     timeidx, method, squeeze,
                                                     cache)
    if geobounds is not None:
        return wrf_proj.cartopy_ylim(geobounds)

    return wrf_proj.cartopy_ylim(native_geobnds)


def ll_points(lat, lon):
    """Return the lower left latitude and longitude point(s).

    This functions extracts the lower left corner points and returns the result
    as either a single :class:`CoordPair` object, or a list of
    :class:`CoordPair` objects.

    This is primarily used for testing or constructing the corner point objects
    from the XLAT and XLONG variables.

    Args:

        lat (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The latitude
            array. Must be at least two dimensions.

        lon (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The
            longitude array. Must be at least two dimensions.

    Returns:

        :class:`wrf.CoordPair` or :obj:`list`: A single :class:`wrf.CoordPair`
        object or a list of :class:`wrf.CoordPair` objects.

    """
    latvals = np.ravel(to_np(lat)[..., 0, 0])
    lonvals = np.ravel(to_np(lon)[..., 0, 0])

    if latvals.shape[0] == 1:
        return CoordPair(lat=float(latvals), lon=float(lonvals))
    else:
        return [CoordPair(lat=latvals[i], lon=lonvals[i])
                for i in py3range(latvals.shape[0])]


def pairs_to_latlon(pairs):
    """Return latitude and longitude arrays from a sequence of \
    :class:`wrf.CoordPair` objects.

    This function converts a sequence of :class:`wrf.CoordPair` objects into
    lists of latitude and longitude points. If the *pairs* argument is a
    single :class:`wrf.CoordPair` object, then a single latitude and
    longitude value is returned.

    Args:

        pairs (:class:`wrf.CoordPair` or sequence): A single
            :class:`wrf.CoordPair` or sequence of :class:`wrf.CoordPair`.

    Returns:

        :obj:`tuple`: A tuple of (lat, lon), where lat and lon are single
        values or lists of values.

    """

    if isinstance(pairs, CoordPair):
        return (pairs.lat, pairs.lon)
    else:
        lats = [pair.lat for pair in pairs]
        lons = [pair.lon for pair in pairs]

        return lats, lons


def is_latlon_pair(pair):
    """Return True if the :class:`wrf.CoordPair` is a lat/lon pair

    Args:

        pair (:class:`wrf.CoordPair`): A single :class:`wrf.CoordPair` object.

    Returns:

        :obj:`bool`: True if the pair is a lat/lon pair.

    """
    if pair is not None:
        return (pair.lat is not None and pair.lon is not None)
    else:
        return False
