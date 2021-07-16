from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict
from collections.abc import Iterable

import wrapt
import numpy as np
import numpy.ma as ma

from .units import do_conversion, check_units, dealias_and_clean_unit
from .util import iter_left_indexes, from_args, to_np, combine_dims
from .py3compat import viewitems, viewvalues, isstr
from .config import xarray_enabled
from .constants import default_fill

if xarray_enabled():
    from xarray import DataArray


def convert_units(unit_type, alg_unit):
    """A decorator that converts the units from the wrapped function's output.

    The desired units are determined from the wrapped function's arguments.

    Args:

        unit_type (:obj:`str`): The unit type.  Choices are: 'wind',
            'pressure', 'temp', or 'height'.

        alg_unit (:obj:`str`): The units returned by the wrapped function,
            which is usually the units returned by the Fortran routine.

    Returns:

        :class:`numpy.ndarray`: The wrapped function's output in the desired
        units.

    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        desired_units = from_args(wrapped, "units", *args, **kwargs)["units"]
        u_cleaned = dealias_and_clean_unit(desired_units)
        check_units(u_cleaned, unit_type)

        # Unit conversion done here
        return do_conversion(wrapped(*args, **kwargs), unit_type,
                             alg_unit, desired_units)

    return func_wrapper


def left_iteration(ref_var_expected_dims,
                   ref_var_right_ndims,
                   insert_dims=None,
                   ref_var_idx=None,
                   ref_var_name=None,
                   ignore_args=None,
                   ignore_kargs=None,
                   outviews="outview",
                   alg_dtype=np.float64,
                   cast_output=True):
    """A decorator to handle iterating over the leftmost dimensions.

    For example, if a wrapped function works with three-dimensional arrays, but
    the variables include a 4th leftmost dimension for 'Time', this decorator
    will iterate over all times, call the 3D Fortran routine, and aggregate the
    results in to a 4D output array.

    It is also important to note that the final output array is allocated
    first, and then views are passed to the wrapped function so that values
    do not need to get copied in to the final output array.

    Args:

        ref_var_expected_dims (:obj:`int`): The number of dimensions that the
            Fortran routine is expecting for the reference variable.

        ref_var_right_ndims (:obj:`int`): The number of dimensions from the
            right to keep for the reference variable when making the output.
            Can also be a :class:`combine_dims` object if the sizes are
            determined from multiple variables.

        insert_dims (sequence of :obj:`int`, optional): A sequence of
            dimensions to insert between the left dimensions (e.g. time) and
            the kept right dimensions. Default is None.

        ref_var_idx (:obj:`int`, optional): The index in the wrapped function's
            positional arguments to be used as the reference variable for
            determining the leftmost dimensions. Must be specified if
            *ref_var_name* is None.  Default is None.

        ref_var_name (:obj:`str`, optional): The keyword argument name for the
            wrapped function's keyword arguments to be used as the reference
            variable for calculating the leftmost dimensions.  Must be
            specified if *ref_var_idx* is None.  Default is None.

        ignore_args (sequence of :obj:`int`): Indexes of any arguments that
            should be ignored when creating the sliced views that are
            passed to the Fortran routine.

        ignore_kargs (sequence of :obj:`str`): Keys of any keyword arguments
            that should be ignored when creating the sliced views that are
            passed to the Fortran routine.

        outviews (:obj:`str` or a sequence): A single key or sequence of keys
            that indicate the wrapped function's keyword argument to use
            as the output variable(s) in the wrapped function.

        alg_dtype (:class:`numpy.dtype` or :obj:`str`): The numpy data type
            used in the wrapped function.

        cast_output (:obj:`bool`): Set to True to cast the wrapped function's
            output to the same type as the reference variable.

    Returns:

        :class:`numpy.ndarray`: The aggregated output array that includes
        all extra leftmost dimensions found in the reference variable.

    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        _ignore_args = ignore_args if ignore_args is not None else ()
        _ignore_kargs = ignore_kargs if ignore_kargs is not None else ()
        _outkeys = [outviews] if isstr(outviews) else outviews

        if ref_var_idx is not None:
            ref_var = args[ref_var_idx]
        else:
            ref_var = kwargs[ref_var_name]

        ref_var_dtype = ref_var.dtype
        ref_var_shape = ref_var.shape
        extra_dim_num = ref_var.ndim - ref_var_expected_dims

        # No special left side iteration, return the function result
        if (extra_dim_num == 0):
            return wrapped(*args, **kwargs)

        # Start by getting the left-most 'extra' dims
        extra_dims = ref_var_shape[0:extra_dim_num]

        mid_dims = () if insert_dims is None else tuple(insert_dims)

        if not isinstance(ref_var_right_ndims, combine_dims):
            right_dims = ref_var_shape[-ref_var_right_ndims:]
        else:
            right_dims = ref_var_right_ndims(*args)

        left_dims = extra_dims

        outdims = left_dims + mid_dims + right_dims

        if "outview" not in kwargs:
            outd = OrderedDict((outkey, np.empty(outdims, alg_dtype))
                               for outkey in _outkeys)

        mask_output = False
        for left_idxs in iter_left_indexes(extra_dims):
            # Make the left indexes plus a single slice object
            # The single slice will handle all the dimensions to
            # the right (e.g. [1,1,:])
            left_and_slice_idxs = left_idxs + (slice(None), )

            # Slice the args if applicable
            new_args = [arg[left_and_slice_idxs]
                        if i not in _ignore_args else arg
                        for i, arg in enumerate(args)]

            # Slice the kwargs if applicable
            new_kargs = {key: (val[left_and_slice_idxs]
                         if key not in _ignore_kargs else val)
                         for key, val in viewitems(kwargs)}

            # Skip the possible empty/missing arrays for the join method
            skip_missing = False
            for arg in new_args:
                try:
                    _ = arg.ndim
                except AttributeError:
                    continue  # Not an array object
                else:
                    arr = to_np(arg)

                try:
                    all_masked = arr.mask.all()
                except AttributeError:
                    pass  # Not a masked array
                else:
                    if all_masked:
                        for output in viewvalues(outd):
                            output[left_and_slice_idxs] = (
                                default_fill(np.float64))
                        skip_missing = True
                        mask_output = True
                        break

            if skip_missing:
                continue

            # Insert the output views if one hasn't been provided
            if "outview" not in new_kargs:
                for outkey, output in viewitems(outd):
                    outview = output[left_and_slice_idxs]
                    new_kargs[outkey] = outview

            result = wrapped(*new_args, **new_kargs)

            # Make sure the result is the same data as what got passed in
            # Can delete this once everything works
            if (result.__array_interface__["data"][0] !=
                    outview.__array_interface__["data"][0]):
                raise RuntimeError("output array was copied")

        if len(outd) == 1:
            output = next(iter(viewvalues(outd)))
        else:
            output = tuple(arr for arr in viewvalues(outd))

        if cast_output:
            if isinstance(output, np.ndarray):
                output = output.astype(ref_var_dtype)
            else:
                output = tuple(arr.astype(ref_var_dtype) for arr in output)

        # Mostly when used with join
        if mask_output:
            if isinstance(output, np.ndarray):
                output = ma.masked_values(output, default_fill(np.float64))
            else:
                output = tuple(ma.masked_values(arr, default_fill(np.float64))
                               for arr in output)

        return output

    return func_wrapper


def cast_type(ref_idx=0, arg_idxs=None, karg_names=None,
              alg_dtype=np.float64, outviews="outview"):
    """A decorator to handle type casting.

    This decorator is used to cast variables to and from the required
    :class:`numpy.dtype` used in the wrapped function.

    Args:

        ref_idx (:obj:`int`, optional): The index in the wrapped function's
            positional arguments to be used as the reference variable for
            determining the :class:`numpy.dtype` to return.  Default is 0.

        arg_idxs (sequence of :obj:`int`, optional): A sequence of indexes in
            the wrapped function's positional arguments that indicate which
            arguments to cast.  Must be specified if *karg_names* is None.
            Default is None.

        karg_names (sequence of :obj:`str`): A sequence of keyword arguments
            in the wrapped function's keyword arguments that indicate the
            arguments to cast.  Must be specified if *arg_idxs* is None.
            Default is None.

        alg_dtype (:class:`numpy.dtype` or :obj:`str`): The numpy data type
            used in the wrapped function.

        outviews (:obj:`str` or a sequence): A single key or sequence of keys
            that indicate the wrapped function's keyword argument to use
            as the output variable(s) in the wrapped function.

    Returns:

        :class:`numpy.ndarray`: The wrapped function's output cast to the
        same :class:`numpy.dtype` as the reference variable.

    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        _arg_idxs = arg_idxs if arg_idxs is not None else ()
        _karg_names = karg_names if karg_names is not None else ()

        # Handle output views if applicable
        _outkeys = [outviews] if isstr(outviews) else outviews
        _outviews = from_args(wrapped, _outkeys, *args, **kwargs)

        has_outview = False
        for outkey in _outkeys:
            _outview = _outviews[outkey]
            if _outview is not None:
                has_outview = True

        orig_type = args[ref_idx].dtype

        new_args = [arg.astype(alg_dtype)
                    if i in _arg_idxs else arg
                    for i, arg in enumerate(args)]

        new_kargs = {key: (val.astype(alg_dtype)
                           if key in _karg_names else val)
                     for key, val in viewitems(kwargs)}

        result = wrapped(*new_args, **new_kargs)

        # Do nothing for supplied output views
        if not has_outview:
            if isinstance(result, np.ndarray):
                if result.dtype == orig_type:
                    return result
                return result.astype(orig_type)
            elif isinstance(result, Iterable):  # got back a sequence of arrays
                return tuple(arr.astype(orig_type)
                             if arr.dtype != orig_type else arr
                             for arr in result)

        return result

    return func_wrapper


def _extract_and_transpose(arg, do_transpose):
    """Return a transposed view of the :class:`numpy.ndarray` inside of a
    :class:`xarray.DataArray` object.

    If the *arg* parameter is not a :class:`xarray.DataArray` object, then
    *arg* is returned.

    Args:

        arg (:class:`xarray.DataArray` or :obj:`object`): Can be any object
            type.

        do_transpose: Set to False to only extract the variable.  When True,
            the extracted array will also be transposed to a Fortran view if
            it is not already Fortran contiguous.

    Returns:

        :class:`numpy.ndarray`: A numpy array.  If *do_transpose* is True,
        the numpy array will also be a Fortran contiguous view.

    """

    if xarray_enabled():
        if isinstance(arg, DataArray):
            arg = to_np(arg)

    if do_transpose:
        if isinstance(arg, np.ndarray):
            if not arg.flags.f_contiguous and arg.ndim > 1:
                return arg.T

    return arg


def extract_and_transpose(do_transpose=True, outviews="outview"):
    """A decorator to extract the data array from a :class:`xarray.DataArray`

    This decorator also transposes the view of the data to Fortran
    contiguous if *do_transpose* is True.

    Args:

        do_transpose: Set to False to only extract the variable.  When True,
            the extracted array will also be transposed to a Fortran view if
            it is not already Fortran contiguous.

        outviews (:obj:`str` or a sequence): A single key or sequence of keys
            that indicate the wrapped function's keyword argument to use
            as the output variable(s) in the wrapped function.

    Returns:

        :class:`numpy.ndarray`: A numpy array.  If *do_transpose* is True,
        the numpy array will also be a Fortran contiguous view.

    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):

        # Handle output views if applicable
        _outkeys = [outviews] if isstr(outviews) else outviews
        _outviews = from_args(wrapped, _outkeys, *args, **kwargs)

        has_outview = False
        for outkey in _outkeys:
            _outview = _outviews[outkey]
            if _outview is not None:
                has_outview = True

        new_args = [_extract_and_transpose(arg, do_transpose) for arg in args]

        new_kargs = {key: _extract_and_transpose(val, do_transpose)
                     for key, val in viewitems(kwargs)}

        result = wrapped(*new_args, **new_kargs)

        # Do nothing for supplied output views
        if has_outview:
            return result

        if isinstance(result, np.ndarray):
            if result.flags.f_contiguous and result.ndim > 1:
                return result.T
        elif isinstance(result, Iterable):
            return tuple(x.T if x.flags.f_contiguous and x.ndim > 1 else x
                         for x in result)

        return result

    return func_wrapper


def check_args(refvaridx, refvarndim, rightdims, stagger=None,
               refstagdim=None):
    """A decorator to check that the wrapped function's arguments are valid.

    An exception is raised when an invalid argument is found.

    Args:

        refvaridx (:obj:`int`): The wrapped function's positional argument
            index to use as the reference variable.

        refvarndim (:obj:`int`): The number of dimensions for the reference
            variable that is expected by the wrapped function.

        rightdims (sequence of :obj:`int`): The expected number of right
            dimensions for each argument.

        stagger (sequence of :obj:`int` or :obj:`None`, optional): The
            dimension that is staggered for each argument in the wrapped
            function.  Use :obj:`None` in the sequence to indicate no
            staggering for that argument.  Default is None.

        refstagdim (:obj:`int`, optional): The staggered dimension for the
            reference variable, if applicable.  Default is None.

    Returns:

        None

    Raises:

        :class:`ValueError`: Raised when an invalid argument is detected.

    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):

        refvar = args[refvaridx]
        try:
            _ndim = refvar.ndim
        except AttributeError:
            raise ValueError("argument {} is not an arraylike "
                             "object".format(refvaridx))
        else:
            extra_dims = refvar.ndim - refvarndim

        # Always use unstaggered as the basis of comparison
        if refstagdim is not None:
            _refshape = list(refvar.shape)
            _refshape[refstagdim] -= 1
            _refshape = tuple(_refshape)
        else:
            _refshape = refvar.shape

        if stagger is None:
            _stagger = [None]*len(rightdims)
        else:
            _stagger = stagger

        for i, ndim in enumerate(rightdims):
            if ndim is None:
                continue

            var = args[i]

            try:
                _ = var.ndim
            except AttributeError:
                raise ValueError("argument {} is not an arraylike "
                                 "object".format(i))

            right_var_ndims = rightdims[i]

            # Check that the number of dims is correct
            if (var.ndim - extra_dims != right_var_ndims):
                raise ValueError("invalid number of dimensions for argument "
                                 "{} (got {}, expected {}).".format(
                                     i,
                                     var.ndim,
                                     right_var_ndims + extra_dims))

            # Add 1 to the reference staggered dim index before doing the check
            if _stagger[i] is not None:
                ref_shape = list(_refshape)
                ref_shape[_stagger[i]] += 1
                ref_shape = tuple(ref_shape)
            else:
                ref_shape = _refshape

            ref_right_sizes = ref_shape[extra_dims:]

            # Check that right dimensions are lined up
            if (var.shape[-right_var_ndims:] !=
                    ref_right_sizes[-right_var_ndims:]):

                raise ValueError("invalid shape for argument "
                                 "{} (got {}, expected {})".format(
                                     i,
                                     var.shape[-right_var_ndims:],
                                     ref_right_sizes[-right_var_ndims:]))

        return wrapped(*args, **kwargs)

    return func_wrapper
