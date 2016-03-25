#from functools import wraps
import wrapt 
from inspect import getargspec
from collections import OrderedDict

import numpy as np
import numpy.ma as ma

from .units import do_conversion, check_units
from .destag import destagger
from .util import (iter_left_indexes, viewitems, extract_vars,
                          combine_with, either)
from .config import xarray_enabled

if xarray_enabled():
    from xarray import DataArray



__all__ = ["convert_units", "handle_left_iter", "uvmet_left_iter", 
           "handle_casting" "copy_and_set_metadata", "set_wind_metadata",
           "set_latlon_metadata", "set_height_metadata"]

# TODO:  In python 3.x, getargspec is deprecated
def from_args(func, argnames, *args, **kwargs):
        """Parses the function args and kargs looking for the desired argument 
        value. Otherwise, the value is taken from the default keyword argument 
        using the arg spec.
        
        """
        if isinstance(argnames, str):
            arglist = [argnames]
        else:
            arglist = argnames
        
        result = {}
        for argname in arglist:
            arg_loc = arg_location(func, argname, args, kwargs)
            
            if arg_loc is not None:
                result[argname] = arg_loc[0][arg_loc[1]] 
            else:
                result[argname] = None
        
        return result

def arg_location(func, argname, args, kwargs):
    """Parses the function args, kargs and signature looking for the desired 
    argument location (either in args, kargs, or argspec.defaults), 
    and returns a tuple of argument_sequence, location.
    
    """

    argspec = getargspec(func)
    
    print argspec
    
    if argname not in argspec.args and argname not in kwargs:
        return None
        
    try:
        result_idx = argspec.args.index(argname)
    except ValueError:
        result_idx = None
    
    result = None
    if (argname in kwargs):
        result = kwargs, argname
    elif (len(args) > result_idx and result_idx is not None):
        result = args, result_idx
    else:    
        arg_idx_from_right = (len(argspec.args) - 
                              argspec.args.index(argname))
        result = list(argspec.defaults), -arg_idx_from_right
        
    return result
  
def convert_units(unit_type, alg_unit):
    """A decorator that applies unit conversion to a function's output array.
    
    Arguments:
    
        - unit_type - the unit type category (wind, pressure, etc)
        - alg_unit - the units that the function returns by default
    
    """
#    def convert_decorator(func):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        desired_units = from_args(wrapped, "units", *args, **kwargs)["units"]
        check_units(desired_units, unit_type)
        
        # Unit conversion done here
        return do_conversion(wrapped(*args, **kwargs), unit_type, 
                             alg_unit, desired_units)
    return func_wrapper
    
#    return convert_decorator


def _calc_out_dims(outvar, left_dims):
    left_dims = [x for x in left_dims]
    right_dims = [x for x in outvar.shape]
    return left_dims + right_dims 

def handle_left_iter(ref_var_expected_dims, ref_var_idx=-1,
                   ref_var_name=None,
                   ignore_args=None, ignore_kargs=None):
    """Decorator to handle iterating over leftmost dimensions when using 
    multiple files and/or multiple times.
    
    Arguments:
        - ref_var_expected_dims - the number of dimensions that the Fortran 
        algorithm is expecting for the reference variable
        - ref_var_idx - the index in args used as the reference variable for 
        calculating leftmost dimensions
        - ref_var_name - the keyword argument name for kargs used as the 
        reference varible for calculating leftmost dimensions
        - alg_out_fixed_dims - additional fixed dimension sizes for the 
        numerical algorithm (e.g. uvmet has a fixed left dimsize of 2)
        - ignore_args - indexes of any arguments which should be passed 
        directly without any slicing
        - ignore_kargs - keys of any keyword arguments which should be passed
        directly without slicing
        - ref_stag_dim - in some cases the reference variable dimensions
        have a staggered dimension that needs to be corrected when calculating
        the output dimensions (e.g. avo)
    
    """
#    def indexing_decorator(func):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        ignore_args = ignore_args if ignore_args is not None else ()
        ignore_kargs = ignore_kargs if ignore_kargs is not None else ()
        
        if ref_var_idx >= 0:
            ref_var = args[ref_var_idx]
        else:
            ref_var = kwargs[ref_var_name]
        
        ref_var_shape = ref_var.shape
        extra_dim_num = ref_var.ndim - ref_var_expected_dims
        
        # No special left side iteration, return the function result
        if (extra_dim_num == 0):
            return wrapped(*args, **kwargs)
        
        # Start by getting the left-most 'extra' dims
        extra_dims = [ref_var_shape[x] for x in xrange(extra_dim_num)]            
        
        out_inited = False
        for left_idxs in iter_left_indexes(extra_dims):
            # Make the left indexes plus a single slice object
            # The single slice will handle all the dimensions to
            # the right (e.g. [1,1,:])
            left_and_slice_idxs = ([x for x in left_idxs] + 
                                   [slice(None, None, None)])
            # Slice the args if applicable
            new_args = [arg[left_and_slice_idxs] 
                        if i not in ignore_args else arg 
                        for i,arg in enumerate(args)]
                
            # Slice the kargs if applicable
            new_kargs = {key:(val[left_and_slice_idxs] 
                         if key not in ignore_kargs else val)
                         for key,val in viewitems(kwargs)}
            
            # Call the numerical routine
            res = wrapped(*new_args, **new_kargs)
            
            if isinstance(res, np.ndarray):
                # Output array
                if not out_inited:
                    outdims = _calc_out_dims(res, extra_dims)
                    if not isinstance(res, ma.MaskedArray):
                        output = np.empty(outdims, ref_var.dtype)
                        masked = False
                    else:
                        output = ma.MaskedArray(
                                        np.zeros(outdims, ref_var.dtype),
                                        mask=np.zeros(outdims, np.bool_),
                                        fill_value=res.fill_value)
                        masked = True
                    
                    out_inited = True 
                
                if not masked:
                    output[left_and_slice_idxs] = res[:]
                else:
                    output.data[left_and_slice_idxs] = res.data[:]
                    output.mask[left_and_slice_idxs] = res.mask[:]
                
            else:   # This should be a list or a tuple (cape)
                if not out_inited:
                    outdims = _calc_out_dims(res[0], extra_dims)
                    if not isinstance(res[0], ma.MaskedArray):
                        output = [np.empty(outdims, ref_var.dtype) 
                                  for i in xrange(len(res))]
                        masked = False
                    else:
                        output = [ma.MaskedArray(
                                    np.zeros(outdims, ref_var.dtype),
                                    mask=np.zeros(outdims, np.bool_),
                                    fill_value=res[0].fill_value) 
                                  for i in xrange(len(res))]
                        masked = True
                    
                    out_inited = True
                
                for i,outarr in enumerate(res):
                    if not masked:
                        output[i][left_and_slice_idxs] = outarr[:]
                    else:
                        output[i].data[left_and_slice_idxs] = outarr.data[:]
                        output[i].mask[left_and_slice_idxs] = outarr.mask[:]
            
        
        return output
    
    return func_wrapper
    
#    return indexing_decorator

def uvmet_left_iter():
    """Decorator to handle iterating over leftmost dimensions when using 
    multiple files and/or multiple times with the uvmet product.
    
    """
    #def indexing_decorator(func):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        u = args[0]
        v = args[1]
        lat = args[2]
        lon = args[3]
        cen_long  = args[4]
        cone = args[5]
        
        if u.ndim == lat.ndim:
            num_right_dims = 2
            is_3d = False
        else:
            num_right_dims = 3
            is_3d = True
        
        is_stag = False
        if ((u.shape[-1] != lat.shape[-1]) or 
            (u.shape[-2] != lat.shape[-2])):
            is_stag = True
        
        if is_3d:
            extra_dim_num = u.ndim - 3
        else:
            extra_dim_num = u.ndim - 2
            
        if is_stag:
            u = destagger(u,-1)
            v = destagger(v,-2)
        
        # No special left side iteration, return the function result
        if (extra_dim_num == 0):
            return wrapped(u,v,lat,lon,cen_long,cone)
        
        # Start by getting the left-most 'extra' dims
        outdims = [u.shape[x] for x in xrange(extra_dim_num)]
        extra_dims = list(outdims) # Copy the left-most dims for iteration
        
        # Append the right-most dimensions
        outdims += [2] # For u/v components
        
        outdims += [u.shape[x] for x in xrange(-num_right_dims,0,1)]
        
        output = np.empty(outdims, u.dtype)
        
        for left_idxs in iter_left_indexes(extra_dims):
            # Make the left indexes plus a single slice object
            # The single slice will handle all the dimensions to
            # the right (e.g. [1,1,:])
            left_and_slice_idxs = ([x for x in left_idxs] + 
                                   [slice(None, None, None)])
                    
            
            new_u = u[left_and_slice_idxs]
            new_v = v[left_and_slice_idxs]
            new_lat = lat[left_and_slice_idxs]
            new_lon = lon[left_and_slice_idxs]
            
            # Call the numerical routine
            res = wrapped(new_u, new_v, new_lat, new_lon, cen_long, cone)
            
            # Note:  The 2D version will return a 3D array with a 1 length
            # dimension.  Numpy is unable to broadcast this without 
            # sqeezing first.
            res = np.squeeze(res) 
            
            output[left_and_slice_idxs] = res[:]
            
        
        return output
    
    return func_wrapper
    
#    return indexing_decorator

def handle_casting(ref_idx=0, arg_idxs=None, karg_names=None, dtype=np.float64):
    """Decorator to handle casting to/from required dtype used in 
    numerical routine.
    
    """
#    def cast_decorator(func):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        arg_idxs = arg_idxs if arg_idxs is not None else ()
        karg_names = karg_names if karg_names is not None else ()
        
        orig_type = args[ref_idx].dtype
        
        new_args = [arg.astype(dtype)
                    if i in arg_idxs else arg 
                    for i,arg in enumerate(args)]
        
        new_kargs = {key:(val.astype(dtype)
                    if key in karg_names else val)
                    for key,val in viewitems()}
        
        res = wrapped(*new_args, **new_kargs) 
        
        if isinstance(res, np.ndarray):
            if res.dtype == orig_type:
                return res
            return res.astype(orig_type)
        else:   # got back a sequence of arrays
            return tuple(arr.astype(orig_type) 
                         if arr.dtype != orig_type else arr
                         for arr in res)
    
    return func_wrapper
    
#    return cast_decorator

def _extract_and_transpose(arg):
    if xarray_enabled():
        if isinstance(arg, DataArray):
            arg = arg.values
    
    if isinstance(arg, np.ndarray):
        if not arg.flags.F_CONTIGUOUS:
            return arg.T
        
    return arg
    
def handle_extract_transpose():
    """Decorator to extract the data array from a DataArray and also
    transposes if the data is not fortran contiguous.
    
    """
#    def extract_transpose_decorator(func):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        new_args = [_extract_and_transpose(arg) for arg in args]
        
        new_kargs = {key:_extract_and_transpose(val) 
                     for key,val in viewitems(kwargs)}
            
        res = wrapped(*new_args, **new_kargs) 
        
        if isinstance(res, np.ndarray):
            if res.flags.F_CONTIGUOUS:
                return res.T
        else:
            return tuple(x.T if x.flags.F_CONTIGUOUS else x for x in res)
        
        return res
    
    return func_wrapper
    
#    return extract_transpose_decorator


def copy_and_set_metadata(copy_varname=None, delete_attrs=None, name=None,
                          dimnames=None, coords=None, **fixed_attrs):
    """Decorator to set the metadata for a WRF method.
    
    A cache is inserted/updated to include the extracted variable that will 
    have its metadata copied. This prevents the variable being extracted more 
    than once.  This extraction can be slow with sequences of large files.
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        if not xarray_enabled():
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfnc", "timeidx", "method", 
                                      "squeeze", "cache", "units"), 
                            *args, **kwargs)
        wrfnc = argvars["wrfnc"]
        timeidx = argvars["timeidx"]
        units = argvars["units"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        cache = argvars["cache"]
        if cache is None:
            cache = {}
        
        
        if (callable(copy_varname)):
            copy_varname = copy_varname(wrfnc)
        
        # Extract the copy_from argument
        var_to_copy = None if cache is None else cache.get(copy_varname, 
                                                           None)

            
        if var_to_copy is None:
            var_to_copy = extract_vars(wrfnc, timeidx, (copy_varname,), 
                                       method, squeeze, cache)[copy_varname]
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[copy_varname] = var_to_copy
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args = list(args)
        new_kargs = dict(kwargs)
        cache_argseq, cache_argloc = arg_location(wrapped, "cache", 
                                                  new_args, new_kargs)
        
        cache_argseq[cache_argloc] = new_cache
        
        result = wrapped(*new_args, **new_kargs)
        
        outname = ""
        outdimnames = list()
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        
        if copy_varname is not None:
            outname = str(var_to_copy.name)
            
            if dimnames is not None:
                if isinstance(combine_with):
                    outdimnames, outcoords = dimnames(copy_varname)
                else:
                    outdimnames = dimnames
                    outcoords = coords
            else:
                outdimnames += var_to_copy.dims
                outcoords.update(var_to_copy.coords)
            
            outattrs.update(var_to_copy.attrs)   
        
        if name is not None:
            outname = name
        
        if units is not None:
            outattrs["units"] = units
            
        for argname, val in viewitems(fixed_attrs):
            outattrs[argname] = val
        
        if delete_attrs is not None:
            for attr in delete_attrs:
                del outattrs[attr]
                
        if isinstance(ma.MaskedArray):
            outattrs["_FillValue"] = result.fill_value
            outattrs["missing_value"] = result.fill_value
        
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
    
    return func_wrapper


def set_wind_metadata(wspd_wdir=False):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        if not xarray_enabled():
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfnc", "timeidx", "units", 
                                      "method", "squeeze", "ten_m", "cache"), 
                          *args, **kwargs)
        wrfnc = argvars["wrfnc"]
        timeidx = argvars["timeidx"]
        units = argvars["units"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        ten_m = argvars["ten_m"]
        cache = argvars["cache"]
        if cache is None:
            cache = {}
        
        lat_var = either("XLAT", "XLAT_M")(wrfnc)
        xlat_var = extract_vars(wrfnc, timeidx, lat_var, 
                                method, squeeze, cache)
        lat = xlat_var[lat_var]
        
        lon_var = either("XLONG", "XLONG_M")
        xlon_var = extract_vars(wrfnc, timeidx, lon_var, 
                                method, squeeze, cache)
        lon = xlon_var[lon_var]
        
        pres_var = either("P", "PRES")
        p_var = extract_vars(wrfnc, timeidx, lon_var, method, squeeze, cache)
        pres = p_var[pres_var]
        
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[lat_var] = lat
        new_cache[lon_var] = lon
        new_cache[pres_var] = pres
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args = list(args)
        new_kargs = dict(kwargs)
        cache_argseq, cache_argloc = arg_location(wrapped, "cache", 
                                                  new_args, new_kargs)
        
        cache_argseq[cache_argloc] = new_cache
        
        result = wrapped(*new_args, **new_kargs)
        
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outname = "uvmet" if not ten_m else "uvmet10"
        
        outdimnames = list(pres.dims)
        outcoords.update(pres.coords)
        outattrs.update(pres.attrs)
        
        if not wspd_wdir:
            outattrs["description"] = ("earth rotated u,v" if not ten_m else 
                                   "10m earth rotated u,v")
            if not ten_m:
                outdimnames.insert(-3, "u_v")
            else:
                outdimnames.insert(-2, "u_v")
        else:
            outattrs["description"] = ("earth rotated wspd,wdir" if not ten_m 
                                       else "10m earth rotated wspd,wdir")
            if not ten_m:
                outdimnames.insert(-3, "wspd_wdir")
            else:
                outdimnames.insert(-2, "wspd_wdir")
            
        if units is not None: 
            outattrs["units"] = units
            
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
        
    return func_wrapper

def set_latlon_metadata(ij=False):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        if not xarray_enabled():
            return wrapped(*args, **kwargs)
        
        res = wrapped(*args, **kwargs)
        
        argnames = ("latitude", "longitude") if not ij else ("i", "j")
        outname = "latlon" if not ij else "ij"
        
        if res.ndim <= 2:
            dimnames = (["lat_lon", "i_j"] if not ij 
                        else ["i_j", "lat_lon"])
        else:
            dimnames = (["lat_lon", "domain", "i_j"] if not ij 
                        else ["i_j", "domain", "lat_lon"])
            
        
        argvars = from_args(wrapped, argnames, *args, **kwargs)
        
        var1 = argvars[argnames[0]]
        var2 = argvars[argnames[1]]
        
        arr1 = np.asarray(var1).ravel()
        arr2 = np.asarray(var2).ravel()
        
        coords = {}
        coords[dimnames[0]] = [x for x in zip(arr1, arr2)]
        
        return DataArray(res, name=outname, dims=dimnames, coords=coords)
    
    return func_wrapper
    
def set_height_metadata(geopt=False):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        if not xarray_enabled():
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfnc", "timeidx", "method", 
                                    "squeeze", "units", "msl", "cache"), 
                          *args, **kwargs)
        wrfnc = argvars["wrfnc"]
        timeidx = argvars["timeidx"]
        units = argvars["units"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        msl = argvars["msl"]
        cache = argvars["cache"]
        if cache is None:
            cache = {}
        
        # For height, either copy the met_em GHT variable or copy and modify
        # pressure (which has the same dims as destaggered height)
        ht_metadata_varname = either("P", "GHT")(wrfnc)
        ht_var = extract_vars(wrfnc, timeidx, ht_metadata_varname, 
                              method, squeeze, cache)
        ht_metadata_var = ht_var[ht_metadata_varname]
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[ht_metadata_var] = ht_metadata_var
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args = list(args)
        new_kargs = dict(kwargs)
        cache_argseq, cache_argloc = arg_location(wrapped, "cache", 
                                                  new_args, new_kargs)
        
        cache_argseq[cache_argloc] = new_cache
        
        result = wrapped(*new_args, **new_kargs)
        
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(ht_metadata_var.dims)
        outcoords.update(ht_metadata_var.coords)
        outattrs.update(ht_metadata_var.attrs)
        
        if geopt:
            outname = "geopt"
            outattrs["units"] = "m2 s-2"
            outattrs["description"] = "full model geopotential"
        else:
            outname = "height" if msl else "height_agl" 
            outattrs["units"] = units
            height_type = "MSL" if msl else "AGL"
            outattrs["description"] = "model height ({})".format(height_type)
        
        
        return DataArray(result, name=outname, 
                         dims=outdimnames, coords=outcoords)
        
        
    
    


