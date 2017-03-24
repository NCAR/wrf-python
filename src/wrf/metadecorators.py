from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
import wrapt 
from collections import OrderedDict

import numpy as np
import numpy.ma as ma

from .extension import _interpline
from .util import (extract_vars, either, from_args, arg_location,
                   is_coordvar, latlon_coordvars, to_np, 
                   from_var, iter_left_indexes)
from .coordpair import CoordPair
from .py3compat import viewkeys, viewitems, py3range, ucode
from .interputils import get_xy_z_params, get_xy, to_xy_coords
from .config import xarray_enabled

if xarray_enabled():
    from xarray import DataArray
    

def copy_and_set_metadata(copy_varname=None, delete_attrs=None, name=None,
                          remove_dims=None, dimnames=None, 
                          coords=None, **fixed_attrs):
    """A decorator that sets the metadata for a wrapped function's output.
    
    Generally, the metadata is copied from the variable specified by 
    *copy_varname*, with other fixed fields set by the other decorator 
    arguments.
    
    The *cache* argument used by most diagnostic routines is supplied by 
    this decorator when the *copy_varname* variable is extracted, in order
    to prevent the variable from being extracted again by the wrapped function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        copy_varname (:obj:`str`, optional): The NetCDF variable name to copy.
            Default is None.
        
        delete_attrs (sequence of :obj:`str`, optional): A sequence of key 
            names to remove from the :attr:`xarray.DataArray.attrs` attribute 
            in the wrapped function output (after being copied from 
            *copy_varname*).  Default is None.
            
        name (:obj:`str`): The name to use for the 
            :attr:`xarray.DataArray.name` attribute.
            
        remove_dims (sequence of :obj:`int`, optional): A sequence of dimension 
            indexes to be removed from the wrapped function output (after being 
            copied from *copy_varname*).  This is useful when the copy 
            variable is three dimensional but the wrapped function output 
            is two dimensional, and you still want to keep the names of 
            the rightmost two dimensions.  Default is None.
        
        dimnames (sequence of :obj:`str`, optional):  A sequence of dimension
            names in order to manually set the 
            :attr:`xarray.DataArray.dims` attribute.  Default is None.
            
        coords (:obj:`dict`): A mapping of coordinate name to coordinate 
            value to manually specify the :attr:`xarray.DataArray.coords` 
            attribute.  Default is None.
            
        **fixed_attrs: These keyword arguments are added to the 
            :attr:`xarray.DataArray.attrs` attribute.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        function output with or without metadata.  If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs): 
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfin", "timeidx", "method", 
                                      "squeeze", "cache", "units", "meta",
                                      "_key"), 
                            *args, **kwargs)
        
        wrfin = argvars["wrfin"]
        timeidx = argvars["timeidx"]
        units = argvars["units"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        cache = argvars["cache"]
        _key = argvars["_key"]
        if cache is None:
            cache = {}
        
        # Note:  can't modify nonlocal var
        if (callable(copy_varname)):
            _copy_varname = copy_varname(wrfin)
        else:
            _copy_varname = copy_varname
        
        # Extract the copy_from argument
        var_to_copy = None if cache is None else cache.get(_copy_varname, 
                                                           None)

            
        if var_to_copy is None:
            var_to_copy = extract_vars(wrfin, timeidx, (_copy_varname,), 
                                       method, squeeze, cache,
                                       meta=True, _key=_key)[_copy_varname]
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[_copy_varname] = var_to_copy
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args, cache_argloc = arg_location(wrapped, "cache", args, kwargs)
        new_args[cache_argloc] = new_cache
        
        result = wrapped(*new_args)
        
        outname = ""
        outdimnames = list()
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        
        if copy_varname is not None:
            outname = ucode(var_to_copy.name)
            
            if dimnames is not None:
                outdimnames = dimnames
                outcoords = coords
            else:
                outdimnames += var_to_copy.dims
                outcoords.update(var_to_copy.coords)
            
            outattrs.update(var_to_copy.attrs)
            
            if remove_dims is not None:
                for dimname in remove_dims:
                    outdimnames.remove(dimname)
                    
                    try:
                        del outcoords[dimname]
                    except KeyError:
                        pass
                     
        
        if name is not None:
            outname = name
        
        if units is not None:
            outattrs["units"] = units
            
        for argname, val in viewitems(fixed_attrs):
            outattrs[argname] = val
        
        if delete_attrs is not None:
            for attr in delete_attrs:
                try:
                    del outattrs[attr]
                except KeyError:
                    pass
                
        if isinstance(result, ma.MaskedArray):
            outattrs["_FillValue"] = result.fill_value
            outattrs["missing_value"] = result.fill_value
        
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
    
    return func_wrapper


def set_wind_metadata(copy_varname, name, description, 
                      wind_ncvar=False, 
                      two_d=False, wspd_wdir=False):
    """A decorator that sets the metadata for a wrapped wind function's output.
    
    This is a special metadata decorator for working with wind functions, which
    include wind extraction routines, uvmet routines, and 
    wind speed / wind direction routines.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        copy_varname (:obj:`str`, optional): The NetCDF variable name to copy.
            Default is None.
            
        name (:obj:`str`): The name to use for the 
            :attr:`xarray.DataArray.name` attribute.
            
        description (:obj:`str`): The description to use for the 'description' 
            key in the :attr:`xarray.DataArray.attrs` attribute.
        
        wind_ncvar (:obj:`bool`, optional): Set to True when the wrapped 
            function is simply extracting a wind variable (U, V, W) from the 
            NetCDF file.  Set to False for other types of wind algorithms 
            (uvmet, wspd_wdir, etc). Default is False.
            
        two_d (:obj:`bool`, optional): Set to True if the wind field is 
            two-dimensional. Set to False for a three-dimensional wind field.
            Default is False.
            
        wspd_wdir (:obj:`bool`): Set to True if the wrapped function is a 
            wind speed / wind direction algorithm.  Otherwise, set to False.
            Default is False.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        wind function output with or without metadata.  If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfin", "timeidx", "units", 
                                      "method", "squeeze", "ten_m", "cache",
                                      "_key"), 
                          *args, **kwargs)
        wrfin = argvars["wrfin"]
        timeidx = argvars["timeidx"]
        units = argvars["units"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        ten_m = argvars["ten_m"]
        cache = argvars["cache"]
        _key = argvars["_key"]
        if cache is None:
            cache = {}
        
        if isinstance(copy_varname, either):
            _copy_varname = copy_varname(wrfin)
        else:
            _copy_varname = copy_varname
        
        copy_var = extract_vars(wrfin, timeidx, _copy_varname, 
                                method, squeeze, cache, 
                                meta=True, _key=_key)[_copy_varname]
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[_copy_varname] = copy_var
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args, cache_argloc = arg_location(wrapped, "cache", args, kwargs)
        new_args[cache_argloc] = new_cache
        
        result = wrapped(*new_args)
        
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        
        outdimnames = list(copy_var.dims)
        #outcoords.update(copy_var.coords)
        outattrs.update(copy_var.attrs)
        
        
        if wind_ncvar:
            outcoords.update(copy_var.coords)
        elif not wspd_wdir:
            if not two_d:
                outdimnames.insert(0, "u_v")
            else:
                outdimnames.insert(0, "u_v")
                outattrs["MemoryOrder"] = "XY"
            outcoords["u_v"] = ["u", "v"]
        else:
            if not two_d:
                outdimnames.insert(0, "wspd_wdir")
            else:
                outdimnames.insert(0, "wspd_wdir")
                outattrs["MemoryOrder"] = "XY"
                
            outcoords["wspd_wdir"] = ["wspd", "wdir"]
        
        if units is not None: 
            outattrs["units"] = units
            
        # xarray doesn't line up coordinate dimensions based on 
        # names, it just remembers the index it originally mapped to.  
        # So, need to rebuild the XLAT, XLONG, coordinates again since the 
        # leftmost index changed.
        if not wind_ncvar:
            for key,dataarray in viewitems(copy_var.coords):
                if is_coordvar(key):
                    outcoords[key] = dataarray.dims, to_np(dataarray)
                elif key == "XTIME":
                    outcoords[key] = dataarray.dims, to_np(dataarray)
                elif key == "Time":
                    outcoords[key] = to_np(dataarray)
                   
        outname = name
        outattrs["description"] = description
        
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
        
    return func_wrapper

def set_cape_metadata(is2d):
    """A decorator that sets the metadata for a wrapped CAPE function's output.
    
    This is a special metadata decorator for working with CAPE functions.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        is2d (:obj:`bool`): Set to True if the wrapped function is for a 
            two-dimensional CAPE routine.  Set to False for a
            three-dimensional CAPE routine.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        CAPE function output with or without metadata.  If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfin", "timeidx", "method", "squeeze", 
                                      "cache", "_key", "missing"), 
                          *args, **kwargs)
        wrfin = argvars["wrfin"]
        timeidx = argvars["timeidx"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        cache = argvars["cache"]
        missing = argvars["missing"]
        _key = argvars["_key"]
        if cache is None:
            cache = {}
        
        _copy_varname = "P"
        copy_var = extract_vars(wrfin, timeidx, _copy_varname, method, squeeze, 
                                cache, meta=True, _key=_key)[_copy_varname]
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[_copy_varname] = copy_var
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args, cache_argloc = arg_location(wrapped, "cache", args, kwargs)
        new_args[cache_argloc] = new_cache
        
        result = wrapped(*new_args)
        
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outattrs.update(copy_var.attrs)
        outdimnames = [None] * result.ndim
        
        if is2d:
            # Right dims
            outdimnames[-2:] = copy_var.dims[-2:]
            # Left dims
            outdimnames[1:-2] = copy_var.dims[0:-3]
            outdimnames[0] = "mcape_mcin_lcl_lfc"
            outattrs["description"] = "mcape ; mcin ; lcl ; lfc"
            outattrs["MemoryOrder"] = "XY"
            outattrs["units"] = "J kg-1 ; J kg-1 ; m ; m"
            outname = "cape_2d"
        else:
            # Right dims
            outdimnames[-3:] = copy_var.dims[-3:]
            # Left dims
            outdimnames[1:-3] = copy_var.dims[0:-3]
            outdimnames[0] = "cape_cin"
            outattrs["description"] = "cape; cin"
            outattrs["units"] = "J kg-1 ; J kg-1"
            outattrs["MemoryOrder"] = "XYZ"
            outname = "cape_3d"
        
        outattrs["_FillValue"] = missing
        outattrs["missing_value"] = missing

        
        # xarray doesn't line up coordinate dimensions based on 
        # names, it just remembers the index it originally mapped to.  
        # So, need to rebuild the XLAT, XLONG, coordinates again since the 
        # leftmost index changed.

        for key,dataarray in viewitems(copy_var.coords):
            if is_coordvar(key):
                outcoords[key] = dataarray.dims, to_np(dataarray)
            elif key == "XTIME":
                outcoords[key] = dataarray.dims, to_np(dataarray)
            elif key == "Time":
                outcoords[key] = to_np(dataarray)
        
        if is2d:
            outcoords["mcape_mcin_lcl_lfc"] = ["mcape", "mcin", "lcl", "lfc"]
        else:
            outcoords["cape_cin"] = ["cape", "cin"]
            
        
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
        
    return func_wrapper


def set_cloudfrac_metadata():
    """A decorator that sets the metadata for a wrapped cloud fraction 
    function's output.
    
    This is a special metadata decorator for working with cloud fraction 
    functions.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        cloud fraction function output with or without metadata.  If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfin", "timeidx", "method", "squeeze", 
                                      "cache", "_key"), 
                          *args, **kwargs)
        wrfin = argvars["wrfin"]
        timeidx = argvars["timeidx"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        cache = argvars["cache"]
        _key = argvars["_key"]
        if cache is None:
            cache = {}
        
        _copy_varname = "P"
        copy_var = extract_vars(wrfin, timeidx, _copy_varname, method, squeeze, 
                                cache, meta=True, _key=_key)[_copy_varname]
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[_copy_varname] = copy_var
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args, cache_argloc = arg_location(wrapped, "cache", args, kwargs)
        new_args[cache_argloc] = new_cache
        
        result = wrapped(*new_args)
        
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outattrs.update(copy_var.attrs)
        outdimnames = [None] * result.ndim
        
        # Right dims
        outdimnames[-2:] = copy_var.dims[-2:]
        # Left dims
        outdimnames[1:-2] = copy_var.dims[0:-3]
        outdimnames[0] = "low_mid_high"
        outattrs["description"] = "low, mid, high clouds"
        outattrs["MemoryOrder"] = "XY"
        outattrs["units"] = "%"
        outname = "cloudfrac"

        # xarray doesn't line up coordinate dimensions based on 
        # names, it just remembers the index it originally mapped to.  
        # So, need to rebuild the XLAT, XLONG, coordinates again since the 
        # leftmost index changed.

        for key,dataarray in viewitems(copy_var.coords):
            if is_coordvar(key):
                outcoords[key] = dataarray.dims, to_np(dataarray)
            elif key == "XTIME":
                outcoords[key] = dataarray.dims, to_np(dataarray)
            elif key == "Time":
                outcoords[key] = to_np(dataarray)
        
        outcoords["low_mid_high"] = ["low", "mid", "high"]
            
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
        
    return func_wrapper

def set_latlon_metadata(xy=False):
    """A decorator that sets the metadata for a wrapped latlon function's 
    output.
    
    This is a special metadata decorator for working with latlon functions.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        xy (:obj:`bool`, optional): Set to True if the wrapped function returns 
            xy values (ll_to_xy).  Set to False if the wrapped function returns 
            latlon values (xy_to_ll).  Default is False.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        latlon function output with or without metadata.  If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        # Set squeeze to False.  Squeeze will be handled in here
        new_args, squeeze_argloc = arg_location(wrapped, "squeeze", args, 
                                                kwargs)
        new_args[squeeze_argloc] = False
        
        result = wrapped(*new_args)
        
        # Want to preserve the input coordinate pair in metadata
        if result.ndim == 1:
            result = result[:, np.newaxis]
        
        argnames = ["x", "y"] if not xy else ["latitude", "longitude"]
        argnames.append("squeeze")
        outname = "latlon" if not xy else "xy"
        
        if result.ndim == 2:
            dimnames = (["lat_lon", "idx"] if not xy else ["x_y", "idx"])
        else:
            dimnames = (["lat_lon", "domain_idx", "idx"] if not xy 
                        else ["x_y", "domain_idx", "idx"])
        
        argvars = from_args(wrapped, argnames, *args, **kwargs)
        
        var1 = argvars[argnames[0]]
        var2 = argvars[argnames[1]]
        squeeze = argvars["squeeze"]
        
        arr1 = np.asarray(var1).ravel()
        arr2 = np.asarray(var2).ravel()
        
        coords = {}
        if not xy:
            coords["xy_coord"] = (dimnames[-1], [CoordPair(x=x[0], y=x[1]) 
                               for x in zip(arr1, arr2)])
            coords[dimnames[0]] = ["lat", "lon"]
        else:
            coords["latlon_coord"] = (dimnames[-1], [CoordPair(lat=x[0], 
                                                               lon=x[1]) 
                               for x in zip(arr1, arr2)])
            coords[dimnames[0]] = ["x", "y"]
        
        da = DataArray(result, name=outname, dims=dimnames, coords=coords)
        
        if squeeze:
            da = da.squeeze()
        
        return da
    
    return func_wrapper
    
def set_height_metadata(geopt=False):
    """A decorator that sets the metadata for a wrapped height function's 
    output.
    
    This is a special metadata decorator for working with height functions.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        geopt (:obj:`bool`, optional): Set to True if the wrapped function 
            returns geopotential.  Set to True if the wrapped function 
            returns geopotential height.  Default is False.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        height function output with or without metadata.  If xarray is 
        enabled and the *meta* parameter is True, then the result will be a 
        :class:`xarray.DataArray` object.  Otherwise, the result will be a 
        :class:`numpy.ndarray` object with no metadata.
        
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfin", "timeidx", "method", 
                                    "squeeze", "units", "msl", "cache",
                                    "_key"), 
                          *args, **kwargs)
        wrfin = argvars["wrfin"]
        timeidx = argvars["timeidx"]
        units = argvars["units"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        msl = argvars["msl"]
        cache = argvars["cache"]
        _key = argvars["_key"]
        
        if cache is None:
            cache = {}
        
        # For height, either copy the met_em GHT variable or copy and modify
        # pressure (which has the same dims as destaggered height)
        ht_metadata_varname = either("P", "GHT")(wrfin)
        ht_var = extract_vars(wrfin, timeidx, ht_metadata_varname, 
                              method, squeeze, cache, meta=True,
                              _key=_key)
        ht_metadata_var = ht_var[ht_metadata_varname]
        
        # Make a copy so we don't modify a user supplied cache
        new_cache = dict(cache) 
        new_cache[ht_metadata_varname] = ht_metadata_var
        
        # Don't modify the original args/kargs.  The args need to be a list
        # so it can be modified.
        new_args, cache_argloc = arg_location(wrapped, "cache", args, kwargs)
        new_args[cache_argloc] = new_cache
        
        result = wrapped(*new_args)
        
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
                         dims=outdimnames, coords=outcoords, attrs=outattrs)
    return func_wrapper

def _set_horiz_meta(wrapped, instance, args, kwargs):  
    """A decorator implementation that sets the metadata for a wrapped 
    horizontal interpolation function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        wrapped: The wrapped function which in turns needs to be called by your 
            wrapper function.
        
        instance: The object to which the wrapped function was bound when it 
            was called.
        
        args: The list of positional arguments supplied when the decorated 
            function was called.
        
        kwargs: The dictionary of keyword arguments supplied when the decorated 
            function was called.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        horiztontal interpolation function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    See Also:
    
        :mod:`wrapt`
        
    """  
    argvars = from_args(wrapped, ("field3d", "vert", "desiredlev", 
                                  "missing"), 
                          *args, **kwargs)  
    
    field3d = argvars["field3d"]
    z = argvars["vert"]
    desiredloc = argvars["desiredlev"]
    missingval = argvars["missing"]
    
    result = wrapped(*args, **kwargs)
    
    # Defaults, in case the data isn't a DataArray
    outname = None
    outdimnames = None
    outcoords = None
    outattrs = None
    
    # Get the vertical level units
    vert_units = None
    if isinstance(z, DataArray):
        vert_units = z.attrs.get("units", None)
    
    # If we have no metadata to start with, only set the level
    levelstr = ("{0} {1}".format(desiredloc, vert_units) 
                if vert_units is not None 
                else "{0}".format(desiredloc))
    
    name_levelstr = ("{0}_{1}".format(desiredloc, vert_units) 
                if vert_units is not None 
                else "{0}".format(desiredloc))
    
    if isinstance(field3d, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(field3d.dims)
        outcoords.update(field3d.coords)
        outdimnames.remove(field3d.dims[-3])
        try:
            del outcoords[field3d.dims[-3]]
        except KeyError:
            pass # xarray 0.9
        outattrs.update(field3d.attrs)
        outname = "{0}_{1}".format(field3d.name, name_levelstr)
        
    else:
        outname = "field3d_{0}".format(name_levelstr)
        outattrs = OrderedDict()
        
    outattrs["level"] = levelstr
    outattrs["missing_value"] = missingval
    outattrs["_FillValue"] = missingval
    
    for key in ("MemoryOrder", "description"):
        try:
            del outattrs[key]
        except KeyError:
            pass
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs)
    
def _set_cross_meta(wrapped, instance, args, kwargs):
    """A decorator implementation that sets the metadata for a wrapped cross \
    section interpolation function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        wrapped: The wrapped function which in turns needs to be called by your 
            wrapper function.
        
        instance: The object to which the wrapped function was bound when it 
            was called.
        
        args: The list of positional arguments supplied when the decorated 
            function was called.
        
        kwargs: The dictionary of keyword arguments supplied when the decorated 
            function was called.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        cross section interpolation function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    See Also:
    
        :mod:`wrapt`
        
    """  
    argvars = from_args(wrapped, ("field3d", "vert", "levels", 
                                  "latlon", "missing", 
                                  "wrfin", "timeidx", "stagger", "projection",
                                  "pivot_point", "angle",
                                  "start_point", "end_point",
                                  "cache"), 
                          *args, **kwargs)  
    
    field3d = argvars["field3d"]
    z = argvars["vert"]
    levels = argvars["levels"]
    inc_latlon = argvars["latlon"]
    missingval = argvars["missing"]
    wrfin = argvars["wrfin"]
    timeidx = argvars["timeidx"]
    stagger = argvars["stagger"]
    projection = argvars["projection"]
    pivot_point = argvars["pivot_point"]
    angle = argvars["angle"]
    start_point = argvars["start_point"]
    end_point = argvars["end_point"]
    cache = argvars["cache"]
    
    start_point_xy = None
    end_point_xy = None
    pivot_point_xy = None
    
    if pivot_point is not None:
        if pivot_point.lat is not None and pivot_point.lon is not None:
            xy_coords = to_xy_coords(pivot_point, wrfin, timeidx, 
                                     stagger, projection)
            pivot_point_xy = (xy_coords.x, xy_coords.y)
        else:
            pivot_point_xy = (pivot_point.x, pivot_point.y)
            
    if start_point is not None and end_point is not None:
        if start_point.lat is not None and start_point.lon is not None:
            xy_coords = to_xy_coords(start_point, wrfin, timeidx, 
                                     stagger, projection)
            start_point_xy = (xy_coords.x, xy_coords.y)
        else:
            start_point_xy = (start_point.x, start_point.y)
            
        if end_point.lat is not None and end_point.lon is not None:
            xy_coords = to_xy_coords(end_point, wrfin, timeidx, 
                                     stagger, projection)
            end_point_xy = (xy_coords.x, xy_coords.y)
        else:
            end_point_xy = (end_point.x, end_point.y)
    
    xy, var2dz, z_var2d = get_xy_z_params(to_np(z), pivot_point_xy, angle,
              start_point_xy, end_point_xy, levels)
    
    # Make a copy so we don't modify a user supplied cache
    if cache is not None:
        new_cache = dict(cache)
    else:
        new_cache = {}
    new_cache["xy"] = xy
    new_cache["var2dz"] = var2dz
    new_cache["z_var2d"] = z_var2d
    
    # Don't modify the original args/kargs.  The args need to be a list
    # so it can be modified.
    new_args, cache_argloc = arg_location(wrapped, "cache", args, kwargs)
    new_args[cache_argloc] = new_cache
        
    result = wrapped(*new_args)
    
    # Defaults, in case the data isn't a DataArray
    outname = None
    outdimnames = None
    outcoords = None
    outattrs = None
    
    # Use XY to set the cross-section metadata
    st_x = xy[0,0]
    st_y = xy[0,1]
    ed_x = xy[-1,0]
    ed_y = xy[-1,1]
    
    cross_str = "({0}, {1}) to ({2}, {3})".format(st_x, st_y, ed_x, ed_y)
    if angle is not None:
        cross_str += " ; center={0} ; angle={1}".format(pivot_point,
                                                        angle)
    
    if isinstance(field3d, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(field3d.dims)
        outcoords.update(field3d.coords)
        for i in py3range(-3,0,1):
            outdimnames.remove(field3d.dims[i])
            try:
                del outcoords[field3d.dims[i]]
            except KeyError:
                pass # Xarray 0.9
        
        
        # Delete any lat,lon coords    
        delkeys = [key for key in viewkeys(outcoords) if is_coordvar(key)]
        for key in delkeys:
            del outcoords[key]
        
        outdimnames.append("vertical")
        outdimnames.append("cross_line_idx")
        outattrs.update(field3d.attrs)
        
        outname = "{0}_cross".format(field3d.name)
        
        for key in ("MemoryOrder",):
            try:
                del outattrs[key]
            except KeyError:
                pass
        
        # Interpolate to get the lat/lon coords, if desired
        if inc_latlon:
            latcoordname, loncoordname = latlon_coordvars(field3d.coords)
            
            if latcoordname is not None and loncoordname is not None:
                latcoord = field3d.coords[latcoordname]
                loncoord = field3d.coords[loncoordname]
                
                if latcoord.ndim == 2:
                    lats = _interpline(latcoord, xy)
                    lons = _interpline(loncoord, xy)
                    
                    outcoords["xy_loc"] = ("cross_line_idx", 
                                           np.asarray(tuple(
                                                CoordPair(x=xy[i,0], y=xy[i,1],
                                                    lat=lats[i], lon=lons[i]) 
                                          for i in py3range(xy.shape[-2])))
                                          )
                # Moving domain
                else:
                    extra_dims = latcoord.shape[0:-2]
                    outdims = extra_dims + xy.shape[-2:-1]
                    
                    latlon_loc = np.empty(outdims, np.object_)
                    for left_dims in iter_left_indexes(extra_dims):
                        idxs = left_dims + (slice(None),)
                        lats = _interpline(latcoord[idxs], xy)
                        lons = _interpline(loncoord[idxs], xy)
                        
                        latlon_loc[idxs] = np.asarray(tuple(
                                            CoordPair(x=xy[i,0], y=xy[i,1],
                                                    lat=lats[i], lon=lons[i]) 
                                            for i in py3range(xy.shape[-2]))
                                            )[:] 
                        
                    
                    extra_dimnames = latcoord.dims[0:-2]
                    loc_dimnames = extra_dimnames + ("cross_line_idx",)
                    outcoords["xy_loc"] = (loc_dimnames, latlon_loc)
                    
            else:
                outcoords["xy_loc"] = ("cross_line_idx", np.asarray(tuple(
                                                CoordPair(xy[i,0], xy[i,1]) 
                                          for i in py3range(xy.shape[-2]))))
            
        else:    
            outcoords["xy_loc"] = ("cross_line_idx", np.asarray(tuple(
                                                CoordPair(xy[i,0], xy[i,1]) 
                                          for i in py3range(xy.shape[-2]))))
        
        outcoords["vertical"] = z_var2d[:]
        
    else:
        outname = "field3d_cross"
        outattrs = OrderedDict()
    
    outattrs["orientation"] = cross_str
    outattrs["missing_value"] = missingval
    outattrs["_FillValue"] = missingval
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs)  
    
    

def _set_line_meta(wrapped, instance, args, kwargs):
    """A decorator implementation that sets the metadata for a wrapped line 
    interpolation function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        wrapped: The wrapped function which in turns needs to be called by your 
            wrapper function.
        
        instance: The object to which the wrapped function was bound when it 
            was called.
        
        args: The list of positional arguments supplied when the decorated 
            function was called.
        
        kwargs: The dictionary of keyword arguments supplied when the decorated 
            function was called.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        line interpolation function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    See Also:
    
        :mod:`wrapt`
        
    """   
    argvars = from_args(wrapped, ("field2d", 
                                  "wrfin", "timeidx", "stagger", "projection",
                                  "pivot_point", "angle",
                                  "start_point", "end_point", "latlon", 
                                  "cache"), 
                          *args, **kwargs)  
    
    field2d = argvars["field2d"]
    wrfin = argvars["wrfin"]
    timeidx = argvars["timeidx"]
    stagger = argvars["stagger"]
    projection = argvars["projection"]
    pivot_point = argvars["pivot_point"]
    angle = argvars["angle"]
    start_point = argvars["start_point"]
    end_point = argvars["end_point"]
    inc_latlon = argvars["latlon"]
    cache = argvars["cache"]
    
    if cache is None:
        cache = {}
        
    start_point_xy = None
    end_point_xy = None
    pivot_point_xy = None
        
    if pivot_point is not None:
        if pivot_point.lat is not None and pivot_point.lon is not None:
            xy_coords = to_xy_coords(pivot_point, wrfin, timeidx, 
                                     stagger, projection)
            pivot_point_xy = (xy_coords.x, xy_coords.y)
        else:
            pivot_point_xy = (pivot_point.x, pivot_point.y)      
    
            
    if start_point is not None and end_point is not None:
        if start_point.lat is not None and start_point.lon is not None:
            xy_coords = to_xy_coords(start_point, wrfin, timeidx, 
                                     stagger, projection)
            start_point_xy = (xy_coords.x, xy_coords.y)
        else:
            start_point_xy = (start_point.x, start_point.y)
            
        if end_point.lat is not None and end_point.lon is not None:
            xy_coords = to_xy_coords(end_point, wrfin, timeidx, 
                                     stagger, projection)
            end_point_xy = (xy_coords.x, xy_coords.y)
        else:
            end_point_xy = (end_point.x, end_point.y)
    
    xy = get_xy(field2d, pivot_point_xy, angle, start_point_xy, end_point_xy)
    
    # Make a copy so we don't modify a user supplied cache
    new_cache = dict(cache) 
    new_cache["xy"] = xy
    
    # Don't modify the original args/kargs.  The args need to be a list
    # so it can be modified.
    new_args, cache_argloc = arg_location(wrapped, "cache", args, kwargs)
    new_args[cache_argloc] = new_cache
        
    result = wrapped(*new_args)
    
    # Defaults, in case the data isn't a DataArray
    outname = None
    outdimnames = None
    outcoords = None
    outattrs = None
    
    # Use XY to set the cross-section metadata
    st_x = xy[0,0]
    st_y = xy[0,1]
    ed_x = xy[-1,0]
    ed_y = xy[-1,1]
    
    cross_str = "({0}, {1}) to ({2}, {3})".format(st_x, st_y, ed_x, ed_y)
    if angle is not None:
        cross_str += " ; center={0} ; angle={1}".format(pivot_point,
                                                        angle)
    
    if isinstance(field2d, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(field2d.dims)
        outcoords.update(field2d.coords)
        for i in py3range(-2,0,1):
            outdimnames.remove(field2d.dims[i])
            try:
                del outcoords[field2d.dims[i]]
            except KeyError:
                pass # xarray 0.9
            
        # Delete any lat,lon coords
        delkeys = [key for key in viewkeys(outcoords) if is_coordvar(key)]
        for key in delkeys:
            del outcoords[key]
        
        outdimnames.append("line_idx")
        outattrs.update(field2d.attrs)
        
        outname = "{0}_line".format(field2d.name)
        
        for key in ("MemoryOrder",):
            try:
                del outattrs[key]
            except KeyError:
                pass
            
        # Interpolate to get the lat/lon coords, if desired
        if inc_latlon:
            latcoordname, loncoordname = latlon_coordvars(field2d.coords)
            
            if latcoordname is not None and loncoordname is not None:
                latcoord = field2d.coords[latcoordname]
                loncoord = field2d.coords[loncoordname]
                
                if latcoord.ndim == 2:
                    lats = _interpline(latcoord, xy)
                    lons = _interpline(loncoord, xy)
                    
                    outcoords["xy_loc"] = ("line_idx",
                                           np.asarray(tuple(
                                                CoordPair(x=xy[i,0], y=xy[i,1],
                                                    lat=lats[i], lon=lons[i]) 
                                           for i in py3range(xy.shape[-2])))
                                           )

                # Moving domain 
                else:
                    extra_dims = latcoord.shape[0:-2]
                    outdims = extra_dims + xy.shape[-2:-1]
                    
                    latlon_loc = np.empty(outdims, np.object_)
                    for left_dims in iter_left_indexes(extra_dims):
                        idxs = left_dims + (slice(None),)
                        lats = _interpline(latcoord[idxs], xy)
                        lons = _interpline(loncoord[idxs], xy)
                        
                        latlon_loc[idxs] = np.asarray(tuple(
                                            CoordPair(x=xy[i,0], y=xy[i,1],
                                                    lat=lats[i], lon=lons[i]) 
                                            for i in py3range(xy.shape[-2]))
                                            )[:] 
                        
                    
                    extra_dimnames = latcoord.dims[0:-2]
                    loc_dimnames = extra_dimnames + ("line_idx",)
                    outcoords["xy_loc"] = (loc_dimnames, latlon_loc)
                    
            else:
                outcoords["xy_loc"] = ("line_idx", np.asarray(tuple(
                                                CoordPair(xy[i,0], xy[i,1]) 
                                          for i in py3range(xy.shape[-2]))))
            
        else:    
            outcoords["xy_loc"] = ("line_idx", np.asarray(tuple(
                                                CoordPair(xy[i,0], xy[i,1]) 
                                          for i in py3range(xy.shape[-2]))))
        
    else:
        outname = "field2d_line"
        outattrs = OrderedDict()
    
    outattrs["orientation"] = cross_str
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs) 
    

def _set_vinterp_meta(wrapped, instance, args, kwargs):
    """A decorator implementation that sets the metadata for a wrapped 
    vertical coordinate interpolation function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        wrapped: The wrapped function which in turns needs to be called by your 
            wrapper function.
        
        instance: The object to which the wrapped function was bound when it 
            was called.
        
        args: The list of positional arguments supplied when the decorated 
            function was called.
        
        kwargs: The dictionary of keyword arguments supplied when the decorated 
            function was called.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        vertical coordinate interpolation function output with or without 
        metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    See Also:
    
        :mod:`wrapt`
        
    """      
    argvars = from_args(wrapped, ("wrfin", "field", "vert_coord", 
                                  "interp_levels", "extrapolate",
                                  "field_type", "log_p",
                                  "timeidx", "method", "squeeze",
                                  "cache"), 
                          *args, **kwargs)  
    
    field = argvars["field"]
    vert_coord = argvars["vert_coord"]
    interp_levels = argvars["interp_levels"]
    field_type = argvars["field_type"]
    
    result = wrapped(*args, **kwargs)
    
    # Defaults, in case the data isn't a DataArray
    outname = None
    outdimnames = None
    outcoords = None
    outattrs = None
    
    
    if isinstance(field, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(field.dims)
        outcoords.update(field.coords)
        
        outdimnames.remove(field.dims[-3])
        try:
            del outcoords[field.dims[-3]]
        except KeyError:
            pass # xarray 0.9
        
        outdimnames.insert(-2, "interp_level")
        outcoords["interp_level"] = interp_levels
        outattrs.update(field.attrs)
        outattrs["vert_interp_type"] = vert_coord
        
        outname = field.name
        
    else:
        outname = field_type
    
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs)  
      
        
def _set_2dxy_meta(wrapped, instance, args, kwargs):
    """A decorator implementation that sets the metadata for a wrapped line 
    cross section interpolation function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        wrapped: The wrapped function which in turns needs to be called by your 
            wrapper function.
        
        instance: The object to which the wrapped function was bound when it 
            was called.
        
        args: The list of positional arguments supplied when the decorated 
            function was called.
        
        kwargs: The dictionary of keyword arguments supplied when the decorated 
            function was called.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        line cross section interpolation function output with or without 
        metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    See Also:
    
        :mod:`wrapt`
        
    """   
    argvars = from_args(wrapped, ("field3d", "xy"), *args, **kwargs)  
    
    field3d = argvars["field3d"]
    xy = argvars["xy"]
    xy = to_np(xy)
    
    result = wrapped(*args, **kwargs)
    
    # Use XY to set the cross-section metadata
    st_x = xy[0,0]
    st_y = xy[0,1]
    ed_x = xy[-1,0]
    ed_y = xy[-1,1]
    
    cross_str = "({0},{1}) to ({2},{3})".format(st_x, st_y, 
                                                ed_x, ed_y)
    
    outname = None
    outdimnames = None
    outcoords = None
    outattrs = None
    
    # Dims are (...,xy,z)
    if isinstance(field3d, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(field3d.dims)
        outcoords.update(field3d.coords)
        
        for i in py3range(-2,0,1):
            try:
                del outcoords[field3d.dims[i]]
            except KeyError:
                pass # xarray 0.9
            outdimnames.remove(field3d.dims[i])
            
        # Need to remove XLAT, XLONG...
        delkeys = (key for key,arr in viewitems(field3d.coords) 
                   if arr.ndim > 1)
        
        for key in delkeys:
            del outcoords[key]
        
        outdimnames.append("line_idx")
        #outattrs.update(field3d.attrs)
        
        desc = field3d.attrs.get("description", None)
        if desc is not None:
            outattrs["description"] = desc
            
        units = field3d.attrs.get("units", None)
        if units is not None:
            outattrs["units"] = units
        
        outname = "{0}_2dxy".format(field3d.name)
        
        outcoords["xy_loc"] = ("line_idx", [CoordPair(xy[i,0], xy[i,1]) 
                           for i in py3range(xy.shape[-2])])
        
        for key in ("MemoryOrder",):
            try:
                del outattrs[key]
            except KeyError:
                pass
        
    else:
        outname = "field3d_2dxy"
    
    outattrs["orientation"] = cross_str
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs) 


def _set_1d_meta(wrapped, instance, args, kwargs):
    """A decorator implementation that sets the metadata for a wrapped 1D 
    interpolation function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        wrapped: The wrapped function which in turns needs to be called by your 
            wrapper function.
        
        instance: The object to which the wrapped function was bound when it 
            was called.
        
        args: The list of positional arguments supplied when the decorated 
            function was called.
        
        kwargs: The dictionary of keyword arguments supplied when the decorated 
            function was called.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        1D interpolation function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    See Also:
    
        :mod:`wrapt`
        
    """   
    argvars = from_args(wrapped, ("field", "z_in", "z_out", "missingval"), 
                          *args, **kwargs)  
    
    field = argvars["field"]
    z_in = argvars["z_in"]
    z_out = argvars["z_out"]
    missingval = argvars["missingval"]
    
    result = wrapped(*args, **kwargs)
    
    outname = None
    outdimnames = None
    outcoords = None
    outattrs = None
    
    # Dims are (...,xy,z)
    if isinstance(field, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(field.dims)
        
        outdimnames.pop(-1)
        
        for name in outdimnames:
            try:
                outcoords[name] = field.coords[name]
            except KeyError:
                continue
        
        outdimnames.append("z")
        outname = "{0}_z".format(field.name)
        outcoords["z"] = z_out
        
        outattrs["_FillValue"] = missingval
        outattrs["missing_value"] = missingval
        
        desc = field.attrs.get("description", None)
        if desc is not None:
            outattrs["description"] = desc
            
        units = field.attrs.get("units", None)
        if units is not None:
            outattrs["units"] = units
        
    else:
        outname = "field_z"
    
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs)
    

def _set_xy_meta(wrapped, instance, args, kwargs):
    """A decorator implementation that sets the metadata for a wrapped xy line 
    interpolation function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        wrapped: The wrapped function which in turns needs to be called by your 
            wrapper function.
        
        instance: The object to which the wrapped function was bound when it 
            was called.
        
        args: The list of positional arguments supplied when the decorated 
            function was called.
        
        kwargs: The dictionary of keyword arguments supplied when the decorated 
            function was called.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        xy line interpolation function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    See Also:
    
        :mod:`wrapt`
        
    """   
    argvars = from_args(wrapped, ("field", "pivot_point", "angle", 
                                  "start_point", "end_point"), 
                        *args, **kwargs)  
    
    field = argvars["field"]
    pivot_point = argvars["pivot_point"]
    angle = argvars["angle"]
    start_point = argvars["start_point"]
    end_point = argvars["end_point"]
    
    result = wrapped(*args, **kwargs)
    
    if isinstance(field, DataArray):
        outname = "{0}_xy".format(field.name)
    else:
        outname = "xy"
    
    outdimnames = ["line_idx", "x_y"]
    outcoords = OrderedDict()
    outattrs = OrderedDict()
    
    outcoords["x_y"] = ["x", "y"]
    
    if pivot_point is not None and angle is not None:
        outattrs["pivot_point"] = pivot_point
        outattrs["angle"] = angle
        
    if start_point is not None and end_point is not None:
        outattrs["start_point"] = start_point
        outattrs["end_point"] = end_point
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs) 
    
    
def set_interp_metadata(interp_type):
    """A decorator that sets the metadata for a wrapped interpolation 
    function.
    
    If the wrapped function's *meta* argument is False, then this decorator 
    returns the wrapped function output without applying the metadata.
    
    Args:
    
        interp_type (:obj:`str`): The type of interpolation routine.  Choices
            are: 'horiz', 'cross', 'line', 'vinterp', '2dxy', '1d', 
            'xy'.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        interpolation function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """   
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        if interp_type == "horiz":
            return _set_horiz_meta(wrapped, instance, args, kwargs)
        elif interp_type == "cross":
            return _set_cross_meta(wrapped, instance, args, kwargs)
        elif interp_type == "line":
            return _set_line_meta(wrapped, instance, args, kwargs)
        elif interp_type == "vinterp":
            return _set_vinterp_meta(wrapped, instance, args, kwargs)
        elif interp_type == "2dxy":
            return _set_2dxy_meta(wrapped, instance, args, kwargs)
        elif interp_type == "1d":
            return _set_1d_meta(wrapped, instance, args, kwargs)
        elif interp_type == "xy":
            return _set_xy_meta(wrapped, instance, args, kwargs)
        
    return func_wrapper


def set_alg_metadata(alg_ndims, refvarname, 
                     refvarndims=None, missingarg=None,
                     stagdim=None, stagsubvar=None,
                     units=None, description=None):
    """A decorator that sets the metadata for a wrapped raw diagnostic 
    function.
    
    Args:
    
        alg_ndims (:obj:`int`): The number of dimensions returned by the 
            wrapped function.
            
        refvarname (:obj:`str`): The wrapped function argument name for the 
            reference variable.
        
        refvarndims (:obj:`int`, optional): The number of right dimensions for 
            the reference variable.  This paramter is required when the 
            wrapped function result has less dimensions than reference. 
            Default is None.
        
        missingarg (:obj:`str`, optional): The wrapped function argument name 
            for the missing value variable.  Default is None.
        
        stagdim (:obj`int`, optional): The staggered dimension for the 
            reference.  This is only needed if the reference variable is 
            on a staggered grid.  Default is None.
        
        stagsubvar (:obj:`str`, optional): The wrapped function argument name 
            to use to supply the unstaggered dimension name for the staggered 
            dimension in the reference. This is needed if *stagdim* is not 
            None. It is primarily used for absolute vorticity. Default is None.
            
        units (:obj:`str`, optional): The units to use if if there is no 
            'units' argument for the wrapped function.  Default is None.
            
        description (:obj:`str`, optional): A description for the wrapped 
            algorithm, which is stored in the :attr:`xarray.DataArray.attrs` 
            attribute under the 'description' key.  Default is None.
            
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        numerical function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
            
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
            
        
        result = wrapped(*args, **kwargs)
        
        outname = wrapped.__name__
        outattrs = OrderedDict()
        
        # Default dimension names
        outdims = ["dim_{}".format(i) for i in py3range(result.ndim)]
        
        if missingarg is not None:
            missingval = from_args(wrapped, (missingarg,), 
                                   *args, **kwargs)[missingarg]
        else:
            missingval = None
        
        if missingval is not None:
            outattrs["_FillValue"] = missingval
            outattrs["missing_value"] = missingval
            result = np.ma.masked_values(result, missingval)
        
        if units is not None:
            if isinstance(units, from_var):
                _units = units(wrapped, *args, **kwargs)
                if _units is not None:
                    outattrs["units"] = _units 
            else:
                outattrs["units"] = units
                
        else:
            # Check for a units argument, if not, just ignore
            _units = from_args(wrapped, ("units",), *args, **kwargs)["units"]
            if _units is not None:
                outattrs["units"] = _units
            
            
        if description is not None:
            if isinstance(description, from_var):
                desc = description(wrapped, *args, **kwargs)
                if desc is not None:
                    outattrs["description"] = desc
            else:
                outattrs["description"] = description
        
        
        # Copy the dimnames from the reference variable, otherwise, use
        # the supplied dimnames
        if refvarname is not None:
            refvar = from_args(wrapped, (refvarname,), 
                               *args, **kwargs)[refvarname]
        else:
            refvar = None
            
        if stagsubvar is not None:
            stagvar = from_args(wrapped, (stagsubvar,), 
                               *args, **kwargs)[stagsubvar]
        else:
            stagvar = None
            
        if isinstance(refvar, DataArray):
            
            # Copy the right dims
            outdims[-alg_ndims:] = refvar.dims[-alg_ndims:]
            
            # Use the stagsubvar if applicable
            if stagvar is not None and stagdim is not None:
                outdims[stagdim] = stagvar.dims[stagdim]
            
            # Left dims 
            if refvarndims is None:
                # Used when result and reference are aligned on right
                if result.ndim > alg_ndims:
                    result_extra = result.ndim - alg_ndims
                    
                    for i in py3range(1, result_extra + 1):
                        idx = -alg_ndims - i
                        if -idx <= refvar.ndim:
                            outdims[idx] = refvar.dims[idx]
                        else:
                            continue
            # When reference and result aren't exactly aligned (slp)
            # (reference is 3D, result is 2D)
            else: 
                ref_extra = refvar.ndim - refvarndims
                ref_left_dimnames = refvar.dims[0:ref_extra]
                
                for i,dimname in enumerate(ref_left_dimnames[::-1], 1):
                    if i <= result.ndim:
                        outdims[-alg_ndims - i] = dimname
                    else:
                        continue
                    
        out = DataArray(result, name=outname, dims=outdims, attrs=outattrs)
        
        return out
    
    return func_wrapper

def set_smooth_metdata():
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs): 
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("field", "passes"), 
                            *args, **kwargs)
        
        field = argvars["field"]
        passes = argvars["passes"]
        
        result = wrapped(*args, **kwargs)
        
        outname = "smooth2d"
        outdimnames = ["dim_{}".format(i) for i in py3range(result.ndim)]
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        
        if isinstance(field, DataArray):
            outname = "smooth_" + ucode(field.name)
            outdimnames = list(field.dims)
            outcoords.update(field.coords)
            outattrs.update(field.attrs)
            outattrs["passes"] = passes
                
        if isinstance(result, ma.MaskedArray):
            outattrs["_FillValue"] = result.fill_value
            outattrs["missing_value"] = result.fill_value
        
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
    
    return func_wrapper
    
    
def set_uvmet_alg_metadata(units=None, description="earth rotated u,v",
                           latarg="lat", windarg="u"):
    """A decorator that sets the metadata for the wrapped raw UVMET diagnostic 
    function.
    
    Args:
            
        units (:obj:`str`, optional): The units to use if if there is no 
            'units' argument for the wrapped function.  Default is None.
            
        description (:obj:`str`, optional): A description for the wrapped 
            algorithm, which is stored in the :attr:`xarray.DataArray.attrs` 
            attribute under the 'description' key.  Default is None.
            
        latarg (:obj:'str`, optional): The wrapped function argument name for 
            latitude.  Default is 'lat'.
            
        windarg (:obj:`str`, optional): The wrapped function argument name for 
            the u wind component.  Default is 'u'.
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        UVMET function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
            
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
            
        result = wrapped(*args, **kwargs)
        
        # Default dimension names
        outdims = ["dim_{}".format(i) for i in py3range(result.ndim)]
        
        outname = "uvmet"
        outattrs = OrderedDict()
        
        if units is not None:
            outattrs["units"] = units
        else:
            _units = from_args(wrapped, ("units",), *args, **kwargs)["units"]
            if _units is not None:
                outattrs["units"] = _units
            
        if description is not None:
            outattrs["description"] = description
        
        latvar = from_args(wrapped, latarg, *args, **kwargs)[latarg]
        uvar = from_args(wrapped, windarg, *args, **kwargs)[windarg]
            
        if isinstance(uvar, DataArray) and isinstance(latvar, DataArray):
            # Right dims come from latvar
            outdims[-2:] = latvar.dims[-2:]
            
            # Left dims come from u-var
            outdims[1:-2] = uvar.dims[0:-2]
            
        # Left-most is always u_v
        outdims[0] = "u_v"
        outcoords = {}
        outcoords["u_v"] = ["u", "v"]
                          
        out = DataArray(result, name=outname, dims=outdims, coords=outcoords,
                        attrs=outattrs)
        
        return out
    
    return func_wrapper


def set_cape_alg_metadata(is2d, copyarg="pres_hpa"):
    """A decorator that sets the metadata for the wrapped raw CAPE diagnostic 
    function.
    
    Args:
            
        is2d (:obj:`bool`): Set to True for the two-dimensional CAPE 
            calculation, False for three-dimensional CAPE.
            
        copyarg (:obj:`str`): The wrapped function argument to use for 
            copying dimension names.  Default is 'pres_hpa'.
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        CAPE function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
            
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
    
        result = wrapped(*args, **kwargs)
        
        argvals = from_args(wrapped, (copyarg,"missing"), *args, **kwargs)
        p = argvals[copyarg]
        missing = argvals["missing"]
        
        # Note: cape_3d supports using only a single column of data
        is1d = p.ndim == 1
        
        # Need to squeeze the right dimensions for 1D cape
        if is1d:
            result = np.squeeze(result)
        
        # Default dimension names
        outdims = ["dim_{}".format(i) for i in py3range(result.ndim)]
        
        outattrs = OrderedDict()
        
        if is2d:
            outname = "cape_2d"
            outattrs["description"] = "mcape ; mcin ; lcl ; lfc"
            outattrs["units"] = "J kg-1 ; J kg-1 ; m ; m"
            outattrs["MemoryOrder"] = "XY"
        else:
            if not is1d:
                outname = "cape_3d"
                outattrs["description"] = "cape; cin"
                outattrs["units"] = "J kg-1 ; J kg-1"
                outattrs["MemoryOrder"] = "XYZ"
            else:
                outname = "cape_column"
                outattrs["description"] = "cape; cin"
                outattrs["units"] = "J kg-1 ; J kg-1"
                outattrs["MemoryOrder"] = "Z"
                
        
        if isinstance(p, DataArray):
            if is2d:
                # Right dims
                outdims[-2:] = p.dims[-2:]
                # Left dims
                outdims[1:-2] = p.dims[0:-3]
                
            else:
                if not is1d:
                    # Right dims
                    outdims[-3:] = p.dims[-3:]
                    # Left dims
                    outdims[1:-3] = p.dims[0:-3]
                else:
                    outdims[1] = p.dims[0]
        
        outcoords = {}     
        # Left-most is always cape_cin or cape_cin_lcl_lfc
        if is2d:
            outdims[0] = "cape_cin_lcl_lfc"
            outcoords["cape_cin_lcl_lfc"] = ["cape", "cin", "lcl", "lfc"]
        else:
            outdims[0] = "cape_cin"
            outcoords["cape_cin"] = ["cape", "cin"]
            
        outattrs["_FillValue"] = missing
        outattrs["missing_value"] = missing
            
        out = DataArray(result, name=outname, dims=outdims, coords=outcoords,
                        attrs=outattrs)
        
        return out
    
    return func_wrapper


def set_cloudfrac_alg_metadata(copyarg="pres"):
    """A decorator that sets the metadata for the wrapped raw cloud fraction 
    diagnostic function.
    
    Args:
            
        copyarg (:obj:`str`): The wrapped function argument to use for 
            copying dimension names.  Default is 'pres'.
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        cloud fraction function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
            
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
            
        result = wrapped(*args, **kwargs)
        
        # Default dimension names
        outdims = ["dim_{}".format(i) for i in py3range(result.ndim)]
        
        outattrs = OrderedDict()
        
        outname = "cloudfrac"
        outattrs["description"] = "low, mid, high clouds"
        outattrs["units"] = "%"
        outattrs["MemoryOrder"] = "XY"
            
        
        p = from_args(wrapped, copyarg, *args, **kwargs)[copyarg]
            
        if isinstance(p, DataArray):
            # Right dims
            outdims[-2:] = p.dims[-2:]
            # Left dims
            outdims[1:-2] = p.dims[0:-3]
                
        
        outcoords = {}     
        # Left-most is always low_mid_high
        outdims[0] = "low_mid_high"
        outcoords["low_mid_high"] = ["low", "mid", "high"]
            
        out = DataArray(result, name=outname, dims=outdims, coords=outcoords,
                        attrs=outattrs)
        
        return out
    
    return func_wrapper


def set_destag_metadata():
    """A decorator that sets the metadata for the wrapped raw destaggering 
    function.
    
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The wrapped 
        destaggering function output with or without metadata.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    """
    
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
            
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
            
        result = wrapped(*args, **kwargs)
        
        # Default dimension names
        outdims = ["dim{}".format(i) for i in py3range(result.ndim)]
        
        destag_args = from_args(wrapped, ("var", "stagger_dim"), 
                                *args, **kwargs)
        var = destag_args["var"]
        destag_dim = destag_args["stagger_dim"] 
            
        if isinstance(var, DataArray):
            
            if var.name is not None:
                outname = "destag_{}".format(var.name)
            else:
                outname = "destag_var"
            
            outattrs = OrderedDict()
            outattrs.update(var.attrs)
            
            outattrs["destag_dim"] = destag_dim
            
            outdims = []
            outdims += var.dims
            
            destag_dim_name = outdims[destag_dim]
            if destag_dim_name.find("_stag") >= 0:
                new_dim_name = destag_dim_name.replace("_stag", "")
            else:
                if destag_dim >= 0:
                    new_dim_name = "dim_{}".format(destag_dim)
                else:
                    dim_num = result.ndim + destag_dim
                    new_dim_name = "dim_{}".format(dim_num)
            
            outdims[destag_dim] = new_dim_name
                          
        out = DataArray(result, name=outname, dims=outdims, attrs=outattrs)
        
        return out
    
    return func_wrapper





