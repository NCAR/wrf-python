from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
import wrapt 
from collections import OrderedDict

import numpy as np
import numpy.ma as ma


from .util import (viewkeys, viewitems, extract_vars, 
                   combine_with, either, from_args, arg_location,
                   _is_coord_var, CoordPair, npvalues, py3range, ucode)
from .interputils import get_xy_z_params, get_xy
from .latlonutils import ij_to_ll, ll_to_ij
from .config import xarray_enabled

if xarray_enabled():
    from xarray import DataArray
    
__all__ = ["copy_and_set_metadata", "set_wind_metadata",
           "set_latlon_metadata", "set_height_metadata",
           "set_interp_metadata"]

def copy_and_set_metadata(copy_varname=None, delete_attrs=None, name=None,
                          remove_dims=None, dimnames=None, 
                          coords=None, **fixed_attrs):
    """Decorator to set the metadata for a WRF method.
    
    A cache is inserted/updated to include the extracted variable that will 
    have its metadata copied. This prevents the variable being extracted more 
    than once.  This extraction can be slow with sequences of large files.
    
    """
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs): 
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        argvars = from_args(wrapped, ("wrfnc", "timeidx", "method", 
                                      "squeeze", "cache", "units", "meta"), 
                            *args, **kwargs)
        
        wrfnc = argvars["wrfnc"]
        timeidx = argvars["timeidx"]
        units = argvars["units"]
        method = argvars["method"]
        squeeze = argvars["squeeze"]
        cache = argvars["cache"]
        if cache is None:
            cache = {}
        
        # Note:  can't modify nonlocal var
        if (callable(copy_varname)):
            _copy_varname = copy_varname(wrfnc)
        else:
            _copy_varname = copy_varname
        
        # Extract the copy_from argument
        var_to_copy = None if cache is None else cache.get(_copy_varname, 
                                                           None)

            
        if var_to_copy is None:
            var_to_copy = extract_vars(wrfnc, timeidx, (_copy_varname,), 
                                       method, squeeze, cache,
                                       meta=True)[_copy_varname]
        
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
                if isinstance(dimnames, combine_with):
                    outdimnames, outcoords = dimnames(var_to_copy)
                else:
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
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
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
        
        if isinstance(copy_varname, either):
            _copy_varname = copy_varname(wrfnc)
        else:
            _copy_varname = copy_varname
        
        copy_var = extract_vars(wrfnc, timeidx, _copy_varname, 
                                method, squeeze, cache, 
                                meta=True)[_copy_varname]
        
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
        outcoords.update(copy_var.coords)
        outattrs.update(copy_var.attrs)
        
        if wind_ncvar:
            pass
        
        elif not wspd_wdir:
            if not two_d:
                outdimnames.insert(-3, "u_v")
            else:
                outdimnames.insert(-2, "u_v")
                outattrs["MemoryOrder"] = "XY"
            outcoords["u_v"] = ["u", "v"]
        else:
            if not two_d:
                outdimnames.insert(-3, "wspd_wdir")
            else:
                outdimnames.insert(-2, "wspd_wdir")
                outattrs["MemoryOrder"] = "XY"
                
            outcoords["wspd_wdir"] = ["wspd", "wdir"]
        
        if units is not None: 
            outattrs["units"] = units
            
        outname = name
        outattrs["description"] = description
        
        return DataArray(result, name=outname, coords=outcoords, 
                       dims=outdimnames, attrs=outattrs)
        
    return func_wrapper

def set_latlon_metadata(ij=False):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
            return wrapped(*args, **kwargs)
        
        res = wrapped(*args, **kwargs)
        
        # Want to preserve the input coordinate pair in metadata
        if res.ndim == 1:
            res = res[np.newaxis, :]
        
        argnames = ["i", "j"] if not ij else ["latitude", "longitude"]
        argnames.append("squeeze")
        outname = "latlon" if not ij else "ij"
        
        if res.ndim == 2:
            dimnames = (["ij", "lat_lon"] if not ij 
                        else ["latlon", "i_j"])
        else:
            dimnames = (["ij", "domain", "lat_lon"] if not ij 
                        else ["latlon", "domain", "i_j"])
        
        argvars = from_args(wrapped, argnames, *args, **kwargs)
        
        var1 = argvars[argnames[0]]
        var2 = argvars[argnames[1]]
        squeeze = argvars["squeeze"]
        
        arr1 = np.asarray(var1).ravel()
        arr2 = np.asarray(var2).ravel()
        
        coords = {}
        if not ij:
            coords["coord_pair"] = (dimnames[0], [CoordPair(i=x[0], j=x[1]) 
                               for x in zip(arr1, arr2)])
            coords[dimnames[-1]] = ["lat", "lon"]
        else:
            coords["coord_pair"] = (dimnames[0], [CoordPair(lat=x[0], lon=x[1]) 
                               for x in zip(arr1, arr2)])
            coords[dimnames[-1]] = ["i", "j"]
        
        da = DataArray(res, name=outname, dims=dimnames, coords=coords)
        
        if squeeze:
            da = da.squeeze()
        
        return da
    
    return func_wrapper
    
def set_height_metadata(geopt=False):
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        do_meta = from_args(wrapped, ("meta",), *args, **kwargs)["meta"]
        
        if do_meta is None:
            do_meta = True
                
        if not xarray_enabled() or not do_meta:
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
                              method, squeeze, cache, meta=True)
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
    argvars = from_args(wrapped, ("field3d", "z", "desiredloc", 
                                  "missingval"), 
                          *args, **kwargs)  
    
    field3d = argvars["field3d"]
    z = argvars["z"]
    desiredloc = argvars["desiredloc"]
    missingval = argvars["missingval"]
    
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
        del outcoords[field3d.dims[-3]]
        outattrs.update(field3d.attrs)
        outname = "{0}_{1}".format(field3d.name, name_levelstr)
        
    else:
        outname = "field3d_{0}".format(levelstr)
        outattrs = OrderedDict()
        
    outattrs["PlotLevelID"] = levelstr
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
    argvars = from_args(wrapped, ("field3d", "z", "missingval", 
                                  "pivot_point", "angle",
                                  "start_point", "end_point",
                                  "cache"), 
                          *args, **kwargs)  
    
    field3d = argvars["field3d"]
    z = argvars["z"]
    missingval = argvars["missingval"]
    pivot_point = argvars["pivot_point"]
    angle = argvars["angle"]
    start_point = argvars["start_point"]
    end_point = argvars["end_point"]
    cache = argvars["cache"]
    
    xy, var2dz, z_var2d = get_xy_z_params(npvalues(z), pivot_point, angle,
              start_point, end_point)
    
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
    
    cross_str = "cross-section: ({0}, {1}) to ({2}, {3})".format(st_x, st_y, 
                                                               ed_x, ed_y)
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
            del outcoords[field3d.dims[i]]
        
        
        # Delete any lat,lon coords
        delkeys = [key for key in viewkeys(outcoords) if _is_coord_var(key)]
        for key in delkeys:
            del outcoords[key]
        
        outdimnames.append("vertical")
        outdimnames.append("xy")
        outattrs.update(field3d.attrs)
        
        outname = "{0}_cross".format(field3d.name)
        
        for key in ("MemoryOrder",):
            try:
                del outattrs[key]
            except KeyError:
                pass
            
        outcoords["xy_loc"] = ("xy", [CoordPair(xy[i,0], xy[i,1]) 
                           for i in py3range(xy.shape[-2])])
        
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
    argvars = from_args(wrapped, ("field2d", "pivot_point", "angle",
                                  "start_point", "end_point", "cache"), 
                          *args, **kwargs)  
    
    field2d = argvars["field2d"]
    pivot_point = argvars["pivot_point"]
    angle = argvars["angle"]
    start_point = argvars["start_point"]
    end_point = argvars["end_point"]
    cache = argvars["cache"]
    
    if cache is None:
        cache = {}
    
    xy = get_xy(field2d, pivot_point, angle, start_point, end_point)
    
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
    
    cross_str = "({0}, {1}) to ({2}, {3})".format(st_x, st_y, 
                                                ed_x, ed_y)
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
            del outcoords[field2d.dims[i]]
            
        # Delete any lat,lon coords
        delkeys = [key for key in viewkeys(outcoords) if _is_coord_var(key)]
        for key in delkeys:
            del outcoords[key]
        
        outdimnames.append("xy")
        outattrs.update(field2d.attrs)
        
        outname = "{0}_line".format(field2d.name)
        
        for key in ("MemoryOrder",):
            try:
                del outattrs[key]
            except KeyError:
                pass
            
        outcoords["xy_loc"] = ("xy", [CoordPair(xy[i,0], xy[i,1]) 
                           for i in py3range(xy.shape[-2])])
        
    else:
        outname = "field2d_line"
        outattrs = OrderedDict()
    
    outattrs["orientation"] = cross_str
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs) 
    

def _set_vinterp_meta(wrapped, instance, args, kwargs):
    argvars = from_args(wrapped, ("wrfnc", "field", "vert_coord", 
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
        del outcoords[field.dims[-3]]
        
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
    argvars = from_args(wrapped, ("field3d", "xy"), 
                          *args, **kwargs)  
    
    field3d = argvars["field3d"]
    xy = argvars["xy"]
    
    result = wrapped(*args, **kwargs)
    
    # Use XY to set the cross-section metadata
    st_x = xy[0,0]
    st_y = xy[0,1]
    ed_x = xy[-1,0]
    ed_y = xy[-1,1]
    
    cross_str = "({0},{1}) to ({2},{3})".format(st_x, st_y, 
                                                ed_x, ed_y)
    
    # Dims are (...,xy,z)
    if isinstance(field3d, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(field3d.dims)
        outcoords.update(field3d.coords)
        for i in py3range(-2,0,1):
            outdimnames.remove(field3d.dims[i])
            del outcoords[field3d.dims[i]]
        
        outdimnames[-2] = "xy"
        outattrs.update(field3d.attrs)
        
        outname = "{0}_xy".format(field3d.name)
        
        outcoords["xy_loc"] = ("xy", [CoordPair(xy[i,0], xy[i,1]) 
                           for i in py3range(xy.shape[-2])])
        
        for key in ("MemoryOrder",):
            try:
                del outattrs[key]
            except KeyError:
                pass
        
    else:
        outname = "field3d_xy"
    
    outattrs["Orientation"] = cross_str
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs) 


def _set_1d_meta(wrapped, instance, args, kwargs):
    argvars = from_args(wrapped, ("v_in", "z_in", "z_out", "missingval"), 
                          *args, **kwargs)  
    
    v_in = argvars["v_in"]
    z_in = argvars["z_in"]
    z_out = argvars["z_out"]
    missingval = argvars["missingval"]
    
    result = wrapped(*args, **kwargs)
    
    # Dims are (...,xy,z)
    if isinstance(v_in, DataArray):
        outcoords = OrderedDict()
        outattrs = OrderedDict()
        outdimnames = list(v_in.dims)
        outcoords.update(v_in.coords)
        
        
        outdimnames.remove(v_in.dims[-1])
        del outcoords[v_in.dims[-1]]
        outdimnames.append("z")
        outname = "{0}_z".format(v_in.name)
        outcoords["z"] = z_out
        
        outattrs.update(v_in.attrs)
        outattrs["_FillValue"] = missingval
        outattrs["missing_value"] = missingval
        
    else:
        outname = "v_in_z"
    
    
    return DataArray(result, name=outname, dims=outdimnames, 
                     coords=outcoords, attrs=outattrs)
    
    
def set_interp_metadata(interp_type):
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
        
    return func_wrapper
