from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import os
from sys import version_info
from copy import copy
from collections import Iterable, Mapping, OrderedDict
from itertools import product
from types import GeneratorType
import datetime as dt
from math import floor, copysign
from inspect import getmodule

try:
    from inspect import signature
except ImportError:
    from inspect import getargspec

try:    
    from inspect import getargvalues
except ImportError:
    from inspect import getgeneratorlocals

import numpy as np
import numpy.ma as ma

from .config import xarray_enabled
from .projection import getproj, NullProjection
from .constants import Constants, ALL_TIMES
from .py3compat import (viewitems, viewkeys, viewvalues, isstr, py2round, 
                        py3range, ucode)

if xarray_enabled():
    from xarray import DataArray
    from pandas import NaT


_COORD_PAIR_MAP = {"XLAT" : ("XLAT", "XLONG"),
                "XLONG" : ("XLAT", "XLONG"),
                "XLAT_M" : ("XLAT_M", "XLONG_M"),
                "XLONG_M" : ("XLAT_M", "XLONG_M"),
                "XLAT_U" : ("XLAT_U", "XLONG_U"),
                "XLONG_U" : ("XLAT_U", "XLONG_U"),
                "XLAT_V" : ("XLAT_V", "XLONG_V"),
                "XLONG_V" : ("XLAT_V", "XLONG_V"),
                "CLAT" : ("CLAT", "CLONG"),
                "CLONG" : ("CLAT", "CLONG")}


_COORD_VARS = ("XLAT", "XLONG", "XLAT_M", "XLONG_M", "XLAT_U", "XLONG_U",
               "XLAT_V", "XLONG_V", "CLAT", "CLONG")

_LAT_COORDS = ("XLAT",  "XLAT_M", "XLAT_U",  "XLAT_V",  "CLAT")

_LON_COORDS = ("XLONG", "XLONG_M", "XLONG_U","XLONG_V", "CLONG")

_TIME_COORD_VARS = ("XTIME",)


def is_time_coord_var(varname):
    return varname in _TIME_COORD_VARS


def get_coord_pairs(varname):
    return _COORD_PAIR_MAP[varname]


def is_multi_time_req(timeidx):
    return timeidx is None


def is_multi_file(wrfnc):
    return (isinstance(wrfnc, Iterable) and not isstr(wrfnc))

def has_time_coord(wrfnc):
    return "XTIME" in wrfnc.variables

def is_mapping(wrfnc):
    return isinstance(wrfnc, Mapping)

def _generator_copy(gen):
    funcname = gen.__name__
    try:
        argvals = getargvalues(gen.gi_frame)
    except NameError:
        argvals = getgeneratorlocals(gen)
    module = getmodule(gen.gi_frame)
    
    if module is not None:
        res = module.get(funcname)(**argvals.locals)
    else:
        # Created in jupyter or the python interpreter
        import __main__
        res = getattr(__main__, funcname)(**argvals.locals)
        
    return res

def test():
    q = [1,2,3]
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

def latlon_coordvars(d):
    lat_coord = None
    lon_coord = None
    
    for name in _LAT_COORDS:
        if name in viewkeys(d):
            lat_coord = name
            break
        
    for name in _LON_COORDS:
        if name in viewkeys(d):
            lon_coord = name
            break
        
    return lat_coord, lon_coord

def is_coordvar(varname):
    return varname in _COORD_VARS

class IterWrapper(Iterable):
    """A wrapper class for generators and custom iterable classes which returns
    a new iterator from the start of the sequence when __iter__ is called"""
    def __init__(self, wrapped):
        self._wrapped = wrapped
        
    def __iter__(self):
        if isinstance(self._wrapped, GeneratorType):
            return _generator_copy(self._wrapped)
        else:
            obj_copy = copy(self._wrapped)
            return obj_copy.__iter__()
            
            
def get_iterable(wrfseq):
    """Returns a resetable iterable object."""
    if not is_multi_file(wrfseq):
        return wrfseq
    else:
        if not is_mapping(wrfseq):
            
            if isinstance(wrfseq, (list, tuple, IterWrapper)):
                return wrfseq
            else:
                return IterWrapper(wrfseq) # generator/custom iterable class
            
        else:
            if isinstance(wrfseq, dict):
                return wrfseq
            else:
                return dict(wrfseq) # generator/custom iterable class
    
# Helper to extract masked arrays from DataArrays that convert to NaN
def npvalues(array_type):
    if not isinstance(array_type, DataArray):
        result = array_type
    else:
        try:
            fill_value = array_type.attrs["_FillValue"]
        except KeyError:
            result = array_type.values
        else:
            result = ma.masked_invalid(array_type.values, copy=False)
            result.set_fill_value(fill_value)
    
    return result

# Helper utilities for metadata
class either(object):
    def __init__(self, *varnames):
        self.varnames = varnames
    
    def __call__(self, wrfnc):
        if is_multi_file(wrfnc):
            if not is_mapping(wrfnc):
                wrfnc = next(iter(wrfnc))
            else:
                entry = wrfnc[next(iter(viewkeys(wrfnc)))]
                return self(entry)
            
        for varname in self.varnames:
            if varname in wrfnc.variables:
                return varname
        
        raise ValueError("{} are not valid variable names".format(
                                                            self.varnames))

class combine_with(object):
    # Remove remove_idx first, then insert_idx is applied to removed set
    def __init__(self, varname, remove_dims=None, insert_before=None, 
                 new_dimnames=None, new_coords=None):
        self.varname = varname
        self.remove_dims = remove_dims
        self.insert_before = insert_before
        self.new_dimnames = new_dimnames if new_dimnames is not None else []
        self.new_coords = (new_coords if new_coords is not None 
                           else OrderedDict())
    
    def __call__(self, var):
        new_dims = list(var.dims)
        new_coords = OrderedDict(var.coords)
        
        if self.remove_dims is not None:
            for dim in self.remove_dims:
                new_dims.remove(dim)
                del new_coords[dim]
        
        if self.insert_before is not None:     
            insert_idx = new_dims.index(self.insert_before)
            new_dims = (new_dims[0:insert_idx] + self.new_dimnames + 
                    new_dims[insert_idx:])
        elif self.new_dimnames is not None:
            new_dims = self.new_dimnames
        
        if self.new_coords is not None:
            new_coords.update(self.new_coords)
        
        return new_dims, new_coords
    
# This should look like:
# [(0, (-3,-2)), # variable 1
# (1, -1)] # variable 2
class combine_dims(object):
    
    def __init__(self, pairs):
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
    def __init__(self, varname, attribute):
        self.varname = varname
        self.attribute = attribute
        
    def __call__(self, wrapped, *args, **kwargs):
        vard = from_args(wrapped, (self.varname,), *args, **kwargs)
        
        var = None
        if vard is not None:
            var = vard[self.varname]
            
        if not isinstance(var, DataArray):
            return None
        
        return var.attrs.get(self.attribute, None)
        

def _corners_moved(wrfnc, first_ll_corner, first_ur_corner, latvar, lonvar):
    lats = wrfnc.variables[latvar]
    lons = wrfnc.variables[lonvar]
    
    # Need to check all times
    for i in py3range(lats.shape[-3]):
        start_idxs = [0]*len(lats.shape) # PyNIO does not support ndim
        start_idxs[-3] = i
        start_idxs = tuple(start_idxs)
        
        end_idxs = [-1]*len(lats.shape)
        end_idxs[-3] = i
        end_idxs = tuple(end_idxs)
        
        if (first_ll_corner[0] != lats[start_idxs] or 
            first_ll_corner[1] != lons[start_idxs] or 
            first_ur_corner[0] != lats[end_idxs] or 
            first_ur_corner[1] != lons[end_idxs]):
            return True
    
    return False


def is_moving_domain(wrfseq, varname=None, latvar=either("XLAT", "XLAT_M"), 
                      lonvar=either("XLONG", "XLONG_M")):
    
    if isinstance(latvar, either):
        latvar = latvar(wrfseq)
        
    if isinstance(lonvar, either):
        lonvar = lonvar(wrfseq)
        
    # In case it's just a single file
    if not is_multi_file(wrfseq):
        wrfseq = [wrfseq]
    
    # Slow, but safe. Compare the corner points to the first item and see
    # any move.  User iterator protocol in case wrfseq is not a list/tuple.
    if not is_mapping(wrfseq):
        wrf_iter = iter(wrfseq)
        first_wrfnc = next(wrf_iter)
    else:
        entry = wrfseq[next(iter(viewkeys(wrfseq)))]
        return is_moving_domain(entry, varname, latvar, lonvar)
    
    # Quick check on pressure coordinates, bypassing the need to search the
    # domain corner points
    try:
        coord_names = getattr(first_wrfnc.variables["P"], 
                              "coordinates").split()
    except KeyError:
        pass
    else:
        if "XTIME" in coord_names:
            return True
        else:
            return False
    
    # The long way of checking all lat/lon corner points
    if varname is not None:
        try:
            coord_names = getattr(first_wrfnc.variables[varname], 
                              "coordinates").split()
        except AttributeError:
            # Variable doesn't have a coordinates attribute, use the 
            # arguments
            lon_coord = lonvar
            lat_coord = latvar
        else:
            # If the XTIME variable is found to be a coordinate variable, 
            # then it's a moving domain file
            try:
                xtime_coord = coord_names[2]
            except IndexError:
                # XTIME is not a coordinate variable, if the variable is in the
                # file, then this is not a moving domain file
                if "XTIME" in first_wrfnc.variables:
                    return False
                
            else:
                # XTIME is a coordinate, so this is a moving domain file
                if xtime_coord == "XTIME":
                    return True
                
            lon_coord = coord_names[0]
            lat_coord = coord_names[1]
    else:
        lon_coord = lonvar
        lat_coord = latvar
    
    # Probably a met_em file, need to search all the files
    lats = first_wrfnc.variables[lat_coord]
    lons = first_wrfnc.variables[lon_coord]
    
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
    
    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            if _corners_moved(wrfnc, ll_corner, ur_corner, 
                              lat_coord, lon_coord):
                return True
    
    return False

def _get_global_attr(wrfnc, attr):
    val = getattr(wrfnc, attr, None)
    
    # PyNIO puts single values in to an array
    if isinstance(val, np.ndarray):
        if len(val) == 1:
            return val[0] 
    return val
        
def extract_global_attrs(wrfnc, attrs):
    if isstr(attrs):
        attrlist = [attrs]
    else:
        attrlist = attrs
        
    multifile = is_multi_file(wrfnc)
    
    if multifile:
        if not is_mapping(wrfnc):
            wrfnc = next(iter(wrfnc))
        else:
            entry = wrfnc[next(iter(viewkeys(wrfnc)))]
            return extract_global_attrs(entry, attrs)
        
    return {attr:_get_global_attr(wrfnc, attr) for attr in attrlist}

def extract_dim(wrfnc, dim):
    if is_multi_file(wrfnc):
        if not is_mapping(wrfnc):
            wrfnc = next(iter(wrfnc))
        else:
            entry = wrfnc[next(iter(viewkeys(wrfnc)))]
            return extract_dim(entry, dim)
    
    d = wrfnc.dimensions[dim]
    if not isinstance(d, int):
        return len(d) #netCDF4
    return d # PyNIO
        
def _combine_dict(wrfdict, varname, timeidx, method, meta):
    """Dictionary combination creates a new left index for each key, then 
    does a cat or join for the list of files for that key"""
    keynames = []
    numkeys = len(wrfdict)  
    
    key_iter = iter(viewkeys(wrfdict))
    first_key = next(key_iter)
    keynames.append(first_key)
    
    is_moving = is_moving_domain(wrfdict, varname)
    
    first_array = _extract_var(wrfdict[first_key], varname, 
                              timeidx, is_moving=is_moving, method=method, 
                              squeeze=False, cache=None, meta=meta)
    
    
    # Create the output data numpy array based on the first array
    outdims = [numkeys]
    outdims += first_array.shape
    outdata = np.empty(outdims, first_array.dtype)
    outdata[0,:] = first_array[:]
    
    idx = 1
    while True:
        try:
            key = next(key_iter)
        except StopIteration:
            break
        else:
            keynames.append(key)
            vardata = _extract_var(wrfdict[key], varname, timeidx, 
                                   is_moving=is_moving, method=method, 
                                   squeeze=False, cache=None, meta=meta)
            
            if outdata.shape[1:] != vardata.shape:
                raise ValueError("data sequences must have the "
                                   "same size for all dictionary keys")
            outdata[idx,:] = npvalues(vardata)[:]
            idx += 1
      
    if xarray_enabled() and meta:
        outname = str(first_array.name)
        # Note: assumes that all entries in dict have same coords
        outcoords = OrderedDict(first_array.coords)
        outdims = ["key"] + list(first_array.dims)
        outcoords["key"] = keynames
        outattrs = OrderedDict(first_array.attrs)
        
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdims, attrs=outattrs)
    else:
        outarr = outdata
        
    return outarr

def _find_coord_names(coords):
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


def _build_data_array(wrfnc, varname, timeidx, is_moving_domain):
    multitime = is_multi_time_req(timeidx)
    time_idx_or_slice = timeidx if not multitime else slice(None)
    var = wrfnc.variables[varname]
    data = var[time_idx_or_slice, :]
    
    # Want to preserve the time dimension
    if not multitime:
        data = data[np.newaxis, :]
    
    attrs = OrderedDict(var.__dict__)
    dimnames = var.dimensions[-data.ndim:]
    
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
            
            if is_moving_domain and has_time_coord(wrfnc):
                time_coord = "XTIME"
                
        else:
            try:
                # met_em files
                stag_attr = getattr(var, "stagger")
            except AttributeError:
                lon_coord = None
                lat_coord = None
            else:
                # For met_em files, use the stagger name to get the lat/lon var
                lat_coord = "XLAT_{}".format(stag_attr)
                lon_coord = "XLONG_{}".format(stag_attr)
    else:
        coord_names = coord_attr.split()
        lon_coord = coord_names[0]
        lat_coord = coord_names[1]
        
        try:
            time_coord = coord_names[2]
        except IndexError:
            time_coord = None
    
    coords = OrderedDict()
    
    # Handle lat/lon coordinates and projection information if available
    if lon_coord is not None and lat_coord is not None:
        lon_coord_var = wrfnc.variables[lon_coord]
        lat_coord_var = wrfnc.variables[lat_coord]
        time_coord_var = (wrfnc.variables[time_coord] 
                          if time_coord is not None 
                          else None)
    
        if multitime:
            if is_moving_domain:
                # Special case with a moving domain in a multi-time file,
                # otherwise the projection parameters don't change
                coords[lon_coord] = lon_coord_var.dimensions, lon_coord_var[:]
                coords[lat_coord] = lat_coord_var.dimensions, lat_coord_var[:]
                
                # Returned lats/lons arrays will have a time dimension, so proj
                # will need to be a list due to moving corner points
                lats, lons, proj_params = get_proj_params(wrfnc, 
                                                          timeidx, 
                                                          varname)
                proj = [getproj(lats=lats[i,:], 
                                lons=lons[i,:],
                                **proj_params) for i in py3range(lats.shape[0])]
                
                if time_coord is not None:
                    coords[time_coord] = (lon_coord_var.dimensions[0], 
                                          time_coord_var[:])
            else:
                coords[lon_coord] = (lon_coord_var.dimensions[1:], 
                                     lon_coord_var[0,:])
                coords[lat_coord] = (lat_coord_var.dimensions[1:], 
                                     lat_coord_var[0,:])
                
                # Domain not moving, so just get the first time
                lats, lons, proj_params = get_proj_params(wrfnc, 0, varname)
                proj = getproj(lats=lats, lons=lons, **proj_params)
        else:
            coords[lon_coord] = (lon_coord_var.dimensions[1:], 
                                 lon_coord_var[timeidx,:])
            coords[lat_coord] = (lat_coord_var.dimensions[1:], 
                                 lat_coord_var[timeidx,:])
            lats, lons, proj_params = get_proj_params(wrfnc, 0, varname)
            proj = getproj(lats=lats, lons=lons, **proj_params)
        
        attrs["projection"] = proj
        
    
    if dimnames[0] == "Time":
        coords[dimnames[0]] = extract_times(wrfnc, timeidx)
    
    data_array = DataArray(data, name=varname, dims=dimnames, coords=coords,
                           attrs=attrs)
    
    
    return data_array


def _find_forward(wrfseq, varname, timeidx, is_moving, meta):

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
                                             is_moving)
                else:
                    return wrfnc.variables[varname][filetimeidx, :]
            else:
                comboidx += numtimes
            
    raise IndexError("invalid time index")


def _find_reverse(wrfseq, varname, timeidx, is_moving, meta):
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
                                             is_moving)
                else:
                    return wrfnc.variables[varname][filetimeidx, :]
            else:
                comboidx += numtimes
            
    raise IndexError("invalid time index")


def _find_arr_for_time(wrfseq, varname, timeidx, is_moving, meta):
    if timeidx >= 0:
        return _find_forward(wrfseq, varname, timeidx, is_moving, meta)
    else:
        return _find_reverse(wrfseq, varname, timeidx, is_moving, meta)
    
# TODO:  implement in C
def _cat_files(wrfseq, varname, timeidx, is_moving, squeeze, meta):
    if is_moving is None:
        is_moving = is_moving_domain(wrfseq, varname)
    
    file_times = extract_times(wrfseq, ALL_TIMES)
    
    multitime = is_multi_time_req(timeidx) 
    
    # For single times, just need to find the ncfile and appropriate 
    # time index, and return that array
    if not multitime:
        return _find_arr_for_time(wrfseq, varname, timeidx, is_moving, meta)

    #time_idx_or_slice = timeidx if not multitime else slice(None)
    
    # If all times are requested, need to build a new array and cat together
    # all of the arrays in the sequence
    wrf_iter = iter(wrfseq)
    
    if xarray_enabled() and meta:
        first_var = _build_data_array(next(wrf_iter), varname, 
                                      ALL_TIMES, is_moving)
    else:
        first_var = (next(wrf_iter)).variables[varname][:]
    
    outdims = [len(file_times)]
    
    # Making a new time dim, so ignore this one
    outdims += first_var.shape[1:]
    
    outdata = np.empty(outdims, first_var.dtype)
    
    numtimes = first_var.shape[0]
    startidx = 0
    endidx = numtimes
    
    outdata[startidx:endidx, :] = first_var[:]
    
    if xarray_enabled() and meta and is_moving:
        latname, lonname, timename = _find_coord_names(first_var.coords)
        outcoorddims = outdims[0:1] + outdims[-2:] 
        
        if latname is not None:
            outlats = np.empty(outcoorddims, first_var.dtype)
            outlats[startidx:endidx, :] = first_var.coords[latname][:]
            
        if lonname is not None:
            outlons = np.empty(outcoorddims, first_var.dtype)
            outlons[startidx:endidx, :] = first_var.coords[lonname][:]
            
        if timename is not None:
            outxtimes = np.empty(outdims[0])
            outxtimes[startidx:endidx] = first_var.coords[timename][:]
            
        # Projections also need to be aggregated
        outprojs = np.empty(outdims[0], np.object)
        
        outprojs[startidx:endidx] = np.asarray(first_var.attrs["projection"],
                                            np.object)[:]
    
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
            
            outdata[startidx:endidx, :] = vardata[:]
            
            if xarray_enabled() and meta and is_moving:
                if latname is not None:
                    latdata = wrfnc.variables[latname][:]
                    outlats[startidx:endidx, :] = latdata[:]
                    
                if lonname is not None:
                    londata = wrfnc.variables[lonname][:]
                    outlons[startidx:endidx, :] = londata[:]
                    
                if timename is not None:
                    xtimedata = wrfnc.variables[timename][:]
                    outxtimes[startidx:endidx] = xtimedata[:]
                    
                lats, lons, proj_params = get_proj_params(wrfnc, 
                                                          ALL_TIMES, 
                                                          varname)
                projs = [getproj(lats=lats[i,:], 
                                lons=lons[i,:],
                                **proj_params) for i in py3range(lats.shape[0])]
                
                outprojs[startidx:endidx] = np.asarray(projs, np.object)[:]
            
            startidx = endidx
    
    if xarray_enabled() and meta:
        outname = ucode(first_var.name)
        outattrs = OrderedDict(first_var.attrs)
        outcoords = OrderedDict(first_var.coords)
        outdimnames = list(first_var.dims)
        
        if "Time" not in outdimnames:
            outdimnames.insert(0, "Time")
        
        if not multitime:
            file_times = [file_times[timeidx]]
        
        outcoords[outdimnames[0]] = file_times
        
        outcoords["datetime"] = outdimnames[0], file_times
        
        # If the domain is moving, need to create the lat/lon/xtime coords
        # since they can't be copied
        if is_moving:
            outlatdims = [outdimnames[0]] + outdimnames[-2:]
            
            if latname is not None:
                outlats = outlats[:]
                outcoords[latname] = outlatdims, outlats
            if lonname is not None:
                outlons = outlons[:]
                outcoords[lonname] = outlatdims, outlons
            if timename is not None:
                outxtimes = outxtimes[:]
                outcoords[timename] = outdimnames[0], outxtimes
                    
            outattrs["projection"] = outprojs[:]
        
        outdata = outdata[:]
            
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdimnames, attrs=outattrs)
        
    else:
        outdata = outdata[:]
        outarr = outdata
        
    return outarr

def _get_numfiles(wrfseq):
    try:
        return len(wrfseq)
    except TypeError:
        wrf_iter = iter(wrfseq)
        return sum(1 for _ in wrf_iter)

# TODO:  implement in C
def _join_files(wrfseq, varname, timeidx, is_moving, meta):
    if is_moving is None:
        is_moving = is_moving_domain(wrfseq, varname)
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
        first_var = _build_data_array(wrfnc, varname, ALL_TIMES, is_moving)
        time_coord = np.full((numfiles, maxtimes), int(NaT), "datetime64[ns]")
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
    outdata = np.full(outdims, Constants.DEFAULT_FILL, first_var.dtype)
    outdata[file_idx, 0:numtimes, :] = first_var[:]
    
    # Create the secondary coordinate arrays
    if xarray_enabled() and meta and is_moving:
        latname, lonname, timename = _find_coord_names(first_var.coords)
        outcoorddims = outdims[0:2] + outdims[-2:] 
        
        if latname is not None:
            outlats = np.full(outcoorddims, Constants.DEFAULT_FILL, 
                              first_var.dtype)
            outlats[file_idx, 0:numtimes, :] = first_var.coords[latname][:]
            
        if lonname is not None:
            outlons = np.full(outcoorddims, Constants.DEFAULT_FILL, 
                              first_var.dtype)
            outlons[file_idx, 0:numtimes, :] = first_var.coords[lonname][:]
            
        if timename is not None:
            outxtimes = np.full(outdims[0:2], Constants.DEFAULT_FILL, 
                                first_var.dtype)
            outxtimes[file_idx, 0:numtimes] = first_var.coords[timename][:]
            
        # Projections also need two dimensions
        outprojs = np.full(outdims[0:2], NullProjection(), np.object)
        
        outprojs[file_idx, 0:numtimes] = np.asarray(
                                            first_var.attrs["projection"],
                                            np.object)[:]
    
    file_idx=1
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
                
            outdata[file_idx, 0:numtimes, :] = outvar[:]
            
            if xarray_enabled() and meta:
                file_times = extract_times(wrfnc, ALL_TIMES)
                time_coord[file_idx, 0:numtimes] = np.asarray(file_times, 
                                                        "datetime64[ns]")[:]
            
            if xarray_enabled() and meta and is_moving:
                if latname is not None:
                    latdata = wrfnc.variables[latname][:]
                    outlats[file_idx, 0:numtimes, :] = latdata[:]
                    
                if lonname is not None:
                    londata = wrfnc.variables[lonname][:]
                    outlons[file_idx, 0:numtimes, :] = londata[:]
                    
                if timename is not None:
                    xtimedata = wrfnc.variables[timename][:]
                    outxtimes[file_idx, 0:numtimes] = xtimedata[:]
                
                lats, lons, proj_params = get_proj_params(wrfnc, 
                                                          ALL_TIMES, 
                                                          varname)
                projs = [getproj(lats=lats[i,:], 
                                lons=lons[i,:],
                                **proj_params) for i in py3range(lats.shape[0])]
                
                outprojs[file_idx, 0:numtimes] = (
                                        np.asarray(projs, np.object)[:])
            
            # Need to update coords here
            file_idx += 1  
    
    # If any of the output files contain less than the max number of times,
    # then a mask array is needed to flag all the missing arrays with 
    # missing values
    if file_times_less_than_max:
        outdata = np.ma.masked_values(outdata, Constants.DEFAULT_FILL)
    
    if xarray_enabled() and meta:
        outname = ucode(first_var.name)
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
            outattrs["_FillValue"] = Constants.DEFAULT_FILL
            outattrs["missing_value"] = Constants.DEFAULT_FILL
            
        # If the domain is moving, need to create the lat/lon/xtime coords
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
            if timename is not None:
                outxtimes = outxtimes[:, time_idx_or_slice]
                if not multitime:
                    outxtimes = outxtimes[:, np.newaxis]
                outcoords[timename] = outdimnames[0:2], outxtimes[:]
            
            if not multitime:
                outattrs["projection"] = outprojs[:, timeidx]
            else:
                outattrs["projection"] = outprojs
            
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

def combine_files(wrfseq, varname, timeidx, is_moving=None,
                  method="cat", squeeze=True, meta=True):
    
    # Handles generators, single files, lists, tuples, custom classes
    wrfseq = get_iterable(wrfseq)
    
    # Dictionary is unique
    if is_mapping(wrfseq):
        outarr = _combine_dict(wrfseq, varname, timeidx, method, meta)
    elif method.lower() == "cat":
        outarr = _cat_files(wrfseq, varname, timeidx, is_moving, 
                            squeeze, meta)
    elif method.lower() == "join":
        outarr = _join_files(wrfseq, varname, timeidx, is_moving, meta)
    else:
        raise ValueError("method must be 'cat' or 'join'")
    
    return outarr.squeeze() if squeeze else outarr


# Cache is a dictionary of already extracted variables
def _extract_var(wrfnc, varname, timeidx, is_moving, 
                 method, squeeze, cache, meta):
    # Mainly used internally so variables don't get extracted multiple times,
    # particularly to copy metadata.  This can be slow.
    if cache is not None:
        try:
            cache_var = cache[varname]
        except KeyError:
            pass
        else:
            if not meta:
                return npvalues(cache_var)
            
            return cache_var
    
    multitime = is_multi_time_req(timeidx)
    multifile = is_multi_file(wrfnc)
    
    if not multifile:
        if xarray_enabled() and meta:
            if is_moving is None:
                is_moving = is_moving_domain(wrfnc, varname)
            result = _build_data_array(wrfnc, varname, timeidx, is_moving)
        else:
            if not multitime:
                result = wrfnc.variables[varname][timeidx,:]
                result = result[np.newaxis, :] # So that no squeeze works
            else:
                result = wrfnc.variables[varname][:]
    else:
        # Squeeze handled in this routine, so just return it
        return combine_files(wrfnc, varname, timeidx, is_moving, 
                             method, squeeze, meta)
        
    return result.squeeze() if squeeze else result


def extract_vars(wrfnc, timeidx, varnames, method="cat", squeeze=True, 
                 cache=None, meta=True):
    if isstr(varnames):
        varlist = [varnames]
    else:
        varlist = varnames
    
    return {var:_extract_var(wrfnc, var, timeidx, None,
                             method, squeeze, cache, meta)
            for var in varlist}

# Python 3 compatability
def npbytes_to_str(var):
    return (bytes(c).decode("utf-8") for c in var[:])


def _make_time(timearr):
    return dt.datetime.strptime("".join(npbytes_to_str(timearr)), 
                                "%Y-%m-%d_%H:%M:%S")

def _file_times(wrfnc, timeidx):
    multitime = is_multi_time_req(timeidx)
    if multitime:
        times = wrfnc.variables["Times"][:,:]
        for i in py3range(times.shape[0]):
            yield _make_time(times[i,:])
    else:
        times = wrfnc.variables["Times"][timeidx,:]
        yield _make_time(times)
        

def extract_times(wrfnc, timeidx):
    multi_file = is_multi_file(wrfnc)
    if not multi_file:
        wrf_list = [wrfnc]
    else:
        wrf_list = wrfnc
    
    return [file_time 
            for wrf_file in wrf_list 
            for file_time in _file_times(wrf_file, timeidx)]        
        
    
def is_standard_wrf_var(wrfnc, var):
    multifile = is_multi_file(wrfnc)
    if multifile:
        if not is_mapping(wrfnc):
            wrfnc = next(iter(wrfnc))
        else:
            entry = wrfnc[next(iter(viewkeys(wrfnc)))]
            return is_standard_wrf_var(entry, var)
                
    return var in wrfnc.variables


def is_staggered(var, wrfnc):
    we = extract_dim(wrfnc, "west_east")
    sn = extract_dim(wrfnc, "south_north")
    bt = extract_dim(wrfnc, "bottom_top")
    
    if (var.shape[-1] != we or var.shape[-2] != sn or var.shape[-3] != bt):
        return True
    
    return False



def get_left_indexes(ref_var, expected_dims):
    """Returns the extra left side dimensions for a variable with an expected
    shape.
    
    For example, if a 2D variable contains an additional left side dimension
    for time, this will return the time dimension size.
    
    """
    extra_dim_num = ref_var.ndim - expected_dims
    
    if (extra_dim_num == 0):
        return []
    
    return tuple([ref_var.shape[x] for x in py3range(extra_dim_num)]) 

def iter_left_indexes(dims):
    """A generator which yields the iteration tuples for a sequence of 
    dimensions sizes.
    
    For example, if an array shape is (3,3), then this will yield:
    
    (0,0), (0,1), (1,0), (1,1)
    
    Arguments:
    
        - dims - a sequence of dimensions sizes (e.g. ndarry.shape)
    
    """
    arg = [py3range(dim) for dim in dims]
    for idxs in product(*arg):
        yield idxs
        
def get_right_slices(var, right_ndims, fixed_val=0):
    extra_dim_num = var.ndim - right_ndims
    if extra_dim_num == 0:
        return [slice(None)] * right_ndims
    
    return tuple([fixed_val]*extra_dim_num + 
                 [slice(None)]*right_ndims)

def get_proj_params(wrfnc, timeidx=0, varname=None):
    proj_params = extract_global_attrs(wrfnc, attrs=("MAP_PROJ", 
                                                "CEN_LAT", "CEN_LON",
                                                "TRUELAT1", "TRUELAT2",
                                                "MOAD_CEN_LAT", "STAND_LON", 
                                                "POLE_LAT", "POLE_LON"))
    multitime = is_multi_time_req(timeidx)
    if not multitime:
        time_idx_or_slice = timeidx
    else:
        time_idx_or_slice = slice(None)
    
    if varname is not None:
        if not is_coordvar(varname):
            coord_names = getattr(wrfnc.variables[varname], 
                                  "coordinates").split()
            lon_coord = coord_names[0]
            lat_coord = coord_names[1]
        else:
            lat_coord, lon_coord = get_coord_pairs(varname)
    else:
        lat_coord, lon_coord = latlon_coordvars(wrfnc.variables)
    
    return (wrfnc.variables[lat_coord][time_idx_or_slice,:],
            wrfnc.variables[lon_coord][time_idx_or_slice,:],
            proj_params)
    

class CoordPair(object):
    def __init__(self, x=None, y=None, i=None, j=None, lat=None, lon=None):
        self.x = x
        self.y = y
        self.i = i
        self.j = j
        self.lat = lat
        self.lon = lon
        
    def __repr__(self):
        args = []
        if self.x is not None:
            args.append("x={}".format(self.x))
            args.append("y={}".format(self.y))
            
        if self.i is not None:
            args.append("i={}".format(self.i))
            args.append("j={}".format(self.j))
        
        if self.lat is not None:
            args.append("lat={}".format(self.lat))
            args.append("lon={}".format(self.lon))
            
        argstr = ", ".join(args)
        
        return "{}({})".format(self.__class__.__name__, argstr)
    
    def __str__(self):
        return self.__repr__()
    
    def xy_str(self, fmt="{:.4f}, {:.4f}"):
        if self.x is None or self.y is None:
            return None
        
        return fmt.format(self.x, self.y)
    
    def latlon_str(self, fmt="{:.4f}, {:.4f}"):
        if self.lat is None or self.lon is None:
            return None
        
        return fmt.format(self.lat, self.lon)
    
    def ij_str(self, fmt="{:.4f}, {:.4f}"):
        if self.i is None or self.j is None:
            return None
        
        return fmt.format(self.i, self.j)
    

def from_args(func, argnames, *args, **kwargs):
    """Parses the function args and kargs looking for the desired argument 
    value. Otherwise, the value is taken from the default keyword argument 
    using the arg spec.
    
    """
    if isstr(argnames):
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

def _args_to_list2(func, args, kwargs):
    argspec = getargspec(func)
    
    # Build the full tuple with defaults filled in
    outargs = [None]*len(argspec.args)
    if argspec.defaults is not None:
        for i,default in enumerate(argspec.defaults[::-1], 1):
            outargs[-i] = default
    
    # Add the supplied args
    for i,arg in enumerate(args):
        outargs[i] = arg
    
    # Fill in the supplied kargs 
    for argname,val in viewitems(kwargs):
        argidx = argspec.args.index(argname)
        outargs[argidx] = val
        
    return outargs

def _args_to_list3(func, args, kwargs):
    sig = signature(func)
    bound = sig.bind(*args, **kwargs)
    bound.apply_defaults()
    
    return [x for x in bound.arguments.values()]
    

def args_to_list(func, args, kwargs):
    """Converts the mixed args/kwargs to a single list of args"""
    if version_info > (3,):
        _args_to_list = _args_to_list3
    else:
        _args_to_list = _args_to_list2
    
    return _args_to_list(func, args, kwargs)
    

def _arg_location2(func, argname, args, kwargs):
    argspec = getargspec(func)
        
    list_args = _args_to_list2(func, args, kwargs)
    
    # Return the new sequence and location
    if argname not in argspec.args and argname not in kwargs:
        return None
    
    result_idx = argspec.args.index(argname)
    
    return list_args, result_idx

def _arg_location3(func, argname, args, kwargs):
    sig = signature(func)
    params = list(sig.parameters.keys())
    
    list_args = _args_to_list3(func, args, kwargs)
    
    try:
        result_idx = params.index(argname) 
    except ValueError:
        return None
        
    return list_args, result_idx
    
    
def arg_location(func, argname, args, kwargs):
    """Parses the function args, kargs and signature looking for the desired 
    argument location (either in args, kargs, or argspec.defaults), 
    and returns a list containing representing all arguments in the 
    correct order with defaults filled in.
    
    """
    if version_info > (3,):
        _arg_location = _arg_location3
    else:
        _arg_location = _arg_location2
        
    return _arg_location(func, argname, args, kwargs)

def psafilepath():
    return os.path.join(os.path.dirname(__file__), "data", "psadilookup.dat")
    
    
    
        


        
    
    
    
    
    
    
    




        
    
    
    