from collections import Iterable, Mapping, OrderedDict
from itertools import product
import datetime as dt
import warnings

import numpy as n

from wrf.var.config import xarray_enabled
from wrf.var.projection import getproj

if xarray_enabled():
    from xarray import DataArray

__all__ = ["extract_vars", "extract_global_attrs", "extract_dim",
           "combine_files", "is_standard_wrf_var", "extract_times",
           "iter_left_indexes", "get_left_indexes", "get_right_slices",
           "is_staggered", "get_proj_params"]

def _is_multi_time(timeidx):
    if timeidx == -1:
        return True
    return False

def _is_multi_file(wrfnc):
    if isinstance(wrfnc, Iterable) and not isinstance(wrfnc, str):
        return True
    return False

def _is_mapping(wrfnc):
    if isinstance(wrfnc, Mapping):
        return True
    return False

def _is_moving_domain(wrfnc, latvar="XLAT", lonvar="XLONG"):
    lat1 = wrfnc.variables[latvar][0,:]
    lat2 = wrfnc.variables[latvar][-1,:]
    lon1 = wrfnc.variables[lonvar][0,:]
    lon2 = wrfnc.variables[lonvar][-1,:]
    
    if (lat1[0,0] != lat2[0,0] or lat1[-1,-1] != lat2[-1,-1] or 
        lon1[0,0] != lon2[0,0] or lon1[-1,-1] != lon2[-1,-1]):
        return True
    
    return False 

def _get_attr(wrfnc, attr):
    val = getattr(wrfnc, attr, None)
    
    # PyNIO puts single values in to an array
    if isinstance(val,n.ndarray):
        if len(val) == 1:
            return val[0] 
    return val
        
def extract_global_attrs(wrfnc, attrs):
    if isinstance(attrs, str):
        attrlist = [attrs]
    else:
        attrlist = attrs
        
    multifile = _is_multi_file(wrfnc)
    
    if multifile:
        if not _is_mapping(wrfnc):
            wrfnc = wrfnc[0]
        else:
            wrfnc = wrfnc[next(wrfnc.iterkeys())]
        
    return {attr:_get_attr(wrfnc, attr) for attr in attrlist}

def extract_dim(wrfnc, dim):
    if _is_multi_file(wrfnc):
        if not _is_mapping(wrfnc):
            wrfnc = wrfnc[0]
        else:
            wrfnc = wrfnc[next(wrfnc.iterkeys())]
    
    d = wrfnc.dimensions[dim]
    if not isinstance(d, int):
        return len(d) #netCDF4
    return d # PyNIO
        

# TODO
def _combine_dict(wrfseq, var, timeidx, method):
    """Dictionary combination creates a new left index for each key, then 
    does a cat or join for the list of files for that key"""
    
    multitime = _is_multi_time(timeidx)
    numfiles = len(wrfseq)  
    
    if not multitime:
        time_idx_or_slice = timeidx
    else:
        time_idx_or_slice = slice(None, None, None)
    
    keys = (list(x for x in xrange(numfiles)) if not _is_mapping(wrfseq) else 
            list(key for key in wrfseq.iterkeys()))
    
    # Check if 
    
    if not xarray_enabled():
        first_var = wrfseq[keys[0]].variables[var][time_idx_or_slice, :]
    else:
        first_var = _build_data_array(wrfseq[keys[0]], var, timeidx)
    
    # Create the output data numpy array based on the first array
    outdims = [numfiles]
    outdims += first_var.shape
    outdata = n.zeros(outdims, first_var.dtype)
    outdata[0,:] = first_var[:]
    
    for idx, key in enumerate(keys[1:], start=1):
        outdata[idx,:] = wrfseq[key].variables[var][time_idx_or_slice, :]
    
    if not xarray_enabled():
        outarr = outdata
    else:
        outname = str(first_var.name)
        outcoords = dict(first_var.coords)
        outdims = ["sequence"] + list(first_var.dims)
        outcoords["sequence"] = keys
        outattrs = dict(first_var.attrs)
        
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdims, attrs=outattrs)
        
    return outarr


def _cat_files(wrfseq, var, timeidx):
    file_times = extract_times(wrfseq, timeidx)
    
    multitime = _is_multi_time(timeidx) 

    time_idx_or_slice = timeidx if not multitime else slice(None, None, None)
    
    first_var = (_build_data_array(wrfseq[0], var, timeidx)
                 if xarray_enabled() else 
                 wrfseq[0].variables[var][time_idx_or_slice, :])
    
    outdims = [len(file_times)]
    
    # Making a new time dim, so ignore this one
    outdims += first_var.shape[1:]
        
    outdata = n.zeros(outdims, first_var.dtype)
    
    numtimes = first_var.shape[0]
    startidx = 0
    endidx = numtimes
    
    outdata[startidx:endidx, :] = first_var[:]
    
    startidx = endidx
    for wrfnc in wrfseq[1:]:
        vardata = wrfnc.variables[var][time_idx_or_slice, :]
        if multitime:
            numtimes = vardata.shape[0]
        else:
            numtimes = 1
            
        endidx = startidx + numtimes
        
        outdata[startidx:endidx, :] = vardata[:]
        
        startidx = endidx
        
    if xarray_enabled():
        # FIXME:  If it's a moving nest, then the coord arrays need to have same
        # time indexes as the whole data set
        outname = str(first_var.name)
        outattrs = dict(first_var.attrs)
        outcoords = dict(first_var.coords)
        outdimnames = list(first_var.dims)
        outcoords[outdimnames[0]] = file_times # New time dimension values
        
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdimnames, attrs=outattrs)
        
    else:
        outarr = outdata
        
    return outarr

def _join_files(wrfseq, var, timeidx):
        
    multitime = _is_multi_time(timeidx)
    numfiles = len(wrfseq)  
    
    if not multitime:
        time_idx_or_slice = timeidx
    else:
        time_idx_or_slice = slice(None, None, None)
    
    if xarray_enabled():
        first_var = _build_data_array(wrfseq[0], var, timeidx)
    else:
        first_var = wrfseq[0].variables[var][time_idx_or_slice, :]
    
    # Create the output data numpy array based on the first array
    outdims = [numfiles]
    outdims += first_var.shape
    outdata = n.zeros(outdims, first_var.dtype)
    
    outdata[0,:] = first_var[:]

    for idx, wrfnc in enumerate(wrfseq[1:], 1):
        print idx
        outdata[idx,:] = wrfnc.variables[var][time_idx_or_slice, :]
        
    if xarray_enabled():
        outname = str(first_var.name)
        outcoords = dict(first_var.coords)
        outattrs = dict(first_var.attrs)
        # New dimensions
        outdimnames = ["file_idx"] + list(first_var.dims)
        outcoords["file_idx"] = [i for i in xrange(numfiles)]
        
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdimnames, attrs=outattrs)
        
    else:
        outarr = outdata
        
    return outarr

def combine_files(wrfseq, var, timeidx, method="cat", squeeze=True):
    # Dictionary is unique
    if _is_mapping(wrfseq):
        outarr = _combine_dict(wrfseq, var, timeidx, method)
    
    if method.lower() == "cat":
        outarr = _cat_files(wrfseq, var, timeidx)
    elif method.lower() == "join":
        outarr = _join_files(wrfseq, var, timeidx)
    else:
        raise ValueError("method must be 'cat' or 'join'")
    
    if squeeze:
        return outarr.squeeze()
    
    return outarr

# Note, always returns the full data set with the time dimension included
def _build_data_array(wrfnc, varname, timeidx):
    multitime = _is_multi_time(timeidx)
    var = wrfnc.variables[varname]
    data = var[:]
    attrs = OrderedDict(var.__dict__)
    dimnames = var.dimensions
    
    # Add the coordinate variables here.
    coord_names = getattr(var, "coordinates").split()
    lon_coord = coord_names[0]
    lat_coord = coord_names[1]
    
    coords = {}
    lon_coord_var = wrfnc.variables[lon_coord]
    lat_coord_var = wrfnc.variables[lat_coord]
    
    if multitime:
        if _is_moving_domain(wrfnc, lat_coord, lon_coord):
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
                            **proj_params) for i in xrange(lats.shape[0])]
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
        
        
    coords[dimnames[0]] = extract_times(wrfnc, timeidx)
    
    attrs["projection"] = proj
    
    data_array = DataArray(data, name=varname, dims=dimnames, coords=coords,
                           attrs=attrs)
    
    return data_array

def _extract_var(wrfnc, varname, timeidx, method, squeeze):
    multitime = _is_multi_time(timeidx)
    multifile = _is_multi_file(wrfnc)
    
    if not multifile:
        if xarray_enabled():
            return _build_data_array(wrfnc, varname, timeidx)
        else:
            if not multitime:
                return wrfnc.variables[varname][timeidx,:]
            else:
                return wrfnc.variables[varname][:]
    else:
        return combine_files(wrfnc, varname, timeidx, method)

def extract_vars(wrfnc, timeidx, varnames, method="cat", squeeze=True):
    if isinstance(varnames, str):
        varlist = [varnames]
    else:
        varlist = varnames
    
    return {var:_extract_var(wrfnc, var, timeidx, method, squeeze)
            for var in varlist}

def _make_time(timearr):
    return dt.datetime.strptime("".join(timearr[:]), "%Y-%m-%d_%H:%M:%S")

def _file_times(wrfnc, timeidx):
    multitime = _is_multi_time(timeidx)
    if multitime:
        times = wrfnc.variables["Times"][:,:]
        for i in xrange(times.shape[0]):
            yield _make_time(times[i,:])
    else:
        times = wrfnc.variables["Times"][timeidx,:]
        yield _make_time(times)
        

def extract_times(wrfnc, timeidx):
    multi_file = _is_multi_file(wrfnc)
    if not multi_file:
        wrf_list = [wrfnc]
    else:
        wrf_list = wrfnc
    
    return [file_time 
            for wrf_file in wrf_list 
            for file_time in _file_times(wrf_file,timeidx) ]        
        
    
def is_standard_wrf_var(wrfnc, var):
    multifile = _is_multi_file(wrfnc)
    if multifile:
        wrfnc = wrfnc[0]
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
    
    return [ref_var.shape[x] for x in xrange(extra_dim_num)] 

def iter_left_indexes(dims):
    """A generator which yields the iteration tuples for a sequence of 
    dimensions sizes.
    
    For example, if an array shape is (3,3), then this will yield:
    
    (0,0), (0,1), (1,0), (1,1)
    
    Arguments:
    
        - dims - a sequence of dimensions sizes (e.g. ndarry.shape)
    
    """
    arg = [xrange(dim) for dim in dims]
    for idxs in product(*arg):
        yield idxs
        
def get_right_slices(var, right_ndims, fixed_val=0):
    extra_dim_num = var.ndim - right_ndims
    if extra_dim_num == 0:
        return [slice(None,None,None)] * right_ndims
    
    return [fixed_val]*extra_dim_num + [slice(None,None,None)]*right_ndims

def get_proj_params(wrfnc, timeidx=0, varname=None):
    proj_params = extract_global_attrs(wrfnc, attrs=("MAP_PROJ", 
                                                "CEN_LAT", "CEN_LON",
                                                "TRUELAT1", "TRUELAT2",
                                                "MOAD_CEN_LAT", "STAND_LON", 
                                                "POLE_LAT", "POLE_LON"))
    multitime = _is_multi_time(timeidx)
    if not multitime:
        time_idx_or_slice = timeidx
    else:
        time_idx_or_slice = slice(None, None, None)
    
    if varname is not None:
        coord_names = getattr(wrfnc.variables[varname], "coordinates").split()
        lon_coord = coord_names[0]
        lat_coord = coord_names[1]
    else:
        lat_coord = "XLAT"
        lon_coord = "XLONG"
    
    return (wrfnc.variables[lat_coord][time_idx_or_slice,:],
            wrfnc.variables[lon_coord][time_idx_or_slice,:],
            proj_params)
        
    
    
    
    
    
    
    




        
    
    
    