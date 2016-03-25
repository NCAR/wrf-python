from collections import Iterable, Mapping, OrderedDict
from itertools import product
import datetime as dt

import numpy as np

from .config import xarray_enabled
from .projection import getproj

if xarray_enabled():
    from xarray import DataArray

__all__ = ["extract_vars", "extract_global_attrs", "extract_dim",
           "combine_files", "is_standard_wrf_var", "extract_times",
           "iter_left_indexes", "get_left_indexes", "get_right_slices",
           "is_staggered", "get_proj_params", "viewitems", "viewkeys",
           "viewvalues"]

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

def _corners_moved(wrfnc, first_ll_corner, first_ur_corner, latvar, lonvar):
    lats = wrfnc.variables[latvar]
    lons = wrfnc.variables[lonvar]
    # Need to check all times
    for i in xrange(lats.shape[-3]):
        start_idxs = [0]*lats.ndim
        start_idxs[-3] = i
        
        end_idxs = [-1]*lats.ndim
        end_idxs[-3] = i
        
        if (first_ll_corner[0] != lats[start_idxs] or 
            first_ll_corner[1] != lons[start_idxs] or 
            first_ur_corner[0] != lats[end_idxs] or 
            first_ur_corner[1] != lons[end_idxs]):
            return True
    
    return False

def _is_moving_domain(wrfseq, varname=None, latvar="XLAT", lonvar="XLONG"):
    # In case it's just a single file
    if not _is_multi_file(wrfseq):
        wrfseq = [wrfseq]
        
    # Slow, but safe. Compare the corner points to the first item and see
    # any move.  User iterator protocol in case wrfseq is not a list/tuple.
    wrf_iter = iter(wrfseq)
    
    first_wrfnc = next(wrf_iter)
    
    if varname is not None:
        coord_names = getattr(first_wrfnc.variables[varname], 
                              "coordinates").split()
        lon_coord = coord_names[0]
        lat_coord = coord_names[1]
    else:
        lon_coord = lonvar
        lat_coord = latvar
    
    lats = first_wrfnc.variables[lat_coord]
    lons = first_wrfnc.variables[lon_coord]
    
    zero_idxs = [0] * first_wrfnc.variables[lat_coord].ndim
    last_idxs = list(zero_idxs)
    last_idxs[-2:] = [-1]*2
    
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
    if isinstance(val,np.ndarray):
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
        
    return {attr:_get_global_attr(wrfnc, attr) for attr in attrlist}

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
        
def _combine_dict(wrfdict, varname, timeidx, method):
    """Dictionary combination creates a new left index for each key, then 
    does a cat or join for the list of files for that key"""
    keynames = []
    numkeys = len(wrfdict)  
    
    key_iter = iter(viewkeys(wrfdict))
    first_key = next(key_iter)
    keynames.append(first_key)
    
    first_array = _extract_var(wrfdict[first_key], varname, 
                              timeidx, method, squeeze=False)
    
    
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
                                   method, squeeze=False)
            
            if outdata.shape[1:] != vardata.shape:
                raise ValueError("data sequences must have the "
                                   "same size for all dictionary keys")
            outdata[idx,:] = vardata.values[:]
            idx += 1
    
    if not xarray_enabled():
        outarr = outdata
    else:
        outname = str(first_array.name)
        # Note: assumes that all entries in dict have same coords
        outcoords = OrderedDict(first_array.coords)
        outdims = ["key"] + list(first_array.dims)
        outcoords["key"] = keynames
        outattrs = OrderedDict(first_array.attrs)
        
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdims, attrs=outattrs)
        
    return outarr

# TODO:  implement in C
def _cat_files(wrfseq, varname, timeidx):
    is_moving = _is_moving_domain(wrfseq, varname)
    
    file_times = extract_times(wrfseq, timeidx)
    
    multitime = _is_multi_time(timeidx) 

    time_idx_or_slice = timeidx if not multitime else slice(None, None, None)
    
    # wrfseq might be a generator
    wrf_iter = iter(wrfseq)
    
    first_var = (_build_data_array(next(wrf_iter), varname, timeidx, is_moving)
                 if xarray_enabled() else 
                 wrfseq[0].variables[varname][time_idx_or_slice, :])
    
    outdims = [len(file_times)]
    
    # Making a new time dim, so ignore this one
    outdims += first_var.shape[1:]
        
    outdata = np.empty(outdims, first_var.dtype)
    
    numtimes = first_var.shape[0]
    startidx = 0
    endidx = numtimes
    
    outdata[startidx:endidx, :] = first_var[:]
    
    startidx = endidx
    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            vardata = wrfnc.variables[varname][time_idx_or_slice, :]
            
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
        outattrs = OrderedDict(first_var.attrs)
        outcoords = OrderedDict(first_var.coords)
        outdimnames = list(first_var.dims)
        outcoords[outdimnames[0]] = file_times # New time dimension values
        
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdimnames, attrs=outattrs)
        
    else:
        outarr = outdata
        
    return outarr

# TODO:  implement in C
def _join_files(wrfseq, varname, timeidx):
    is_moving = _is_moving_domain(wrfseq, varname)
    multitime = _is_multi_time(timeidx)
    numfiles = len(wrfseq)  
    
    if not multitime:
        time_idx_or_slice = timeidx
    else:
        time_idx_or_slice = slice(None, None, None)
    
    # wrfseq might be a generator
    wrf_iter = iter(wrfseq)
    
    if xarray_enabled():
        first_var = _build_data_array(next(wrf_iter), varname, 
                                      timeidx, is_moving)
    else:
        first_var = (next(wrf_iter)).variables[varname][time_idx_or_slice, :]
    
    # Create the output data numpy array based on the first array
    outdims = [numfiles]
    outdims += first_var.shape
    outdata = np.empty(outdims, first_var.dtype)
    
    outdata[0,:] = first_var[:]
    
    idx=1
    while True:
        try:
            wrfnc = next(wrf_iter)
        except StopIteration:
            break
        else:
            outdata[idx,:] = wrfnc.variables[varname][time_idx_or_slice, :]
            idx += 1
        
    if xarray_enabled():
        outname = str(first_var.name)
        outcoords = OrderedDict(first_var.coords)
        outattrs = OrderedDict(first_var.attrs)
        # New dimensions
        outdimnames = ["file_idx"] + list(first_var.dims)
        outcoords["file_idx"] = [i for i in xrange(numfiles)]
        
        outarr = DataArray(outdata, name=outname, coords=outcoords, 
                           dims=outdimnames, attrs=outattrs)
        
    else:
        outarr = outdata
        
    return outarr

def combine_files(wrfseq, varname, timeidx, method="cat", squeeze=True):
    # Dictionary is unique
    if _is_mapping(wrfseq):
        outarr = _combine_dict(wrfseq, varname, timeidx, method)
    elif method.lower() == "cat":
        outarr = _cat_files(wrfseq, varname, timeidx)
    elif method.lower() == "join":
        outarr = _join_files(wrfseq, varname, timeidx)
    else:
        raise ValueError("method must be 'cat' or 'join'")
    
    return outarr.squeeze() if squeeze else outarr

# Note, always returns the full data set with the time dimension included
def _build_data_array(wrfnc, varname, timeidx, is_moving_domain):
    multitime = _is_multi_time(timeidx)
    var = wrfnc.variables[varname]
    data = var[:]
    attrs = OrderedDict(var.__dict__)
    dimnames = var.dimensions
    
    # WRF variables will have a coordinates attribute.  MET_EM files have 
    # a stagger attribute which indicates the coordinate variable.
    try:
        # WRF files
        coord_attr = getattr(var, "coordinates")
    except KeyError:
        try:
            # met_em files
            stag_attr = getattr(var, "stagger")
        except KeyError:
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
    
    coords = OrderedDict()
    
    # Handle lat/lon coordinates and projection information if available
    if lon_coord is not None and lat_coord is not None:
        lon_coord_var = wrfnc.variables[lon_coord]
        lat_coord_var = wrfnc.variables[lat_coord]
    
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
        
        attrs["projection"] = proj
        
    
    if dimnames[0] == "Time":
        coords[dimnames[0]] = extract_times(wrfnc, timeidx)
    
    data_array = DataArray(data, name=varname, dims=dimnames, coords=coords,
                           attrs=attrs)
    
    return data_array

# Cache is a dictionary of already extracted variables
def _extract_var(wrfnc, varname, timeidx, method, squeeze, cache):
    # Mainly used internally so variables don't get extracted multiple times,
    # particularly to copy metadata.  This can be slow.
    if cache is not None:
        try:
            return cache[varname]
        except KeyError:
            pass
    
    is_moving = _is_moving_domain(wrfnc, varname)
    multitime = _is_multi_time(timeidx)
    multifile = _is_multi_file(wrfnc)
    
    if not multifile:
        if xarray_enabled():
            result = _build_data_array(wrfnc, varname, timeidx, is_moving)
        else:
            if not multitime:
                result = wrfnc.variables[varname][timeidx,:]
            else:
                result = wrfnc.variables[varname][:]
    else:
        # Squeeze handled in this routine, so just return it
        return combine_files(wrfnc, varname, timeidx, method, squeeze)
        
    return result.squeeze() if squeeze else result

def extract_vars(wrfnc, timeidx, varnames, method="cat", squeeze=True, 
                 cache=None):
    if isinstance(varnames, str):
        varlist = [varnames]
    else:
        varlist = varnames
    
    return {var:_extract_var(wrfnc, var, timeidx, method, squeeze, cache)
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
    
# Dictionary python 2-3 compatibility stuff
def viewitems(d):
    func = getattr(d, "viewitems", None)
    if func is None:
        func = d.items
    return func()


def viewkeys(d):
    func = getattr(d, "viewkeys", None)
    if func is None:
        func = d.keys
    return func()


def viewvalues(d):
    func = getattr(d, "viewvalues", None)
    if func is None:
        func = d.values
    return func()

# Helper utilities for metadata
class either(object):
    def __init__(self, *varnames):
        self.varnames = varnames
    
    def __call__(self, wrfnc):
        if _is_multi_file(wrfnc):
            wrfnc = next(iter(wrfnc))
            
        for varname in self.varnames:
            if varname in wrfnc:
                return varname
        
        raise ValueError("{} are not valid variable names".format(
                                                            self.varnames))

class combine_with:
    # Remove remove_idx first, then insert_idx is applied to removed set
    def __init__(self, varname, remove_dims=None, insert_before=None, 
                 new_dimnames=None, new_coords=None):
        self.varname = varname
        self.remove_dims = remove_dims
        self.insert_before = insert_before
        self.new_dimnames = new_dimnames
        self.new_coords = new_coords
    
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
        
        


        
    
    
    
    
    
    
    




        
    
    
    