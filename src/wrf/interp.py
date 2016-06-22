from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np
import numpy.ma as ma

# from .extension import (interpz3d, interp2dxy, interp1d,
#                         smooth2d, monotonic, vintrp, computevertcross,
#                         computeinterpline)

from .extension import (_interpz3d, _interp2dxy, _interp1d, _vertcross,
                        _interpline, _smooth2d, monotonic, vintrp)

from .metadecorators import set_interp_metadata
from .util import extract_vars, is_staggered
from .interputils import get_xy, get_xy_z_params
from .constants import Constants, ConversionFactors
from .terrain import get_terrain
from .geoht import get_height
from .temp import get_theta, get_temp, get_eth
from .pressure import get_pressure


#  Note:  Extension decorator is good enough to handle left dims
@set_interp_metadata("horiz")
def interplevel(field3d, z, desiredlev, missing=Constants.DEFAULT_FILL, 
                meta=True):
    """Interpolates a three-dimensional field specified pressure or 
    height level.

    Parameters
    ----------
    field3d : `xarray.DataArray` or `numpy.ndarray`
        A three-dimensional field.
    
    z : `xarray.DataArray` or `numpy.ndarray`
        A three-dimensional array for the vertical coordinate, typically 
        pressure or height.  
        
    desiredlev : float
        The desired vertical level.  Must be in the same units as the `z`
        parameter.
        
    missing : float
        The fill value to use for the output.  
        `Default is wrf.Constants.DEFAULT_FILL`.
        
    meta : {True, False}
        Set to False to disable metadata and return `numpy.ndarray` instead of 
        `xarray.DataArray`.  Default is True.
        
    Returns
    -------
    out : 'xarray.DataArray` or `numpy.ndarray`
        Returns the interpolated variable.  If xarray is enabled and 
        the meta parameter is True, then the result will be an 
        `xarray.DataArray` object.  Otherwise, the result will be a 
        `numpy.ndarray` object with no metadata.
    
    """
    r1 = _interpz3d(field3d, z, desiredlev, missing)
    masked_r1 = ma.masked_values (r1, missing)
    
    return masked_r1

@set_interp_metadata("cross")
def vertcross(field3d, z, missing=Constants.DEFAULT_FILL, 
              pivot_point=None, angle=None,
              start_point=None, end_point=None,
              include_latlon=False, 
              cache=None, meta=True):
    """Return the vertical cross section for a 3D field, interpolated 
    to a verical plane defined by a horizontal line.
    
    The horizontal line is defined by either including the 
    `pivot_point` and `angle` parameters, or the `start_point` and 
    `end_point` parameters.
    
    Parameters
    ----------
    field3d : `xarray.DataArray` or `numpy.ndarray`
        A three-dimensional field.
        
    z : `xarray.DataArray` or `numpy.ndarray`
        A three-dimensional array for the vertical coordinate, typically 
        pressure or height
    
    pivot_point : tuple or list, optional
        A tuple or list with two entries, in the form of [x, y] 
        (or [west_east, south_north]), which indicates the x,y location 
        through which the plane will pass.  Must also specify `angle`.
    
    angle : float, optional
        Only valid for cross sections where a plane will be plotted through 
        a given point on the model domain. 0.0 represents a S-N cross section 
        and 90.0 a W-E cross section. 
        
    start_point : tuple or list, optional
        A tuple or list with two entries, in the form of [x, y] 
        (or [west_east, south_north]), which indicates the start x,y location 
        through which the plane will pass.
        
    end_point : tuple or list, optional
        A tuple or list with two entries, in the form of [x, y] 
        (or [west_east, south_north]), which indicates the end x,y location 
        through which the plane will pass.
        
    include_latlon : {True, False}, optional
        Set to True to also interpolate the two-dimensional latitude and 
        longitude coordinates along the same horizontal line and include
        this information in the metadata (if enabled).  This can be 
        helpful for plotting.  Default is False.
        
    cache : dict, optional
        A dictionary of (varname, numpy.ndarray) pairs which can be used to 
        supply pre-extracted NetCDF variables to the computational routines.  
        This can be used to prevent the repeated variable extraction from large
        sequences of data files.  Default is None.
        
    meta : {True, False}, optional
        Set to False to disable metadata and return `numpy.ndarray` instead of 
        `xarray.DataArray`.  Default is True.
        
    Returns
    -------
    out : 'xarray.DataArray` or `numpy.ndarray`
        Returns the interpolated variable.  If xarray is enabled and 
        the meta parameter is True, then the result will be an 
        `xarray.DataArray` object.  Otherwise, the result will be a 
        `numpy.ndarray` object with no metadata.
        
    """
    
    try:
        xy = cache["xy"]
        var2dz = cache["var2dz"]
        z_var2d = cache["z_var2d"]
    except (KeyError, TypeError):
        xy, var2dz, z_var2d = get_xy_z_params(z, pivot_point, angle,
                                              start_point, end_point)
        
    res = _vertcross(field3d, xy, var2dz, z_var2d, missing)
    
    return ma.masked_values(res, missing)


@set_interp_metadata("line")
def interpline(field2d, pivot_point=None, 
               angle=None, start_point=None,
               end_point=None, include_latlon=False, 
               cache=None, meta=True):
    """Return the two-dimensional field interpolated along a line.
    
    Parameters
    ----------
    field2d : `xarray.DataArray` or `numpy.ndarray`
        A two-dimensional field.
    
    pivot_point : tuple or list, optional
        A tuple or list with two entries, in the form of [x, y] 
        (or [west_east, south_north]), which indicates the x,y location 
        through which the plane will pass.  Must also specify `angle`.
    
    angle : float, optional
        Only valid for cross sections where a plane will be plotted through 
        a given point on the model domain. 0.0 represents a S-N cross section 
        and 90.0 a W-E cross section. 
        
    start_point : tuple or list, optional
        A tuple or list with two entries, in the form of [x, y] 
        (or [west_east, south_north]), which indicates the start x,y location 
        through which the plane will pass.
        
    end_point : tuple or list, optional
        A tuple or list with two entries, in the form of [x, y] 
        (or [west_east, south_north]), which indicates the end x,y location 
        through which the plane will pass.
        
    include_latlon : {True, False}, optional
        Set to True to also interpolate the two-dimensional latitude and 
        longitude coordinates along the same horizontal line and include
        this information in the metadata (if enabled).  This can be 
        helpful for plotting.  Default is False.
        
    cache : dict, optional
        A dictionary of (varname, numpy.ndarray) pairs which can be used to 
        supply pre-extracted NetCDF variables to the computational routines.  
        This can be used to prevent the repeated variable extraction from large
        sequences of data files.  Default is None.
        
    meta : {True, False}, optional
        Set to False to disable metadata and return `numpy.ndarray` instead of 
        `xarray.DataArray`.  Default is True.
    
    
    Returns
    -------
    out : 'xarray.DataArray` or `numpy.ndarray`
        Returns the interpolated variable.  If xarray is enabled and 
        the meta parameter is True, then the result will be an 
        `xarray.DataArray` object.  Otherwise, the result will be a 
        `numpy.ndarray` object with no metadata.
        
    """
    
    try:
        xy = cache["xy"]
    except (KeyError, TypeError):
        xy = get_xy(field2d, pivot_point, angle, start_point, end_point)
        
    return _interpline(field2d, xy)


@set_interp_metadata("vinterp")
def vinterp(wrfnc, field, vert_coord, interp_levels, extrapolate=False, 
            field_type=None, log_p=False, timeidx=0, method="cat", 
            squeeze=True, cache=None, meta=True):
    """Return the field vertically interpolated to the given the type of 
    surface and a set of new levels.
    
    Parameters
    ----------
    wrfnc : `netCD4F.Dataset`, `Nio.NioFile`, or a sequence
        Input WRF ARW NetCDF data as a `netCDF4.Dataset`, `Nio.NioFile` or an 
        iterable sequence of the aforementioned types.
        
    field2d : `xarray.DataArray` or `numpy.ndarray`
        A two-dimensional field.
        
    vert_coord : {"pressure", "pres", "p", "ght_msl", 
                  "ght_agl", "theta", "th", "theta-e", 
                  "thetae", "eth"}
        A string indicating the vertical coordinate type to interpolate to.
        Valid strings are: 
            * "pressure", "pres", "p": pressure [hPa]
            * "ght_msl": grid point height msl [km]
            * "ght_agl": grid point height agl [km]
            * "theta", "th": potential temperature [K]
            * "theta-e", "thetae", "eth": equivalent potential temperature [K]
            
    interp_levels : sequence
        A 1D sequence of vertical levels to interpolate to.
        
    extrapolate : {True, False}, optional
        Set to True to extrapolate values below ground.  Default is False.
        
    field_type : {"none", "pressure", "pres", "p", "z", "tc", "tk", "theta", 
                  "th", "theta-e", "thetae", "eth", "ght"}, optional
        The type of field. Default is None.
        
    log_p : {True, False}
        Use the log of the pressure for interpolation instead of just pressure.
        Default is False.
        
    timeidx : int, optional
        The time index to use when extracting auxiallary variables used in 
        the interpolation.  This value must be set to match the same value
        used when the `field` variable was extracted.  Default is 0.
        
    method : {'cat', 'join'}, optional
        The aggregation method to use for sequences, either 'cat' or 'join'.  
        'cat' combines the data along the Time index. 'join' is creates a new 
        index for the file.  This must be set to the same method used when 
        extracting the `field` variable. The default is 'cat'.
        
    squeeze : {True, False}, optional
        Set to False to prevent dimensions with a size of 1 from being removed
        from the shape of the output.  Default is True.
        
    cache : dict, optional
        A dictionary of (varname, ndarray) which can be used to supply 
        pre-extracted NetCDF variables to the computational routines.  This can 
        be used to prevent the repeated variable extraction from large
        sequences of data files.  Default is None.
        
    meta : {True, False}, optional
        Set to False to disable metadata and return `numpy.ndarray` instead of 
        `xarray.DataArray`.  Default is True.
    
    
    Returns
    -------
    out : 'xarray.DataArray` or `numpy.ndarray`
        Returns the interpolated variable.  If xarray is enabled and 
        the meta parameter is True, then the result will be an 
        `xarray.DataArray` object.  Otherwise, the result will be a 
        `numpy.ndarray` object with no metadata.
        
    """
    
    # Remove case sensitivity
    field_type = field_type.lower() if field_type is not None else "none"
    vert_coord = vert_coord.lower() if vert_coord is not None else "none"
        
    valid_coords = ("pressure", "pres", "p", "ght_msl", 
                    "ght_agl", "theta", "th", "theta-e", "thetae", "eth")
    
    valid_field_types = ("none", "pressure", "pres", "p", "z",
                         "tc", "tk", "theta", "th", "theta-e", "thetae", 
                         "eth", "ght")
    
    icase_lookup = {"none" : 0,
                    "p" : 1,
                    "pres" : 1,
                    "pressure" : 1,
                    "z" : 2,
                    "ght" : 2,
                    "tc" : 3, 
                    "tk" : 4,
                    "theta" : 5,
                    "th" : 5,
                    "theta-e" : 6,
                    "thetae" : 6,
                    "eth" : 6}
    
    # These constants match what's in the fortran code.  
    rgas    = 287.04     #J/K/kg
    ussalr  = .0065      # deg C per m, avg lapse rate
    sclht   = rgas*256./9.81
    
    # interp_levels might be a list or tuple, make a numpy array
    if not isinstance(interp_levels, np.ndarray):
        interp_levels = np.asarray(interp_levels, np.float64)
        
    # TODO: Check if field is staggered
    if is_staggered(field, wrfnc):
        raise RuntimeError("Please unstagger field in the vertical")
    
    # Check for valid coord
    if vert_coord not in valid_coords:
        raise RuntimeError("'%s' is not a valid vertical "
                           "coordinate type" %  vert_coord)
    
    # Check for valid field type
    if field_type not in valid_field_types:
        raise RuntimeError("'%s' is not a valid field type" % field_type)
    
    log_p_int = 1 if log_p else 0
    
    icase = 0
    extrap = 0
    
    if extrapolate:
        extrap = 1
        icase = icase_lookup[field_type]
    
    # Extract vriables
    #timeidx = -1 # Should this be an argument?
    ncvars = extract_vars(wrfnc, timeidx, ("PSFC", "QVAPOR", "F"), 
                          method, squeeze, cache, meta=False)
    
    sfp = ncvars["PSFC"] * ConversionFactors.PA_TO_HPA
    qv = ncvars["QVAPOR"]
    coriolis = ncvars["F"]
    
    terht = get_terrain(wrfnc, timeidx, units="m", 
                        method=method, squeeze=squeeze, cache=cache)
    t = get_theta(wrfnc, timeidx,  units="k", 
                  method=method, squeeze=squeeze, cache=cache)
    tk = get_temp(wrfnc, timeidx, units="k",  
                  method=method, squeeze=squeeze, cache=cache)
    p = get_pressure(wrfnc, timeidx, units="pa",  
                     method=method, squeeze=squeeze, cache=cache)
    ght = get_height(wrfnc, timeidx, msl=True, units="m", 
                     method=method, squeeze=squeeze, cache=cache)
    ht_agl = get_height(wrfnc, timeidx, msl=False, units="m",
                        method=method, squeeze=squeeze, cache=cache)
    
    smsfp = _smooth2d(sfp, 3)        

    vcor = 0
    
    if vert_coord in ("pressure", "pres", "p"):
        vcor = 1
        vcord_array = p * ConversionFactors.PA_TO_HPA
        
    elif vert_coord == "ght_msl":
        vcor = 2
        vcord_array = np.exp(-ght/sclht)
        
    elif vert_coord == "ght_agl":
        vcor = 3
        vcord_array = np.exp(-ht_agl/sclht)
    
    elif vert_coord in ("theta", "th"):
        vcor = 4
        idir = 1
        icorsw = 0
        delta = 0.01
        
        p_hpa = p * ConversionFactors.PA_TO_HPA
        
        vcord_array = monotonic(t, p_hpa, coriolis, idir, delta, icorsw)
        
        # We only extrapolate temperature fields below ground 
        # if we are interpolating to pressure or height vertical surfaces.
        
        icase = 0 
        
    elif vert_coord in ("theta-e", "thetae", "eth"):
        vcor = 5
        icorsw = 0
        idir = 1
        delta = 0.01
        
        eth = get_eth(wrfnc, timeidx)
        
        p_hpa = p * ConversionFactors.PA_TO_HPA
        
        vcord_array = monotonic(eth, p_hpa, coriolis, idir, delta, icorsw)
        # We only extrapolate temperature fields below ground if we are
        # interpolating to pressure or height vertical surfaces
        icase = 0
    
    # Set the missing value
    if isinstance(field, ma.MaskedArray):
        missing = field.fill_value
    else:
        missing = Constants.DEFAULT_FILL
    
    if (field.shape != p.shape):
        raise ValueError("'field' shape does not match other variable shapes. "
                         "Verify that the 'timeidx' parameter matches the "
                         "same value used when extracting the 'field' "
                         "variable.")
            
    res = vintrp(field, p, tk, qv, ght, terht, sfp, smsfp,
                 vcord_array, interp_levels,
                 icase, extrap, vcor, log_p_int, missing)
    
    return ma.masked_values(res, missing)


    
    
    
    
