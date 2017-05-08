from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np
import numpy.ma as ma

from .constants import Constants
from .extension import (_interpz3d, _interp2dxy, _interp1d, _slp, _tk, _td, 
                        _rh, _uvmet, _smooth2d, _cape, _cloudfrac, _ctt, _dbz,
                        _srhel, _udhel, _avo, _pvo, _eth, _wetbulb, _tv, 
                        _omega, _pw)
from .decorators import convert_units
from .metadecorators import (set_alg_metadata, set_uvmet_alg_metadata, 
                             set_interp_metadata, set_cape_alg_metadata,
                             set_cloudfrac_alg_metadata,
                             set_smooth_metdata)
from .interputils import get_xy

@set_interp_metadata("xy")
def xy(field, pivot_point=None, angle=None, start_point=None, end_point=None,
       meta=True):
    """Return the x,y points for a line within a two-dimensional grid.
    
    This function is primarily used to obtain the x,y points when making a 
    cross section.
    
    Args:
    
        field (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            field with at least two dimensions.
                    
        pivot_point (:obj:`tuple` or :obj:`list`, optional): A 
            :obj:`tuple` or :obj:`list` with two entries, 
            in the form of [x, y] (or [west_east, south_north]), which 
            indicates the x,y location through which the plane will pass.  
            Must also specify `angle`.
        
        angle (:obj:`float`, optional): Only valid for cross sections where 
            a plane will be plotted through 
            a given point on the model domain. 0.0 represents a S-N cross 
            section.  90.0 is a W-E cross section. 
            
        start_point (:obj:`tuple` or :obj:`list`, optional): A 
            :obj:`tuple` or :obj:`list` with two entries, in the form of 
            [x, y] (or [west_east, south_north]), which indicates the start 
            x,y location through which the plane will pass.
            
        end_point (:obj:`tuple` or :obj:`list`, optional): A 
            :obj:`tuple` or :obj:`list` with two entries, in the form of 
            [x, y] (or [west_east, south_north]), which indicates the end x,y 
            location through which the plane will pass.
            
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: An array of 
        x,y points, which has shape num_points x 2.  
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    Examples:
    
        Example 1: Using Pivot Point and Angle
        
        .. code-block:: python
            
            from wrf import getvar, xy
            from netCDF4 import Dataset
        
            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
            field = wrf.getvar(wrfnc, "slp")
            
            # Use the center of the grid
            pivot = (field.shape[-1]/2.0, field.shape[-2]/2.0)
            
            # West-East
            angle = 90.0
            
            xy_points = xy(field, pivot_point=pivot, angle=angle)
            
        Example 2: Using Start Point and End Point
        
        .. code-block:: python
        
            from wrf import getvar, xy
            from netCDF4 import Dataset
        
            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
            field = wrf.getvar(wrfnc, "slp")
            
            # Make a diagonal of lower left to upper right
            start = (0, 0)
            end = (-1, -1)
            
            xy_points = xy(field, start_point=start, end_point=end)
            
    
    """
    return get_xy(field, pivot_point, angle, start_point, end_point)
    

@set_interp_metadata("1d")
def interp1d(field, z_in, z_out, missing=Constants.DEFAULT_FILL, 
             meta=True):
    """Return the linear interpolation of a one-dimensional variable.
    
    This function is typically used to interpolate a variable in a vertical
    column, but the coordinate system need not be a vertical coordinate 
    system.  Multiple interpolation points may be specified in the *z_out* 
    parameter.
    
    Args:
    
        field (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A
            one-dimensional field.  Metadata for *field* is only copied 
            to the output if *field* is a :class:`xarray.DataArray` object.
            
        z_in (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The 
            one-dimensional coordinates associated with *field* (usually the 
            vertical coordinates, either height or pressure).
            
        z_out (:class:`xarray.DataArray`, :class:`numpy.ndarray`): A 
            one-dimensional array of *z_in* coordinate points to interpolate 
            to.  Must be the same type as *z_in*. 
            
        missing (:obj:`float`, optional): The fill value to use for the 
            output.  Default is :data:`wrf.Constants.DEFAULT_FILL`.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
        
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: An array with the
        same dimensionality as *z_out* containing the interpolated values.
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
        
    Examples:
    
        Example 1:  Calculate the 850 hPa and 500 hPa values at location \
        x,y = (100,200)
        
        .. code-block:: python
            
            import numpy as np
            from wrf import getvar, interp1d
            from netCDF4 import Dataset
            
            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
            
            # Get a 1D vertical column for pressure at location x,y = 100,200
            p_1d = wrf.getvar(wrfnc, "pres", units="hPa")[:,200,100]
            
            # Get a 1D vertical column for height at location 100,200
            ht_1d = wrf.getvar(wrfnc, "z", units="dm")[:,200,100]
            
            # Want the heights (in decameters) at 850, 500 hPa
            levels = np.asarray([850., 500.])
            
            # Get the 850 hPa and 500 hPa values at location 100,200.
            interp_vals = interp1d(p_1d, ht_1d, levels)
    
    """
    return _interp1d(field, z_in, z_out, missing)


@set_interp_metadata("2dxy")
def interp2dxy(field3d, xy, meta=True):
    """Return a cross section for a three-dimensional field.
    
    The returned array will hold the vertical cross section data along the 
    line described by *xy*.
    
    This method differs from :meth:`wrf.vertcross` in that it will return 
    all vertical levels found in *field3d*.  :meth:`wrf.vertcross` includes 
    an additional interpolation to set the output to a fixed number of 
    vertical levels.  Also, a :class:`numpy.ma.MaskedArray` is not created 
    and this routine should be considered as low-level access to the underlying 
    Fortran routine.
    
    See Also:
    
        :meth:`wrf.xy`, :meth:`wrf.vertcross`
        
    Args:
    
        field3d (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The  
            array to interpolate with at least three dimensions, whose 
            rightmost dimensions are nz x ny x nx.
            
        xy (:class:`xarray.DataArray` or :class:`numpy.ndarray`): An array 
            of one less dimension than *field3d*, whose rightmost dimensions 
            are nxy x 2. This array holds the x,y pairs of a line across the 
            model domain. The requested vertical cross section will be 
            extracted from *field3d* along this line.
            
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
            
            
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: An array 
        containing the vertical cross section along the line *xy*.  The 
        returned dimensions will be the same as *xy*, but with the rightmost
        dimensions being nz x nxy. If xarray is enabled and the *meta* 
        parameter is True, then the result will be a :class:`xarray.DataArray` 
        object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
        
    Examples:
    
        Example 1:  Calculate the vertical cross section for RH for a diagonal 
        line from the lower left to the upper right of the domain.
        
        .. code-block:: python
            
            from wrf import getvar, xy, interp2dxy
            from netCDF4 import Dataset
            
            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
            
            rh = getvar(wrfnc, "rh")
            start = (0, 0)
            end = (-1, -1)
            xy_line = xy(rh, start_point=start, end_point=end)
            
            vert_cross = interp2dxy(rh, xy_line)
    
    """
    return _interp2dxy(field3d, xy)


@set_interp_metadata("horiz")
def interpz3d(field3d, vert, desiredlev, missing=Constants.DEFAULT_FILL,
              meta=True):
    """Return the field interpolated to a specified pressure or height level.
    
    This function is roughly equivalent to :meth:`interplevel`, but does not 
    handle multi-product diagnostics (uvmet, cape_3d, etc) that contain an 
    additional leftmost dimension for the product type.  Also, a 
    :class:`numpy.ma.MaskedArray` is not created and this routine should
    be considered as low-level access to the underlying Fortran routine.
    
    See Also:
    
        :meth:`wrf.interplevel`
        
    Args:
    
        field3d (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            three-dimensional field to interpolate, with the rightmost 
            dimensions of nz x ny x nx.
    
        vert (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            three-dimensional array for the vertical coordinate, typically 
            pressure or height.  This array must have the same dimensionality
            as *field3d*.
        
        desiredlev (:obj:`float`): The desired vertical level.  
            Must be in the same units as the *vert* parameter.
        
        missing (:obj:`float`): The fill value to use for the output.  
            Default is :data:`wrf.Constants.DEFAULT_FILL`.
        
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        interpolated variable.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    Example:
    
        Example 1: Interpolate Geopotential Height to 500 hPa
    
        .. code-block:: python
        
            from netCDF4 import Dataset
            from wrf import getvar, interpz3d
        
            wrfin = Dataset("wrfout_d02_2010-06-13_21:00:00")
            
            p = getvar(wrfin, "pressure")
            ht = getvar(wrfin, "z", units="dm")
            
            ht_500 = interpz3d(ht, p, 500.0)
    
    """
    return _interpz3d(field3d, vert, desiredlev, missing)


@set_alg_metadata(2, "pres", refvarndims=3, description="sea level pressure")
@convert_units("pressure", "hpa")
def slp(height, tkel, pres, qv, meta=True, units="hPa"):
    """Return the sea level pressure.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args:
    
        height (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Geopotential height in [m] with the rightmost dimensions being 
            bottom_top x south_north x west_east.
            
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *height*.  
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa] with 
            the same dimensionality as *height*.  
            
            Note:    
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as *height*.
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'slp'.  Default 
            is 'hPa'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        sea level pressure.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.temp`, :meth:`wrf.tk`
    
    """
    return _slp(height, tkel, pres, qv)


@set_alg_metadata(3, "pres", description="temperature")
@convert_units("temp", "k")
def tk(pres, theta, meta=True, units="K"):
    """Return the temperature.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa] with at least
            three dimensions.  The rightmost dimensions are bottom_top x 
            south_north x west_east.  
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        theta (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Potential
            temperature (perturbation plus reference temperature) in [K] with
            the same dimensionality as *pres*.
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'temp'.  Default 
            is 'K'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        temperature in the specified units.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.tk`
    
    """
    
    return _tk(pres, theta)


@set_alg_metadata(3, "pres", description="dew point temperature")
@convert_units("temp", "c")
def td(pres, qv, meta=True, units="degC"):
    """Return the dewpoint temperature.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa] with at least
            three dimensions.  The rightmost dimensions are bottom_top x 
            south_north x west_east.  
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as *pres*.
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'dp'.  Default 
            is 'degC'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        dewpoint temperature.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.rh`
    
    """
    return _td(pres, qv)


@set_alg_metadata(3, "pres", description="relative humidity", units=None)
def rh(qv, pres, tkel, meta=True):
    """Return the relative humidity.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with at least three dimensions. The 
            rightmost dimensions are bottom_top x south_north x west_east. 
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa] with the 
            same dimensionality as *qv*.  
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *qv*.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        relative humidity.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.td`
    
    """
    return _rh(qv, pres, tkel)


@set_uvmet_alg_metadata(latarg="lat", windarg="u")
@convert_units("wind", "m s-1")
def uvmet(u, v, lat, lon, cen_long, cone, meta=True, units="m s-1"):
    """Return the u,v components of the wind rotated to earth coordinates.
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain U 
        - return_val[1,...] will contain V
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        u (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The u 
            component of the wind [m s-1]. This variable can be staggered or 
            unstaggered, but must be at least two dimensions. If staggered,
            the rightmost dimensions are south_north x west east.
            
            If staggered, the rightmost dimensions are south_north x 
            west_east_stag.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used. 
                
        v (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The v 
            component of the wind [m s-1]. This variable can be staggered or 
            unstaggered, but must be at least two dimensions. If staggered,
            the rightmost dimensions are south_north x west east.
            
            If staggered, the rightmost dimensions are south_north_stag x 
            west_east.
            
        lat (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The 
            latitude array.  
            
            This array can either be:
            
                - two-dimensional of size south_north x west_east.
                - multi-dimensional with the same number of dimensions as *u* 
                  and *v*, but with rightmost dimensions south_north x 
                  west_east and the same leftmost dimensions as *u* and *v*
                - multi-dimensional with one fewer dimensions as *u* and *v*, 
                  with rightmost dimensions south_north x west_east and the same 
                  leftmost dimensions as *u* and *v*, minus the 
                  third-from-the-right dimension of *u* and *v*.
                
            Note:
            
                This variable must also be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used. 
                
        lon (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The 
            longitude array.  
            
            This array can either be:
            
                - two-dimensional of size south_north x west_east.
                - multi-dimensional with the same number of dimensions as *u* 
                  and *v*, but with rightmost dimensions south_north x 
                  west_east and the same leftmost dimensions as *u* and *v*
                - multi-dimensional with one fewer dimensions as *u* and *v*, 
                  with rightmost dimensions south_north x west_east and the same 
                  leftmost dimensions as *u* and *v*, minus the 
                  third-from-the-right dimension of *u* and *v*.
                
        
        cen_long (:obj:`float`): The standard longitude for the map projection.
        
        cone (:obj:`float`): The cone factor used for the map project. If the 
            projection is not a conic projection, the *cone* is simply 1.0. 
            For conic projections, the cone factor is given by:
            
            .. code-block:: python
            
                if((fabs(true_lat1 - true_lat2) > 0.1) and 
                (fabs(true_lat2 - 90.) > 0.1)): 
                    cone = (log(cos(true_lat1*radians_per_degree)) 
                    - log(cos(true_lat2*radians_per_degree)))
                    
                    cone = (cone / 
                    (log(tan((45.-fabs(true_lat1/2.))*radians_per_degree)) 
                    - log(tan((45.-fabs(true_lat2/2.))*radians_per_degree)))) 
                else:
                    cone = sin(fabs(true_lat1)*radians_per_degree)

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'uvmet'.  Default 
            is 'm s-1'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        u,v components of the wind rotated to earth coordinates.  The leftmost 
        dimension size is 2, for u and v. If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`
    
    """
    return _uvmet(u, v, lat, lon, cen_long, cone)


@set_smooth_metdata()
def smooth2d(field, passes, meta=True):
    """Return the field smoothed.
    
    This routine does not modify the original data.
    
    Args:
    
        field (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The field
            to smooth, which must be at least two dimensions.  Missing/fill 
            values will be ignored as long as the type is either a 
            :class:`numpy.ma.MaskedArray or a :class:`xarray.DataArray` with 
            a *_FillValue* attribute.
            
        passes (:obj:`int`): The number of smoothing passes.
        
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
        
    Returns:
    
        :class:`xarray.DataArray`, :class:`numpy.ma.MaskedArray` or \
        :class:`numpy.ndarray`): The smoothed field. If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be either a :class:`numpy.ndarray` or a :class:`numpy.ma.MaskedArray` 
        depending on the type for *field*.
        
    
    """
    return _smooth2d(field, passes)


@set_cape_alg_metadata(is2d=True, copyarg="pres_hpa")
def cape_2d(pres_hpa, tkel, qv, height, terrain, psfc_hpa, ter_follow, 
            missing=Constants.DEFAULT_FILL, meta=True):
    """Return the two-dimensional CAPE, CIN, LCL, and LFC.
    
    This function calculates the maximum convective available potential 
    energy (CAPE), maximum convective inhibition (CIN), 
    lifted condensation level (LCL), and level of free convection (LFC). This 
    function uses the RIP [Read/Interpolate/plot] code to calculate 
    potential energy (CAPE) and convective inhibition 
    (CIN) [J kg-1] only for the parcel with max theta-e 
    in the column (i.e. something akin to Colman's MCAPE). CAPE is defined as 
    the accumulated buoyant energy from the level of free convection (LFC) to 
    the equilibrium level (EL). CIN is defined as the accumulated negative 
    buoyant energy from the parcel starting point to the LFC. The word 'parcel' 
    here refers to a 500 meter deep parcel, with actual temperature and 
    moisture averaged over that depth. 
    
    The leftmost dimension of the returned array represents four different 
    quantities:
        
        - return_val[0,...] will contain CAPE [J kg-1]
        - return_val[1,...] will contain CIN [J kg-1]
        - return_val[2,...] will contain LCL [m]
        - return_val[3,...] will contain LFC [m]
    
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
            
        pres_hpa (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [hPa] with at 
            least three dimensions. The rightmost dimensions can be 
            top_bottom x south_north x west_east or bottom_top x south_north x
            west_east.
            
            Note:
            
                The units for *pres_hpa* are [hPa].
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *pres_hpa*.
            
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as *pres_hpa*.
            
        height (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Geopotential height in [m] with the same dimensionality as 
            *pres_hpa*.
            
        terrain (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Terrain height in [m].  This is at least a two-dimensional array 
            with the same dimensionality as *pres_hpa*, excluding the vertical 
            (bottom_top/top_bottom) dimension.
            
        psfc_hpa (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            The surface pressure in [hPa].  This is at least a two-dimensional 
            array with the same dimensionality as *pres_hpa*, excluding the 
            vertical (bottom_top/top_bottom) dimension.
            
            Note:
            
                The units for *psfc_hpa* are [hPa].
                
        ter_follow (:obj:`bool`): A boolean that should be set to True if the 
            data uses terrain following coordinates (WRF data).  Set to 
            False for pressure level data.
            
        missing (:obj:`float`, optional): The fill value to use for the 
            output.  Default is :data:`wrf.Constants.DEFAULT_FILL`.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        cape, cin, lcl, and lfc values as an array whose 
        leftmost dimension is 4 (0=CAPE, 1=CIN, 2=LCL, 3=LFC) .  If xarray is 
        enabled and the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.cape_3d`
    
    """
    
    if isinstance(ter_follow, bool):
        ter_follow = 1 if ter_follow else 0
    
    i3dflag = 0
    cape_cin = _cape(pres_hpa, tkel, qv, height, terrain, psfc_hpa, 
                     missing, i3dflag, ter_follow)
    
    left_dims = cape_cin.shape[1:-3]
    right_dims = cape_cin.shape[-2:]
    
    resdim = (4,) + left_dims + right_dims
    
    # Make a new output array for the result
    result = np.zeros(resdim, cape_cin.dtype)
    
    # Cape 2D output is not flipped in the vertical, so index from the 
    # end
    result[0,...,:,:] = cape_cin[0,...,-1,:,:]
    result[1,...,:,:] = cape_cin[1,...,-1,:,:]
    result[2,...,:,:] = cape_cin[1,...,-2,:,:]
    result[3,...,:,:] = cape_cin[1,...,-3,:,:]
    
    return ma.masked_values(result, missing)


@set_cape_alg_metadata(is2d=False, copyarg="pres_hpa")
def cape_3d(pres_hpa, tkel, qv, height, terrain, psfc_hpa, ter_follow, 
            missing=Constants.DEFAULT_FILL, meta=True):
    """Return the three-dimensional CAPE and CIN.
    
    This function calculates the maximum convective available potential 
    energy (CAPE) and maximum convective inhibition (CIN).  This 
    function uses the RIP [Read/Interpolate/plot] code to calculate 
    potential energy (CAPE) and convective inhibition 
    (CIN) [J kg-1] for every grid point in the entire 3D domain 
    (treating each grid point as a parcel).
    
    The leftmost dimension of the returned array represents two different 
    quantities:
        
        - return_val[0,...] will contain CAPE [J kg-1]
        - return_val[1,...] will contain CIN [J kg-1]
    
    This function also supports computing CAPE along a single vertical 
    column.  In this mode, the *pres_hpa*, *tkel*, *qv* and *height* arguments
    must be one-dimensional vertical columns, and the *terrain* and 
    *psfc_hpa* arguments must be scalar values 
    (:obj:`float`, :class:`numpy.float32` or :class:`numpy.float64`).
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
            
        pres_hpa (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [hPa] with at 
            least three dimensions when operating on a grid of values. The 
            rightmost dimensions can be top_bottom x south_north x west_east 
            or bottom_top x south_north x west_east.  
            When operating on only a single column of values, the vertical 
            column can be bottom_top or top_bottom.  In this case, *terrain* 
            and *psfc_hpa* must be scalars.
            
            Note:
            
                The units for *pres_hpa* are [hPa].
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *pres_hpa*.
            
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as *pres_hpa*.
            
        height (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Geopotential height in [m] with the same dimensionality as 
            *pres_hpa*.
            
        terrain (:class:`xarray.DataArray`, :class:`numpy.ndarray`, \
            or a scalar): Terrain height in [m].  When operating on a grid of 
            values, this argument is at least a two-dimensional array 
            with the same dimensionality as *pres_hpa*, excluding the vertical 
            (bottom_top/top_bottom) dimension.  When operating on a single 
            vertical column, this argument must be a scalar (:obj:`float`, 
            :class:`numpy.float32`, or :class:`numpy.float64`).
            
        psfc_hpa (:class:`xarray.DataArray`, :class:`numpy.ndarray`, \
            or a scalar): Surface pressure in [hPa].  When operating on a 
            grid of values, this argument is at least a two-dimensional array 
            with the same dimensionality as *pres_hpa*, excluding the vertical 
            (bottom_top/top_bottom) dimension.  When operating on a single 
            vertical column, this argument must be a scalar (:obj:`float`, 
            :class:`numpy.float32`, or :class:`numpy.float64`).
            
            Note:
            
                The units for *psfc_hpa* are [hPa].
                
        ter_follow (:obj:`bool`): A boolean that should be set to True if the 
            data uses terrain following coordinates (WRF data).  Set to 
            False for pressure level data.
            
        missing (:obj:`float`, optional): The fill value to use for the 
            output.  Default is :data:`wrf.Constants.DEFAULT_FILL`.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        CAPE and CIN as an array whose 
        leftmost dimension is 2 (0=CAPE, 1=CIN).  If xarray is 
        enabled and the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.cape_2d`
    
    """
    
    if isinstance(ter_follow, bool):
        ter_follow = 1 if ter_follow else 0
    
    i3dflag = 1
    cape_cin = _cape(pres_hpa, tkel, qv, height, terrain, psfc_hpa, 
                     missing, i3dflag, ter_follow)
    
    return ma.masked_values(cape_cin, missing)


@set_cloudfrac_alg_metadata(copyarg="pres")
def cloudfrac(pres, relh, meta=True):
    """Return the cloud fraction.
    
    The leftmost dimension of the returned array represents three different 
    quantities:
        
        - return_val[0,...] will contain LOW level cloud fraction
        - return_val[1,...] will contain MID level cloud fraction
        - return_val[2,...] will contain HIGH level cloud fraction
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args:  
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa], with the 
            rightmost dimensions as bottom_top x south_north x west_east
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        relh (:class:`xarray.DataArray` or :class:`numpy.ndarray`) Relative 
            humidity with the same dimensionality as *pres*

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        cloud fraction array whose leftmost dimension is 3 (LOW=0, MID=1, 
        HIGH=2).  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.rh`
    
    """
    return _cloudfrac(pres, relh)


@set_alg_metadata(2, "pres_hpa", refvarndims=3, 
                  description="cloud top temperature")
@convert_units("temp", "c")
def ctt(pres_hpa, tkel, qv, qcld, height, terrain, qice=None, meta=True,
        units="degC"):
    """Return the cloud top temperature.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args:  
            
        pres_hpa (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [hPa], with the 
            rightmost dimensions as bottom_top x south_north x west_east
            
            Note:
            
                The units for *psfc_hpa* are [hPa].
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
        
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *pres_hpa*.
            
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as *pres_hpa*.
            
        qcld (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Cloud water 
            vapor mixing ratio in [kg/kg] with the same dimensionality as 
            *pres_hpa*.
            
        height (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Geopotential height in [m] with the same dimensionality as 
            *pres_hpa*.
            
        terrain (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Terrain height in [m].  This is at least a two-dimensional array 
            with the same dimensionality as *pres_hpa*, excluding the vertical 
            (bottom_top/top_bottom) dimension.
            
        qice (:class:`xarray.DataArray` or :class:`numpy.ndarray`, optional): 
            Ice mixing ratio in [kg/kg] with the same dimensionality as 
            *pres_hpa*.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'ctt'.  Default 
            is 'degC'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        cloud top temperature.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.cloudfrac`
    
    """
    
    # Qice and QCLD need to be in g/kg
    if qice is None:
        qice = np.zeros(qv.shape, qv.dtype)
        haveqci = 0
    else:
        haveqci = 1 if qice.any() else 0
    
    return _ctt(pres_hpa, tkel, qice, qcld, qv, height, terrain, haveqci)
    


@set_alg_metadata(3, "pres", units="dBZ",
                  description="radar reflectivity")
def dbz(pres, tkel, qv, qr, qs=None, qg=None, use_varint=False, 
        use_liqskin=False, meta=True):
    """Return the simulated radar reflectivity.
    
    This function computes equivalent reflectivity factor [dBZ] at each 
    model grid point assuming spherical particles of constant density, 
    with exponential size distributions. This function is based on 
    "dbzcalc.f" in RIP. 
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args:  
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa], with the 
            rightmost dimensions as bottom_top x south_north x west_east
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
        
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *pres*.
            
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as *pres*.
            
        qr (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Rain water 
            vapor mixing ratio in [kg/kg] with the same dimensionality as 
            *pres*.
            
        qs (:class:`xarray.DataArray` or :class:`numpy.ndarray`, optional): 
            Snow mixing ratio in [kg/kg] with the same dimensionality as 
            *pres*.
            
        qg (:class:`xarray.DataArray` or :class:`numpy.ndarray`, optional): 
            Graupel mixing ratio in [kg/kg] with the same dimensionality as 
            *pres*.
            
        use_varint (:obj:`bool`, optional): When set to False, 
            the intercept parameters are assumed constant 
            (as in MM5's Reisner-2 bulk microphysical scheme). 
            When set to True, the variable intercept 
            parameters are used as in the more recent version of Reisner-2 
            (based on Thompson, Rasmussen, and Manning, 2004, Monthly weather 
            Review, Vol. 132, No. 2, pp. 519-542.).  
            
        use_liqskin (:obj:`bool`, optional): When set to True, frozen particles 
            that are at a temperature above freezing are assumed to scatter 
            as a liquid particle.  Set to False to disable.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        simulated radar reflectivity.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`
    
    """
    
    if qs is None:
        qs = np.zeros(qv.shape, qv.dtype)
        
    if qg is None:
        qg = np.zeros(qv.shape, qv.dtype)
    
    sn0 = 1 if qs.any() else 0
    ivarint = 1 if use_varint else 0
    iliqskin = 1 if use_liqskin else 0
    
    return _dbz(pres, tkel, qv, qr, qs, qg, sn0, ivarint, iliqskin)


@set_alg_metadata(2, "terrain", units="m2 s-2",
                  description="storm relative helicity")
def srhel(u, v, height, terrain, top=3000.0, meta=True):
    """Return the storm relative helicity.
    
    This function calculates storm relative helicity from WRF ARW output. 
    SRH (Storm Relative Helicity) is a measure of the potential for cyclonic 
    updraft rotation in right-moving supercells, and is calculated for the 
    lowest 1-km and 3-km layers above ground level. There is no clear threshold 
    value for SRH when forecasting supercells, since the formation of 
    supercells appears to be related more strongly to the deeper layer 
    vertical shear. Larger values of 0-3 km SRH (greater than 250 m2 s-2) 
    and 0-1 km SRH (greater than 100 m2 s-2), however, do suggest an 
    increased threat of tornadoes with supercells. For SRH, larger values are 
    generally better, but there are no clear "boundaries" between non-tornadic 
    and significant tornadic supercells.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        u (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The u 
            component of the wind that must have at least three dimensions.  
            The rightmost dimensions are bottom_top x south_north x west_east.
                
        v (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The v 
            component of the wind with the same dimensionality as *u*. 
            
        height (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Geopotential height in [m] with the same dimensionality as 
            *u*.
            
        terrain (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Terrain height in [m].  This is at least a two-dimensional array 
            with the same dimensionality as *u*, excluding the bottom_top 
            dimension.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used. 
                
        top (:obj:`float`):  The height of the layer below which helicity is 
            calculated (meters above ground level).

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        storm relative helicity.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.udhel`
        
    """
    
    # u, v get swapped in vertical
    _u = np.ascontiguousarray(u[...,::-1,:,:])
    _v = np.ascontiguousarray(v[...,::-1,:,:])
    _height = np.ascontiguousarray(height[...,::-1,:,:])
    
    return _srhel(_u, _v, _height, terrain, top)


@set_alg_metadata(2, "u", refvarndims=3, units="m2 s-2",
                  description="updraft helicity")
def udhel(zstag, mapfct, u, v, wstag, dx, dy, bottom=2000.0, top=5000.0, 
          meta=True):
    """Return the updraft helicity.
    
    This function calculates updraft helicity to detect 
    rotating updrafts. The formula follows Kain et al., 2008, Wea. and 
    Forecasting, 931-952, but this version has controls for the limits of 
    integration, *bottom* to *top*, in m AGL. Kain et al used 2000 to 5000 m. 
    The expected range is 25 to 250 m-2/s-2. Keith Brewster, CAPS/Univ. of 
    Oklahoma ; March, 2010
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        zstag (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Geopotential height in [m] that is at least three dimensions with
            a staggered vertical dimension.  The rightmost dimensions are 
            bottom_top_stag x south_north x west_east.
            
        mapfct (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The map 
            scale factor on the mass grid. An array of at least 
            two dimensions, whose rightmost two dimensions must be 
            south_north x west_east. If this array is more than two dimensions, 
            they must be the same as *zstag*'s leftmost dimensions.
    
        u (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The u 
            component of the wind [m s-1] whose rightmost three dimensions 
            must be bottom_top x south_north x west_east. The leftmost 
            dimensions must be the same as zp's leftmost dimensions.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used. 
                
        v (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The v 
            component of the wind [m s-1] whose rightmost three dimensions 
            must be bottom_top x south_north x west_east. The leftmost 
            dimensions must be the same as *zstag*'s leftmost dimensions. 
            
        wstag (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The z
            component of the wind [m s-1] with the same dimensionality as 
            *zstag*.
            
        dx (:obj:`float`): The distance between x grid points.
        
        dy (:obj:`float`): The distance between y grid points.
        
        bottom (:obj:`float`, optional): The bottom limit of integration. 
            Default is 2000.0.
            
        top (:obj:`float`, optional): The upper limit of integration. 
            Default is 5000.0.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        updraft helicity.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.srhel`
        
    """
    return _udhel(zstag, mapfct, u, v, wstag, dx, dy, bottom, top)


# Requires both u an v for dimnames
@set_alg_metadata(3, "ustag", units="10-5 s-1",
                  stagdim=-1, stagsubvar="vstag",
                  description="absolute vorticity")
def avo(ustag, vstag, msfu, msfv, msfm, cor, dx, dy, meta=True):
    """Return the absolute vorticity.
    
    This function returns absolute vorticity [10-5 s-1], which is the sum of 
    the relative vorticity at each grid point and the Coriolis parameter 
    at the latitude.

    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        ustag (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            The u component of the wind in [m s-1] that is at least three 
            dimensions with a staggered west_east dimension.  The rightmost 
            dimensions are bottom_top x south_north x west_east_stag.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used. 
            
        vstag (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            The v component of the wind in [m s-1] that is at least three 
            dimensions with a staggered south_north dimension.  The rightmost 
            dimensions are bottom_top x south_north_stag x west_east.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used. 
            
        msfu (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The map 
            scale factor on the u-grid that is at least 
            two dimensions, whose rightmost two dimensions must be 
            the same as *ustag*. If this array contains more than two 
            dimensions, they must be the same as *ustag* and *vstag*'s leftmost 
            dimensions.
            
        msfv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The map 
            scale factor on the v-grid that is at least 
            two dimensions, whose rightmost two dimensions must be 
            the same as *vstag*. If this array contains more than two 
            dimensions, they must be the same as *ustag* and *vstag*'s leftmost 
            dimensions.
            
        msfm (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The map 
            scale factor on the mass grid that is at least 
            two dimensions, whose rightmost two dimensions must be 
            south_north x west_east. If this array contains more than two 
            dimensions, they must be the same as *ustag* and *vstag*'s leftmost 
            dimensions.
            
        cor (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The Coriolis 
            sine latitude array that is at least 
            two dimensions, whose dimensions must be the same as *msfm*.
            
        dx (:obj:`float`): The distance between x grid points.
        
        dy (:obj:`float`): The distance between y grid points.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        absolute vorticity.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.pvo`
        
    """
    return _avo(ustag, vstag, msfu, msfv, msfm, cor, dx, dy)


@set_alg_metadata(3, "theta", units="PVU",
                  description="potential vorticity")
def pvo(ustag, vstag, theta, pres, msfu, msfv, msfm, cor, dx, dy, meta=True):
    """Return the potential vorticity.
    
    This function calculates the potential vorticity [PVU] at each grid point.

    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        ustag (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            The u component of the wind in [m s-1] that is at least three 
            dimensions with a staggered west_east dimension.  The rightmost 
            dimensions are bottom_top x south_north x west_east_stag.
             
            
        vstag (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            The v component of the wind in [m s-1] that is at least three 
            dimensions with a staggered south_north dimension.  The rightmost 
            dimensions are bottom_top x south_north_stag x west_east.
            
        theta (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The 
            potential temperature field [K] whose rightmost dimensions are 
            bottom_top x south_north x west_east and whose leftmost dimensions
            are the same as *ustag*.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa], with the 
            same dimensions as *theta*.
            
        msfu (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The map 
            scale factor on the u-grid that is at least 
            two dimensions, whose rightmost two dimensions must be 
            the same as *ustag*. If this array contains more than two 
            dimensions, they must be the same as *ustag* and *vstag*'s leftmost 
            dimensions.
            
        msfv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The map 
            scale factor on the v-grid that is at least 
            two dimensions, whose rightmost two dimensions must be 
            the same as *vstag*. If this array contains more than two 
            dimensions, they must be the same as *ustag* and *vstag*'s leftmost 
            dimensions.
            
        msfm (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The map 
            scale factor on the mass grid that is at least 
            two dimensions, whose rightmost two dimensions must be 
            south_north x west_east. If this array contains more than two 
            dimensions, they must be the same as *ustag* and *vstag*'s leftmost 
            dimensions.
            
        cor (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The Coriolis 
            sine latitude array that is at least 
            two dimensions, whose dimensions must be the same as *msfm*.
            
        dx (:obj:`float`): The distance between x grid points.
        
        dy (:obj:`float`): The distance between y grid points.

        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        potential vorticity.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.avo`
        
    """
    return _pvo(ustag, vstag, theta, pres, msfu, msfv, msfm, cor, dx, dy)


@set_alg_metadata(3, "qv",
                  description="equivalent potential temperature")
@convert_units("temp", "k")
def eth(qv, tkel, pres, meta=True, units="K"):
    """Return the equivalent potential temperature.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] that is at least three dimensions, with 
            the rightmost dimensions of bottom_top x south_north x west_east.
            
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *qv*.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa] with the 
            same dimensionality as *qv*. 
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'eth'.  Default 
            is 'K'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        equivalent potential temperature.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.temp`, :meth:`wrf.wetbulb`, 
        :meth:`tvirtual`
    
    """
    
    return _eth(qv, tkel, pres)


@set_alg_metadata(3, "pres",
                  description="wetbulb temperature")
@convert_units("temp", "k")
def wetbulb(pres, tkel, qv, meta=True, units="K"):
    """Return the wetbulb temperature.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa], with the 
            rightmost dimensions as bottom_top x south_north x west_east
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
                
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with same dimensionality as *pres*.
    
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as 
            *pres*
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'twb'.  Default 
            is 'K'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        wetbulb temperature.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.temp`, :meth:`wrf.eth`, 
        :meth:`tvirtual`
    
    """
    return _wetbulb(pres, tkel, qv)


@set_alg_metadata(3, "tkel",
                  description="virtual temperature")
@convert_units("temp", "k")
def tvirtual(tkel, qv, meta=True, units="K"):
    """Return the virtual temperature.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with the rightmost dimensions as bottom_top x south_north 
            x west_east.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as 
            *tkel*
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        units (:obj:`str`): The desired units.  Refer to the :meth:`getvar` 
            product table for a list of available units for 'tv'.  Default 
            is 'K'.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        virtual temperature.  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`wrf.temp`, :meth:`wrf.eth`, 
        :meth:`wetbulb`
    
    """
    return _tv(tkel, qv)


@set_alg_metadata(3, "qv", units="Pa s-1",
                  description="omega")
def omega(qv, tkel, w, pres, meta=True):
    """Return omega.
    
    This function calculates omega (dp/dt) [Pa s-1].
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
        
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the rightmost dimensions as 
            bottom_top x south_north x west_east.
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
            
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with the same dimensionality as *qv*.
            
        w (:class:`xarray.DataArray` or :class:`numpy.ndarray`): The vertical 
            velocity [m s-1] with the same dimensionality as *qv*.
            
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa] with the 
            same dimensionality as *qv*. 
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: Omega.  
        If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`, :meth:`uvmet`
    
    """
    return _omega(qv, tkel, w, pres)


@set_alg_metadata(2, "pres", refvarndims=3, units="kg m-2",
                  description="precipitable water")
def pw(pres, tkel, qv, height, meta=True):
    """Return the precipitable water.
    
    This is the raw computational algorithm and does not extract any variables 
    from WRF output files.  Use :meth:`wrf.getvar` to both extract and compute
    diagnostic variables.
    
    Args: 
    
        pres (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Full 
            pressure (perturbation + base state pressure) in [Pa], with the 
            rightmost dimensions as bottom_top x south_north x west_east
            
            Note:
            
                This variable must be 
                supplied as a :class:`xarray.DataArray` in order to copy the 
                dimension names to the output.  Otherwise, default names will
                be used.
                
        tkel (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Temperature
            in [K] with the same dimensionality as *pres*.
        
        qv (:class:`xarray.DataArray` or :class:`numpy.ndarray`): Water vapor 
            mixing ratio in [kg/kg] with the same dimensionality as *pres*
            
        height (:class:`xarray.DataArray` or :class:`numpy.ndarray`): 
            Geopotential height in [m] with the same dimensionality as 
            *pres*.
            
        meta (:obj:`bool`): Set to False to disable metadata and return 
            :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
    
    Warning:
    
        The input arrays must not contain any missing/fill values or 
        :data:`numpy.nan` values.
        
    Returns:
        
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The precipitable 
        water [kg m-2].  If xarray is enabled and 
        the *meta* parameter is True, then the result will be an 
        :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    See Also:
    
        :meth:`wrf.getvar`
    
    """
    tv = _tv(tkel, qv)
    
    return _pw(pres, tv, qv, height)

    
    
