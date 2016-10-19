from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from math import floor, ceil

import numpy as np

from .extension import _interp2dxy
from .py3compat import py3range


def to_positive_idxs(shape, coord):
    """Return the positive index values.
    
    This function converts negative index values to positive index values.
    
    Args:
    
        shape (indexable sequence): The array shape.
        
        coord (indexable sequence): The coordinate pair for x and y.
        
    Returns:
        
        :obj:`list`: The coordinate values with all positive indexes.
        
    """
    if (coord[-2] >= 0 and coord[-1] >= 0):
        return coord
    
    return [x if (x >= 0) else shape[-i-1]+x for (i,x) in enumerate(coord)]


def _calc_xy(xdim, ydim, pivot_point=None, angle=None, 
           start_point=None, end_point=None):
    """Return the x,y points for the horizontal cross section line.
    
    Args:
    
        xdim (:obj:`int`): The x-dimension size.
        
        ydim (:obj:`int`): The y-dimension size.
        
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
            
    Returns:
    
        :class:`np.ndarray`: A two-dimensional array with the left index 
        representing each point along the line, and the rightmost dimension 
        having two values for the x and y coordinates [0=X, 1=Y].
    
    """
    # Have a pivot point with an angle to find cross section
    if pivot_point is not None and angle is not None:
        xp = pivot_point[-2]
        yp = pivot_point[-1]
        
        if xp >= xdim or yp >= ydim:
            raise ValueError("pivot point {} is outside of domain "
                             "with shape {}".format(pivot_point,
                                                         (xdim, ydim)))
        
        if (angle > 315.0 or angle < 45.0 
            or ((angle > 135.0) and (angle < 225.0))):
            
            #x = y*slope + intercept
            slope = -(360.-angle)/45.
            if( angle < 45. ):
                slope = angle/45.
            if( angle > 135.):
                slope = (angle-180.)/45.
            
            intercept = xp - yp*slope
            
            # find intersections with domain boundaries
            y0 = 0.
            x0 = y0*slope + intercept
            
            if( x0 < 0.):  # intersect outside of left boundary
                x0 = 0.
                y0 =  (x0 - intercept)/slope
            if( x0 > xdim-1):  #intersect outside of right boundary
                x0 = xdim-1
                y0 =  (x0 - intercept)/slope
            y1 = ydim-1.  #need to make sure this will be a float?
            x1 = y1*slope + intercept
            
            if( x1 < 0.):  # intersect outside of left boundary
                x1 = 0.
                y1 =  (x1 - intercept)/slope
            
            if( x1 > xdim-1):  # intersect outside of right boundary
                x1 = xdim-1
                y1 =  (x1 - intercept)/slope
        else:
            #  y = x*slope + intercept
            slope = (90.-angle)/45.
            if( angle > 225. ):
                slope = (270.-angle)/45.
            intercept = yp - xp*slope

            #find intersections with domain boundaries
            x0 = 0.
            y0 = x0*slope + intercept
            
            if( y0 < 0.):  # intersect outside of bottom boundary
                y0 = 0.
                x0 =  (y0 - intercept)/slope
            
            if( y0 > ydim-1):  # intersect outside of top boundary
                y0 = ydim-1
                x0 =  (y0 - intercept)/slope
            
            x1 = xdim-1.  #  need to make sure this will be a float?
            y1 = x1*slope + intercept
            
            if( y1 < 0.):  # intersect outside of bottom boundary
                y1 = 0.
                x1 =  (y1 - intercept)/slope
            
            if( y1 > ydim-1):# intersect outside of top boundary
                y1 = ydim-1
                x1 =  (y1 - intercept)/slope
    elif start_point is not None and end_point is not None:
        x0 = start_point[-2]
        y0 = start_point[-1]
        x1 = end_point[-2]
        y1 = end_point[-1]
        
        if x0 >= xdim or y0 >= ydim:
            raise ValueError("start_point {} is outside of domain "
                             "with shape {}".format(start_point, (xdim, ydim)))
            
        if x1 >= xdim or y1 >= ydim:
            raise ValueError("end_point {} is outside of domain "
                             "with shape {}".format(end_point, (xdim, ydim)))
        
        if ( x1 > xdim-1 ): 
            x1 = xdim
        if ( y1 > ydim-1): 
            y1 = ydim
    else:
        raise ValueError("invalid start/end or pivot/angle arguments")
    
    dx = x1 - x0
    dy = y1 - y0
    distance = (dx*dx + dy*dy)**0.5
    npts = int(distance)
    dxy = distance/npts
    
    xy = np.zeros((npts,2), "float")

    dx = dx/npts
    dy = dy/npts
    
    for i in py3range(npts):
        xy[i,0] = x0 + i*dx
        xy[i,1] = y0 + i*dy
        
    return xy


def get_xy_z_params(z, pivot_point=None, angle=None,
                    start_point=None, end_point=None,
                    levels=None):
    """Return the cross section parameters.
    
    This function returns the xy horizontal cross section line coordinates, 
    the xy x z vertical values interpolated along the xy cross section 
    line, and the fixed vertical levels to be used by the cross section
    algorithm (at ~1% increments for the minimum to maximum vertical 
    span).
    
    Args:
    
        z (:class:`numpy.ndarray`): The vertical coordinate, whose rightmost 
            dimensions are bottom_top x south_north x west_east.
            
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
            
        levels (sequence): A sequence of :obj:`float` for the desired
            vertical levels in the output array.  If None, a fixed set of 
            vertical levels is provided.  Default is None.
            
    Returns:
    
        :obj:`tuple`:  A tuple containing the xy horizontal cross section  
        coordinates, the vertical values interpolated along the xy cross 
        section line, and the fixed vertical levels used by the 
        cross section algorithm at ~1% increments of minimum to maximum 
        vertical span.
    
    """
    
    xy = get_xy(z, pivot_point, angle, start_point, end_point)
    
    # Interp z
    var2dz = _interp2dxy(z, xy)
    
    extra_dim_num = z.ndim - 3
    idx1 = tuple([0]*extra_dim_num + [0,0])
    idx2 = tuple([0]*extra_dim_num + [-1,0])
    
    if levels is None:
        #  interp to constant z grid
        if(var2dz[idx1] > var2dz[idx2]):  # monotonically decreasing coordinate
            z_max = floor(np.amax(z)/10) * 10     # bottom value
            z_min = ceil(np.amin(z)/10) * 10      # top value
            dz = 10
            nlevels = int((z_max - z_min)/dz)
            z_var2d = np.zeros((nlevels), dtype=z.dtype)
            z_var2d[0] = z_max
            dz = -dz
        else:
            z_max = np.amax(z)
            z_min = 0.
            dz = 0.01*z_max
            nlevels = int(z_max/dz)
            z_var2d = np.zeros((nlevels), dtype=z.dtype)
            z_var2d[0] = z_min
    else:
        z_var2d = np.asarray(levels)
    
    for i in py3range(1,nlevels):
        z_var2d[i] = z_var2d[0] + i*dz
        
    return xy, var2dz, z_var2d


def get_xy(var, pivot_point=None, angle=None, 
           start_point=None, end_point=None):
    """Return the x,y points for the horizontal cross section line.
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A variable
            that contains a :attr:`shape` attribute.
        
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
            
    Returns:
    
        :class:`np.ndarray`: A two-dimensional array with the left index 
        representing each point along the line, and the rightmost dimension 
        having two values for the x and y coordinates [0=X, 1=Y].
    
    """
    if pivot_point is not None:
        pos_pivot = to_positive_idxs(var.shape[-2:], pivot_point)
    else:
        pos_pivot = pivot_point
        
    if start_point is not None:
        pos_start = to_positive_idxs(var.shape[-2:], start_point)
    else:
        pos_start = start_point
    
    if end_point is not None:
        pos_end = to_positive_idxs(var.shape[-2:], end_point)
    else:
        pos_end = start_point
        
    xdim = var.shape[-1]
    ydim = var.shape[-2]
    
    xy = _calc_xy(xdim, ydim, pos_pivot, angle, pos_start, pos_end)
    
    return xy
