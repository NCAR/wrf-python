from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from math import floor, ceil

import numpy as np

from .extension import interp2dxy

__all__ = ["to_positive_idxs", "calc_xy", "get_xy_z_params", "get_xy"]

def to_positive_idxs(shape, coord):
    if (coord[-2] >= 0 and coord[-1] >= 0):
        return coord
    
    return [x if (x >= 0) else shape[i]+x for (i,x) in enumerate(coord) ]

def calc_xy(xdim, ydim, pivot_point=None, angle=None, 
           start_point=None, end_point=None):
    """Returns the x,y points for the horizontal cross section line.
    
    xdim - maximum x-dimension
    ydim - maximum y-dimension
    pivot_point - a pivot point of (south_north, west_east) 
                  (must be used with angle)
    angle - the angle through the pivot point in degrees
    start_point - a start_point sequence of [south_north1, west_east1]
    end_point - an end point sequence of [south_north2, west_east2]
    
    """ 
    
    # Have a pivot point with an angle to find cross section
    if pivot_point is not None and angle is not None:
        xp = pivot_point[-1]
        yp = pivot_point[-2]
        
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
        x0 = start_point[-1]
        y0 = start_point[-2]
        x1 = end_point[-1]
        y1 = end_point[-2]
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
    
    for i in xrange(npts):
        xy[i,0] = x0 + i*dx
        xy[i,1] = y0 + i*dy
        
    return xy

def get_xy_z_params(z, pivot_point=None, angle=None,
                    start_point=None, end_point=None):
    
    xy = get_xy(z, pivot_point, angle, start_point, end_point)
    
    # Interp z
    var2dz = interp2dxy(z, xy)
    
    extra_dim_num = z.ndim - 3
    idx1 = tuple([0]*extra_dim_num + [0,0])
    idx2 = tuple([0]*extra_dim_num + [-1,0])
    
    #  interp to constant z grid
    if(var2dz[idx1] > var2dz[idx2]):  # monotonically decreasing coordinate
        z_max = floor(np.amax(z) / 10) * 10     # bottom value
        z_min = ceil(np.amin(z) / 10) * 10      # top value
        dz = 10
        nlevels = int((z_max-z_min) / dz)
        z_var2d = np.zeros((nlevels), dtype=z.dtype)
        z_var2d[0] = z_max
        dz = -dz
    else:
        z_max = np.amax(z)
        z_min = 0.
        dz = 0.01 * z_max
        nlevels = int(z_max / dz)
        z_var2d = np.zeros((nlevels), dtype=z.dtype)
        z_var2d[0] = z_min
    
    for i in xrange(1,nlevels):
        z_var2d[i] = z_var2d[0] + i*dz
        
    return xy, var2dz, z_var2d

def get_xy(var, pivot_point=None, angle=None, 
           start_point=None, end_point=None):
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
    
    xy = calc_xy(xdim, ydim, pos_pivot, angle, pos_start, pos_end)
    
    return xy
