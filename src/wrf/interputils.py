from __future__ import (absolute_import, division, print_function)

from math import floor, ceil

import numpy as np

from .extension import _interp2dxy
from .py3compat import py3range
from .coordpair import CoordPair
from .constants import Constants, ProjectionTypes
from .latlonutils import _ll_to_xy
from .util import pairs_to_latlon


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

    return [x if (x >= 0) else shape[-i-1]+x for (i, x) in enumerate(coord)]


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

            slope = -(360.-angle)/45.
            if(angle < 45.):
                slope = angle/45.
            if(angle > 135.):
                slope = (angle-180.)/45.

            intercept = xp - yp*slope

            # find intersections with domain boundaries
            y0 = 0.
            x0 = y0*slope + intercept

            if(x0 < 0.):  # intersect outside of left boundary
                x0 = 0.
                y0 = (x0 - intercept)/slope
            if(x0 > xdim-1):  # intersect outside of right boundary
                x0 = xdim-1
                y0 = (x0 - intercept)/slope
            y1 = ydim-1.  # need to make sure this will be a float?
            x1 = y1*slope + intercept

            if(x1 < 0.):  # intersect outside of left boundary
                x1 = 0.
                y1 = (x1 - intercept)/slope

            if(x1 > xdim-1):  # intersect outside of right boundary
                x1 = xdim-1
                y1 = (x1 - intercept)/slope
        else:
            #  y = x*slope + intercept
            slope = (90.-angle)/45.
            if (angle > 225.):
                slope = (270.-angle)/45.
            intercept = yp - xp*slope

            # Find intersections with domain boundaries
            x0 = 0.
            y0 = x0*slope + intercept

            if (y0 < 0.):  # intersect outside of bottom boundary
                y0 = 0.
                x0 = (y0 - intercept)/slope

            if (y0 > ydim-1):  # intersect outside of top boundary
                y0 = ydim-1
                x0 = (y0 - intercept)/slope

            x1 = xdim-1.  # need to make sure this will be a float?
            y1 = x1*slope + intercept

            if (y1 < 0.):  # intersect outside of bottom boundary
                y1 = 0.
                x1 = (y1 - intercept)/slope

            if (y1 > ydim-1):  # intersect outside of top boundary
                y1 = ydim-1
                x1 = (y1 - intercept)/slope
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
    else:
        raise ValueError("invalid start/end or pivot/angle arguments")

    dx = x1 - x0
    dy = y1 - y0
    distance = (dx*dx + dy*dy)**0.5
    npts = int(distance) + 1

    xy = np.zeros((npts, 2), "float")

    dx = dx/(npts-1)
    dy = dy/(npts-1)

    for i in py3range(npts):
        xy[i, 0] = x0 + i*dx
        xy[i, 1] = y0 + i*dy

    return xy


def get_xy_z_params(z, pivot_point=None, angle=None,
                    start_point=None, end_point=None,
                    levels=None, autolevels=100):
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

        autolevels(:obj:`int`, optional): The number of evenly spaced
            automatically chosen vertical levels to use when *levels*
            is None. Default is 100.

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
    idx1 = tuple([0]*extra_dim_num + [0, 0])
    idx2 = tuple([0]*extra_dim_num + [-1, 0])

    if levels is None:
        #  interp to constant z grid
        if(var2dz[idx1] > var2dz[idx2]):  # monotonically decreasing coordinate
            z_max = floor(np.amax(z)/10) * 10     # bottom value
            z_min = ceil(np.amin(z)/10) * 10      # top value
            dz = (1.0/autolevels) * (z_max - z_min)
            z_var2d = np.zeros((autolevels), dtype=z.dtype)
            z_var2d[0] = z_max
            dz = -dz
        else:
            z_max = np.amax(z)
            z_min = 0.
            dz = (1.0/autolevels)*z_max
            z_var2d = np.zeros((autolevels), dtype=z.dtype)
            z_var2d[0] = z_min

        for i in py3range(1, autolevels):
            z_var2d[i] = z_var2d[0] + i*dz
    else:
        z_var2d = np.asarray(levels, z.dtype)

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


def to_xy_coords(pairs, wrfin=None, timeidx=0, stagger=None, projection=None,
                 ll_point=None):
    """Return the coordinate pairs in grid space.

    This function converts latitude,longitude coordinate pairs to
    x,y coordinate pairs.

    Args:

        pairs (:class:`CoordPair` or sequence): A single coordinate pair or
            a sequence of coordinate pairs to be converted.

        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable, optional): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types. This is used
            to obtain the map projection when using latitude,longitude
            coordinates. Should not be used when working with x,y
            coordinates. Default is None.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index when obtaining map boundary information
            from moving nests. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. Only required when
            *wrfin* is specified and the nest is moving.  Default is 0.

        stagger (:obj:`str`): If using latitude, longitude coordinate pairs
            for *start_point*, *end_point*, or *pivot_point*,
            set the appropriate grid staggering type for *field2d*. By default,
            the mass grid is used.  The options are:

                - 'm': Use the mass grid (default).
                - 'u': Use the same staggered grid as the u wind component,
                  which has a staggered west_east (x) dimension.
                - 'v': Use the same staggered grid as the v wind component,
                  which has a staggered south_north (y) dimension.

        projection (:class:`wrf.WrfProj`, optional): The map
            projection object to use when working with latitude, longitude
            coordinates, and must be specified if *wrfin* is None. Default
            is None.

        ll_point (:class:`wrf.CoordPair`, sequence of :class:`wrf.CoordPair`, \
        optional): The lower left latitude, longitude point for your domain,
            and must be specified
            if *wrfin* is None. If the domain is a moving nest, this should be
            a sequence of :class:`wrf.CoordPair`. Default is None.

    Returns:

        :class:`wrf.CoordPair` or sequence: The coordinate pair(s) in
        x,y grid coordinates.

    """

    if (wrfin is None and (projection is None or ll_point is None)):
        raise ValueError("'wrfin' parameter or "
                         "'projection' and 'll_point' parameters "
                         "are required")

    lat, lon = pairs_to_latlon(pairs)

    if wrfin is not None:
        xy_vals = _ll_to_xy(lat, lon, wrfin=wrfin, timeidx=timeidx,
                            squeeze=True, meta=False, stagger=stagger,
                            as_int=True)

    else:
        map_proj = projection.map_proj

        if map_proj == ProjectionTypes.LAT_LON:
            pole_lat = projection.pole_lat
            pole_lon = projection.pole_lon
            latinc = ((projection.dy*360.0)/2.0 /
                      Constants.PI/Constants.WRF_EARTH_RADIUS)
            loninc = ((projection.dx*360.0)/2.0 /
                      Constants.PI/Constants.WRF_EARTH_RADIUS)
        else:
            pole_lat = 90.0
            pole_lon = 0.0
            latinc = 0.0
            loninc = 0.0

        ll_lat, ll_lon = pairs_to_latlon(ll_point)
        xy_vals = _ll_to_xy(lat, lon, meta=False, squeeze=True,
                            as_int=True,
                            map_proj=projection.map_proj,
                            truelat1=projection.truelat1,
                            truelat2=projection.truelat2,
                            stand_lon=projection.stand_lon,
                            ref_lat=ll_lat,
                            ref_lon=ll_lon,
                            pole_lat=pole_lat,
                            pole_lon=pole_lon,
                            known_x=0,
                            known_y=0,
                            dx=projection.dx,
                            dy=projection.dy,
                            latinc=latinc,
                            loninc=loninc)

    xy_vals = xy_vals.squeeze()

    if xy_vals.ndim == 1:
        return CoordPair(x=xy_vals[0], y=xy_vals[1])
    else:
        return [CoordPair(x=xy_vals[0, i], y=xy_vals[1, i])
                for i in py3range(xy_vals.shape[1])]
