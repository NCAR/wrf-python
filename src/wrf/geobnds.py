from __future__ import (absolute_import, division, print_function)

from .coordpair import CoordPair


class GeoBounds(object):
    """A class that stores the geographic boundaries.

    Currently, only corner points are used, specified as the bottom left and
    top right corners.  Users can specify the corner points directly, or
    specify two-dimensional latitude and longitude arrays and the corner points
    will be extracted from them.

    Attributes:

        bottom_left (:class:`wrf.CoordPair`): The bottom left coordinate.
        top_right (:class:`wrf.CoordPair`): The top right coordinate.

    """
    def __init__(self, bottom_left=None, top_right=None, lats=None, lons=None):
        """ Initialize a :class:`wrf.GeoBounds` object.

        Args:

            bottom_left (:class:`wrf.CoordPair`, optional): The lower left
                corner. Must also specify *top_right* if used.
                Default is None.

            top_right (:class:`wrf.CoordPair`, optional): The upper right
                corner. Must also specify *bottom_left* if used.
                Default is None.

            lats (:class:`numpy.ndarray`, optional): An array of at least
                two dimensions containing all of the latitude values.  Must
                also specify *lons* if used.  Default is None.

            lons (:class:`numpy.ndarray`, optional): An array of at least
                two dimensions containing all of the longitude values.  Must
                also specify *lats* if used.  Default is None.

        """
        if bottom_left is not None and top_right is not None:
            self.bottom_left = bottom_left
            self.top_right = top_right

            # Make sure the users set lat/lon coordinates
            if self.bottom_left.lat is None:
                raise ValueError("'bottom_left' parameter does not contain a "
                                 "'lat' attribute")
            if self.bottom_left.lon is None:
                raise ValueError("'bottom_left' parameter does not contain a"
                                 "'lon' attribute")
            if self.top_right.lat is None:
                raise ValueError("'top_right' parameter does not contain a"
                                 "'lat' attribute")
            if self.top_right.lon is None:
                raise ValueError("'top_right' parameter does not contain a"
                                 "'lon' attribute")
        elif lats is not None and lons is not None:
            self.bottom_left = CoordPair(lat=lats[0, 0], lon=lons[0, 0])
            self.top_right = CoordPair(lat=lats[-1, -1], lon=lons[-1, -1])
        else:
            raise ValueError("must specify either 'bottom_top' and "
                             "'top_right' parameters "
                             "or 'lats' and 'lons' parameters")

    def __repr__(self):
        argstr = "{}, {}".format(repr(self.bottom_left),
                                 repr(self.top_right))

        return "{}({})".format(self.__class__.__name__, argstr)


class NullGeoBounds(GeoBounds):
    """An emtpy :class:`wrf.GeoBounds` subclass.

    This is used for initializing arrays of :class:`wrf.GeoBounds`, in
    particular when working with moving domains and variables combined with the
    'join' method.

    """
    def __init__(self):
        """ Initialize a :class:`wrf.NullGeoBounds` object."""
        pass

    def __repr__(self):
        return "{}()".format(self.__class__.__name__)
