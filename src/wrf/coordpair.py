from __future__ import (absolute_import, division, print_function)

from .py3compat import py2round


def _binary_operator(operator):
    """Function wrapper for binary operators.

    Args:

        operator (method): The operator to wrap.

    Returns:

        method: An implementation for the *operator* type.

    """
    def func(self, other):
        """Operator implementation.

        Operator action is performed across the same class attributes when
        the *other* object is a :class:`CoordPair`. If the *other* object is
        a scalar value, the operator action is performed across all
        attributes with the scalar value.

        Args:

            other (:class:`CoordPair` or scalar): A separate :class:`CoordPair`
                object or scalar value.

        Returns:

            :class:`CoordPair`: A new :class:`CoordPair` object that is the
            result of the operator action.

        """
        if isinstance(other, CoordPair):
            args = [None if getattr(self, attr) is None or
                    getattr(other, attr) is None else
                    getattr(getattr(self, attr), operator)(getattr(other,
                                                                   attr))
                    for attr in ("x", "y", "lat", "lon")]
        else:
            args = [None if getattr(self, attr) is None
                    else getattr(getattr(self, attr), operator)(other)
                    for attr in ("x", "y", "lat", "lon")]

        return CoordPair(*args)

    return func


def _unary_operator(operator):
    """Function wrapper for unary operators.

    Args:

        operator (method): The operator to wrap.

    Returns:

        method: An implementation for the *operator* type.

    """
    def func(self):
        """Operator implementation.

        Operator action is performed across all class attributes.

        Returns:

            :class:`CoordPair`: A new :class:`CoordPair` object that is the
            result of the operator action.

        """
        args = [None if getattr(self, attr) is None
                else getattr(getattr(self, attr), operator)()
                for attr in ("x", "y", "lat", "lon")]

        return CoordPair(*args)

    return func


def _cmp_operator(operator):
    """Function wrapper for comparison operators.

    Args:

        operator (method): The operator to wrap.

    Returns:

        method: An implementation for the *operator* type.

    """

    def func(self, other):
        """Operator implementation.

        Performs a comparison operation across all of the same class
        attributes, and returns True if all these operations are True.

        Returns:

            :obj:`boot`: Returns True if all comparisons across class
            attributes returns True, otherwise False.

        """
        vals = [getattr(getattr(self, attr), operator)(getattr(other, attr))
                for attr in ("x", "y", "lat", "lon")
                if getattr(self, attr) is not None]

        return all(vals)

    return func


class CoordPair(object):
    """A class that stores (x, y) and/or (latitude, longitude)
    coordinate pairs.

    Most math operators are supplied.  When the other operand is a
    :class:`CoordPair`, the operation is performed with the same attribute.
    When a math operation uses a scalar as the other operand, the
    operation is applied across all attributes.

    Attributes:

        x (:obj:`float`): The x-coordinate.
        y (:obj:`float`): The y-coordinate.
        lat (:obj:`float`): The latitude coordinate.
        lon (:obj:`float`): The longitude coordinate.


    """
    def __init__(self, x=None, y=None, lat=None, lon=None):
        """Initialize a :class:`CoordPair` object.

        Args:

            x (:obj:`float`, optional): The x-coordinate.
            y (:obj:`float`, optional): The y-coordinate.
            lat (:obj:`float`, optional): The latitude coordinate.
            lon (:obj:`float`, optional): The longitude coordinate.


        """
        self.x = x
        self.y = y
        self.lat = lat
        self.lon = lon

    def __repr__(self):
        args = []
        if self.x is not None:
            args.append("x={}".format(self.x))
            args.append("y={}".format(self.y))

        if self.lat is not None:
            args.append("lat={}".format(self.lat))
            args.append("lon={}".format(self.lon))

        argstr = ", ".join(args)

        return "{}({})".format(self.__class__.__name__, argstr)

    def __str__(self):
        return self.__repr__()

    def xy_str(self, fmt="{:.4f}, {:.4f}"):
        """Return a :obj:`str` for the (x,y) coordinate pair.

        Args:

            fmt (:obj:`str`): The format string.  Default is '{:.4f}, {:.4f}'

        Returns:

            :obj:`str`: A string for the (x,y) coordinate pair

        """
        if self.x is None or self.y is None:
            return None

        return fmt.format(self.x, self.y)

    def latlon_str(self, fmt="{:.4f}, {:.4f}"):
        """Return a :obj:`str` for the (latitude, longitude) coordinate pair.

        Args:

            fmt (:obj:`str`): The format string.  Default is '{:.4f}, {:.4f}'

        Returns:

            :obj:`str`: A string for the (latitude, longitude) coordinate pair

        """
        if self.lat is None or self.lon is None:
            return None

        return fmt.format(self.lat, self.lon)

    def __round__(self, ndigits=None):
        """Return a new :class:`CoordPair` object with all coordinate values
        rounded to the nearest integer.

        Args:

            ndigits (:obj:`int`): The number of digits.

        Returns:

            :class:`CoordPair`: A CoordPair object.

        """
        args = [None if getattr(self, attr) is None
                else py2round(getattr(self, attr), ndigits)
                for attr in ("x", "y", "lat", "lon")]

        return CoordPair(*args)

    def __pow__(self, other, modulo=None):
        if isinstance(other, CoordPair):
            args = [None if getattr(self, attr) is None or
                    getattr(other, attr) is None
                    else getattr(getattr(self, attr), "__pow__")(
                        getattr(other, attr), modulo)
                    for attr in ("x", "y", "lat", "lon")]
        else:
            args = [None if getattr(self, attr) is None
                    else getattr(getattr(self, attr), "__pow__")(other, modulo)
                    for attr in ("x", "y", "lat", "lon")]

        return CoordPair(*args)

    def __rpow__(self, other):
        return self.__pow__(other)


for operator in ("__add__", "__divmod__", "__floordiv__", "__mod__",
                 "__mul__", "__sub__", "__truediv__", "__radd__",
                 "__rdivmod__", "__rsub__", "__rmul__", "__rtruediv__",
                 "__rfloordiv__", "__rmod__"):
    setattr(CoordPair, operator, _binary_operator(operator))


for operator in ("__neg__", "__pos__", "__abs__", "__invert__"):
    setattr(CoordPair, operator, _unary_operator(operator))


for operator in ("__lt__", "__le__", "__eq__", "__ne__", "__gt__", "__ge__"):
    setattr(CoordPair, operator, _cmp_operator(operator))
