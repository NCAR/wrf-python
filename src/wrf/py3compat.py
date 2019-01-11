from __future__ import (absolute_import, division, print_function)

from sys import version_info
from math import floor, copysign


# Dictionary python 2-3 compatibility stuff
def viewitems(d):
    """Return either the items or viewitems method for a dictionary.

    Args:

        d (:obj:`dict`): A dictionary.

    Returns:

        view method: Either the items or viewitems method.

    """
    func = getattr(d, "viewitems", None)
    if func is None:
        func = d.items
    return func()


def viewkeys(d):
    """Return either the keys or viewkeys method for a dictionary.

    Args:

        d (:obj:`dict`): A dictionary.

    Returns:

        view method: Either the keys or viewkeys method.

    """
    func = getattr(d, "viewkeys", None)
    if func is None:
        func = d.keys
    return func()


def viewvalues(d):
    """Return either the values or viewvalues method for a dictionary.

    Args:

        d (:obj:`dict`): A dictionary.

    Returns:

        view method: Either the values or viewvalues method.

    """
    func = getattr(d, "viewvalues", None)
    if func is None:
        func = d.values
    return func()


def isstr(s):
    """Return True if the object is a string type.

    Args:

        s (string): A string (str, unicode, bytes).

    Returns:

        :obj:`bool`: True if the object is a type of string. Otherwise, False.

    """
    try:
        return isinstance(s, basestring)
    except NameError:
        return isinstance(s, str)


# Python 2 rounding behavior
def _round2(x, d=None):
    """Return the result of Python 2.x rounding, which is to round the number
    to the nearest integer.

    Python 3.x uses banker's rounding, which is not applicable for nearest
    neighbor approaches with grid boxes.

    Args:

        x (:obj:`float`): A number, usually a float.

        d (:obj:`int`, optional): The number of digits.  Default is None,
            which indicates the nearest integer.

    Returns:

        :obj:`float`: The rounded number.

    """
    d = 0 if d is None else d
    p = 10 ** d
    return float(floor((x * p) + copysign(0.5, x)))/p


def py2round(x, d=None):
    """Return the result of Python 2.x rounding, which is to round the number
    to the nearest integer.

    Python 3.x uses banker's rounding, which is not applicable for nearest
    neighbor approaches with grid boxes.

    Args:

        x (:obj:`float`): A number, usually a float.

        d (:obj:`int`, optional): The number of digits.  Default is None,
            which indicates the nearest integer.

    Returns:

        :obj:`float`: The rounded number.

    """
    if version_info >= (3,):
        return _round2(x, d)

    return round(x, d)


def py3range(*args):
    """Return the equivalent of the range function in Python 3.x.

    For Python 2.x, this is the same as the xrange function.

    Args:

        *args: The function arguments for range or xrange.

    Returns:

        iterable: An iterable sequence.

    """
    if version_info >= (3,):
        return range(*args)

    return xrange(*args)


def ucode(*args, **kwargs):
    """Return a Python 3.x unicode string.

    For Python 2.x, this is accomplished by using the unicode function.

    Args:

        *args: The function positional arguments for str or unicode.

        **kwargs: The function keyword arguments for str or unicode.

    Returns:

        string: A unicode string.

    """
    if version_info >= (3, ):
        return str(*args, **kwargs)

    return unicode(*args, **kwargs)
