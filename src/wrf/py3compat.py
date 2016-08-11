from sys import version_info

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

def isstr(s):
    try:
        return isinstance(s, basestring)
    except NameError:
        return isinstance(s, str)


# Python 2 rounding behavior  
def _round2(x, d=0):
    p = 10 ** d
    return float(floor((x * p) + copysign(0.5, x)))/p


def py2round(x, d=0):
    if version_info >= (3,):
        return _round2(x, d)
    
    return round(x, d)


def py3range(*args):
    if version_info >= (3,):
        return range(*args)
    
    return xrange(*args)


def ucode(*args, **kwargs):
    if version_info >= (3, ):
        return str(*args, **kwargs)
    
    return unicode(*args, **kwargs)