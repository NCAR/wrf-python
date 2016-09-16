from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .py3compat import py2round
        

def _binary_operator(operator):
    def func(self, other):
        if isinstance(other, CoordPair):
            args = [
            None if getattr(self, attr) is None or getattr(other, attr) is None 
            else getattr(getattr(self, attr), operator)(getattr(other, attr))      
            for attr in ("x", "y", "lat", "lon")]
        else:
            args = [
            None if getattr(self, attr) is None 
            else getattr(getattr(self, attr), operator)(other)      
            for attr in ("x", "y", "lat", "lon")]
            
        return CoordPair(*args)
    
    return func


def _unary_operator(operator):
    def func(self):
        args = [None if getattr(self, attr) is None
                else getattr(getattr(self, attr), operator)()   
                for attr in ("x", "y", "lat", "lon")]
        
        return CoordPair(*args)
    
    return func


def _cmp_operator(operator):
    def func(self, other):
        vals = [getattr(getattr(self, attr), operator)(getattr(other, attr))
                for attr in ("x", "y", "lat", "lon") 
                if getattr(self, attr) is not None]
                
        return all(vals)
    
    return func
  
    
class CoordPair(object):
    def __init__(self, x=None, y=None, lat=None, lon=None):
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
        if self.x is None or self.y is None:
            return None
        
        return fmt.format(self.x, self.y)
    
    
    def latlon_str(self, fmt="{:.4f}, {:.4f}"):
        if self.lat is None or self.lon is None:
            return None
        
        return fmt.format(self.lat, self.lon)
     
       
    def __round__(self, d=None):
        args = [None if getattr(self, attr) is None
                else py2round(getattr(self, attr), d)  
                for attr in ("x", "y", "lat", "lon")]
        
        return CoordPair(*args)
    
    
    def __pow__(self, other, modulo=None):
        if isinstance(other, CoordPair):
            args = [
            None if getattr(self, attr) is None or getattr(other, attr) is None 
            else getattr(getattr(self, attr), "__pow__")(getattr(other, attr),
                                                        modulo)      
            for attr in ("x", "y", "lat", "lon")]
        else:
            args = [
            None if getattr(self, attr) is None 
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
    

    