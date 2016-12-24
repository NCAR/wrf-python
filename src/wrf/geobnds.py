from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .coordpair import CoordPair

class GeoBounds(object):
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
        elif lats is not None and lons is not None:
            self.bottom_left = CoordPair(lat=lats[0,0], lon=lons[0,0])
            self.top_right = CoordPair(lat=lats[-1,-1], lon=lons[-1,-1])
        else:
            raise ValueError("invalid corner point arguments")
        
    def __repr__(self):
        argstr = "{}, {}".format(repr(self.bottom_left), 
                                 repr(self.top_right))
        
        return "{}({})".format(self.__class__.__name__, argstr)

 
class NullGeoBounds(GeoBounds):
    def __init__(self):
        pass
    
    def __repr__(self):
        return "{}()".format(self.__class__.__name__)
        
    
    
                