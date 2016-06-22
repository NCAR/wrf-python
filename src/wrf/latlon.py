from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .util import (extract_vars)
from .latlonutils import (_lat_varname, _lon_varname, ll_to_ij, ij_to_ll)
from .metadecorators import set_latlon_metadata


def get_lat(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True,
            stagger=None):
    
    varname = _lat_varname(wrfnc, stagger)
    lat_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           meta)
    
    return lat_var[varname]

        
def get_lon(wrfnc, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True,
            stagger=None):
    
    varname = _lon_varname(wrfnc, stagger)
    lon_var = extract_vars(wrfnc, timeidx, varname, method, squeeze, cache,
                           meta)
    
    return lon_var[varname]


# TODO:  Rename these!
@set_latlon_metadata(ij=True) 
def get_ij(wrfnc, latitude, longitude, timeidx=0,
           stagger=None, method="cat", squeeze=True, cache=None, meta=True):
    return ll_to_ij(wrfnc, latitude, longitude, timeidx, stagger, 
                    method, squeeze, cache)
    
    
@set_latlon_metadata(ij=False) 
def get_ll(wrfnc, i, j, timeidx=0,
           stagger=None, method="cat", squeeze=True, cache=None, meta=True):
    return ij_to_ll(wrfnc, i, j, timeidx, stagger, 
                    method, squeeze, cache)
    
    
    