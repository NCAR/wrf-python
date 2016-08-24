from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

try:
    from xarray import DataArray
except ImportError:
    _xarray_enabled = False
else:
    _xarray_enabled = True
    
try:
    from cartopy import crs
except ImportError:
    _cartopy_enabled = False
else:
    _cartopy_enabled = True
    
try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    _basemap_enabled = False
else:
    _basemap_enabled = True
    
try:
    from Ngl import Resources
except ImportError:
    _pyngl_enabled = False
else:
    _pyngl_enabled = True
    
_cache_size = 5

def xarray_enabled():
    global _xarray_enabled
    return _xarray_enabled

def disable_xarray():
    global _xarray_enabled
    _xarray_enabled = False
    
def enable_xarray():
    global _xarray_enabled
    _xarray_enabled = True
    
def cartopy_enabled():
    global _cartopy_enabled
    return _cartopy_enabled

def enable_cartopy():
    global _cartopy_enabled
    _cartopy_enabled = True
    
def disable_cartopy():
    global _cartopy_enabled
    _cartopy_enabled = True
    
def basemap_enabled():
    global _basemap_enabled
    return _basemap_enabled

def enable_basemap():
    global _basemap_enabled
    _basemap_enabled = True
    
def disable_basemap():
    global _basemap_enabled
    _basemap_enabled = True
    
def pyngl_enabled():
    global _pyngl_enabled
    return _pyngl_enabled

def enable_pyngl():
    global _pyngl_enabled
    _pyngl_enabled = True
    
def disable_pyngl():
    global _pyngl_enabled
    _pyngl_enabled = True
    
def set_cache_size(size):
    global _cache_size
    _cache_size = size
    
def get_cache_size():
    return int(_cache_size)
    


