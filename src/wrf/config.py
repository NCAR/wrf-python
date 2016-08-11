from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

try:
    from xarray import DataArray
except ImportError:
    _XARRAY_ENABLED = False
else:
    _XARRAY_ENABLED = True
    
try:
    from cartopy import crs
except ImportError:
    _CARTOPY_ENABLED = False
else:
    _CARTOPY_ENABLED = True
    
try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    _BASEMAP_ENABLED = False
else:
    _BASEMAP_ENABLED = True
    
try:
    from Ngl import Resources
except ImportError:
    _PYNGL_ENABLED = False
else:
    _PYNGL_ENABLED = True

def xarray_enabled():
    global _XARRAY_ENABLED
    return _XARRAY_ENABLED

def disable_xarray():
    global _XARRAY_ENABLED
    _XARRAY_ENABLED = False
    
def enable_xarray():
    global _XARRAY_ENABLED
    _XARRAY_ENABLED = True
    
def cartopy_enabled():
    global _CARTOPY_ENABLED
    return _CARTOPY_ENABLED

def enable_cartopy():
    global _CARTOPY_ENABLED
    _CARTOPY_ENABLED = True
    
def disable_cartopy():
    global _CARTOPY_ENABLED
    _CARTOPY_ENABLED = True
    
def basemap_enabled():
    global _BASEMAP_ENABLED
    return _BASEMAP_ENABLED

def enable_basemap():
    global _BASEMAP_ENABLED
    _BASEMAP_ENABLED = True
    
def disable_basemap():
    global _BASEMAP_ENABLED
    _BASEMAP_ENABLED = True
    
def pyngl_enabled():
    global _PYNGL_ENABLED
    return _PYNGL_ENABLED

def enable_pyngl():
    global _PYNGL_ENABLED
    _PYNGL_ENABLED = True
    
def disable_pyngl():
    global _PYNGL_ENABLED
    _PYNGL_ENABLED = True
