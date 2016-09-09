from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from threading import local

_local_config = local()
_local_config.xarray_enabled = True
_local_config.cartopy_enabled = True
_local_config.basemap_enabled = True
_local_config.pyngl_enabled = True
_local_config.cache_size = 20

try:
    from xarray import DataArray
except ImportError:
    _local_config.xarray_enabled = False
    
try:
    from cartopy import crs
except ImportError:
    _local_config.cartopy_enabled = False
    
try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    _local_config.basemap_enabled = False
    
try:
    from Ngl import Resources
except ImportError:
    _local_config.pyngl_enabled = False


def xarray_enabled():
    global _local_config
    return _local_config.xarray_enabled


def disable_xarray():
    global _local_config
    _local_config.xarray_enabled = False
    
    
def enable_xarray():
    global _local_config
    _local_config.xarray_enabled = True
    
    
def cartopy_enabled():
    global _local_config
    return _local_config.cartopy_enabled


def enable_cartopy():
    global _local_config
    _local_config.cartopy_enabled = True
    
    
def disable_cartopy():
    global _local_config
    _local_config.cartopy_enabled = True
    
    
def basemap_enabled():
    global _local_config
    return _local_config.basemap_enabled


def enable_basemap():
    global _local_config
    _local_config.basemap_enabled = True
    
    
def disable_basemap():
    global _local_config
    _local_config.basemap_enabled = True
    
def pyngl_enabled():
    global _local_config
    return _local_config.pyngl_enabled


def enable_pyngl():
    global _local_config
    _local_config.pyngl_enabled = True
    
    
def disable_pyngl():
    global _local_config
    _local_config.pyngl_enabled = True
    
    
def set_cache_size(size):
    global _local_config
    _local_config.cache_size = size
    
    
def get_cache_size():
    global _local_config
    return int(_local_config.cache_size)
    


