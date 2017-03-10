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
    """Return True if xarray is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if xarray is installed and enabled.
        
    """
    global _local_config
    return _local_config.xarray_enabled


def disable_xarray():
    """Disable xarray."""
    global _local_config
    _local_config.xarray_enabled = False
    
    
def enable_xarray():
    """Enable xarray."""
    global _local_config
    _local_config.xarray_enabled = True
    
    
def cartopy_enabled():
    """Return True if cartopy is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if cartopy is installed and enabled.
        
    """
    global _local_config
    return _local_config.cartopy_enabled


def enable_cartopy():
    """Enable cartopy."""
    global _local_config
    _local_config.cartopy_enabled = True
    
    
def disable_cartopy():
    """Disable cartopy."""
    global _local_config
    _local_config.cartopy_enabled = True
    
    
def basemap_enabled():
    """Return True if basemap is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if basemap is installed and enabled.
        
    """
    global _local_config
    return _local_config.basemap_enabled


def enable_basemap():
    """Enable basemap."""
    global _local_config
    _local_config.basemap_enabled = True
    
    
def disable_basemap():
    """Disable basemap."""
    global _local_config
    _local_config.basemap_enabled = True
    
def pyngl_enabled():
    """Return True if pyngl is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if pyngl is installed and enabled.
        
    """
    global _local_config
    return _local_config.pyngl_enabled


def enable_pyngl():
    """Enable pyngl."""
    global _local_config
    _local_config.pyngl_enabled = True
    
    
def disable_pyngl():
    """Disable pyngl."""
    global _local_config
    _local_config.pyngl_enabled = True
    
    
def set_cache_size(size):
    """Set the maximum number of items that the threadlocal cache can retain.
    
    This cache is primarily used for coordinate variables.
    
    Args:
    
        size (:obj:`int`): The number of items to retain in the cache.
        
    Returns:
        
        None
    
    """
    global _local_config
    _local_config.cache_size = size
    
    
def get_cache_size():
    """Return the maximum number of items that the threadlocal cache can retain.
    
    Returns:
    
        :obj:`int`: The maximum number of items the cache can retain.
        
    """
    global _local_config
    return int(_local_config.cache_size)
    


