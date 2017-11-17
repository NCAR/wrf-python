from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from threading import local
import wrapt

_local_config = local()

def _init_local():
    global _local_config
    
    _local_config.xarray_enabled = True
    _local_config.cartopy_enabled = True
    _local_config.basemap_enabled = True
    _local_config.pyngl_enabled = True
    _local_config.cache_size = 20
    _local_config.initialized = True
    
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
        

# Initialize the main thread's configuration
_init_local()


def init_local():
    """A decorator that initializes thread local data if necessary."""
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        global _local_config
        try:
            init = _local_config.init
        except AttributeError:
            _init_local()
        else:
            if not init:
                _init_local() 
                
        return wrapped(*args, **kwargs)
        
    return func_wrapper


@init_local()
def xarray_enabled():
    """Return True if xarray is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if xarray is installed and enabled.
        
    """
    global _local_config
    return _local_config.xarray_enabled


@init_local()
def disable_xarray():
    """Disable xarray."""
    global _local_config
    _local_config.xarray_enabled = False
    

@init_local()
def enable_xarray():
    """Enable xarray."""
    global _local_config
    _local_config.xarray_enabled = True
    

@init_local()
def cartopy_enabled():
    """Return True if cartopy is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if cartopy is installed and enabled.
        
    """
    global _local_config
    return _local_config.cartopy_enabled


@init_local()
def enable_cartopy():
    """Enable cartopy."""
    global _local_config
    _local_config.cartopy_enabled = True
    

@init_local()
def disable_cartopy():
    """Disable cartopy."""
    global _local_config
    _local_config.cartopy_enabled = True
    

@init_local()
def basemap_enabled():
    """Return True if basemap is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if basemap is installed and enabled.
        
    """
    global _local_config
    return _local_config.basemap_enabled


@init_local()
def enable_basemap():
    """Enable basemap."""
    global _local_config
    _local_config.basemap_enabled = True
    

@init_local()
def disable_basemap():
    """Disable basemap."""
    global _local_config
    _local_config.basemap_enabled = True


@init_local()
def pyngl_enabled():
    """Return True if pyngl is installed and enabled.
    
    Returns:
    
        :obj:`bool`: True if pyngl is installed and enabled.
        
    """
    global _local_config
    return _local_config.pyngl_enabled


@init_local()
def enable_pyngl():
    """Enable pyngl."""
    global _local_config
    _local_config.pyngl_enabled = True
    

@init_local()  
def disable_pyngl():
    """Disable pyngl."""
    global _local_config
    _local_config.pyngl_enabled = True
    

@init_local() 
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
    

@init_local()  
def get_cache_size():
    """Return the maximum number of items that the threadlocal cache can retain.
    
    Returns:
    
        :obj:`int`: The maximum number of items the cache can retain.
        
    """
    global _local_config
    return int(_local_config.cache_size)
    


