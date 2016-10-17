from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from threading import local
from collections import OrderedDict

from .config import get_cache_size

_local_storage = local()

def cache_item(key, product, value):
    """Store an item in the cache.
    
    The cache should be viewed as two nested dictionaries.  The outer key is 
    usually the id for the sequence where the cached item was generated.  The 
    inner key is the product type.
    
    Storing a cached item behaves like:
    
        cache[key][product] = value
    
    The cache is thread local, so stored items are only available in 
    the thread that cached them.
    
    Args:
    
        key (:obj:`int`): The outer dictionary cache key, which is typically 
            the id of the sequence where the cached item was generated.
            
        product (:obj:`str`): The inner dictionary cache key, which is a 
            string for the product type.
        
        value (:obj:`object`): The object to store in the cache.
    
    Returns:
    
        None.
        
    See Also:
    
        :meth:`get_cached_item`
        
    """
    global _local_storage
    
    if key is None:
        return
    
    try:
        cache = _local_storage.cache
    except AttributeError:
        _local_storage.cache = OrderedDict()
        cache = _local_storage.cache
    
    try:
        _ = cache[key]
    except KeyError:
        if len(cache) >= get_cache_size():
            cache.popitem(last=False) # Remove the oldest dataset
        
        cache[key] = OrderedDict()
        
    cache[key][product] = value
    
    
def get_cached_item(key, product):
    """Return an item from the cache.
    
    The cache should be viewed as two nested dictionaries.  The outer key is 
    usually the id for the sequence where the cached item was generated.  The 
    inner key is the product type.
    
    Retrieving a cached item behaves like:
    
        value = cache[key][product]
    
    The cache is thread local, so stored items are only available in 
    the thread that cached them.
    
    Args:
    
        key (:obj:`int`): The outer dictionary cache key, which is typically 
            the id of the sequence where the cached item was generated.
            
        product (:obj:`str`): The inner dictionary cache key, which is a 
            string for the product type.
    
    Returns:
    
        :obj:`object`: The cached object.
        
    See Also:
    
        :meth:`cache_item`
        
    """
    if key is None:
        return None
    
    cache = getattr(_local_storage, "cache", None)
    
    if cache is None:
        return None
    
    prod_dict = cache.get(key, None)
    
    if prod_dict is None:
        return None
    
    result = prod_dict.get(product, None)
        
    return result

def _get_cache():
    """Return the cache.
    
    This is primarily used for testing.
    
    Returns:
    
        :class:`threading.local`
    
    """
    return getattr(_local_storage, "cache", None)
    
    
    
        
    