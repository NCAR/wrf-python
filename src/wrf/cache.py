from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from threading import local
from collections import OrderedDict

from .config import get_cache_size

_local_storage = local()

def cache_item(key, product, value):
    global _local_storage
    
    try:
        cache = _local_storage.cache
    except AttributeError:
        _local_storage.cache = OrderedDict()
        cache = _local_storage.cache
    
    try:
        prod_dict = cache[key]
    except KeyError:
        if len(cache) >= get_cache_size():
            cache.popitem(last=False) # Remove the oldest dataset
        
        cache[key] = OrderedDict()
        prod_dict = cache[key]
        
    cache[key][product] = value
    
    
def get_cached_item(key, product):
    cache = getattr(_local_storage, "cache", None)
    
    if cache is None:
        return None
    
    prod_dict = cache.get(key, None)
    
    if prod_dict is None:
        return None
    
    return prod_dict.get(product, None)


def _get_cache():
    return getattr(_local_storage, "cache", None)
    
    
    
        
    