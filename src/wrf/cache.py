from __future__ import (absolute_import, division, print_function)

from threading import local
from collections import OrderedDict

from .py3compat import py3range
from .config import get_cache_size

_local_storage = local()


def _shrink_cache():
    """Shrink the cache if applicable.

    This only applies if a user has modified the cache size, otherwise it
    just returns.

    Returns:

        None

    """
    global _local_storage

    try:
        cache = _local_storage.cache
    except AttributeError:
        return

    diff = len(cache) - get_cache_size()

    if diff > 0:
        for _ in py3range(diff):
            cache.popitem(last=False)


def cache_item(key, product, value):
    """Store an item in the threadlocal cache.

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

    _shrink_cache()

    if key is None or get_cache_size() == 0:
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
            cache.popitem(last=False)  # Remove the oldest dataset

        cache[key] = OrderedDict()

    cache[key][product] = value


def get_cached_item(key, product):
    """Return an item from the threadlocal cache.

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
    global _local_storage

    _shrink_cache()

    if key is None or get_cache_size == 0:
        return None

    cache = getattr(_local_storage, "cache", None)

    if cache is None:
        return None

    if len(cache) == 0:
        return None

    prod_dict = cache.get(key, None)

    if prod_dict is None:
        return None

    result = prod_dict.get(product, None)

    return result


def _get_cache():
    """Return the threadlocal cache.

    This is primarily used for testing.

    Returns:

        :class:`threading.local`

    """
    global _local_storage

    _shrink_cache()
    return getattr(_local_storage, "cache", None)
