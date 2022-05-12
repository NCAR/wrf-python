from __future__ import (absolute_import, division, print_function)

from threading import local
import wrapt

from ._wrffortran import (fomp_enabled, fomp_set_num_threads,
                          fomp_set_schedule, fomp_set_dynamic,
                          omp_constants)

_local_config = local()


def _try_enable_xarray():
    global _local_config
    _local_config.xarray_enabled = True
    try:
        from xarray import DataArray
    except ImportError:
        _local_config.xarray_enabled = False


def _try_enable_cartopy():
    global _local_config
    _local_config.cartopy_enabled = True
    try:
        from cartopy import crs
    except ImportError:
        _local_config.cartopy_enabled = False


def _try_enable_basemap():
    global _local_config
    _local_config.basemap_enabled = True
    try:
        from mpl_toolkits.basemap import Basemap
    except ImportError:
        _local_config.basemap_enabled = False


def _try_enable_pyngl():
    global _local_config
    _local_config.pyngl_enabled = True
    try:
        from Ngl import Resources
    except ImportError:
        _local_config.pyngl_enabled = False


def _init_local():
    global _local_config

    _try_enable_xarray()
    _try_enable_cartopy()
    _try_enable_basemap()
    _try_enable_pyngl()

    _local_config.cache_size = 20
    _local_config.initialized = True


# Initialize the main thread's configuration
_init_local()


def init_local():
    """A decorator that initializes thread local data if necessary."""
    @wrapt.decorator
    def func_wrapper(wrapped, instance, args, kwargs):
        global _local_config
        try:
            initialized = _local_config.initialized
        except AttributeError:
            _init_local()
        else:
            if not initialized:
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
def enable_xarray():
    """Enable xarray if it is installed."""
    _try_enable_xarray()


@init_local()
def disable_xarray():
    """Disable xarray."""
    global _local_config
    _local_config.xarray_enabled = False


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
    """Enable cartopy if it is installed."""
    _try_enable_cartopy()


@init_local()
def disable_cartopy():
    """Disable cartopy."""
    global _local_config
    _local_config.cartopy_enabled = False


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
    """Enable basemap if it is installed."""
    _try_enable_basemap()


@init_local()
def disable_basemap():
    """Disable basemap."""
    global _local_config
    _local_config.basemap_enabled = False


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
    """Enable pyngl if it is installed."""
    _try_enable_pyngl()


@init_local()
def disable_pyngl():
    """Disable pyngl."""
    global _local_config
    _local_config.pyngl_enabled = False


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
    """Return the maximum number of items that the threadlocal cache can
    retain.

    Returns:

        :obj:`int`: The maximum number of items the cache can retain.

    """
    global _local_config
    return int(_local_config.cache_size)


def omp_enabled():
    """Return True if OpenMP is enabled.

    OpenMP is only enabled if compiled with OpenMP features.

    Returns:

        :obj:`bool`: True if OpenMP is enabled, otherwise False.

    """

    return True if fomp_enabled() else False


# Set OpenMP to use 1 thread, static scheduler, and no dynamic
# Note: Using the raw extension functions here to prevent possible
# circular import problems in the future.
fomp_set_num_threads(1)
fomp_set_schedule(omp_constants.fomp_sched_static, 0)
fomp_set_dynamic(False)
