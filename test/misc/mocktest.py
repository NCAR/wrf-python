import sys
import os

try:
    from unittest.mock import MagicMock
except ImportError:
    from mock import Mock as MagicMock


class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return Mock()


MOCK_MODULES = ["numpy", "numpy.ma", "xarray", "cartopy",
                "pandas", "matplotlib", "netCDF4", "mpl_toolkits.basemap",
                "wrf._wrffortran"]
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

consts = {"DEFAULT_FILL": 9.9692099683868690E36,
          "DEFAULT_FILL_INT8": -127,
          "DEFAULT_FILL_INT16": -32767,
          "DEFAULT_FILL_INT32": -2147483647,
          "DEFAULT_FILL_INT64": -9223372036854775806,
          "DEFAULT_FILL_FLOAT": 9.9692099683868690E36,
          "DEFAULT_FILL_DOUBLE": 9.9692099683868690E36,
          "fomp_sched_static": 1,
          "fomp_sched_dynamic": 2,
          "fomp_sched_guided": 3,
          "fomp_sched_auto": 4}


class MockWrfConstants(object):
    def __init__(self):
        self.__dict__ = consts


def mock_asscalar(val):
    return float(val)


sys.modules["wrf._wrffortran"].wrf_constants = MockWrfConstants()
sys.modules["wrf._wrffortran"].omp_constants = MockWrfConstants()

sys.modules["numpy"].asscalar = mock_asscalar

try:
    import wrf
except ImportError:
    pass

print(wrf.get_coord_pairs.__doc__)
