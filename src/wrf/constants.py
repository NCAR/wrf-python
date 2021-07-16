from __future__ import (absolute_import, division, print_function)

from sys import version_info
import struct
import numpy as np

from .py3compat import viewitems
from ._wrffortran import wrf_constants, omp_constants

#: Indicates that all times should be used in a diagnostic routine.
ALL_TIMES = None


class Constants(object):
    pass


for key, val in viewitems(wrf_constants.__dict__):
    setattr(Constants, key.upper(), np.array(val).item())

OMP_SCHED_STATIC = omp_constants.fomp_sched_static
OMP_SCHED_DYNAMIC = omp_constants.fomp_sched_dynamic
OMP_SCHED_GUIDED = omp_constants.fomp_sched_guided
OMP_SCHED_AUTO = omp_constants.fomp_sched_auto


class ConversionFactors(object):
    PA_TO_HPA = .01
    PA_TO_TORR = 760.0/101325.0
    PA_TO_MMHG = PA_TO_TORR * 1.000000142466321
    PA_TO_ATM = 1.0 / 1.01325E5
    MPS_TO_KTS = 1.94384
    MPS_TO_KMPH = 3.60
    MPS_TO_MPH = 2.23694
    MPS_TO_FPS = 3.28084
    M_TO_KM = 1.0/1000.0
    M_TO_DM = 1.0/10.0
    M_TO_FT = 3.28084
    M_TO_MILES = .000621371


class ProjectionTypes(object):
    ZERO = 0
    LAMBERT_CONFORMAL = 1
    POLAR_STEREOGRAPHIC = 2
    MERCATOR = 3
    LAT_LON = 6


# Create the default fill mapping based on type.
_DEFAULT_FILL_MAP = {None: Constants.DEFAULT_FILL,
                     np.dtype(np.bool_): False,
                     np.dtype(np.intc): Constants.DEFAULT_FILL_INT32,
                     np.dtype(np.int8): Constants.DEFAULT_FILL_INT8,
                     np.dtype(np.uint8): 255,
                     np.dtype(np.int16): Constants.DEFAULT_FILL_INT16,
                     np.dtype(np.uint16): 65535,
                     np.dtype(np.int32): Constants.DEFAULT_FILL_INT32,
                     np.dtype(np.uint32): 4294967295,
                     np.dtype(np.int64): Constants.DEFAULT_FILL_INT64,
                     np.dtype(np.uint64): 18446744073709551614,
                     np.dtype(np.float_): Constants.DEFAULT_FILL_DOUBLE,
                     np.dtype(np.float32): Constants.DEFAULT_FILL_FLOAT,
                     np.dtype(np.float64): Constants.DEFAULT_FILL_DOUBLE
                     }

if version_info >= (3, ):
    _DEFAULT_FILL_MAP[np.int_] = Constants.DEFAULT_FILL_INT64
else:
    _DEFAULT_FILL_MAP[np.int_] = Constants.DEFAULT_FILL_INT32

if (struct.calcsize("P") == 8):
    _DEFAULT_FILL_MAP[np.intp] = Constants.DEFAULT_FILL_INT64
else:
    _DEFAULT_FILL_MAP[np.intp] = Constants.DEFAULT_FILL_INT32


# Add the integers based on python 2.x or 3.x
def default_fill(dtype=None):
    dt = np.dtype(dtype) if dtype is not None else None
    return _DEFAULT_FILL_MAP.get(dt, Constants.DEFAULT_FILL)
