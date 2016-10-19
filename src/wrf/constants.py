from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

from .py3compat import viewitems
from wrf._wrffortran import wrf_constants

#: Indicates that all times should be used in a diagnostic routine.
ALL_TIMES = None

class Constants(object):
    pass
    
for key,val in viewitems(wrf_constants.__dict__):
    setattr(Constants, key.upper(), np.asscalar(val))

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
    
    