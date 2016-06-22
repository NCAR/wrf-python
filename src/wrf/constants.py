from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

ALL_TIMES = None

class Constants(object):
    R = 287.06
    CP = 1005.0
    G = 9.81
    TCK0 = 273.15
    T_BASE = 300.0 # In WRF the base temperature is always 300 (not var T00)
    PI = 3.141592653589793
    DEFAULT_FILL = 9.9692099683868690e+36 # Value is from netcdf.h
    WRF_EARTH_RADIUS = 6370000.
    ERRLEN = 512

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
    
    