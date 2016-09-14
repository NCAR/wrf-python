from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants, ConversionFactors


# Handles unit conversions that only differ by multiplication factors
def _apply_conv_fact(var, vartype, var_unit, dest_unit):
    if var_unit == dest_unit:
        return var
    
    # Note, case where var_unit and dest_unit are base unit, should be 
    # handled above
    if var_unit == _BASE_UNITS[vartype]:
        return var*(_CONV_FACTORS[vartype]["to_dest"][dest_unit])
    else:
        if dest_unit == _BASE_UNITS[vartype]:
            return var*(_CONV_FACTORS[vartype]["to_base"][var_unit])
        else:
            return var*(_CONV_FACTORS[vartype]["to_base"][var_unit] * 
                    _CONV_FACTORS[vartype]["to_dest"][dest_unit])


def _to_kelvin(var, var_unit):
    if var_unit == "c":
        return var + Constants.CELKEL
    elif var_unit == "f":
        return (var - 32.0) * (5.0/9.0) + Constants.CELKEL
    
    
def _k_to_c(var):
    return var - Constants.CELKEL


def _k_to_f(var):
    return 1.8 * _k_to_c(var) + 32.0


# Temperature is a more complicated operation so requires functions
def _apply_temp_conv(var, var_unit, dest_unit):
    if dest_unit == var_unit:
        return var
    
    if var_unit != _BASE_UNITS["temp"]:
        tk = _to_kelvin(var, var_unit)
        if dest_unit == _BASE_UNITS["temp"]:
            return tk
        else:
            return (_TEMP_CONV_METHODS[dest_unit])(tk)
    else:
        return (_TEMP_CONV_METHODS[dest_unit])(var)
    

_UNIT_ALIASES = {"mps" : "m s-1",
                 "m/s" : "m s-1",
                 "ms-1" : "m s-1",
                 "meters_per_second" : "m s-1",
                 "metres_per_second" : "m s-1",
                 "knots" : "kt",
                 "knot" : "kt",
                 "kts" : "kt",
                 "kn" : "kt",
                 "miles_per_hour" : "mi h-1",
                 "mih-1" : "mi h-1",
                 "mph" : "mi h-1",
                 "mi/h" : "mi h-1",
                 "kmph" : "km h-1",
                 "kmh-1" : "km h-1",
                 "km/h" : "km h-1",
                 "kilometers_per_hour" : "km h-1",
                 "kilometres_per_hour" : "km h-1",
                 "ft/s" : "ft s-1",
                 "ft/sec" : "ft s-1",
                 "fps" : "ft s-1",
                 "fs-1" : "ft s-1",
                 "feet_per_second" : "ft s-1",
                 
                 "pascal" : "pa",
                 "pascals" : "pa",
                 "hecto_pascal" : "hpa",
                 "hecto_pascals" : "hpa",
                 "millibar" : "mb",
                 "millibars" : "mb",
                 "mbar" : "mb",
                 
                 "kelvin" : "k",
                 "degree_kelvin" : "k",
                 "degrees_kelvin" : "k",
                 "degree_k" : "k",
                 "degrees_k" : "k",
                 "degreek" : "k",
                 "degreesk" : "k",
                 "degk" : "k",
                 "degsk" : "k",
                 "deg_k" : "k",
                 "degs_k" : "k",
                 "deg k" : "k",
                 "degs k" : "k",
                 
                 "celsius" : "c",
                 "degree_celsius" : "c",
                 "degrees_celsius" : "c",
                 "degree_c" : "c",
                 "degrees_c" : "c",
                 "degreec" : "c",
                 "degreesc" : "c",
                 "degc" : "c",
                 "degsc" : "c",
                 "deg_c" : "c",
                 "degs_c" : "c",
                 "deg c" : "c",
                 "degs c" : "c",

                 "fahrenheit" : "f",
                 "degree_fahrenheit" : "f",
                 "degrees_fahrenheit" : "f",
                 "degree_f" : "f",
                 "degrees_f" : "f",
                 "degreef" : "f",
                 "degreesf" : "f",
                 "degf" : "f",
                 "degsf" : "f",
                 "deg_f" : "f",
                 "degs_f" : "f",
                 "deg f" : "f",
                 "degs f" : "f",
                 
                 "meter" : "m",
                 "meters" : "m",
                 "metre" : "m",
                 "metres" : "m",
                 "kilometer" : "km",
                 "kilometers" : "km",
                 "dekameter" : "dm",
                 "dekameters" : "dm",
                 "decameter" : "dm",
                 "decameters" : "dm",
                 "dekametre" : "dm",
                 "dekametres" : "dm",
                 "decametre" : "dm",
                 "decametres" : "dm",
                 "dam" : "dm",
                 "dkm" : "dm",
                 "feet" : "ft",
                 "foot" : "ft",
                 "mile" : "mi",
                 "miles" : "mi"
    
    }

_VALID_UNITS = {"wind" : ["m s-1", "kt", "mi h-1", "km h-1", "ft s-1"],
          "pressure" : ["pa", "hpa", "mb", "torr", "mmhg", "atm"],
          "temp" : ["k", "f", "c"],
          "height" : ["m", "km", "dm", "ft", "mi"]
          }

_WIND_BASE_FACTORS = {"kt" : ConversionFactors.MPS_TO_KTS,
                "km h-1" : ConversionFactors.MPS_TO_KMPH,
                "mi h-1" : ConversionFactors.MPS_TO_MPH,
                "ft s-1" : ConversionFactors.MPS_TO_FPS
            }

_WIND_TOBASE_FACTORS = {"kt" : 1.0/ConversionFactors.MPS_TO_KTS,
                        "km h-1" : 1.0/ConversionFactors.MPS_TO_KMPH,
                        "mi h-1" : 1.0/ConversionFactors.MPS_TO_MPH,
                        "ft s-1" : 1.0/ConversionFactors.MPS_TO_FPS
                        }

_PRES_BASE_FACTORS = {"hpa" : ConversionFactors.PA_TO_HPA,
                      "mb" : ConversionFactors.PA_TO_HPA,
                      "torr" : ConversionFactors.PA_TO_TORR,
                      "mmhg" : ConversionFactors.PA_TO_MMHG,
                      "atm" : ConversionFactors.PA_TO_ATM
                      }

_PRES_TOBASE_FACTORS = {"hpa" : 1.0/ConversionFactors.PA_TO_HPA,
                        "mb" : 1.0/ConversionFactors.PA_TO_HPA,
                        "torr" : 1.0/ConversionFactors.PA_TO_TORR,
                        "mmhg" : 1.0/ConversionFactors.PA_TO_MMHG,
                        "atm" : 1.0/ConversionFactors.PA_TO_ATM
                        }

_HEIGHT_BASE_FACTORS = {"km" : ConversionFactors.M_TO_KM,
                        "dm" : ConversionFactors.M_TO_DM,
                        "ft" : ConversionFactors.M_TO_FT,
                        "mi" : ConversionFactors.M_TO_MILES
                        }

_HEIGHT_TOBASE_FACTORS = {"km" : 1.0/ConversionFactors.M_TO_KM,
                        "dm" : 1.0/ConversionFactors.M_TO_DM,
                        "ft" : 1.0/ConversionFactors.M_TO_FT,
                        "mi" : 1.0/ConversionFactors.M_TO_MILES
                        }

_BASE_UNITS = {"wind" : "m s-1",
               "pressure" : "pa",
               "temp" : "k",
               "height" : "m"
               }

_CONV_FACTORS = {"wind" : {"to_dest" : _WIND_BASE_FACTORS,
                           "to_base" : _WIND_TOBASE_FACTORS},
                 "pressure" : {"to_dest" : _PRES_BASE_FACTORS,
                               "to_base" : _PRES_TOBASE_FACTORS},
                 "height" : {"to_dest" : _HEIGHT_BASE_FACTORS,
                             "to_base" : _HEIGHT_TOBASE_FACTORS}            
                 }

_TEMP_CONV_METHODS = {"c" : _k_to_c,
                      "f" : _k_to_f
                      }

def dealias_and_clean_unit(unit):
    cleaned_unit = " ".join(unit.lower().split())
    dealiased = _UNIT_ALIASES.get(cleaned_unit, None)
    
    return cleaned_unit if dealiased is None else dealiased
    

def check_units(unit, unit_type):
    u_cleaned = dealias_and_clean_unit(unit)
    if u_cleaned not in _VALID_UNITS[unit_type]:
        raise ValueError("invalid unit type '%s'" % unit)


def do_conversion(var, vartype, var_unit, dest_unit):
    u_cleaned = dealias_and_clean_unit(dest_unit) 
    if vartype != "temp":
        return _apply_conv_fact(var, vartype, var_unit.lower(), u_cleaned)
    else:
        return _apply_temp_conv(var, var_unit.lower(), u_cleaned)
    


    
    
