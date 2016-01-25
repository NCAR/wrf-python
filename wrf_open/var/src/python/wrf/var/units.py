
from wrf.var.constants import Constants, ConversionFactors

__all__ = ["check_units", "do_conversion"]

# Handles unit conversions that only differ by multiplication factors
def _apply_conv_fact(var, vartype, var_unit, dest_unit):
    if var_unit == dest_unit:
        return var
    
    # Note, case where var_unit and dest_unit are base unit, should be 
    # handled above
    if var_unit == _BASE_UNITS[vartype]:
        return var * _CONV_FACTORS[vartype]["to_dest"][dest_unit]
    else:
        if dest_unit == _BASE_UNITS[vartype]:
            return var*(_CONV_FACTORS[vartype]["to_base"][var_unit])
        else:
            return var*(_CONV_FACTORS[vartype]["to_base"][var_unit] * 
                    _CONV_FACTORS[vartype]["to_dest"][dest_unit])

def _to_celsius(var, var_unit):
    if var_unit == "k":
        return var - Constants.TCK0
    elif var_unit == "f":
        return (var - 32.0) * (5.0/9.0)
    
def _c_to_k(var):
    return var + Constants.TCK0

def _c_to_f(var):
    return 1.8 * var + 32.0

# Temperature is a more complicated operation so requres functions
def _apply_temp_conv(var, var_unit, dest_unit):
    if dest_unit == var_unit:
        return var
    
    if var_unit != _BASE_UNITS["temp"]:
        tc = _to_celsius(var, var_unit)
        if dest_unit == _BASE_UNITS["temp"]:
            return tc
        else:
            return (_TEMP_CONV_METHODS[dest_unit])(tc)
    else:
        return (_TEMP_CONV_METHODS[dest_unit])(var)

_VALID_UNITS = {"wind" : ["mps", "kts", "mph", "kmph", "fps"],
          "pressure" : ["pa", "hpa", "mb", "torr", "mmhg", "atm"],
          "temp" : ["k", "f", "c"],
          "height" : ["m", "km", "dm", "ft", "miles"]
          }

_WIND_BASE_FACTORS = {"kts" : ConversionFactors.MPS_TO_KTS,
                "kmph" : ConversionFactors.MPS_TO_KMPH,
                "mph" : ConversionFactors.MPS_TO_MPH,
                "fps" : ConversionFactors.MPS_TO_FPS
            }

_WIND_TOBASE_FACTORS = {"kts" : 1.0/ConversionFactors.MPS_TO_KTS,
                        "kmph" : 1.0/ConversionFactors.MPS_TO_KMPH,
                        "mph" : 1.0/ConversionFactors.MPS_TO_MPH,
                        "fps" : 1.0/ConversionFactors.MPS_TO_FPS
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
                        "miles" : ConversionFactors.M_TO_MILES
                        }

_HEIGHT_TOBASE_FACTORS = {"km" : 1.0/ConversionFactors.M_TO_KM,
                        "dm" : 1.0/ConversionFactors.M_TO_DM,
                        "ft" : 1.0/ConversionFactors.M_TO_FT,
                        "miles" : 1.0/ConversionFactors.M_TO_MILES
                        
                        }

_BASE_UNITS = {"wind" : "mps",
               "pressure" : "pa",
               "temp" : "c",
               "height" : "m"
               }

_CONV_FACTORS = {"wind" : {"to_dest" : _WIND_BASE_FACTORS,
                           "to_base" : _WIND_TOBASE_FACTORS},
                 "pressure" : {"to_dest" : _PRES_BASE_FACTORS,
                               "to_base" : _PRES_TOBASE_FACTORS},
                 "height" : {"to_dest" : _HEIGHT_BASE_FACTORS,
                             "to_base" : _HEIGHT_TOBASE_FACTORS}            
                 }

_TEMP_CONV_METHODS = {"k" : _c_to_k,
                      "f" : _c_to_f
                      }

def check_units(unit, unit_type):
    unitl = unit.lower()
    if unitl not in _VALID_UNITS[unit_type]:
        raise ValueError("invalid unit type '%s'" % unit)

def do_conversion(var, vartype, var_unit, dest_unit):
    if vartype != "temp":
        return _apply_conv_fact(var, vartype, 
                                var_unit.lower(), dest_unit.lower())
    else:
        return _apply_temp_conv(var, var_unit.lower(), dest_unit.lower())
    


    
    
