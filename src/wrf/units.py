from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants, ConversionFactors


def _apply_conv_fact(var, vartype, var_unit, dest_unit):
    """Return the variable converted to different units using a conversion 
    factor.
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            variable.
            
        vartype (:obj:`str`): The type of variable.  Choices are: 'wind', 
            'pressure', 'temp', or 'height'.
            
        var_unit (:obj:`str`): The variable's current units.
        
        dest_unit (:obj:`str`): The desired units.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The variable in
        the desired units.  
    
    """
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
    """Return the variable in Kelvin.
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            variable.
            
        var_unit (:obj:`str`): The variable's current units.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The variable in 
        Kelvin.
    
    """
    if var_unit == "c":
        return var + Constants.CELKEL
    elif var_unit == "f":
        return (var - 32.0) * (5.0/9.0) + Constants.CELKEL
    
    
def _k_to_c(var):
    """Return the variable in Celsius.
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            variable in units of Kelvin.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The variable in 
        Celsius.
    
    """
    return var - Constants.CELKEL


def _k_to_f(var):
    """Return the variable in Fahrenheit.
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            variable in units of Kelvin.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The variable in 
        Fahrenheit.
    
    """
    return 1.8 * _k_to_c(var) + 32.0


def _apply_temp_conv(var, var_unit, dest_unit):
    """Return the variable converted to different units using a temperature
    conversion algorithm.
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            variable.
            
        var_unit (:obj:`str`): The variable's current units.
        
        dest_unit (:obj:`str`): The desired units.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The variable in
        the desired units.  
    
    """
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
    

# A mapping of unit names to their dictionary key names
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

# A mapping of unit types to the avaible units
_VALID_UNITS = {"wind" : ["m s-1", "kt", "mi h-1", "km h-1", "ft s-1"],
          "pressure" : ["pa", "hpa", "mb", "torr", "mmhg", "atm"],
          "temp" : ["k", "f", "c"],
          "height" : ["m", "km", "dm", "ft", "mi"]
          }

# Conversion factor map for wind from base units
_WIND_BASE_FACTORS = {"kt" : ConversionFactors.MPS_TO_KTS,
                "km h-1" : ConversionFactors.MPS_TO_KMPH,
                "mi h-1" : ConversionFactors.MPS_TO_MPH,
                "ft s-1" : ConversionFactors.MPS_TO_FPS
            }

# Conversion factor map to base units
_WIND_TOBASE_FACTORS = {"kt" : 1.0/ConversionFactors.MPS_TO_KTS,
                        "km h-1" : 1.0/ConversionFactors.MPS_TO_KMPH,
                        "mi h-1" : 1.0/ConversionFactors.MPS_TO_MPH,
                        "ft s-1" : 1.0/ConversionFactors.MPS_TO_FPS
                        }

# Conversion factor map for pressure from base units
_PRES_BASE_FACTORS = {"hpa" : ConversionFactors.PA_TO_HPA,
                      "mb" : ConversionFactors.PA_TO_HPA,
                      "torr" : ConversionFactors.PA_TO_TORR,
                      "mmhg" : ConversionFactors.PA_TO_MMHG,
                      "atm" : ConversionFactors.PA_TO_ATM
                      }

# Conversion factor map for pressure to base units
_PRES_TOBASE_FACTORS = {"hpa" : 1.0/ConversionFactors.PA_TO_HPA,
                        "mb" : 1.0/ConversionFactors.PA_TO_HPA,
                        "torr" : 1.0/ConversionFactors.PA_TO_TORR,
                        "mmhg" : 1.0/ConversionFactors.PA_TO_MMHG,
                        "atm" : 1.0/ConversionFactors.PA_TO_ATM
                        }

# Conversion factor map for height from base units
_HEIGHT_BASE_FACTORS = {"km" : ConversionFactors.M_TO_KM,
                        "dm" : ConversionFactors.M_TO_DM,
                        "ft" : ConversionFactors.M_TO_FT,
                        "mi" : ConversionFactors.M_TO_MILES
                        }

# Conversion factor map for height to base units
_HEIGHT_TOBASE_FACTORS = {"km" : 1.0/ConversionFactors.M_TO_KM,
                        "dm" : 1.0/ConversionFactors.M_TO_DM,
                        "ft" : 1.0/ConversionFactors.M_TO_FT,
                        "mi" : 1.0/ConversionFactors.M_TO_MILES
                        }

# Mapping of unit type to base unit type
_BASE_UNITS = {"wind" : "m s-1",
               "pressure" : "pa",
               "temp" : "k",
               "height" : "m"
               }

# A mapping of unit type to a mapping of to/from base conversion factors
_CONV_FACTORS = {"wind" : {"to_dest" : _WIND_BASE_FACTORS,
                           "to_base" : _WIND_TOBASE_FACTORS},
                 "pressure" : {"to_dest" : _PRES_BASE_FACTORS,
                               "to_base" : _PRES_TOBASE_FACTORS},
                 "height" : {"to_dest" : _HEIGHT_BASE_FACTORS,
                             "to_base" : _HEIGHT_TOBASE_FACTORS}            
                 }

# A mapping of temperature type to the conversion function
_TEMP_CONV_METHODS = {"c" : _k_to_c,
                      "f" : _k_to_f
                      }

def dealias_and_clean_unit(unit):
    """Return the properly cleaned and dealiased unit name.
    
    Args:
    
        unit (:obj:`str`): The unit name.
        
    Returns:
    
        :obj:`str`: A unit name suitable for dictionary key lookups.
    
    """
    cleaned_unit = " ".join(unit.lower().split())
    dealiased = _UNIT_ALIASES.get(cleaned_unit, None)
    
    return cleaned_unit if dealiased is None else dealiased
    

def check_units(unit, unit_type):
    """Raise an exception if the unit name is invalid.
    
    Args:
        
        unit (:obj:`str`):  The unit name.
        
        unit_type (:obj:`str`): The type of unit.
        
    Returns:
    
        None
        
    Raises:
    
        :class:`ValueError`: Raised when the unit name is invalid.
    
    """
    u_cleaned = dealias_and_clean_unit(unit)
    if u_cleaned not in _VALID_UNITS[unit_type]:
        raise ValueError("invalid unit type '%s'" % unit)


def do_conversion(var, vartype, var_unit, dest_unit):
    """Return the variable converted to different units.
    
    Args:
    
        var (:class:`xarray.DataArray` or :class:`numpy.ndarray`): A 
            variable.
            
        vartype (:obj:`str`): The type of variable.  Choices are: 'wind', 
            'pressure', 'temp', or 'height'.
            
        var_unit (:obj:`str`): The variable's current units.
        
        dest_unit (:obj:`str`): The desired units.
        
    Returns:
    
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The variable in
        the desired units.  
    
    """
    u_cleaned = dealias_and_clean_unit(dest_unit) 
    if vartype != "temp":
        return _apply_conv_fact(var, vartype, var_unit.lower(), u_cleaned)
    else:
        return _apply_temp_conv(var, var_unit.lower(), u_cleaned)
    


    
    
