from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .util import get_iterable, is_standard_wrf_var, extract_vars, viewkeys
from .cape import get_2dcape, get_3dcape
from .ctt import get_ctt
from .dbz import get_dbz, get_max_dbz
from .dewpoint import get_dp, get_dp_2m
from .geoht import get_geopt, get_height
from .helicity import get_srh, get_uh
from .latlon import get_lat, get_lon
from .omega import get_omega
from .pressure import get_pressure, get_pressure_hpa
from .pw import get_pw
from .rh import get_rh, get_rh_2m
from .slp import get_slp
from .temp import get_tc, get_eth, get_temp, get_theta, get_tk, get_tv, get_tw
from .terrain import get_terrain
from .uvmet import (get_uvmet, get_uvmet10, get_uvmet10_wspd_wdir, 
                    get_uvmet_wspd_wdir)
from .vorticity import get_avo, get_pvo
from .wind import (get_destag_wspd_wdir, get_destag_wspd_wdir10, 
                   get_u_destag, get_v_destag, get_w_destag)
from .times import get_times
from .cloudfrac import get_cloudfrac


# func is the function to call.  kargs are required arguments that should 
# not be altered by the user
_FUNC_MAP = {"cape2d" : get_2dcape,
             "cape3d" : get_3dcape,
             "dbz" : get_dbz,
             "maxdbz" : get_max_dbz,
             "dp" : get_dp,
             "dp2m" : get_dp_2m, 
             "height" : get_height,
             "geopt" : get_geopt,
             "srh" : get_srh,
             "uhel" : get_uh,
             "omega" : get_omega,
             "pw" : get_pw,
             "rh" : get_rh, 
             "rh2m" : get_rh_2m, 
             "slp" : get_slp, 
             "theta" : get_theta,
             "temp" : get_temp, 
             "tk" : get_tk,
             "tc" : get_tc,
             "theta_e" : get_eth, 
             "tv" : get_tv, 
             "twb" : get_tw, 
             "terrain" : get_terrain,
             "times" : get_times, 
             "uvmet" : get_uvmet, 
             "uvmet10" : get_uvmet10,
             "avo" : get_avo,
             "pvo" : get_pvo,
             "ua" : get_u_destag, 
             "va" : get_v_destag,
             "wa" : get_w_destag,
             "lat" : get_lat,
             "lon" : get_lon,
             "pressure" : get_pressure_hpa,
             "pres" : get_pressure,
             "wspd_wdir" : get_destag_wspd_wdir,
             "wspd_wdir10" : get_destag_wspd_wdir10,
             "wspd_wdir_uvmet" : get_uvmet_wspd_wdir,
             "wspd_wdir_uvmet10" : get_uvmet10_wspd_wdir,
             "ctt" : get_ctt,
             "cloudfrac" : get_cloudfrac
             }
  
_VALID_KARGS = {"cape2d" : ["missing"],
             "cape3d" : ["missing"],
             "dbz" : ["do_variant", "do_liqskin"],
             "maxdbz" : ["do_variant", "do_liqskin"],
             "dp" : ["units"],
             "dp2m" : ["units"], 
             "height" : ["msl", "units"],
             "geopt" : [],
             "srh" : ["top"],
             "uhel" : ["bottom", "top"],
             "omega" : [],
             "pw" : [],
             "rh" : [], 
             "rh2m" : [], 
             "slp" : ["units"], 
             "temp" : ["units"], 
             "tk" : [],
             "tc" : [],
             "theta" : ["units"],
             "theta_e" : ["units"], 
             "tv" : ["units"], 
             "twb" : ["units"], 
             "terrain" : ["units"],
             "times" : [], 
             "uvmet" : ["units"], 
             "uvmet10" : ["units"],
             "avo" : [],
             "pvo" : [],
             "ua" : ["units"], 
             "va" : ["units"],
             "wa" : ["units"],
             "lat" : [],
             "lon" : [],
             "pres" : ["units"],
             "pressure" : ["units"],
             "wspd_wdir" : ["units"],
             "wspd_wdir10" : ["units"],
             "wspd_wdir_uvmet" : ["units"],
             "wspd_wdir_uvmet10" : ["units"],  
             "ctt" : [],
             "cloudfrac" : [],
             "default" : []
            }
  
_ALIASES = {"cape_2d" : "cape2d",
            "cape_3d" : "cape3d",
            "eth" : "theta_e",
            "mdbz" : "maxdbz",
            "geopotential" : "geopt",
            "helicity" : "srh",
            "latitude" : "lat",
            "longitude" : "lon",
            "omg" : "omega",
            "p" : "pres",
            "rh2" : "rh2m",
            "z": "height",
            "ter" : "terrain",
            "updraft_helicity" : "uhel",
            "td" : "dp",
            "td2" : "dp2m"
            }
  
class ArgumentError(Exception):
    def __init__(self, msg):
        self.msg = msg
          
    def __str__(self):
        return self.msg
  
def _undo_alias(alias):
    actual = _ALIASES.get(alias, None)
    if actual is None:
        return alias
    else:
        return actual
  
def _check_kargs(var, kargs):
    for arg in viewkeys(kargs):
        if arg not in _VALID_KARGS[var]:
            raise ArgumentError("'%s' is an invalid keyword "
                          "argument for '%s'" % (arg, var))
              
  
def getvar(wrfnc, varname, timeidx=0, 
           method="cat", squeeze=True, cache=None, meta=True, 
           **kargs):
    
    """Returns basic diagnostics from the WRF ARW model output.
    
    Below is a table of available diagnostics.
    
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | Variable Name      | Description                                                   | Units               | Additional Keyword Arguments                                                                  |
    +====================+===============================================================+=====================+===============================================================================================+
    | avo                | Absolute Vorticity                                            | 10-5 s-1            |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | eth/theta_e        | Equivalent Potential Temperature                              | K                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | cape_2d            | 2D cape (mcape/mcin/lcl/lfc)                                  | J/kg / J/kg / m / m | missing: Fill value for output only (float)                                                   |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | cape_3d            | 3D cape and cin                                               | J/kg                | missing: Fill value for output only (float)                                                   |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | ctt                | Cloud Top Temperature                                         | C                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | cloudfrac          | Cloud Fraction                                                | %                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | dbz                | Reflectivity                                                  | dBz                 | do_variant: Set to True to enable variant calculation. Default is False.                      |
    |                    |                                                               |                     | do_liqskin : Set to True to enable liquid skin calculation. Default is False.                 |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | mdbz               | Maximum Reflectivity                                          | dBz                 | do_variant: Set to True to enable variant calculation. Default is False.                      |
    |                    |                                                               |                     | do_liqskin: Set to True to enable liquid skin calculation. Default is False.                  |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | geopt/geopotential | Full Model Geopotential                                       | m2 s-2              |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | helicity           | Storm Relative Helicity                                       | m-2/s-2             | top: The top level for the calculation in meters (float). Default is 3000.0.                  |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | lat                | Latitude                                                      | decimal degrees     |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | lon                | Longitude                                                     | decimal degrees     |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | omg/omega          | Omega                                                         | Pa/s                |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | p/pres             | Full Model Pressure                                           | Pa                  |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | pressure           | Full Model Pressure                                           | hPa                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | pvo                | Potential Vorticity                                           | PVU                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | pw                 | Precipitable Water                                            | kg m-2              |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | rh2                | 2m Relative Humidity                                          | %                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | slp                | Sea Level Pressure                                            | hPa                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | ter                | Model Terrain Height                                          | m                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | td2                | 2m Dew Point Temperature                                      | C                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | td                 | Dew Point Temperature                                         | C                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | tc                 | Temperature                                                   | C                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | th/theta           | Potential Temperature                                         | K                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | tk                 | Temperature                                                   | K                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | times              | Times in the File or Sequence                                 |                     |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | tv                 | Virtual Temperature                                           | K                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | twb                | Wet Bulb Temperature                                          | K                   |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | updraft_helicity   | Updraft Helicity                                              | m-2/s-2             | bottom: The bottom level for the calculation in meters (float). Default is 2000.0.            |                                                              
    |                    |                                                               |                     | top: The top level for the calculation in meters (float). Default is 5000.0.                  |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | ua                 | U-component of Wind on Mass Points                            | m/s                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | va                 | V-component of Wind on Mass Points                            | m/s                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | wa                 | W-component of Wind on Mass Points                            | m/s                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | uvmet10            | 10 m U and V Components of Wind Rotated to Earth Coordinates  | m/s                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | uvmet              | U and V Components of Wind Rotated to Earth Coordinates       | m/s                 |                                                                                               |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+
    | z/height           | Full Model Height                                             | m                   | msl: Set to False to return AGL values. Otherwise, MSL.  Default is True.                     |
    +--------------------+---------------------------------------------------------------+---------------------+-----------------------------------------------------------------------------------------------+


    Parameters
    ----------
    wrfnc : `netCD4F.Dataset`, `Nio.NioFile`, or a sequence
        Input WRF ARW NetCDF data as a `netCDF4.Dataset`, `Nio.NioFile` or an 
        iterable sequence of the aforementioned types.
    
    varname : str
        The variable name.
        
    timeidx : int or `wrf.ALL_TIMES`, optional
        The desired time index.  This value can be a positive integer, 
        negative integer, or `wrf.ALL_TIMES` (an alias for None) to return 
        all times in the file or sequence.  The default is 0.
        
    method : {'cat', 'join'}, optional
        The aggregation method to use for sequences, either 'cat' or 'join'.  
        'cat' combines the data along the Time index.  'join' is creates a new 
        index for the file.  The default is 'cat'.
        
    squeeze : {True, False}, optional
        Set to False to prevent dimensions with a size of 1 from being removed
        from the shape of the output.  Default is True.
        
    cache : dict, optional
        A dictionary of (varname, ndarray) which can be used to supply 
        pre-extracted NetCDF variables to the computational routines.  This can 
        be used to prevent the repeated variable extraction from large
        sequences of data files.  Default is None.
        
    meta : {True, False}, optional
        Set to False to disable metadata and return `numpy.ndarray` instead of 
        `xarray.DataArray`.  Default is True.
        
   
    Returns
    -------
    `xarray.DataArray` or `numpy.ndarray`
        Returns the specified diagnostics output.  If xarray is enabled and 
        the meta parameter is True, then the result will be an 
        `xarray.DataArray` object.  Otherwise, the result will be a 
        `numpy.ndarray` object with no metadata.
    
    """
    
    wrfnc = get_iterable(wrfnc)
    
    if is_standard_wrf_var(wrfnc, varname):
        return extract_vars(wrfnc, timeidx, varname, 
                            method, squeeze, cache, meta)[varname]
      
    actual_var = _undo_alias(varname)
    if actual_var not in _VALID_KARGS:
        raise ArgumentError("'%s' is not a valid variable name" % (varname))
      
    _check_kargs(actual_var, kargs)
    return _FUNC_MAP[actual_var](wrfnc, timeidx, 
                                 method, squeeze, cache, meta, **kargs)
    
