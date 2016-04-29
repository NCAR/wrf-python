

from .util import _get_iterable, is_standard_wrf_var, extract_vars
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

__all__ = ["getvar"]

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
             "ctt" : get_ctt
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
             "wspddir" : ["units"],
             "wspddir_uvmet" : ["units"],
             "wspddir_uvmet10" : ["units"],  
             "ctt" : [],
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
    for arg in kargs.iterkeys():
        if arg not in _VALID_KARGS[var]:
            raise ArgumentError("'%s' is an invalid keyword "
                          "argument for '%s'" % (arg, var))
              
  
def getvar(wrfnc, var, timeidx=0, 
           method="cat", squeeze=True, cache=None, meta=True, 
           **kargs):
    
    wrfnc = _get_iterable(wrfnc)
    
    if is_standard_wrf_var(wrfnc, var):
        return extract_vars(wrfnc, timeidx, var, 
                            method, squeeze, cache, meta)[var]
      
    actual_var = _undo_alias(var)
    if actual_var not in _VALID_KARGS:
        raise ArgumentError("'%s' is not a valid variable name" % (var))
      
    _check_kargs(actual_var, kargs)
    return _FUNC_MAP[actual_var](wrfnc, timeidx, 
                                 method, squeeze, cache, meta, **kargs)
