import warnings

from . import config
from .config import *
from . import extension
from .extension import *
from . import util
from .util import *
from . import cape
from .cape import *
from . import constants
from .constants import *
from . import ctt
from .ctt import *
from . import dbz
from .dbz import *
from . import destag
from .destag import *
from . import dewpoint
from .dewpoint import *
from . import geoht
from .geoht import *
from . import helicity
from .helicity import *
from . import interp
from .interp import *
from . import latlon
from .latlon import *
from . import omega
from .omega import *
from . import precip
from .precip import *
from . import pressure
from .pressure import *
from . import psadlookup
from .psadlookup import *
from . import pw
from .pw import *
from . import rh
from .rh import *
from . import slp
from .slp import *
from . import temp
from .temp import *
from . import terrain
from .terrain import *
from . import uvmet
from .uvmet import *
from . import vorticity
from .vorticity import *
from . import wind
from .wind import *
from . import times
from .times import *
from . import units
from .units import *
from . import projection
from .projection import *

__all__ = ["getvar"]
__all__.extend(config.__all__)
__all__.extend( extension.__all__)
__all__.extend(util.__all__)
__all__.extend(cape.__all__)
__all__.extend(constants.__all__)
__all__.extend(ctt.__all__)
__all__.extend(dbz.__all__)
__all__.extend(destag.__all__)
__all__.extend(dewpoint.__all__)
__all__.extend(geoht.__all__)
__all__.extend(helicity.__all__)
__all__.extend(interp.__all__)
__all__.extend(latlon.__all__)
__all__.extend(omega.__all__)
__all__.extend(precip.__all__)
__all__.extend(psadlookup.__all__)
__all__.extend(pw.__all__)
__all__.extend(rh.__all__)
__all__.extend(slp.__all__)
__all__.extend(temp.__all__)
__all__.extend(terrain.__all__)
__all__.extend(uvmet.__all__)
__all__.extend(vorticity.__all__)
__all__.extend(wind.__all__)
__all__.extend(times.__all__)
__all__.extend(pressure.__all__)
__all__.extend(units.__all__)
__all__.extend(projection.__all__)

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
           method="cat", squeeze=True, cache=None, 
           **kargs):
    if is_standard_wrf_var(wrfnc, var):
        return extract_vars(wrfnc, timeidx, var, method, squeeze, cache)[var]
    
    actual_var = _undo_alias(var)
    if actual_var not in _VALID_KARGS:
        raise ArgumentError("'%s' is not a valid variable name" % (var))
    
    _check_kargs(actual_var, kargs)
    return _FUNC_MAP[actual_var](wrfnc,timeidx,**kargs)
    

    
    
    
    