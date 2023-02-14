from __future__ import (absolute_import, division, print_function)

from .util import (get_iterable, is_standard_wrf_var, extract_vars, viewkeys,
                   get_id)
from .g_cape import (get_2dcape, get_3dcape, get_cape2d_only,
                     get_cin2d_only, get_lcl, get_lfc, get_3dcape_only,
                     get_3dcin_only)
from .g_ctt import get_ctt
from .g_dbz import get_dbz, get_max_dbz
from .g_dewpoint import get_dp, get_dp_2m
from .g_geoht import (get_geopt, get_height, get_stag_geopt, get_stag_height,
                      get_height_agl)
from .g_helicity import get_srh, get_uh
from .g_latlon import get_lat, get_lon
from .g_omega import get_omega
from .g_pressure import get_pressure, get_pressure_hpa
from .g_pw import get_pw
from .g_rh import get_rh, get_rh_2m
from .g_slp import get_slp
from .g_temp import (get_tc, get_eth, get_temp, get_theta, get_tk, get_tv,
                     get_tw)
from .g_terrain import get_terrain
from .g_uvmet import (get_uvmet, get_uvmet10, get_uvmet10_wspd_wdir,
                      get_uvmet_wspd_wdir, get_uvmet_wspd, get_uvmet_wdir,
                      get_uvmet10_wspd, get_uvmet10_wdir)
from .g_vorticity import get_avo, get_pvo
from .g_wind import (get_destag_wspd_wdir, get_destag_wspd_wdir10,
                     get_u_destag, get_v_destag, get_w_destag,
                     get_destag_wspd, get_destag_wdir, get_destag_wspd10,
                     get_destag_wdir10)
from .g_times import get_times, get_xtimes
from .g_cloudfrac import (get_cloudfrac, get_low_cloudfrac, get_mid_cloudfrac,
                          get_high_cloudfrac)


# func is the function to call.  kargs are required arguments that should
# not be altered by the user
_FUNC_MAP = {"cape2d": get_2dcape,
             "cape3d": get_3dcape,
             "dbz": get_dbz,
             "maxdbz": get_max_dbz,
             "dp": get_dp,
             "dp2m": get_dp_2m,
             "height": get_height,
             "geopt": get_geopt,
             "srh": get_srh,
             "uhel": get_uh,
             "omega": get_omega,
             "pw": get_pw,
             "rh": get_rh,
             "rh2m": get_rh_2m,
             "slp": get_slp,
             "theta": get_theta,
             "temp": get_temp,
             "tk": get_tk,
             "tc": get_tc,
             "theta_e": get_eth,
             "tv": get_tv,
             "twb": get_tw,
             "terrain": get_terrain,
             "times": get_times,
             "xtimes": get_xtimes,
             "uvmet": get_uvmet,
             "uvmet10": get_uvmet10,
             "avo": get_avo,
             "pvo": get_pvo,
             "ua": get_u_destag,
             "va": get_v_destag,
             "wa": get_w_destag,
             "lat": get_lat,
             "lon": get_lon,
             "pressure": get_pressure_hpa,
             "pres": get_pressure,
             "wspd_wdir": get_destag_wspd_wdir,
             "wspd_wdir10": get_destag_wspd_wdir10,
             "uvmet_wspd_wdir": get_uvmet_wspd_wdir,
             "uvmet10_wspd_wdir": get_uvmet10_wspd_wdir,
             "ctt": get_ctt,
             "cloudfrac": get_cloudfrac,
             "geopt_stag": get_stag_geopt,
             "zstag": get_stag_height,
             "height_agl": get_height_agl,
             # Diagnostics below are extracted from multi-product diagnostics
             "cape2d_only": get_cape2d_only,
             "cin2d_only": get_cin2d_only,
             "lcl": get_lcl,
             "lfc": get_lfc,
             "cape3d_only": get_3dcape_only,
             "cin3d_only": get_3dcin_only,
             "uvmet_wspd": get_uvmet_wspd,
             "uvmet_wdir": get_uvmet_wdir,
             "uvmet10_wspd": get_uvmet10_wspd,
             "uvmet10_wdir": get_uvmet10_wdir,
             "wspd": get_destag_wspd,
             "wdir": get_destag_wdir,
             "wspd10": get_destag_wspd10,
             "wdir10": get_destag_wdir10,
             "low_cloudfrac": get_low_cloudfrac,
             "mid_cloudfrac": get_mid_cloudfrac,
             "high_cloudfrac": get_high_cloudfrac
             }

_VALID_KARGS = {"cape2d": ["missing"],
                "cape3d": ["missing"],
                "dbz": ["do_variant", "do_liqskin"],
                "maxdbz": ["do_variant", "do_liqskin"],
                "dp": ["units"],
                "dp2m": ["units"],
                "height": ["msl", "units"],
                "geopt": [],
                "srh": ["top"],
                "uhel": ["bottom", "top"],
                "omega": [],
                "pw": [],
                "rh": [],
                "rh2m": [],
                "slp": ["units"],
                "temp": ["units"],
                "tk": [],
                "tc": [],
                "theta": ["units"],
                "theta_e": ["units"],
                "tv": ["units"],
                "twb": ["units"],
                "terrain": ["units"],
                "times": [],
                "xtimes": [],
                "uvmet": ["units"],
                "uvmet10": ["units"],
                "avo": [],
                "pvo": [],
                "ua": ["units"],
                "va": ["units"],
                "wa": ["units"],
                "lat": [],
                "lon": [],
                "pres": ["units"],
                "pressure": ["units"],
                "wspd_wdir": ["units"],
                "wspd_wdir10": ["units"],
                "uvmet_wspd_wdir": ["units"],
                "uvmet10_wspd_wdir": ["units"],
                "ctt": ["fill_nocloud", "missing", "opt_thresh", "units"],
                "cloudfrac": ["vert_type", "low_thresh",
                              "mid_thresh", "high_thresh"],
                "geopt_stag": [],
                "zstag": ["msl", "units"],
                "height_agl": ["units"],
                "cape2d_only": ["missing"],
                "cin2d_only": ["missing"],
                "lcl": ["missing"],
                "lfc": ["missing"],
                "cape3d_only": ["missing"],
                "cin3d_only": ["missing"],
                "uvmet_wspd": ["units"],
                "uvmet_wdir": ["units"],
                "uvmet10_wspd": ["units"],
                "uvmet10_wdir": ["units"],
                "wspd": ["units"],
                "wdir": ["units"],
                "wspd10": ["units"],
                "wdir10": ["units"],
                "low_cloudfrac": ["vert_type", "low_thresh",
                                  "mid_thresh", "high_thresh"],
                "mid_cloudfrac": ["vert_type", "low_thresh",
                                  "mid_thresh", "high_thresh"],
                "high_cloudfrac": ["vert_type", "low_thresh",
                                   "mid_thresh", "high_thresh"],
                "default": []
                }

_ALIASES = {"cape_2d": "cape2d",
            "cape_3d": "cape3d",
            "eth": "theta_e",
            "mdbz": "maxdbz",
            "geopotential": "geopt",
            "helicity": "srh",
            "latitude": "lat",
            "longitude": "lon",
            "omg": "omega",
            "p": "pres",
            "rh2": "rh2m",
            "z": "height",
            "ter": "terrain",
            "updraft_helicity": "uhel",
            "td": "dp",
            "td2": "dp2m",
            "cfrac": "cloudfrac",
            "wspd_wdir_uvmet": "uvmet_wspd_wdir",
            "wspd_wdir_uvmet10": "uvmet10_wspd_wdir",
            "th": "theta",
            "low_cfrac": "low_cloudfrac",
            "mid_cfrac": "mid_cloudfrac",
            "high_cfrac": "high_cloudfrac",
            "wspd_uvmet": "uvmet_wspd",
            "wdir_uvmet": "uvmet_wdir",
            "wspd_uvmet10": "uvmet10_wspd",
            "wdir_uvmet10": "uvmet10_wdir",
            "mcape": "cape2d_only",
            "mcin": "cin2d_only"
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
            if var != "default":
                raise ValueError("'{}' is an invalid keyword "
                                 "argument for '{}'".format(arg, var))
            else:
                raise ValueError("'{}' is an invalid keyword "
                                 "argument".format(arg))


def getvar(wrfin, varname, timeidx=0,
           method="cat", squeeze=True, cache=None, meta=True,
           **kwargs):

    """Returns basic diagnostics from the WRF ARW model output.

    A table of all available diagnostics is below.

    .. include:: ../../_templates/product_table.txt

    Args:
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
            or an iterable sequence of the aforementioned types.

        varname (:obj:`str`) : The variable name.

        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
            desired time index. This value can be a positive integer,
            negative integer, or
            :data:`wrf.ALL_TIMES` (an alias for None) to return
            all times in the file or sequence. The default is 0.

        method (:obj:`str`, optional): The aggregation method to use for
            sequences.  Must be either 'cat' or 'join'.
            'cat' combines the data along the Time dimension.
            'join' creates a new dimension for the file index.
            The default is 'cat'.

        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
            with a size of 1 from being automatically removed from the shape
            of the output. Default is True.

        cache (:obj:`dict`, optional): A dictionary of (varname, ndarray)
            that can be used to supply pre-extracted NetCDF variables to the
            computational routines.  It is primarily used for internal
            purposes, but can also be used to improve performance by
            eliminating the need to repeatedly extract the same variables
            used in multiple diagnostics calculations, particularly when using
            large sequences of files.
            Default is None.

        meta (:obj:`bool`, optional): Set to False to disable metadata and
            return :class:`numpy.ndarray` instead of
            :class:`xarray.DataArray`.  Default is True.

        **kwargs: Optional keyword arguments for certain diagnostics.
            See table above.


    Returns:

        :class:`xarray.DataArray` or :class:`numpy.ndarray`: If xarray is
        enabled and the *meta* parameter is True, then the result will be a
        :class:`xarray.DataArray` object.  Otherwise, the result will be a
        :class:`numpy.ndarray` object with no metadata.


    Raises:
        :class:`ValueError`: Raised when an invalid diagnostic type or
            keyword argument is passed to the routine.
        :class:`FortranError`: Raised when a problem occurs during a Fortran
            calculation.

    See Also:

        :class:`numpy.ndarray`, :class:`xarray.DataArray`


    Examples:
        Using netCDF4

        .. code-block:: python

            from netCDF4 import Dataset
            from wrf import getvar

            wrfnc = Dataset("wrfout_d02_2010-06-13_21:00:00")
            slp = getvar(wrfnc, "slp")

        Using PyNIO

        .. code-block:: python

            from Nio import open_file
            from wrf import getvar

            wrfnc = open_file("wrfout_d02_2010-06-13_21:00:00"+".nc", "r")
            slp = getvar(wrfnc, "slp")

        Using Iterables:

        .. code-block:: python

            import os
            from netCDF4 import Dataset
            from wrf import getvar

            filedir = "/path/to/wrf/files"
            wrfin = [Dataset(f) for f in os.listdir(filedir)
                    if f.startswith("wrfout_d02_")]

            uvmet = getvar(wrfin, "uvmet", timeidx=3, units="kt")


    """

    _key = get_id(wrfin)

    wrfin = get_iterable(wrfin)

    if is_standard_wrf_var(wrfin, varname) and varname != "Times":
        _check_kargs("default", kwargs)
        return extract_vars(wrfin, timeidx, varname,
                            method, squeeze, cache, meta, _key)[varname]
    elif varname == "Times":
        varname = "times"  # Diverting to the get_times routine

    actual_var = _undo_alias(varname)
    if actual_var not in _VALID_KARGS:
        raise ValueError("'{}' is not a valid variable name".format(varname))

    _check_kargs(actual_var, kwargs)

    return _FUNC_MAP[actual_var](wrfin, timeidx, method, squeeze, cache,
                                 meta, _key, **kwargs)
