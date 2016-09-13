from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np
import numpy.ma as ma

from .constants import Constants
from .extension import (_interpz3d, _interp2dxy, _interp1d, _slp, _tk, _td, 
                        _rh, _uvmet, _smooth2d, _cape, _cloudfrac, _ctt, _dbz,
                        _srhel, _udhel, _avo, _pvo, _eth, _wetbulb, _tv, 
                        _omega, _pw)
from .util import from_var
from .metadecorators import (set_alg_metadata, set_uvmet_alg_metadata, 
                             set_interp_metadata, set_cape_alg_metadata,
                             set_cloudfrac_alg_metadata)
from .interputils import get_xy

@set_interp_metadata("xy")
def xy(field, pivot_point=None, angle=None, start_point=None, end_point=None,
       meta=True):
    return get_xy(field, pivot_point, angle, start_point, end_point)
    

@set_interp_metadata("1d")
def interp1d(field, z_in, z_out, missingval=Constants.DEFAULT_FILL, 
             meta=True):
    return _interp1d(field, z_in, z_out, missingval)


@set_interp_metadata("2dxy")
def interp2dxy(field3d, xy, meta=True):
    return _interp2dxy(field3d, xy)


@set_interp_metadata("horiz")
def interpz3d(field3d, z, desiredloc, missingval=Constants.DEFAULT_FILL,
              meta=True):
    return _interpz3d(field3d, z, desiredloc, missingval)


@set_alg_metadata(2, "pres", refvarndims=3, units="hpa",
                  description="sea level pressure")
def slp(height, tkel, pres, qv, meta=True):
    return _slp(height, tkel, pres, qv)


@set_alg_metadata(3, "pres", units="K",
                  description="temperature")
def tk(pres, theta, meta=True):
    return _tk(pres, theta)


@set_alg_metadata(3, "pres", units="degC",
                  description="dew point temperature")
def td(pres, qv, meta=True):
    return _td(pres, qv)


@set_alg_metadata(3, "pres", 
                  description="relative humidity", units=None)
def rh(qv, pres, tkel, meta=True):
    return _rh(qv, pres, tkel)


@set_uvmet_alg_metadata(latarg="lat", windarg="u")
def uvmet(u, v, lat, lon, cen_long, cone, meta=True):
    return _uvmet(u, v, lat, lon, cen_long, cone)


@set_alg_metadata(2, "field", 
                  description=from_var("field", "description"), 
                  units=from_var("field", "units"))
def smooth2d(field, passes, meta=True):
    return _smooth2d(field, passes)


@set_cape_alg_metadata(is2d=True, copyarg="pres_hpa")
def cape_2d(pres_hpa, tkel, qvapor, height, terrain, psfc_hpa, ter_follow, 
            missing=Constants.DEFAULT_FILL, meta=True):
    
    if isinstance(ter_follow, bool):
        ter_follow = 1 if ter_follow else 0
    
    i3dflag = 0
    cape_cin = _cape(pres_hpa, tkel, qvapor, height, terrain, psfc_hpa, 
                     missing, i3dflag, ter_follow)
    
    left_dims = cape_cin.shape[1:-3]
    right_dims = cape_cin.shape[-2:]
    
    resdim = (4,) + left_dims + right_dims
    
    # Make a new output array for the result
    result = np.zeros(resdim, cape_cin.dtype)
    
    # Cape 2D output is not flipped in the vertical, so index from the 
    # end
    result[0,...,:,:] = cape_cin[0,...,-1,:,:]
    result[1,...,:,:] = cape_cin[1,...,-1,:,:]
    result[2,...,:,:] = cape_cin[1,...,-2,:,:]
    result[3,...,:,:] = cape_cin[1,...,-3,:,:]
    
    return ma.masked_values(result, missing)


@set_cape_alg_metadata(is2d=False, copyarg="pres_hpa")
def cape_3d(pres_hpa, tkel, qvapor, height, terrain, psfc_hpa, ter_follow, 
            missing=Constants.DEFAULT_FILL, meta=True):
    
    if isinstance(ter_follow, bool):
        ter_follow = 1 if ter_follow else 0
    
    i3dflag = 1
    cape_cin = _cape(pres_hpa, tkel, qvapor, height, terrain, psfc_hpa, 
                     missing, i3dflag, ter_follow)
    
    return ma.masked_values(cape_cin, missing)


@set_cloudfrac_alg_metadata(copyarg="pres")
def cloudfrac(pres, relh, meta=True):
    return _cloudfrac(pres, relh)


@set_alg_metadata(2, "pres_hpa", refvarndims=3, units="degC",
                  description="cloud top temperature")
def ctt(pres_hpa, tkel, qv, qcld, height, terrain, qice=None, meta=True):
    
    # Qice and QCLD need to be in g/kg
    if qice is None:
        qice = np.zeros(qv.shape, qv.dtype)
        haveqci = 0
    else:
        haveqci = 1 if qice.any() else 0
    
    return _ctt(pres_hpa, tkel, qice, qcld, qv, height, terrain, haveqci)
    


@set_alg_metadata(3, "pres", units="dBZ",
                  description="radar reflectivity")
def dbz(pres, tkel, qv, qr, qs=None, qg=None, use_varint=False, use_liqskin=False, 
        meta=True):
    
    if qs is None:
        qs = np.zeros(qv.shape, qv.dtype)
        
    if qg is None:
        qg = np.zeros(qv.shape, qv.dtype)
    
    sn0 = 1 if qs.any() else 0
    ivarint = 1 if use_varint else 0
    iliqskin = 1 if use_liqskin else 0
    
    return _dbz(pres, tkel, qv, qr, qs, qg, sn0, ivarint, iliqskin)


@set_alg_metadata(2, "terrain", units="m-2/s-2",
                  description="storm relative helicity")
def srhel(u, v, z, terrain, top=3000.0, meta=True):
    # u, v get swapped in vertical
    _u = np.ascontiguousarray(u[...,::-1,:,:])
    _v = np.ascontiguousarray(v[...,::-1,:,:])
    _z = np.ascontiguousarray(z[...,::-1,:,:])
    
    return _srhel(_u, _v, _z, terrain, top)


@set_alg_metadata(2, "u", refvarndims=3, units="m-2/s-2",
                  description="updraft helicity")
def udhel(zstag, mapfct, u, v, wstag, dx, dy, bottom=2000.0, top=5000.0, 
          meta=True):
    return _udhel(zstag, mapfct, u, v, wstag, dx, dy, bottom, top)


# Requires both u an v for dimnames
@set_alg_metadata(3, "ustag", units="10-5 s-1",
                  stagdim=-1, stagsubvar="vstag",
                  description="absolute vorticity")
def avo(ustag, vstag, msfu, msfv, msfm, cor, dx, dy, meta=True):
    return _avo(ustag, vstag, msfu, msfv, msfm, cor, dx, dy)


@set_alg_metadata(3, "tkel", units="PVU",
                  description="potential vorticity")
def pvo(ustag, vstag, tkel, pres, msfu, msfv, msfm, cor, dx, dy, meta=True):
    return _pvo(ustag, vstag, tkel, pres, msfu, msfv, msfm, cor, dx, dy)


@set_alg_metadata(3, "qv", units="K",
                  description="equivalent potential temperature")
def eth(qv, tkel, pres, meta=True):
    return _eth(qv, tkel, pres)


@set_alg_metadata(3, "pres", units="K",
                  description="wetbulb temperature")
def wetbulb(pres, tkel, qv, meta=True):
    return _wetbulb(pres, tkel, qv)


@set_alg_metadata(3, "tkel", units="K",
                  description="virtual temperature")
def tvirtual(tkel, qv, meta=True):
    return _tv(tkel, qv)


@set_alg_metadata(3, "qv", units="Pa/s",
                  description="omega")
def omega(qv, tkel, w, pres, meta=True):
    return _omega(qv, tkel, w, pres)


@set_alg_metadata(2, "pres", refvarndims=3, units="kg m-2",
                  description="precipitable water")
def pw(pres, tkel, qv, height, meta=True):
    tv = _tv(tkel, qv)
    
    return _pw(pres, tv, qv, height)

    
    
