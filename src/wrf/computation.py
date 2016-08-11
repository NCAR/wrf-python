from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np
import numpy.ma as ma

from .constants import Constants
from .extension import (_interpz3d, _interp2dxy, _interp1d, _slp, _tk, _td, 
                        _rh, _uvmet, _smooth2d, _cape, _cloudfrac, _ctt, _dbz,
                        _srhel, _udhel, _eth, _wetbulb, _tv, _omega)
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
def interp1d(v_in, z_in, z_out, missingval=Constants.DEFAULT_FILL, 
             meta=True):
    return _interp1d(v_in, z_in, z_out, missingval)


@set_interp_metadata("2dxy")
def interp2dxy(field3d, xy, meta=True):
    return _interp2dxy(field3d, xy)


@set_interp_metadata("horiz")
def interpz3d(field3d, z, desiredloc, missingval=Constants.DEFAULT_FILL,
              meta=True):
    return _interpz3d(field3d, z, desiredloc, missingval)


@set_alg_metadata(2, "z", refvarndims=3, units="hpa",
                  description="sea level pressure")
def slp(z, t, p, q, meta=True):
    return _slp(z, t, p, q)


@set_alg_metadata(3, "pressure", units="K",
                  description="temperature")
def tk(pressure, theta, meta=True):
    return _tk(pressure, theta)


@set_alg_metadata(3, "pressure", units="degC",
                  description="dew point temperature")
def td(pressure, qv_in, meta=True):
    return _td(pressure, qv_in)


@set_alg_metadata(3, "pressure", 
                  description="relative humidity", units=None)
def rh(qv, q, t, meta=True):
    return _rh(qv, q, t, meta)


@set_uvmet_alg_metadata()
def uvmet(u, v, lat, lon, cen_long, cone, meta=True):
    return _uvmet(u, v, lat, lon, cen_long, cone)


@set_alg_metadata(2, "field", 
                  description=from_var("field", "description"), 
                  units=from_var("field", "units"))
def smooth2d(field, passes, meta=True):
    return _smooth2d(field, passes)


@set_cape_alg_metadata(is2d=True)
def cape_2d(p_hpa, tk, qvapor, height, terrain, psfc_hpa, ter_follow, 
            missing=Constants.DEFAULT_FILL, meta=True):
    
    if isinstance(ter_follow, bool):
        ter_follow = 1 if ter_follow else 0
    
    i3dflag = 0
    cape_cin = _cape(p_hpa, tk, qvapor, height, terrain, psfc_hpa, 
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


@set_cape_alg_metadata(is2d=False)
def cape_3d(p_hpa, tk, qvapor, height, terrain, psfc_hpa, ter_follow, 
            missing=Constants.DEFAULT_FILL, meta=True):
    
    if isinstance(ter_follow, bool):
        ter_follow = 1 if ter_follow else 0
    
    i3dflag = 1
    cape_cin = _cape(p_hpa, tk, qvapor, height, terrain, psfc_hpa, 
                     missing, i3dflag, ter_follow)
    
    return ma.masked_values(cape_cin, missing)


@set_cloudfrac_alg_metadata()
def cloudfrac(p, rh, meta=True):
    return _cloudfrac(p, rh)


@set_alg_metadata(2, "p_hpa", refvarndims=3, units="degC",
                  description="cloud top temperature")
def ctt(p_hpa, tk, qv, qcld, ght, ter, qice=None, meta=True):
    
    # Qice and QCLD need to be in g/kg
    if qice is None:
        qice = n.zeros(qv.shape, qv.dtype)
        haveqci = 0
    else:
        haveqci = 1 if qice.any() else 0
    
    ctt = _ctt(p_hpa, tk, qice, qcld, qv, ght, ter, haveqci)
    return _ctt(p, rh)


@set_alg_metadata(3, "p", units="dBZ",
                  description="radar reflectivity")
def dbz(p, tk, qv, qr, ivarint, iliqskin, qs=None, qg=None, meta=True):
    
    if qs is None:
        qs = np.zeros(qv.shape, qv.dtype)
        
    if qg is None:
        qg = np.zeros(qv.shape, qv.dtype)
    
    sn0 = 1 if qs.any() else 0
    ivarint = 1 if do_varint else 0
    iliqskin = 1 if do_liqskin else 0
    
    return _dbz(p, tk, qv, qr, qs, qg, sn0, ivarint, iliqskin)


@set_alg_metadata(2, "ter", units="m-2/s-2",
                  description="storm relative helicity")
def srhel(u, v, z, ter, top):
    return _srhel(u, v, z, ter, top=3000.0)


@set_alg_metadata(2, "u", refvarndims=3, units="m-2/s-2",
                  description="updraft helicity")
def udhel(zstag, mapfct, u, v, wstag, dx, dy, bottom=2000.0, top=5000.0):
    return _udhel(zstag, mapfct, u, v, wstag, dx, dy, bottom, top)


# Requires both u an v for dimnames
@set_alg_metadata(3, "u", units="10-5 s-1",
                  stagdim=-1, stagsubvar="vstag",
                  description="absolute vorticity")
def avo(ustag, vstag, msfu, msfv, msfm, cor, dx, dy):
    return _avo(ustag, vstag, msfu, msfv, msfm, cor, dx, dy)


@set_alg_metadata(3, "t", units="PVU",
                  description="potential vorticity")
def pvo(ustag, vstag, t, p, msfu, msfv, msfm, cor, dx, dy):
    return _pvo(ustag, vstag, t, p, msfu, msfv, msfm, cor, dx, dy)


@set_alg_metadata(3, "qv", units="K",
                  description="equivalent potential temperature")
def eth(qv, tk, p):
    return _eth(qv, tk, p)


@set_alg_metadata(3, "p", units="K",
                  description="wetbulb temperature")
def wetbulb(p, tk, qv):
    return _wetbulb(p, tk, qv)


@set_alg_metadata(3, "tk", units="K",
                  description="virtual temperature")
def tvirtual(tk, qv):
    return _tv(tk, qv)


@set_alg_metadata(3, "qv", units="Pa/s",
                  description="omega")
def omega(qv, tk, w, p):
    return _omega(qv, tk, w, p)


@set_alg_metadata(2, "p", refvarndims=3, units="kg m-2",
                  description="precipitable water")
def pw(p, tk, qv, ht):
    tv = _tv(tk, qv)
    return _pw(full_p, tv, qv, ht)

    
    
