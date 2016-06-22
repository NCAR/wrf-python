
from .constants import Constants
from .extension import (_interpz3d, _interp2dxy, _interp1d, _slp, _tk, _td, 
                        _rh, _uvmet, _smooth2d)
from .util import from_var
from .metadecorators import (set_alg_metadata, set_uvmet_alg_metadata, 
                             set_interp_metadata)
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


@set_alg_metadata(2, refvarname="z", refvarndims=3, units="hpa",
                  description="sea level pressure")
def slp(z, t, p, q, meta=True):
    return _slp(z, t, p, q)


@set_alg_metadata(3, refvarname="pressure", units="K",
                  description="temperature")
def tk(pressure, theta, meta=True):
    return _tk(pressure, theta)


@set_alg_metadata(3, refvarname="pressure", units="degC",
                  description="dew point temperature")
def td(pressure, qv_in, meta=True):
    return _td(pressure, qv_in)


@set_alg_metadata(3, refvarname="pressure", 
                  description="relative humidity", units=None)
def rh(qv, q, t, meta=True):
    return _rh(qv, q, t, meta)


@set_uvmet_alg_metadata()
def uvmet(u, v, lat, lon, cen_long, cone, meta=True):
    return _uvmet(u, v, lat, lon, cen_long, cone)


@set_alg_metadata(2,
                  refvarname="field", 
                  description=from_var("field", "description"), 
                  units=from_var("field", "units"))
def smooth2d(field, passes, meta=True):
    return _smooth2d(field, passes)
