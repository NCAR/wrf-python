from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

from .constants import Constants

from ._wrffortran import (dcomputetk, dinterp3dz, dinterp2dxy, dinterp1d,
                          dcomputeseaprs, dfilter2d, dcomputerh, dcomputeuvmet,
                          dcomputetd, dcapecalc3d, dcloudfrac, wrfcttcalc, 
                          calcdbz, dcalrelhl, dcalcuh, dcomputepv, 
                          dcomputeabsvort, dlltoij, dijtoll, deqthecalc,
                          omgcalc, virtual_temp, wetbulbcalc, dcomputepw,
                          wrf_monotonic, wrf_vintrp, dcomputewspd, 
                          dcomputewdir)

from .decorators import (left_iteration, cast_type, 
                         extract_and_transpose, check_args)
from .util import combine_dims, npbytes_to_str, psafilepath
from .py3compat import py3range
from .specialdec import (uvmet_left_iter, cape_left_iter, 
                         cloudfrac_left_iter, check_cape_args)

class DiagnosticError(Exception):
    """Raised when an error occurs in a diagnostic routine."""
    def __init__(self, message=None):
        """Initialize a :class:`wrf.DiagnosticError` objection.
        
        Args:
        
            message (:obj:`str`): The error message.
        
        """
        self._msg = message
        
    def __str__(self):
        return self._msg
    
    def __call__(self, message):
        """Callable method to make the exception object raise itself.
        
        This allows the exception to be thrown from inside Fortran routines 
        by using f2py's callback mechanism.  This is no longer used within 
        wrf-python, but may be useful to other users.
        
        See Also:
        
            `f2py doc <http://docs.scipy.org/doc/numpy-1.11.0/f2py/>`_
        
        """
        raise self.__class__(message)

# The routines below are thin wrappers around the Fortran functions.  These 
# are not meant to be called by end users.  Use the public API instead for 
# that purpose.

# IMPORTANT!  Unless otherwise noted, all variables used in the routines 
# below assume that Fortran-ordered views are being used.  This allows
# f2py to pass the array pointers directly to the Fortran routine.

@check_args(0, 3, (3, 3, None, None))
@left_iteration(3, 2, ref_var_idx=0, ignore_args=(2,3))
@cast_type(arg_idxs=(0,1))
@extract_and_transpose()
def _interpz3d(field3d, z, desiredloc, missingval, outview=None):
    """Wrapper for dinterp3dz.
    
    Located in wrf_user.f90.
    
    """
    if outview is None:
        outview = np.empty(field3d.shape[0:2], np.float64, order="F")
        
    result = dinterp3dz(field3d, 
                        outview,
                        z, 
                        desiredloc, 
                        missingval)
    return result


@check_args(0, 3, (3,))
@left_iteration(3, combine_dims([(0,-3),(1,-2)]), ref_var_idx=0, 
                  ignore_args=(1,))
@cast_type(arg_idxs=(0,1))
@extract_and_transpose()
def _interp2dxy(field3d, xy, outview=None):
    """Wrapper for dinterp2dxy.
    
    Located in wrf_user.f90.
    
    """
    if outview is None:
        outview = np.empty((xy.shape[-1], field3d.shape[-1]), np.float64, 
                           order="F")
    
    result = dinterp2dxy(field3d,
                         outview,
                         xy)
    return result


@check_args(0, 1, (1,1,None,None))
@left_iteration(1, combine_dims([(2,0)]), ref_var_idx=0, ignore_args=(2,3))
@cast_type(arg_idxs=(0,1,2))
@extract_and_transpose()
def _interp1d(v_in, z_in, z_out, missingval, outview=None):
    """Wrapper for dinterp1d.
    
    Located in wrf_user.f90.
    
    """
    if outview is None:
        outview = np.empty_like(z_out)

    result = dinterp1d(v_in,
                       outview,
                       z_in,
                       z_out,
                       missingval)
    
    return result


@left_iteration(3, combine_dims([(3,0), (1,0)]), 
                  ref_var_idx=0, ignore_args=(1,3,4))
@cast_type(arg_idxs=(0,))
@extract_and_transpose(do_transpose=False)
def _vertcross(field3d, xy, var2dz, z_var2d, missingval, outview=None):
    """Return the vertical cross section.
    
    This routine was originally written in scripted NCL code and doesn't 
    directly wrap a Fortran routine.
    
    Located in WRFUserARW.ncl.
    
    """
    # Note:  This is using C-ordering
    if outview is None:
        outview = np.empty((z_var2d.shape[0], xy.shape[0]), dtype=var2dz.dtype)
        
    var2dtmp = _interp2dxy(field3d, xy)
    
    for i in py3range(xy.shape[0]):
        outview[:,i] = _interp1d(var2dtmp[:,i], var2dz[:,i], z_var2d, 
                                 missingval)
    
    return outview


@left_iteration(2, combine_dims([(1,0)]), ref_var_idx=0, ignore_args=(1,))
@cast_type(arg_idxs=(0,))
@extract_and_transpose(do_transpose=False)
def _interpline(field2d, xy, outview=None):
    """Return the two-dimensional field interpolated to a line.
    
    This routine was originally written in scripted NCL code and doesn't 
    directly wrap a Fortran routine.
    
    Located in WRFUserARW.ncl.
    
    """
    # Note:  This is using C-ordering
    if outview is None:
        outview = np.empty(xy.shape[0], dtype=field2d.dtype)

    tmp_shape = (1,) + field2d.shape
    var2dtmp = np.empty(tmp_shape, field2d.dtype)
    var2dtmp[0,:,:] = field2d[:,:]
    
    var1dtmp = _interp2dxy(var2dtmp, xy)
    
    outview[:] = var1dtmp[0, :]
    
    return outview


@check_args(0, 3, (3,3,3,3))                         
@left_iteration(3, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0,1,2,3))
@extract_and_transpose()
def _slp(z, t, p, q, outview=None):
    """Wrapper for dcomputeseaprs.
    
    Located in wrf_user.f90.
    
    """
    t_surf = np.zeros(z.shape[0:2], np.float64, order="F")
    t_sea_level = np.zeros(z.shape[0:2], np.float64, order="F")
    level = np.zeros(z.shape[0:2], np.int32, order="F")
    
    if outview is None:
        outview = np.empty(z.shape[0:2], np.float64, order="F")
        
    errstat = np.array(0)
    errmsg = np.zeros(Constants.ERRLEN, "c")
    
    result = dcomputeseaprs(z,
                            t,
                            p,
                            q,
                            outview,
                            t_sea_level,
                            t_surf,
                            level,
                            errstat=errstat,
                            errmsg=errmsg)
    
    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())
    
    return result


@check_args(0, 3, (3,3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0,1))
@extract_and_transpose()
def _tk(pressure, theta, outview=None):
    """Wrapper for dcomputetk.
    
    Located in wrf_user.f90.
    
    """
    # No need to transpose here since operations on 1D array
    shape = pressure.shape
    if outview is None: 
        outview = np.empty_like(pressure)
    result = dcomputetk(outview.ravel(order="A"),
                     pressure.ravel(order="A"), 
                     theta.ravel(order="A"))
    result = np.reshape(result, shape, order="F")
    
    return result 


@check_args(0, 2, (2,2))
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0,1))
@extract_and_transpose()
def _td(pressure, qv_in, outview=None):
    """Wrapper for dcomputetd.
    
    Located in wrf_user.f90.
    
    """
    shape = pressure.shape
    if outview is None:
        outview = np.empty_like(pressure)
    
    result = dcomputetd(outview.ravel(order="A"),
                        pressure.ravel(order="A"), 
                        qv_in.ravel(order="A"))
    result = np.reshape(result, shape, order="F")
    
    return result


@check_args(0, 2, (2,2,2))
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0,1,2))
@extract_and_transpose()
def _rh(qv, q, t, outview=None):
    """Wrapper for dcomputerh.
    
    Located in wrf_user.f90.
    
    """
    shape = qv.shape
    if outview is None:
        outview = np.empty_like(qv)
    result = dcomputerh(qv.ravel(order="A"),
                      q.ravel(order="A"),
                      t.ravel(order="A"),
                      outview.ravel(order="A"))
    result = np.reshape(result, shape, order="F")
    
    return result


# Note:  combining the -3 and -2 dimensions from u, then the -1 dimension 
# from v
@check_args(0, 3, (3,3,2,2,2,2), stagger=(-1,-2,-1,-2,None,None),
            refstagdim=-1)
@left_iteration(3, combine_dims([(0, (-3,-2)), 
                                   (1, (-1,))]), 
                  ref_var_idx=0, ignore_args=(6,7))
@cast_type(arg_idxs=(0,1,2,3,4,5))
@extract_and_transpose()
def _avo(u, v, msfu, msfv, msfm, cor, dx, dy, outview=None):
    """Wrapper for dcomputeabsvort.
    
    Located in wrf_pvo.f90.
    
    """
    if outview is None:
        outshape = (v.shape[0],) + u.shape[1:]
        outview = np.empty(outshape, np.float64, order="F")
    
    result = dcomputeabsvort(outview,
                             u,
                             v,
                             msfu,
                             msfv,
                             msfm,
                             cor,
                             dx,
                             dy)
    
    return result


@check_args(0, 3, (3,3,3,3,2,2,2,2), stagger=(-1,-2,None,None,-1,-2,None, 
                                              None),
            refstagdim=-1)
@left_iteration(3, 3, ref_var_idx=2, ignore_args=(8,9))
@cast_type(arg_idxs=(0,1,2,3,4,5,6,7))
@extract_and_transpose()
def _pvo(u, v, theta, prs, msfu, msfv, msfm, cor, dx, dy, outview=None):
    """Wrapper for dcomputepv.
    
    Located in wrf_pvo.f90.
    
    """
    if outview is None:
        outview = np.empty_like(prs)
        
    result = dcomputepv(outview,
                        u,
                        v,
                        theta,
                        prs,
                        msfu,
                        msfv,
                        msfm,
                        cor,
                        dx,
                        dy)
    
    return result


@check_args(0, 3, (3,3,3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0,1,2))
@extract_and_transpose()
def _eth(qv, tk, p, outview=None):
    """Wrapper for deqthecalc.
    
    Located in eqthecalc.f90.
    
    """
    if outview is None:
        outview = np.empty_like(qv)
    
    result = deqthecalc(qv,
                        tk,
                        p,
                        outview)
    
    return result


@uvmet_left_iter()
@cast_type(arg_idxs=(0,1,2,3))
@extract_and_transpose()
def _uvmet(u, v, lat, lon, cen_long, cone, isstag=0, has_missing=False, 
           umissing=Constants.DEFAULT_FILL, vmissing=Constants.DEFAULT_FILL, 
           uvmetmissing=Constants.DEFAULT_FILL, outview=None):
    """Wrapper for dcomputeuvmet.
    
    Located in wrf_user.f90.
    
    """
    longca = np.zeros(lat.shape[0:2], np.float64, order="F")
    longcb = np.zeros(lon.shape[0:2], np.float64, order="F")
    rpd = Constants.PI/180.
        
    if outview is None:
        outdims = u.shape + (2,)
        outview = np.empty(outdims, np.float64, order="F")
    
    result = dcomputeuvmet(u,
                           v,
                           outview,
                           longca,
                           longcb,
                           lon,
                           lat,
                           cen_long,
                           cone,
                           rpd,
                           isstag, 
                           has_missing,
                           umissing,
                           vmissing,
                           uvmetmissing)
    
    return result


@check_args(0, 3, (3,3,3,3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0,1,2,3))
@extract_and_transpose()
def _omega(qv, tk, w, p, outview=None):
    """Wrapper for omgcalc.
    
    Located in wrf_rip_phys_routines.f90.
    
    """
    if outview is None:
        outview = np.empty_like(qv)
    
    result = omgcalc(qv,
                     tk,
                     w,
                     p,
                     outview)
    
    return result


@check_args(0, 3, (3,3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0,1))
@extract_and_transpose()
def _tv(tk, qv, outview=None):
    """Wrapper for virtual_temp.
    
    Located in wrf_rip_phys_routines.f90.
    
    """
    if outview is None:
        outview = np.empty_like(tk)
        
    result = virtual_temp(tk,
                          qv,
                          outview)
    
    return result


@check_args(0, 3, (3,3,3))                         
@left_iteration(3, 3, ref_var_idx=0, ignore_args=(3,))
@cast_type(arg_idxs=(0,1,2))
@extract_and_transpose()
def _wetbulb(p, tk, qv,  psafile=psafilepath(), outview=None):
    """Wrapper for wetbulbcalc.
    
    Located in wrf_rip_phys_routines.f90.
    
    """
    if outview is None:
        outview = np.empty_like(p)
    
    errstat = np.array(0)
    errmsg = np.zeros(Constants.ERRLEN, "c")
    
    result = wetbulbcalc(p,
                         tk,
                         qv,
                         outview,
                         psafile,
                         errstat,
                         errmsg)
    
    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())
    
    return result


@check_args(0, 3, (3,3,3,2))                         
@left_iteration(3, 2, ref_var_idx=0, ignore_args=(4,))
@cast_type(arg_idxs=(0,1,2,3))
@extract_and_transpose()
def _srhel(u, v, z, ter, top, outview=None):
    """Wrapper for dcalrelhl.
    
    Located in wrf_relhl.f90.
    
    """
    if outview is None:
        outview = np.empty_like(ter)
        
    result = dcalrelhl(u, 
                       v, 
                       z, 
                       ter, 
                       top,
                       outview)
    
    return result


@check_args(2, 3, (3,2,3,3,3), stagger=(-3,None,None,None,-3))                         
@left_iteration(3, 2, ref_var_idx=2, ignore_args=(5,6,7,8))
@cast_type(arg_idxs=(0,1,2,3,4))
@extract_and_transpose()
def _udhel(zstag, mapfct, u, v, wstag, dx, dy, bottom, top, outview=None):
    """Wrapper for dcalcuh.
    
    Located in calc_uh.f90.
    
    """
    if outview is None:
        outview = np.empty_like(mapfct)
        
    tem1 = np.zeros((u.shape[0], u.shape[1], u.shape[2]), np.float64, 
                    order="F")
    tem2 = np.zeros((u.shape[0], u.shape[1], u.shape[2]), np.float64, 
                    order="F")
    
    result = dcalcuh(zstag, 
                     mapfct, 
                     dx, 
                     dy, 
                     bottom, 
                     top, 
                     u,
                     v, 
                     wstag, 
                     outview, 
                     tem1, 
                     tem2)
    
    return result


@check_args(0, 3, (3,3,3,3), stagger=(None, None, None, -3))                         
@left_iteration(3, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0,1,2,3))
@extract_and_transpose()
def _pw(p, tv, qv, ht, outview=None):
    """Wrapper for dcomputepw.
    
    Located in wrf_pw.f90.
    
    """
    if outview is None:
        outview = np.empty(p.shape[0:2], p.dtype, order="F")
        
    result = dcomputepw(p,
                        tv,
                        qv,
                        ht,
                        outview)
    
    return result


@check_args(0, 3, (3,3,3,3,3,3))                         
@left_iteration(3, 3, ref_var_idx=0, ignore_args=(6,7,8))
@cast_type(arg_idxs=(0,1,2,3,4,5))
@extract_and_transpose()
def _dbz(p, tk, qv, qr, qs, qg, sn0, ivarint, iliqskin, outview=None):
    """Wrapper for calcdbz.
    
    Located in wrf_user_dbz.f90.
    
    """
    if outview is None:
        outview = np.empty_like(p)
        
    result = calcdbz(p,
                     tk,
                     qv,
                     qr,
                     qs,
                     qg,
                     sn0,
                     ivarint,
                     iliqskin,
                     outview)
    
    return result


@check_cape_args()
@cape_left_iter()
@cast_type(arg_idxs=(0,1,2,3,4,5), outviews=("capeview", "cinview"))
@extract_and_transpose(outviews=("capeview", "cinview"))
def _cape(p_hpa, tk, qv, ht, ter, sfp, missing, i3dflag, ter_follow,
          psafile=psafilepath(), capeview=None, cinview=None):
    """Wrapper for dcapecalc3d.
    
    Located in rip_cape.f90.
    
    """
    if capeview is None:
        capeview = np.zeros(p_hpa.shape[0:3], p_hpa.dtype, order="F")
        
    if cinview is None:
        cinview = np.zeros(p_hpa.shape[0:3], p_hpa.dtype, order="F")
        
    errstat = np.array(0)
    errmsg = np.zeros(Constants.ERRLEN, "c")
    
    # note that p_hpa, tk, qv, and ht have the vertical flipped
    result = dcapecalc3d(p_hpa,
                         tk,
                         qv,
                         ht,
                         ter,
                         sfp,
                         capeview,
                         cinview,
                         missing,
                         i3dflag,
                         ter_follow,
                         psafile,
                         errstat,
                         errmsg)
    
    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())
    
    return result

@check_args(0, 3, (3,3))
@cloudfrac_left_iter()
@cast_type(arg_idxs=(0, 1), outviews=("lowview", "medview", "highview"))
@extract_and_transpose(outviews=("lowview", "medview", "hightview"))
def _cloudfrac(p, rh, lowview=None, medview=None, highview=None):
    """Wrapper for dcloudfrace.
    
    Located in wrf_cloud_fracf.f90.
    
    """
    if lowview is None:
        lowview = np.zeros(p.shape[0:2], p.dtype, order="F")
        
    if medview is None:
        medview = np.zeros(p.shape[0:2], p.dtype, order="F")
        
    if highview is None:
        highview = np.zeros(p.shape[0:2], p.dtype, order="F")
    
    result = dcloudfrac(p, 
                        rh, 
                        lowview, 
                        medview, 
                        highview)
    
    return result


def _lltoxy(map_proj, truelat1, truelat2, stdlon,
            lat1, lon1, pole_lat, pole_lon,
            known_x, known_y, dx, dy, latinc, loninc, lat, lon,
            outview=None):
    """Wrapper for dlltoij.
    
    Located in wrf_user_latlon_routines.f90.
    
    """
    if outview is None:
        outview = np.zeros((2), dtype=np.float64, order="F")
    
    errstat = np.array(0)
    errmsg = np.zeros(Constants.ERRLEN, "c")
    
    result = dlltoij(map_proj,
                     truelat1,
                     truelat2,
                     stdlon,
                     lat1,
                     lon1,
                     pole_lat,
                     pole_lon,
                     known_x,
                     known_y,
                     dx,
                     dy,
                     latinc,
                     loninc,
                     lat,
                     lon,
                     outview,
                     errstat,
                     errmsg)
    
    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())
    
    return result


def _xytoll(map_proj, truelat1, truelat2, stdlon, lat1, lon1,
             pole_lat, pole_lon, known_x, known_y, dx, dy, latinc,
             loninc, x, y, outview=None):
    """Wrapper for dijtoll.
    
    Located in wrf_user_latlon_routines.f90.
    
    """
    if outview is None:
        outview = np.zeros((2), dtype=np.float64, order="F")
    
    errstat = np.array(0)
    errmsg = np.zeros(Constants.ERRLEN, "c")
    
    result = dijtoll(map_proj,
                     truelat1,
                     truelat2,
                     stdlon,
                     lat1,
                     lon1,
                     pole_lat,
                     pole_lon,
                     known_x,
                     known_y,
                     dx,
                     dy,
                     latinc,
                     loninc,
                     x,
                     y,
                     outview,
                     errstat,
                     errmsg)
    
    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())
    
    return result


@check_args(0, 3, (3,3,3,3,3,3,2))                         
@left_iteration(3, 2, ref_var_idx=0, ignore_args=(7,))
@cast_type(arg_idxs=(0,1,2,3,4,5,6))
@extract_and_transpose()
def _ctt(p_hpa, tk, qice, qcld, qv, ght, ter, haveqci, outview=None):
    """Wrapper for wrfcttcalc.
    
    Located in wrf_fctt.f90.
    
    """
    if outview is None:
        outview = np.empty_like(ter)
        
    result = wrfcttcalc(p_hpa,
                        tk,
                        qice,
                        qcld,
                        qv,
                        ght,
                        ter,
                        outview,
                        haveqci)
    
    return result


@check_args(0, 2, (2,))
@left_iteration(2, 2, ref_var_idx=0, ignore_args=(1,))
@cast_type(arg_idxs=(0,))
@extract_and_transpose()
def _smooth2d(field, passes, outview=None):
    """Wrapper for dfilter2d.
    
    Located in wrf_user.f90.
    
    """
    # Unlike NCL, this routine will not modify the values in place, but 
    # copies the original data before modifying it.
    
    if isinstance(field, np.ma.MaskedArray):
        missing = field.fill_value
    else:
        missing = Constants.DEFAULT_FILL
    
    if outview is None:
        outview = field.copy(order="A")
    else:
        outview[:] = field[:]
        
    field_tmp = np.zeros(outview.shape, outview.dtype, order="F") 

    dfilter2d(outview, 
              field_tmp,               
              passes,
              missing)
    
    return outview


@check_args(0, 3, (3,3,2))    
@left_iteration(3, 3, ref_var_idx=0, ignore_args=(3,4,5))
@cast_type(arg_idxs=(0,1,2))
@extract_and_transpose()
def _monotonic(var, lvprs, coriolis, idir, delta, icorsw, outview=None):
    """Wrapper for wrf_monotonic.
    
    Located in wrf_vinterp.f90.
    
    """
    # If icorsw is not 0, then the input variable might get modified by the 
    # fortran routine.  We don't want this, so make a copy and pass that on.
    var = var.copy(order="A") if icorsw != 0 else var
        
    if outview is None:
        outview = np.empty_like(var)
    
    result = wrf_monotonic(outview,
                           var,
                           lvprs,
                           coriolis,
                           idir,
                           delta,
                           icorsw)
    
    return result


# Output shape is interp_levels.shape + field.shape[-2:]
@check_args(0, 3, (3,3,3,3,3,2,2,2,3))    
@left_iteration(3, combine_dims([(9, (-1,)),
                                   (0, (-2,-1))]), 
                  ref_var_idx=0, ignore_args=(9,10,11,12,13,14))
@cast_type(arg_idxs=(0,1,2,3,4,5,6,7,8,9))
@extract_and_transpose()
def _vintrp(field, pres, tk, qvp, ght, terrain, sfp, smsfp,
            vcarray, interp_levels, icase, extrap, vcor, logp,
            missing, outview=None):
    """Wrapper for wrf_vintrp.
    
    Located in wrf_vinterp.f90.
    
    """
    if outview is None:
        outdims = field.shape[0:2] + interp_levels.shape
        outview = np.empty(outdims, field.dtype, order="F")
        
    errstat = np.array(0)
    errmsg = np.zeros(Constants.ERRLEN, "c")
    
    result = wrf_vintrp(field,
                        outview,
                        pres,
                        tk,
                        qvp,
                        ght,
                        terrain,
                        sfp,
                        smsfp,
                        vcarray,
                        interp_levels,
                        icase,
                        extrap,
                        vcor,
                        logp,
                        missing,
                        errstat,
                        errmsg)
    
    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())
    
    return result

@check_args(0, 2, (2,2))                         
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0,1))
@extract_and_transpose()
def _wspd(u, v, outview=None):
    """Wrapper for dcomputewspd.
    
    Located in wrf_wind.f90.
    
    """
    if outview is None:
        outview = np.empty_like(u)
        
    result = dcomputewspd(outview,
                          u,
                          v)
    
    return result


@check_args(0, 2, (2,2))                         
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0,1))
@extract_and_transpose()
def _wdir(u, v, outview=None):
    """Wrapper for dcomputewdir.
    
    Located in wrf_wind.f90.
    
    """
    if outview is None:
        outview = np.empty_like(u)
        
    result = dcomputewdir(outview,
                          u,
                          v)
    
    return result
    

