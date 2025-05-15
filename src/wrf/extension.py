from __future__ import (absolute_import, division, print_function)

import numpy as np

from .constants import Constants, default_fill

from wrf._wrffortran import (dcomputetk, dinterp3dz, dinterp2dxy, dinterp1d,
                             dcomputeseaprs, dfilter2d, dcomputerh,
                             dcomputeuvmet, dcomputetd, dcapecalc2d,
                             dcapecalc3d, dcloudfrac2, wrfcttcalc, calcdbz,
                             dcalrelhl, dcalcuh, dcomputepv, dcomputeabsvort,
                             dlltoij, dijtoll, deqthecalc, omgcalc,
                             virtual_temp, wetbulbcalc, dcomputepw,
                             wrf_monotonic, wrf_vintrp, dcomputewspd,
                             dcomputewdir, dinterp3dz_2dlev,
                             fomp_set_num_threads, fomp_get_num_threads,
                             fomp_get_max_threads, fomp_get_thread_num,
                             fomp_get_num_procs, fomp_in_parallel,
                             fomp_set_dynamic, fomp_get_dynamic,
                             fomp_set_nested, fomp_get_nested,
                             fomp_set_schedule, fomp_get_schedule,
                             fomp_get_thread_limit, fomp_set_max_active_levels,
                             fomp_get_max_active_levels, fomp_get_level,
                             fomp_get_ancestor_thread_num, fomp_get_team_size,
                             fomp_get_active_level, fomp_in_final,
                             fomp_init_lock, fomp_init_nest_lock,
                             fomp_destroy_lock, fomp_destroy_nest_lock,
                             fomp_set_lock, fomp_set_nest_lock,
                             fomp_unset_lock, fomp_unset_nest_lock,
                             fomp_test_lock, fomp_test_nest_lock,
                             fomp_get_wtime, fomp_get_wtick, fomp_enabled)

from .decorators import (left_iteration, cast_type,
                         extract_and_transpose, check_args)
from .util import combine_dims, npbytes_to_str, psafilepath
from .py3compat import py3range
from .specialdec import (uvmet_left_iter, cape_left_iter,
                         cloudfrac_left_iter, check_cape_args,
                         interplevel_left_iter, check_interplevel_args)


class DiagnosticError(Exception):
    """Raised when an error occurs in a diagnostic routine."""
    def __init__(self, message=None):
        """Initialize a :class:`wrf.DiagnosticError` object.

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
@check_interplevel_args(is2dlev=False)
@interplevel_left_iter(is2dlev=False)
@cast_type(arg_idxs=(0, 1, 2))
@extract_and_transpose()
def _interpz3d(field3d, z, desiredloc, missingval, outview=None):
    """Wrapper for dinterp3dz.

    Located in wrf_user.f90.

    """
    if outview is None:
        outshape = field3d.shape[0:2] + desiredloc.shape
        outview = np.empty(outshape, np.float64, order="F")

    result = dinterp3dz(field3d,
                        outview,
                        z,
                        desiredloc,
                        missingval)
    return result


@check_interplevel_args(is2dlev=True)
@interplevel_left_iter(is2dlev=True)
@cast_type(arg_idxs=(0, 1, 2))
@extract_and_transpose()
def _interpz3d_lev2d(field3d, z, lev2d, missingval, outview=None):
    """Wrapper for dinterp3dz.

    Located in wrf_user.f90.

    """
    if outview is None:
        outview = np.empty(field3d.shape[0:2], np.float64, order="F")

    result = dinterp3dz_2dlev(field3d,
                              outview,
                              z,
                              lev2d,
                              missingval)
    return result


@check_args(0, 3, (3, ))
@left_iteration(3, combine_dims([(0, -3), (1, -2)]), ref_var_idx=0,
                ignore_args=(1, ))
@cast_type(arg_idxs=(0, 1))
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


@check_args(0, 1, (1, 1, None, None))
@left_iteration(1, combine_dims([(2, 0)]), ref_var_idx=0, ignore_args=(2, 3))
@cast_type(arg_idxs=(0, 1, 2))
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


@left_iteration(3, combine_dims([(3, 0), (1, 0)]),
                ref_var_idx=0, ignore_args=(1, 3, 4))
@cast_type(arg_idxs=(0, ))
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
        outview[:, i] = _interp1d(var2dtmp[:, i], var2dz[:, i], z_var2d,
                                  missingval)

    return outview


@left_iteration(2, combine_dims([(1, 0)]), ref_var_idx=0, ignore_args=(1, ))
@cast_type(arg_idxs=(0, ))
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
    var2dtmp[0, :, :] = field2d[:, :]

    var1dtmp = _interp2dxy(var2dtmp, xy)

    outview[:] = var1dtmp[0, :]

    return outview


@check_args(0, 3, (3, 3, 3, 3))
@left_iteration(3, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1, 2, 3))
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


@check_args(0, 3, (3, 3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1))
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


@check_args(0, 2, (2, 2))
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1))
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


@check_args(0, 2, (2, 2, 2))
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1, 2))
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
@check_args(0, 3, (3, 3, 2, 2, 2, 2), stagger=(-1, -2, -1, -2, None, None),
            refstagdim=-1)
@left_iteration(3, combine_dims([(0, (-3, -2)),
                                (1, (-1, ))]),
                ref_var_idx=0, ignore_args=(6, 7))
@cast_type(arg_idxs=(0, 1, 2, 3, 4, 5))
@extract_and_transpose()
def _avo(u, v, msfu, msfv, msfm, cor, dx, dy, outview=None):
    """Wrapper for dcomputeabsvort.

    Located in wrf_pvo.f90.

    """
    if outview is None:
        outshape = (v.shape[0], ) + u.shape[1:]
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


@check_args(0, 3, (3, 3, 3, 3, 2, 2, 2, 2),
            stagger=(-1, -2, None, None, -1, -2, None, None), refstagdim=-1)
@left_iteration(3, 3, ref_var_idx=2, ignore_args=(8, 9))
@cast_type(arg_idxs=(0, 1, 2, 3, 4, 5, 6, 7))
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


@check_args(0, 3, (3, 3, 3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1, 2))
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
@cast_type(arg_idxs=(0, 1, 2, 3))
@extract_and_transpose()
def _uvmet(u, v, lat, lon, cen_long, cone, isstag=0, has_missing=False,
           umissing=default_fill(np.float64),
           vmissing=default_fill(np.float64),
           uvmetmissing=default_fill(np.float64),
           outview=None):
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


@check_args(0, 3, (3, 3, 3, 3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1, 2, 3))
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


@check_args(0, 3, (3, 3))
@left_iteration(3, 3, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1))
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


@check_args(0, 3, (3, 3, 3))
@left_iteration(3, 3, ref_var_idx=0, ignore_args=(3,))
@cast_type(arg_idxs=(0, 1, 2))
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


@check_args(0, 3, (3, 3, 3, 2, 2))
@left_iteration(3, 2, ref_var_idx=0, ignore_args=(5, ))
@cast_type(arg_idxs=(0, 1, 2, 3, 4))
@extract_and_transpose()
def _srhel(u, v, z, ter, lats, top, outview=None):
    """Wrapper for dcalrelhl.

    Located in wrf_relhl.f90.

    """
    if outview is None:
        outview = np.empty_like(ter)

    result = dcalrelhl(u,
                       v,
                       z,
                       ter,
                       lats,
                       top,
                       outview)

    return result


@check_args(2, 3, (3, 2, 3, 3, 3), stagger=(-3, None, None, None, -3))
@left_iteration(3, 2, ref_var_idx=2, ignore_args=(5, 6, 7, 8))
@cast_type(arg_idxs=(0, 1, 2, 3, 4))
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


@check_args(0, 3, (3, 3, 3, 3), stagger=(None, None, None, -3))
@left_iteration(3, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1, 2, 3))
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


@check_args(0, 3, (3, 3, 3, 3, 3, 3))
@left_iteration(3, 3, ref_var_idx=0, ignore_args=(6, 7, 8))
@cast_type(arg_idxs=(0, 1, 2, 3, 4, 5))
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
@cast_type(arg_idxs=(0, 1, 2, 3, 4, 5), outviews=("capeview", "cinview"))
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

    if i3dflag:
        cape_routine = dcapecalc3d
    else:
        cape_routine = dcapecalc2d

    # Work arrays
    k_left_shape = (p_hpa.shape[2], p_hpa.shape[0], p_hpa.shape[1])
    prsf = np.empty(k_left_shape, np.float64, order="F")
    prs_new = np.empty(k_left_shape, np.float64, order="F")
    tmk_new = np.empty(k_left_shape, np.float64, order="F")
    qvp_new = np.empty(k_left_shape, np.float64, order="F")
    ght_new = np.empty(k_left_shape, np.float64, order="F")

    # note that p_hpa, tk, qv, and ht have the vertical flipped
    result = cape_routine(p_hpa,
                          tk,
                          qv,
                          ht,
                          ter,
                          sfp,
                          capeview,
                          cinview,
                          prsf,
                          prs_new,
                          tmk_new,
                          qvp_new,
                          ght_new,
                          missing,
                          ter_follow,
                          psafile,
                          errstat,
                          errmsg)

    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())

    return result


@check_args(0, 3, (3, 3))
@cloudfrac_left_iter()
@cast_type(arg_idxs=(0, 1), outviews=("lowview", "midview", "highview"))
@extract_and_transpose(outviews=("lowview", "midview", "highview"))
def _cloudfrac(vert, rh, vert_inc_w_height, low_thresh, mid_thresh,
               high_thresh, missing, lowview=None, midview=None,
               highview=None):
    """Wrapper for dcloudfrac2.

    Located in wrf_cloud_fracf.f90.

    """
    if lowview is None:
        lowview = np.zeros(vert.shape[0:2], vert.dtype, order="F")

    if midview is None:
        midview = np.zeros(vert.shape[0:2], vert.dtype, order="F")

    if highview is None:
        highview = np.zeros(vert.shape[0:2], vert.dtype, order="F")

    result = dcloudfrac2(vert,
                         rh,
                         vert_inc_w_height,
                         low_thresh,
                         mid_thresh,
                         high_thresh,
                         missing,
                         lowview,
                         midview,
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


@check_args(0, 3, (3, 3, 3, 3, 3, 3, 2))
@left_iteration(3, 2, ref_var_idx=0, ignore_args=(7, 8, 9, 10))
@cast_type(arg_idxs=(0, 1, 2, 3, 4, 5, 6))
@extract_and_transpose()
def _ctt(p_hpa, tk, qice, qcld, qv, ght, ter, haveqci, fill_nocloud,
         missing, opt_thresh, outview=None):
    """Wrapper for wrfcttcalc.

    Located in wrf_fctt.f90.

    """
    if outview is None:
        outview = np.empty_like(ter)

    pf = np.empty(p_hpa.shape[0:3], np.float64, order="F")

    result = wrfcttcalc(p_hpa,
                        tk,
                        qice,
                        qcld,
                        qv,
                        ght,
                        ter,
                        outview,
                        pf,
                        haveqci,
                        fill_nocloud,
                        missing,
                        opt_thresh)

    return result


@check_args(0, 2, (2, ))
@left_iteration(2, 2, ref_var_idx=0, ignore_args=(1, 2))
@cast_type(arg_idxs=(0, ))
@extract_and_transpose()
def _smooth2d(field, passes, cenweight, outview=None):
    """Wrapper for dfilter2d.

    Located in wrf_user.f90.

    """
    # Unlike NCL, this routine will not modify the values in place, but
    # copies the original data before modifying it.

    if isinstance(field, np.ma.MaskedArray):
        missing = field.fill_value
    else:
        missing = default_fill(np.float64)

    if outview is None:
        outview = field.copy(order="A")
    else:
        outview[:] = field[:]

    field_tmp = np.zeros(outview.shape, outview.dtype, order="F")

    dfilter2d(outview,
              field_tmp,
              passes,
              missing,
              cenweight)

    return outview


@check_args(0, 3, (3, 3, 2))
@left_iteration(3, 3, ref_var_idx=0, ignore_args=(3, 4, 5))
@cast_type(arg_idxs=(0, 1, 2))
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
@check_args(0, 3, (3, 3, 3, 3, 3, 2, 2, 2, 3))
@left_iteration(3, combine_dims([(9, (-1, )),
                                 (0, (-2, -1))]),
                ref_var_idx=0, ignore_args=(9, 10, 11, 12, 13, 14))
@cast_type(arg_idxs=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
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

    tempout = np.zeros(field.shape[0:2], np.float64, order="F")

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
                        tempout,
                        missing,
                        errstat,
                        errmsg)

    if int(errstat) != 0:
        raise DiagnosticError("".join(npbytes_to_str(errmsg)).strip())

    return result


@check_args(0, 2, (2, 2))
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1))
@extract_and_transpose()
def _wspd(u, v, outview=None):
    """Wrapper for dcomputewspd.

    Located in wrf_wind.f90.

    """
    shape = u.shape
    if outview is None:
        outview = np.empty_like(u)

    result = dcomputewspd(outview.ravel(order="A"),
                          u.ravel(order="A"),
                          v.ravel(order="A"))

    result = np.reshape(result, shape, order="F")

    return result


@check_args(0, 2, (2, 2))
@left_iteration(2, 2, ref_var_idx=0)
@cast_type(arg_idxs=(0, 1))
@extract_and_transpose()
def _wdir(u, v, outview=None):
    """Wrapper for dcomputewdir.

    Located in wrf_wind.f90.

    """
    shape = u.shape
    if outview is None:
        outview = np.empty_like(u)

    result = dcomputewdir(outview.ravel(order="A"),
                          u.ravel(order="A"),
                          v.ravel(order="A"))

    result = np.reshape(result, shape, order="F")

    return result


# OpenMP runtime wrappers

def omp_set_num_threads(num_threads):
    """Specify the number of threads to use.

    The omp_set_num_threads routine affects the number of threads to be used
    for subsequent parallel regions that do not specify a num_threads
    clause, by setting the value of the first element of the nthreads-var
    ICV of the current task.

    Args:

        num_threads (a positive :obj:`int`): The number of threads.  Must be
            positive.

    Returns:

        None.

    """
    if num_threads < 0 and fomp_enabled():
        raise ValueError("'num_threads' must be a positive integer.")

    fomp_set_num_threads(num_threads)


def omp_get_num_threads():
    """Return the number of threads in the current team.

    The omp_get_num_threads routine returns the number of threads in the
    team executing the parallel region to which the routine region binds.
    If called from the sequential part of a program, this routine returns 1.

    Note:

        This function always returns 1 when called from within Python.

    Returns:

        :obj:`int`: The number of threads in the current team.

    See Also:

        :meth:`wrf.omp_get_max_threads`, :meth:`wrf.omp_set_num_threads`

    """
    return fomp_get_num_threads()


def omp_get_max_threads():
    """Return the maximum number of threads that can be used in a parallel \
    region.

    The omp_get_max_threads routine returns an upper bound on the number of
    threads that could be used to form a new team if a parallel construct
    without a num_threads clause were encountered after execution returns from
    this routine.

    Returns:

        :obj:`int`: The number of threads in the current team.

    See Also:

        :meth:`wrf.omp_set_num_threads`

    """
    return fomp_get_max_threads()


def omp_get_thread_num():
    """Return the thread number, within the current team, of the \
    calling thread.

    The omp_get_thread_num routine returns the thread number of the calling
    thread, within the team executing the parallel region to which the routine
    region binds. The thread number is an integer between 0 and one less than
    the value returned by omp_get_num_threads, inclusive. The thread number of
    the master thread of the team is 0. The routine returns 0 if it is called
    from the sequential part of a program.

    Note:

        This function always returns 0 when called from within Python.

    Returns:

        :obj:`int`: The thread number.

    See Also:

        :meth:`wrf.omp_get_num_procs`

    """
    return fomp_get_thread_num()


def omp_get_num_procs():
    """Return the number of processors on the device.

    The omp_get_num_procs routine returns the number of processors that are
    available to the device at the time the routine is called. This value may
    change between the time that it is determined by the omp_get_num_procs
    routine and the time that it is read in the calling context due to system
    actions outside the control of the OpenMP implementation.

    Returns:

        :obj:`int`: The number of processors.

    """
    return fomp_get_num_procs()


def omp_in_parallel():
    """Return 1 if the active-levels-var ICV is greater than zero; \
    otherwise, return 0.

    The effect of the omp_in_parallel routine is to return 1 if the current
    task is enclosed by an active parallel region, and the parallel region is
    enclosed by the outermost initial task region on the device; otherwise it
    returns 0.

    Note:

        This function always returns 0 when called from within Python.

    Returns:

        :obj:`int`: Returns 1 if the active-levels-var ICV is greater than
        zero.  Otherwise, it returns 0.

    """
    return fomp_in_parallel()


def omp_set_dynamic(dynamic_threads):
    """Enable or disable dynamic adjustment of the number of threads \
    available for the execution of subsequent parallel regions by setting the \
    value of the dyn-var ICV.

    For implementations that support dynamic adjustment of the number of
    threads, if the argument to omp_set_dynamic evaluates to True, dynamic
    adjustment is enabled for the current task; otherwise, dynamic adjustment
    is disabled for the current task. For implementations that do not support
    dynamic adjustment of the number of threads this routine has no effect:
    the value of dyn-var remains false.

    Args:

        dynamic_threads (:obj:`bool`): Set to True to support the dynamic
            adjustment of the number of threads.  Otherwise, set to False.

    Returns:

        None.

    See Also:

        :meth:`wrf.omp_get_dynamic`

    """
    fomp_set_dynamic(dynamic_threads)


def omp_get_dynamic():
    """Return the value of the dyn-var ICV, which determines whether \
    dynamic adjustment of the number of threads is enabled or disabled.

    This routine returns 1 if dynamic adjustment of the number of threads
    is enabled for the current task; it returns 0, otherwise. If an
    implementation does not support dynamic adjustment of the
    number of threads, then this routine always returns 0.

    Returns:

        :obj:`int`: Returns 1 if dynamic thread adjustment is enabled, 0
        if disabled.

    See Also:

        :meth:`wrf.omp_set_dynamic`

    """
    return fomp_get_dynamic()


def omp_set_nested(nested):
    """Enable or disable nested parallelism, by setting the nest-var ICV

    For implementations that support nested parallelism, if the argument to
    omp_set_nested evaluates to True, nested parallelism is enabled for the
    current task; otherwise, nested parallelism is disabled for the current
    task. For implementations that do not support nested parallelism, this
    routine has no effect: the value of nest-var remains False.

    Args:

        dynamic_threads (:obj:`bool`): Set to True to support nested
            parallelism, otherwise False.

    Returns:

        None

    See Also:

        :meth:`wrf.omp_get_nested`

    """
    fomp_set_nested(nested)


def omp_get_nested():
    """Return the value of the nest-var ICV, which determines if nested \
    parallelism is enabled or disabled

    This routine returns 1 if nested parallelism is enabled for the current
    task; it returns 0, otherwise. If an implementation does not support
    nested parallelism, this routine always returns 0.

    Returns:

        :obj:`int`: Returns 1 if nested parallelism is enabled, otherwise 0.

    See Also:

        :meth:`wrf.omp_set_nested`

    """
    return fomp_get_nested()


def omp_set_schedule(kind, modifier=0):
    """Set the schedule that is applied when *runtime* is used as \
    schedule kind, by setting the value of the run-sched-var ICV.

    The effect of this routine is to set the value of the run-sched-var ICV
    of the current task to the values specified in the two arguments. The
    schedule is set to the schedule type specified by the first argument kind.
    It can be any of the standard schedule types or any other implementation
    specific one. For the schedule types static, dynamic, and guided the
    chunk_size is set to the value of the second argument, or to the default
    chunk_size if the value of the second argument is less than 1; for the
    schedule type auto the second argument has no meaning; for implementation
    specific schedule types, the values and associated meanings of the second
    argument are implementation defined.

    Args:

        kind (:obj:`int`):  Must be :data:`wrf.OMP_SCHED_STATIC`,
            :data:`wrf.OMP_SCHED_DYNAMIC`, :data:`wrf.OMP_SCHED_GUIDED`,
            or :data:`wrf.OMP_SCHED_AUTO`.

        modifier(:obj:`int`): An implementation specific value, depending on
            the choice for *kind*. This parameter is alternatively named
            chunk_size in some OpenMP documentation. Default is 0, which
            means the OpenMP implementation will use its default value.

    Returns:

        None

    See Also:

        :meth:`wrf.omp_get_schedule`

    """
    fomp_set_schedule(kind, modifier)


def omp_get_schedule():
    """Return the schedule that is applied when the runtime schedule is used.

    This routine returns the run-sched-var ICV in the task to which the routine
    binds. The first item is the schedule kind, which will be one of
    :data:`wrf.OMP_SCHED_STATIC`, :data:`wrf.OMP_SCHED_DYNAMIC`,
    :data:`wrf.OMP_SCHED_GUIDED`, or :data:`wrf.OMP_SCHED_AUTO`.  The second
    item returned is the modifier, which is often named chunk_size in
    OpenMP documentation.

    Returns:

        :obj:`tuple`: The first item is an :obj:`int` for the schedule *kind*.
        The second items is :obj:`int` for the *modifier* (chunk_size).

    See Also:

        :meth:`wrf.omp_set_schedule`

    """
    return fomp_get_schedule()


def omp_get_thread_limit():
    """Return the maximum number of OpenMP threads available to participate \
    in the current contention group.

    The omp_get_thread_limit routine returns the value of the thread-limit-var
    ICV.

    Returns:

        :obj:`int`: The number of OpenMP threads available to participate in
        the current contention group.

    See Also:

        :meth:`wrf.omp_get_max_threads`

    """
    return fomp_get_thread_limit()


def omp_set_max_active_levels(max_levels):
    """Limit the number of nested active parallel regions on the device, \
    by setting the max-active-levels-var ICV.

    The effect of this routine is to set the value of the max-active-levels-var
    ICV to the value specified in the argument. If the number of parallel
    levels requested exceeds the number of levels of parallelism supported by
    the implementation, the value of the max-active-levels-var ICV will be set
    to the number of parallel levels supported by the implementation. This
    routine has the described effect only when called from a sequential part
    of the program. When called from within an explicit parallel region, the
    effect of this routine is implementation defined.

    Args:

        max_levels (:obj:`int`):  The maximum number of nested active parallel
            regions.

    Returns:

        None.

    See Also:

        :meth:`wrf.omp_get_max_active_levels`

    """
    fomp_set_max_active_levels(max_levels)


def omp_get_max_active_levels():
    """Return the value of the max-active-levels-var ICV, which determines \
    the maximum number of nested active parallel regions on the device

    The omp_get_max_active_levels routine returns the value of the
    max-active-levels-var ICV, which determines the maximum number of nested
    active parallel regions on the device.

    Returns:

        :obj:`int`: The maximum number of nested active parallel regions.

    See Also:

        :meth:`wrf.omp_set_max_active_levels`

    """
    return fomp_get_max_active_levels()


def omp_get_level():
    """Return the value of the levels-var ICV.

    The effect of the omp_get_level routine is to return the number of nested
    parallel regions (whether active or inactive) enclosing the current task
    such that all of the parallel regions are enclosed by the outermost initial
    task region on the current device.

    Returns:

        :obj:`int`: The number of nested parallel regions.

    See Also:

        :meth:`wrf.omp_get_active_level`

    """
    return fomp_get_level()


def omp_get_ancestor_thread_num(level):
    """Return, for a given nested level of the current thread, the thread \
    number of the ancestor of the current thread.

    The omp_get_ancestor_thread_num routine returns the thread number of the
    ancestor at a given nest level of the current thread or the thread number
    of the current thread. If the requested nest level is outside the range of
    0 and the nest level of the current thread, as returned by the
    omp_get_level routine, the routine returns -1.

    Args:

        level (:obj:`int`):  The nested level of the current thread.

    Returns:

        :obj:`int`: The thread number of the ancestor at a given nest level
            of the current thread.

    See Also:

        :meth:`wrf.omp_get_max_active_levels`, :meth:`wrf.omp_get_level`

    """
    return fomp_get_ancestor_thread_num(level)


def omp_get_team_size(level):
    """Return, for a given nested level of the current thread, the size \
    of the thread team to which the ancestor or the current thread belongs

    The omp_get_team_size routine returns the size of the thread team to which
    the ancestor or the current thread belongs. If the requested nested level
    is outside the range of 0 and the nested level of the current thread, as
    returned by the omp_get_level routine, the routine returns -1. Inactive
    parallel regions are regarded like active parallel regions executed with
    one thread.

    Args:

        level (:obj:`int`):  The nested level of the current thread.

    Returns:

        :obj:`int`: The size of the thread team.

    See Also:

        :meth:`wrf.omp_get_ancestor_thread_num`

    """
    return fomp_get_team_size(level)


def omp_get_active_level():
    """Return the value of the active-level-vars ICV.

    The effect of the omp_get_active_level routine is to return the number of
    nested, active parallel regions enclosing the current task such that all
    of the parallel regions are enclosed by the outermost initial task region
    on the current device.

    Returns:

        :obj:`int`: The number of nested activate parallel regions.

    See Also:

        :meth:`wrf.omp_get_team_size`

    """
    return fomp_get_active_level()


def omp_in_final():
    """Return 1 (True) if the routine is executed in a final task region; \
    otherwise, it returns 0 (False).

    Returns:

        :obj:`int`: Return 1 if the routine is executed in a final task
        region, 0 otherwise.

    """
    return fomp_in_final()


def omp_init_lock():
    """Initialize a simple OpenMP lock.

    Returns:

        :obj:`int`:  An integer representing the lock.

    See Also:

        :meth:`wrf.omp_init_nest_lock`, :meth:`wrf.omp_destroy_lock`

    """
    return fomp_init_lock()


def omp_init_nest_lock():
    """Initialize a nestable OpenMP lock.

    Returns:

        :obj:`int`:  An integer representing the nestable lock.

    See Also:

        :meth:`wrf.omp_init_lock`

    """
    return fomp_init_nest_lock()


def omp_destroy_lock(svar):
    """Destroy a simple OpenMP lock.

    This sets the lock to an uninitialized state.

    Args:

        svar (:obj:`int`):  An integer representing the lock.

    See Also:

        :meth:`wrf.omp_destroy_nest_lock`, :meth:`wrf.omp_init_lock`

    """
    fomp_destroy_lock(svar)


def omp_destroy_nest_lock(nvar):
    """Destroy a nestable OpenMP lock.

    This sets the lock to an uninitialized state.

    Args:

        nvar (:obj:`int`):  An integer representing the nestable lock.

    See Also:

        :meth:`wrf.omp_destroy_lock`, :meth:`wrf.omp_init_nest_lock`

    """
    fomp_destroy_nest_lock(nvar)


def omp_set_lock(svar):
    """Set a simple OpenMP lock.

    Args:

        svar (:obj:`int`):  An integer representing the lock.

    See Also:

        :meth:`wrf.omp_unset_lock`, :meth:`wrf.omp_set_nest_lock`

    """
    fomp_set_lock(svar)


def omp_set_nest_lock(nvar):
    """Set a nestable OpenMP lock.

    Args:

        nvar (:obj:`int`):  An integer representing the nestable lock.

    See Also:

        :meth:`wrf.omp_unset_nest_lock`, :meth:`wrf.omp_set_lock`

    """
    fomp_set_nest_lock(nvar)


def omp_unset_lock(svar):
    """Unset a simple OpenMP lock.

    Args:

        svar (:obj:`int`):  An integer representing the simple lock.

    See Also:

        :meth:`wrf.omp_unset_nest_lock`, :meth:`wrf.omp_set_lock`

    """
    fomp_unset_lock(svar)


def omp_unset_nest_lock(nvar):
    """Unset a nestable OpenMP lock.

    Args:

        nvar (:obj:`int`):  An integer representing the nestable lock.

    See Also:

        :meth:`wrf.omp_set_nest_lock`, :meth:`wrf.omp_unset_lock`

    """
    fomp_unset_nest_lock(nvar)


def omp_test_lock(svar):
    """Test a simple OpenMP lock.

    This method attempts to set the lock, but does not suspend execution.

    Args:

        svar (:obj:`int`):  An integer representing the simple lock.

    Returns:

        :obj:`int`:  Returns 1 (True) if the lock is successfully set,
        otherwise 0 (False).

    See Also:

        :meth:`wrf.test_nest_lock`, :meth:`wrf.omp_set_lock`

    """
    return fomp_test_lock(svar)


def omp_test_nest_lock(nvar):
    """Test a nestable OpenMP lock.

    This method attempts to set the lock, but does not suspend execution.

    Args:

        nvar (:obj:`int`):  An integer representing the simple lock.

    Returns:

        :obj:`int`:  Returns the nesting count if successful,
        otherwise 0 (False).

    See Also:

        :meth:`wrf.test_lock`, :meth:`wrf.omp_set_nest_lock`

    """
    return fomp_test_nest_lock(nvar)


def omp_get_wtime():
    """Return elapsed wall clock time in seconds.

    The omp_get_wtime routine returns a value equal to the elapsed wall clock
    time in seconds since some "time in the past". The actual
    "time in the past" is arbitrary, but it is guaranteed not to change during
    the execution of the application program. The time returned is a
    "per-thread time", so it is not required to be globally consistent across
    all threads participating in an application.

    Returns:

        :obj:`float`:  Returns the wall clock time in seconds.

    See Also:

        :meth:`wrf.omp_get_wtick`

    """
    return fomp_get_wtime()


def omp_get_wtick():
    """Return the precision of the timer used by :meth:`wrf.omp_get_wtime`.

    The omp_get_wtick routine returns a value equal to the number of
    seconds between successive clock ticks of the timer used by
    :meth:`wrf.omp_get_wtime`.

    Returns:

        :obj:`float`:  Returns the precision of the timer.

    See Also:

        :meth:`wrf.omp_get_wtime`

    """
    return fomp_get_wtick()
