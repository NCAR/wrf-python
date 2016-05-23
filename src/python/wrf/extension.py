from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

from .constants import Constants
from .psadlookup import get_lookup_tables
from ._wrfext import (f_interpz3d, f_interp2dxy,f_interp1d,
                     f_computeslp, f_computetk, f_computetd, f_computerh, 
                     f_computeabsvort,f_computepvo, f_computeeth, 
                     f_computeuvmet, 
                     f_computeomega, f_computetv, f_computewetbulb,
                     f_computesrh, f_computeuh, f_computepw, f_computedbz,
                     f_lltoij, f_ijtoll, f_converteta, f_computectt,
                     f_monotonic, f_filter2d, f_vintrp)
from ._wrfcape import f_computecape
from .decorators import (handle_left_iter, handle_casting, 
                         handle_extract_transpose)
from .util import range2
from .uvdecorator import uvmet_left_iter

__all__ = ["FortranException", "computeslp", "computetk", "computetd", 
           "computerh", "computeavo", "computepvo", "computeeth", 
           "computeuvmet","computeomega", "computetv", "computesrh", 
           "computeuh", "computepw","computedbz","computecape", 
           "computeij", "computell", "computeeta", "computectt",
           "interp2dxy", "interpz3d", "interp1d", "computeinterpline",
           "computevertcross"]

class FortranException(Exception):
    def __call__(self, message):
        raise self.__class__(message)


@handle_left_iter(3,0, ignore_args=(2,3))
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def interpz3d(field3d, z, desiredloc, missingval):
    res = f_interpz3d(field3d, 
                      z, 
                      desiredloc, 
                      missingval)
    return res

@handle_left_iter(3,0, ignore_args=(1,))
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def interp2dxy(field3d, xy):
    res = f_interp2dxy(field3d, 
                       xy)
    return res

@handle_left_iter(1, 0, ignore_args=(2,3))
@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def interp1d(v_in, z_in, z_out, missingval):
    res = f_interp1d(v_in,
                     z_in,
                     z_out,
                     missingval)
    
    return res

@handle_left_iter(3, 0, ignore_args=(1,4,3))
@handle_casting(arg_idxs=(0,))
@handle_extract_transpose(do_transpose=False)
def computevertcross(field3d, xy, var2dz, z_var2d, missingval):
    var2d = np.empty((z_var2d.size, xy.shape[0]), dtype=var2dz.dtype)
    var2dtmp = interp2dxy(field3d, xy)
    
    for i in range2(xy.shape[0]):
        var2d[:,i] = interp1d(var2dtmp[:,i], var2dz[:,i], z_var2d, missingval)
    
    return var2d

@handle_left_iter(2, 0, ignore_args=(1,))
@handle_casting(arg_idxs=(0,))
@handle_extract_transpose(do_transpose=False)
def computeinterpline(field2d, xy):

    tmp_shape = [1] + [x for x in field2d.shape]
    var2dtmp = np.empty(tmp_shape, field2d.dtype)
    var2dtmp[0,:,:] = field2d[:,:]
    
    var1dtmp = interp2dxy(var2dtmp, xy)
    
    return var1dtmp[0,:]

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2,3))
@handle_extract_transpose()
def computeslp(z, t, p, q):
    
    t_surf = np.zeros(z.shape[0:2], np.float64, order="F")
    t_sea_level = np.zeros(z.shape[0:2], np.float64, order="F")
    level = np.zeros(z.shape[0:2], np.int32, order="F")
    
    res = f_computeslp(z, 
                       t, 
                       p, 
                       q,
                       t_sea_level,
                       t_surf, 
                       level,
                       FortranException())
    
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def computetk(pres, theta):
    # No need to transpose here since operations on 1D array
    shape = pres.shape
    res = f_computetk(pres.ravel(order="A"), 
                      theta.ravel(order="A"))
    res = np.reshape(res, shape, order="F")
    
    return res

# Note: No left iteration decorator needed with 1D arrays
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def computetd(pressure, qv_in):
    shape = pressure.shape
    res = f_computetd(pressure.ravel(order="A"), 
                      qv_in.ravel(order="A"))
    res = np.reshape(res, shape, order="F")
    return res

# Note:  No decorator needed with 1D arrays
@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def computerh(qv, q, t):
    shape = qv.shape
    res = f_computerh(qv.ravel(order="A"),
                      q.ravel(order="A"),
                      t.ravel(order="A"))
    res = np.reshape(res, shape, order="F")
    return res

@handle_left_iter(3,0, ignore_args=(6,7))
@handle_casting(arg_idxs=(0,1,2,3,4,5))
@handle_extract_transpose()
def computeavo(u, v, msfu, msfv, msfm, cor, dx, dy):
    res = f_computeabsvort(u,
                           v,
                           msfu,
                           msfv,
                           msfm,
                           cor,
                           dx,
                           dy)
    
    return res

@handle_left_iter(3,2, ignore_args=(8,9))
@handle_casting(arg_idxs=(0,1,2,3,4,5,6,7))
@handle_extract_transpose()
def computepvo(u, v, theta, prs, msfu, msfv, msfm, cor, dx, dy):
    
    res = f_computepvo(u,
                       v,
                       theta,
                       prs,
                       msfu,
                       msfv,
                       msfm,
                       cor,
                       dx,
                       dy)
    
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def computeeth(qv, tk, p):
    
    res = f_computeeth(qv,
                       tk,
                       p)
    
    return res

@uvmet_left_iter()
@handle_casting(arg_idxs=(0,1,2,3))
@handle_extract_transpose()
def computeuvmet(u, v, lat, lon, cen_long, cone):
    longca = np.zeros((lat.shape[0], lat.shape[1]), np.float64, order="F")
    longcb = np.zeros((lon.shape[0], lon.shape[1]), np.float64, order="F")
    rpd = Constants.PI/180.
    
    
    # Make the 2D array a 3D array with 1 dimension
    if u.ndim < 3:
        u = u.reshape((u.shape[0], u.shape[1], 1), order="F")
        v = v.reshape((v.shape[0], v.shape[1], 1), order="F")

    # istag will always be false since winds are destaggered already
    # Missing values don't appear to be used, so setting to 0
    res = f_computeuvmet(u,
               v,
               longca,
               longcb,
               lon,
               lat,
               cen_long,
               cone,
               rpd,
               0,
               False,
               0,
               0,
               0)

    
    return np.squeeze(res)

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2,3))
@handle_extract_transpose()
def computeomega(qv, tk, w, p):
    
    res = f_computeomega(qv,
                    tk,
                    w,
                    p)
    
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def computetv(tk, qv):
    res = f_computetv(tk,
                      qv)
    
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def computewetbulb(p, tk, qv):
    PSADITHTE, PSADIPRS, PSADITMK = get_lookup_tables()
    PSADITMK = PSADITMK.T
    
    res = f_computewetbulb(p,
                     tk,
                     qv,
                     PSADITHTE,
                     PSADIPRS,
                     PSADITMK,
                     FortranException())
    
    return res

@handle_left_iter(3,0, ignore_args=(4,))
@handle_casting(arg_idxs=(0,1,2,3))
@handle_extract_transpose()
def computesrh(u, v, z, ter, top):

    res = f_computesrh(u, 
                       v, 
                       z, 
                       ter, 
                       top)
    
    return res

@handle_left_iter(3,2, ignore_args=(5,6,7,8))
@handle_casting(arg_idxs=(0,1,2,3,4))
@handle_extract_transpose()
def computeuh(zp, mapfct, u, v, wstag, dx, dy, bottom, top):
    
    tem1 = np.zeros((u.shape[0], u.shape[1], u.shape[2]), np.float64, 
                    order="F")
    tem2 = np.zeros((u.shape[0], u.shape[1], u.shape[2]), np.float64, 
                    order="F")
    
    res = f_computeuh(zp,
                      mapfct,
                      dx,
                      dy,
                      bottom,
                      top,
                      u,
                      v,
                      wstag,
                      tem1,
                      tem2)
    
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2,3))
@handle_extract_transpose()
def computepw(p, tv, qv, ht):

    zdiff = np.zeros((p.shape[0], p.shape[1]), np.float64, order="F")
    res = f_computepw(p,
                      tv,
                      qv,
                      ht,
                      zdiff)
    
    return res

@handle_left_iter(3,0, ignore_args=(6,7,8))
@handle_casting(arg_idxs=(0,1,2,3,4,5))
@handle_extract_transpose()
def computedbz(p, tk, qv, qr, qs, qg, sn0, ivarint, iliqskin):
    
    res = f_computedbz(p,
                       tk,
                       qv,
                       qr,
                       qs,
                       qg,
                       sn0,
                       ivarint,
                       iliqskin)
    
    return res

@handle_left_iter(3,0,ignore_args=(6,7,8))
@handle_casting(arg_idxs=(0,1,2,3,4,5))
@handle_extract_transpose()
def computecape(p_hpa, tk, qv, ht, ter, sfp, missing, i3dflag, ter_follow):
    flip_cape = np.zeros((p_hpa.shape[0], p_hpa.shape[1], p_hpa.shape[2]), 
                        np.float64, order="F")
    flip_cin = np.zeros((p_hpa.shape[0], p_hpa.shape[1], p_hpa.shape[2]), 
                       np.float64, order="F")
    PSADITHTE, PSADIPRS, PSADITMK = get_lookup_tables()
    PSADITMK = PSADITMK.T
    
    # The fortran routine needs pressure to be ascending in z-direction, 
    # along with tk,qv,and ht.
    # The extra mumbo-jumbo is so that the view created by numpy is fortran
    # contiguous.  'ascontiguousarray' only works in C ordering, hence the 
    # extra transposes.
    flip_p = np.ascontiguousarray(p_hpa[:,:,::-1].T).T
    flip_tk = np.ascontiguousarray(tk[:,:,::-1].T).T
    flip_qv = np.ascontiguousarray(qv[:,:,::-1].T).T
    flip_ht = np.ascontiguousarray(ht[:,:,::-1].T).T
    
    f_computecape(flip_p,
                  flip_tk,
                  flip_qv,
                  flip_ht,
                  ter,
                  sfp,
                  flip_cape,
                  flip_cin,
                  PSADITHTE,
                  PSADIPRS,
                  PSADITMK,
                  missing,
                  i3dflag,
                  ter_follow,
                  FortranException())
    
    # Need to flip the vertical back to decending pressure with height.
    cape = np.ascontiguousarray(flip_cape[:,:,::-1].T).T
    cin = np.ascontiguousarray(flip_cin[:,:,::-1].T).T

    return (cape, cin)

def computeij(map_proj, truelat1, truelat2, stdlon,
               lat1, lon1, pole_lat, pole_lon,
               knowni, knownj, dx, latinc, loninc, lat, lon):
    
    res = f_lltoij(map_proj,
                   truelat1,
                   truelat2,
                   stdlon,
                   lat1,
                   lon1,
                   pole_lat,
                   pole_lon,
                   knowni,
                   knownj,
                   dx,
                   latinc,
                   loninc,
                   lat,
                   lon,
                   FortranException())
    
    return res

def computell(map_proj, truelat1, truelat2, stdlon, lat1, lon1,
             pole_lat, pole_lon, knowni, knownj, dx, latinc,
             loninc, i, j):
    
    res = f_ijtoll(map_proj,
                   truelat1,
                   truelat2,
                   stdlon,
                   lat1,
                   lon1,
                   pole_lat,
                   pole_lon,
                   knowni,
                   knownj,
                   dx,
                   latinc,
                   loninc,
                   i,
                   j,
                   FortranException())
    
    return res

@handle_left_iter(3,0, ignore_args=(3,))
@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def computeeta(full_t, znu, psfc, ptop):
    pcalc = np.zeros(full_t.shape, np.float64, order="F")
    mean_t = np.zeros(full_t.shape, np.float64, order="F")
    temp_t = np.zeros(full_t.shape, np.float64, order="F")
    
    res = f_converteta(full_t, 
                       znu, 
                       psfc, 
                       ptop, 
                       pcalc, 
                       mean_t, 
                       temp_t)
    
    return res

@handle_left_iter(3,0,ignore_args=(7,))
@handle_casting(arg_idxs=(0,1,2,3,4,5,6))
@handle_extract_transpose()
def computectt(p_hpa, tk, qice, qcld, qv, ght, ter, haveqci):
    res = f_computectt(p_hpa,
                    tk,
                    qice,
                    qcld,
                    qv,
                    ght,
                    ter,
                    haveqci)
    
    return res

@handle_left_iter(2,0,ignore_args=(1,))
@handle_casting(arg_idxs=(0,))
@handle_extract_transpose()
def smooth2d(field, passes):
    # Unlike NCL, this routine will not modify the values in place, but 
    # copies the original data before modifying it.  This allows the decorator
    # to work properly and also behaves like the other methods.
    
    if isinstance(field, np.ma.MaskedArray):
        missing = field.fill_value
    else:
        missing = Constants.DEFAULT_FILL
    
    field_copy = field.copy(order="A")
    field_tmp = np.zeros(field_copy.shape, field_copy.dtype, order="F")  
    
    f_filter2d(field_copy, 
               field_tmp, 
               missing,
               passes)
    
    # Don't transpose here since the fortran routine is not returning an
    # array.  It's only modifying the existing array.
    return field_copy
    
@handle_left_iter(3,0,ignore_args=(3,4,5))
@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def monotonic(var, lvprs, coriolis, idir, delta, icorsw):
    res = f_monotonic(var,
                      lvprs,
                      coriolis,
                      idir,
                      delta,
                      icorsw)
    
    return res

@handle_left_iter(3,0,ignore_args=(9,10,11,12,13,14))
@handle_casting(arg_idxs=(0,1,2,3,4,5,6,7,8,9))
@handle_extract_transpose()
def vintrp(field, pres, tk, qvp, ght, terrain, sfp, smsfp,
           vcarray, interp_levels, icase, extrap, vcor, logp,
           missing):
    
    res = f_vintrp(field,
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
             missing)
    
    return res
    

    
    
    
    
    
    

