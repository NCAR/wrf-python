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
from .decorators import (handle_left_iter, uvmet_left_iter, 
                                handle_casting, handle_extract_transpose)

__all__ = ["FortranException", "computeslp", "computetk", "computetd", 
           "computerh", "computeavo", "computepvo", "computeeth", 
           "computeuvmet","computeomega", "computetv", "computesrh", 
           "computeuh", "computepw","computedbz","computecape", 
           "computeij", "computell", "computeeta", "computectt"]

class FortranException(Exception):
    def __call__(self, message):
        raise self.__class__(message)

@handle_left_iter(3,0, ignore_args=(2,3))
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def interpz3d(data3d, zdata, desiredloc, missingval):
    res = f_interpz3d(data3d, 
                      zdata, 
                      desiredloc, 
                      missingval)
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def interp2dxy(data3d,xy):
    res = f_interp2dxy(data3d,
                       xy)
    return res

@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def interp1d(v_in, z_in, z_out, missingval):
    res = f_interp1d(v_in,
                     z_in,
                     z_out,
                     missingval)
    
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2,3))
@handle_extract_transpose()
def computeslp(z, t, p, q):
    t_surf = np.zeros((z.shape[-2], z.shape[-1]), np.float64, order="F")
    t_sea_level = np.zeros((z.shape[-2], z.shape[-1]), np.float64, order="F")
    level = np.zeros((z.shape[-2], z.shape[-1]), np.int32, order="F")
    
    res = f_computeslp(z, 
                       t, 
                       p, 
                       q,
                       t_sea_level, # Should come in with fortran ordering
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
    res = np.reshape(res, shape)
    return res

# Note: No left iteration decorator needed with 1D arrays
@handle_casting(arg_idxs=(0,1))
@handle_extract_transpose()
def computetd(pressure, qv_in):
    shape = pressure.shape
    res = f_computetd(pressure.ravel(order="A"), 
                      qv_in.ravel(order="A"))
    res = np.reshape(res, shape)
    return res

# Note:  No decorator needed with 1D arrays
@handle_casting(arg_idxs=(0,1,2))
@handle_extract_transpose()
def computerh(qv, q, t):
    shape = qv.shape
    res = f_computerh(qv.ravel(order="A"),
                      q.ravel(order="A"),
                      t.ravel(order="A"))
    res = np.reshape(res, shape)
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
    longca = np.zeros((lat.shape[-2], lat.shape[-1]), np.float64, order="F")
    longcb = np.zeros((lon.shape[-2], lon.shape[-1]), np.float64, order="F")
    rpd = Constants.PI/180.
    
    
    # Make the 2D array a 3D array with 1 dimension
    if u.ndim < 3:
        u = u.reshape((u.shape[-2], u.shape[-1], 1), order="F")
        v = v.reshape((v.shape[-2], v.shape[-1], 1), order="F")

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
    
    #return res.T
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
    
    tem1 = np.zeros((u.shape[-3],u.shape[-2],u.shape[-1]), np.float64, 
                    order="F")
    tem2 = np.zeros((u.shape[-3],u.shape[-2],u.shape[-1]), np.float64, 
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
    # Note, dim -3 is height, we only want y and x
    zdiff = np.zeros((p.shape[-2], p.shape[-1]), np.float64, order="F")
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
    flip_cape = np.zeros((p_hpa.shape[-3],p_hpa.shape[-2],p_hpa.shape[-1]), 
                        np.float64, order="F")
    flip_cin = np.zeros((p_hpa.shape[-3],p_hpa.shape[-2],p_hpa.shape[-1]), 
                       np.float64, order="F")
    PSADITHTE, PSADIPRS, PSADITMK = get_lookup_tables()
    PSADITMK = PSADITMK.T
    
    # The fortran routine needs pressure to be ascending in z-direction, 
    # along with tk,qv,and ht.
    flip_p = p_hpa[::-1,:,:]
    flip_tk = tk[::-1,:,:]
    flip_qv = qv[::-1,:,:]
    flip_ht = ht[::-1,:,:]
    
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
    
    # Don't need to transpose since we only passed a view to fortran
    # Remember to flip cape and cin back to descending p coordinates
    cape = flip_cape[::-1,:,:]
    cin = flip_cin[::-1,:,:]
    
    return (cape, cin)

def computeij(map_proj, truelat1, truelat2, stdlon,
               lat1, lon1, pole_lat, pole_lon,
               knowni, knownj, dx, latinc, loninc, lat, lon):
    
    res = f_lltoij(map_proj,truelat1,truelat2,stdlon,
                   lat1,lon1,pole_lat,pole_lon,
                   knowni,knownj,dx,latinc,loninc,lat,lon,
                   FortranException())
    
    return res

def computell(map_proj, truelat1, truelat2, stdlon, lat1, lon1,
             pole_lat, pole_lon, knowni, knownj, dx, latinc,
             loninc, i, j):
    
    res = f_ijtoll(map_proj,truelat1,truelat2,stdlon,lat1,lon1,
                   pole_lat,pole_lon,knowni,knownj,dx,latinc,
                   loninc,i,j,FortranException())
    
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
    
    field_copy = field.copy()
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
    

    
    
    
    
    
    

