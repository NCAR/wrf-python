import numpy as n
import numpy.ma as ma

from wrf.var.constants import Constants
from wrf.var.psadlookup import get_lookup_tables
from wrf.var._wrfext import (f_interpz3d, f_interp2dxy,f_interp1d,
                     f_computeslp, f_computetk, f_computetd, f_computerh, 
                     f_computeabsvort,f_computepvo, f_computeeth, 
                     f_computeuvmet, 
                     f_computeomega, f_computetv, f_computewetbulb,
                     f_computesrh, f_computeuh, f_computepw, f_computedbz,
                     f_lltoij, f_ijtoll, f_converteta, f_computectt)
from wrf.var._wrfcape import f_computecape
from wrf.var.decorators import handle_left_iter

__all__ = ["FortranException", "computeslp", "computetk", "computetd", 
           "computerh", "computeavo", "computepvo", "computeeth", 
           "computeuvmet","computeomega", "computetv", "computesrh", 
           "computeuh", "computepw","computedbz","computecape", 
           "computeij", "computell", "computeeta", "computectt"]

class FortranException(Exception):
    def __call__(self, message):
        raise self.__class__(message)

def interpz3d(data3d,zdata,desiredloc,missingval):
    res = f_interpz3d(data3d.astype("float64").T, 
                      zdata.astype("float64").T, 
                      desiredloc, 
                      missingval)
    return res.astype("float32").T

def interpz2d(data3d,xy):
    res = f_interp2dxy(data3d.astype("float64").T,
                       xy.astype("float64").T)
    # Note: Fortran routine does not support missing values, so no masking
    return res.astype("float32").T

def interp1d(v_in,z_in,z_out,missingval):
    res = f_interp1d(v_in.astype("float64"),
                     z_in.astype("float64"),
                     z_out.astype("float64"),
                     missingval)
    
    return res.astype("float32")

@handle_left_iter(2,3,0)
def computeslp(z,t,p,q):
    t_surf = n.zeros((z.shape[-2], z.shape[-1]), "float64")
    t_sea_level = n.zeros((z.shape[-2], z.shape[-1]), "float64")
    level = n.zeros((z.shape[-2], z.shape[-1]), "int32")
    
    res = f_computeslp(z.astype("float64").T, 
                       t.astype("float64").T, 
                       p.astype("float64").T, 
                       q.astype("float64").T,
                       t_sea_level.T, 
                       t_surf.T, 
                       level.T,
                       FortranException())
    
    return res.astype("float32").T

@handle_left_iter(3,3,0)
def computetk(pres, theta):
    # No need to transpose here since operations on 1D array
    shape = pres.shape
    res = f_computetk(pres.astype("float64").flatten("A"), 
                      theta.astype("float64").flatten("A"))
    res = n.reshape(res, shape, "A")
    return res.astype("float32")

# Note: No decorator needed with 1D arrays
def computetd(pressure,qv_in):
    shape = pressure.shape
    res = f_computetd(pressure.astype("float64").flatten("A"), qv_in.astype("float64").flatten("A"))
    res = n.reshape(res, shape, "A")
    return res.astype("float32")

# Note:  No decorator needed with 1D arrays
def computerh(qv,q,t):
    shape = qv.shape
    res = f_computerh(qv.astype("float64").flatten("A"),
                      q.astype("float64").flatten("A"),
                      t.astype("float64").flatten("A"))
    res = n.reshape(res, shape, "A")
    return res.astype("float32")

@handle_left_iter(3,3,0, ignore_args=(6,7))
def computeavo(u,v,msfu,msfv,msfm,cor,dx,dy):
    res = f_computeabsvort(u.astype("float64").T,
                           v.astype("float64").T,
                           msfu.astype("float64").T,
                           msfv.astype("float64").T,
                           msfm.astype("float64").T,
                           cor.astype("float64").T,
                           dx,
                           dy)
    
    return res.astype("float32").T

@handle_left_iter(3,3,0, ignore_args=(8,9))
def computepvo(u,v,theta,prs,msfu,msfv,msfm,cor,dx,dy):
    
    res = f_computepvo(u.astype("float64").T,
                       v.astype("float64").T,
                       theta.astype("float64").T,
                       prs.astype("float64").T,
                       msfu.astype("float64").T,
                       msfv.astype("float64").T,
                       msfm.astype("float64").T,
                       cor.astype("float64").T,
                       dx,
                       dy)
    
    return res.astype("float32").T

@handle_left_iter(3,3,0)
def computeeth(qv, tk, p):
    
    res = f_computeeth(qv.astype("float64").T,
                       tk.astype("float64").T,
                       p.astype("float64").T)
    
    return res.astype("float32").T

@handle_left_iter(4,2,2, ignore_args=(4,5), alg_out_fixed_dims=(2,))
def computeuvmet(u,v,lat,lon,cen_long,cone):
    longca = n.zeros((lat.shape[-2], lat.shape[-1]), "float64")
    longcb = n.zeros((lon.shape[-2], lon.shape[-1]), "float64")
    rpd = Constants.PI/180.
    
    
    # Make the 2D array a 3D array with 1 dimension
    if u.ndim != 3:
        u = u.reshape((1,u.shape[-2], u.shape[-1]))
        v = v.reshape((1,v.shape[-2], v.shape[-1]))

    # istag will always be false since winds are destaggered already
    # Missing values don't appear to be used, so setting to 0
    res = f_computeuvmet(u.astype("float64").T,
               v.astype("float64").T,
               longca.T,
               longcb.T,
               lon.astype("float64").T,
               lat.astype("float64").T,
               cen_long,
               cone,
               rpd,
               0,
               False,
               0,
               0,
               0)

    
    return res.astype("float32").T

@handle_left_iter(3,3,0)
def computeomega(qv, tk, w, p):
    
    res = f_computeomega(qv.astype("float64").T,
                    tk.astype("float64").T,
                    w.astype("float64").T,
                    p.astype("float64").T)
    
    #return res.T
    return res.astype("float32").T

@handle_left_iter(3,3,0)
def computetv(tk,qv):
    res = f_computetv(tk.astype("float64").T,
                      qv.astype("float64").T)
    
    return res.astype("float32").T

@handle_left_iter(3,3,0)
def computewetbulb(p,tk,qv):
    PSADITHTE, PSADIPRS, PSADITMK = get_lookup_tables()
    
    res = f_computewetbulb(p.astype("float64").T,
                     tk.astype("float64").T,
                     qv.astype("float64").T,
                     PSADITHTE,
                     PSADIPRS,
                     PSADITMK.T,
                     FortranException())
    
    return res.astype("float32").T

@handle_left_iter(2,3,0, ignore_args=(4,))
def computesrh(u, v, z, ter, top):
    
    res = f_computesrh(u.astype("float64").T, 
                       v.astype("float64").T, 
                       z.astype("float64").T, 
                       ter.astype("float64").T, 
                       top)
    
    return res.astype("float32").T

@handle_left_iter(2,3,2, ignore_args=(5,6,7,8))
def computeuh(zp, mapfct, u, v, wstag, dx, dy, bottom, top):
    
    tem1 = n.zeros((u.shape[-3],u.shape[-2],u.shape[-1]), "float64")
    tem2 = n.zeros((u.shape[-3],u.shape[-2],u.shape[-1]), "float64")
    
    res = f_computeuh(zp.astype("float64").T,
                      mapfct.astype("float64").T,
                      dx,
                      dy,
                      bottom,
                      top,
                      u.astype("float64").T,
                      v.astype("float64").T,
                      wstag.astype("float64").T,
                      tem1.T,
                      tem2.T)
    
    return res.astype("float32").T

@handle_left_iter(2,3,0)
def computepw(p,tv,qv,ht):
    # Note, dim 0 is height, we only want y and x
    zdiff = n.zeros((p.shape[-2], p.shape[-1]), "float64")
    res = f_computepw(p.astype("float64").T,
                      tv.astype("float64").T,
                      qv.astype("float64").T,
                      ht.astype("float64").T,
                      zdiff.T)
    
    return res.astype("float32").T

@handle_left_iter(3,3,0, ignore_args=(6,7,8))
def computedbz(p,tk,qv,qr,qs,qg,sn0,ivarint,iliqskin):
    
    res = f_computedbz(p.astype("float64").T,
                       tk.astype("float64").T,
                       qv.astype("float64").T,
                       qr.astype("float64").T,
                       qs.astype("float64").T,
                       qg.astype("float64").T,
                       sn0,
                       ivarint,
                       iliqskin)
    
    return res.astype("float32").T

@handle_left_iter(3,3,0,ignore_args=(6,7,8))
def computecape(p_hpa,tk,qv,ht,ter,sfp,missing,i3dflag,ter_follow):
    flip_cape = n.zeros((p_hpa.shape[-3],p_hpa.shape[-2],p_hpa.shape[-1]), "float64")
    flip_cin = n.zeros((p_hpa.shape[-3],p_hpa.shape[-2],p_hpa.shape[-1]), "float64")
    PSADITHTE, PSADIPRS, PSADITMK = get_lookup_tables()
    
    # The fortran routine needs pressure to be ascending in z-direction, 
    # along with tk,qv,and ht.
    flip_p = p_hpa[::-1,:,:]
    flip_tk = tk[::-1,:,:]
    flip_qv = qv[::-1,:,:]
    flip_ht = ht[::-1,:,:]
    
    f_computecape(flip_p.astype("float64").T,
                  flip_tk.astype("float64").T,
                  flip_qv.astype("float64").T,
                  flip_ht.astype("float64").T,
                  ter.astype("float64").T,
                  sfp.astype("float64").T,
                  flip_cape.T,
                  flip_cin.T,
                  PSADITHTE,
                  PSADIPRS,
                  PSADITMK.T,
                  missing,
                  i3dflag,
                  ter_follow,
                  FortranException())
    
    # Don't need to transpose since we only passed a view to fortran
    cape = flip_cape.astype("float32")
    cin = flip_cin.astype("float32")
    # Remember to flip cape and cin back to descending p coordinates
    return (cape[::-1,:,:],cin[::-1,:,:])

# TODO:  This should handle lists of coords
def computeij(map_proj,truelat1,truelat2,stdlon,
               lat1,lon1,pole_lat,pole_lon,
               knowni,knownj,dx,latinc,loninc,lat,lon):
    
    res = f_lltoij(map_proj,truelat1,truelat2,stdlon,
                   lat1,lon1,pole_lat,pole_lon,
                   knowni,knownj,dx,latinc,loninc,lat,lon,
                   FortranException())
    
    return res[0],res[1]

# TODO:  This should handle lists of coords
def computell(map_proj,truelat1,truelat2,stdlon,lat1,lon1,
             pole_lat,pole_lon,knowni,knownj,dx,latinc,
             loninc,i,j):
    
    res = f_ijtoll(map_proj,truelat1,truelat2,stdlon,lat1,lon1,
                   pole_lat,pole_lon,knowni,knownj,dx,latinc,
                   loninc,i,j,FortranException())
    
    # Want lon,lat
    return res[1],res[0]

@handle_left_iter(3,3,0, ignore_args=(3,))
def computeeta(full_t, znu, psfc, ptop):
    pcalc = n.zeros(full_t.shape, "float64")
    mean_t = n.zeros(full_t.shape, "float64")
    temp_t = n.zeros(full_t.shape, "float64")
    
    res = f_converteta(full_t.astype("float64").T, 
                       znu.astype("float64"), 
                       psfc.astype("float64").T, 
                       ptop, 
                       pcalc.T, 
                       mean_t.T, 
                       temp_t.T)
    
    return res.astype("float32").T

@handle_left_iter(2,3,0,ignore_args=(7,))
def computectt(p_hpa,tk,qice,qcld,qv,ght,ter,haveqci):
    res = f_computectt(p_hpa.astype("float64").T,
                    tk.astype("float64").T,
                    qice.astype("float64").T,
                    qcld.astype("float64").T,
                    qv.astype("float64").T,
                    ght.astype("float64").T,
                    ter.astype("float64").T,
                    haveqci)
    
    return res.astype("float32").T

    
    
    
    
    
    

