import numpy as n

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
from wrf.var.decorators import (handle_left_iter, uvmet_left_iter, 
                                handle_casting)

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
def interpz3d(data3d,zdata,desiredloc,missingval):
    res = f_interpz3d(data3d.T, 
                      zdata.T, 
                      desiredloc, 
                      missingval)
    return res.T

@handle_left_iter(2,0)
@handle_casting(arg_idxs=(0,1))
def interp2dxy(data3d,xy):
    res = f_interp2dxy(data3d.T,
                       xy.T)
    return res.T

@handle_casting(arg_idxs=(0,1,2))
def interp1d(v_in,z_in,z_out,missingval):
    res = f_interp1d(v_in,
                     z_in,
                     z_out,
                     missingval)
    
    return res

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2,3))
def computeslp(z,t,p,q):
    t_surf = n.zeros((z.shape[-2], z.shape[-1]), "float64")
    t_sea_level = n.zeros((z.shape[-2], z.shape[-1]), "float64")
    level = n.zeros((z.shape[-2], z.shape[-1]), "int32")
    
    res = f_computeslp(z.T, 
                       t.T, 
                       p.T, 
                       q.T,
                       t_sea_level.T, 
                       t_surf.T, 
                       level.T,
                       FortranException())
    
    return res.T

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1))
def computetk(pres, theta):
    # No need to transpose here since operations on 1D array
    shape = pres.shape
    res = f_computetk(pres.flatten("A"), 
                      theta.flatten("A"))
    res = n.reshape(res, shape, "A")
    return res

# Note: No left iteration decorator needed with 1D arrays
@handle_casting(arg_idxs=(0,1))
def computetd(pressure,qv_in):
    shape = pressure.shape
    res = f_computetd(pressure.flatten("A"), 
                      qv_in.flatten("A"))
    res = n.reshape(res, shape, "A")
    return res

# Note:  No decorator needed with 1D arrays
@handle_casting(arg_idxs=(0,1,2))
def computerh(qv,q,t):
    shape = qv.shape
    res = f_computerh(qv.flatten("A"),
                      q.flatten("A"),
                      t.flatten("A"))
    res = n.reshape(res, shape, "A")
    return res

@handle_left_iter(3,0, ignore_args=(6,7))
@handle_casting(arg_idxs=(0,1,2,3,4,5))
def computeavo(u,v,msfu,msfv,msfm,cor,dx,dy):
    res = f_computeabsvort(u.T,
                           v.T,
                           msfu.T,
                           msfv.T,
                           msfm.T,
                           cor.T,
                           dx,
                           dy)
    
    return res.T

@handle_left_iter(3,2, ignore_args=(8,9))
@handle_casting(arg_idxs=(0,1,2,3,4,5,6,7))
def computepvo(u,v,theta,prs,msfu,msfv,msfm,cor,dx,dy):
    
    res = f_computepvo(u.T,
                       v.T,
                       theta.T,
                       prs.T,
                       msfu.T,
                       msfv.T,
                       msfm.T,
                       cor.T,
                       dx,
                       dy)
    
    return res.T

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2))
def computeeth(qv, tk, p):
    
    res = f_computeeth(qv.T,
                       tk.T,
                       p.T)
    
    return res.T

@uvmet_left_iter()
@handle_casting(arg_idxs=(0,1,2,3))
def computeuvmet(u,v,lat,lon,cen_long,cone):
    longca = n.zeros((lat.shape[-2], lat.shape[-1]), "float64")
    longcb = n.zeros((lon.shape[-2], lon.shape[-1]), "float64")
    rpd = Constants.PI/180.
    
    
    # Make the 2D array a 3D array with 1 dimension
    if u.ndim < 3:
        u = u.reshape((1,u.shape[-2], u.shape[-1]))
        v = v.reshape((1,v.shape[-2], v.shape[-1]))

    # istag will always be false since winds are destaggered already
    # Missing values don't appear to be used, so setting to 0
    res = f_computeuvmet(u.T,
               v.T,
               longca.T,
               longcb.T,
               lon.T,
               lat.T,
               cen_long,
               cone,
               rpd,
               0,
               False,
               0,
               0,
               0)

    
    return n.squeeze(res.T)

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2,3))
def computeomega(qv, tk, w, p):
    
    res = f_computeomega(qv.T,
                    tk.T,
                    w.T,
                    p.T)
    
    #return res.T
    return res.T

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1))
def computetv(tk,qv):
    res = f_computetv(tk.T,
                      qv.T)
    
    return res.T

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2))
def computewetbulb(p,tk,qv):
    PSADITHTE, PSADIPRS, PSADITMK = get_lookup_tables()
    
    res = f_computewetbulb(p.T,
                     tk.T,
                     qv.T,
                     PSADITHTE,
                     PSADIPRS,
                     PSADITMK.T,
                     FortranException())
    
    return res.T

@handle_left_iter(3,0, ignore_args=(4,))
@handle_casting(arg_idxs=(0,1,2,3))
def computesrh(u, v, z, ter, top):

    res = f_computesrh(u.T, 
                       v.T, 
                       z.T, 
                       ter.T, 
                       top)
    
    return res.T

@handle_left_iter(3,2, ignore_args=(5,6,7,8))
@handle_casting(arg_idxs=(0,1,2,3,4))
def computeuh(zp, mapfct, u, v, wstag, dx, dy, bottom, top):
    
    tem1 = n.zeros((u.shape[-3],u.shape[-2],u.shape[-1]), "float64")
    tem2 = n.zeros((u.shape[-3],u.shape[-2],u.shape[-1]), "float64")
    
    res = f_computeuh(zp.T,
                      mapfct.T,
                      dx,
                      dy,
                      bottom,
                      top,
                      u.T,
                      v.T,
                      wstag.T,
                      tem1.T,
                      tem2.T)
    
    return res.T

@handle_left_iter(3,0)
@handle_casting(arg_idxs=(0,1,2,3))
def computepw(p,tv,qv,ht):
    # Note, dim -3 is height, we only want y and x
    zdiff = n.zeros((p.shape[-2], p.shape[-1]), "float64")
    res = f_computepw(p.T,
                      tv.T,
                      qv.T,
                      ht.T,
                      zdiff.T)
    
    return res.T

@handle_left_iter(3,0, ignore_args=(6,7,8))
@handle_casting(arg_idxs=(0,1,2,3,4,5))
def computedbz(p,tk,qv,qr,qs,qg,sn0,ivarint,iliqskin):
    
    res = f_computedbz(p.T,
                       tk.T,
                       qv.T,
                       qr.T,
                       qs.T,
                       qg.T,
                       sn0,
                       ivarint,
                       iliqskin)
    
    return res.T

@handle_left_iter(3,0,ignore_args=(6,7,8))
@handle_casting(arg_idxs=(0,1,2,3,4,5))
def computecape(p_hpa,tk,qv,ht,ter,sfp,missing,i3dflag,ter_follow):
    flip_cape = n.zeros((p_hpa.shape[-3],p_hpa.shape[-2],p_hpa.shape[-1]), 
                        "float64")
    flip_cin = n.zeros((p_hpa.shape[-3],p_hpa.shape[-2],p_hpa.shape[-1]), 
                       "float64")
    PSADITHTE, PSADIPRS, PSADITMK = get_lookup_tables()
    
    # The fortran routine needs pressure to be ascending in z-direction, 
    # along with tk,qv,and ht.
    flip_p = p_hpa[::-1,:,:]
    flip_tk = tk[::-1,:,:]
    flip_qv = qv[::-1,:,:]
    flip_ht = ht[::-1,:,:]
    
    f_computecape(flip_p.T,
                  flip_tk.T,
                  flip_qv.T,
                  flip_ht.T,
                  ter.T,
                  sfp.T,
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
    # Remember to flip cape and cin back to descending p coordinates
    cape = flip_cape[::-1,:,:]
    cin = flip_cin[::-1,:,:]
    
    
    return (cape, cin)

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

@handle_left_iter(3,0, ignore_args=(3,))
@handle_casting(arg_idxs=(0,1,2))
def computeeta(full_t, znu, psfc, ptop):
    pcalc = n.zeros(full_t.shape, "float64")
    mean_t = n.zeros(full_t.shape, "float64")
    temp_t = n.zeros(full_t.shape, "float64")
    
    res = f_converteta(full_t.T, 
                       znu.T, 
                       psfc.T, 
                       ptop, 
                       pcalc.T, 
                       mean_t.T, 
                       temp_t.T)
    
    return res.T

@handle_left_iter(3,0,ignore_args=(7,))
@handle_casting(arg_idxs=(0,1,2,3,4,5,6))
def computectt(p_hpa,tk,qice,qcld,qv,ght,ter,haveqci):
    res = f_computectt(p_hpa.T,
                    tk.T,
                    qice.T,
                    qcld.T,
                    qv.T,
                    ght.T,
                    ter.T,
                    haveqci)
    
    return res.T

    
    
    
    
    
    

