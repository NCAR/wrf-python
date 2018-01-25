from math import fabs, log, tan, sin, cos

import unittest as ut
import numpy.testing as nt 
import numpy as np
import numpy.ma as ma
import os, sys
import subprocess

from netCDF4 import Dataset as nc

from wrf import *

NCL_EXE = "/Users/ladwig/nclbuild/6.3.0/bin/ncl"
TEST_FILE = "/Users/ladwig/Documents/wrf_files/wrfout_d01_2010-06-13_21:00:00"
OUT_NC_FILE = "/tmp/wrftest.nc"
NCFILE = nc(TEST_FILE)
NCGROUP = [NCFILE, NCFILE, NCFILE]

# Python 3
if sys.version_info > (3,):
    xrange = range


ROUTINE_MAP = {"avo" : avo, 
          "eth" : eth, 
          "cape_2d" : cape_2d, 
          "cape_3d" : cape_3d, 
          "ctt" : ctt, 
          "dbz" : dbz, 
          "helicity" : srhel, 
          "omg" : omega, 
          "pvo" : pvo, 
          "pw" : pw,  
          "rh" : rh, 
          "slp" : slp,  
          "td" : td, 
          "tk" : tk, 
          "tv" : tvirtual, 
          "twb" : wetbulb, 
          "updraft_helicity" : udhel, 
          "uvmet" : uvmet,  
          "cloudfrac" : cloudfrac}

class ProjectionError(RuntimeError):
    pass

def get_args(varname, wrfnc, timeidx, method, squeeze):
    if varname == "avo":
        ncvars = extract_vars(wrfnc, timeidx, ("U", "V", "MAPFAC_U",
                                           "MAPFAC_V", "MAPFAC_M",
                                           "F"),
                          method, squeeze, cache=None, meta=True)
    
        attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
        u = ncvars["U"]
        v = ncvars["V"]
        msfu = ncvars["MAPFAC_U"]
        msfv = ncvars["MAPFAC_V"]
        msfm = ncvars["MAPFAC_M"]
        cor = ncvars["F"]
        
        dx = attrs["DX"]
        dy = attrs["DY"]
        
        return (u, v, msfu, msfv, msfm, cor, dx, dy)
    
    if varname == "pvo":
        ncvars = extract_vars(wrfnc, timeidx, ("U", "V", "T", "P",
                                           "PB", "MAPFAC_U",
                                           "MAPFAC_V", "MAPFAC_M",
                                           "F"),
                          method, squeeze, cache=None, meta=True)
        attrs = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
        
        u = ncvars["U"]
        v = ncvars["V"]
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        msfu = ncvars["MAPFAC_U"]
        msfv = ncvars["MAPFAC_V"]
        msfm = ncvars["MAPFAC_M"]
        cor = ncvars["F"]
        
        dx = attrs["DX"]
        dy = attrs["DY"]
        
        full_t = t + 300
        full_p = p + pb
        
        return (u, v, full_t, full_p, msfu, msfv, msfm, cor, dx, dy)
    
    if varname == "eth":
        varnames=("T", "P", "PB", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
    
        full_t = t + Constants.T_BASE
        full_p = p + pb
        tkel = tk(full_p, full_t, meta=False)
    
        return (qv, tkel, full_p)
    
    if varname == "cape_2d":
        varnames = ("T", "P", "PB", "QVAPOR", "PH","PHB", "HGT", "PSFC")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        ter = ncvars["HGT"]
        psfc = ncvars["PSFC"]
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        tkel = tk(full_p, full_t, meta=False)
        
        geopt = ph + phb
        geopt_unstag = destagger(geopt, -3)
        z = geopt_unstag/Constants.G
        
        # Convert pressure to hPa
        p_hpa = ConversionFactors.PA_TO_HPA * full_p
        psfc_hpa = ConversionFactors.PA_TO_HPA * psfc 
        
        i3dflag = 0
        ter_follow = 1
        
        return (p_hpa, tkel, qv, z, ter, psfc_hpa, ter_follow)
    
    if varname == "cape_3d":
        varnames = ("T", "P", "PB", "QVAPOR", "PH", "PHB", "HGT", "PSFC")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        ter = ncvars["HGT"]
        psfc = ncvars["PSFC"]
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        tkel = tk(full_p, full_t, meta=False)
        
        geopt = ph + phb
        geopt_unstag = destagger(geopt, -3)
        z = geopt_unstag/Constants.G
        
        # Convert pressure to hPa
        p_hpa = ConversionFactors.PA_TO_HPA * full_p
        psfc_hpa = ConversionFactors.PA_TO_HPA * psfc 
        
        i3dflag = 1
        ter_follow = 1
        
        return (p_hpa, tkel, qv, z, ter, psfc_hpa, ter_follow)
    
    if varname == "ctt":
        varnames = ("T", "P", "PB", "PH", "PHB", "HGT", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        ter = ncvars["HGT"]
        qv = ncvars["QVAPOR"] * 1000.0 # g/kg
        
        haveqci = 1
        try:
            icevars = extract_vars(wrfnc, timeidx, "QICE", 
                                   method, squeeze, cache=None, meta=False)
        except KeyError:
            qice = np.zeros(qv.shape, qv.dtype)
            haveqci = 0
        else:
            qice = icevars["QICE"] * 1000.0 #g/kg
        
        try:
            cldvars = extract_vars(wrfnc, timeidx, "QCLOUD", 
                                   method, squeeze, cache=None, meta=False)
        except KeyError:
            raise RuntimeError("'QCLOUD' not found in NetCDF file")
        else:
            qcld = cldvars["QCLOUD"] * 1000.0 #g/kg
        
        full_p = p + pb
        p_hpa = full_p * ConversionFactors.PA_TO_HPA
        full_t = t + Constants.T_BASE
        tkel = tk(full_p, full_t, meta=False)
        
        geopt = ph + phb
        geopt_unstag = destagger(geopt, -3)
        ght = geopt_unstag / Constants.G
        
        return (p_hpa, tkel, qv, qcld, ght, ter, qice)
    
    if varname == "dbz":
        varnames = ("T", "P", "PB", "QVAPOR", "QRAIN")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
        qr = ncvars["QRAIN"]
        
        try:
            snowvars = extract_vars(wrfnc, timeidx, "QSNOW", 
                                    method, squeeze, cache=None, meta=False)
        except KeyError:
            qs = np.zeros(qv.shape, qv.dtype)
        else:
            qs = snowvars["QSNOW"]
        
        try:
            graupvars = extract_vars(wrfnc, timeidx, "QGRAUP", 
                                     method, squeeze, cache=None, meta=False)
        except KeyError:
            qg = np.zeros(qv.shape, qv.dtype)
        else:
            qg = graupvars["QGRAUP"]
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        tkel = tk(full_p, full_t, meta=False)
        
        return (full_p, tkel, qv, qr, qs, qg)
    
    if varname == "helicity":
        # Top can either be 3000 or 1000 (for 0-1 srh or 0-3 srh)
        
        ncvars = extract_vars(wrfnc, timeidx, ("HGT", "PH", "PHB"),
                              method, squeeze, cache=None, meta=True)
        
        ter = ncvars["HGT"]
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        
        # As coded in NCL, but not sure this is possible
        varname = "U"
        u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, 
                              cache=None, meta=False)
        u = destagger(u_vars[varname], -1) 
        
        varname = "V"
        v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, 
                              cache=None, meta=False)
        v = destagger(v_vars[varname], -2)
    
        geopt = ph + phb
        geopt_unstag = destagger(geopt, -3)
        
        z = geopt_unstag / Constants.G
        
        return (u, v, z, ter)
    
    if varname == "updraft_helicity":
        ncvars = extract_vars(wrfnc, timeidx, ("W", "PH", "PHB", "MAPFAC_M"),
                          method, squeeze, cache=None, meta=True)
    
        wstag = ncvars["W"]
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        mapfct = ncvars["MAPFAC_M"]
        
        attrs  = extract_global_attrs(wrfnc, attrs=("DX", "DY"))
        dx = attrs["DX"]
        dy = attrs["DY"]
        
        # As coded in NCL, but not sure this is possible
        varname = "U"
        u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, 
                              cache=None, meta=True)
        u = destagger(u_vars[varname], -1, meta=True) 
        
        varname = "V"
        v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, 
                              cache=None, meta=True)
        v = destagger(v_vars[varname], -2, meta=True) 
        
        zstag = ph + phb
        
        return (zstag, mapfct, u, v, wstag, dx, dy)
    
    if varname == "omg":
        varnames=("T", "P", "W", "PB", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        w = ncvars["W"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
        
        wa = destagger(w, -3)
        full_t = t + Constants.T_BASE
        full_p = p + pb
        tkel = tk(full_p, full_t, meta=False)
        
        return (qv, tkel, wa, full_p)
    
    if varname == "pw":
        varnames=("T", "P", "PB", "PH", "PHB", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
    
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        qv = ncvars["QVAPOR"]
    
        # Change this to use real virtual temperature!
        full_p = p + pb
        ht = (ph + phb)/Constants.G
        full_t =  t + Constants.T_BASE
    
        tkel = tk(full_p, full_t, meta=False)
        
        return (full_p, tkel, qv, ht)
    
    if varname == "rh":
        varnames=("T", "P", "PB", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qvapor = to_np(ncvars["QVAPOR"])
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        qvapor[qvapor < 0] = 0
        tkel = tk(full_p, full_t, meta=False)
        
        return (qvapor, full_p, tkel)
    
    if varname == "slp":
        varnames=("T", "P", "PB", "QVAPOR", "PH", "PHB")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
    
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qvapor = to_np(ncvars["QVAPOR"])
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        qvapor[qvapor < 0] = 0.
        
        full_ph = (ph + phb) / Constants.G
        
        destag_ph = destagger(full_ph, -3)
        
        tkel = tk(full_p, full_t, meta=False)
        
        return (destag_ph, tkel, full_p, qvapor)
    
    if varname == "td":
        varnames=("P", "PB", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        
        p = ncvars["P"]
        pb = ncvars["PB"]
        qvapor = to_np(ncvars["QVAPOR"])
        
        # Algorithm requires hPa
        full_p = .01*(p + pb)
        qvapor[qvapor < 0] = 0
        
        return (full_p, qvapor)
    
    if varname == "tk":
        varnames=("T", "P", "PB")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        
        return (full_p, full_t)
    
    if varname == "tv":
        varnames=("T", "P", "PB", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        tkel = tk(full_p, full_t)
        
        return (tkel, qv)
    
    if varname == "twb":
        varnames=("T", "P", "PB", "QVAPOR")
        ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, 
                              cache=None, meta=True)
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        
        tkel = tk(full_p, full_t)
        
        return (full_p, tkel, qv)
    
    if varname == "uvmet":
        varname = "U"
        u_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, 
                              cache=None, meta=True)
        
        u = destagger(u_vars[varname], -1, meta=True)
        
        varname = "V"
        v_vars = extract_vars(wrfnc, timeidx, varname, method, squeeze, 
                              cache=None, meta=True)
        v = destagger(v_vars[varname], -2, meta=True)
        
        map_proj_attrs = extract_global_attrs(wrfnc, attrs="MAP_PROJ")
        map_proj = map_proj_attrs["MAP_PROJ"]
        
        if map_proj in (0,3,6):
            raise ProjectionError("Map projection does not need rotation")
        elif map_proj in (1,2):
            lat_attrs = extract_global_attrs(wrfnc, attrs=("TRUELAT1",
                                                           "TRUELAT2"))
            radians_per_degree = Constants.PI/180.0
            # Rotation needed for Lambert and Polar Stereographic
            true_lat1 = lat_attrs["TRUELAT1"]
            true_lat2 = lat_attrs["TRUELAT2"]
            
            try:
                lon_attrs = extract_global_attrs(wrfnc, attrs="STAND_LON")
            except AttributeError:
                try:
                    cen_lon_attrs = extract_global_attrs(wrfnc, attrs="CEN_LON")
                except AttributeError:
                    raise RuntimeError("longitude attributes not found in NetCDF")
                else:
                    cen_lon = cen_lon_attrs["CEN_LON"]
            else:
                cen_lon = lon_attrs["STAND_LON"]
            
            
            varname = "XLAT"
            xlat_var = extract_vars(wrfnc, timeidx, varname, 
                                    method, squeeze, cache=None, meta=True)
            lat = xlat_var[varname]
            
            varname = "XLONG"
            xlon_var = extract_vars(wrfnc, timeidx, varname, 
                                    method, squeeze, cache=None, meta=True)
            lon = xlon_var[varname]
            
            if map_proj == 1:
                if((fabs(true_lat1 - true_lat2) > 0.1) and
                        (fabs(true_lat2 - 90.) > 0.1)): 
                    cone = (log(cos(true_lat1*radians_per_degree)) 
                        - log(cos(true_lat2*radians_per_degree)))
                    cone = (cone / 
                            (log(tan((45.-fabs(true_lat1/2.))*radians_per_degree)) 
                        - log(tan((45.-fabs(true_lat2/2.))*radians_per_degree)))) 
                else:
                    cone = sin(fabs(true_lat1)*radians_per_degree)
            else:
                cone = 1
                
        return (u, v, lat, lon, cen_lon, cone)
    
    if varname == "cloudfrac":
        from wrf.g_geoht import get_height
        vars = extract_vars(wrfnc, timeidx, ("P", "PB", "QVAPOR", "T"), 
                          method, squeeze, cache=None, meta=True)
        
        p = vars["P"]
        pb = vars["PB"]
        qv = vars["QVAPOR"]
        t = vars["T"]
        
        geoht_agl = get_height(wrfnc, timeidx, method, squeeze, 
                               cache=None, meta=True, msl=False)
        
        full_p = p + pb
        full_t = t + Constants.T_BASE
        
        tkel = tk(full_p, full_t)
        relh = rh(qv, full_p, tkel)
        
        return (geoht_agl, relh, 1, 300., 2000., 6000.)
        
        
class WRFVarsTest(ut.TestCase):
    longMessage = True
    
def make_func(varname, wrfnc, timeidx, method, squeeze, meta):
    def func(self):
        
        try:
            args = get_args(varname, wrfnc, timeidx, method, squeeze)
        except ProjectionError: # Don't fail for this
            return
        
        routine = ROUTINE_MAP[varname]
        
        kwargs = {"meta" : meta}
        result = routine(*args, **kwargs)
        
        ref = getvar(wrfnc, varname, timeidx, method, squeeze, cache=None, 
                     meta=meta)
        
        nt.assert_allclose(to_np(result), to_np(ref))
        
        if meta:
            self.assertEqual(result.dims, ref.dims)
        
    return func


def test_cape3d_1d(wrfnc):
    
    def func(self):
        varnames = ("T", "P", "PB", "QVAPOR", "PH", "PHB", "HGT", "PSFC")
        ncvars = extract_vars(wrfnc, 0, varnames, method="cat", squeeze=True, 
                              cache=None, meta=True)
        
        t = ncvars["T"]
        p = ncvars["P"]
        pb = ncvars["PB"]
        qv = ncvars["QVAPOR"]
        ph = ncvars["PH"]
        phb = ncvars["PHB"]
        ter = ncvars["HGT"]
        psfc = ncvars["PSFC"]
        
        col_idxs = (slice(None), t.shape[-2]//2, t.shape[-1]//2)
        
        t = t[col_idxs]
        p = p[col_idxs]
        pb = pb[col_idxs]
        qv = qv[col_idxs]
        ph = ph[col_idxs]
        phb = phb[col_idxs]
        ter = float(ter[col_idxs[1:]])
        psfc = float(psfc[col_idxs[1:]])
        
        full_t = t + Constants.T_BASE
        full_p = p + pb
        tkel = tk(full_p, full_t, meta=False)
        
        geopt = ph + phb
        geopt_unstag = destagger(geopt, -1)
        z = geopt_unstag/Constants.G
        
        # Convert pressure to hPa
        p_hpa = ConversionFactors.PA_TO_HPA * full_p
        psfc_hpa = ConversionFactors.PA_TO_HPA * psfc 
        
        i3dflag = 1
        ter_follow = 1
        
        result = cape_3d(p_hpa, tkel, qv, z, ter, psfc_hpa, ter_follow)
    
        ref = getvar(wrfnc, "cape_3d")

        ref = ref[(slice(None),) + col_idxs]
        
        nt.assert_allclose(to_np(result), to_np(ref))
    
    return func
        
        
if __name__ == "__main__":
    varnames = ["avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
            "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
            "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
            "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
            "wa", "uvmet10", "uvmet", "z", "cloudfrac"]
    
    #varnames = ["helicity"]
    varnames=["avo", "pvo", "eth", "dbz", "helicity", "updraft_helicity",
              "omg", "pw", "rh", "slp", "td", "tk", "tv", "twb", "uvmet",
              "cloudfrac"]
    
    omp_set_num_threads(omp_get_num_procs()-1)
    omp_set_schedule(OMP_SCHED_STATIC, 0)
    omp_set_dynamic(False)
    
    # Turn this one off when not needed, since it's slow
    #varnames += ["cape_2d", "cape_3d"]
    
    for varname in varnames:
        for i,wrfnc in enumerate((NCFILE, NCGROUP)):
            for j,timeidx in enumerate((0, ALL_TIMES)):
                for method in ("cat", "join"):
                    for squeeze in (True, False):
                        for meta in (True, False):
                            func = make_func(varname, wrfnc, timeidx, method, 
                                                  squeeze, meta)
                            ncname = "single" if i == 0 else "multi"
                            timename = "t0" if j == 0 else "all"
                            squeeze_name = "squeeze" if squeeze else "nosqueeze"
                            meta_name = "meta" if meta else "nometa"
                            test_name = "test_{}_{}_{}_{}_{}_{}".format(varname, 
                                                        ncname, timename, method,
                                                        squeeze_name, meta_name)
    
                            setattr(WRFVarsTest, test_name, func)
    
    func = test_cape3d_1d(wrfnc)
    setattr(WRFVarsTest, "test_cape3d_1d", func)
    
    
    ut.main()
    
    