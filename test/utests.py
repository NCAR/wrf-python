import unittest as ut
import numpy.testing as nt 
import numpy as np
import numpy.ma as ma
import os, sys
import subprocess
import glob

from wrf import (getvar, interplevel, interpline, vertcross, vinterp,
                 disable_xarray, xarray_enabled, to_np,
                 xy_to_ll, ll_to_xy, xy_to_ll_proj, ll_to_xy_proj,
                 extract_global_attrs, viewitems, CoordPair, ll_points)
from wrf.util import is_multi_file

NCL_EXE = "/Users/ladwig/miniconda2/envs/ncl_build/bin/ncl"
NCARG_ROOT = "/Users/ladwig/miniconda2/envs/ncl_build"
#TEST_FILE = "/Users/ladwig/Documents/wrf_files/wrfout_d01_2010-06-13_21:00:00"
DIRS = ["/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/moving_nest",
        "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/static_nest"]
PATTERN = "wrfout_d02_*"
REF_NC_FILES = ["/tmp/wrftest_moving.nc", "/tmp/wrftest_static.nc"]
NEST = ["moving", "static"]

# Python 3
if sys.version_info > (3,):
    xrange = range

def setUpModule():
    #ncarg_root = os.environ.get("NCARG_ROOT", None)
    #if ncarg_root is None:
    #    raise RuntimeError("$NCARG_ROOT environment variable not set")
    
    os.environ["NCARG_ROOT"] = NCARG_ROOT
    os.environ["NCARG_NCARG"] = os.path.join(NCARG_ROOT, "lib", "ncarg")
    os.environ["OMP_NUM_THREADS"] = "4"
     
        
    this_path = os.path.realpath(__file__)
    ncl_script = os.path.join(os.path.dirname(this_path),
                              "ncl_get_var.ncl")
    
    for dir,outfile in zip(DIRS, REF_NC_FILES):
        cmd = "%s %s 'dir=\"%s\"' 'pattern=\"%s\"'  'out_file=\"%s\"'" % (
            NCL_EXE,
            ncl_script,
            dir,
            PATTERN,
            outfile)
    
        print(cmd)
    
        if not os.path.exists(outfile):
            status = subprocess.call(cmd, shell=True)
            if (status != 0):
                raise RuntimeError("NCL script failed. Could not set up test.")

# Using helpful information at: 
# http://eli.thegreenplace.net/2014/04/02/dynamically-generating-python-test-cases
def make_test(varname, dir, pattern, referent, multi=False, pynio=False):
    def test(self):

        try:
            from netCDF4 import Dataset as NetCDF
        except:
            pass
        
        try:
            import Nio
        except:
            pass
        
        timeidx = 0 if not multi else None
        pat = os.path.join(dir, pattern)
        wrf_files = glob.glob(pat)
        wrf_files.sort()
        
        wrfin = []
        for fname in wrf_files:
            if not pynio:
                f = NetCDF(fname)
                try:
                    f.set_always_mask(False)
                except:
                    pass
                wrfin.append(f)
            else:
                if not fname.endswith(".nc"):
                    _fname = fname + ".nc"
                else:
                    _fname = fname
                f = Nio.open_file(_fname)
                wrfin.append(f)

        refnc = NetCDF(referent)
        try:
            refnc.set_auto_mask(False)
        except:
            pass
        
        # These have a left index that defines the product type
        multiproduct = varname in ("uvmet", "uvmet10", "cape_2d", "cape_3d", 
                                   "cfrac")
        multi2d = ("uvmet10", "cape_2d", "cfrac")
        multi3d = ("uvmet", "cape_3d")
        
        # These varnames don't have NCL functions to test against
        ignore_referent = ("zstag", "geopt_stag")
        
        if varname not in ignore_referent:
            if not multi:
                if varname in multi2d:
                    ref_vals = refnc.variables[varname][...,0,:,:]
                elif varname in multi3d:
                    ref_vals = refnc.variables[varname][...,0,:,:,:]
                else:
                    ref_vals = refnc.variables[varname][0,:]
            else:
                ref_vals = refnc.variables[varname][:]
        
        if (varname == "tc"):
            my_vals = getvar(wrfin, "temp", timeidx=timeidx, units="c")
            tol = 1/100.
            atol = .1 # Note:  NCL uses 273.16 as conversion for some reason
            nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
        elif (varname == "height_agl"):
            # Change the vert_type to height_agl when NCL gets updated.
            my_vals = getvar(wrfin, "z", timeidx=timeidx, msl=False)
            tol = 1/100.
            atol = .1 # Note:  NCL uses 273.16 as conversion for some reason
            nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
        elif (varname == "cfrac"):
            # Change the vert_type to height_agl when NCL gets updated.
            my_vals = getvar(wrfin, "cfrac", timeidx=timeidx)
            tol = 1/100.
            atol = .1 # Note:  NCL uses 273.16 as conversion for some reason
            nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
        elif (varname == "pw"):
            my_vals = getvar(wrfin, "pw", timeidx=timeidx)
            tol = .5/100.0
            atol = 0 # NCL uses different constants and doesn't use same
                     # handrolled virtual temp in method
            try:
                nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
            except AssertionError:
                print (np.amax(np.abs(to_np(my_vals) - ref_vals)))
                raise
        elif (varname == "cape_2d"):
            cape_2d = getvar(wrfin, varname, timeidx=timeidx)
            tol = 0/100.
            atol = 200.0
            # Let's only compare CAPE values until the F90 changes are 
            # merged back in to NCL.  The modifications to the R and CP
            # changes TK enough that non-lifting parcels could lift, thus
            # causing wildly different values in LCL
            nt.assert_allclose(to_np(cape_2d[0,:]), ref_vals[0,:], tol, atol)
        elif (varname == "cape_3d"):
            cape_3d = getvar(wrfin, varname, timeidx=timeidx)
            # Changing the R and CP constants, while keeping TK within
            # 2%, can lead to some big changes in CAPE.  Tolerances 
            # have been set wide when comparing the with the original
            # NCL.  Change back when the F90 code is merged back with 
            # NCL
            tol = 0/100.
            atol = 200.0
            
            #print np.amax(np.abs(to_np(cape_3d[0,:]) - ref_vals[0,:]))
            nt.assert_allclose(to_np(cape_3d), ref_vals, tol, atol)
        elif (varname == "zstag" or varname == "geopt_stag"):
            v = getvar(wrfin, varname, timeidx=timeidx)
            # For now, only make sure it runs without crashing since no NCL
            # to compare with yet.
        else:
            my_vals = getvar(wrfin, varname, timeidx=timeidx)
            tol = 2/100.
            atol = 0.1
            #print (np.amax(np.abs(to_np(my_vals) - ref_vals)))
            try:
                nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
            except:
                absdiff = np.abs(to_np(my_vals) - ref_vals)
                maxdiff = np.amax(absdiff)
                print(maxdiff)
                print(np.argwhere(absdiff == maxdiff))
                
                raise
    
    
    return test

def _get_refvals(referent, varname, multi):
    try:
        from netCDF4 import Dataset as NetCDF
    except:
        pass
        
    refnc = NetCDF(referent)
    try:
        pass
        #refnc.set_auto_mask(False)
    except:
        pass
    
    multi2d = ("uvmet10", "cape_2d", "cfrac")
    multi3d = ("uvmet", "cape_3d")
    interpline = ("t2_line","t2_line2", "t2_line3")
    
    if not multi:
        if varname in multi2d:
            ref_vals = refnc.variables[varname][...,0,:,:]
        elif varname in multi3d:
            ref_vals = refnc.variables[varname][...,0,:,:,:]
        else:
            v = refnc.variables[varname][:]
            if v.ndim == 2:
                if varname in interpline:
                    ref_vals = v[0,:]
                else:
                    ref_vals = v
            else:
                ref_vals = v[0,:]
    else:
        ref_vals = refnc.variables[varname][:]
                
    return ref_vals

def make_interp_test(varname, dir, pattern, referent, multi=False, 
                     pynio=False):
    def test(self):
        try:
            from netCDF4 import Dataset as NetCDF
        except:
            pass
        
        try:
            import Nio
        except:
            pass
        
        timeidx = 0 if not multi else None
        pat = os.path.join(dir, pattern)
        wrf_files = glob.glob(pat)
        wrf_files.sort()
        
        wrfin = []
        for fname in wrf_files:
            if not pynio:
                f = NetCDF(fname)
                try:
                    f.set_always_mask(False)
                except:
                    pass
                wrfin.append(f)
            else:
                if not fname.endswith(".nc"):
                    _fname = fname + ".nc"
                else:
                    _fname = fname
                f = Nio.open_file(_fname)
                wrfin.append(f)
        
        if (varname == "interplevel"):
            ref_ht_500 = _get_refvals(referent, "z_500", multi)
            ref_p_5000 = _get_refvals(referent, "p_5000", multi)
            ref_ht_multi = _get_refvals(referent, "z_multi", multi)
            ref_p_multi = _get_refvals(referent, "p_multi", multi)
            
            ref_ht2_500 = _get_refvals(referent, "z2_500", multi)
            ref_p2_5000 = _get_refvals(referent, "p2_5000", multi)
            ref_ht2_multi = _get_refvals(referent, "z2_multi", multi)
            ref_p2_multi = _get_refvals(referent, "p2_multi", multi)
            
            ref_p_lev2d = _get_refvals(referent, "p_lev2d", multi)
            
            hts = getvar(wrfin, "z", timeidx=timeidx)
            p = getvar(wrfin, "pressure", timeidx=timeidx)
            wspd_wdir = getvar(wrfin, "wspd_wdir", timeidx=timeidx)
            
            # Make sure the numpy versions work first
            hts_500 = interplevel(to_np(hts), to_np(p), 500)
            hts_500 = interplevel(hts, p, 500)
            
            # Note: the '*2*' versions in the reference are testing 
            # against the new version of interplevel in NCL
            nt.assert_allclose(to_np(hts_500), ref_ht_500)
            nt.assert_allclose(to_np(hts_500), ref_ht2_500)
            
            # Make sure the numpy versions work first
            p_5000 = interplevel(to_np(p), to_np(hts), 5000)
            p_5000 = interplevel(p, hts, 5000)
            
            
            nt.assert_allclose(to_np(p_5000), ref_p_5000)
            nt.assert_allclose(to_np(p_5000), ref_p2_5000)
            
            hts_multi= interplevel(to_np(hts), to_np(p), 
                                   [1000., 850., 500., 250.])
            hts_multi = interplevel(hts, p, [1000., 850., 500., 250.])
            
            nt.assert_allclose(to_np(hts_multi), ref_ht_multi)
            nt.assert_allclose(to_np(hts_multi), ref_ht2_multi)
            
            p_multi= interplevel(to_np(p), to_np(hts), 
                                   [500., 2500., 5000., 10000. ])
            p_multi = interplevel(p, hts, [500., 2500., 5000., 10000. ])
            
            nt.assert_allclose(to_np(p_multi), ref_p_multi)
            nt.assert_allclose(to_np(p_multi), ref_p2_multi)
            
            pblh = getvar(wrfin, "PBLH", timeidx=timeidx)
            p_lev2d = interplevel(to_np(p), to_np(hts), to_np(pblh))
            p_lev2d = interplevel(p, hts, pblh)
            nt.assert_allclose(to_np(p_lev2d), ref_p_lev2d)
            
            # Just make sure these run below
            wspd_wdir_500 = interplevel(to_np(wspd_wdir), to_np(p), 500)
            wspd_wdir_500 = interplevel(wspd_wdir, p, 500)
            #print(wspd_wdir_500)
            
            wspd_wdir_multi= interplevel(to_np(wspd_wdir), 
                                         to_np(p), [1000,500,250])
            wdpd_wdir_multi = interplevel(wspd_wdir, p, [1000,500,250])
            
            
            wspd_wdir_pblh = interplevel(to_np(wspd_wdir), to_np(hts), pblh)
            wspd_wdir_pblh = interplevel(wspd_wdir, hts, pblh)
            
            if multi:
                wspd_wdir_pblh_2 = interplevel(to_np(wspd_wdir), 
                                         to_np(hts), pblh[0,:])
                wspd_wdir_pblh_2 = interplevel(wspd_wdir, hts, pblh[0,:])
                
                # Since PBLH doesn't change in this case, it should match
                # the 0 time from previous computation. Note that this 
                # only works when the data has 1 time step that is repeated.
                # If you use a different case with multiple times, 
                # this will probably fail.
                nt.assert_allclose(to_np(wspd_wdir_pblh_2[:,0,:]), 
                                   to_np(wspd_wdir_pblh[:,0,:]))
                
                nt.assert_allclose(to_np(wspd_wdir_pblh_2[:,-1,:]), 
                                   to_np(wspd_wdir_pblh[:,0,:]))
                
            
        elif (varname == "vertcross"):
            ref_ht_cross = _get_refvals(referent, "ht_cross", multi)
            ref_p_cross = _get_refvals(referent, "p_cross", multi)
            ref_ht_vertcross1 = _get_refvals(referent, "ht_vertcross1", 
                                            multi)
            ref_ht_vertcross2 = _get_refvals(referent, "ht_vertcross2", 
                                            multi)
            ref_ht_vertcross3 = _get_refvals(referent, "ht_vertcross3", 
                                            multi)
            
            hts = getvar(wrfin, "z", timeidx=timeidx)
            p = getvar(wrfin, "pressure", timeidx=timeidx)
            
            pivot_point = CoordPair(hts.shape[-1] / 2, hts.shape[-2] / 2) 
            
            # Beginning in wrf-python 1.3, users can select number of levels.
            # Originally, for pressure, dz was 10, so let's back calculate
            # the number of levels.
            p_max = np.floor(np.amax(p)/10) * 10     # bottom value
            p_min = np.ceil(np.amin(p)/10) * 10      # top value
            
            p_autolevels = int((p_max - p_min) /10)
            
            # Make sure the numpy versions work first
            
            ht_cross = vertcross(to_np(hts), to_np(p), 
                                 pivot_point=pivot_point, angle=90.,
                                 autolevels=p_autolevels)
            
            ht_cross = vertcross(hts, p, pivot_point=pivot_point, angle=90.,
                                 autolevels=p_autolevels)
            
            nt.assert_allclose(to_np(ht_cross), ref_ht_cross, rtol=.01)
            
            lats = hts.coords["XLAT"]
            lons = hts.coords["XLONG"]
                
            # Test the manual projection method with lat/lon
            # Only do this for the non-multi case, since the domain 
            # might be moving
            if not multi:
                if lats.ndim > 2: # moving nest
                    lats = lats[0,:]
                    lons = lons[0,:]
                    
                ll_point = ll_points(lats, lons)
            
                pivot = CoordPair(lat=lats[int(lats.shape[-2]/2), 
                                       int(lats.shape[-1]/2)],
                                  lon=lons[int(lons.shape[-2]/2), 
                                       int(lons.shape[-1]/2)])
            
                v1 = vertcross(hts,p,wrfin=wrfin,pivot_point=pivot_point,
                           angle=90.0)
                v2 = vertcross(hts,p,projection=hts.attrs["projection"], 
                          ll_point=ll_point,
                          pivot_point=pivot_point, angle=90.)
            
                nt.assert_allclose(to_np(v1), to_np(v2), rtol=.01)
            
            # Test opposite
            
            p_cross1 = vertcross(p,hts,pivot_point=pivot_point, angle=90.0)
            
            nt.assert_allclose(to_np(p_cross1), 
                               ref_p_cross, 
                               rtol=.01)
            # Test point to point
            start_point = CoordPair(0, hts.shape[-2]/2)
            end_point = CoordPair(-1, hts.shape[-2]/2)
          
             
            p_cross2 = vertcross(p,hts,start_point=start_point, 
                                end_point=end_point)
            
            nt.assert_allclose(to_np(p_cross1), 
                               to_np(p_cross2))
            
            # Check the new vertcross routine
            pivot_point = CoordPair(hts.shape[-1] / 2, hts.shape[-2] / 2) 
            ht_cross = vertcross(hts, p, 
                                 pivot_point=pivot_point, angle=90.,
                                 latlon=True)
            
            nt.assert_allclose(to_np(ht_cross), 
                               to_np(ref_ht_vertcross1), atol=.01)
            
            
            
            levels = [1000., 850., 700., 500., 250.]
            ht_cross = vertcross(hts, p, 
                                 pivot_point=pivot_point, angle=90.,
                                 levels=levels, latlon=True)
            
            nt.assert_allclose(to_np(ht_cross), 
                               to_np(ref_ht_vertcross2), atol=.01)
            
            idxs = (0, slice(None)) if lats.ndim > 2 else (slice(None),)
            
            start_lat = np.amin(lats[idxs]) + .25*(np.amax(lats[idxs]) - np.amin(lats[idxs]))
            end_lat = np.amin(lats[idxs]) + .65*(np.amax(lats[idxs]) - np.amin(lats[idxs]))
            
            start_lon = np.amin(lons[idxs]) + .25*(np.amax(lons[idxs]) - np.amin(lons[idxs]))
            end_lon = np.amin(lons[idxs]) + .65*(np.amax(lons[idxs]) - np.amin(lons[idxs]))
            
            start_point = CoordPair(lat=start_lat, lon=start_lon)
            end_point = CoordPair(lat=end_lat, lon=end_lon)
            
            ll_point = ll_points(lats[idxs], lons[idxs])
            
            ht_cross = vertcross(hts, p, 
                                 start_point=start_point, 
                                 end_point=end_point, 
                                 projection=hts.attrs["projection"], 
                                 ll_point=ll_point, 
                                 latlon=True,
                                 autolevels=1000)
            
            
            nt.assert_allclose(to_np(ht_cross), 
                               to_np(ref_ht_vertcross3),
                               rtol=.01)
            
            if multi:
                ntimes = hts.shape[0]
                
                for t in range(ntimes):
                    hts = getvar(wrfin, "z", timeidx=t)
                    p = getvar(wrfin, "pressure", timeidx=t)
                    
                    ht_cross = vertcross(hts, p, 
                                 start_point=start_point, 
                                 end_point=end_point, 
                                 wrfin=wrfin,
                                 timeidx=t,
                                 latlon=True,
                                 autolevels=1000)
                    
                    refname = "ht_vertcross_t{}".format(t)
                    ref_ht_vertcross = _get_refvals(referent, refname, False)
                    
                    nt.assert_allclose(to_np(ht_cross), 
                               to_np(ref_ht_vertcross),rtol=.02)
            
              
        elif (varname == "interpline"):
            
            ref_t2_line = _get_refvals(referent, "t2_line", multi)
            ref_t2_line2 = _get_refvals(referent, "t2_line2", multi)
            ref_t2_line3 = _get_refvals(referent, "t2_line3", multi)
            
            t2 = getvar(wrfin, "T2", timeidx=timeidx)
            pivot_point = CoordPair(t2.shape[-1] / 2, t2.shape[-2] / 2)
            
            # Make sure the numpy version works
            t2_line1 = interpline(to_np(t2), pivot_point=pivot_point, 
                                  angle=90.0)
            t2_line1 = interpline(t2, pivot_point=pivot_point, angle=90.0)
            
            nt.assert_allclose(to_np(t2_line1), ref_t2_line)
            
            # Test the new NCL wrf_user_interplevel result
            nt.assert_allclose(to_np(t2_line1), ref_t2_line2)
            
            # Test the manual projection method with lat/lon
            lats = t2.coords["XLAT"]
            lons = t2.coords["XLONG"]
            if multi:
                if lats.ndim > 2: # moving nest
                    lats = lats[0,:]
                    lons = lons[0,:]
            
            ll_point = ll_points(lats, lons)
            
            pivot = CoordPair(lat=lats[int(lats.shape[-2]/2), 
                                       int(lats.shape[-1]/2)],
                              lon=lons[int(lons.shape[-2]/2), 
                                       int(lons.shape[-1]/2)])

            l1 = interpline(t2,wrfin=wrfin,pivot_point=pivot_point,
                           angle=90.0)
            
            l2 = interpline(t2,projection=t2.attrs["projection"], 
                          ll_point=ll_point, 
                          pivot_point=pivot_point, angle=90.)
            nt.assert_allclose(to_np(l1), to_np(l2), rtol=.01)
            
            # Test point to point
            start_point = CoordPair(0, t2.shape[-2]/2)
            end_point = CoordPair(-1, t2.shape[-2]/2)
            
            t2_line2 = interpline(t2, start_point=start_point, 
                                  end_point=end_point)
            
            nt.assert_allclose(to_np(t2_line1), to_np(t2_line2))
            
            # Now test the start/end with lat/lon points
            
            start_lat = float(np.amin(lats) + .25*(np.amax(lats) 
                                                   - np.amin(lats)))
            end_lat = float(np.amin(lats) + .65*(np.amax(lats) 
                                                 - np.amin(lats)))
            
            start_lon = float(np.amin(lons) + .25*(np.amax(lons) 
                                                   - np.amin(lons)))
            end_lon = float(np.amin(lons) + .65*(np.amax(lons) 
                                                 - np.amin(lons)))
            
            start_point = CoordPair(lat=start_lat, lon=start_lon)
            end_point = CoordPair(lat=end_lat, lon=end_lon)
            
            t2_line3 = interpline(t2,wrfin=wrfin,timeidx=0,
                                  start_point=start_point,
                                  end_point=end_point,latlon=True)
            
            
            nt.assert_allclose(to_np(t2_line3), ref_t2_line3, rtol=.01)
            
            # Test all time steps
            if multi:
                refnc = NetCDF(referent)
                ntimes = t2.shape[0]
                
                for t in range(ntimes):
                    t2 = getvar(wrfin, "T2", timeidx=t)
                    
                    line = interpline(t2,wrfin=wrfin,timeidx=t,
                                  start_point=start_point,
                                  end_point=end_point,latlon=True)
                    
                    refname = "t2_line_t{}".format(t)
                    refline = refnc.variables[refname][:]
                    
                    nt.assert_allclose(to_np(line), 
                               to_np(refline),rtol=.005)
            
                refnc.close()
            
            # Test NCLs single time case
            if not multi:
                refnc = NetCDF(referent)
                ref_t2_line4 = refnc.variables["t2_line4"][:]
                
                t2 = getvar(wrfin, "T2", timeidx=0)
                line = interpline(t2,wrfin=wrfin,timeidx=0,
                                  start_point=start_point,
                                  end_point=end_point,latlon=True)
                
                nt.assert_allclose(to_np(line), 
                               to_np(ref_t2_line4),rtol=.005)
                refnc.close()

        elif (varname == "vinterp"):
            # Tk to theta
            fld_tk_theta = _get_refvals(referent, "fld_tk_theta", multi)
            fld_tk_theta = np.squeeze(fld_tk_theta)
            
            tk = getvar(wrfin, "temp", timeidx=timeidx, units="k")
            
            interp_levels = [200,300,500,1000]
            
            # Make sure the numpy version works
            field = vinterp(wrfin, 
                            field=to_np(tk), 
                            vert_coord="theta", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk",
                            timeidx=timeidx, 
                            log_p=True)
            
            field = vinterp(wrfin, 
                            field=tk, 
                            vert_coord="theta", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk",
                            timeidx=timeidx, 
                            log_p=True)
            
            tol = 5/100.
            atol = 0.0001
            
            field = np.squeeze(field)
            #print (np.amax(np.abs(to_np(field) - fld_tk_theta)))
            nt.assert_allclose(to_np(field), fld_tk_theta, tol, atol)
            
            # Tk to theta-e
            fld_tk_theta_e = _get_refvals(referent, "fld_tk_theta_e", multi)
            fld_tk_theta_e = np.squeeze(fld_tk_theta_e)
            
            interp_levels = [200,300,500,1000]
            
            field = vinterp(wrfin, 
                            field=tk, 
                            vert_coord="theta-e", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk", 
                            timeidx=timeidx, 
                            log_p=True)
            
            tol = 3/100.
            atol = 50.0001
            
            field = np.squeeze(field)
            #print (np.amax(np.abs(to_np(field) - fld_tk_theta_e)/fld_tk_theta_e)*100)
            nt.assert_allclose(to_np(field), fld_tk_theta_e, tol, atol)
            
            # Tk to pressure
            fld_tk_pres = _get_refvals(referent, "fld_tk_pres", multi)
            fld_tk_pres = np.squeeze(fld_tk_pres)
            
            interp_levels = [850,500]
            
            field = vinterp(wrfin, 
                            field=tk, 
                            vert_coord="pressure", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk",
                            timeidx=timeidx,  
                            log_p=True)
            
            field = np.squeeze(field)
            
            #print (np.amax(np.abs(to_np(field) - fld_tk_pres)))
            nt.assert_allclose(to_np(field), fld_tk_pres, tol, atol)
            
            # Tk to geoht_msl
            fld_tk_ght_msl = _get_refvals(referent, "fld_tk_ght_msl", multi)
            fld_tk_ght_msl = np.squeeze(fld_tk_ght_msl)
            interp_levels = [1,2]
            
            field = vinterp(wrfin, 
                            field=tk, 
                            vert_coord="ght_msl", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk",
                            timeidx=timeidx,  
                            log_p=True)
            
            field = np.squeeze(field)
            #print (np.amax(np.abs(to_np(field) - fld_tk_ght_msl)))
            nt.assert_allclose(to_np(field), fld_tk_ght_msl, tol, atol)
            
            # Tk to geoht_agl
            fld_tk_ght_agl = _get_refvals(referent, "fld_tk_ght_agl", multi)
            fld_tk_ght_agl = np.squeeze(fld_tk_ght_agl)
            interp_levels = [1,2]
            
            field = vinterp(wrfin, 
                            field=tk, 
                            vert_coord="ght_agl", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk", 
                            timeidx=timeidx, 
                            log_p=True)
            
            field = np.squeeze(field)
            #print (np.amax(np.abs(to_np(field) - fld_tk_ght_agl)))
            nt.assert_allclose(to_np(field), fld_tk_ght_agl, tol, atol)
            
            # Hgt to pressure
            fld_ht_pres = _get_refvals(referent, "fld_ht_pres", multi)
            fld_ht_pres = np.squeeze(fld_ht_pres)
            
            z = getvar(wrfin, "height", timeidx=timeidx, units="m")
            interp_levels = [500,50]
            field = vinterp(wrfin, 
                            field=z, 
                            vert_coord="pressure", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="ght", 
                            timeidx=timeidx, 
                            log_p=True)
            
            field = np.squeeze(field)
            #print (np.amax(np.abs(to_np(field) - fld_ht_pres)))
            nt.assert_allclose(to_np(field), fld_ht_pres, tol, atol)
            
            # Pressure to theta
            fld_pres_theta = _get_refvals(referent, "fld_pres_theta", multi)
            fld_pres_theta = np.squeeze(fld_pres_theta)
            
            p = getvar(wrfin, "pressure", timeidx=timeidx)
            interp_levels = [200,300,500,1000]
            field = vinterp(wrfin, 
                            field=p, 
                            vert_coord="theta", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="pressure", 
                            timeidx=timeidx, 
                            log_p=True)
            
            field = np.squeeze(field)
            #print (np.amax(np.abs(to_np(field) - fld_pres_theta)))
            nt.assert_allclose(to_np(field), fld_pres_theta, tol, atol)
            
            # Theta-e to pres
            fld_thetae_pres = _get_refvals(referent, "fld_thetae_pres", multi)
            fld_thetae_pres = np.squeeze(fld_thetae_pres)
            
            eth = getvar(wrfin, "eth", timeidx=timeidx)
            interp_levels = [850,500,5]
            field = vinterp(wrfin, 
                            field=eth, 
                            vert_coord="pressure", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="theta-e", 
                            timeidx=timeidx, 
                            log_p=True)
            
            field = np.squeeze(field)
            #print (np.amax(np.abs(to_np(field) - fld_thetae_pres)))
            nt.assert_allclose(to_np(field), fld_thetae_pres, tol, atol)
    
    return test

def extract_proj_params(wrfnc, timeidx=0):
    attrs = extract_global_attrs(wrfnc, ("MAP_PROJ", "TRUELAT1", "TRUELAT2",
                                         "STAND_LON", "POLE_LAT", "POLE_LON", 
                                         "DX", "DY"))
    
    result = {key.lower(): val for key,val in viewitems(attrs)}
    
    _timeidx = timeidx
    if is_multi_file(wrfnc):
        wrfnc0 = wrfnc[0]
        num_times_per_file = len(wrfnc0.dimensions["Time"])
        file_idx = timeidx // num_times_per_file
        _timeidx = timeidx % num_times_per_file
        
        wrfnc = wrfnc[file_idx]
    
    result["known_x"] = 0
    result["known_y"] = 0
    result["ref_lat"] = wrfnc.variables["XLAT"][_timeidx,0,0]
    result["ref_lon"] = wrfnc.variables["XLONG"][_timeidx,0,0]
    
    return result

def make_latlon_test(testid, dir, pattern, referent, single, 
                     multi=False, pynio=False):
    def test(self):
        try:
            from netCDF4 import Dataset as NetCDF
        except:
            pass
        
        try:
            import Nio
        except:
            pass
        
        timeidx = 0 if not multi else None
        pat = os.path.join(dir, pattern)
        wrf_files = glob.glob(pat)
        wrf_files.sort()
        
        refnc = NetCDF(referent)
        try:
            refnc.set_always_mask(False)
        except:
            pass
        
        wrfin = []
        for fname in wrf_files:
            if not pynio:
                f = NetCDF(fname)
                try:
                    f.set_auto_mask(False)
                except:
                    pass
                wrfin.append(f)
            else:
                if not fname.endswith(".nc"):
                    _fname = fname + ".nc"
                else:
                    _fname = fname
                f = Nio.open_file(_fname)
                wrfin.append(f)
                
        if testid == "xy":
            
            # Lats/Lons taken from NCL script, just hard-coding for now
            lats = [22.0, 25.0, 27.0]
            lons = [-90.0, -87.5, -83.75]
            
            # Just call with a single lat/lon
            if single:
                timeidx = 8
                ref_vals = refnc.variables["xy2"][:]
                
                xy = ll_to_xy(wrfin, lats[0], lons[0], timeidx=timeidx,
                              as_int=True)
                ref = ref_vals[:,0]
                
                nt.assert_allclose(to_np(xy), ref)
                
                # Next make sure the 'proj' version works
                projparams = extract_proj_params(wrfin, timeidx=timeidx)
                xy_proj = ll_to_xy_proj(lats[0], lons[0],
                                        as_int=True, 
                                        **projparams)
                
                nt.assert_allclose(to_np(xy_proj), to_np(xy))
                
            
            else:
                ref_vals = refnc.variables["xy1"][:]
                xy = ll_to_xy(wrfin, lats, lons, timeidx=None, as_int=False)
                
                ref = ref_vals[:]
                
                nt.assert_allclose(to_np(xy), ref)
                
                if xy.ndim > 2:
                    # Moving nest
                    is_moving = True
                    numtimes = xy.shape[-2]
                else:
                    is_moving = False
                    numtimes = 1
                
                for tidx in range(9):
                
                    # Next make sure the 'proj' version works
                    projparams = extract_proj_params(wrfin, timeidx=tidx)
                    xy_proj = ll_to_xy_proj(lats, lons, as_int=False,
                                            **projparams)
                    
                    if is_moving:
                        idxs = (slice(None), tidx, slice(None))
                    else:
                        idxs = (slice(None),)
                
                    nt.assert_allclose(to_np(xy_proj), to_np(xy[idxs]))
                
        else:
             # i_s, j_s taken from NCL script, just hard-coding for now
             # NCL uses 1-based indexing for this, so need to subtract 1
            x_s = np.asarray([10, 50, 90], int)
            y_s = np.asarray([10, 50, 90], int)
            
            if single:
                timeidx = 8
                ref_vals = refnc.variables["ll2"][:]
                ll = xy_to_ll(wrfin, x_s[0], y_s[0], timeidx=timeidx)
                ref = ref_vals[::-1,0]
                
                nt.assert_allclose(to_np(ll), ref)
                
                # Next make sure the 'proj' version works
                projparams = extract_proj_params(wrfin, timeidx=8)
                ll_proj = xy_to_ll_proj(x_s[0], y_s[0], **projparams)
                
                nt.assert_allclose(to_np(ll_proj), to_np(ll))
                 
                
            else:
                ref_vals = refnc.variables["ll1"][:]
                ll = xy_to_ll(wrfin, x_s, y_s, timeidx=None)
                ref = ref_vals[::-1,:]
                
                nt.assert_allclose(to_np(ll), ref)
                
                if ll.ndim > 2:
                    # Moving nest
                    is_moving = True
                    numtimes = ll.shape[-2]
                else:
                    is_moving = False
                    numtimes = 1
                
                for tidx in range(numtimes):
                    # Next make sure the 'proj' version works
                    projparams = extract_proj_params(wrfin, timeidx=tidx)
                    ll_proj = xy_to_ll_proj(x_s, y_s, **projparams)
                    
                    if is_moving:
                        idxs = (slice(None), tidx, slice(None))
                    else:
                        idxs = (slice(None),)
                        
                    nt.assert_allclose(to_np(ll_proj), to_np(ll[idxs]))
        
    return test

class WRFVarsTest(ut.TestCase):
    longMessage = True
    
class WRFInterpTest(ut.TestCase):
    longMessage = True
    
class WRFLatLonTest(ut.TestCase):
    longMessage = True
        

if __name__ == "__main__":
    from wrf import (omp_set_num_threads, omp_set_schedule, omp_get_schedule, 
                     omp_set_dynamic, omp_get_num_procs, OMP_SCHED_STATIC)
    omp_set_num_threads(omp_get_num_procs()//2)
    omp_set_schedule(OMP_SCHED_STATIC, 0)
    omp_set_dynamic(False)
    
    ignore_vars = []  # Not testable yet
    wrf_vars = ["avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
                "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
                "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
                "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
                "wa", "uvmet10", "uvmet", "z", "cfrac", "zstag", "geopt_stag",
                "height_agl"]
    interp_methods = ["interplevel", "vertcross", "interpline", "vinterp"]
    latlon_tests = ["xy", "ll"]
    
    for dir, ref_nc_file, nest in zip(DIRS, REF_NC_FILES, NEST):
        try:
            import netCDF4
        except ImportError:
            pass
        else:
            for var in wrf_vars:
                if var in ignore_vars:
                    continue
                
                test_func1 = make_test(var, dir, PATTERN, ref_nc_file)
                test_func2 = make_test(var, dir, PATTERN, ref_nc_file, multi=True)
                setattr(WRFVarsTest, 'test_{0}_{1}'.format(nest,var), test_func1)
                setattr(WRFVarsTest, 'test_{0}_multi_{1}'.format(nest,var), test_func2)
                
            for method in interp_methods:
                test_interp_func1 = make_interp_test(method, dir, PATTERN, 
                                                     ref_nc_file)
                test_interp_func2 = make_interp_test(method, dir, PATTERN, 
                                                     ref_nc_file, multi=True)
                setattr(WRFInterpTest, 'test_{0}_{1}'.format(nest,method), 
                        test_interp_func1)
                setattr(WRFInterpTest, 'test_{0}_multi_{1}'.format(nest,method), 
                        test_interp_func2)
            
            for testid in latlon_tests:
                for single in (True, False):
                    for multi in (True, False):
                        test_ll_func = make_latlon_test(testid, dir, PATTERN, 
                                                        ref_nc_file, 
                                                        single=single, 
                                                        multi=multi, 
                                                        pynio=False)
                        multistr = "" if not multi else "_multi"
                        singlestr = "_nosingle" if not single else "_single"
                        test_name = "test_{}_{}{}{}".format(nest, testid, singlestr, 
                                                           multistr)
                        setattr(WRFLatLonTest, test_name, test_ll_func)
                    
        try:
            import PyNIO
        except ImportError:
            pass
        else:
            for var in wrf_vars:
                if var in ignore_vars:
                    continue
                
                test_func1 = make_test(var, dir, PATTERN, ref_nc_file, pynio=True)
                test_func2 = make_test(var, dir, PATTERN, ref_nc_file, multi=True,
                                       pynio=True)
                setattr(WRFVarsTest, 'test_pynio_{0}_{1}'.format(nest,var), test_func1)
                setattr(WRFVarsTest, 'test_pynio_{0}_multi_{1}'.format(nest,var), 
                        test_func2)
                
            for method in interp_methods:
                test_interp_func1 = make_interp_test(method, dir, PATTERN, 
                                                     ref_nc_file)
                test_interp_func2 = make_interp_test(method, dir, PATTERN, 
                                                     ref_nc_file, multi=True)
                setattr(WRFInterpTest, 'test_pynio_{0}_{1}'.format(nest,method), 
                        test_interp_func1)
                setattr(WRFInterpTest, 'test_pynio_{0}_multi_{1}'.format(nest,method), 
                        test_interp_func2)
                
            for testid in latlon_tests:
                for single in (True, False):
                    for multi in (True, False):
                        test_ll_func = make_latlon_test(testid, dir, PATTERN, 
                                                        ref_nc_file, 
                                                        single=single, 
                                                        multi=multi, 
                                                        pynio=False)
                        multistr = "" if not multi else "_multi"
                        singlestr = "_nosingle" if not single else "_single"
                        test_name = "test_pynio_{}_{}{}{}".format(nest, testid, 
                                                                  singlestr, 
                                                                  multistr)
                        setattr(WRFLatLonTest, test_name, test_ll_func)
     
    ut.main()
    