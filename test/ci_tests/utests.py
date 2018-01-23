import unittest as ut
import numpy.testing as nt 
import numpy as np
import numpy.ma as ma
import os, sys
import subprocess

from wrf import (getvar, interplevel, interpline, vertcross, vinterp,
                 disable_xarray, xarray_enabled, to_np,
                 xy_to_ll, ll_to_xy, xy_to_ll_proj, ll_to_xy_proj,
                 extract_global_attrs, viewitems, CoordPair, 
                 omp_get_num_procs, omp_set_num_threads)
from wrf.util import is_multi_file

TEST_FILE = "ci_test_file.nc"
REF_FILE = "ci_result_file.nc"


# Python 3
if sys.version_info > (3,):
    xrange = range

# Using helpful information at: 
# http://eli.thegreenplace.net/2014/04/02/dynamically-generating-python-test-cases
def make_test(varname, wrf_in, referent, multi=False, repeat=3, pynio=False):
    def test(self):
        
        from netCDF4 import Dataset as NetCDF
        timeidx = 0
        in_wrfnc = NetCDF(wrf_in)
        refnc = NetCDF(referent)
        
        # These have a left index that defines the product type
        multiproduct = varname in ("uvmet", "uvmet10", "cape_2d", "cape_3d", 
                                   "cfrac")
        
        
        ref_vals = refnc.variables[varname][:]
        
        if (varname == "tc"):
            my_vals = getvar(in_wrfnc, "temp", timeidx=timeidx, units="c")
            tol = 1/100.
            atol = .1 # Note:  NCL uses 273.16 as conversion for some reason
            nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
        elif (varname == "pw"):
            my_vals = getvar(in_wrfnc, "pw", timeidx=timeidx)
            tol = .5/100.0
            atol = 0 # NCL uses different constants and doesn't use same
                     # handrolled virtual temp in method
            nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
        elif (varname == "cape_2d"):
            cape_2d = getvar(in_wrfnc, varname, timeidx=timeidx)
            tol = 0/100.
            atol = 200.0
            # Let's only compare CAPE values until the F90 changes are 
            # merged back in to NCL.  The modifications to the R and CP
            # changes TK enough that non-lifting parcels could lift, thus
            # causing wildly different values in LCL
            nt.assert_allclose(to_np(cape_2d[0,:]), ref_vals[0,:], tol, atol)
        elif (varname == "cape_3d"):
            cape_3d = getvar(in_wrfnc, varname, timeidx=timeidx)
            # Changing the R and CP constants, while keeping TK within
            # 2%, can lead to some big changes in CAPE.  Tolerances 
            # have been set wide when comparing the with the original
            # NCL.  Change back when the F90 code is merged back with 
            # NCL
            tol = 0/100.
            atol = 200.0
            
            #print np.amax(np.abs(to_np(cape_3d[0,:]) - ref_vals[0,:]))
            nt.assert_allclose(to_np(cape_3d), ref_vals, tol, atol)
        else:
            my_vals = getvar(in_wrfnc, varname, timeidx=timeidx)
            tol = 2/100.
            atol = 0.1
            #print (np.amax(np.abs(to_np(my_vals) - ref_vals)))
            nt.assert_allclose(to_np(my_vals), ref_vals, tol, atol)
    
    
    return test

def _get_refvals(referent, varname, repeat, multi):
    from netCDF4 import Dataset as NetCDF
        
    refnc = NetCDF(referent)
    
    ref_vals = refnc.variables[varname][:]
                
    return ref_vals

def make_interp_test(varname, wrf_in, referent, multi=False, 
                     repeat=3, pynio=False):
    def test(self):
        from netCDF4 import Dataset as NetCDF
        
        timeidx = 0
        in_wrfnc = NetCDF(wrf_in)
        
        if (varname == "interplevel"):
            ref_ht_850 = _get_refvals(referent, "interplevel", repeat, multi)
            hts = getvar(in_wrfnc, "z", timeidx=timeidx)
            p = getvar(in_wrfnc, "pressure", timeidx=timeidx)

            # Check that it works with numpy arrays
            hts_850 = interplevel(to_np(hts), p, 850)
            #print (hts_850)
            hts_850 = interplevel(hts, p, 850)
            
            nt.assert_allclose(to_np(hts_850), ref_ht_850)
            
        elif (varname == "vertcross"):
            ref_ht_cross = _get_refvals(referent, "vertcross", repeat, multi)
            
            hts = getvar(in_wrfnc, "z", timeidx=timeidx)
            p = getvar(in_wrfnc, "pressure", timeidx=timeidx)
            
            pivot_point = CoordPair(hts.shape[-1] // 2, hts.shape[-2] // 2) 
            # Check that it works with numpy arrays
            ht_cross = vertcross(to_np(hts), to_np(p), 
                                 pivot_point=pivot_point, angle=90.)
            #print (ht_cross)
            ht_cross = vertcross(hts, p, pivot_point=pivot_point, angle=90.)

            nt.assert_allclose(to_np(ht_cross), ref_ht_cross, rtol=.01)
            
              
        elif (varname == "interpline"):
            
            ref_t2_line = _get_refvals(referent, "interpline", repeat, multi)
            
            t2 = getvar(in_wrfnc, "T2", timeidx=timeidx)
            pivot_point = CoordPair(t2.shape[-1] // 2, t2.shape[-2] // 2)
            
            # Check that it works with numpy arrays
            t2_line1 = interpline(to_np(t2), pivot_point=pivot_point, 
                                  angle=90.0)
            #print (t2_line1)
            t2_line1 = interpline(t2, pivot_point=pivot_point, angle=90.0)
            
            nt.assert_allclose(to_np(t2_line1), ref_t2_line)
            
        elif (varname == "vinterp"):
            # Tk to theta
            fld_tk_theta = _get_refvals(referent, "vinterp", repeat, multi)
            fld_tk_theta = np.squeeze(fld_tk_theta)
            
            tk = getvar(in_wrfnc, "temp", timeidx=timeidx, units="k")
            
            interp_levels = [200,300,500,1000]
            
            # Check that it works with numpy arrays
            field = vinterp(in_wrfnc, 
                            field=to_np(tk), 
                            vert_coord="theta", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk",
                            timeidx=timeidx, 
                            log_p=True)
            #print (field)
            
            field = vinterp(in_wrfnc, 
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
            
            nt.assert_allclose(to_np(field), fld_tk_theta, tol, atol)
            
    return test


def make_latlon_test(testid, wrf_in, referent, single, multi=False, repeat=3, 
                     pynio=False):
    def test(self):
        from netCDF4 import Dataset as NetCDF
        
        timeidx = 0
        in_wrfnc = NetCDF(wrf_in)
 
        refnc = NetCDF(referent)
                
        if testid == "xy":
            # Since this domain is not moving, the reference values are the 
            # same whether there are multiple or single files
            ref_vals = refnc.variables["xy"][:]
            # Lats/Lons taken from NCL script, just hard-coding for now
            lats = [-55, -60, -65]
            lons = [25, 30, 35]
            
            xy = ll_to_xy(in_wrfnc, lats[0], lons[0])
                
            nt.assert_allclose(to_np(xy), ref_vals)
                
                
        else:
            # Since this domain is not moving, the reference values are the 
            # same whether there are multiple or single files
            ref_vals = refnc.variables["ll"][:]
            
             # i_s, j_s taken from NCL script, just hard-coding for now
             # NCL uses 1-based indexing for this, so need to subtract 1
            i_s = np.asarray([10, 100, 150], int) - 1
            j_s = np.asarray([10, 100, 150], int) - 1
            
            ll = xy_to_ll(in_wrfnc, i_s[0], j_s[0])
            
            nt.assert_allclose(to_np(ll), ref_vals)
                
        
    return test

class WRFVarsTest(ut.TestCase):
    longMessage = True
    
class WRFInterpTest(ut.TestCase):
    longMessage = True
    
class WRFLatLonTest(ut.TestCase):
    longMessage = True
        

if __name__ == "__main__":
    ignore_vars = []  # Not testable yet
    wrf_vars = ["avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
                "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
                "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
                "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
                "wa", "uvmet10", "uvmet", "z", "cfrac", "zstag"]
    interp_methods = ["interplevel", "vertcross", "interpline", "vinterp"]
    latlon_tests = ["xy", "ll"]
    
    import netCDF4
    for var in wrf_vars:
        if var in ignore_vars:
            continue
        
        test_func1 = make_test(var, TEST_FILE, REF_FILE)
        setattr(WRFVarsTest, 'test_{0}'.format(var), test_func1)
        
    for method in interp_methods:
        test_interp_func1 = make_interp_test(method, TEST_FILE, 
                                             REF_FILE)
        setattr(WRFInterpTest, 'test_{0}'.format(method), 
                test_interp_func1)
    
    for testid in latlon_tests:
        for single in (True,):
            for multi in (False,):
                test_ll_func = make_latlon_test(testid, TEST_FILE, 
                                                REF_FILE, 
                                                single=single, multi=multi, 
                                                repeat=3, pynio=False)
                multistr = "" if not multi else "_multi"
                singlestr = "_nosingle" if not single else "_single"
                test_name = "test_{}{}{}".format(testid, singlestr, 
                                                   multistr)
                setattr(WRFLatLonTest, test_name, test_ll_func)
                
     
    ut.main()
    