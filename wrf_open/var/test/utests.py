import unittest as ut
import numpy.testing as nt 
import numpy as n
import numpy.ma as ma
import os, sys
import subprocess

from wrf.var import getvar, interplevel, interpline, vertcross, vinterp

NCL_EXE = "/Users/ladwig/nclbuild/6.3.0/bin/ncl"
TEST_FILE = "/Users/ladwig/Documents/wrf_files/wrfout_d01_2010-06-13_21:00:00"
OUT_NC_FILE = "/tmp/wrftest.nc"

def setUpModule():
    ncarg_root = os.environ.get("NCARG_ROOT", None)
    if ncarg_root is None:
        raise RuntimeError("$NCARG_ROOT environment variable not set")
        
    
    this_path = os.path.realpath(__file__)
    ncl_script = os.path.join(os.path.dirname(this_path),
                              "ncl_get_var.ncl")
    ncfile = TEST_FILE + ".nc" # NCL requires extension
    
    # This needs to be set when PyNIO is installed, since PyNIOs data does
    # not contain the dat file for the CAPE calcluations
    os.environ["NCARG_NCARG"] = os.path.join(os.environ["NCARG_ROOT"],
                                             "lib", "ncarg")
    cmd = "%s %s 'in_file=\"%s\"' 'out_file=\"%s\"'" % (NCL_EXE,
                                                        ncl_script,
                                                        ncfile,
                                                        OUT_NC_FILE)
    
    #print cmd
    
    if not os.path.exists(OUT_NC_FILE):
        status = subprocess.call(cmd, shell=True)
        if (status != 0):
            raise RuntimeError("NCL script failed.  Could not set up test.")

# Using helpful information at: 
# http://eli.thegreenplace.net/2014/04/02/dynamically-generating-python-test-cases
def make_test(varname, wrf_in, referent, multi=False, repeat=3, pynio=False):
    def test(self):
        try:
            from netCDF4 import Dataset as NetCDF
        except:
            pass
        
        try:
            from PyNIO import Nio
        except:
            pass
        
        if not multi:
            if not pynio:
                in_wrfnc = NetCDF(wrf_in)
            else:
                # Note: Python doesn't like it if you reassign an outer scoped
                # variable (wrf_in in this case)
                if not wrf_in.endswith(".nc"):
                    wrf_file = wrf_in + ".nc"
                else:
                    wrf_file = wrf_in
                in_wrfnc = Nio.open_file(wrf_file)
        else:
            if not pynio:
                nc = NetCDF(wrf_in)
                in_wrfnc = [nc for i in xrange(repeat)]
            else:
                if not wrf_in.endswith(".nc"):
                    wrf_file = wrf_in + ".nc"
                else:
                    wrf_file = wrf_in
                nc = Nio.open_file(wrf_file)
                in_wrfnc = [nc for i in xrange(repeat)]
        
        refnc = NetCDF(referent)
        
        if not multi:
            ref_vals = refnc.variables[varname][:]
        else:
            data = refnc.variables[varname][:]
            new_dims = [repeat] + [x for x in data.shape]
            masked=False
            if (isinstance(data, ma.core.MaskedArray)):
                masked=True
            
            if not masked:
                ref_vals = n.zeros(new_dims, data.dtype)
            else:
                ref_vals = ma.asarray(n.zeros(new_dims, data.dtype))
                
            for i in xrange(repeat):
                ref_vals[i,:] = data[:]
                
                if masked:
                    ref_vals.mask[i,:] = data.mask[:]
        
        if (varname == "tc"):
            my_vals = getvar(in_wrfnc, "temp", units="c")
            tol = 0
            atol = .1 # Note:  NCL uses 273.16 as conversion for some reason
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "pw"):
            my_vals = getvar(in_wrfnc, "pw")
            tol = .5/100.0
            atol = 0 # NCL uses different constants and doesn't use same
                     # handrolled virtual temp in method
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "cape_2d"):
            cape_2d = getvar(in_wrfnc, varname)
            tol = 0/100.    # Not sure why different, F77 vs F90?
            atol = .1
            nt.assert_allclose(cape_2d, ref_vals, tol, atol)
        elif (varname == "cape_3d"):
            cape_3d = getvar(in_wrfnc, varname)
            tol = 0/100.  # Not sure why different, F77 vs F90?
            atol = .01
            nt.assert_allclose(cape_3d, ref_vals, tol, atol)
        else:
            my_vals = getvar(in_wrfnc, varname)
            tol = 0/100.
            atol = 0.0001
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
    
    
    return test

def _get_refvals(referent, varname, repeat, multi):
    try:
        from netCDF4 import Dataset as NetCDF
    except:
        pass
        
    refnc = NetCDF(referent)
    
    if not multi:
        ref_vals = refnc.variables[varname][:]
    else:
        data = refnc.variables[varname][:]
        new_dims = [repeat] + [x for x in data.shape]
        masked=False
        if (isinstance(data, ma.core.MaskedArray)):
            masked=True
          
        if not masked:
            ref_vals = n.zeros(new_dims, data.dtype)
        else:
            ref_vals = ma.asarray(n.zeros(new_dims, data.dtype))
              
        for i in xrange(repeat):
            ref_vals[i,:] = data[:]
              
            if masked:
                ref_vals.mask[i,:] = data.mask[:]
                
    return ref_vals

def make_interp_test(varname, wrf_in, referent, multi=False, 
                     repeat=3, pynio=False):
    def test(self):
        try:
            from netCDF4 import Dataset as NetCDF
        except:
            pass
        
        try:
            from PyNIO import Nio
        except:
            pass
        
        if not multi:
            if not pynio:
                in_wrfnc = NetCDF(wrf_in)
            else:
                # Note: Python doesn't like it if you reassign an outer scoped
                # variable (wrf_in in this case)
                if not wrf_in.endswith(".nc"):
                    wrf_file = wrf_in + ".nc"
                else:
                    wrf_file = wrf_in
                in_wrfnc = Nio.open_file(wrf_file)
        else:
            if not pynio:
                nc = NetCDF(wrf_in)
                in_wrfnc = [nc for i in xrange(repeat)]
            else:
                if not wrf_in.endswith(".nc"):
                    wrf_file = wrf_in + ".nc"
                else:
                    wrf_file = wrf_in
                nc = Nio.open_file(wrf_file)
                in_wrfnc = [nc for i in xrange(repeat)]
        
        if (varname == "interplevel"):
            ref_ht_500 = _get_refvals(referent, "z_500", repeat, multi)
            hts = getvar(in_wrfnc, "z")
            p = getvar(in_wrfnc, "pressure")
            hts_500 = interplevel(hts, p, 500)
            
            nt.assert_allclose(hts_500, ref_ht_500)
            
        elif (varname == "vertcross"):
            ref_ht_cross = _get_refvals(referent, "ht_cross", repeat, multi)
            ref_p_cross = _get_refvals(referent, "p_cross", repeat, multi)
            
            hts = getvar(in_wrfnc, "z")
            p = getvar(in_wrfnc, "pressure")
            
            pivot_point = (hts.shape[-2] / 2, hts.shape[-1] / 2) 
            ht_cross = vertcross(hts, p, pivot_point=pivot_point,angle=90.)
            
            nt.assert_allclose(ht_cross, ref_ht_cross, rtol=.01)
            
            # Test opposite
            p_cross1 = vertcross(p,hts,pivot_point=pivot_point, angle=90.0)

            nt.assert_allclose(p_cross1, 
                               ref_p_cross, 
                               rtol=.01)
            # Test point to point
            start_point = (hts.shape[-2]/2, 0)
            end_point = (hts.shape[-2]/2, -1)
            
            p_cross2 = vertcross(p,hts,start_point=start_point, 
                                end_point=end_point)
            
            nt.assert_allclose(p_cross1, p_cross2)
            
        elif (varname == "interpline"):
            ref_t2_line = _get_refvals(referent, "t2_line", repeat, multi)
            
            t2 = getvar(in_wrfnc, "T2")
            pivot_point = (t2.shape[-2] / 2, t2.shape[-1] / 2) 
            
            t2_line1 = interpline(t2, pivot_point=pivot_point, angle=90.0)
            
            nt.assert_allclose(t2_line1, ref_t2_line)
            
            # Test point to point
            start_point = (t2.shape[-2]/2, 0)
            end_point = (t2.shape[-2]/2, -1)
            
            t2_line2 = interpline(t2, start_point=start_point, 
                                  end_point=end_point)
            
            nt.assert_allclose(t2_line1, t2_line2)
        elif (varname == "vinterp"):
            # Tk to theta
            fld_tk_theta = _get_refvals(referent, "fld_tk_theta", repeat, multi)
            fld_tk_theta = n.squeeze(fld_tk_theta)
            
            tk = getvar(in_wrfnc, "temp", units="k")
            
            interp_levels = [200,300,500,1000]
            
            field = vinterp(in_wrfnc, 
                            field=tk, 
                            vert_coord="theta", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk", 
                            log_p=True)
            
            tol = 0/100.
            atol = 0.0001
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_tk_theta, tol, atol)
            
            # Tk to theta-e
            fld_tk_theta_e = _get_refvals(referent, "fld_tk_theta_e", repeat, multi)
            fld_tk_theta_e = n.squeeze(fld_tk_theta_e)
            
            interp_levels = [200,300,500,1000]
            
            field = vinterp(in_wrfnc, 
                            field=tk, 
                            vert_coord="theta-e", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk", 
                            log_p=True)
            
            tol = 0/100.
            atol = 0.0001
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_tk_theta_e, tol, atol)
            
            # Tk to pressure
            fld_tk_pres = _get_refvals(referent, "fld_tk_pres", repeat, multi)
            fld_tk_pres = n.squeeze(fld_tk_pres)
            
            interp_levels = [850,500]
            
            field = vinterp(in_wrfnc, 
                            field=tk, 
                            vert_coord="pressure", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk", 
                            log_p=True)
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_tk_pres, tol, atol)
            
            # Tk to geoht_msl
            fld_tk_ght_msl = _get_refvals(referent, "fld_tk_ght_msl", repeat, multi)
            fld_tk_ght_msl = n.squeeze(fld_tk_ght_msl)
            interp_levels = [1,2]
            
            field = vinterp(in_wrfnc, 
                            field=tk, 
                            vert_coord="ght_msl", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk", 
                            log_p=True)
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_tk_ght_msl, tol, atol)
            
            # Tk to geoht_agl
            fld_tk_ght_agl = _get_refvals(referent, "fld_tk_ght_agl", repeat, multi)
            fld_tk_ght_agl = n.squeeze(fld_tk_ght_agl)
            interp_levels = [1,2]
            
            field = vinterp(in_wrfnc, 
                            field=tk, 
                            vert_coord="ght_agl", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="tk", 
                            log_p=True)
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_tk_ght_agl, tol, atol)
            
            # Hgt to pressure
            fld_ht_pres = _get_refvals(referent, "fld_ht_pres", repeat, multi)
            fld_ht_pres = n.squeeze(fld_ht_pres)
            
            z = getvar(in_wrfnc, "height", units="m")
            interp_levels = [500,50]
            field = vinterp(in_wrfnc, 
                            field=z, 
                            vert_coord="pressure", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="ght", 
                            log_p=True)
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_ht_pres, tol, atol)
            
            # Pressure to theta
            fld_pres_theta = _get_refvals(referent, "fld_pres_theta", repeat, multi)
            fld_pres_theta = n.squeeze(fld_pres_theta)
            
            p = getvar(in_wrfnc, "pressure")
            interp_levels = [200,300,500,1000]
            field = vinterp(in_wrfnc, 
                            field=p, 
                            vert_coord="theta", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="pressure", 
                            log_p=True)
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_pres_theta, tol, atol)
            
            # Theta-e to pres
            fld_thetae_pres = _get_refvals(referent, "fld_thetae_pres", repeat, multi)
            fld_thetae_pres = n.squeeze(fld_thetae_pres)
            
            eth = getvar(in_wrfnc, "eth")
            interp_levels = [850,500,5]
            field = vinterp(in_wrfnc, 
                            field=eth, 
                            vert_coord="pressure", 
                            interp_levels=interp_levels, 
                            extrapolate=True, 
                            field_type="theta-e", 
                            log_p=True)
            
            field = n.squeeze(field)
            nt.assert_allclose(field, fld_thetae_pres, tol, atol)
    
    return test

class WRFVarsTest(ut.TestCase):
    longMessage = True
    
class WRFInterpTest(ut.TestCase):
    longMessage = True
        

if __name__ == "__main__":
    ignore_vars = []  # Not testable yet
    wrf_vars = ["avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
                "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
                "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
                "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
                "wa", "uvmet10", "uvmet", "z", "ctt", "cape_2d", "cape_3d"]
    interp_methods = ["interplevel", "vertcross", "interpline", "vinterp"]
    
    try:
        import netCDF4
    except ImportError:
        pass
    else:
        for var in wrf_vars:
            if var in ignore_vars:
                continue
            
            test_func1 = make_test(var, TEST_FILE, OUT_NC_FILE)
            test_func2 = make_test(var, TEST_FILE, OUT_NC_FILE, multi=True)
            setattr(WRFVarsTest, 'test_{0}'.format(var), test_func1)
            setattr(WRFVarsTest, 'test_multi_{0}'.format(var), test_func2)
            
        for method in interp_methods:
            test_interp_func1 = make_interp_test(method, TEST_FILE, 
                                                 OUT_NC_FILE)
            test_interp_func2 = make_interp_test(method, TEST_FILE, 
                                                 OUT_NC_FILE, multi=True)
            setattr(WRFInterpTest, 'test_{0}'.format(method), 
                    test_interp_func1)
            setattr(WRFInterpTest, 'test_multi_{0}'.format(method), 
                    test_interp_func2)
        
    try:
        import PyNIO
    except ImportError:
        pass
    else:
        for var in wrf_vars:
            if var in ignore_vars:
                continue
            
            test_func1 = make_test(var, TEST_FILE, OUT_NC_FILE,pynio=True)
            test_func2 = make_test(var, TEST_FILE, OUT_NC_FILE, multi=True,
                                   pynio=True)
            setattr(WRFVarsTest, 'test_pynio_{0}'.format(var), test_func1)
            setattr(WRFVarsTest, 'test_pynio_multi_{0}'.format(var), 
                    test_func2)
            
        for method in interp_methods:
            test_interp_func1 = make_interp_test(method, TEST_FILE, 
                                                 OUT_NC_FILE)
            test_interp_func2 = make_interp_test(method, TEST_FILE, 
                                                 OUT_NC_FILE, multi=True)
            setattr(WRFInterpTest, 'test_pynio_{0}'.format(method), 
                    test_interp_func1)
            setattr(WRFInterpTest, 'test_pynio_multi_{0}'.format(method), 
                    test_interp_func2)
        
    ut.main()
    