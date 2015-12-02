import unittest as ut
import numpy.testing as nt 
import numpy as n
import os, sys
import subprocess

import Ngl, Nio
from wrf.var import getvar

# This should work with Nio
from netCDF4 import Dataset as NetCDF
from boto.ec2.instancestatus import Status

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
def make_test(varname, wrf_in, referent):
    def test(self):
        in_wrfnc = NetCDF(wrf_in)
        refnc = NetCDF(referent)
        
        ref_vals = refnc.variables[varname][...]
        
        if (varname == "tc"):
            my_vals = getvar(in_wrfnc, "temp", units="c")
            tol = 0
            atol = .1
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "tk"):
            my_vals = getvar(in_wrfnc, "temp", units="k")
            tol = 0
            atol = .1
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "td"):
            my_vals = getvar(in_wrfnc, "td", units="c")
            tol = 0
            atol = .1
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "pressure"):
            my_vals = getvar(in_wrfnc, varname, units="hpa")
            tol = 1/100.
            atol = 0
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "p"):
            my_vals = getvar(in_wrfnc, varname, units="pa")
            tol = 1/100.
            atol = 0
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "slp"):
            my_vals = getvar(in_wrfnc, varname, units="hpa")
            tol = 2/100.
            atol = 0
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
            
        elif (varname == "uvmet"):
            my_vals = getvar(in_wrfnc, varname)
            tol = 1/100.
            atol = .5
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "uvmet10"):
            my_vals = getvar(in_wrfnc, varname)
            tol = 1/100.
            atol = .5
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
            
        elif (varname == "omg"):
            my_vals = getvar(in_wrfnc, varname)
            tol = 2/100.
            atol = 0
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "ctt"):
            my_vals = getvar(in_wrfnc, varname)
            tol = 1/100.
            atol = 0
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
        elif (varname == "cape_2d"):
            mcape, mcin, lcl, lfc = getvar(in_wrfnc, varname)
            tol = 1/100.
            atol = 0
            nt.assert_allclose(mcape, ref_vals[0,:,:], tol, atol)
            nt.assert_allclose(mcin, ref_vals[1,:,:], tol, atol)
            nt.assert_allclose(lcl, ref_vals[2,:,:], tol, atol)
            nt.assert_allclose(lfc, ref_vals[3,:,:], tol, atol)
        elif (varname == "cape_3d"):
            cape, cin = getvar(in_wrfnc, varname)
            tol = 1/100.
            atol = 0
            nt.assert_allclose(cape, ref_vals[0,:,:], tol, atol)
            nt.assert_allclose(cin, ref_vals[1,:,:], tol, atol)
            
            
        else:
            my_vals = getvar(in_wrfnc, varname)
            tol = 1/100.
            atol = 0
            nt.assert_allclose(my_vals, ref_vals, tol, atol)
    
    
    return test

class WRFVarsTest(ut.TestCase):
    longMessage = True
        

if __name__ == "__main__":
    ignore_vars = []  # Not testable yet
    wrf_vars = ["avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
                "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
                "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
                "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
                "wa", "uvmet10", "uvmet", "z", "ctt", "cape_2d", "cape_3d"]
    
    for var in wrf_vars:
        if var in ignore_vars:
            continue
        
        test_func = make_test(var, TEST_FILE, OUT_NC_FILE)
        setattr(WRFVarsTest, 'test_{0}'.format(var), test_func)
    
    
    ut.main()
    