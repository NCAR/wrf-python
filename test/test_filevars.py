import unittest as ut
import numpy.testing as nt 
import numpy as np
import numpy.ma as ma
import os, sys
import subprocess

from wrf import getvar, ALL_TIMES

TEST_DIR = "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi"
TEST_FILENAMES = ["wrfout_d02_2005-08-28_00:00:00",
              "wrfout_d02_2005-08-28_12:00:00",
              "wrfout_d02_2005-08-29_00:00:00"]
TEST_FILES = [os.path.join(TEST_DIR, x) for x in TEST_FILENAMES]

# Python 3
if sys.version_info > (3,):
    xrange = range
    

class WRFFileVarsTest(ut.TestCase):
    longMessage = True

def make_test(ncfiles, varname):
    def test(self):
        t1 = getvar(ncfiles, varname, 0)
        t2 = getvar(ncfiles, varname, 0, meta=False)
        t3 = getvar(ncfiles, varname, ALL_TIMES)
        t4 = getvar(ncfiles, varname, ALL_TIMES, meta=False)
        t5 = getvar(ncfiles, varname, ALL_TIMES, method="join")
        
    return test
        

if __name__ == "__main__":
    from netCDF4 import Dataset
    
    ncfiles = [Dataset(x) for x in TEST_FILES]
    
    file_vars = ncfiles[0].variables.keys()
    
    ignore_vars = []
    
    for var in file_vars:
        if var in ignore_vars:
            continue
        
        test_func1 = make_test(ncfiles, var)
        setattr(WRFFileVarsTest, 'test_{0}'.format(var), test_func1)
     
    ut.main()