import sys
import unittest as ut

import numpy.testing as nt 
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset as nc

from wrf import getvar

TEST_FILE = "/Users/ladwig/Documents/wrf_files/wrfout_d01_2010-06-13_21:00:00"

# Python 3
if sys.version_info > (3,):
    xrange = range
    

class TestUnits(ut.TestCase):
    longMessage = True
    
    def test_temp_units(self):
        wrfnc = nc(TEST_FILE)
        
        for units in ("K", "degC", "degF"):
            var = getvar(wrfnc, "temp", units=units)
            
            self.assertEqual(var.attrs["units"], units)
    
    def test_wind_units(self):
        wrfnc = nc(TEST_FILE)
        
        for units in ("m s-1", "kt", "mi h-1", "km h-1", "ft s-1"):
            var = getvar(wrfnc, "uvmet", units=units)
            
            self.assertEqual(var.attrs["units"], units)
            
    def test_pres_units(self):
        wrfnc = nc(TEST_FILE)
        
        for units in ("Pa", "hPa", "mb", "torr", "mmHg", "atm"):
            var = getvar(wrfnc, "slp", units=units)
            
            self.assertEqual(var.attrs["units"], units)
            
    def test_height_units(self):
        wrfnc = nc(TEST_FILE)
        
        for units in ("m", "km", "dam", "ft", "mi"):
            var = getvar(wrfnc, "z", units=units)
            
            self.assertEqual(var.attrs["units"], units)
        
            

if __name__ == "__main__":
    ut.main()