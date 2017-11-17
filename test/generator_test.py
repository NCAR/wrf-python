from __future__ import (absolute_import, division, print_function, unicode_literals)

from wrf import getvar
from netCDF4 import Dataset as nc
#ncfile = nc("/Users/ladwig/Documents/wrf_files/wrfout_d01_2016-02-25_18_00_00")
ncfile = nc("/Users/ladwig/Documents/wrf_files/wrfout_d01_2016-10-07_00_00_00")

def gen_seq():
    wrfseq = [ncfile, ncfile, ncfile]
    for wrf in wrfseq:
        yield wrf
        
p_gen = getvar(gen_seq(), "P", method="join")

print(p_gen)
del p_gen

