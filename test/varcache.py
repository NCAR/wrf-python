from __future__ import print_function

import time
from netCDF4 import Dataset
from wrf import getvar, ALL_TIMES, extract_vars

wrf_filenames = ["/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-28_00:00:00",
                 "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-28_12:00:00", 
                 "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-29_00:00:00"]

wrfin = [Dataset(x) for x in wrf_filenames]

my_cache = extract_vars(wrfin, ALL_TIMES, ("P", "PB", "PH", "PHB", "T", "QVAPOR", "HGT", "U", "V", "W", "PSFC"))

start = time.time()
for var in ("avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
            "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
            "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
            "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
            "wa", "uvmet10", "uvmet", "z", "cfrac", "zstag", "geopt_stag"):
    v = getvar(wrfin, var, ALL_TIMES)
end = time.time()

print ("Time taken without variable cache: ", (end-start))

start = time.time()
for var in ("avo", "eth",  "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
            "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
            "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
            "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
            "wa", "uvmet10", "uvmet", "z", "cfrac", "zstag", "geopt_stag"):
    v = getvar(wrfin, var, ALL_TIMES, cache=my_cache)
end = time.time()

print ("Time taken with variable cache: ", (end-start))

