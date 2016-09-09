import time
from wrf.cache import _get_cache
from wrf import getvar
from netCDF4 import Dataset as nc

#a = nc("/Users/ladwig/Documents/wrf_files/wrf_vortex_single/wrfout_d02_2005-08-28_00:00:00")
#b = nc("/Users/ladwig/Documents/wrf_files/wrf_vortex_single/wrfout_d02_2005-08-28_03:00:00")
a = nc("/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-28_00:00:00")
b = nc("/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-28_12:00:00")
q = {"outoutoutout" : {"outoutout" : {"outout" : {"out1" : {"blah" : [a,b], "blah2" : [a,b]}, "out2" : {"blah" : [a,b], "blah2" : [a,b]} } } } }

t1 = time.time()
c = getvar(q, "rh", method="cat", timeidx=None, squeeze=True)
t2 = time.time()
print (c)
print ("time taken: {}".format((t2-t1)*1000.))

t1 = time.time()
c = getvar(q, "rh", method="cat", timeidx=None, squeeze=False)
t2 = time.time()
print (c)
print ("time taken: {}".format((t2-t1)*1000.))

t1 = time.time()
c = getvar(q, "rh", method="cat", timeidx=1, squeeze=True)
t2 = time.time()
print (c)
print ("time taken: {}".format((t2-t1)*1000.))

t1 = time.time()
c = getvar(q, "rh", method="cat", timeidx=1, squeeze=False)
t2 = time.time()
print(c)
print ("time taken: {}".format((t2-t1)*1000.))

t1 = time.time()
c = getvar(q, "rh", method="join", timeidx=None, squeeze=True)
t2 = time.time()
print (c)
print ("time taken: {}".format((t2-t1)*1000.))

t1 = time.time()
c = getvar(q, "rh", method="join", timeidx=None, squeeze=False)
t2 = time.time()
print(c)
print ("time taken: {}".format((t2-t1)*1000.))

t1 = time.time()
c = getvar(q, "rh", method="join", timeidx=1, squeeze=True)
t2 = time.time()
print (c)
print ("time taken: {}".format((t2-t1)*1000.))

t1 = time.time()
c = getvar(q, "rh", method="join", timeidx=1, squeeze=False)
t2 = time.time()
print (c)
print ("time taken: {}".format((t2-t1)*1000.))
