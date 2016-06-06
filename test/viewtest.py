import numpy as np

import wrf._wrffortran

a = np.ones((3,3,3))
b = np.zeros((3,3,3,3))
errstat = np.array(0)
errmsg = np.zeros(512, "c")


for i in xrange(2):
    outview = b[i,:]
    outview = outview.T
    q = wrf._wrffortran.testfunc(a,outview,errstat=errstat,errstr=errmsg)
    q[1,1,1] = 100


print errstat
print b
str_bytes = (bytes(c).decode("utf-8") for c in errmsg[:])
print "".join(str_bytes).strip()