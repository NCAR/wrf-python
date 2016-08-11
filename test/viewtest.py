import numpy as np

import wrf._wrffortran
errlen = int(wrf._wrffortran.constants.errlen)


a = np.ones((3,3,3))
b = np.zeros((3,3,3,3))
c = np.zeros(errlen, "c")
errstat = np.array(0)
errmsg = np.zeros(errlen, "c")

c[:] = "Test String"


for i in xrange(2):
    outview = b[i,:]
    outview = outview.T
    q = wrf._wrffortran.testfunc(a,outview,c,errstat=errstat,errmsg=errmsg)
    print errstat



str_bytes = (bytes(char).decode("utf-8") for char in errmsg[:])
print repr(errmsg)
print "".join(str_bytes).strip()