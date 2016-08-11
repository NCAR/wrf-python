from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants
from .destag import destagger
#from .extension import computeomega,computetk
from .extension import _omega, _tk
from .util import extract_vars
from .metadecorators import copy_and_set_metadata


@copy_and_set_metadata(copy_varname="T", name="omega", 
                       description="omega",
                       units="Pa/s")
def get_omega(wrfnc, timeidx=0, method="cat", squeeze=True, cache=None,
              meta=True):
    varnames=("T", "P", "W", "PB", "QVAPOR")
    ncvars = extract_vars(wrfnc, timeidx, varnames, method, squeeze, cache,
                          meta=False)
    t = ncvars["T"]
    p = ncvars["P"]
    w = ncvars["W"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    
    wa = destagger(w, -3)
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = _tk(full_p, full_t)
    
    omega = _omega(qv, tk, wa, full_p)
    
    return omega
    