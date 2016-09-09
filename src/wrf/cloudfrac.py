from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .constants import Constants
from .extension import _tk, _rh, _cloudfrac
from .metadecorators import set_cloudfrac_metadata
from .util import extract_vars

@set_cloudfrac_metadata()
def get_cloudfrac(wrfnc, timeidx=0, method="cat", squeeze=True, 
                 cache=None, meta=True, _key=None):
    
    vars = extract_vars(wrfnc, timeidx, ("P", "PB", "QVAPOR", "T"), 
                          method, squeeze, cache, meta=False,
                          _key=_key)
    
    p = vars["P"]
    pb = vars["PB"]
    qv = vars["QVAPOR"]
    t = vars["T"]
    
    full_p = p + pb
    full_t = t + Constants.T_BASE
    
    tk = _tk(full_p, full_t)
    rh = _rh(qv, full_p, tk)
    
    return _cloudfrac(full_p, rh)
