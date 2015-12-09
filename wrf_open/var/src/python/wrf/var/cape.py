import numpy as n
import numpy.ma as ma

from wrf.var.extension import computetk,computecape
from wrf.var.destagger import destagger
from wrf.var.constants import Constants, ConversionFactors
from wrf.var.util import extract_vars, hold_dim_fixed

__all__ = ["get_2dcape", "get_3dcape"]

def get_2dcape(wrfnc, missing=-999999.0, timeidx=0):
    """Return the 2d fields of cape, cin, lcl, and lfc"""
    ncvars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "QVAPOR", "PH",
                                              "PHB", "HGT", "PSFC"))
    
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    ter = ncvars["HGT"]
    psfc = ncvars["PSFC"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, 0)
    z = geopt_unstag/Constants.G
    
    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc 
    
    i3dflag = 0
    ter_follow = 1
    
    cape_res,cin_res = computecape(p_hpa,tk,qv,z,ter,psfc_hpa,
                                   missing,i3dflag,ter_follow)
     
    cape = hold_dim_fixed(cape_res, -3, 0)
    cin = hold_dim_fixed(cin_res, -3, 0)
    lcl = hold_dim_fixed(cin_res, -3, 1)
    lfc = hold_dim_fixed(cin_res, -3, 2)
    
    resdim = [4]
    resdim += cape.shape
    
    # Make a new output array for the result
    res = n.zeros((resdim), cape.dtype)
    
    res[0,:] = cape[:]
    res[1,:] = cin[:]
    res[2,:] = lcl[:]
    res[3,:] = lfc[:]
    
    return ma.masked_values(res,missing)

def get_3dcape(wrfnc, missing=-999999.0, timeidx=0):
    """Return the 3d fields of cape and cin"""
    
    ncvars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "QVAPOR", 
                                                "PH", "PHB", "HGT", "PSFC"))
    t = ncvars["T"]
    p = ncvars["P"]
    pb = ncvars["PB"]
    qv = ncvars["QVAPOR"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    ter = ncvars["HGT"]
    psfc = ncvars["PSFC"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, 0)
    z = geopt_unstag/Constants.G
    
    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc 
    
    i3dflag = 1
    ter_follow = 1
    
    cape,cin = computecape(p_hpa,tk,qv,z,ter,psfc_hpa,
                           missing,i3dflag,ter_follow)
    
    # Make a new output array for the result
    resdim = [2]
    resdim += cape.shape
    res = n.zeros((resdim), cape.dtype)
    
    res[0,:] = cape[:]
    res[1,:] = cin[:]
    
    return ma.masked_values(res, missing)
    
    
    