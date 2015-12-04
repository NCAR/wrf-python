import numpy.ma as ma

from wrf.var.extension import computetk,computecape
from wrf.var.destagger import destagger
from wrf.var.constants import Constants, ConversionFactors
from wrf.var.util import extract_vars

__all__ = ["get_2dcape", "get_3dcape"]

def get_2dcape(wrfnc, missing=-999999.0, timeidx=0):
    """Return the 2d fields of cape, cin, lcl, and lfc"""
    vars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "QVAPOR", "PH",
                                              "PHB", "HGT", "PSFC"))
    
    t = vars["T"]
    p = vars["P"]
    pb = vars["PB"]
    qv = vars["QVAPOR"]
    ph = vars["PH"]
    phb = vars["PHB"]
    ter = vars["HGT"]
    psfc = vars["PSFC"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, 0)
    z = geopt_unstag/Constants.G
    
    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc # This may be the bug in NCL, as they pass this in 
                          # has Pa, but other pressure is hPa.  Converting to 
                          # hPa here.
    
    i3dflag = 0
    ter_follow = 1
    
    cape_res,cin_res = computecape(p_hpa,tk,qv,z,ter,psfc_hpa,
                                   missing,i3dflag,ter_follow)
    
    cape = cape_res[0,:,:]
    cin = cin_res[0,:,:]
    lcl = cin_res[1,:,:]
    lfc = cin_res[2,:,:]
    
    return (ma.masked_values(cape,missing), 
            ma.masked_values(cin,missing), 
            ma.masked_values(lcl,missing), 
            ma.masked_values(lfc,missing))

def get_3dcape(wrfnc, missing=-999999.0, timeidx=0):
    """Return the 3d fields of cape and cin"""
    
    vars = extract_vars(wrfnc, timeidx, vars=("T", "P", "PB", "QVAPOR", "PH",
                                              "PHB", "HGT", "PSFC"))
    t = vars["T"]
    p = vars["P"]
    pb = vars["PB"]
    qv = vars["QVAPOR"]
    ph = vars["PH"]
    phb = vars["PHB"]
    ter = vars["HGT"]
    psfc = vars["PSFC"]
    
    full_t = t + Constants.T_BASE
    full_p = p + pb
    tk = computetk(full_p, full_t)
    
    geopt = ph + phb
    geopt_unstag = destagger(geopt, 0)
    z = geopt_unstag/Constants.G
    
    # Convert pressure to hPa
    p_hpa = ConversionFactors.PA_TO_HPA * full_p
    psfc_hpa = ConversionFactors.PA_TO_HPA * psfc # This may be the bug in NCL, as they pass this in 
                          # has Pa, but other pressure is hPa.  Converting to 
                          # hPa here.
    
    i3dflag = 1
    ter_follow = 1
    
    cape,cin = computecape(p_hpa,tk,qv,z,ter,psfc_hpa,
                           missing,i3dflag,ter_follow)
    return (ma.masked_values(cape, missing), 
            ma.masked_values(cin, missing))
    
    
    