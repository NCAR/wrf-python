import numpy as n

from wrf.var.extension import computeeta
from wrf.var.constants import Constants
from wrf.var.decorators import convert_units
from wrf.var.util import extract_vars

#__all__ = ["convert_eta"]
__all__ = []
# A useful utility, but should probably just use geopotential height when 
# plotting for AGL levels

# Eta definition (nu):
# nu = (P - Ptop) / (Psfc - Ptop)

# def convert_eta(wrfnc, p_or_z="ht", timeidx=0):
#     if (p_or_z.lower() == "height" or p_or_z.lower() == "ht" 
#             or p_or_z.lower() == "h"):
#         return_z = True
#     elif (p_or_z.lower() == "p" or p_or_z.lower() == "pres" 
#           or p_or_z.lower() == "pressure"):
#         return_z = False
#     
#     R = Constants.R
#     G = Constants.G
#     CP = Constants.CP
#     
#     # Keeping the slice notation to show the dimensions
#     # Note: Not sure if T00 should be used (290) or the usual hard-coded 300 for base
#     # theta
#     height_data = wrfnc.variables["HGT"][timeidx,:,:]
#     znu_data = wrfnc.variables["ZNU"][timeidx,:]
#     #t00_data = wrfnc.variables["T00"][timeidx]
#     psfc_data = wrfnc.variables["PSFC"][timeidx,:,:]
#     ptop_data = wrfnc.variables["P_TOP"][timeidx]
#     pth_data = wrfnc.variables["T"][timeidx,:,:,:] # Pert potential temp
#     
#     pcalc_data = n.zeros(pth_data.shape, dtype=n.float32)
#     mean_t_data = n.zeros(pth_data.shape, dtype=n.float32)
#     temp_data = n.zeros(pth_data.shape, dtype=n.float32)
#     z_data = n.zeros(pth_data.shape, dtype=n.float32)
#     
#     #theta_data = pth_data + t00_data
#     theta_data = pth_data + Constants.T_BASE
#     
#     for k in xrange(znu_data.shape[0]):
#         pcalc_data[k,:,:] = znu_data[k]*(psfc_data[:,:] - (ptop_data)) + (ptop_data)
#          
#     # Potential temperature:
#     # theta = T * (Po/P)^(R/CP)
#     #
#     # Hypsometric equation:
#     # h = z2-z1 = R*Tbar/G * ln(p1/p2)
#     # where z1, p1 are the surface
#     if return_z:
#         for k in xrange(znu_data.shape[0]):
#             temp_data[k,:,:] = (theta_data[k,:,:]) / ((100000.0 / (pcalc_data[k,:,:]))**(R/CP))      
#             mean_t_data[k,:,:] = n.mean(temp_data[0:k+1,:,:], axis=0)
#             z_data[k,:,:] = ((R*mean_t_data[k,:,:]/G) * n.log(psfc_data[:,:]/pcalc_data[k,:,:]))
#         
#         return z_data
#     else:
#         return pcalc_data * .01

# def convert_eta(wrfnc, units="m", msl=False, timeidx=0):
#     check_units(units, "height")
#     hgt = wrfnc.variables["HGT"][timeidx,:,:]
#     znu = wrfnc.variables["ZNU"][timeidx,:]
#     psfc = wrfnc.variables["PSFC"][timeidx,:,:]
#     ptop = wrfnc.variables["P_TOP"][timeidx]
#     t = wrfnc.variables["T"][timeidx,:,:,:]
#     
#     full_theta = t + Constants.T_BASE
#     
#     z = computeeta(full_theta, znu, psfc, ptop)
#     
#     if not msl:
#         return convert_units(z, "height", "m", units)
#     else:
#         return convert_units(z + hgt, "height", "m", units)


