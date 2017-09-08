from __future__ import print_function, division

import os
import argparse

import numpy as np
from netCDF4 import Dataset

from wrf import (getvar, interplevel, interpline, vertcross, vinterp, py2round,
                 CoordPair, ll_to_xy, xy_to_ll, to_np)

VARS_TO_KEEP = ("XLAT", "XLONG", "XLAT_U", "XLAT_V", "XLONG_U", "XLONG_V",
                "U", "V", "W", "PH", "PHB", "T", "P", "PB", "Q2", "T2",
                "PSFC", "U10", "V10", "XTIME", "QVAPOR", "QCLOUD", 
                "QGRAUP", "QRAIN", "QSNOW", "MAPFAC_M", "MAPFAC_U",
                "MAPFAC_V", "F", "HGT", "RAINC", "RAINSH", "RAINNC")

DIMS_TO_TRIM = ("west_east", "south_north", "bottom_top", "bottom_top_stag",
                "west_east_stag", "south_north_stag")

WRF_DIAGS = ["avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
             "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
             "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
             "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
             "wa", "uvmet10", "uvmet", "z", "cfrac"]

INTERP_METHS = ["interplevel", "vertcross", "interpline", "vinterp"]

LATLON_METHS = ["xy", "ll"]


def copy_and_reduce(opts):
    infilename = opts.filename
    outfilename = os.path.expanduser(
        os.path.join(opts.outdir, "ci_test_file.nc"))
    with Dataset(infilename) as infile, Dataset(outfilename, "w") as outfile:
    
        # Copy the global attributes
        outfile.setncatts(infile.__dict__)
        
        # Copy Dimensions
        for name, dimension in infile.dimensions.iteritems():
            if name in DIMS_TO_TRIM:
                if name.find("_stag") > 0:
                    dimsize = (len(dimension) / 2 
                               if len(dimension) % 2 == 0 
                               else (len(dimension) + 1) / 2)
                else:
                    dimsize = (len(dimension) / 2 
                               if len(dimension) % 2 == 0 
                               else (len(dimension) - 1) / 2)
            else:
                dimsize = len(dimension)
                
            outfile.createDimension(name, dimsize)
        
        # Copy Variables  
        for name, variable in infile.variables.iteritems():
            if name not in VARS_TO_KEEP: 
                continue
        
            outvar = outfile.createVariable(name, variable.datatype, 
                                       variable.dimensions, zlib=True)
            
            in_slices = tuple(slice(0, dimsize) for dimsize in outvar.shape)
            
            outvar[:] = variable[in_slices]
            outvar.setncatts(infile.variables[name].__dict__)

       
def add_to_ncfile(outfile, var, varname):
    dim_d = dict(zip(var.dims, var.shape))
            
    for dimname, dimsize in dim_d.items():
        if dimname not in outfile.dimensions:
            outfile.createDimension(dimname, dimsize)
    
    fill_value = None
    if "_FillValue" in var.attrs:
        fill_value = var.attrs["_FillValue"]
        
    ncvar = outfile.createVariable(varname, var.dtype, var.dims,
                                   zlib=True, fill_value=fill_value)
    
    var_vals = to_np(var)
    ncvar[:] = var_vals[:]


def make_result_file(opts):
    infilename = os.path.expanduser(
        os.path.join(opts.outdir, "ci_test_file.nc"))
    outfilename = os.path.expanduser(
        os.path.join(opts.outdir, "ci_result_file.nc"))
    
    with Dataset(infilename) as infile, Dataset(outfilename, "w") as outfile:
    
        for varname in WRF_DIAGS:
            var = getvar(infile, varname)
            add_to_ncfile(outfile, var, varname)
            
        for interptype in INTERP_METHS:
            if interptype == "interpline":
                hts = getvar(infile, "z")
                p = getvar(infile, "pressure")
                hts_850 = interplevel(hts, p, 850.)
                
                add_to_ncfile(outfile, hts_850, "interplevel")
                
            if interptype == "vertcross":
                
                hts = getvar(infile, "z")
                p = getvar(infile, "pressure")
                
                pivot_point = CoordPair(hts.shape[-1] // 2, hts.shape[-2] // 2) 
                ht_cross = vertcross(hts, p, pivot_point=pivot_point, 
                                     angle=90.)
    
                add_to_ncfile(outfile, ht_cross, "vertcross")
                
            if interptype == "interpline":
            
                t2 = getvar(infile, "T2")
                pivot_point = CoordPair(t2.shape[-1] // 2, t2.shape[-2] // 2)
                
                t2_line = interpline(t2, pivot_point=pivot_point, angle=90.0)
                
                add_to_ncfile(outfile, t2_line, "interpline")
                
            if interptype == "vinterp":
                
                tk = getvar(infile, "temp", units="k")
                
                interp_levels = [200,300,500,1000]
                
                field = vinterp(infile, 
                                field=tk, 
                                vert_coord="theta", 
                                interp_levels=interp_levels, 
                                extrapolate=True, 
                                field_type="tk",
                                log_p=True)
                
                add_to_ncfile(outfile, field, "vinterp")
                
        for latlonmeth in LATLON_METHS:
            if latlonmeth == "xy":
                # Hardcoded values from other unit tests
                lats = [-55, -60, -65]
                lons = [25, 30, 35]
                
                xy = ll_to_xy(infile, lats[0], lons[0])
                add_to_ncfile(outfile, xy, "xy")
            else:
                # Hardcoded from other unit tests
                i_s = np.asarray([10, 100, 150], int) - 1
                j_s = np.asarray([10, 100, 150], int) - 1
                
                ll = xy_to_ll(infile, i_s[0], j_s[0])
                add_to_ncfile(outfile, ll, "ll")

def main(opts):
    copy_and_reduce(opts)
    make_result_file(opts)
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Generate conda test files "
                                     "for unit testing.")
    parser.add_argument("-f", "--filename", required=True, 
                        help="the WRF test file")
    parser.add_argument("-o", "--outdir", required=True,
                        help="the location for the output files")
    opts = parser.parse_args()
    
    main(opts)