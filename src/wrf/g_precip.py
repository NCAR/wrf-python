from __future__ import (absolute_import, division, print_function)

from .util import extract_vars


def get_accum_precip(wrfin, timeidx=0):
    ncvars = extract_vars(wrfin, timeidx, varnames=("RAINC", "RAINNC"))
    rainc = ncvars["RAINC"]
    rainnc = ncvars["RAINNC"]

    rainsum = rainc + rainnc

    return rainsum


def get_precip_diff(wrfin1, wrfin2, timeidx=0):
    vars1 = extract_vars(wrfin1, timeidx, varnames=("RAINC", "RAINNC"))
    vars2 = extract_vars(wrfin2, timeidx, varnames=("RAINC", "RAINNC"))
    rainc1 = vars1["RAINC"]
    rainnc1 = vars1["RAINNC"]

    rainc2 = vars2["RAINC"]
    rainnc2 = vars2["RAINNC"]

    rainsum1 = rainc1 + rainnc1
    rainsum2 = rainc2 + rainnc2

    return (rainsum1 - rainsum2)

# TODO:  Handle bucket flipping
