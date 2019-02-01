import unittest as ut
import numpy.testing as nt
import numpy as np
import numpy.ma as ma
import os
import sys
import subprocess

from wrf import (getvar, interplevel, interpline, vertcross, vinterp,
                 disable_xarray, xarray_enabled, to_np,
                 xy_to_ll, ll_to_xy, xy_to_ll_proj, ll_to_xy_proj,
                 extract_global_attrs, viewitems, CoordPair, ll_points)
from wrf.util import is_multi_file

IN_DIR = "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/moving_nest"
TEST_FILES = [os.path.join(IN_DIR, "wrfout_d02_2005-08-28_00:00:00"),
              os.path.join(IN_DIR, "wrfout_d02_2005-08-28_12:00:00"),
              os.path.join(IN_DIR, "wrfout_d02_2005-08-29_00:00:00")]


def wrfin_gen(wrf_in):
    for x in wrf_in:
        yield x


class wrf_in_iter_class(object):
    def __init__(self, wrf_in):
        self._wrf_in = wrf_in
        self._total = len(wrf_in)
        self._i = 0

    def __iter__(self):
        return self

    def next(self):
        if self._i >= self._total:
            raise StopIteration
        else:
            result = self._wrf_in[self._i]
            self._i += 1
            return result

    # Python 3
    def __next__(self):
        return self.next()


def make_test(varname, wrf_in):
    def test(self):
        x = getvar(wrf_in, varname, 0)
        xa = getvar(wrf_in, varname, None)

    return test


def make_interp_test(varname, wrf_in):
    def test(self):

        # Only testing vinterp since other interpolation just use variables
        if (varname == "vinterp"):
            for timeidx in (0, None):
                eth = getvar(wrf_in, "eth", timeidx=timeidx)
                interp_levels = [850, 500, 5]
                field = vinterp(wrf_in,
                                field=eth,
                                vert_coord="pressure",
                                interp_levels=interp_levels,
                                extrapolate=True,
                                field_type="theta-e",
                                timeidx=timeidx,
                                log_p=True)
        else:
            pass

    return test


def make_latlon_test(testid, wrf_in):
    def test(self):
        if testid == "xy_out":
            # Lats/Lons taken from NCL script, just hard-coding for now
            lats = [-55, -60, -65]
            lons = [25, 30, 35]
            xy = ll_to_xy(wrf_in, lats, lons, timeidx=0)
            xy = ll_to_xy(wrf_in, lats, lons, timeidx=None)
        else:
            i_s = np.asarray([10, 100, 150], int)
            j_s = np.asarray([10, 100, 150], int)
            ll = xy_to_ll(wrf_in, i_s, j_s, timeidx=0)
            ll = xy_to_ll(wrf_in, i_s, j_s, timeidx=None)

    return test


class WRFVarsTest(ut.TestCase):
    longMessage = True


class WRFInterpTest(ut.TestCase):
    longMessage = True


class WRFLatLonTest(ut.TestCase):
    longMessage = True


if __name__ == "__main__":
    from wrf import (omp_set_num_threads, omp_set_schedule, omp_get_schedule,
                     omp_set_dynamic, omp_get_num_procs, OMP_SCHED_STATIC)
    omp_set_num_threads(omp_get_num_procs()-1)
    omp_set_schedule(OMP_SCHED_STATIC, 0)
    omp_set_dynamic(False)

    ignore_vars = []  # Not testable yet
    wrf_vars = ["avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz",
                "geopt", "helicity", "lat", "lon", "omg", "p", "pressure",
                "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc",
                "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va",
                "wa", "uvmet10", "uvmet", "z", "cfrac", "zstag", "geopt_stag"]
    interp_methods = ["interplevel", "vertcross", "interpline", "vinterp"]
    latlon_tests = ["xy_out", "ll_out"]

    for nc_lib in ("netcdf4", "pynio", "scipy"):
        if nc_lib == "netcdf4":
            try:
                from netCDF4 import Dataset
            except ImportError:
                print("netcdf4-python not installed")
                continue
            else:
                test_in = [Dataset(x) for x in TEST_FILES]
        elif nc_lib == "pynio":
            try:
                from Nio import open_file
            except ImportError:
                print("PyNIO not installed")
                continue
            else:
                test_in = [open_file(x + ".nc", "r") for x in TEST_FILES]
        elif nc_lib == "scipy":
            try:
                from scipy.io.netcdf import netcdf_file
            except ImportError:
                print("scipy.io.netcdf not installed")
            else:
                test_in = [netcdf_file(x, mmap=False) for x in TEST_FILES]

        input0 = test_in[0]
        input1 = test_in
        input2 = tuple(input1)
        input3 = wrfin_gen(test_in)
        input4 = wrf_in_iter_class(test_in)
        input5 = {"A": input1,
                  "B": input2}
        input6 = {"A": {"AA": input1},
                  "B": {"BB": input2}}

        for i, input in enumerate((input0, input1, input2,
                                  input3, input4, input5, input6)):
            for var in wrf_vars:
                if var in ignore_vars:
                    continue
                test_func1 = make_test(var, input)

                setattr(WRFVarsTest, "test_{0}_input{1}_{2}".format(
                    nc_lib, i, var), test_func1)

            for method in interp_methods:
                test_interp_func1 = make_interp_test(method, input)
                setattr(WRFInterpTest, "test_{0}_input{1}_{2}".format(
                    nc_lib, i, method), test_interp_func1)

            for testid in latlon_tests:
                test_ll_func = make_latlon_test(testid, input)
                test_name = "test_{0}_input{1}_{2}".format(nc_lib, i, testid)
                setattr(WRFLatLonTest, test_name, test_ll_func)

    ut.main()
