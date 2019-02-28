from __future__ import print_function, division

import os
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, ALL_TIMES, to_np
import xarray

filename_list = ["/Users/ladwig/Documents/wrf_files/"
                 "wrf_vortex_single/wrfout_d02_2005-08-28_00:00:00",
                 "/Users/ladwig/Documents/wrf_files/wrf_vortex_single/"
                 "wrfout_d02_2005-08-28_03:00:00",
                 "/Users/ladwig/Documents/wrf_files/wrf_vortex_single/"
                 "wrfout_d02_2005-08-28_06:00:00",
                 "/Users/ladwig/Documents/wrf_files/wrf_vortex_single/"
                 "wrfout_d02_2005-08-28_09:00:00"]

result_shape = (4, 1, 96, 96)

# Let's get the first time so we can copy the metadata later
f = Dataset(filename_list[0])
# By setting squeeze to False, you'll get all the dimension names.
z1 = getvar(f, "T2", 0, squeeze=False)
xlat = getvar(f, "XLAT", 0)
xlong = getvar(f, "XLONG", 0)


z_final = np.empty(result_shape, np.float32)

# Modify this number if using more than 1 time per file
times_per_file = 1

data_times = []
xtimes = []
for timeidx in range(result_shape[0]):
    # Compute the file index and the time index inside the file
    fileidx = timeidx // times_per_file
    file_timeidx = timeidx % times_per_file

    f = Dataset(filename_list[fileidx])
    z = getvar(f, "T2", file_timeidx)
    t = getvar(f, "Times", file_timeidx)
    xt = getvar(f, "xtimes", file_timeidx)
    data_times.append(to_np(t))
    xtimes.append(to_np(xt))
    z_final[timeidx, :] = z[:]
    f.close()

# Let's make the metadata. Dimension names should copy easily if you set
# sqeeze to False, otherwise you can just set them yourself is a tuple of
# dimension names. Since you wanted
# to keep the bottom_top dimension for this 2D variable (which is normally
# removed), I'm doing this manually.
z_dims = ["Time", "bottom_top", "south_north", "west_east"]

# Xarray doesn't copy coordinates easily (it always complains about shape
# mismatches), so do this manually
z_coords = {}
z_coords["Time"] = data_times
z_coords["XTIME"] = ("Time",), xtimes
z_coords["XLAT"] = ("south_north", "west_east"), xlat
z_coords["XLONG"] = ("south_north", "west_east"), xlong
z_name = "T2"

# Attributes copy nicely
z_attrs = {}
z_attrs.update(z1.attrs)

z_with_meta = xarray.DataArray(z_final, coords=z_coords, dims=z_dims,
                               attrs=z_attrs, name=z_name)
