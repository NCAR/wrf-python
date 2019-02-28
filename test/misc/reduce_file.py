from netCDF4 import Dataset as NetCDF

f = "/Users/ladwig/Documents/wrf_files/wrfout_d01_2016-02-25_18_00_00"
outfilename = "/Users/ladwig/Documents/wrf_files/rotated_pole_test.nc"

in_nc = NetCDF(f, mode='r', format="NETCDF3_CLASSIC")
out_nc = NetCDF(outfilename, mode='w', format="NETCDF3_CLASSIC")

# Copy Global Attributes
for att_name in in_nc.ncattrs():
    setattr(out_nc, att_name, getattr(in_nc, att_name))

# Copy Dimensions, but modify the vertical dimensions
for dimname, dim in in_nc.dimensions.iteritems():
    out_nc.createDimension(dimname, len(dim))

# Copy Variables and their Attributes, using the reduced vertical dimension
for varname, var in in_nc.variables.iteritems():
    if varname in ("T2", "XLAT", "XLONG", "XTIME"):
        datatype = var.datatype
        dimensions = var.dimensions
        shape = var.shape

        new_shape = shape

        new_var = out_nc.createVariable(varname, datatype, dimensions)

        new_var[:] = var[:]

        for att_name in var.ncattrs():
            setattr(new_var, att_name, getattr(var, att_name))

in_nc.close()
out_nc.close()
