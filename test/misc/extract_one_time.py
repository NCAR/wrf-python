
from netCDF4 import Dataset

VARS_TO_KEEP = ("Times", "XLAT", "XLONG", "XLAT_U", "XLAT_V", "XLONG_U",
                "XLONG_V", "U", "V", "W", "PH", "PHB", "T", "P", "PB", "Q2",
                "T2", "PSFC", "U10", "V10", "XTIME", "QVAPOR", "QCLOUD",
                "QGRAUP", "QRAIN", "QSNOW", "QICE", "MAPFAC_M", "MAPFAC_U",
                "MAPFAC_V", "F", "HGT", "RAINC", "RAINSH", "RAINNC",
                "I_RAINC", "I_RAINNC")

INPUT_FILE = "wrfout_d02_2005-08-28_00:00:00"
OUTPUT_FILE = "wrfout_problem_file"
DESIRED_TIME_INDEX = 0

with Dataset(INPUT_FILE) as infile, Dataset(OUTPUT_FILE, mode="w") as outfile:
    # Copy the global attributes
    outfile.setncatts(infile.__dict__)

    # Copy Dimensions
    for name, dimension in infile.dimensions.items():
        if name != "Time":
            dimsize = len(dimension)
            outfile.createDimension(name, dimsize)
        else:
            outfile.createDimension(name, 1)

    # Copy Variables
    for name, variable in infile.variables.iteritems():
        if name not in VARS_TO_KEEP:
            continue

        outvar = outfile.createVariable(name, variable.datatype,
                                        variable.dimensions, zlib=True)

        if len(variable.dimensions) > 1:
            outvar[:] = variable[DESIRED_TIME_INDEX, :]
        else:
            outvar[:] = variable[DESIRED_TIME_INDEX]

        outvar.setncatts(variable.__dict__)
