import os
import setuptools
import socket

if not socket.gethostname().startswith("cheyenne"):
    import numpy.distutils.core
else:
    import chey_intel
    import numpy.distutils.core
    import numpy.distutils.fcompiler.intel

    numpy.distutils.fcompiler.intel.IntelFCompiler = chey_intel.IntelFCompiler
    numpy.distutils.fcompiler.intel.IntelVisualFCompiler = (
        chey_intel.IntelVisualFCompiler)
    numpy.distutils.fcompiler.intel.IntelItaniumFCompiler = (
        chey_intel.IntelItaniumFCompiler)
    numpy.distutils.fcompiler.intel.IntelItaniumVisualFCompiler = (
        chey_intel.IntelItaniumVisualFCompiler)
    numpy.distutils.fcompiler.intel.IntelEM64VisualFCompiler = (
        chey_intel.IntelEM64VisualFCompiler)
    numpy.distutils.fcompiler.intel.IntelEM64TFCompiler = (
        chey_intel.IntelEM64TFCompiler)

ext1 = numpy.distutils.core.Extension(
    name="wrf._wrffortran",
    sources=["fortran/wrf_constants.f90",
             "fortran/wrf_testfunc.f90",
             "fortran/wrf_user.f90",
             "fortran/rip_cape.f90",
             "fortran/wrf_cloud_fracf.f90",
             "fortran/wrf_fctt.f90",
             "fortran/wrf_user_dbz.f90",
             "fortran/wrf_relhl.f90",
             "fortran/calc_uh.f90",
             "fortran/wrf_user_latlon_routines.f90",
             "fortran/wrf_pvo.f90",
             "fortran/eqthecalc.f90",
             "fortran/wrf_rip_phys_routines.f90",
             "fortran/wrf_pw.f90",
             "fortran/wrf_vinterp.f90",
             "fortran/wrf_wind.f90",
             "fortran/omp.f90"]
    )

on_rtd = os.environ.get("READTHEDOCS", None) == "True"
# on_rtd=True
if on_rtd:
    requirements = ["mock; python_version < 3.3"]
    ext_modules = []

else:
    # Place install_requires into the text file "requirements.txt"
    with open("requirements.txt") as f2:
        requirements = f2.read().strip().splitlines()

    ext_modules = [ext1]

numpy.distutils.core.setup(
    install_requires=requirements,
    ext_modules=ext_modules,
    scripts=[]
)
