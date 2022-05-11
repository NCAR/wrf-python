import os
import sys
import setuptools
import socket

# Bootstrap a numpy installation before trying to import it.
import importlib
try:
    importlib.util.find_spec('numpy')
except ImportError:
    import subprocess
    subprocess.call([sys.executable, '-m', 'pip', 'install', 'numpy'])

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

with open("src/wrf/version.py") as f:
    exec(f.read())

on_rtd = os.environ.get("READTHEDOCS", None) == "True"
# on_rtd=True
if on_rtd:
    if sys.version_info < (3, 3):
        requirements = ["mock"]  # for python2 and python < 3.3
    else:
        requirements = []  # for >= python3.3
    ext_modules = []

else:
    # Place install_requires into the text file "requirements.txt"
    with open("requirements.txt") as f2:
        requirements = f2.read().strip().splitlines()

        # if sys.version_info < (3,3):
        #     requirements.append("mock")
    ext_modules = [ext1]

numpy.distutils.core.setup(
    name='wrf-python',
    author="Bill Ladwig",
    maintainer="GeoCAT",
    maintainer_email="geocat@ucar.edu",
    description="Diagnostic and interpolation routines for WRF-ARW data.",
    long_description=("A collection of diagnostic and interpolation "
                      "routines to be used with WRF-ARW data.\n\n"
                      "GitHub Repository:\n\n"
                      "https://github.com/NCAR/wrf-python\n\n"
                      "Documentation:\n\n"
                      "http://wrf-python.rtfd.org\n"),
    url="https://github.com/NCAR/wrf-python",
    version=__version__,
    package_dir={"": "src"},
    keywords=["python", "wrf-python", "wrf", "forecast", "model",
              "weather research and forecasting", "interpolation",
              "plotting", "plots", "meteorology", "nwp",
              "numerical weather prediction", "diagnostic",
              "science", "numpy"],
    python_requires='>=3.7',
    install_requires=requirements,
    classifiers=["Development Status :: 5 - Production/Stable",
                 "Intended Audience :: Science/Research",
                 "Intended Audience :: Developers",
                 "License :: OSI Approved :: Apache Software License",
                 "Programming Language :: Fortran",
                 "Programming Language :: Python :: 3.7",
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: 3.9',
                 "Topic :: Scientific/Engineering :: Atmospheric Science",
                 "Topic :: Software Development",
                 "Operating System :: POSIX",
                 "Operating System :: Unix",
                 "Operating System :: MacOS",
                 "Operating System :: Microsoft :: Windows"],
    platforms=["any"],
    license="Apache License 2.0",
    packages=setuptools.find_packages("src"),
    ext_modules=ext_modules,
    download_url="http://python.org/pypi/wrf-python",
    package_data={"wrf": ["data/psadilookup.dat"]},
    scripts=[]
)
