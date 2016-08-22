import setuptools
import numpy.distutils.core
 

ext1 = numpy.distutils.core.Extension(
    name = "wrf._wrffortran",
    sources = ["fortran/wrf_constants.f90",
               "fortran/wrf_testfunc.f90",
               "fortran/wrf_user.f90",
               "fortran/rip_cape.f90",
               "fortran/cloud_fracf.f90",
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
               "fortran/wrffortran.pyf"]
    )

with open('src/wrf/version.py') as f: 
    exec(f.read())

requirements = [
    "numpy>=1.9.0",
    "wrapt>=1.10"
    ] 

numpy.distutils.core.setup( 
    author = "Bill Ladwig",
    author_email = "ladwig@ucar.edu",
    description = "Diagnostic and interpolation routines for WRF-ARW data",
    long_description = """
    A collection of diagnostics and interpolation routines to be used with 
    WRF-ARW data.
    """,
    url = "https://github.com/NCAR/wrf-python",
    keywords = ["python", "wrf-python", "wrf", "forecast", "model",
                "weather research and forecasting", "interpolation", 
                "plotting", "plots", "meteorology", "nwp", 
                "numerical weather predication", "diagnostic"],
    install_requires = requirements,
    classifiers = ["Development Status :: 3 - Alpha",
                "Intended Audience :: Science/Research",
                "Intended Audience :: Developers",
                "License :: OSI Approved",
                "Programming Language :: Fortran",
                "Programming Language :: Python :: 2.7",
                "Programming Language :: Python :: 3.4",
                "Programming Language :: Python :: 3.5",
                "Topic :: Scientific/Engineering :: Atmospheric Science",
                "Topic :: Software Development",
                "Operating System :: POSIX",
                "Operating System :: Unix",
                "Operating System :: MacOS",
                "Operating System :: Microsoft :: Windows"],
    name = "wrf",
    version =  __version__,
    packages = setuptools.find_packages("src"),
    ext_modules = [ext1],
    package_dir = {"" : "src"},
    #namespace_packages=["wrf"],
    # Note:  If this doesn't work, you need to add the file to MANIFEST
    package_data={"wrf" : ["data/psadilookup.dat"]},
    scripts=[]
)  
