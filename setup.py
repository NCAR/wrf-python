import setuptools
import numpy.distutils.core
 
#ext1 = numpy.distutils.core.Extension(
#    name = "wrf._wrfext",
#    sources = ["src/wrf/wrfext.f90",
#               "src/wrf/wrfext.pyf"]
#    )

#ext2 = numpy.distutils.core.Extension(
#    name = "wrf._wrfcape",
#    sources = ["src/wrf/wrfcape.f90",
#               "src/wrf/wrfcape.pyf"]
#    )

ext1 = numpy.distutils.core.Extension(
    name = "wrf._wrffortran",
    sources = ["fortran/constants.f90",
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

numpy.distutils.core.setup( 
    name = "wrf",
    version = "0.0.1",        
    packages = setuptools.find_packages("src"),   
    ext_modules = [ext1],
    package_dir={"" : "src"},
    #namespace_packages=["wrf"],
    # Note:  If this doesn't work, you need to add the file to MANIFEST
    package_data={"wrf" : ["data/psadilookup.dat"]},
    scripts=[]
)  
