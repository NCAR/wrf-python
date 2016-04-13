import setuptools
import numpy.distutils.core
 
ext1 = numpy.distutils.core.Extension(
    name = "wrf.var._wrfext",
    sources = ["src/python/wrf/var/wrfext.f90",
               "src/python/wrf/var/wrfext.pyf"]
    )

ext2 = numpy.distutils.core.Extension(
    name = "wrf.var._wrfcape",
    sources = ["src/python/wrf/var/wrfcape.f90",
               "src/python/wrf/var/wrfcape.pyf"]
    )

numpy.distutils.core.setup( 
    name = "wrf.var",
    version = "0.0.1",        
    packages = setuptools.find_packages("src/python"),   
    ext_modules = [ext1,ext2],
    package_dir={"":"src/python"},
    namespace_packages=["wrf"],
    scripts=[],
)  
