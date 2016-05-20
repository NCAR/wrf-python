import setuptools
import numpy.distutils.core
 
ext1 = numpy.distutils.core.Extension(
    name = "wrf._wrfext",
    sources = ["src/python/wrf/wrfext.f90",
               "src/python/wrf/wrfext.pyf"]
    )

ext2 = numpy.distutils.core.Extension(
    name = "wrf._wrfcape",
    sources = ["src/python/wrf/wrfcape.f90",
               "src/python/wrf/wrfcape.pyf"]
    )

numpy.distutils.core.setup( 
    name = "wrf",
    version = "0.0.1",        
    packages = setuptools.find_packages("src/python"),   
    ext_modules = [ext1,ext2],
    package_dir={"":"src/python"},
    #namespace_packages=["wrf"],
    scripts=[],
)  
