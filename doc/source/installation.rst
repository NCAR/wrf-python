Installation
============

Required Dependencies
----------------------

    - Python 2.7, 3.4, or 3.5+
    - numpy (1.11 or later; 1.14 required to build on Windows)
    - wrapt (1.10 or later)
    - setuptools (38.0 or later) 


Highly Recommended Packages
----------------------------

    - xarray (0.7.0 or later)
    - PyNIO (1.4.3 or later)
    - netCDF4-python (1.2.0 or later)


Plotting Packages
-------------------------

    - PyNGL (1.4.3 or later)
    - matplotlib (1.4.3 or later)
        - cartopy (0.13 or later)
        - basemap (1.0.8 or later)


Installing via Conda
---------------------

The easiest way to install wrf-python is using 
`Conda <http://conda.pydata.org/docs/>`_::

    conda install -c conda-forge wrf-python
    
.. note::

   If you use conda to install wrf-python on a supercomputer like 
   Yellowstone or Cheyenne, we recommend that you do not load any python 
   related modules via the 'module load' command. The packages installed 
   by the 'module load' system will not play nicely with packages installed 
   via conda.
   
   Further, some systems will install python packages to a ~/.local directory, 
   which will be found by the miniconda python interpreter and cause various 
   import problems.  If you have a ~/.local directory, we strongly suggest 
   renaming it (mv ~/.local ~/.local_backup).
    

Installing on Yellowstone
----------------------------

On Yellowstone, wrf-python can also be installed using the module load system, 
if this is preferred over using conda.

Unfortunately, because wrf-python requires newer dependencies, it is not 
available using the 'all-python-libs' module, so many of the dependencies 
need to be manually installed (most are for xarray).

Also, make sure you are running in the gnu/4.8.2 compiler environment or 
you will get import errors for a missing libquadmath library when you 
go to import wrf-python.  

To install::

    module load gnu/4.8.2 or module swap intel gnu/4.8.2
    module load python/2.7.7
    module load numpy/1.11.0 wrapt/1.10.10 scipy/0.17.1 bottleneck/1.1.0 numexpr/2.6.0 pyside/1.1.2 matplotlib/1.5.1 pandas/0.18.1 netcdf4python/1.2.4 xarray/0.8.2
    module load wrf-python/1.0.1


Installing via Source Code
--------------------------

Installation via source code will require a Fortran and C compiler in order 
to run f2py.  You can get them
`here <https://gcc.gnu.org/wiki/GFortranBinaries>`_.

The source code is available via github:

https://github.com/NCAR/wrf-python

Or PyPI:

https://pypi.python.org/pypi/wrf-python

To install, if you do not need OpenMP support, change your working directory 
to the wrf-python source directory and run::

    $ pip install .
    
Beginning with wrf-python 1.1, OpenMP is supported, but preprocessing the 
ompgen.F90 file is required, which also requires running f2py to 
build the .pyf file. To simplify this process, you can use the build scripts in 
the *build_scripts* directory as a guide, or just call the script directly.

Below is a sample from a build script for GNU compiler with OpenMP enabled:

.. code-block:: none

   cd ../fortran/build_help
   
   gfortran -o sizes -fopenmp omp_sizes.f90
   
   python sub_sizes.py

   cd ..
   
   gfortran -E ompgen.F90 -fopenmp -cpp -o omp.f90
   
   f2py *.f90 -m _wrffortran -h wrffortran.pyf --overwrite-signature
   
   cd ..

   python setup.py clean --all
   
   python setup.py config_fc --f90flags="-mtune=generic -fopenmp" build_ext --libraries="gomp" build
   
   pip install .

Beginning with numpy 1.14, f2py extensions can now be built using the MSVC 
compiler and mingw gfortran compiler. Numpy 1.14 is required to build 
wrf-python for Python 3.5+. 

.. note::

   If you are building on a supercomputer and receiving linker related 
   errors (e.g. missing symbols, undefined references, etc), you probably 
   need to unset the LDFLAGS environment variable. System administrators on 
   supercomputing systems often define LDFLAGS for you so that you don't need 
   to worry about where libraries like NetCDF are installed. Unfortunately, 
   this can cause problems with the numpy.distutils build system. To fix, 
   using the build command from above::
   
       $ unset LDFLAGS python setup.py config_fc --f90flags="-mtune=generic -fopenmp" build_ext --libraries="gomp" build
       

