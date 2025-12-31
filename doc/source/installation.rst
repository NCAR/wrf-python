Installation
============

Installing via Conda
---------------------

The easiest way to install wrf-python is using 
`Conda <https://docs.conda.io>`_::

    conda install conda-forge::wrf-python
    
.. note::

   If you use Conda to install WRF-Python on a supercomputer like 
   Derecho or Casper, we recommend that you do not load any Python 
   related modules via the 'module load' command. The packages installed 
   by the 'module load' system may not play nicely with packages installed 
   via Conda.
   
   Further, some systems will install Python packages to a ~/.local directory, 
   which will be found by the Conda Python interpreter and cause various 
   import problems.  If you have a ~/.local directory, we strongly suggest 
   renaming it (mv ~/.local ~/.local_backup).
    

WRF-Python on NSF NCAR HPC
--------------------------

WRF-Python is included in the `NCAR Python Library <https://ncar-hpc-docs.readthedocs.io/en/latest/environment-and-software/user-environment/package-managers/conda/#the-ncar-python-library>`_ on the NSF NCAR HPC systems or can be installed using a supported package manager of your choice.


Installing via Source Code
--------------------------

Installation via source code will require a Fortran and C compiler in order 
to run F2PY.

The source code is available on GitHub:

https://github.com/NCAR/wrf-python

To install, if you do not need OpenMP support, change your working directory 
to the wrf-python source directory and run::

    $ pip install .
    
Beginning with WRF-Python 1.1, OpenMP is supported, but preprocessing the 
ompgen.F90 file is required, which also requires running F2PY to 
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

Beginning with NumPy 1.14, F2PY extensions can now be built using the MSVC 
compiler and MinGW GFortran compiler. NumPy 1.14 is required to build 
WRF-Python for Python 3.5+. 

.. note::

   If you are building on a supercomputer and receiving linker related 
   errors (e.g. missing symbols, undefined references, etc), you probably 
   need to unset the LDFLAGS environment variable. System administrators on 
   supercomputing systems often define LDFLAGS for you so that you don't need 
   to worry about where libraries like NetCDF are installed. Unfortunately, 
   this can cause problems with the build system. To fix, 
   using the build command from above::
   
       $ unset LDFLAGS python setup.py config_fc --f90flags="-mtune=generic -fopenmp" build_ext --libraries="gomp" build
       

