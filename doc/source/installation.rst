Installation
============

Required Dependencies
----------------------

    - Python 2.7, 3.4, or 3.5
    - numpy (1.9 or later)
    - wrapt (1.10 or later)


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

To install, change to the wrf-python directory and run::

    $ pip install .

Note that building on Win64 with Python 3.5+ and the mingw-64 compiler
is very difficult, due to incompatibilities with the runtime libraries and 
lack of support from numpy's distutils. Improved support for these 
configurations, along with numpy distutils support, should take place this 
year.  But for now, visual studio and the intel compiler may be required.  
Otherwise, Python 2.7 or Python 3.4 is recommended. 
