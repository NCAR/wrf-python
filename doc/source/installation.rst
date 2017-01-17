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

    $ conda install -c conda-forge wrf-python
    
    
Installing via Source Code
--------------------------

Installation via source code will require a Fortran and C compiler in order 
to run f2py.  You can get them
`here <https://gcc.gnu.org/wiki/GFortranBinaries>`_.

The source code is available via github:

https://github.com/NCAR/wrf-python

To install, change to the wrf-python directory and run::

    $ pip install .


