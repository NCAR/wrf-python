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

To install, change your working directory 
to the wrf-python source directory and run::

    $ pip install .
 
