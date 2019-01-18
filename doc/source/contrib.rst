.. _contrib_guide:

Contributor Guide
=================================

.. note::

   This contributor guide is written for wrf-python v1.3.x. In the 
   not-too-distant future, wrf-python will undergo a significant refactoring 
   to remove the wrapt decorators (which don't serialize for dask), but the 
   concepts will remain the same as described below.

  
Ways to Contribute
-----------------------------

Users are encouraged to contribute various ways. This includes:

- Submitting a bug report
- Submitting bug fixes
- Submitting new Fortran computational routines
- Submitting new Python computational routines
- Submitting fully wrapped computational routines


Getting the source code
------------------------------

The source code is available on GitHub:

    https://github.com/NCAR/wrf-python

To checkout the code::

    git clone https://github.com/NCAR/wrf-python    


Git Flow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This project follows the GitFlow Workflow, which you can read about here if it
is new to you:

https://leanpub.com/git-flow/read

When you first clone the repository, by default you will be on the 'develop' 
branch, which is what you should use for your development. 


Pull Requests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to submit changes, you must use GitHub to issue a pull request. 


Overview of WRF-Python Internals
----------------------------------

WRF-Python is a collection of diagnostic and interpolation routines for WRF-ARW
data. The API consists of a handful of functions 



