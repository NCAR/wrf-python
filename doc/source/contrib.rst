.. _contrib_guide:

Contributor Guide
=================================

.. note::

   This contributor guide is written for wrf-python v1.3.x. In the 
   not-too-distant future, wrf-python will undergo a significant refactoring 
   to remove the wrapt decorators (which don't serialize for dask), but the 
   concepts will remain similar to what is described in :ref:`internals`.


Introduction
-----------------------------

Thank you for your interest in contributing to the WRF-Python project. 
WRF-Python is made up of a very small amount of developers, tasked with 
supporting more than one project, so we rely on outside contributions 
to help keep the project moving forward.

The guidelines below help to ensure that the developers and outside 
collaborators remain on the same page regarding contributions. 


Source Code Location
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

For external contributors, this isn't important, other than making you aware 
that when you first clone the repository, you will be on the 
**develop** branch, which is what you should use for your development. 

Since you will be submitting pull requests for your contributions, you don't 
really need to know much about GitFlow, other than making sure that you 
are not developing off of the master branch.

  
Ways to Contribute
-----------------------------

Users are encouraged to contribute various ways. This includes:

- Submitting a bug report
- Submitting bug fixes
- Submitting new Fortran computational routines
- Submitting new Python computational routines
- Submitting fully wrapped computational routines
- Fixing documentation errors
- Creating new examples in the documentation (e.g. plotting examples)


Ground Rules
------------------------------

Please follow the code of conduct.

- Each pull request should be for a logical collection of changes. You can 
  submit multiple bug fixes in a single pull request if the bugs are related. 
  Otherwise, please submit seperate pull requests.
- Do not commit changes to files that are unrelated to your bug fix 
  (e.g. .gitignore).
- The pull request and code review process is not immediate, so please be 
  patient. 
  

Submitting Bug Reports
-----------------------------

Submitting bug reports is the easiest way to contribute. You will need to 
create an account on GitHub to submit a report.

1. Go to the issues page here: 

   https://github.com/NCAR/wrf-python/issues
   
2. Check to see if an issue has already been created for the problem that 
   you are having.
   
3. If an issue already exists for your problem, feel free to add any 
   additional information to the issue conversation.
   
4. If there is not an issue created yet for your problem, use the 
   "New Issue" button to start your new issue.
   
5. Please provide as much information as you can for the issue. Please supply 
   your version of WRF-Python you are using and which platform you are 
   using (e.g. conda-forge build on OSX). Supply a code snippet if you 
   are doing something more detailed than simply calling :meth:`wrf.getvar`.
   
6. If you are getting a crash (e.g. segmentation fault), we will most likely 
   need to see your data file if we cannot reproduce the problem here. 
   See :ref:`submitting_files`.


Setting Up Your Development Environment
---------------------------------------------

We recommend using the `conda <https://conda.io/en/latest/>`_ 
package manager for your Python environments. Our recommended setup for 
contributing is:

- Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
- Install git on your system if it is not already there (install XCode command 
  line tools on a Mac or git bash on Windows)
- Login to your GitHub account and make a fork of the
  `WRF-Python <https://github.com/ncar/wrf-python>`_ repository by clicking 
  the **Fork** button.
- Clone your fork of the WRF-Python repository (in terminal on Mac/Linux or 
  git shell/ GUI on Windows) in the location you'd like to keep it.
  
  .. code::
   
     git clone https://github.com/your-user-name/wrf-python.git

- Navigate to that folder in the terminal or in Anaconda Prompt if you're 
  on Windows.
  
  .. code::
  
     cd wrf-python
     
- Connect your repository to the NCAR WRF-Python repository. 

  .. code::
  
     git remote add ncar https://github.com/ncar/wrf-python.git
     
- To create the development environment, you'll need to run the appropriate 
  command below for your operating system.
  
  OSX:
  
  .. code::
  
     conda env create -f osx.yml
     
  Linux:
  
  .. code::
  
     conda env create -f linux.yml
     
  Win64:
  
  .. code::
  
     conda env create -f win64.yml
     
  Note: For Win64, you will also need VS2015 installed on your system.
  
- Activate your conda environment.

  .. code::
  
     conda activate develop
     
- CD to the build_scripts directory.

  .. code::
  
     cd build_scripts
     
- Build and install WRF-Python.

  OSX/Linux:
  
  .. code::
  
     sh gnu_omp.sh
     
  Windows:
  
     ./win_msvc_mingw_omp.bat
     
- The previous step will build and install WRF-Python in to the 'develop' 
  environment. If you make changes and want to rebuild, uninstall WRF-Python 
  by running:
  
  .. code::
  
     pip uninstall wrf-python
     
  Now follow the previous step to rebuild.


Pull Requests
--------------------------

In order to submit changes, you must use GitHub to issue a pull request. Your 
pull requests should be made against the **develop** branch, since we are 
following GitFlow for this project. 


Code Style
--------------------------

Python Contributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

The Python code in WRF-Python follows the 
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ coding standard. All 
Python code submitted must pass the PEP8 checks performed by the 
`pycodestyle <https://pycodestyle.readthedocs.io/en/latest/>`_ code 
style guide utility. The utility must pass without any errors or warnings.
For a tool to help automate some of the mundane formatting corrections (e.g. 
whitespace characters in blank lines, etc.), try the 
`autopep8 <https://pypi.org/project/autopep8/0.8/>`_ utility.


Fortran Contributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WRF-Python is a Fortran friendly project and we appreciate your contributions. 
However, we are only accepting Fortran 90 contributions, so you must 
convert any F77 code to F90 before contributing.

Although there is no formal style guide for Fortran contributions, Fortran 
code should look similar to other Fortran code in the WRF-Python *fortran* 
directory. 

A summary of style notes is below:

- Fortran 90 only.
- Use 4 spaces for indentation, not tabs.
- Use all capital letters for Fortran key words (e.g. IF, DO, REAL, INTENT)
- Use all capital letters for Fortran intrinsics (e.g. MAX, MIN, SUM)
- Use all capital letters for any PARAMETER constants.
- Use all lowercase letters for variables with '_' separting words 
  (snake case).
- Use all lowercase letters for functions and subroutines with '_' separting 
  words (snake case).
- Declare your REAL variables as REAL(KIND=8), unless you really need 4-byte
  REALs for a specific reason.
- Do not allocate any memory in your Fortran routine (e.g work arrays). We 
  will use numpy arrays to manage all memory. Instead, declare your work 
  array (or dynamic array) as an INOUT argument in your function 
  signature.
- Avoid submitting code that uses global variables (other than for read only 
  constants). All Fortran contributions must be threadsafe and have no side 
  effects.
- Place any computational constants in the wrf_constants module found in 
  wrf_constants.f90 and use "USE wrf_constants, ONLY : YOUR_CONSTANT" 
  declaration in your function.
- Please do not redefine constants already declared in 
  wrf_constants.f90 (e.g. G, RD, RV, etc). Although the WRF model itself 
  does not adhere to this, we are trying to be consistent with the constants 
  used throughout this project.
- Do not put any STOP statements in your code to deal with errors. STOP
  statements will bring down the entire Python interpreter with it. Instead, 
  add *errstat* and *errmsg* arguments to your function signature to tell 
  Python about the error so it can throw an exception. See WETBULBCALC
  in wrf_rip_phys_routines.f90 for how this is handled. 
- Don't worry about adding OpenMP directives to your code if you are 
  unfamiliar OpenMP, but feel free to do so if you are already familiar.








