.. _contrib_guide:

Contributor Guide
=================================

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
are not developing off of the main branch.

  
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

Please follow the `Code of Conduct <https://github.com/NCAR/wrf-python/blob/develop/CODE_OF_CONDUCT.md>`_.

- Please create an issue on GitHub for any pull request you wish to submit, 
  except for documentation issues.
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
   

Submitting Fortran Computational Routines
--------------------------------------------

If you have Fortran computational routines that you'd like to contribute,
but don't know how to wrap them in to Python, please follow the instructions
below.

1. Only Fortran 90 code will be accepted, so please port your F77 code to 
   F90.
   
2. Follow the :ref:`fortranstyle`.

3. Please only submit routines relevant to WRF-Python (e.g. diagnostics, 
   interpolation). General purpose climate/meteorology should go in to the 
   SkyLab project (a project providing similar functionality as 
   NCL).
   
4. If you are unsure if you should contribute your Fortran code, make an 
   issue on GitHub and we can begin a discussion there. 
   
5. Place your code in the fortran/contrib directory in the WRF-Python 
   source tree.
   
6. Document your code with a text file that has the same name as your Fortran 
   file, but ending in .rst. This file should placed with your F90 code 
   in the fortran/contrib directory. Your documentation can use 
   restructured text formatting, or just plain text. This documentation 
   will be used for the docstring when Python wrappers are made.

7. If you are unable to provide any type of test for your routine, please 
   ensure that your documentation describes what your computation 
   should produce. You can submit auxiallary documentation and/or images for 
   this purpose if needed.


Submitting Python Computational Routines
---------------------------------------------

If you would like to submit a computational routine written in Python, but 
don't know how to integrate it with the rest of WRF-Python's internals 
(e.g. left indexing, arg checking, etc), feel free to 
submit the pure Python routine. Below is the guide for submitting pure 
Python routines.

1. These routines should be placed in src/wrf/contrib.py. These algorithms 
   will not be imported in to WRF-Python's default namespace.
   
2. Follow the :ref:`pythonstyle`. 
  
2. Write your computation as dimension unaware as possible. For example, 
   adding pressure and perturbation pressure is simply P + PB.
   
3. If dimensionality is needed, then write for the minimum dimensionality 
   required to make the computation for one time step (if applicable). For 
   example, if you're computing CAPE, then you should use three dimensions for 
   your algorithm, and we will handle the looping over all times.
   
4. Document your routine by creating a docstring that follows Google docstring 
   format (see `Sphinx Napoleon <https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html#google-vs-numpy>`_).
   
5. If you are unable to provide a test for this function, please provide 
   additional documentation (or images) to show what this function should 
   produce.


Submitting Fully Wrapped Computational Routines
---------------------------------------------------

Submitting a fully wrapped computational routines is the fastest way to get 
your contributation released. However, it requires the most effort on your 
part. 

1. Read the :ref:`internals` guide. This will show you how to wrap your
   routine.

2. Follow the :ref:`fortranstyle` and :ref:`pythonstyle`.

3. You should create your contribution in the WRF-Python source tree as if 
   you were one of the core developers of it. This means:
   
   - Your Fortran code (if applicable) should be placed in the *fortran*
     directory.
   
   - Update the "ext1 = numpy.distutils.core.Extension" section of *setup.py* 
     to include your new Fortran source (if applicable).
     
   - Update *extension.py* to create the Python wrapper that calls your 
     Fortran function. This must include the appropriate function decorators
     for handling argument checking, leftmost dimension indexing, etc. as 
     described in :ref:`internals`.
     
   - If the current function decorators do not cover your specific needs, 
     place your custom decorator in *specialdec.py*.  Most of the decorators 
     in specialdec.py are used for products that contain multiple outputs like 
     cape_2d, but you can use it for other purposes.
     
   - If your function is pure python, you can create a new module for it, 
     or place it in another module with similar functionality. For example, 
     if your routine is a new interpolation routine, then it should go 
     in interp.py. Remember to apply the same type of decorators as 
     done with Fortran extensions (checking args, leftmost dimension 
     indexing, etc).
     
   - Create a 'getter' routine which is responsible for extracting the 
     required variables from a WRF file and calling your computational 
     routine. This is what will be called by :meth:`wrf.getvar`. 
     This function should be placed in a new python module with the prefix 
     'g\_' (i.e. g_yourdiagnostic.py).
     
   - Decorate your getter routine with an appropriate metadata handling 
     decorator. If you need to make a custom decorator for the metadata, 
     place it in *metadecorators.py*. 
     
   - Update the mappings in *routines.py* to map your diagnostic name to your 
     function, and to declare any keyword arguments that your function 
     needs aside from the usual wrfin, varname, timeidx, method, 
     squeeze, cache, and meta.
     
   - If you would like to make your routine available as a raw computation,
     you will need to place an additional thin wrapper in *computation.py*. 
     This thin wrapper must be decorated with an appropriate metadata decorator 
     found in *metadecorators.py* (usually set_alg_metadata). If you need to 
     write your own custom metadata decorator, write it in *metadecorators.py*.
  
   - You must provide a docstring for every function you create using 
     Google docstring format (see `Sphinx Napoleon <https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html#google-vs-numpy>`_).
    
   - You must provide a test for your function. See :ref:`testing`.
   

Fixing Documentation Errors
--------------------------------------

1. Documenation is made with Sphinx using restructured text.

2. Python docstrings follow `Google docstring <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_ format.

3. Documentation can be found in the *doc* directory, along with the 
   docstrings contained within the Python code.
   
4. For documentation fixes, you can just submit a pull request with the 
   appropriate corrections already made.


Creating New Examples
--------------------------------------

1. Examples are made with Sphinx using restructured text.

2. Examples are currently found in the *doc* directory, mostly within the 
   *basic_usage.rst* and *plot.rst* files. Feel free to contribute more 
   examples to these files.
   
3. Unless you are drastically changing the documentation structure, you can 
   submit a pull request with your examples without creating a GitHub 
   issue. If you are making a large change, or are unsure about it, then 
   go ahead and create a GitHub issue to discuss with the developers.
   

.. _dev_setup:

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
  
     conda env create -f build_envs/Darwin.yml
     
  Linux:
  
  .. code::
  
     conda env create -f build_envs/Linux.yml
     
  Win64:
  
  .. code::
  
     conda env create -f build_envs/Win64.yml
     
  Note: For Win64, you will also need VS2015 installed on your system.
  
- Activate your conda environment.

  .. code::
  
     conda activate wrf_python_build
     
- CD to the build_scripts directory.

  .. code::
  
     cd build_scripts
     
- Build and install WRF-Python.

  OSX/Linux:
  
  .. code::
  
     sh gnu_omp.sh
     
  Windows:
  
  .. code::
  
     ./win_msvc_mingw_omp.bat
     
- The previous step will build and install WRF-Python in to the 'develop' 
  environment. If you make changes and want to rebuild, uninstall WRF-Python 
  by running:
  
  .. code::
  
     pip uninstall wrf-python
     
  Now follow the previous step to rebuild.



Code Style
--------------------------

.. _pythonstyle:

Python Style Guide
^^^^^^^^^^^^^^^^^^^^^^^^^^

The Python code in WRF-Python follows the 
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ coding standard. All 
Python code submitted must pass the PEP8 checks performed by the 
`pycodestyle <https://pycodestyle.readthedocs.io/en/latest/>`_ code 
style guide utility. The utility must pass without any errors or warnings.
For a tool to help automate some of the mundane formatting corrections (e.g. 
whitespace characters in blank lines, etc.), try the 
`autopep8 <https://pypi.org/project/autopep8/0.8/>`_ utility.


.. _fortranstyle:

Fortran Style Guide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At this time, we are only accepting Fortran 90 contributions, so you must 
convert any F77 code to F90 before contributing.

Although there is no formal style guide for Fortran contributions, Fortran 
code should look similar to other Fortran code in the WRF-Python *fortran* 
directory. 

A summary of style notes is below:

- Fortran 90 only.
- Use 4 spaces for indentation, not tabs.
- Use all capital letters for Fortran key words (e.g. IF, DO, REAL, INTENT)
- Use all capital letters for Fortran intrinsics (e.g. MAX, MIN, SUM)
- Use all capital letters for any constants declared as PARAMETER (e.g. RD).
- Use all lowercase letters for variables with '_' separting words 
  (snake case).
- Use all lowercase letters for functions and subroutines using '_' to 
  separate words (snake case).
- Declare your REAL variables as REAL(KIND=8), unless you really need 4-byte
  REALs for a specific reason.
- Do not allocate any memory in your Fortran routine (e.g work arrays). We 
  will use numpy arrays to manage all memory. Instead, declare your work 
  array (or dynamic array) as an INTENT(INOUT) argument in your function 
  signature.
- Avoid submitting code that uses global variables (other than for read only 
  constants). All Fortran contributions must be threadsafe and have no side 
  effects.
- Place any computational constants in the wrf_constants module found in 
  *wrf_constants.f90* and put a "USE wrf_constants, ONLY : YOUR_CONSTANT,..." 
  declaration in your function.
- Please do not redefine constants already declared in 
  wrf_constants.f90 (e.g. G, RD, RV, etc). Although the WRF model itself 
  does not adhere to this, we are trying to be consistent with the constants 
  used throughout WRF-Python.
- Do not put any STOP statements in your code to handle errors. STOP
  statements will bring the entire Python interpreter down with it. Instead, 
  use *errstat* and *errmsg* arguments to your function signature to tell 
  Python about the error so it can throw an exception. See WETBULBCALC
  in *wrf_rip_phys_routines.f90* for how this is handled. 
- If you know how to use OpenMP directives, feel free to add them to your 
  routine, but this is not required.


Pull Requests
--------------------------

In order to submit changes, you must use GitHub to issue a pull request. Your 
pull requests should be made against the **develop** branch, since we are 
following GitFlow for this project. 

Following a pull request, automated continuous integration tools will be 
run to ensure that your code follows the PEP 8 style guide, and verifies that
a basic suite of unit tests run. 

If your pull request is for a bug fix to an existing computational routine, 
then the automated unit tests will probably fail due to the new values. This 
is not a problem, but be sure to indicate to the developers in your GitHub 
issue that the tests will need to be updated.

.. testing_::

Tests
---------------------------

Once you have submitted your contribution, we need a way to test your 
code. Currently, most of WRF-Python's tests are written to ensure that 
WRF-Python produces the same result as the NCAR Command Language (NCL), which 
is where the code was originally derived. However, this isn't applicable for 
new contributions and bug fixes, since there is nothing to test against for 
new contributions and bug fixes might change the numerical result. So, we have 
some recommendations below for how you can create your own tests.

Sample Data
^^^^^^^^^^^^^^^^^^^

You can download sample data for Hurricane Katrina here: <contact developers>
This data includes both a moving nest and a static nest version. You should 
create your tests with this data set (both static and moving nests), unless 
you are unable to reproduce a particular problem with it.

Supplying Data
^^^^^^^^^^^^^^^^^^^^^^

If you need to supply us data for your test (due to a bug) please provide us a 
link to the file or upload it using :ref:`submitting_files`. 
Unless the data is very small, do not add it to the GitHub 
repository.

If you can demonstrate the problem/solution with a minimal set of hand created 
values, then you can use that in your test.


Guidelines
^^^^^^^^^^^^^^^^^^^

The following are guidelines for testing you contributions. The developers are 
aware that some issues have unique needs, so you can use the GitHub 
issue related to your contribution to discuss with developers.

1. New computations must work for both moving nests and static nests. 
   Generally this is not an issue unless your data makes use of lat/lon 
   information (e.g. cross sections with lat/lon line definitions).

2. WRF-Python's tests can be found in the *test* directory.
   
3. WRF-Python's tests were written using the standard *unittest* package,
   along with numpy's test package for the assert fuctions. One 
   reason for this is that many of the tests are dynamically generated, and 
   other testing frameworks can't find the tests when generated this way.
   If you need to use another test framework, that's fine, just let us know 
   in your GitHub issue.
   
4. Place your test in the test/contrib directory.

5. For new contributions, images may be sufficient to show that your 
   code is working. Please discuss with the developers in you GitHub issue.
   
6. For bug related issues, try to create a case that demonstrates the problem, 
   and demonstrates the fix. If your problem is a crash, then proving that 
   your code runs without crashing should be sufficient. 
   
7. You might need some creativity here.





