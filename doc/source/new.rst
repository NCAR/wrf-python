What's New
===========

Releases
-------------

v1.0.5
^^^^^^^^^^^^^^

- Release 1.0.5
- Reduced the CI test file sizes by half.  


v1.0.4
^^^^^^^^^^^^^^

- Release 1.0.4
- Fix warnings with CI tests which were caused by fill values being written 
  as NaN to the NetCDF result file.
- Added the __eq__ operator to the WrfProj projection base class.
- Fixed array order issue when using the raw CAPE routine with 1D arrays.

v1.0.3
^^^^^^^^^^^^^^

- Relase 1.0.3
- Fixed an issue with the cartopy Mercator subclass where the xlimits were 
  being calculated to the same value (or very close), causing blank plots.

v1.0.2
^^^^^^^^^^^^^^

- Release 1.0.2
- Fixed issue with the wspd_wdir product types when sequences of files are 
  used.


v1.0.1
^^^^^^^^^^^^^

- Release 1.0.1
- Fixed issue with initialization of PolarStereographic and LatLon map 
  projection objects.
- Fixed issue where XTIME could be included in the coordinate list of a 
  variable, but the actual XTIME variable could be missing.  NCL allows this,
  so wrf-python should as well.
  

v1.0.0
^^^^^^^^^^^^^

- Release 1.0.0.
- Fixed issue with not being able to set the thread-local coordinate cache to 
  0 to disable it.  Also, the cache will now correctly resize itself when 
  the size is reduced to less than its current setting.
- Fixed an issue with the '0000-00-00 00:00:00' time used in geo_em files 
  causing crashes due to the invalid time.  The time is now set to 
  numpy.datetime64('NaT').
- Fixed issue with wrf.cape_3d not working correctly with a single 
  column of data.


Beta Releases
--------------

v1.0b3
^^^^^^^^^^^^^

- Beta release 3.
- Improvements made for conda-forge integration testing.
- Fixed an incorrectly initialized variable issue with vinterp.  This issue 
  mainly impacts the unit tests for continuous integration testing with 
  conda-forge, since the data set used for these tests is heavily cropped.
- Back-ported the inspect.BoundArguments.apply_defaults so that Python 3.4
  works.  Windows users that want to try out wrf-python with Python 3.4
  can use the bladwig conda channel to get it.

v1.0b2
^^^^^^^^^^^^^^

- Beta release 2.
- xarray 0.9 no longer includes default index dimensions in the coordinate 
  mappings.  This was causing a crash in the routines that cause a reduction
  in dimension shape, mainly the interpolation routines.  This has been 
  fixed.
- Documentation updated to show the new output from xarray.

v1.0b1
^^^^^^^^^^^^^

- Beta release 1.
- Added more packaging boilerplate.
- Note:  Currently unable to build with Python 3.5 on Windows, due to
  issues with distutils, numpy distutils, and mingw compiler.  Will attempt
  to find a workaround before the next release. Windows users should use 
  Python 2.7 or Python 3.4 for now.


----------------

Alpha Releases
----------------

v1.0a3
^^^^^^^^^^^^

- Alpha release 3.
- Added docstrings.
- The mapping API has changed.
    - The projection attributes are no longer arrays for moving domains.
    - Utility functions have been added for extracting geobounds.  It is now 
      easier to get map projection objects from sliced variables.
    - Utility functions have been added for getting cartopy, basemap, and pyngl
      objects.
    - Users should no longer need to use xarray attributes directly
- Now uses CoordPair for cross sections so that lat/lon can be used instead of 
  raw x,y grid coordinates.
- Renamed npvalues to to_np which is more intuitive.
- Fixed issue with generator expressions.
- Renamed some functions and arguments.


-------------

  
Known Issues
--------------

v1.0.0
^^^^^^^^

- Currently unable to build on Windows with Python 3.5+ using open source 
  mingw compiler.  The mingwpy project is working on resolving the 
  incompatibilities between mingw and Visual Studio 2015 that was used to 
  build Python 3.5+.  Numpy 1.13 also has improved f2py support for 
  Python 3.5+ on Windows, so this will be revisited when it is released.
  


