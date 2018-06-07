What's New
===========

Releases
-------------

v1.2.0
^^^^^^^^^^^^^^

- Release 1.2.0
- Previous versions of wrf-python promoted the strings used in xarray (e.g. 
  name, attributes) to Unicode strings for Python 2.7. This caused problems 
  when porting examples for PyNGL to use wrf-python in Python 3.x. All strings 
  are now the native string type for the Python version being used. While this 
  change should be transparent to most users, any users that worked with the 
  xarray name or attribute values on Python 2.7 may run in to string related 
  errors, so we've decided to bump the major version number. 


v1.1.3
^^^^^^^^^^^^^^

- Release 1.1.3
- Fixed/Enhanced the cloud top temperature diagnostic.
   - Optical depth was not being calculated correctly when 
     cloud ice mixing ratio was not available.
   - Fixed an indexing bug that caused crashes on Windows, but should have been 
     crashing on all platforms.
   - Users can now specify if they want cloud free regions to use fill values,
     rather than the default behavior of using the surface temperature.
   - Users can now specify the optical depth required to trigger the cloud
     top temperature calculation. However, the default value of 1.0 should be 
     sufficient for most users.
- Added 'th' alias for the theta product.
- Fixed a crash issue related to updraft helicity when a dictionary is 
  used as the input.
- Dictionary inputs now work correctly with xy_to_ll and ll_to_xy.
- The cape_2d diagnostic can now work with a single column of data, just like 
  cape_3d.
  

v1.1.2
^^^^^^^^^^^^^^

- Release 1.1.2
- Fix OpenMP directive issue with cloud top temperature.


v1.1.1
^^^^^^^^^^^^^^

- Release 1.1.1
- Added script for building on Cheyenne with maxed out Intel settings, which 
  also required a patch for numpy.distutils.
- Fixed a few unicode characters hiding in a docstring that were causing 
  problems on Cheyenne, and also building the docs with Sphinx on Python 2.x.
- Fix issue with np.amax not working with xarray on Cheyenne, causing an error
  with the mdbz product.
- Fix cape_2d private variable bug when running with multiple CPUs.


v1.1.0
^^^^^^^^^^^^^^

- Release 1.1.0
- Computational routines now support multiple cores using OpenMP.  See 
  :ref:`using_omp` for details on how to use this new feature.
- The CAPE routines should be noticeably faster, even in the single threaded 
  case (thank you supreethms1809!).
- :meth:`wrf.getvar` now works correctly with non-gridded NetCDF variables
- The cloud fraction diagnostic has changed:
   - Users can now select their own cloud threshold levels, and can choose 
     between a vertical coordinate defined as height (AGL), height (MSL), or 
     pressure. 
   - The default vertical coordinate type has been changed to be height (AGL). 
     This ensures that clouds appear over mountainous regions. If you need 
     the old behavior, set the *vert_type* argument to 'pressure'.
   - Fixed a bug involving the cloud threshold search algorithm, where if the 
     surface was higher than the threshold for a cloud level, the algorithm
     would use whatever was there before (uninitialized variable bug). This 
     caused some interesting visualization issues when plotted.  Now, whenever 
     the surface is above a cloud level threshold, a fill value is used to 
     indicate that data is unavailable for that location.
- The cartopy object for LambertConformal should now work correctly in the 
  southern hemisphere.
- Fixed a bug with the PolarStereographic projection missing a geobounds 
  argument (thank you hanschen!).
- Renamed the modules containing the 'get_product' routines used 
  by :meth:`wrf.getvar` to avoid naming conflicts with the raw computational 
  routine names. Users should be using :meth:`wrf.getvar` instead of these 
  routines, but for those that imported the 'get_product' routines 
  directly, you will need to modify your code.
- Fixed a uniqueness issue with the internal coordinate cache that was causing
  crashes when input data is changed to a different file in a jupyter notebook 
  cell.
- Added code to better support building wheels on Windows (thank you letmaik!)
- Improved support for scipy.io.netcdf objects. 
- Added a new 'zstag' diagnostic that returns the height values for the 
  vertically staggered grid.
- A DOI is now available for wrf-python. Please cite wrf-python if you are 
  using it for your research. (See :ref:`citation`)
- Fixed issue with vertcross and interpline not working correctly when a 
  projection object is used. Users will now have to supply the lower left 
  latitude and longitude corner point.
- Beginning with numpy 1.14, wrf-python can be built using the MSVC 
  compiler with gfortran. WRF-Python can now be built for Python 3.5+ on 
  services like AppVeyor.


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
  


