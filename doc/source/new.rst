What's New
===========

Releases
-------------

v1.3.2 (February 2019)
^^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.3.2
- Coordinate name index positions are no longer assumed and are searched 
  instead. Some users use xarray to rewrite WRF output files, and xarray 
  might reorder the coordinate name positions.
- Fixed a segfault issue with CAPE when more than 150 vertical levels are 
  used (e.g. LES runs).
- setup.py will now bootstrap the numpy installation (thanks bbonenfant!).


v1.3.1 (January 2019)
^^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.3.1
- Added the ompgen template file to the manifest.


v1.3.0 (December 2018)
^^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.3.0
- Fixed FutureWarning issue with destag routine (thank you honnorat!)
- Fixed computational problems with updraft_helicity, and values are no longer 
  scaled by 1000.
- Removed version constraints for wrapt and setuptools.
- Fixed xarray being a hard dependency.
- Fixed unit issues with vinterp when pressure values are extracted below 
  ground. Also added support for height fields in km and pressure fields in 
  hPa. The documentation has been improved.
- Fixed the smooth2d routine so that it actually works. It never worked 
  correctly before (nor did it work in NCL). Users can now specify the 
  center weight of the kernel and the documentation has been updated to 
  describe how it works.
- Fixed the storm relative helicity algorithm so that it works in the southern
  hemisphere. The raw algorithm now requires latitude input if used 
  in the southern hemisphere, otherwise the northern hemisphere is assumed.
- Fixed an issue with the latest version of cartopy 0.17 (thanks honnorat!)
- Fixed an issue where invalid keyword arguments weren't throwing errors when 
  extracting standard WRF variables.
- Fixed minor issues related to moving nests when using line interpolation and 
  vertical cross sections. It is still an error to request all times when 
  using lat/lon coordinates with a moving nest, but otherwise knows how to 
  run when all times are requested. This never really worked quite right.
- Removed the pyf file from setup.py since it is generated via the build
  system.
- Added an autolevels parameter for the vertical cross section so that users 
  can specify the number of vertical levels to use if they don't want to 
  specify them manually.
- The interplevel routine has been improved. Users can now specify a single 
  level, multiple levels, or a 2D array (e.g. PBLH) to interpolate to. 
  Performance has been improved when interpolating a multiple product 
  field like wspd_wdir.
- Products that produce multiple outputs can now have the outputs requested 
  individually. See :ref:`subdiagnostic-table` for a list of what is available.
- Much of this version of wrf-python has been back ported to NCL in the 
  upcoming 6.6.0 release. The diagnostics should produce the same results 
  in both packages.
- Now released under the Apache 2.0 license.



v1.2.0 (May 2018)
^^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.2.0
- Previous versions of wrf-python promoted the strings used in xarray (e.g. 
  name, attributes) to Unicode strings for Python 2.7. This caused problems 
  when porting examples for PyNGL to use wrf-python in Python 3.x. All strings 
  are now the native string type for the Python version being used. While this 
  change should be transparent to most users, any users that worked with the 
  xarray name or attribute values on Python 2.7 may run in to string related 
  errors, so we've decided to bump the major version number. 


v1.1.3 (March 2018)
^^^^^^^^^^^^^^^^^^^^^^^^^

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
  

v1.1.2 (February 2018)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.1.2
- Fix OpenMP directive issue with cloud top temperature.


v1.1.1 (February 2018)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.1.1
- Added script for building on Cheyenne with maxed out Intel settings, which 
  also required a patch for numpy.distutils.
- Fixed a few unicode characters hiding in a docstring that were causing 
  problems on Cheyenne, and also building the docs with Sphinx on Python 2.x.
- Fix issue with np.amax not working with xarray on Cheyenne, causing an error
  with the mdbz product.
- Fix cape_2d private variable bug when running with multiple CPUs.


v1.1.0 (January 2018)
^^^^^^^^^^^^^^^^^^^^^^^^^

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


v1.0.5 (September 2017)
^^^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.0.5
- Reduced the CI test file sizes by half.  


v1.0.4 (September 2017)
^^^^^^^^^^^^^^^^^^^^^^^^

- Release 1.0.4
- Fix warnings with CI tests which were caused by fill values being written 
  as NaN to the NetCDF result file.
- Added the __eq__ operator to the WrfProj projection base class.
- Fixed array order issue when using the raw CAPE routine with 1D arrays.


v1.0.3 (June 2017)
^^^^^^^^^^^^^^^^^^^^^

- Relase 1.0.3
- Fixed an issue with the cartopy Mercator subclass where the xlimits were 
  being calculated to the same value (or very close), causing blank plots.


v1.0.2 (May 2017)
^^^^^^^^^^^^^^^^^^^^^

- Release 1.0.2
- Fixed issue with the wspd_wdir product types when sequences of files are 
  used.


v1.0.1 (March 2017)
^^^^^^^^^^^^^^^^^^^^^

- Release 1.0.1
- Fixed issue with initialization of PolarStereographic and LatLon map 
  projection objects.
- Fixed issue where XTIME could be included in the coordinate list of a 
  variable, but the actual XTIME variable could be missing.  NCL allows this,
  so wrf-python should as well.
  

v1.0.0 (March 2017)
^^^^^^^^^^^^^^^^^^^^^

- Release 1.0.0.
- Fixed issue with not being able to set the thread-local coordinate cache to 
  0 to disable it.  Also, the cache will now correctly resize itself when 
  the size is reduced to less than its current setting.
- Fixed an issue with the '0000-00-00 00:00:00' time used in geo_em files 
  causing crashes due to the invalid time.  The time is now set to 
  numpy.datetime64('NaT').
- Fixed issue with wrf.cape_3d not working correctly with a single 
  column of data.


  


