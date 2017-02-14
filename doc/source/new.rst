What's New
===========

v1.0b2
----------

- Beta release 2.
- xarray 0.9 no longer includes default index dimensions in the coordinate 
  mappings.  This was causing a crash in the routines that cause a reduction
  in dimension shape, mainly the interpolation routines.  This has been 
  fixed.
- Documentation updated to show the new output from xarray.

v1.0b1
----------

- Beta release 1.
- Added more packaging boilerplate.
- Note:  Currently unable to build with Python 3.5 on Windows, due to
  issues with distutils, numpy distutils, and mingw compiler.  Will attempt
  to find a workaround before the next release. Windows users should use 
  Python 2.7 or Python 3.4 for now.

v1.0a3
-----------

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

  
  


