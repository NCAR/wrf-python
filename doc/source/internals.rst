.. _internals:

WRF-Python Internals
========================================

WRF-Python is a collection of diagnostic and interpolation routines for 
WRF-ARW data. The API is kept to a minimal set of functions, since we've found
this to be the easiest to teach to new programmers, students, and scientists. 
Future plans do include adopting the Pangeo xarray/dask model for more 
advanced programmers, but is not currently supported as of this user guide.

A typical use case for a WRF-Python user is to:

1) Open a WRF data file (or sequence of files) using NetCDF4-python or PyNIO.
2) Compute a WRF diagnostic using :meth:`wrf.getvar`.
3) Performing other computations using methods outside of WRF-Python.
4) Creating a plot of the output using matplotlib (basemap or cartopy) or 
   PyNGL.
   
The purpose of this guide is to explain the internals of item 2 so that 
users can help contribute or support the computational diagnostics.


Overview of a :meth:`wrf.getvar` Diagnostic Computation
---------------------------------------------------------------

A diagnostic computed using the :meth:`wrf.getvar` function consists of the 
following steps:

1) Using the diagnostic string, call the appropriate 'get' function. This 
   step occurs in the :met:`wrf.getvar` routine in routines.py. 
2) Extract the required variables from the NetCDF data file (or files).
3) Compute the diagnostic using a wrapped Fortran, C, or Python routine.
4) Convert to the desired units if applicable.
5) If desired, set the metadata and return the result as an 
   :class:`xarray.DataArray`, or return a :class:`numpy.ndarray` if no 
   metadata is desired.
   
In the source directory, the :meth:`wrf.getvar` 'get' routines have a 
"g_" prefix for the naming convention (the "g" stands for "get", but didn't 
want to cause namespace conflicts with functions already named with "get" in 
the title). 

The unit conversion is handled by a wrapt decorator that can be found in 
decorators.py. The setting of the metadata is handled using a wrapt decorator, 
which can be found in the metadecorators.py file.


Overview of Compiled Computational Routines
---------------------------------------------------------

Currently, the compiled computational routines are written in Fortran 
90 and exposed the Python using f2py. The routines have been aquired over 
decades, originated from NCL's Fortran77 codebase or other tools like RIP 
(Read Interpolate Plot), and do not necessarily conform to a common 
programming mindset (e.g. some use 1D arrays, 2D arrays, etc).

The raw Fortran routines are compiled in to the :mod:`wrf._wrffortran` 
extension module, but are not particularly useful for applications in that 
raw form. These routines are imported in the extention.py module, where 
additional functionality is added to make the routines more user friendly.

The common behavior for a fully exported Fortran routine in extension.py 
is:

1) Verify that the arguments passed in are valid in shape. While f2py does this 
   as well, the errors thrown by f2py are confusing to users, so this step 
   helps provide better error messages.

2) Allocate the ouput array based on the output shape of the algorithm, 
   number of "leftmost" dimensions, and size of the data.
   
3) Iterate over the leftmost dimensions and compute output for argument 
   data slices that are of the same dimensionality as the compiled algorithm. 
   For example, if the compiled algorithm is written for two dimensional data, 
   but your data is four dimensional, you have two leftmost dimensions.
   
4) Cast the argument arrays in to the type used in the 
   compiled routine (usually for WRF data, the conversion is from 4-byte float 
   to 8-byte double).
   
5) Extract the argument arrays out of xarray in to numpy arrays 
   (if applicable) and transpose them in to Fortran ordering. Note that this 
   does not actually do any copying of the data, it simply reorders the shape 
   tuple for the data and sets the Fortran ordering flag. This allows data
   pointers from the output array to be directly passed through f2py so that 
   copying is not required from the result in to the output array.
   
The steps described above are handled in :mod:`wrapt` decorators that can be 
found in decorators.py. For some routines that produce multiple outputs or have 
atypical behavior, the special case decorators are located in specialdec.py. 



