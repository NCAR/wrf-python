How To Use
============

Basic Usage
----------------

.. _diagnostic-usage:

Computing Diagnostic Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The primary use for the :meth:`wrf.getvar` function is to return diagnostic 
variables that require a calculation, since WRF does not produce these variables
natively. These diagnostics include CAPE, storm relative helicity, 
omega, sea level pressure, etc. A table of all available diagnostics can be 
found here: :ref:`diagnostic-table`.

In the example below, sea level pressure is calculated and printed.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the Sea Level Pressure
    slp = getvar(ncfile, "slp")
    
    print(slp)
    
Result: 

.. code-block:: none

    <xarray.DataArray u'slp' (south_north: 1059, west_east: 1799)>
    array([[ 1012.22033691,  1012.29815674,  1012.24786377, ...,
         1010.13201904,  1009.93231201,  1010.06707764],
       [ 1012.43286133,  1012.44476318,  1012.33666992, ...,
         1010.1072998 ,  1010.10845947,  1010.04760742],
       [ 1012.39544678,  1012.38085938,  1012.41705322, ...,
         1010.22937012,  1010.05596924,  1010.02679443],
       ..., 
        [ 1009.0423584 ,  1009.06921387,  1008.98779297, ...,
         1019.19281006,  1019.14434814,  1019.1105957 ],
       [ 1009.22485352,  1009.07513428,  1008.98638916, ...,
         1019.07189941,  1019.04266357,  1019.0612793 ],
       [ 1009.18896484,  1009.1071167 ,  1008.97979736, ...,
         1018.91778564,  1018.95684814,  1019.04748535]], dtype=float32) 
    Coordinates:
        XLONG        (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT         (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
        Time         datetime64[ns] 2016-10-07
      * south_north  (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east    (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
    Attributes:
        FieldType: 104
        MemoryOrder: XY
        description: sea level pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
    
.. _extract_ncvars:

Extracting WRF NetCDF Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to computing diagnostic variables (see :ref:`diagnostic-usage`), 
the :meth:`wrf.getvar` function can be used to extract regular WRF-ARW output 
NetCDF variables.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    p = getvar(ncfile, "P")
    
    print(p)

Result:

.. code-block:: none    

    <xarray.DataArray u'P' (bottom_top: 50, south_north: 1059, west_east: 1799)>
    array([[[  1.21753906e+03,   1.22532031e+03,   1.22030469e+03, ...,
           1.00760156e+03,   9.87640625e+02,   1.00111719e+03],
        [  1.23877344e+03,   1.24004688e+03,   1.22926562e+03, ...,
           1.00519531e+03,   1.00529688e+03,   9.99171875e+02],
        [  1.23503906e+03,   1.23367188e+03,   1.23731250e+03, ...,
           1.01739844e+03,   1.00005469e+03,   9.97093750e+02],
        ..., 
        [  1.77978516e+00,   1.77050781e+00,   1.79003906e+00, ...,
           4.22949219e+00,   4.25659180e+00,   4.13647461e+00],
        [  1.73291016e+00,   1.76879883e+00,   1.77978516e+00, ...,
           4.24047852e+00,   4.24707031e+00,   4.13549805e+00],
        [  1.71533203e+00,   1.65722656e+00,   1.67480469e+00, ...,
           4.06884766e+00,   4.03637695e+00,   4.04785156e+00]]], dtype=float32)
    Coordinates:
        XLONG        (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT         (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
        Time         datetime64[ns] 2016-10-07
      * bottom_top   (bottom_top) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * south_north  (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east    (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
           
                    
Disabling xarray and metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes you just want a regular numpy array and don't care about metadata.  
This is often the case when you are working with compiled extensions.  Metadata 
can be disabled in one of two ways.

#. disable xarray completely
#. set the *meta* function parameter to False.
    
The example below illustrates both.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, disable_xarray
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Disable xarray completely
    disable_xarray()
    p_no_meta = getvar(ncfile, "P")
    print (type(p_no_meta))
    enable_xarray()
    
    # Disable by using the meta parameter
    p_no_meta = getvar(ncfile, "P", meta=False)
    print (type(p_no_meta))
    
Result:

.. code-block:: none

    <type 'numpy.ndarray'>
    <type 'numpy.ndarray'>

Extracting a Numpy Array from a DataArray
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you need to convert an :class:`xarray.DataArray` to a :class:`numpy.ndarray`,
wrf-python provides the :meth:`wrf.npvalues` function for this purpose. Although
an :class:`xarray.DataArary` object already contains the 
:attr:`xarray.DataArray.values` attribute to extract the Numpy array, there is a 
problem when working with compiled extensions. The behavior for xarray (and pandas) 
is to convert missing/fill values to NaN, which may cause crashes when working
with compiled extensions.  Also, some existing code may be designed to work with 
:class:`numpy.ma.MaskedArray`, and numpy arrays with NaN may not work with it.

The :meth:`wrf.npvalues` function does the following:

#. If no missing/fill values are used, :meth:`wrf.npvalues` simply returns the 
   :attr:`xarray.DataArray.values` attribute.

#. If missing/fill values are used, then :meth:`wrf.npvalues` replaces the NaN
   values with the _FillValue found in the :attr:`xarray.DataArray.attrs` 
   attribute (required) and a :class:`numpy.ma.MaskedArray` is returned.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the Sea Level Pressure
    cape_3d = getvar(ncfile, "cape_3d")
    
    cape_3d_ndarray = npvalues(cape_3d)
    
    print(type(cape_3d_ndarray))


Result:

.. code-block:: none

    <class 'numpy.ma.core.MaskedArray'>


Sequences of Files
----------------------
    
Combining Multiple Files Using the 'cat' Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 'cat' (concatenate) method aggregates all files in the sequence along the 
'Time' dimension, which will be the leftmost dimension for the output array.  
To include all of the times, in all of the files, in the output array, set the 
*timeidx* parameter to :data:`wrf.ALL_TIMES` (an alias for None).  If a single 
value is specified for *timeidx*, then the time index is assumed to be taken from 
the concatenation of all times for all files.

It is import to note that no sorting is performed in the :meth:`wrf.getvar` 
routine, so all files in the sequence must be sorted prior to calling this 
function.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, ALL_TIMES
    
    
    # Creating a simple test list with three timesteps
    wrflist = [Dataset("wrfout_d01_2016-10-07_00_00_00"), 
               Dataset("wrfout_d01_2016-10-07_01_00_00"), 
               Dataset("wrfout_d01_2016-10-07_02_00_00")]
    
    # Extract the 'P' variable for all times          
    p_cat = getvar(wrflist, "P", timeidx=ALL_TIMES, method="cat")
    
    print(p_cat)

Result:

.. code-block:: none

    <xarray.DataArray u'P' (Time: 3, bottom_top: 50, south_north: 1059, west_east: 1799)>
    array([[[[  1.21753906e+03,   1.22532031e+03,   1.22030469e+03, ...,
            1.00760156e+03,   9.87640625e+02,   1.00111719e+03],
         [  1.23877344e+03,   1.24004688e+03,   1.22926562e+03, ...,
            1.00519531e+03,   1.00529688e+03,   9.99171875e+02],
         [  1.23503906e+03,   1.23367188e+03,   1.23731250e+03, ...,
            1.01739844e+03,   1.00005469e+03,   9.97093750e+02],
         ..., 
         [  1.77978516e+00,   1.77050781e+00,   1.79003906e+00, ...,
            4.22949219e+00,   4.25659180e+00,   4.13647461e+00],
         [  1.73291016e+00,   1.76879883e+00,   1.77978516e+00, ...,
            4.24047852e+00,   4.24707031e+00,   4.13549805e+00],
         [  1.71533203e+00,   1.65722656e+00,   1.67480469e+00, ...,
            4.06884766e+00,   4.03637695e+00,   4.04785156e+00]]]], dtype=float32)
    Coordinates:
        XLONG        (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT         (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
      * Time         (Time) datetime64[ns] 2016-10-07 2016-10-07 2016-10-07
      * bottom_top   (bottom_top) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * south_north  (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east    (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
        datetime     (Time) datetime64[ns] 2016-10-07T00:00:00 ...
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        

Combining Multiple Files Using the 'join' Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 'join' method combines a sequence of files by adding a new leftmost 
dimension for the file/sequence index. In situations where there are multiple 
files with multiple times, and the last file contains less times than the 
previous files, the remaining arrays will be arrays filled with missing values.  
There are checks in place within the wrf-python algorithms to look for these missing
arrays, but be careful when calling compiled routines outside of wrf-python.
    
In most cases, *timeidx* parameter should be set to :data:`wrf.ALL_TIMES`.  If 
a *timeidx* value is specified, then this time index is used when extracting the 
variable from each file.  In cases where there are multiple files with multiple 
time steps, this is probably nonsensical, since the nth time index for each 
file represents a different time.
    
In general, join is rarely used, so the concatenate method should be used 
for most cases. 

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, ALL_TIMES
    
    
    # Creating a simple test list with three timesteps
    wrflist = [Dataset("wrfout_d01_2016-10-07_00_00_00"), 
               Dataset("wrfout_d01_2016-10-07_01_00_00"), 
               Dataset("wrfout_d01_2016-10-07_02_00_00")]
    
    # Extract the 'P' variable for all times          
    p_join = getvar(wrflist, "P", timeidx=ALL_TIMES, method="join")
    
    print(p_join)
    
Result:

.. code-block:: none

    <xarray.DataArray u'P' (file: 3, bottom_top: 50, south_north: 1059, west_east: 1799)>
    array([[[[  1.21753906e+03,   1.22532031e+03,   1.22030469e+03, ...,
            1.00760156e+03,   9.87640625e+02,   1.00111719e+03],
         [  1.23877344e+03,   1.24004688e+03,   1.22926562e+03, ...,
            1.00519531e+03,   1.00529688e+03,   9.99171875e+02],
         [  1.23503906e+03,   1.23367188e+03,   1.23731250e+03, ...,
            1.01739844e+03,   1.00005469e+03,   9.97093750e+02],
         ..., 
         [  1.77978516e+00,   1.77050781e+00,   1.79003906e+00, ...,
            4.22949219e+00,   4.25659180e+00,   4.13647461e+00],
         [  1.73291016e+00,   1.76879883e+00,   1.77978516e+00, ...,
            4.24047852e+00,   4.24707031e+00,   4.13549805e+00],
         [  1.71533203e+00,   1.65722656e+00,   1.67480469e+00, ...,
            4.06884766e+00,   4.03637695e+00,   4.04785156e+00]]]], dtype=float32)
    Coordinates:
        XLONG        (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT         (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
      * bottom_top   (bottom_top) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * south_north  (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east    (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * file         (file) int64 0 1 2
        datetime     (file) datetime64[ns] 2016-10-07T00:00:00 ...
        Time         int64 0
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
    
                        
Note how the 'Time' dimension was replaced with the 'file' dimension, due to the 
numpy's automatic squeezing of the single 'Time' dimension. To maintain the 
'Time' dimension, set the *squeeze* parameter to False.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, ALL_TIMES
    
    
    # Creating a simple test list with three timesteps
    wrflist = [Dataset("wrfout_d01_2016-10-07_00_00_00"), 
               Dataset("wrfout_d01_2016-10-07_01_00_00"), 
               Dataset("wrfout_d01_2016-10-07_02_00_00")]
    
    # Extract the 'P' variable for all times          
    p_join = getvar(wrflist, "P", timeidx=ALL_TIMES, method="join", squeeze=False)
    
    print(p_join)
    
Result

.. code-block:: none

    <xarray.DataArray u'P' (file: 3, Time: 1, bottom_top: 50, south_north: 1059, west_east: 1799)>
    array([[[[[  1.21753906e+03,   1.22532031e+03,   1.22030469e+03, ...,
             1.00760156e+03,   9.87640625e+02,   1.00111719e+03],
          [  1.23877344e+03,   1.24004688e+03,   1.22926562e+03, ...,
             1.00519531e+03,   1.00529688e+03,   9.99171875e+02],
          [  1.23503906e+03,   1.23367188e+03,   1.23731250e+03, ...,
             1.01739844e+03,   1.00005469e+03,   9.97093750e+02],
          ..., 
          [  1.77978516e+00,   1.77050781e+00,   1.79003906e+00, ...,
             4.22949219e+00,   4.25659180e+00,   4.13647461e+00],
          [  1.73291016e+00,   1.76879883e+00,   1.77978516e+00, ...,
             4.24047852e+00,   4.24707031e+00,   4.13549805e+00],
          [  1.71533203e+00,   1.65722656e+00,   1.67480469e+00, ...,
             4.06884766e+00,   4.03637695e+00,   4.04785156e+00]]]]], dtype=float32)
    Coordinates:
        XLONG        (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT         (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
      * bottom_top   (bottom_top) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * south_north  (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east    (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * file         (file) int64 0 1 2
        datetime     (file, Time) datetime64[ns] 2016-10-07T00:00:00 ...
      * Time         (Time) int64 0
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)

                   
Dictionaries of WRF File Sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dictionaries can also be used as input to the :meth:`wrf.getvar` functions.
This can be useful when working with ensembles.  However, all WRF files in the 
dictionary must have the same dimensions.  The result is an array where the 
leftmost dimension is the keys from the dictionary.  Nested dictionaries 
are allowed.

The *method* argument is used to describe how each sequence in the dictionary 
will be combined.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, ALL_TIMES
    
    wrf_dict = {"ens1" : [Dataset("ens1/wrfout_d01_2016-10-07_00_00_00"), 
                          Dataset("ens1/wrfout_d01_2016-10-07_01_00_00"), 
                          Dataset("ens1/wrfout_d01_2016-10-07_02_00_00")],
                "ens2" : [Dataset("ens2/wrfout_d01_2016-10-07_00_00_00"), 
                          Dataset("ens2/wrfout_d01_2016-10-07_01_00_00"), 
                          Dataset("ens2/wrfout_d01_2016-10-07_02_00_00")]
                }
    
    p = getvar(wrf_dict, "P", timeidx=ALL_TIMES)
    
    print(p)
    
Result:

.. code-block:: none

    <xarray.DataArray 'P' (key_0: 2, Time: 2, bottom_top: 50, south_north: 1059, west_east: 1799)>
    array([[[[[  1.21753906e+03,   1.22532031e+03,   1.22030469e+03, ...,
             1.00760156e+03,   9.87640625e+02,   1.00111719e+03],
          [  1.23877344e+03,   1.24004688e+03,   1.22926562e+03, ...,
             1.00519531e+03,   1.00529688e+03,   9.99171875e+02],
          [  1.23503906e+03,   1.23367188e+03,   1.23731250e+03, ...,
             1.01739844e+03,   1.00005469e+03,   9.97093750e+02],
          ..., 
          [  1.77978516e+00,   1.77050781e+00,   1.79003906e+00, ...,
             4.22949219e+00,   4.25659180e+00,   4.13647461e+00],
          [  1.73291016e+00,   1.76879883e+00,   1.77978516e+00, ...,
             4.24047852e+00,   4.24707031e+00,   4.13549805e+00],
          [  1.71533203e+00,   1.65722656e+00,   1.67480469e+00, ...,
             4.06884766e+00,   4.03637695e+00,   4.04785156e+00]]]]], dtype=float32)
    Coordinates:
        XLONG        (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT         (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
      * Time         (Time) datetime64[ns] 2016-10-07T00:00:00 ...
      * bottom_top   (bottom_top) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * south_north  (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east    (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
        datetime     (Time) datetime64[ns] 2016-10-07T00:00:00 ...
      * key_0        (key_0) <U6 u'ens1' u'ens2'
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
                        
Interpolation Routines
--------------------------

Interpolating to a Horizontal Level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`wrf.interplevel` function is used to interpolate a 3D field to 
a specific horizontal level, usually pressure or height.

.. code-block:: python
    
    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, interplevel
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Extract the Geopotential Height and Pressure (hPa) fields
    z = getvar(ncfile, "z")
    p = getvar(ncfile, "pressure")
    
    # Compute the 500 MB Geopotential Height
    ht_500mb = interplevel(z, p, 500.)
    
    print(ht_500mb)

Result:

.. code-block:: none

    <xarray.DataArray u'height_500_hPa' (south_north: 1059, west_east: 1799)>
    array([[ 5882.16992188,  5881.87939453,  5881.81005859, ...,
         5890.14501953,  5890.23583984,  5890.33349609],
       [ 5882.71777344,  5882.17529297,  5882.1171875 , ...,
         5890.37695312,  5890.38525391,  5890.27978516],
       [ 5883.32177734,  5882.47119141,  5882.34130859, ...,
         5890.48339844,  5890.42871094,  5890.17724609],
       ..., 
       [ 5581.45800781,  5580.46826172,  5579.32617188, ...,
         5788.93554688,  5788.70507812,  5788.64453125],
       [ 5580.32714844,  5579.51611328,  5578.34863281, ...,
         5788.15869141,  5787.87304688,  5787.65527344],
       [ 5579.64404297,  5578.30957031,  5576.98632812, ...,
         5787.19384766,  5787.10888672,  5787.06933594]], dtype=float32)
    Coordinates:
        XLONG        (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT         (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
        Time         datetime64[ns] 2016-10-07
      * south_north  (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east    (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
    Attributes:
        FieldType: 104
        units: m
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        level: 500 hPa
        missing_value: 9.96920996839e+36
        _FillValue: 9.96920996839e+36


.. _vert_cross_interp:

Vertical Cross Sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`wrf.vertcross` function is used to create vertical cross sections.  
To define a cross section, a start point and an end point needs to be specified.  
Alternatively, a pivot point and an angle may be used.  The start point, 
end point, and pivot point are specified using a :class:`wrf.CoordPair` object,
and coordinates can either be in grid (x,y) coordinates or (latitude,longitude) 
coordinates. When using (latitude,longitude) coordinates, a NetCDF file object or 
a :class:`wrf.WrfProj` object must be provided.

The vertical levels can also be specified using the *levels* parameter.  If 
not specified, then approximately 100 levels will be chosen in 1% increments.

Example Using Start Point and End Point
*****************************************

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, vertcross, CoordPair
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the geopotential height (m) and pressure (hPa).
    z = getvar(ncfile, "z")
    p = getvar(ncfile, "pressure")
    
    # Define a start point and end point in grid coordinates
    start_point = CoordPair(x=0, y=(z.shape[-2]-1)//2)
    end_point = CoordPair(x=-1, y=(z.shape[-2]-1)//2)
    
    # Calculate the vertical cross section.  By setting latlon to True, this 
    # also calculates the latitude and longitude coordinates along the cross 
    # section line and adds them to the 'xy_loc' metadata to help with plotting.
    p_vert = vertcross(p, z, start_point=start_point, end_point=end_point, latlon=True)
    
    print(p_vert)
    
Result:

.. code-block:: none

    <xarray.DataArray u'pressure_cross' (vertical: 100, idx: 1798)>
    array([[          nan,           nan,           nan, ...,           nan,
                  nan,           nan],
       [ 989.66168213,  989.66802979,  989.66351318, ...,  988.05737305,
         987.99151611,  987.96917725],
       [ 959.49450684,  959.50109863,  959.50030518, ...,  958.96948242,
         958.92980957,  958.89294434],
       ..., 
       [  24.28092003,   24.27359581,   24.27034378, ...,   24.24800491,
          24.2486496 ,   24.24947357],
       [  23.2868309 ,   23.27933884,   23.27607918, ...,   23.25231361,
          23.2530098 ,   23.25384521],
       [          nan,           nan,           nan, ...,           nan,
                  nan,           nan]], dtype=float32)
    Coordinates:
        Time      datetime64[ns] 2016-10-07
        xy_loc    (idx) object CoordPair(x=0.0, y=529.0, lat=34.5279502869, lon=-127.398925781) ...
      * vertical  (vertical) float32 0.0 261.828 523.656 785.484 1047.31 1309.14 ...
      * idx       (idx) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 ...
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        orientation: (0.0, 529.0) to (1797.0, 529.0)
        missing_value: 9.96920996839e+36
        _FillValue: 9.96920996839e+36
    
    
Example Using Pivot Point and Angle
*************************************

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, vertcross, CoordPair  
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    # Get the geopotential height (m) and pressure (hPa).
    z = getvar(ncfile, "z")
    p = getvar(ncfile, "pressure")
    
    # Define a pivot point and angle in grid coordinates, with the 
    # pivot point being the center of the grid.
    pivot_point = CoordPair(x=(z.shape[-1]-1)//2, y=(z.shape[-2]-1)//2) 
    angle = 90.0
    
    # Calculate the vertical cross section.  By setting latlon to True, this 
    # also calculates the latitude and longitude coordinates along the line
    # and adds them to the metadata to help with plotting labels.
    p_vert = vertcross(p, z, pivot_point=pivot_point, angle=angle, latlon=True)
    
    print (p_vert)
    
Result:

.. code-block:: none

    <xarray.DataArray u'pressure_cross' (vertical: 100, idx: 1798)>
    array([[          nan,           nan,           nan, ...,           nan,
                  nan,           nan],
       [ 989.66168213,  989.66802979,  989.66351318, ...,  988.05737305,
         987.99151611,  987.96917725],
       [ 959.49450684,  959.50109863,  959.50030518, ...,  958.96948242,
         958.92980957,  958.89294434],
       ..., 
       [  24.28092003,   24.27359581,   24.27034378, ...,   24.24800491,
          24.2486496 ,   24.24947357],
       [  23.2868309 ,   23.27933884,   23.27607918, ...,   23.25231361,
          23.2530098 ,   23.25384521],
       [          nan,           nan,           nan, ...,           nan,
                  nan,           nan]], dtype=float32)
    Coordinates:
        Time      datetime64[ns] 2016-10-07
        xy_loc    (idx) object CoordPair(x=0.0, y=529.0, lat=34.5279502869, lon=-127.398925781) ...
      * vertical  (vertical) float32 0.0 261.828 523.656 785.484 1047.31 1309.14 ...
      * idx       (idx) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 ...
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        orientation: (0.0, 529.0) to (1797.0, 529.0) ; center=CoordPair(x=899.0, y=529.0) ; angle=90.0
        missing_value: 9.96920996839e+36
        _FillValue: 9.96920996839e+36

    
Example Using Lat/Lon Coordinates
*************************************

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, vertcross, CoordPair  
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    # Get the geopotential height (m) and pressure (hPa).
    z = getvar(ncfile, "z")
    p = getvar(ncfile, "pressure")
    lats = getvar(ncfile, "lat")
    lons = getvar(ncfile, "lon")
    
    # Making the same horizontal line, but with lats/lons
    start_lat = lats[(lats.shape[-2]-1)//2, 0]
    end_lat = lats[(lats.shape[-2]-1)//2, -1]
    start_lon = lons[(lats.shape[-2]-1)//2, 0]
    end_lon = lons[(lats.shape[-2]-1)//2, -1]
    
    # Cross section line using start_point and end_point. 
    start_point = CoordPair(lat=start_lat, lon=start_lon)
    end_point = CoordPair(lat=end_lat, lon=end_lon)
    
    # When using lat/lon coordinates, you must supply a netcdf file object, or a 
    # projection object.
    p_vert = vertcross(p, z, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True)
    print(p_vert)
    
Result:

.. code-block:: none

    <xarray.DataArray u'pressure_cross' (vertical: 100, idx: 1798)>
    array([[          nan,           nan,           nan, ...,           nan,
                  nan,           nan],
       [ 989.66168213,  989.66802979,  989.66351318, ...,  988.05737305,
         987.99151611,  987.96917725],
       [ 959.49450684,  959.50109863,  959.50030518, ...,  958.96948242,
         958.92980957,  958.89294434],
       ..., 
       [  24.28092003,   24.27359581,   24.27034378, ...,   24.24800491,
          24.2486496 ,   24.24947357],
       [  23.2868309 ,   23.27933884,   23.27607918, ...,   23.25231361,
          23.2530098 ,   23.25384521],
       [          nan,           nan,           nan, ...,           nan,
                  nan,           nan]], dtype=float32)
    Coordinates:
        Time      datetime64[ns] 2016-10-07
        xy_loc    (idx) object CoordPair(x=0.0, y=529.0, lat=34.5279502869, lon=-127.398925781) ...
      * vertical  (vertical) float32 0.0 261.828 523.656 785.484 1047.31 1309.14 ...
      * idx       (idx) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 ...
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        orientation: (0.0, 529.0) to (1797.0, 529.0)
        missing_value: 9.96920996839e+36
        _FillValue: 9.96920996839e+36


Example Using Specified Vertical Levels
*****************************************

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, vertcross, CoordPair  
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    # Get the geopotential height (m) and pressure (hPa).
    z = getvar(ncfile, "z")
    p = getvar(ncfile, "pressure")
    lats = getvar(ncfile, "lat")
    lons = getvar(ncfile, "lon")
    
    # Making the same horizontal line, but with lats/lons
    start_lat = lats[(lats.shape[-2]-1)//2, 0]
    end_lat = lats[(lats.shape[-2]-1)//2, -1]
    start_lon = lons[(lats.shape[-2]-1)//2, 0]
    end_lon = lons[(lats.shape[-2]-1)//2, -1]
    
    # Pressure using start_point and end_point.  These were obtained using 
    start_point = CoordPair(lat=start_lat, lon=start_lon)
    end_point = CoordPair(lat=end_lat, lon=end_lon)
    
    # Specify vertical levels
    levels = [1000., 2000., 3000.]
    
    # Calculate the cross section
    p_vert = vertcross(p, z, wrfin=ncfile, levels=levels, start_point=start_point, end_point=end_point, latlon=True)
    
    print(p_vert)
    
Result:

.. code-block:: none

    <xarray.DataArray u'pressure_cross' (vertical: 3, idx: 1798)>
    array([[ 906.375     ,  906.38043213,  906.39367676, ...,  907.6661377 ,
             907.63006592,  907.59191895],
           [ 804.24737549,  804.26885986,  804.28076172, ...,  806.98632812,
             806.95556641,  806.92608643],
           [ 713.24578857,  713.2722168 ,  713.27886963, ...,  716.09594727,
             716.06610107,  716.03503418]], dtype=float32)
    Coordinates:
        Time      datetime64[ns] 2016-10-07
        xy_loc    (idx) object CoordPair(x=0.0, y=529.0, lat=34.5279502869, lon=-127.398925781) ...
      * vertical  (vertical) float32 1000.0 2000.0 3000.0
      * idx       (idx) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 ...
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        orientation: (0.0, 529.0) to (1797.0, 529.0)
        missing_value: 9.96920996839e+36
        _FillValue: 9.96920996839e+36


Interpolating Two-Dimensional Fields to a Line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two-dimensional fields can be interpolated along a line, in a manner similar to 
the vertical cross section (see :ref:`vert_cross_interp`), using the 
:meth:`wrf.interpline` function. To define the line 
to interpolate along, a start point and an end point needs to be specified.  
Alternatively, a pivot point and an angle may be used.  The start point, 
end point, and pivot point are specified using a :class:`wrf.CoordPair` object,
and coordinates can either be in grid (x,y) coordinates or (latitude,longitude) 
coordinates.  When using (latitude,longitude) coordinates, a NetCDF file object or 
a :class:`wrf.WrfProj` object must also be provided.

Example Using Start Point and End Point
*****************************************

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, interpline, CoordPair  
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    # Get the 2m temperature
    t2 = getvar(ncfile, "T2")
    
    # Create a south-north line in the center of the domain using 
    # start point and end point
    start_point = CoordPair(x=(t2.shape[-1]-1)//2, y=0)
    end_point = CoordPair(x=(t2.shape[-1]-1)//2, y=-1)
    
    # Calculate the vertical cross section.  By setting latlon to True, this 
    # also calculates the latitude and longitude coordinates along the line
    # and adds them to the metadata to help with plotting labels.
    t2_line = interpline(t2, start_point=start_point, end_point=end_point, latlon=True)
    
    print(t2_line, "\n")
    
Result:

.. code-block:: none

    <xarray.DataArray u'T2_line' (line_idx: 1058)>
    array([ 302.07214355,  302.08505249,  302.08688354, ...,  279.18557739,
            279.1998291 ,  279.23132324], dtype=float32)
    Coordinates:
        Time      datetime64[ns] 2016-10-07
        xy_loc    (line_idx) object CoordPair(x=899.0, y=0.0, lat=24.3645858765, lon=-97.5) ...
      * line_idx  (line_idx) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
    Attributes:
        FieldType: 104
        description: TEMP at 2 M
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        orientation: (899.0, 0.0) to (899.0, 1057.0) 


Example Using Pivot Point and Angle
*****************************************

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, interpline, CoordPair  
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    # Get the 2m temperature
    t2 = getvar(ncfile, "T2")
    
    # Create a south-north line using pivot point and angle
    pivot_point = CoordPair((t2.shape[-1]-1)//2, (t2.shape[-2]-1)//2) 
    angle = 0.0
    
    # Calculate the vertical cross section.  By setting latlon to True, this 
    # also calculates the latitude and longitude coordinates along the line
    # and adds them to the metadata to help with plotting labels.
    t2_line = interpline(t2, start_point=start_point, end_point=end_point, latlon=True)
    
    print(t2_line, "\n")
    
Result:

.. code-block:: none

    <xarray.DataArray u'T2_line' (line_idx: 1058)>
    array([ 302.07214355,  302.08505249,  302.08688354, ...,  279.18557739,
            279.1998291 ,  279.23132324], dtype=float32)
    Coordinates:
        Time      datetime64[ns] 2016-10-07
        xy_loc    (line_idx) object CoordPair(x=899.0, y=0.0, lat=24.3645858765, lon=-97.5) ...
      * line_idx  (line_idx) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
    Attributes:
        FieldType: 104
        description: TEMP at 2 M
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        orientation: (899.0, 0.0) to (899.0, 1057.0) ; center=CoordPair(x=899, y=529) ; angle=0.0 

        
Example Using Lat/Lon Coordinates
*************************************

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, interpline, CoordPair  
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    t2 = getvar(ncfile, "T2")
    lats = getvar(ncfile, "lat")
    lons = getvar(ncfile, "lon")
    
    # Select the latitude,longitude points for a vertical line through 
    # the center of the domain.
    start_lat = lats[0, (lats.shape[-1]-1)//2]
    end_lat = lats[-1, (lats.shape[-1]-1)//2]
    start_lon = lons[0, (lons.shape[-1]-1)//2]
    end_lon = lons[-1, (lons.shape[-1]-1)//2]
    
    # Create the CoordPairs
    start_point = CoordPair(lat=start_lat, lon=start_lon)
    end_point = CoordPair(lat=end_lat, lon=end_lon)
    
    # Calculate the vertical cross section.  By setting latlon to True, this 
    # also calculates the latitude and longitude coordinates along the line
    # and adds them to the metadata to help with plotting labels.
    t2_line = interpline(t2, wrfin=ncfile, start_point=start_point, end_point=end_point, latlon=True)
    
    print (t2_line)

Result:

.. code-block:: none

    <xarray.DataArray u'T2_line' (line_idx: 1058)>
    array([ 302.07214355,  302.08505249,  302.08688354, ...,  279.18557739,
            279.1998291 ,  279.23132324], dtype=float32)
    Coordinates:
        Time      datetime64[ns] 2016-10-07
        xy_loc    (line_idx) object CoordPair(x=899.0, y=0.0, lat=24.3645858765, lon=-97.5) ...
      * line_idx  (line_idx) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
    Attributes:
        FieldType: 104
        description: TEMP at 2 M
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        orientation: (899.0, 0.0) to (899.0, 1057.0)    
            

Interpolating a 3D Field to a Surface Type 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`wrf.vinterp` is used to interpolate a field to a type of surface.  
The available surfaces are pressure, geopotential height, theta, and theta-e. 
The surface levels to interpolate also need to be specified.

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, interpline, CoordPair  
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    # Interpolate tk to theta-e levels                
    interp_levels = [200, 300, 500, 1000]
    
    interp_field = vinterp(ncfile, 
                   field=tk, 
                   vert_coord="eth", 
                   interp_levels=interp_levels, 
                   extrapolate=True, 
                   field_type="tk", 
                   log_p=True)
                    
    print(interp_field)
    
Result:

.. code-block:: none

    <xarray.DataArray u'temp' (interp_level: 4, south_north: 1059, west_east: 1799)>
    array([[[ 296.12872314,  296.1166687 ,  296.08905029, ...,  301.71026611,
              301.67956543,  301.67791748],
            [ 296.11352539,  295.95581055,  295.91555786, ...,  301.63052368,
              301.62905884,  301.65887451],
            [ 296.07556152,  295.91577148,  295.88214111, ...,  301.61499023,
              301.60287476,  301.63961792],
            ..., 
            [ 219.11134338,  219.08581543,  219.08602905, ...,  218.29879761,
              218.30923462,  218.3787384 ],
            [ 219.09260559,  219.07765198,  219.08340454, ...,  218.2855072 ,
              218.30444336,  218.37931824],
            [ 219.07936096,  219.08181763,  219.10089111, ...,  218.31173706,
              218.34288025,  218.3687439 ]]], dtype=float32)
    Coordinates:
        XLONG         (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT          (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
        Time          datetime64[ns] 2016-10-07
      * south_north   (south_north) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
      * west_east     (west_east) int64 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...
      * interp_level  (interp_level) int64 200 300 500 1000
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: temperature
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(bottom_left=(21.138123, -122.71953), 
                    top_right=(47.843636, -60.901367), stand_lon=-97.5, 
                    moad_cen_lat=38.5000038147, truelat1=38.5, truelat2=38.5, 
                    pole_lat=90.0, pole_lon=0.0)
        vert_interp_type: eth

            
Lat/Lon <-> XY Routines
--------------------------

wrf-python includes a set of routines for converting back and forth between 
latitude,longitude space and x,y space.  The methods are :meth:`wrf.xy_to_ll`,
:meth:`wrf.xy_to_ll_proj`, :meth:`wrf.ll_to_xy`, :meth:`wrf.ll_to_xy_proj`. 
The *latitude*, *longitude*, *x*, and *y* parameters to these methods 
can contain sequences if multiple points are desired to be converted.

Example With Single Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, interpline, CoordPair, xy_to_ll, ll_to_xy
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    lat_lon = xy_to_ll(ncfile, 400, 200)
    
    print(lat_lon)
    
    x_y = ll_to_xy(ncfile, lat_lon[0], lat_lon[1])
    
    print (x_y)
    
Result:

.. code-block:: none

    <xarray.DataArray u'latlon' (lat_lon: 2)>
    array([  28.55816408, -112.67827617])
    Coordinates:
      * lat_lon   (lat_lon) <U3 u'lat' u'lon'
        xy_coord  object CoordPair(x=400, y=200)
        idx       int64 0
        
        
    <xarray.DataArray u'xy' (x_y: 2)>
    array([400, 200])
    Coordinates:
        latlon_coord  object CoordPair(lat=28.5581640822, lon=-112.678276173)
      * x_y           (x_y) <U1 u'x' u'y'
        idx           int64 0
    
    
Example With Multiple Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset
    from wrf import getvar, interpline, CoordPair, xy_to_ll, ll_to_xy
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    lat_lon = xy_to_ll(ncfile, [400,105], [200,205])
    
    print(lat_lon)
    
    x_y = ll_to_xy(ncfile, lat_lon[0,:], lat_lon[1,:])
    
    print (x_y)
    
Result:

.. code-block:: none

    <xarray.DataArray u'latlon' (lat_lon: 2, idx: 2)>
    array([[  28.55816408,   27.03835783],
           [-112.67827617, -121.36392174]])
    Coordinates:
      * lat_lon   (lat_lon) <U3 u'lat' u'lon'
        xy_coord  (idx) object CoordPair(x=400, y=200) CoordPair(x=105, y=205)
      * idx       (idx) int64 0 1
        
        
    <xarray.DataArray u'xy' (x_y: 2, idx: 2)>
    array([[400, 105],
           [200, 205]])
    Coordinates:
        latlon_coord  (idx) object CoordPair(lat=28.5581640822, lon=-112.678276173) ...
      * x_y           (x_y) <U1 u'x' u'y'
      * idx           (idx) int64 0 1


