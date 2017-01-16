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

    from __future__ import print_function

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

    from __future__ import print_function

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

    from __future__ import print_function

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
wrf-python provides the :meth:`wrf.to_np` function for this purpose. Although
an :class:`xarray.DataArary` object already contains the 
:attr:`xarray.DataArray.values` attribute to extract the Numpy array, there is a 
problem when working with compiled extensions. The behavior for xarray (and pandas) 
is to convert missing/fill values to NaN, which may cause crashes when working
with compiled extensions.  Also, some existing code may be designed to work with 
:class:`numpy.ma.MaskedArray`, and numpy arrays with NaN may not work with it.

The :meth:`wrf.to_np` function does the following:

#. If no missing/fill values are used, :meth:`wrf.to_np` simply returns the 
   :attr:`xarray.DataArray.values` attribute.

#. If missing/fill values are used, then :meth:`wrf.to_np` replaces the NaN
   values with the _FillValue found in the :attr:`xarray.DataArray.attrs` 
   attribute (required) and a :class:`numpy.ma.MaskedArray` is returned.

.. code-block:: python

    from __future__ import print_function

    from netCDF4 import Dataset
    from wrf import getvar
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the 3D CAPE, which contains missing values
    cape_3d = getvar(ncfile, "cape_3d")
    
    # Since there are missing values, this should return a MaskedArray
    cape_3d_ndarray = to_np(cape_3d)
    
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

    from __future__ import print_function
    
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

    from __future__ import print_function

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

    from __future__ import print_function

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
    
    from __future__ import print_function

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
    
    from __future__ import print_function
    
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

    from __future__ import print_function, division

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

    from __future__ import print_function, division

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

    from __future__ import print_function, division

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

    from __future__ import print_function, division

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

    from __future__ import print_function, division

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

    from __future__ import print_function, division

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
    t2_line = interpline(t2, pivot_point=pivot_point, angle=angle, latlon=True)
    
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

    from __future__ import print_function, division

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

    from __future__ import print_function

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

    from __future__ import print_function

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

    from __future__ import print_function
    
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


Mapping Helper Routines
-------------------------

wrf-python includes several routines to assist with plotting, primarily for 
obtaining the mapping object used for cartopy, basemap, and PyNGL.  For all 
three plotting systems, the mapping object can be determined directly from 
a variable when using xarray, or can be obtained from the WRF output file(s) 
if xarray is turned off.  

Also included are utilities for extracting the geographic boundaries 
directly from xarray variables.  This can be useful in situations where you 
only want to work with subsets (slices) of a large domain, but don't want to 
define the map projection over the subset region.


Cartopy Example Using a Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we're going to extract the cartopy mapping object from a
diagnostic variable (slp), the lat,lon coordinates, and the geographic 
boundaries.  Next, we're going to take a subset of the diagnostic variable 
and extract the geographic boundaries.  Some of the variables 
will be printed for demonstration.

.. code-block:: python

    from __future__ import print_function
    
    from netCDF4 import Dataset
    from wrf import getvar, get_cartopy, latlon_coords, geo_bounds

    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Use SLP for the example variable
    slp = getvar(ncfile, "slp")
    
    # Get the cartopy mapping object
    cart_proj = get_cartopy(slp)
    
    print (cart_proj)
    
    # Get the latitude and longitude coordinate.  This is usually needed for plotting.
    lats, lons = latlon_coords(slp)
    
    # Get the geobounds for the SLP variable
    bounds = geo_bounds(slp)
    
    print (bounds)
    
    # Get the geographic boundaries for a subset of the domain
    slp_subset = slp[150:250, 150:250]
    slp_subset_bounds = geo_bounds(slp_subset)
    
    print (slp_subset_bounds)


Result:

.. code-block:: none

    <cartopy.crs.LambertConformal object at 0x115374290>
    GeoBounds(CoordPair(lat=25.9246292114, lon=-119.675048828), CoordPair(lat=29.0761833191, lon=-117.46484375))
    GeoBounds(CoordPair(lat=25.9246292114, lon=-119.675048828), CoordPair(lat=29.0761833191, lon=-117.46484375))


Cartopy Example Using WRF Output Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, the cartopy mapping object and geographic boundaries 
will be extracted directly from the netcdf variable.

.. code-block:: python

    from __future__ import print_function
    
    from netCDF4 import Dataset
    from wrf import get_cartopy, geo_bounds
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the cartopy mapping object from the netcdf file
    cart_proj = get_cartopy(wrfin=ncfile)
    
    print (cart_proj)
    
    # Get the geobounds from the netcdf file (by default, uses XLAT, XLONG)
    # You can supply a variable name to get the staggered boundaries
    bounds = geo_bounds(wrfin=ncfile)
    
    print (bounds)
    
Result:

.. code-block:: none

    <cartopy.crs.LambertConformal object at 0x11d3be650>
    GeoBounds(CoordPair(lat=21.1381225586, lon=-122.719528198), CoordPair(lat=47.8436355591, lon=-60.9013671875))
    

Basemap Example Using a Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we're going to extract the basemap mapping object from a
diagnostic variable (slp), the lat,lon coordinates, and the geographic 
boundaries.  Next, we're going to take a subset of the diagnostic variable 
and extract the geographic boundaries.  Some of the variables will be 
printed for demonstration.

.. code-block:: python

    from __future__ import print_function
    
    from netCDF4 import Dataset
    from wrf import getvar, get_basemap, latlon_coords, geo_bounds

    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    slp = getvar(ncfile, "slp")
    
    # Get the basemap mapping object
    bm = get_basemap(slp)
    
    print (bm)
    
    # Get the latitude and longitude coordinate.  This is usually needed for plotting.
    lats, lons = latlon_coords(slp)
    
    # Get the geobounds for the SLP variable
    bounds = geo_bounds(slp)
    
    print(bounds)
    
    # Get the geographic boundaries for a subset of the domain
    slp_subset = slp[150:250, 150:250]
    slp_subset_bounds = geo_bounds(slp_subset)
    
    print (slp_subset_bounds)

Result:

.. code-block:: none

    <mpl_toolkits.basemap.Basemap object at 0x114d65650>
    GeoBounds(CoordPair(lat=21.1381225586, lon=-122.719528198), CoordPair(lat=47.8436355591, lon=-60.9013671875))
    GeoBounds(CoordPair(lat=25.9246292114, lon=-119.675048828), CoordPair(lat=29.0761833191, lon=-117.46484375)


Basemap Example Using WRF Output Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, the basemap mapping object and geographic boundaries 
will be extracted directly from the netcdf variable.

.. code-block:: python

    from __future__ import print_function
    
    from netCDF4 import Dataset
    from wrf import get_basemap, geo_bounds
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the basemap object from the netcdf file
    bm = get_basemap(wrfin=ncfile)
    
    print (bm)
    
    # Get the geographic boundaries from the netcdf file
    bounds = geo_bounds(wrfin=ncfile)
    
    print (bounds)
    
Result:

.. code-block:: none

    <mpl_toolkits.basemap.Basemap object at 0x125bb4750>
    GeoBounds(CoordPair(lat=21.1381225586, lon=-122.719528198), CoordPair(lat=47.8436355591, lon=-60.9013671875))
  
    
PyNGL Example Using a Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we're going to extract the basemap mapping object from a
diagnostic variable (slp), the lat,lon coordinates, and the geographic 
boundaries.  Next, we're going to take a subset of the diagnostic variable 
and extract the geographic boundaries.  Some of the variables will be 
printed for demonstration.

.. code-block:: python

    from __future__ import print_function
    
    from netCDF4 import Dataset
    from wrf import getvar, get_pyngl, latlon_coords, geo_bounds

    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Use SLP as the example variable
    slp = getvar(ncfile, "slp")
    
    # Get the pyngl resources from the variable
    pyngl_resources = get_pyngl(slp)
    
    print (pyngl_resources)
    
    # Get the latitude and longitude coordinate.  This is needed for plotting.
    lats, lons = latlon_coords(slp)
    
    # Get the geobounds from the SLP variable
    bounds = geo_bounds(slp)
    
    print(bounds)
    
    # Get the geographic boundaries for a subset of the domain
    slp_subset = slp[150:250, 150:250]
    slp_subset_bounds = geo_bounds(slp_subset)
    
    print (slp_subset_bounds)

Result:

.. code-block:: none

    <Ngl.Resources instance at 0x114cabbd8>
    GeoBounds(CoordPair(lat=21.1381225586, lon=-122.719528198), CoordPair(lat=47.8436355591, lon=-60.9013671875))
    GeoBounds(CoordPair(lat=25.9246292114, lon=-119.675048828), CoordPair(lat=29.0761833191, lon=-117.46484375))


PyNGL Example Using WRF Output Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, the basemap mapping object and geographic boundaries 
will be extracted directly from the netcdf variable.

.. code-block:: python

    from __future__ import print_function
    
    from netCDF4 import Dataset
    from wrf import get_pyngl, geo_bounds
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the pyngl resources from the netcdf file
    pyngl_resources = get_pyngl(wrfin=ncfile)
    
    print (pyngl_resources)
    
    # Get the geographic boundaries from the netcdf file
    bounds = geo_bounds(wrfin=ncfile)
    
    print (bounds)
    
Result:

.. code-block:: none

    <Ngl.Resources instance at 0x115391f80>
    GeoBounds(CoordPair(lat=21.1381225586, lon=-122.719528198), CoordPair(lat=47.8436355591, lon=-60.9013671875))

   
Moving Nests
^^^^^^^^^^^^^^^^^^^^

When a domain nest is moving, the domain boundaries become a function of time when 
combining the files using the 'cat' method.  When using 'join', the domain boundaries
become a function of both file and time. As a result, the methods that 
depend on geographic boundaries (:meth:`wrf.geo_bounds`, :meth:`wrf.get_basemap`, etc)
will return arrays of objects rather than a single object when multiple times 
and/or files are detected in the underlying coordinate variables.  

An exception is :meth:`wrf.get_cartopy`, which contains no geographic 
boundary information in the mapping object.  Instead, the 
:meth:`wrf.cartopy_xlim` and :meth:`wrf.cartopy_ylim` methods can be used to 
get the array of matplotlib axes boundaries (returned in the axes projection 
coordinates).

Geographic Boundaries with Moving Nest Example
***************************************************

In this example, the geographic boundaries are extracted from a sequence 
of files that use a moving nest.  The result will be an array of 
:class:`wrf.GeoBounds` objects.

.. code-block:: python

    from __future__ import print_function
    
    from glob import glob
    from netCDF4 import Dataset as nc
    
    from wrf import getvar, ALL_TIMES, geo_bounds 
    
    # Get all the domain 02 files
    wrf_filenames = glob("wrf_files/wrf_vortex_multi/wrfout_d02_*")
    ncfiles = [nc(x) for x in wrf_filenames]
    
    # SLP is the example variable and includes all times
    slp = getvar(ncfiles, "slp", timeidx=ALL_TIMES)
    
    # Get the geographic boundaries
    bounds = geo_bounds(slp)
    print (bounds)

Result:

.. code-block:: none

    [ GeoBounds(CoordPair(lat=21.3020038605, lon=-90.5740585327), CoordPair(lat=29.0274410248, lon=-82.0291671753))
     GeoBounds(CoordPair(lat=21.3020038605, lon=-90.3042221069), CoordPair(lat=29.0274410248, lon=-81.7593231201))
     GeoBounds(CoordPair(lat=21.3020038605, lon=-90.8438949585), CoordPair(lat=29.0274410248, lon=-82.2990036011))
     GeoBounds(CoordPair(lat=21.3020038605, lon=-91.1137390137), CoordPair(lat=29.0274410248, lon=-82.5688400269))
     GeoBounds(CoordPair(lat=21.8039493561, lon=-91.6534042358), CoordPair(lat=29.4982528687, lon=-83.1085205078))
     GeoBounds(CoordPair(lat=22.0542640686, lon=-92.193107605), CoordPair(lat=29.7328338623, lon=-83.6481933594))
     GeoBounds(CoordPair(lat=22.5535621643, lon=-92.7327728271), CoordPair(lat=30.2003688812, lon=-84.1878738403))
     GeoBounds(CoordPair(lat=22.8025398254, lon=-93.0026092529), CoordPair(lat=30.4333114624, lon=-84.4577102661))
     GeoBounds(CoordPair(lat=23.0510597229, lon=-93.2724456787), CoordPair(lat=30.665681839, lon=-84.7275543213))]


Cartopy Mapping with Moving Nest Example
********************************************

In this example, a cartopy mapping object is extracted from a variable
that uses a moving nest.  Since cartopy objects do not include geographic 
boundary information, only a single cartopy object is returned.  However, 
if the axes xlimits and ylimits are desired, the :meth:`wrf.cartopy_xlim` and 
:meth:`wrf.cartopy_ylim` functions can be used to obtain the array of 
moving boundaries in the axes projected coordinate space.

.. code-block:: python
    
    from __future__ import print_function
    
    from glob import glob
    from netCDF4 import Dataset as nc
    
    from wrf import getvar, ALL_TIMES, get_cartopy, cartopy_xlim, cartopy_ylim 
    
    # Get all of the domain 02 WRF output files
    wrf_filenames = glob("wrf_files/wrf_vortex_multi/wrfout_d02_*")
    ncfiles = [nc(x) for x in wrf_filenames]
    
    # Use SLP as the example variable and include all times
    slp = getvar(ncfiles, "slp", timeidx=ALL_TIMES)
    
    # Get the cartopy mapping object
    cart_proj = get_cartopy(slp)
    print (cart_proj)
    print ("\n")
    
    # Get the array of axes x-limits
    xlims = cartopy_xlim(slp)
    print (xlims)
    print ("\n")
    
    # Get the array of axes y-limits
    ylims = cartopy_ylim(slp)
    print (ylims)


Result:

.. code-block:: none

    <wrf.projection.MercatorWithLatTS object at 0x13893c9b0>
    
    [[-174999.8505754546, 774999.5806103835]
     [-145000.11853874932, 805000.1608638937]
     [-204999.58261215844, 744999.8485736783]
     [-235000.16286567, 715000.1165369744]
     [-294998.77872227144, 654999.804246759]
     [-355001.6356629085, 595000.34017335]
     [-415000.25151950994, 535000.0278831345]
     [-444999.98355621524, 505000.29584642925]
     [-474999.7155929191, 474999.7155929177]]
    
    [[2424828.507236154, 3374828.14098255]
     [2424828.507236154, 3374828.14098255]
     [2424828.507236154, 3374828.14098255]
     [2424828.507236154, 3374828.14098255]
     [2484829.1182174017, 3434828.972518358]
     [2514829.1041871575, 3464828.196283651]
     [2574829.0041584675, 3524828.8880928173]
     [2604829.1786526926, 3554829.5610342724]
     [2634828.9016262344, 3584828.016406863]]


Basemap Mapping with Moving Nest Example
*******************************************

In this example, basemap objects are extracted from a variable that uses a moving 
nest.  An array of basemap objects is returned because the 
basemap object includes geographic boundary information.  

.. code-block:: python
    
    from __future__ import print_function
    
    from glob import glob
    from netCDF4 import Dataset as nc
    
    from wrf import getvar, ALL_TIMES, get_basemap 
    
    # Get all of the domain 02 WRF output files
    wrf_filenames = glob("wrf_files/wrf_vortex_multi/wrfout_d02_*")
    ncfiles = [nc(x) for x in wrf_filenames]
    
    # Use SLP as the reference variable and include all times
    slp = getvar(ncfiles, "slp", timeidx=ALL_TIMES)
    
    # Get the array of basemap objects
    bm = get_basemap(slp)
    print (bm)
    print ("\n")
    
Result:

.. code-block:: none

    [<mpl_toolkits.basemap.Basemap object at 0x1327bc510>
     <mpl_toolkits.basemap.Basemap object at 0x115a9a790>
     <mpl_toolkits.basemap.Basemap object at 0x115a9a750>
     <mpl_toolkits.basemap.Basemap object at 0x115a9a7d0>
     <mpl_toolkits.basemap.Basemap object at 0x115a9a850>
     <mpl_toolkits.basemap.Basemap object at 0x115a9a8d0>
     <mpl_toolkits.basemap.Basemap object at 0x115a9a950>
     <mpl_toolkits.basemap.Basemap object at 0x115a9a9d0>
     <mpl_toolkits.basemap.Basemap object at 0x115a9aa50>]
    

PyNGL Mapping with Moving Nest Example
*****************************************

In this example, pyngl resource objects are extracted from a variable that uses 
a moving nest.  An array of pyngl resource objects is returned because the 
pyngl object includes geographic boundary information.

.. code-block:: python
    
    from __future__ import print_function
    
    from glob import glob
    from netCDF4 import Dataset as nc
    
    from wrf import getvar, ALL_TIMES, get_pyngl 
    
    # Get the domain 02 WRF output files
    wrf_filenames = glob("wrf_files/wrf_vortex_multi/wrfout_d02_*")
    ncfiles = [nc(x) for x in wrf_filenames]
    
    # Use SLP as the example variable and include all times
    slp = getvar(ncfiles, "slp", timeidx=ALL_TIMES)
    
    # Get the array of pyngl resource objects
    bm = get_pyngl(slp)
    print (bm)
    print ("\n")
    
Result:

.. code-block:: none

    [<Ngl.Resources instance at 0x140cd30e0>
     <Ngl.Resources instance at 0x11d3187a0>
     <Ngl.Resources instance at 0x11d3185a8>
     <Ngl.Resources instance at 0x11d3188c0>
     <Ngl.Resources instance at 0x11d318878>
     <Ngl.Resources instance at 0x11d3183f8>
     <Ngl.Resources instance at 0x11d318950>
     <Ngl.Resources instance at 0x11d318a70>
     <Ngl.Resources instance at 0x11d318710>]
    
