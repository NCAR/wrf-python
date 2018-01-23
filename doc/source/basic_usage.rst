How To Use
============

Introduction
---------------

The API for wrf-python can be summarized as a variable computation/extraction
routine, several interpolation routines, and a few plotting helper utilities. 
The API is kept as simple as possible to help minimize the 
learning curve for new programmers, students, and scientists. In the future, 
we plan to extend xarray for programmers desiring a more object oriented API, 
but this remains a work in progress.

The five most commonly used routines can be summarized as:

- :meth:`wrf.getvar` - Extracts WRF-ARW NetCDF variables and 
  computes diagnostic variables that WRF does not compute (e.g. storm 
  relative helicity). This is the routine that you will use most often.
  
- :meth:`wrf.interplevel` - Interpolates a three-dimensional field to a 
  horizontal plane at a specified level using simple (fast) linear 
  interpolation (e.g. 850 hPa temperature).
  
- :meth:`wrf.vertcross` - Interpolates a three-dimensional field to a vertical 
  plane through a user-specified horizontal line (i.e. a cross section).
  
- :meth:`wrf.interpline` - Interpolates a two-dimensional field to a 
  user-specified line.
  
- :meth:`wrf.vinterp` - Interpolates a three-dimensional field to 
  user-specified  'surface' levels (e.g. theta-e levels). This is a smarter, 
  albeit slower, version of :meth:`wrf.interplevel`. 

Basic Usage
----------------

.. _diagnostic-usage:

Computing Diagnostic Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The primary use for the :meth:`wrf.getvar` function is to return diagnostic 
variables that require a calculation, since WRF does not produce these 
variables natively. These diagnostics include CAPE, storm relative helicity, 
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
    array([[ 1012.220337,  1012.298157,  1012.247864, ...,  1010.132019,
             1009.932312,  1010.067078],
           [ 1012.432861,  1012.444763,  1012.33667 , ...,  1010.1073  ,
             1010.108459,  1010.047607],
           [ 1012.395447,  1012.380859,  1012.417053, ...,  1010.22937 ,
             1010.055969,  1010.026794],
           ..., 
           [ 1009.042358,  1009.069214,  1008.987793, ...,  1019.19281 ,
             1019.144348,  1019.110596],
           [ 1009.224854,  1009.075134,  1008.986389, ...,  1019.071899,
             1019.042664,  1019.061279],
           [ 1009.188965,  1009.107117,  1008.979797, ...,  1018.917786,
             1018.956848,  1019.047485]], dtype=float32)
    Coordinates:
        XLONG    (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT     (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
        Time     datetime64[ns] 2016-10-07
    Dimensions without coordinates: south_north, west_east
    Attributes:
        FieldType: 104
        MemoryOrder: XY
        description: sea level pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)

    
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
    array([[[  1.217539e+03,   1.225320e+03, ...,   9.876406e+02,   1.001117e+03],
            [  1.238773e+03,   1.240047e+03, ...,   1.005297e+03,   9.991719e+02],
            ..., 
            [  9.208594e+02,   9.059141e+02, ...,   1.902922e+03,   1.904805e+03],
            [  9.172734e+02,   9.091094e+02, ...,   1.894375e+03,   1.903422e+03]],
    
           [[  1.219562e+03,   1.210273e+03, ...,   9.973984e+02,   9.907891e+02],
            [  1.224578e+03,   1.223508e+03, ...,   9.985547e+02,   9.921172e+02],
            ..., 
            [  9.012734e+02,   9.052031e+02, ...,   1.897766e+03,   1.894500e+03],
            [  9.137500e+02,   9.071719e+02, ...,   1.893273e+03,   1.893664e+03]],
    
           ..., 
           [[  7.233154e+00,   7.224121e+00, ...,   3.627930e+00,   3.613770e+00],
            [  7.226318e+00,   7.358154e+00, ...,   3.725098e+00,   3.634033e+00],
            ..., 
            [  5.354248e+00,   5.406006e+00, ...,   1.282715e+01,   1.264844e+01],
            [  5.295410e+00,   5.177490e+00, ...,   1.256274e+01,   1.257642e+01]],
    
           [[  2.362061e+00,   2.376221e+00, ...,   1.151367e+00,   1.156982e+00],
            [  2.342529e+00,   2.403809e+00, ...,   1.198486e+00,   1.155273e+00],
            ..., 
            [  1.732910e+00,   1.768799e+00, ...,   4.247070e+00,   4.135498e+00],
            [  1.715332e+00,   1.657227e+00, ...,   4.036377e+00,   4.047852e+00]]], dtype=float32)
    Coordinates:
        XLONG    (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT     (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
        Time     datetime64[ns] 2016-10-07
    Dimensions without coordinates: bottom_top, south_north, west_east
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
           
                    
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
        datetime     (Time) datetime64[ns] 2016-10-07T00:00:00 ...
    Dimensions without coordinates: bottom_top, south_north, west_east
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
        

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
    array([[[[  1.217539e+03, ...,   1.001117e+03],
             ..., 
             [  9.172734e+02, ...,   1.903422e+03]],
            ..., 
            [[  2.362061e+00, ...,   1.156982e+00],
             ..., 
             [  1.715332e+00, ...,   4.047852e+00]]]], dtype=float32)
    Coordinates:
        XLONG     (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT      (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
      * file      (file) int64 0 1 2
        datetime  (file) datetime64[ns] 2016-10-07 ...
    Dimensions without coordinates: bottom_top, south_north, west_east
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
    
                        
Note how the 'Time' dimension was replaced with the 'file' dimension, due to  
numpy's automatic squeezing of the single element 'Time' dimension. To maintain 
the 'Time' dimension, set the *squeeze* parameter to False.

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
    array([[[[[  1.217539e+03, ...,   1.001117e+03],
          ..., 
          [  9.172734e+02, ...,   1.903422e+03]],

         ..., 
         [[  2.362061e+00, ...,   1.156982e+00],
          ..., 
          [  1.715332e+00, ...,   4.047852e+00]]]]], dtype=float32)
    Coordinates:
        XLONG     (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT      (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
      * file      (file) int64 0 1 2
        datetime  (file, Time) datetime64[ns] 2016-10-07 2016-10-07 2016-10-07
    Dimensions without coordinates: Time, bottom_top, south_north, west_east
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)

                   
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
    array([[[[[  1.217539e+03, ...,   1.001117e+03],
              ..., 
              [  9.172734e+02, ...,   1.903422e+03]],
    
             ..., 
             [[  2.362061e+00, ...,   1.156982e+00],
              ..., 
              [  1.715332e+00, ...,   4.047852e+00]]]]], dtype=float32)
    Coordinates:
        XLONG     (south_north, west_east) float32 -122.72 -122.693 -122.666 ...
        XLAT      (south_north, west_east) float32 21.1381 21.1451 21.1521 ...
      * Time      (Time) datetime64[ns] 2016-10-07 ...
        datetime  (Time) datetime64[ns] 2016-10-07 ...
      * key_0     (key_0) <U6 u'label1' u'label2'
    Dimensions without coordinates: bottom_top, south_north, west_east
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: perturbation pressure
        units: Pa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
                        
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
    Dimensions without coordinates: south_north, west_east
    Attributes:
        FieldType: 104
        units: m
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    Dimensions without coordinates: idx
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    Dimensions without coordinates: idx
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    
    # When using lat/lon coordinates, you must supply a WRF netcdf file object, 
    # or a projection object with the lower left latitude and lower left 
    # longitude points.
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
    Dimensions without coordinates: idx
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    Dimensions without coordinates: idx
    Attributes:
        FieldType: 104
        description: pressure
        units: hPa
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    Dimensions without coordinates: line_idx
    Attributes:
        FieldType: 104
        description: TEMP at 2 M
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    Dimensions without coordinates: line_idx
    Attributes:
        FieldType: 104
        description: TEMP at 2 M
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    
    # Calculate the interpolated line.  To use latitude and longitude points, 
    # you must supply a WRF NetCDF file object, or a projection object along 
    # with the lower left latitude and lower left longitude points. 
    # Also, by setting latlon to True, this routine will 
    # also calculate the latitude and longitude coordinates along the line
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
    Dimensions without coordinates: line_idx
    Attributes:
        FieldType: 104
        description: TEMP at 2 M
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
        orientation: (899.0, 0.0) to (899.0, 1057.0)    
            

Interpolating a 3D Field to a Surface Type 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`wrf.vinterp` is used to interpolate a field to a type of surface.  
The available surfaces are pressure, geopotential height, theta, and theta-e. 
The surface levels to interpolate also need to be specified.

.. code-block:: python

    from __future__ import print_function

    from netCDF4 import Dataset
    from wrf import getvar, vinterp 
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")  
    
    tk = getvar(ncfile, "tk")
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
      * interp_level  (interp_level) int64 200 300 500 1000
    Dimensions without coordinates: south_north, west_east
    Attributes:
        FieldType: 104
        MemoryOrder: XYZ
        description: temperature
        units: K
        stagger: 
        coordinates: XLONG XLAT
        projection: LambertConformal(stand_lon=-97.5, moad_cen_lat=38.5000038147, 
                                     truelat1=38.5, truelat2=38.5, pole_lat=90.0, 
                                     pole_lon=0.0)
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
    Dimensions without coordinates: idx
        
        
    <xarray.DataArray u'xy' (x_y: 2)>
    array([400, 200])
    Coordinates:
        latlon_coord  object CoordPair(lat=28.5581640822, lon=-112.678276173)
      * x_y           (x_y) <U1 u'x' u'y'
    Dimensions without coordinates: idx
    
    
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
    Dimensions without coordinates: idx
        
        
    <xarray.DataArray u'xy' (x_y: 2, idx: 2)>
    array([[400, 105],
           [200, 205]])
    Coordinates:
        latlon_coord  (idx) object CoordPair(lat=28.5581640822, lon=-112.678276173) ...
      * x_y           (x_y) <U1 u'x' u'y'
    Dimensions without coordinates: idx


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
    
.. _using_omp:

Using OpenMP
-------------------------

Beginning in version 1.1, the Fortran computational routines in wrf-python make 
use of OpenMP directives. OpenMP enables the calculations to use multiple CPU 
cores, which can improve performance. In order to use OpenMP features, 
wrf-python has to be compiled with OpenMP enabled (most pre-built binary 
installations will have this enabled).

The Fortran computational routines have all been built using runtime 
scheduling, instead of compile time scheduling, so that the user can choose the 
scheduler type within their Python application. By default, the scheduling 
type is set to :data:`wrf.OMP_SCHED_STATIC` using only 1 CPU core, so 
wrf-python will behave similarly to the non-OpenMP built versions. For the most 
part, the difference between the scheduling types is minimal, with the exception 
being the :data:`wrf.OMP_SCHED_DYNAMIC` scheduler that is much slower due to 
the additional overhead associated with it. For new users, using the default 
scheduler should be sufficient.


Verifying that OpenMP is Enabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To take advantage of the performance improvements offered by OpenMP, wrf-python
needs to have been compiled with OpenMP features enabled. The example below 
shows how you can determine if OpenMP is enabled in your build of wrf-python.

.. code-block:: python

   from __future__ import print_function

   from wrf import omp_enabled

   print(omp_enabled())


Result:

.. code-block:: none

   True


Determining the Number of Available Processors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The example below shows how you can get the maximum number of processors 
that are available on your system.

.. code-block:: python

   from __future__ import print_function

   from wrf import omp_get_num_procs

   print(omp_get_num_procs())


Result:

.. code-block:: none

   8


Specifying the Number of Threads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To enable multicore support via OpenMP, specifying the maximum number 
of OpenMP threads (i.e. CPU cores) is the only step that you need to take.  

In the example below, :meth:`wrf.omp_set_num_threads` is used to set the 
maximum number of threads to use, and :meth:`wrf.omp_get_max_threads` is used 
to retrieve (and print) the maximum number of threads used.

.. note::

   Although there is an OpenMP routine named :meth:`wrf.omp_get_num_threads`, 
   this routine will always return 1 when called from the sequential part of 
   the program. Use :meth:`wrf.omp_get_max_threads` to return the value set by 
   :meth:`wrf.omp_set_num_threads`.

.. code-block:: python

   from __future__ import print_function

   from wrf import omp_set_num_threads, omp_get_max_threads

   omp_set_num_threads(4)
   
   print (omp_get_max_threads())


Result:

.. code-block:: none

   4
   
Setting a Different Scheduler Type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When an OpenMP directive is encountered in the Fortran code, a scheduler is 
used to determine how the work is divided among the threads. All of the 
Fortran routines are compiled to use a 'runtime' scheduler, which indicates 
that the scheduler type (from the four listed below) is to be chosen at 
runtime (i.e. inside a Python script)

By default, the scheduler chosen is the :data:`wrf.OMP_SCHED_STATIC` scheduler,
which should be sufficient for most users. However, OpenMP and wrf-python 
include the following options for the scheduler type:

- :data:`wrf.OMP_SCHED_STATIC`
- :data:`wrf.OMP_SCHED_DYNAMIC`
- :data:`wrf.OMP_SCHED_GUIDED`
- :data:`wrf.OMP_SCHED_AUTO`

Refer to the 
`OpenMP Specification <http://www.openmp.org/wp-content/uploads/openmp-4.5.pdf>`_.
for more information about these scheduler types. In local testing, 
:data:`wrf.OMP_SCHED_GUIDED` produced the best results, but 
differences between :data:`wrf.OMP_SCHED_STATIC`, 
:data:`wrf.OMP_SCHED_GUIDED`, and 
:data:`wrf.OMP_SCHED_AUTO` were minor. However, 
:data:`wrf.OMP_SCHED_DYNAMIC` produced noticeably slower results 
due to the overhead of using a dynamic scheduler.

When setting a scheduler type, the :meth:`wrf.omp_set_schedule` takes two 
arguments.  The first is the scheduler type (one from the list above), and the 
second optional argument is a modifier, which is usually referred as the chunk 
size. If the modifier/chunk_size is set to 0, then the OpenMP default 
implementation is used. For :data:`wrf.OMP_SCHED_AUTO`, the 
modifier is ignored.

If you are new to OpenMP and all this sounds confusing, don't worry about 
setting a scheduler type.  The default static scheduler will be good enough.

In the example below, the scheduler type is set to 
:data:`wrf.OMP_SCHED_GUIDED` and uses the default chunk size of 0. The 
scheduler type is then read back using :meth:`wrf.omp_get_schedule` 
and printed.

.. code-block:: python

   from __future__ import print_function

   from wrf import omp_set_schedule, omp_get_schedule, OMP_SCHED_GUIDED

   omp_set_schedule(OMP_SCHED_GUIDED, 0)

   sched, modifier = omp_get_schedule()

   print(sched, modifier)


Result:

.. code-block:: none

   3 1
   
Notice that the printed scheduler type (*sched* variable) is set to a 
value of 3, which is the actual integer constant value for the 
:data:`wrf.OMP_SCHED_GUIDED` scheduler type. The *modifier* is returned as a 
value of 1, which is different than the 0 that was supplied to the 
:meth:`wrf.omp_set_schedule` routine. This is because the 0 tells OpenMP to use 
its own default value for the scheduler, which is 1 for this type of scheduler.

.. _performance:

Performance Tips
--------------------

Memory Issues with :data:`wrf.ALL_TIMES` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The use of :data:`wrf.ALL_TIMES` for the *timeidx* parameter to 
:meth:`wrf.getvar` is convenient for computing diagnostic variables across 
multiple files/times, but there is something that users should be aware of. 
When :data:`wrf.ALL_TIMES` is set as the *timeidx* argument, all arrays used 
in the computation are extracted for all times before the computation 
is started. This can cause serious memory issues on smaller hardware systems 
like laptops.  

In this example, the user wants to use a data set that is 289 x 39 x 300 x 300
and compute z for the entire data set. The user is using a laptop with 
8 GB of memory.  

.. code-block:: python

   from netCDF4 import Dataset
   from wrf import getvar, ALL_TIMES
   
   file_list = [Dataset("/path/to/file1"), Dataset("/path/to/file2"),...]
   z = getvar(file_list, "z", ALL_TIMES)
   
Five hours later, the computation finished. What happened?

In wrf-python, all of the computational routines use 8-byte REAL variables so 
that both the 4-byte and 8-byte version of WRF output can be used. The 
calculation for z extracts three variables (P, PHB, and HGT) and returns a 
fourth array (RESULT). The RESULT will get cut in half to 4-byte REALs 
after the computation, but needs an 8-byte REAL when the result is computed. 

Let's look at the approximate amount memory needed:

**P**: 289 x 39 x 300 x 300 x 8 = 8,115,120,000 bytes (~8 GB!)    

**PHB**: 289 x 39 x 300 x 300 x 8 = 8,115,120,000 bytes (~8 GB!)    

**HGT**: 289 x 300 x 300 x 8 = 208,080,000 (~208 MB)

**RESULT**: 289 x 39 x 300 x 300 x 8 = 8,115,120,000 bytes (~8 GB!)

Yikes! So, in order to do this calculation using :data:`wrf.ALL_TIMES` as 
the *timeidx*, over 24.2 GB are needed for this one calculation. When the 
laptop runs out of memory, it begins using the hard drive for swap memory, 
which runs hundreds of times slower than real memory.

To fix this situation, it is better to allocate the output array yourself and 
run the calculation for each time step in a loop 
("loop-and-fill"). The required memory requirements change to:

(Note: only need to store the result in a 4-byte REAL)

**FINAL_RESULT**: 289 x 39 x 300 x 300 x 4 = 4,057560,000 bytes (~4 GB)

(Note: the numbers below are for each loop iteration)

**P**: 39 x 300 x 300 x 8 = 28,080,000 bytes (~28 MB)

**PHB**: 39 x 300 x 300 x 8 = 28,080,000 bytes (~28 MB)

**HGT**: 300 x 300 x 8 = 720,000 bytes (720 KB)

**RESULT**: 39 x 300 x 300 x 8 = 28,080,000 bytes (~28 MB)

Since the memory for the computation is deleted after each 
loop iteration, the total memory usage drops to approximately 4.1 GB.

The moral of the story is that you need to make sure that your system has 
enough memory to extract everything it needs up front if you want to use 
:data:`wrf.ALL_TIMES`, otherwise it is better to "loop-and-fill" yourself.

Here is an example of the "loop-and-fill" technique:

.. code-block:: python
    
   from __future__ import print_function, division
   
   import numpy as np
   from netCDF4 import Dataset
   from wrf import getvar, ALL_TIMES
   
   filename_list = ["/path/to/file1", "/path/to/file2",...]
   
   # Result shape (hard coded for this example)
   result_shape = (289, 39, 300, 300)
   
   # Only need 4-byte floats
   z_final = np.empty(result_shape, np.float32)
   
   # Modify this number if using more than 1 time per file
   times_per_file = 1
   
   for timeidx in range(result_shape[0]):
       # Compute the file index and the time index inside the file
       fileidx = timeidx // times_per_file
       file_timeidx = timeidx % times_per_file
       
       f = Dataset(filename_list[fileidx])  
       z = getvar(f, "z", file_timeidx)
       
       z_final[timeidx,:] = z[:]
       f.close()

      
The *cache* Argument for :meth:`wrf.getvar`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     

If you have read through the documentation, you may have noticed that the 
:meth:`wrf.getvar` routine contains a *cache* argument. What is this for?

Internally, if metadata is turned on, a variable is extracted from the NetCDF 
file and its metadata is copied to form the result's metadata. Often this 
variable is one of the computation's function arguments, so rather than 
spend time extracting the variable again for the computation, it is placed 
in a cache (dictionary) and passed on to the computational function.

What isn't widely known is that this cache argument can also be supplied by 
end users wishing to speed up their application. This can be useful in 
situations where numerous calculations are being performed on the same 
data set. For many algorithms, the time needed to extract the arrays from the 
NetCDF file is on par with the time needed to perform the calculation. If you 
are computing numerous diagnostics, extracting the variables up front allows 
you to only pay this extraction penalty once, rather than inside of each call 
to :meth:`wrf.getvar`.

The cache is nothing more than a dictionary where each key is the variable 
name (e.g. "P") and the value is the :class:`xarray.DataArray` or 
:class:`numpy.ndarray` variable. Creating the cache dictionary is easy, 
since the :meth:`wrf.extract_vars` routine returns a dictionary for a 
sequence of variables. 

.. note:: 

   The *timeidx* parameter supplied to :meth:`extract_vars` 
   must be the same *timeidx* parameter that you plan to use for 
   :meth:`wrf.getvar`. Otherwise, it will crash with dimension mismatch errors.

Some common variables that you can use to create an effective cache are: P, PB, 
PH, PHB, T, QVAPOR, HGT, PSFC, U, V, W.

Below is an example showing the same computation done with and without the 
cache. The execution time is printed. The hardware used is a 2.8 GHz Intel Core 
i7, which contains 4 CPU cores with 2 hyper threads (8 total threads). This 
will be interpreted as 8 CPUs for OpenMP.

.. code-block:: python

   from __future__ import print_function

   import time
   from netCDF4 import Dataset
   from wrf import getvar, ALL_TIMES, extract_vars

   # The first two files contain four times, the last file contains only one.
   wrf_filenames = ["/path/to/wrfout_d02_2005-08-28_00:00:00",
                    "/path/to/wrfout_d02_2005-08-28_12:00:00", 
                    "/path/to/wrfout_d02_2005-08-29_00:00:00"]

   wrfin = [Dataset(x) for x in wrf_filenames]

   start = time.time()
   my_cache = extract_vars(wrfin, ALL_TIMES, ("P", "PSFC", "PB", "PH", "PHB", 
                                              "T", "QVAPOR", "HGT", "U", "V", 
                                              "W"))
   end = time.time()
   print ("Time taken to build cache: ", (end-start), "s")
    
   vars = ("avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
           "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
           "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
           "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
           "wa", "uvmet10", "uvmet", "z", "cfrac", "zstag", "geopt_stag")
   
   # No cache
   start = time.time()
   for var in vars:
       v = getvar(wrfin, var, ALL_TIMES)
   end = time.time()
   no_cache_time = (end-start)

   print ("Time taken without variable cache: ", no_cache_time, "s")

   # With a cache
   start = time.time()
   for var in vars:
       v = getvar(wrfin, var, ALL_TIMES, cache=my_cache)
   end = time.time()
   cache_time = (end-start)

   print ("Time taken with variable cache: ", cache_time, "s")

   improvement = ((no_cache_time-cache_time)/no_cache_time) * 100 
   print ("The cache decreased computation time by: ", improvement, "%")


Result:

.. code-block:: none

   Time taken to build cache:  0.28154706955 s
   Time taken without variable cache:  11.0905270576 s
   Time taken with variable cache:  8.25931215286 s
   The cache decreased computation time by:  25.5282268378 %
   
By removing the repeated extraction of common variables in the 
:meth:`wrf.getvar` routine, for the single threaded case, the computation 
time has been reduced by 25.5% in this particular example.

Things get more interesting when OpenMP is turned on and set to use the 
maximum number of processors (in this case 8 threads are used).  

.. code-block:: python

   from __future__ import print_function

   import time
   from netCDF4 import Dataset
   from wrf import (getvar, ALL_TIMES, extract_vars, 
                    omp_set_num_threads, omp_get_num_procs)
   
   # The first two files contain four times, the last file contains only one.
   wrf_filenames = ["/path/to/wrfout_d02_2005-08-28_00:00:00",
                    "/path/to/wrfout_d02_2005-08-28_12:00:00", 
                    "/path/to/wrfout_d02_2005-08-29_00:00:00"]

   wrfin = [Dataset(x) for x in wrf_filenames]

   start = time.time()
   my_cache = extract_vars(wrfin, ALL_TIMES, ("P", "PSFC", "PB", "PH", "PHB", 
                                              "T", "QVAPOR", "HGT", "U", "V", 
                                              "W"))
   end = time.time()
   print ("Time taken to build cache: ", (end-start), "s")

   omp_set_num_threads(omp_get_num_procs())
   
   vars = ("avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", 
           "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", 
           "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", 
           "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", 
           "wa", "uvmet10", "uvmet", "z", "cfrac", "zstag", "geopt_stag")
   
   # No cache
   start = time.time()
   for var in vars:
       v = getvar(wrfin, var, ALL_TIMES)
   end = time.time()
   no_cache_time = (end-start)

   print ("Time taken without variable cache: ", no_cache_time, "s")

   # With a cache
   start = time.time()
   for var in vars:
       v = getvar(wrfin, var, ALL_TIMES, cache=my_cache)
   end = time.time()
   cache_time = (end-start)

   print ("Time taken with variable cache: ", cache_time, "s")

   improvement = ((no_cache_time-cache_time)/no_cache_time) * 100 
   print ("The cache decreased computation time by: ", improvement, "%")

Result:

.. code-block:: none

   Time taken to build cache:  0.2700548172 s
   Time taken without variable cache:  6.02652812004 s
   Time taken with variable cache:  3.27777099609 s
   The cache decreased computation time by:  45.6109565772 %
   
In this example, 4 CPU cores (8 total threads) are used. When the cache is 
used, the computation time drops by 45%, so almost half the time was spent 
simply extracting variables from the NetCDF file. When compared to the 
11.09 s needed to compute the single threaded case with no variable cache, the
computation time drops by roughly 70% (compared to 45% with 8 threads but 
no cache). 

In summary, if you are computing a lot of diagnostic variables, consider using
the *cache* argument to improve performance, particularly if you want to 
maximize your multithreaded performance with OpenMP.
