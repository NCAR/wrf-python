User API
=============

Routines
------------------

Diagnostic Routine
^^^^^^^^^^^^^^^^^^^^^^^

The routine below is the primary routine for extracting variables from a 
WRF-ARW NetCDF file (or sequence of files) and performing diagnostic 
calculations.  

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   wrf.getvar
   
   
Interpolation Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are the primary routines used for performing interpolation 
calculations.  

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.interplevel
   wrf.vertcross
   wrf.interpline
   wrf.vinterp
   
Lat-Lon to/from XY Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are used for converting back and forth between xy-grid 
space and latitude-longitude space.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.ll_to_xy
   wrf.xy_to_ll
   wrf.ll_to_xy_proj
   wrf.xy_to_ll_proj
   
   
Grid Destaggering Routine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The routine below is used to convert a variable on a staggered grid to the 
unstaggered grid.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.destagger


Numpy Extraction Routine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The routine below is used to extract a :class:`numpy.ndarray` from a 
:class:`xarray.DataArray`.  This routine must be used before passing 
the array object to a compiled extension.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.to_np
    
Variable Extraction Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are primarily used internally by :meth:`wrf.getvar`, but 
some users may find them useful to manually extract variables from a 
WRF NetCDF file (or a sequence of NetCDF files).

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

    wrf.extract_vars
    wrf.combine_files
    wrf.extract_dim
    wrf.extract_global_attrs
    wrf.extract_times
    
    
Plotting Helper Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are used to assist with plotting.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
    wrf.geo_bounds
    wrf.latlon_coords
    wrf.get_cartopy
    wrf.get_basemap
    wrf.get_pyngl
    wrf.cartopy_xlim
    wrf.cartopy_ylim

Raw Diagnostic Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below can be used when working with variables that are not 
contained in a WRF-ARW NetCDF file.  They can also be used with non-WRF data.
However, if you are working with WRF-ARW NetCDF files, 
use :meth:`wrf.getvar` instead.

Keep in mind that these routines were developed for WRF-ARW, so your mileage 
may vary when working with non-WRF data.  Also, the vast majority of these 
routines do not allow for missing values in any of the input arrays, so make 
sure they are removed before calling these routines.


.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.xy
   wrf.interp1d
   wrf.interp2dxy
   wrf.interpz3d
   wrf.slp
   wrf.tk
   wrf.td
   wrf.rh
   wrf.uvmet
   wrf.smooth2d
   wrf.cape_2d
   wrf.cape_3d
   wrf.cloudfrac
   wrf.ctt
   wrf.dbz
   wrf.srhel
   wrf.udhel
   wrf.avo
   wrf.pvo
   wrf.eth
   wrf.wetbulb
   wrf.tvirtual
   wrf.omega
   wrf.pw

OpenMP Runtime Library Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are the OpenMP runtime libraries that have been wrapped 
for wrf-python. The entire library (OpenMP 3.x) has been wrapped, but many of 
the routines are only useful inside of an OpenMP thread, so they aren't useful
from inside the Python interpreter. Also, the Fortran code in wrf-python is 
fairly simple in terms of threading, so features like nested threads aren't 
used. The documentation below is split in to the useful OpenMP functions and 
the less useful functions.

The documentation for each routine was taken directly from the 
`OpenMP Specification <http://www.openmp.org/wp-content/uploads/openmp-4.5.pdf>`_. 
Read the specification for more details about these routines.
 
Useful OpenMP Routines
*****************************

The routines below are useful when called from within a Python program. These 
routines handle setting the number of threads, setting up the scheduler, 
and timing.

It is also important to note that the OpenMP directives within the Fortran 
code all specify a runtime scheduler. This means that the user can control 
the type of scheduling to use from within their Python application by using the 
routines below.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.omp_enabled
   wrf.omp_set_num_threads
   wrf.omp_get_max_threads
   wrf.omp_get_num_procs
   wrf.omp_set_dynamic
   wrf.omp_get_dynamic
   wrf.omp_set_schedule
   wrf.omp_get_schedule
   wrf.omp_get_thread_limit
   wrf.omp_get_wtime
   wrf.omp_get_wtick
   
Less Useful OpenMP Routines
*******************************

The routines below are less useful because wrf-python does not use nested 
parallelism and some of the routines are only applicable when called from 
within an OpenMP thread.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.omp_get_num_threads
   wrf.omp_get_thread_num
   wrf.omp_in_parallel
   wrf.omp_set_nested
   wrf.omp_get_nested
   wrf.omp_set_max_active_levels
   wrf.omp_get_max_active_levels
   wrf.omp_get_level
   wrf.omp_get_ancestor_thread_num
   wrf.omp_get_team_size
   wrf.omp_get_active_level
   wrf.omp_in_final
   wrf.omp_init_lock
   wrf.omp_init_nest_lock
   wrf.omp_destroy_lock
   wrf.omp_destroy_nest_lock
   wrf.omp_set_lock
   wrf.omp_set_nest_lock
   wrf.omp_unset_lock
   wrf.omp_unset_nest_lock
   wrf.omp_test_lock
   wrf.omp_test_nest_lock

Configuration Routines
^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are used to configure wrf-python by enabling or 
disabling third party packages.  For the most part, these settings are 
configured automatically based on the presence of a third party package.  
However, disabling xarray can be useful when you want to turn off all metadata 
in one place.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.xarray_enabled
   wrf.disable_xarray
   wrf.enable_xarray
   wrf.cartopy_enabled
   wrf.disable_cartopy
   wrf.enable_cartopy
   wrf.basemap_enabled
   wrf.disable_basemap
   wrf.enable_basemap
   wrf.pyngl_enabled
   wrf.enable_pyngl
   wrf.disable_pyngl
   wrf.set_cache_size
   wrf.get_cache_size
   wrf.omp_enabled
   

Miscellaneous Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are primarily used internally, but some users may find 
them helpful for other purposes.  

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.is_time_coord_var
   wrf.get_coord_pairs
   wrf.is_multi_time_req
   wrf.is_multi_file
   wrf.has_time_coord
   wrf.is_mapping
   wrf.latlon_coordvars
   wrf.is_coordvar
   wrf.get_iterable
   wrf.is_moving_domain
   wrf.npbytes_to_str
   wrf.is_standard_wrf_var
   wrf.is_staggered
   wrf.get_left_indexes
   wrf.iter_left_indexes
   wrf.get_right_slices
   wrf.get_proj_params
   wrf.psafilepath
   wrf.get_id
   wrf.getproj
   wrf.cache_item
   wrf.get_cached_item
   wrf.ll_points
   wrf.pairs_to_latlon
 
 
------------------------

 
Classes
----------------------

Exceptions
^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   wrf.DiagnosticError
   

CoordPair Class
^^^^^^^^^^^^^^^^^^^^^^^

The class below is used for storing coordinate metadata from routines that 
use a single point for an (x, y) or (lat, lon) location. 

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.CoordPair
   
CoordPair Methods
************************

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.CoordPair.latlon_str
   wrf.CoordPair.xy_str
   
GeoBounds Class
^^^^^^^^^^^^^^^^^^^^^^^

The class below is used for specifying geographic boundaries. 

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.GeoBounds
   
Projection Classes
^^^^^^^^^^^^^^^^^^^^^^^^

The classes below are used to hold the projection information in the 
'projection' entry within a :attr:`xarray.DataArray.attrs` attribute.

Projection Base Class
******************************

The base class for all map projection types.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.WrfProj
   
Projection Base Class Methods
***********************************

The class methods for all projection types.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   wrf.WrfProj.basemap
   wrf.WrfProj.cartopy
   wrf.WrfProj.cartopy_xlim
   wrf.WrfProj.cartopy_ylim
   wrf.WrfProj.pyngl
   wrf.WrfProj.cf
   wrf.WrfProj.proj4
   
   
Projection Subclasses
*****************************

See :class:`wrf.WrfProj` for methods and attributes.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   wrf.NullProjection
   wrf.LambertConformal
   wrf.Mercator
   wrf.PolarStereographic
   wrf.LatLon
   wrf.RotatedLatLon
   
   
   
   
