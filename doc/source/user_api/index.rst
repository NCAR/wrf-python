User API
=============

Diagnostic Routine
-------------------

The routine below is the primary routine for extracting variables from a 
WRF-ARW NetCDF file (or sequence of files) and performing diagnostic 
calculations.  

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   wrf.getvar
   
   
Interpolation Routines
----------------------

The routines below are the primary routines used for performing interpolation 
calculations.  

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.interplevel
   wrf.vertcross
   wrf.interpline
   wrf.vinterp


Numpy Extraction Routine
--------------------------

The routine below is used to extract a :class:`numpy.ndarray` from a 
:class:`xarray.DataArray`.  This routine must be used before passing 
the array object to a compiled extension.  Otherwise, unusually crashes
may occur.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.npvalues
   
    
Variable Extraction Routines
-----------------------------

The routines below are primarily used internally by :meth:`wrf.getvar`, but 
some users may find them useful to manually extract variables from a 
WRF NetCDF file (or a sequence of NetCDF files).

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

    wrf.extract_vars
    wrf.combine_files
    wrf.extract_dim
    wrf.extract_extract_global_attrs
    wrf.extract_times
    
   
   
Raw Diagnostic Routines
------------------------

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
   wrf.temp
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


CoordPair Class
----------------------

The class below is used for coordinate metadata. 

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.CoordPair
   

Miscellaneous Routines
-----------------------

The routines below are primarily used internally, but some users may find 
them helpful.  

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
   
   
   
   