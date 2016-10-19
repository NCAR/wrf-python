Internal API
=============

Algorithm Decorators
--------------------

The decorators below are used for performing common operations related to  
diagnostic and interpolation calculations.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   wrf.convert_units
   wrf.left_iteration
   wrf.cast_type
   wrf.extract_and_transpose
   wrf.check_args
   wrf.uvmet_left_iter
   wrf.cape_left_iter
   wrf.cloudfrac_left_iter

  
Metadata Decorators
--------------------

The decorators below are used for performing common operations related to 
setting metadata.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.copy_and_set_metadata
   wrf.set_wind_metadata
   wrf.set_cape_metadata
   wrf.set_cloudfrac_metadata
   wrf.set_latlon_metadata
   wrf.set_height_metadata
   wrf.set_interp_metadata
   wrf.set_alg_metadata
   wrf.set_uvmet_alg_metadata
   wrf.set_cape_alg_metadata
   wrf.set_cloudfrac_alg_metadata
   wrf.set_destag_metadata
   
   
Decorator Utilities
--------------------

The routines below are used within decorators.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.either
   wrf.combine_dims
   wrf.from_var
   wrf.from_args
   wrf.args_to_list
   wrf.arg_location
   
Miscellaneous Classes
----------------------


.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.IterWrapper
   
   
