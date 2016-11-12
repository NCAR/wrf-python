Internal API
=============

Routines
-------------

Extraction and Diagnostic Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The routines below are called internally by :meth:`wrf.getvar`.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.cape.get_2dcape
   wrf.cape.get_3dcape
   wrf.cloudfrac.get_cloudfrac
   wrf.ctt.get_ctt
   wrf.dbz.get_dbz
   wrf.dbz.get_max_dbz
   wrf.dewpoint.get_dp
   wrf.dewpoint.get_dp_2m
   wrf.geoht.get_geopt
   wrf.geoht.get_height
   wrf.helicity.get_srh
   wrf.helicity.get_uh
   wrf.omega.get_omega
   wrf.pressure.get_pressure
   wrf.pressure.get_pressure_hpa
   wrf.pw.get_pw
   wrf.rh.get_rh
   wrf.rh.get_rh_2m
   wrf.slp.get_slp
   wrf.temp.get_theta
   wrf.temp.get_temp
   wrf.temp.get_eth
   wrf.temp.get_tv
   wrf.temp.get_tw
   wrf.temp.get_tk
   wrf.temp.get_tc
   wrf.terrain.get_terrain
   wrf.times.get_times
   wrf.times.get_xtimes
   wrf.uvmet.get_uvmet
   wrf.uvmet.get_uvmet10
   wrf.uvmet.get_uvmet_wspd_wdir
   wrf.uvmet.get_uvmet10_wspd_wdir
   wrf.vorticity.get_avo
   wrf.vorticity.get_pvo
   wrf.wind.get_u_destag
   wrf.wind.get_v_destag
   wrf.wind.get_w_destag
   wrf.wind.get_destag_wspd_wdir
   wrf.wind.get_destag_wspd_wdir10
   
-------------------------

Decorators
----------------


Algorithm Decorators
^^^^^^^^^^^^^^^^^^^^^^^^

The decorators below are used for performing common operations related to  
diagnostic and interpolation calculations.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/

   wrf.decorators.convert_units
   wrf.decorators.left_iteration
   wrf.decorators.cast_type
   wrf.decorators.extract_and_transpose
   wrf.decorators.check_args
   wrf.specialdec.uvmet_left_iter
   wrf.specialdec.cape_left_iter
   wrf.specialdec.cloudfrac_left_iter

  
Metadata Decorators
^^^^^^^^^^^^^^^^^^^^^^

The decorators below are used for performing common operations related to 
setting metadata.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.metadecorators.copy_and_set_metadata
   wrf.metadecorators.set_wind_metadata
   wrf.metadecorators.set_cape_metadata
   wrf.metadecorators.set_cloudfrac_metadata
   wrf.metadecorators.set_latlon_metadata
   wrf.metadecorators.set_height_metadata
   wrf.metadecorators.set_interp_metadata
   wrf.metadecorators.set_alg_metadata
   wrf.metadecorators.set_uvmet_alg_metadata
   wrf.metadecorators.set_cape_alg_metadata
   wrf.metadecorators.set_cloudfrac_alg_metadata
   wrf.metadecorators.set_destag_metadata
   
   
Decorator Utilities
^^^^^^^^^^^^^^^^^^^^^^^

The routines below are used within the decorators.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.either
   wrf.combine_dims
   wrf.from_var
   wrf.from_args
   wrf.args_to_list
   wrf.arg_location

  
------------------------

Classes
-----------------------
  
Iterable Wrapper Class
^^^^^^^^^^^^^^^^^^^^^^^

The class below is an Iterable wrapper class and provides an __iter__ function 
that always returns the beginning of the sequence, regardless of the 
Iterable type.

.. autosummary::
   :nosignatures:
   :toctree: ./generated/
   
   wrf.IterWrapper
   
   
