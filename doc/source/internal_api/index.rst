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
   
   wrf.g_cape.get_2dcape
   wrf.g_cape.get_3dcape
   wrf.g_cloudfrac.get_cloudfrac
   wrf.g_ctt.get_ctt
   wrf.g_dbz.get_dbz
   wrf.g_dbz.get_max_dbz
   wrf.g_dewpoint.get_dp
   wrf.g_dewpoint.get_dp_2m
   wrf.g_geoht.get_geopt
   wrf.g_geoht.get_height
   wrf.g_geoht.get_height_agl
   wrf.g_geoht.get_stag_geopt
   wrf.g_geoht.get_stag_height
   wrf.g_helicity.get_srh
   wrf.g_helicity.get_uh
   wrf.g_omega.get_omega
   wrf.g_pressure.get_pressure
   wrf.g_pressure.get_pressure_hpa
   wrf.g_pw.get_pw
   wrf.g_rh.get_rh
   wrf.g_rh.get_rh_2m
   wrf.g_slp.get_slp
   wrf.g_temp.get_theta
   wrf.g_temp.get_temp
   wrf.g_temp.get_eth
   wrf.g_temp.get_tv
   wrf.g_temp.get_tw
   wrf.g_temp.get_tk
   wrf.g_temp.get_tc
   wrf.g_terrain.get_terrain
   wrf.g_times.get_times
   wrf.g_times.get_xtimes
   wrf.g_uvmet.get_uvmet
   wrf.g_uvmet.get_uvmet10
   wrf.g_uvmet.get_uvmet_wspd_wdir
   wrf.g_uvmet.get_uvmet10_wspd_wdir
   wrf.g_vorticity.get_avo
   wrf.g_vorticity.get_pvo
   wrf.g_wind.get_u_destag
   wrf.g_wind.get_v_destag
   wrf.g_wind.get_w_destag
   wrf.g_wind.get_destag_wspd_wdir
   wrf.g_wind.get_destag_wspd_wdir10
   wrf.g_wind.get_destag_wspd
   wrf.g_wind.get_destag_wdir
   wrf.g_wind.get_destag_wspd10
   wrf.g_wind.get_destag_wdir10
   wrf.g_uvmet.get_uvmet_wspd
   wrf.g_uvmet.get_uvmet_wdir
   wrf.g_uvmet.get_uvmet10_wspd
   wrf.g_uvmet.get_uvmet10_wdir
   wrf.g_cloudfrac.get_low_cloudfrac
   wrf.g_cloudfrac.get_mid_cloudfrac
   wrf.g_cloudfrac.get_high_cloudfrac
   wrf.g_cape.get_cape2d_only
   wrf.g_cape.get_cin2d_only
   wrf.g_cape.get_lcl
   wrf.g_cape.get_lfc
   wrf.g_cape.get_3dcape_only
   wrf.g_cape.get_3dcin_only
   
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
   
   
