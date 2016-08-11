undef("set_mp_wrf_map_resources")
function set_mp_wrf_map_resources(in_file[1]:file,opt_args[1]:logical)   

begin
;
    opts = opt_args      ; Make a copy of the resource list

; Set some resources depending on what kind of map projection is 
; chosen.
;
;   MAP_PROJ = 0 : "CylindricalEquidistant"
;   MAP_PROJ = 1 : "LambertConformal"
;   MAP_PROJ = 2 : "Stereographic"
;   MAP_PROJ = 3 : "Mercator"
;   MAP_PROJ = 6 : "Lat/Lon"

    if(isatt(in_file,"MAP_PROJ"))

;   CylindricalEquidistant
      if(in_file@MAP_PROJ .eq. 0)
        projection          = "CylindricalEquidistant"
        opts@mpProjection = projection
        opts@mpGridSpacingF = 45
        opts@mpCenterLatF   = get_res_value_keep(opts, "mpCenterLatF", 0.0)
        if(isatt(in_file,"STAND_LON"))
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if

;   LambertConformal projection
      if(in_file@MAP_PROJ .eq. 1)
        projection               = "LambertConformal"
        opts@mpProjection = projection
        opts@mpLambertParallel1F = get_res_value_keep(opts, "mpLambertParallel1F",in_file@TRUELAT1)
        opts@mpLambertParallel2F = get_res_value_keep(opts, "mpLambertParallel2F",in_file@TRUELAT2)
        if(isatt(in_file,"STAND_LON"))
          opts@mpLambertMeridianF  = get_res_value_keep(opts, "mpLambertMeridianF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            opts@mpLambertMeridianF  = get_res_value_keep(opts, "mpLambertMeridianF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if
  
;   Stereographic projection
      if(in_file@MAP_PROJ .eq. 2)
        projection          = "Stereographic"
        opts@mpProjection = projection
        opts@mpCenterLatF   = get_res_value_keep(opts, "mpCenterLatF", in_file@CEN_LAT)
        if(isatt(in_file,"STAND_LON"))
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if
  
;   Mercator projection
      if(in_file@MAP_PROJ .eq. 3)
        projection          = "Mercator"
        opts@mpProjection = projection
        opts@mpCenterLatF   = get_res_value_keep(opts, "mpCenterLatF", 0.0)
        if(isatt(in_file,"STAND_LON"))
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if
  
;    global WRF CylindricalEquidistant
      if(in_file@MAP_PROJ .eq. 6)
        projection          = "CylindricalEquidistant"
        opts@mpProjection = projection
        opts@mpGridSpacingF = 45

        if (isatt(in_file,"POLE_LON") .and. isatt(in_file,"POLE_LAT") .and. isatt(in_file,"STAND_LON")) then

          if (in_file@POLE_LON .eq. 0 .and. in_file@POLE_LAT .eq. 90) then
            ; not rotated

            opts@mpCenterLatF   = get_res_value_keep(opts, "mpCenterLatF", 0.0)
            opts@mpCenterLonF   = get_res_value_keep(opts, "mpCenterLonF",180 - in_file@STAND_LON)

          else 
            ; rotated

            southern = False ; default to northern hemisphere
            if (in_file@POLE_LON .eq. 0.0) then
              southern = True
            else if (in_file@POLE_LON .ne. 180) then
              if (isatt(in_file,"CEN_LAT") .and. in_file@CEN_LAT .lt. 0.0) then
                southern = True  ; probably but not necessarily true -- no way to tell for sure
              end if
            end if
            end if

            if (.not. southern) then
              opts@mpCenterLatF = get_res_value_keep(opts, "mpCenterLatF", 90.0 - in_file@POLE_LAT) 
              opts@mpCenterLonF = get_res_value_keep(opts, "mpCenterLonF", -in_file@STAND_LON) 
            else 
              opts@mpCenterLatF = get_res_value_keep(opts, "mpCenterLatF", in_file@POLE_LAT - 90) 
              opts@mpCenterLonF = get_res_value_keep(opts, "mpCenterLonF", 180 - in_file@STAND_LON) 
            end if

          end if

        else if (isatt(in_file,"ref_lon") .and. isatt(in_file,"ref_lat")) then
          ;; this is definitely true for NMM grids but unlikely for ARW grids especially if ref_x and ref_y are set
          opts@mpCenterLatF = get_res_value_keep(opts, "mpCenterLatF", in_file@REF_LAT) 
          opts@mpCenterLonF = get_res_value_keep(opts, "mpCenterLonF", in_file@REF_LON) 

        else if (isatt(in_file,"cen_lat") .and. isatt(in_file,"cen_lon")) then
          ;; these usually specifiy the center of the coarse domain --- not necessarily the center of the projection
          opts@mpCenterLatF  = get_res_value_keep(opts, "mpCenterLatF",in_file@CEN_LAT)
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@CEN_LON)

        else 
          ;; default values for global grid
          opts@mpCenterLatF  = get_res_value_keep(opts, "mpCenterLatF", 0.0)
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF", 180.0)

        end if
        end if
        end if
      end if
  
    end if

  return(opts)                                     ; Return.

end


undef("wrf_map_resources")
function wrf_map_resources(in_file[1]:file,map_args[1]:logical)   
local lat, lon, x1, x2, y1, y2, dims, ii, jj, southern
begin
;
; This function sets resources for a WRF map plot, basing the projection on
; the MAP_PROJ attribute in the given file. It's intended to be callable
; by users who need to set mpXXXX resources for other plotting scripts.
;

; Set some resources depending on what kind of map projection is 
; chosen.
;
;   MAP_PROJ = 0 : "CylindricalEquidistant"
;   MAP_PROJ = 1 : "LambertConformal"
;   MAP_PROJ = 2 : "Stereographic"
;   MAP_PROJ = 3 : "Mercator"
;   MAP_PROJ = 6 : "Lat/Lon"

    if(isatt(in_file,"MAP_PROJ"))

;   CylindricalEquidistant
      if(in_file@MAP_PROJ .eq. 0)
        map_args@mpProjection          = "CylindricalEquidistant"
        map_args@mpGridSpacingF = 45
        map_args@mpCenterLatF   = get_res_value_keep(map_args, "mpCenterLatF", 0.0)
        if(isatt(in_file,"STAND_LON"))
          map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if

;   LambertConformal projection
      if(in_file@MAP_PROJ .eq. 1)
        map_args@mpProjection               = "LambertConformal"
        map_args@mpLambertParallel1F = get_res_value_keep(map_args, "mpLambertParallel1F",in_file@TRUELAT1)
        map_args@mpLambertParallel2F = get_res_value_keep(map_args, "mpLambertParallel2F",in_file@TRUELAT2)
        if(isatt(in_file,"STAND_LON"))
          map_args@mpLambertMeridianF  = get_res_value_keep(map_args, "mpLambertMeridianF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            map_args@mpLambertMeridianF  = get_res_value_keep(map_args, "mpLambertMeridianF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if
  
;   Stereographic projection
      if(in_file@MAP_PROJ .eq. 2)
        map_args@mpProjection          = "Stereographic"
        map_args@mpCenterLatF   = get_res_value_keep(map_args, "mpCenterLatF", in_file@CEN_LAT)
        if(isatt(in_file,"STAND_LON"))
          map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if
  
;   Mercator projection
      if(in_file@MAP_PROJ .eq. 3)
        map_args@mpProjection          = "Mercator"
        map_args@mpCenterLatF   = get_res_value_keep(map_args, "mpCenterLatF", 0.0)
        if(isatt(in_file,"STAND_LON"))
          map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@STAND_LON)
        else
          if(isatt(in_file,"CEN_LON"))
            map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@CEN_LON)
          else
           print("ERROR: Found neither STAND_LON or CEN_LON in file")
          end if
        end if
      end if
  
;    global WRF CylindricalEquidistant
      if(in_file@MAP_PROJ .eq. 6)
        print ("YES, THIS WORKED")
        projection          = "CylindricalEquidistant"
        map_args@mpProjection = projection
        map_args@mpGridSpacingF = 45
        
        ;; according to the docs if POLE_LON is 0 then the projection center is in the southern hemisphere
        ;; if POLE_LON is 180 the projection center is in the northern hemisphere
        ;; otherwise you can't tell for sure -- CEN_LAT does not have to be the projection center but hopefully
        ;; it is in the same hemisphere. The same is true for REF_LAT except that if REF_Y is specified REF_LAT might
        ;; be in a corner or somewhere else and therefore it is even less reliable
        ;; 

        if (isatt(in_file,"POLE_LON") .and. isatt(in_file,"POLE_LAT") .and. isatt(in_file,"STAND_LON")) then

          if (in_file@POLE_LON .eq. 0 .and. in_file@POLE_LAT .eq. 90) then
            ; not rotated

            map_args@mpCenterLatF   = get_res_value_keep(map_args, "mpCenterLatF", 0.0)
            map_args@mpCenterLonF   = get_res_value_keep(map_args, "mpCenterLonF",180 - in_file@STAND_LON)

          else 
            ; rotated

            southern = False ; default to northern hemisphere
            if (in_file@POLE_LON .eq. 0.0) then
              southern = True
            else if (in_file@POLE_LON .ne. 180) then
              if (isatt(in_file,"CEN_LAT") .and. in_file@CEN_LAT .lt. 0.0) then
                southern = True  ; probably but not necessarily true -- no way to tell for sure
              end if
            end if
            end if

            if (.not. southern) then
              map_args@mpCenterLatF = get_res_value_keep(map_args, "mpCenterLatF", 90.0 - in_file@POLE_LAT) 
              map_args@mpCenterLonF = get_res_value_keep(map_args, "mpCenterLonF", -in_file@STAND_LON) 
            else 
              map_args@mpCenterLatF = get_res_value_keep(map_args, "mpCenterLatF", in_file@POLE_LAT - 90) 
              map_args@mpCenterLonF = get_res_value_keep(map_args, "mpCenterLonF", 180 - in_file@STAND_LON) 
            end if

          end if

        else if (isatt(in_file,"ref_lon") .and. isatt(in_file,"ref_lat")) then
          ;; this is definitely true for NMM grids but unlikely for ARW grids especially if ref_x and ref_y are set
          map_args@mpCenterLatF = get_res_value_keep(map_args, "mpCenterLatF", in_file@REF_LAT) 
          map_args@mpCenterLonF = get_res_value_keep(map_args, "mpCenterLonF", in_file@REF_LON) 

        else if (isatt(in_file,"cen_lat") .and. isatt(in_file,"cen_lon")) then
          ;; these usually specifiy the center of the coarse domain --- not necessarily the center of the projection
          map_args@mpCenterLatF  = get_res_value_keep(map_args, "mpCenterLatF",in_file@CEN_LAT)
          map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@CEN_LON)

        else 
          ;; default values for global grid
          map_args@mpCenterLatF  = get_res_value_keep(map_args, "mpCenterLatF", 0.0)
          map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF", 180.0)

        end if
        end if
        end if

      end if

    else
  
      return(map_args)
  
    end if

    map_args@mpNestTime = get_res_value_keep(map_args, "mpNestTime",0)

      if(isfilevar(in_file,"XLAT"))
        lat = in_file->XLAT(map_args@mpNestTime,:,:)
        lon = in_file->XLONG(map_args@mpNestTime,:,:)
      else
        lat = in_file->XLAT_M(map_args@mpNestTime,:,:)
        lon = in_file->XLONG_M(map_args@mpNestTime,:,:)
      end if
      dims = dimsizes(lat)

      do ii = 0, dims(0)-1
      do jj = 0, dims(1)-1
        if ( lon(ii,jj) .lt. 0.0) then
          lon(ii,jj) = lon(ii,jj) + 360.
        end if
      end do
      end do

      map_args@start_lat = lat(0,0)
      map_args@start_lon = lon(0,0)
      map_args@end_lat   = lat(dims(0)-1,dims(1)-1)
      map_args@end_lon   = lon(dims(0)-1,dims(1)-1)
      
      ; end_lon must be greater than start_lon, or errors are thrown
      if (map_args@end_lon .le. map_args@start_lon) then
          map_args@end_lon = map_args@end_lon + 360.0
      end if


; Set some resources common to all map projections.
      map_args = set_mp_resources(map_args)

      if ( isatt(map_args,"ZoomIn") .and. map_args@ZoomIn ) then
        y1 = 0
        x1 = 0
        y2 = dims(0)-1
        x2 = dims(1)-1
        if ( isatt(map_args,"Ystart") ) then
          y1 = map_args@Ystart
          delete(map_args@Ystart)
        end if
        if ( isatt(map_args,"Xstart") ) then
          x1 = map_args@Xstart
          delete(map_args@Xstart)
        end if
        if ( isatt(map_args,"Yend") ) then
          if ( map_args@Yend .le. y2 ) then
            y2 = map_args@Yend
          end if
          delete(map_args@Yend)
        end if
        if ( isatt(map_args,"Xend") ) then
          if ( map_args@Xend .le. x2 ) then
            x2 = map_args@Xend
          end if
          delete(map_args@Xend)
        end if

        map_args@mpLeftCornerLatF      = lat(y1,x1)
        map_args@mpLeftCornerLonF      = lon(y1,x1)
        map_args@mpRightCornerLatF     = lat(y2,x2)
        map_args@mpRightCornerLonF     = lon(y2,x2)
        
        if ( map_args@mpRightCornerLonF .lt. 0.0 ) then
          map_args@mpRightCornerLonF  = map_args@mpRightCornerLonF + 360.0
        end if 
        
        if ( map_args@mpRightCornerLonF .le. map_args@mpRightCornerLonF ) then
          map_args@mpRightCornerLonF  = map_args@mpRightCornerLonF + 360.0
        end if

        delete(map_args@ZoomIn)
      end if

      return(map_args)
end