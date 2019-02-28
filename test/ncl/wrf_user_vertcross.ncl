;--------------------------------------------------------------------------------

undef("wrf_user_ll_to_xy")
function wrf_user_ll_to_xy( file_handle, longitude:numeric, latitude:numeric, \
                            opts_args:logical )
                            
; This is the same as wrf_user_ll_to_ij, but returns 0-based indexes

begin
;
; As of NCL V6.0.0, wrf_user_ll_to_ij can now handle a file 
; or a list of files.
;
  if(typeof(file_handle).eq."file") then
    ISFILE = True
    nc_file = file_handle
  elseif(typeof(file_handle).eq."list") then
    ISFILE = False
    nc_file = file_handle[0]
  else
    print("wrf_user_ll_to_xy: Error: the first argument must be a file or a list of files opened with addfile or addfiles")
    return
  end if

  opts = opts_args
  useT  = get_res_value(opts,"useTime",0)
  returnI= get_res_value(opts,"returnInt",True)

  res = True
  res@MAP_PROJ  = nc_file@MAP_PROJ
  res@TRUELAT1  = nc_file@TRUELAT1
  res@TRUELAT2  = nc_file@TRUELAT2
  res@STAND_LON = nc_file@STAND_LON
  res@DX        = nc_file@DX
  res@DY        = nc_file@DY

  if (res@MAP_PROJ .eq. 6) then
    res@POLE_LAT  = nc_file@POLE_LAT
    res@POLE_LON  = nc_file@POLE_LON
    res@LATINC = (res@DY*360.)/2.0/3.141592653589793/6370000.
    res@LONINC = (res@DX*360.)/2.0/3.141592653589793/6370000.
  else
    res@POLE_LAT = 90.0
    res@POLE_LON =  0.0
    res@LATINC = 0.0
    res@LONINC = 0.0
  end if

  if(isfilevar(nc_file,"XLAT"))
    if(ISFILE) then
      XLAT  = nc_file->XLAT(useT,:,:)
      XLONG = nc_file->XLONG(useT,:,:)
    else
      XLAT  = file_handle[:]->XLAT
      XLONG = file_handle[:]->XLONG
    end if
  else
    if(ISFILE) then
      XLAT  = nc_file->XLAT_M(useT,:,:)
      XLONG = nc_file->XLONG_M(useT,:,:)
    else
      XLAT  = file_handle[:]->XLAT_M
      XLONG = file_handle[:]->XLONG_M
    end if
  end if


  if(dimsizes(dimsizes(XLAT)).eq.2) then
; Rank 2
    res@REF_LAT = XLAT(0,0)
    res@REF_LON = XLONG(0,0)
  else
; Rank 3
    res@REF_LAT = XLAT(useT,0,0)
    res@REF_LON = XLONG(useT,0,0)
  end if
  res@KNOWNI  = 1.0
  res@KNOWNJ  = 1.0

  loc = wrf_ll_to_ij (longitude, latitude, res)
  loc = loc - 1
  
  if (dimsizes(dimsizes(loc)) .eq. 1) then
      loc!0 = "x_y"
  elseif (dimsizes(dimsizes(loc)) .eq. 2) then
      loc!0 = "x_y"
      loc!1 = "idx"
  else ; Not currently supported
      loc!0 = "x_y"
      loc!1 = "domain_idx"
      loc!2 = "idx"
  end if
  
  if ( returnI  ) then
    loci = new(dimsizes(loc),integer)
    ;loci@_FillValue = default_fillvalue("integer")   ; was -999
    loci = tointeger(loc + .5)
    loci!0 = loc!0
    return(loci)
  else
    return(loc)
  end if


end

;--------------------------------------------------------------------------------

undef("wrf_user_xy_to_ll")
function wrf_user_xy_to_ll( file_handle, x:numeric, y:numeric, \
                            opts_args:logical )

begin
;
; As of NCL V6.0.0, wrf_user_ll_to_ij can now handle a file 
; or a list of files.
;
  if(typeof(file_handle).eq."file") then
    ISFILE = True
    nc_file = file_handle
  elseif(typeof(file_handle).eq."list") then
    ISFILE = False
    nc_file = file_handle[0]
  else
    print("wrf_user_xy_to_ll: Error: the first argument must be a file or a list of files opened with addfile or addfiles")
    return
  end if

  opts = opts_args
  useT  = get_res_value(opts,"useTime",0)

  res = True
  res@MAP_PROJ  = nc_file@MAP_PROJ
  res@TRUELAT1  = nc_file@TRUELAT1
  res@TRUELAT2  = nc_file@TRUELAT2
  res@STAND_LON = nc_file@STAND_LON
  res@DX        = nc_file@DX
  res@DY        = nc_file@DY

  if (res@MAP_PROJ .eq. 6) then
    res@POLE_LAT  = nc_file@POLE_LAT
    res@POLE_LON  = nc_file@POLE_LON
    res@LATINC = (res@DY*360.)/2.0/3.141592653589793/6370000.
    res@LONINC = (res@DX*360.)/2.0/3.141592653589793/6370000.
  else
    res@POLE_LAT = 90.0
    res@POLE_LON =  0.0
    res@LATINC = 0.0
    res@LONINC = 0.0
  end if


  if(isfilevar(nc_file,"XLAT")) then
    if(ISFILE) then
      XLAT  = nc_file->XLAT(useT,:,:)
      XLONG = nc_file->XLONG(useT,:,:)
    else
      XLAT  = file_handle[:]->XLAT
      XLONG = file_handle[:]->XLONG
    end if
  else
    if(ISFILE) then
      XLAT  = nc_file->XLAT_M(useT,:,:)
      XLONG = nc_file->XLONG_M(useT,:,:)
    else
      XLAT  = file_handle[:]->XLAT_M
      XLONG = file_handle[:]->XLONG_M
    end if
  end if

  if(dimsizes(dimsizes(XLAT)).eq.2) then
; Rank 2
    res@REF_LAT = XLAT(0,0)
    res@REF_LON = XLONG(0,0)
  else
; Rank 3
    res@REF_LAT = XLAT(useT,0,0)
    res@REF_LON = XLONG(useT,0,0)
  end if
  res@KNOWNI  = 1.0
  res@KNOWNJ  = 1.0
    
  ; Convert to 1-based indexes for Fortran
  new_x = x + 1
  new_y = y + 1
  
  loc = wrf_ij_to_ll (new_x,new_y,res)
  
  if (dimsizes(dimsizes(loc)) .eq. 1) then
      loc!0 = "lon_lat"
  elseif (dimsizes(dimsizes(loc)) .eq. 2) then
      loc!0 = "lon_lat"
      loc!1 = "idx"
  else ; Not currently supported
      loc!0 = "lon_lat"
      loc!1 = "domain_idx"
      loc!2 = "idx"
  end if

  return(loc)


end

;--------------------------------------------------------------------------------
 
undef("wrf_user_vertcross")
function wrf_user_vertcross(var3d:numeric, z_in:numeric, \
                           loc_param:numeric, opts:logical )

; var3d      - 3d field to interpolate (all input fields must be unstaggered)
; z_in       - interpolate to this field (either p/z)
; loc_param  - an array of 4 values representing the start point and end point 
;              for the cross section (start_x, start_y, end_x, end_y) OR a single 
;              point when opt@use_pivot is True representing the pivot point. 
;              The values can be in grid coordinates or lat/lon coordinates 
;              (start_x = start_lon, start_y = start_lat, ...). If using 
;              lat/lon coordinates, then opt@latlon must be True.
; opts       - optional arguments
;   use_pivot   - set to True to indicate that loc_param and angle are used, 
;                 otherwise loc_param is set to 4 values to indicate a start and
;                 end point 
;   angle       - an angle for vertical plots - 90 represent a WE cross section, 
;                 ignored if use_pivot is False.
;   levels      - the vertical levels to use in the same units as z_in. Set to 
;                 False to automatically generate the number of levels specified 
;                 by autolevels.
;   latlon      - set to True if the values in loc_param are latitude and longitude
;                 values rather than grid values
;   file_handle - must be set to a file handle when latlon is True or 
;                 linecoords is True, otherwise this is ignored.
;   timeidx     - the time index to use for moving nests when latlon is True. Set
;                 to 0 if the nest is not moving.
;   linecoords  - set to True to include the latitude and longitude coordinates
;                 for the cross section line in the output attributes.
;   autolevels  - set to the desired number of levels when levels are 
;                 selected automatically (default 100).

begin

     use_pivot = get_res_value(opts, "use_pivot", False)
     angle = get_res_value(opts, "angle", 0.0)
     levels = get_res_value(opts, "levels", new(1,integer))
     latlon = get_res_value(opts, "latlon", False)
     file_handle = get_res_value(opts, "file_handle", 0)
     timeidx = get_res_value(opts, "timeidx", 0)
     linecoords = get_res_value(opts, "linecoords", False)
     nlevels = get_res_value(opts, "autolevels", 100)
     
     dims = dimsizes(var3d)
     nd = dimsizes(dims)
    
     dimX = dims(nd-1)
     dimY = dims(nd-2)
     dimZ = dims(nd-3)
     
     if ( nd .eq. 4 ) then
       z = z_in(0,:,:,:)
     else
       z = z_in
     end if

; Convert latlon to xy coordinates

    if (use_pivot) then 
        if (latlon) then
            opt = True
            opt@returnInt = True
            opt@useTime = timeidx
            ij := wrf_user_ll_to_xy(file_handle, loc_param(0), loc_param(1), opt)
            start_x = ij(0)
            start_y = ij(1)
        else
            start_x = loc_param(0)
            start_y = loc_param(1)
        end if
    else
        if (latlon) then
            opt = True
            opt@returnInt = True
            opt@useTime = timeidx
            ij := wrf_user_ll_to_xy(file_handle, (/ loc_param(0), loc_param(2) /), (/ loc_param(1), loc_param(3) /), opt)
            start_x = ij(0,0)
            start_y = ij(1,0)
            end_x = ij(0,1)
            end_y = ij(1,1)
        else
            start_x = loc_param(0)
            start_y = loc_param(1)
            end_x = loc_param(2)
            end_y = loc_param(3)
        end if
    end if 
 
; get the lat/lons along the cross section line if requested
     
; set the cross section line coordinates if requested
     if (linecoords) then 
     
         latname = "XLAT" 
         lonname = "XLONG" 
         if(.not. isfilevar(file_handle,"XLAT")) then
           if(isfilevar(file_handle,"XLAT_M")) then
             latname = "XLAT_M"
             lonname = "XLONG_M"
           end if
         end if
         
         latvar  = _get_wrf_var(file_handle, latname, timeidx)
         lonvar = _get_wrf_var(file_handle, lonname, timeidx)
         
        if (use_pivot) then
            loc := (/start_x, start_y/)
            linelats = wrf_user_intrp2d(latvar, loc, angle, False)
            linelons = wrf_user_intrp2d(lonvar, loc, angle, False)
        else
            loc := (/start_x, start_y, end_x, end_y /)
            linelats = wrf_user_intrp2d(latvar, loc, angle, True)
            linelons = wrf_user_intrp2d(lonvar, loc, angle, True)
        end if
        
     end if

; set vertical cross section
; Note for wrf_user_set_xy, opt is False when pivot and angle used.
     if (use_pivot) then
         xy = wrf_user_set_xy( z, start_x, start_y, \ ; assumes 0-based indexing in v6.5.0
                                0.0, 0.0, angle, False )
     
     else
         xy = wrf_user_set_xy( z, start_x, start_y, \    ; assumes 0-based indexing in v6.5.0
                                end_x, end_y, \
                                angle, True)
         
     end if
     xp = dimsizes(xy)
  
  
; first we interp z
     var2dz   = wrf_interp_2d_xy( z, xy)

;  interp to constant z grid
     if (all(ismissing(levels))) then
         if(var2dz(0,0) .gt. var2dz(1,0) ) then  ; monotonically decreasing coordinate
            z_max = floor(max(z)/10)*10     ; bottom value
            z_min = ceil(min(z)/10)*10      ; top value
            dz = (1.0/nlevels) * (z_max - z_min)
            ;nlevels = tointeger( (z_max-z_min)/dz)
            z_var2d = new( (/nlevels/), typeof(z))
            z_var2d(0) = z_max
            dz = -dz
         else
            z_max = max(z)
            z_min = 0.
            dz = (1.0/nlevels) * z_max
            ;nlevels = tointeger( z_max/dz )
            z_var2d = new( (/nlevels/), typeof(z))
            z_var2d(0) = z_min
         end if
      
         do i=1, nlevels-1
            z_var2d(i) = z_var2d(0)+i*dz
         end do
     else
         z_var2d = levels
         nlevels = dimsizes(z_var2d)
     end if
  
; interp the variable
     if ( dimsizes(dims) .eq. 4 ) then
       var2d = new( (/dims(0), nlevels, xp(0)/), typeof(var2dz))
       do it = 0,dims(0)-1
         var2dtmp = wrf_interp_2d_xy( var3d(it,:,:,:), xy)
         do i=0,xp(0)-1
            var2d(it,:,i) = wrf_interp_1d( var2dtmp(:,i), var2dz(:,i), z_var2d)
         end do
       end do
       var2d!0 = var3d!0
       var2d!1 = "vertical"
       var2d!2 = "cross_line_idx"
     else
       var2d = new( (/nlevels, xp(0)/), typeof(var2dz))
       var2dtmp = wrf_interp_2d_xy( var3d, xy)
       do i=0,xp(0)-1
          var2d(:,i) = wrf_interp_1d( var2dtmp(:,i), var2dz(:,i), z_var2d)
       end do
       var2d!0 = "vertical"
       var2d!1 = "cross_line_idx"
     end if

     st_x = tointeger(xy(0,0)) ; + 1 (removed 1-based indexing in 6.5.0) 
     st_y = tointeger(xy(0,1)) ; + 1
     ed_x = tointeger(xy(xp(0)-1,0)) ; + 1
     ed_y = tointeger(xy(xp(0)-1,1)) ; + 1
     if (.not. use_pivot) then
       var2d@Orientation = "Cross-Section: (" + \
                            st_x + "," + st_y + ") to (" + \
                            ed_x + "," + ed_y + ")"
     else
       var2d@Orientation = "Cross-Section: (" + \
                            st_x + "," + st_y + ") to (" + \
                            ed_x + "," + ed_y + ") ; center=(" + \
                            start_x + "," + start_y + \
                            ") ; angle=" + angle
     end if
     
     if (linecoords) then 
        var2d@lats = linelats
        var2d@lons = linelons
     end if
     
     var2d&vertical = z_var2d

     return(var2d)

end
