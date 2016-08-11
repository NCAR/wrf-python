;--------------------------------------------------------------------------------
; function wrf_user_set_xy( var:numeric, xp:numeric, yp:numeric, x1:numeric, \
;                           y1:numeric,angle:numeric ,opts )
; function wrf_user_intrp3d( var3d:numeric, z:numeric, plot_type:string, \
;                            loc_param:numeric, angle:numeric, opts:logical )
; function wrf_user_intrp2d( var2d:numeric, \
;                            loc_param:numeric, angle:numeric, opts:logical )
; function wrf_user_vert_interp(file_handle,field:float,\
;                               vert_coordinate[1]:string, \
        ;                       interp_levels[*]:numeric,opts[1]:logical)
; function wrf_user_getvar( file_handle, variable:string, time:integer )
; function wrf_user_list_times( file_handle )
; function wrf_user_ll_to_ij( nc_file:file, longitude:numeric, latitude:numeric, \
;                             opts:logical )
; function wrf_user_ij_to_ll( nc_file:file, i:numeric, j:numeric \
;                             opts:logical )
; function wrf_contour(nc_file:file,wks[1]: graphic, data[*][*]:numeric, \
;                      opt_args[1]:logical)
; function wrf_vector(nc_file:file,wks[1]: graphic, data_u[*][*]:numeric, \
;                     data_v[*][*]:numeric, opt_args[1]:logical)
; function wrf_map_resources(in_file[1]:file,opt_args[1]:logical)   
; function wrf_map(wks[1]:graphic,in_file[1]:file,opt_args[1]:logical)   
; function wrf_map_overlays(in_file[1]:file,wks:graphic,base[1]:graphic, \
;                           plots[*]:graphic,opt_arg[1]:logical,mp_arg[1]:logical)
; function wrf_overlays(in_file[1]:file,wks:graphic, plots[*]:graphic, \
;                       opt_arg[1]:logical)
; function wrf_user_unstagger( varin:numeric, unstagDim:string )
;
; function wrf_cloud_fraction(file_handle,out_arrays[*][*][*]:numeric,time:integer)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Thse functions are still experimental
; function wrf_wps_dom(wks[1]:graphic,mpres[1]:logical,lnres[1]:logical,txres[1]:logical)   
; function wrf_contour_ps(nc_file:file,wks[1]: graphic, data[*][*]:numeric, \
;                         opt_args[1]:logical)
; function wrf_vector_ps(nc_file:file,wks[1]: graphic, \
;                        data_u[*][*]:numeric, data_v[*][*]:numeric, \
;                        opt_args[1]:logical)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; These functions/procedures are obsolete as of version 5.0.1 of NCL.
; Do not use them. 
; Use wrf_overlays instead of wrf_overlay
; Use wrf_map_overlays instead of wrf_map_overlay
; Use wrf_user_ll_to_ij instead of wrf_user_latlon_to_ij
;
; function wrf_map_zoom(wks[1]:graphic,in_file[1]:file,opt_args[1]:logical, \
;                       y1:integer,y2:integer,x1:integer,x2:integer)
; procedure wrf_map_overlay(wks:graphic,base[1]:graphic, \
;                           plots[*]:graphic, \
;                           opt_arg[1]:logical)
; procedure wrf_overlay(wks:graphic, plots[*]:graphic, \
;                           opt_arg[1]:logical)
; function wrf_user_latlon_to_ij( nc_file:file, latitude:numeric, 
;                                 longitude:numeric )
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; function add_white_space(str:string,maxlen:integer)
; procedure print_opts(opts_name,opts,debug)
; procedure print_header(icount:integer,debug)
;
;--------------------------------------------------------------------------------
 external  fortran_clouds "./cloudsf.so"
;--------------------------------------------------------------------------------
undef("wrf_user_set_xy")
function wrf_user_set_xy( var:numeric, xp:numeric, yp:numeric, x1:numeric, \
                          y1:numeric, angle:numeric, opts )

;  mass coordinate version of ncl user routines

local dims,x,y,slope,intercept,distance,dx,dy,dxy,npts,xy

begin

; find intersection of line and domain boundaries

  dims = dimsizes(var)

  if (.not. opts) then        ; We have a pivot point and location and 
			      ; need to calculate the start and end points of 
			      ; the cross section

     if ((angle .gt. 315.) .or. (angle .lt. 45.) .or. \
         ((angle .gt. 135.) .and. (angle .lt. 225.)) )  then

;  x = y*slope + intercept

          slope = -(360.-angle)/45.
          if( angle .lt. 45. ) then
            slope = angle/45.
          end if
          if( angle .gt. 135.) then
            slope = (angle-180.)/45.
          end if
          intercept = xp - yp*slope

; find intersections with domain boundaries

          y0 = 0.
          x0 = y0*slope + intercept
   
          if( x0 .lt. 0.) then  ; intersect outside of left boundary
            x0 = 0.
            y0 =  (x0 - intercept)/slope
          end if
          if( x0 .gt. dims(2)-1) then  ; intersect outside of right boundary
            x0 = dims(2)-1
            y0 =  (x0 - intercept)/slope
          end if
   
          y1 = dims(1)-1.  ;  need to make sure this will be a float?
          x1 = y1*slope + intercept
   
          if( x1 .lt. 0.) then  ; intersect outside of left boundary
            x1 = 0.
            y1 =  (x1 - intercept)/slope
          end if
          if( x1 .gt. dims(2)-1) then  ; intersect outside of right boundary
            x1 = dims(2)-1
            y1 =  (x1 - intercept)/slope
          end if

     else

;  y = x*slope + intercept

          slope = (90.-angle)/45.
          if( angle .gt. 225. ) then
            slope = (270.-angle)/45.
          end if
          intercept = yp - xp*slope

; find intersections with domain boundaries

          x0 = 0.
          y0 = x0*slope + intercept
   
          if( y0 .lt. 0.) then  ; intersect outside of bottom boundary
            y0 = 0.
            x0 =  (y0 - intercept)/slope
          end if
          if( y0 .gt. dims(1)-1) then  ; intersect outside of top boundary
            y0 = dims(1)-1
            x0 =  (y0 - intercept)/slope
          end if

          x1 = dims(2)-1.  ;  need to make sure this will be a float?
          y1 = x1*slope + intercept
   
          if( y1 .lt. 0.) then  ; intersect outside of bottom boundary
            y1 = 0.
            x1 =  (y1 - intercept)/slope
          end if
          if( y1 .gt. dims(1)-1) then  ; intersect outside of top boundary
            y1 = dims(1)-1
            x1 =  (y1 - intercept)/slope
          end if
   
     end if  ; we have beginning and ending points

  end if


  if (opts) then   	; We have a specified start and end point
    x0 = xp
    y0 = yp
    if ( x1 .gt. dims(2)-1 ) then 
      x1 = dims(2)
    end if
    if ( y1 .gt. dims(1)-1 ) then 
      y1 = dims(1)
    end if
  end if

  dx = x1 - x0
  dy = y1 - y0
  distance = (dx*dx + dy*dy)^0.5
  npts = tointeger(distance)
  dxy = new(1,typeof(distance))
  dxy = distance/npts

  xy = new((/ npts, 2 /),typeof(x1))

  dx = dx/npts
  dy = dy/npts

  do i=0,npts-1
    xy(i,0) = x0 + i*dx
    xy(i,1) = y0 + i*dy
  end do

; print(xy)
  return(xy)

end

;--------------------------------------------------------------------------------
 
undef("wrf_user_intrp3d")
function wrf_user_intrp3d( var3d:numeric, z_in:numeric, \
                           plot_type:string, \
                           loc_param:numeric, angle:numeric, opts:logical )

; var3d      - 3d field to interpolate (all input fields must be unstaggered)
; z_in       - interpolate to this field (either p/z)
; plot_type  - interpolate horizontally "h", or vertically "v"
; loc_param  - level(s) for horizontal plots (eg. 500hPa ; 3000m - scalar), 
;              plane for vertical plots (2 values representing an xy point 
;              on the model domain through which the vertical plane will pass
;              OR 4 values specifying start and end values
; angle      - 0.0 for horizontal plots, and 
;              an angle for vertical plots - 90 represent a WE cross section
; opts         Used IF opts is TRUE, else use loc_param and angle to determine crosssection

begin

  if(plot_type .eq. "h" ) then   ;  horizontal cross section needed

     dimL = dimsizes(loc_param)

     dims = dimsizes(var3d)
     nd = dimsizes(dims)

     dimX = dims(nd-1)
     dimY = dims(nd-2)
     dimZ = dims(nd-3)
     dim4 = 1
     dim5 = 1
     if ( nd .eq. 4 ) then
       dim4 = dims(nd-4)
     end if
     if ( nd .eq. 5 ) then
       dim4 = dims(nd-4)
       dim5 = dims(nd-5)
     end if

     var3  = new ( (/ dim5, dim4, dimZ, dimY, dimX /) , typeof(var3d) )
     z     = new ( (/ dim5, dim4, dimZ, dimY, dimX /) , typeof(var3d) )
     var2d = new ( (/ dim5, dim4, dimL, dimY, dimX /) , typeof(var3d) )

     if ( nd .eq. 5 ) then
       var3 = var3d
          z =  z_in
     end if
     if ( nd .eq. 4 ) then
       var3(0,:,:,:,:) = var3d(:,:,:,:)
          z(0,:,:,:,:) =  z_in(:,:,:,:)
     end if
     if ( nd .eq. 3 ) then
       var3(0,0,:,:,:) = var3d(:,:,:)
          z(0,0,:,:,:) =  z_in(:,:,:)
     end if


     if ( z(0,0,dimZ-1,0,0) .lt. z(0,0,dimZ-2,0,0) ) then
       ; We must be interpolating to pressure
       ; This routine needs input field and level in hPa - lets make sure of this
       if ( max(z) .gt. 2000. ) then
         ; looks like we have Pa as input - make this hPa
         z = z * 0.01
       end if
       if ( loc_param(0) .gt. 2000. ) then
         ; looks like the input was specified in Pa - change this
         loc_param = loc_param * 0.01
       end if
     end if

     do il = 0,dimL-1
       var = wrf_interp_3d_z(var3,z,loc_param(il))
       var2d(:,:,il,:,:) = var(:,:,:,:)
     end do

     copy_VarAtts(var3d,var3)
     if(isatt(var3,"description")) then
       delete_VarAtts(var3,(/"description"/))
     end if
     if(isatt(var3,"units")) then
       delete_VarAtts(var3,(/"units"/))
     end if
     if(isatt(var3,"MemoryOrder")) then
       delete_VarAtts(var3,(/"MemoryOrder"/))
     end if
     if(isatt(var3,"_FillValue")) then
       delete_VarAtts(var3,(/"_FillValue"/))
     end if
     copy_VarAtts(var3,var2d)

     nn = nd-2
     var2d!nn = "plevs"

     if ( dimL .gt. 1 ) then
       if ( nd .eq. 5 ) then
         return( var2d )
       end if
       if ( nd .eq. 4 ) then
         return( var2d(0,:,:,:,:) )
       end if
       if ( nd .eq. 3 ) then
         return( var2d(0,0,:,:,:) )
       end if
     else
       if ( z(0,0,dimZ-1,0,0) .lt. z(0,0,dimZ-2,0,0) ) then
          var2d@PlotLevelID = loc_param + " hPa"
       else
          var2d@PlotLevelID = .001*loc_param + " km"
       end if
       if ( nd .eq. 5 ) then
         return( var2d(:,:,0,:,:) )
       end if
       if ( nd .eq. 4 ) then
         return( var2d(0,:,0,:,:) )
       end if
       if ( nd .eq. 3 ) then
         return( var2d(0,0,0,:,:) )
       end if
     end if


  end if




  if(plot_type .eq. "v" ) then   ;  vertical cross section needed

     dims = dimsizes(var3d)
     nd = dimsizes(dims)

     dimX = dims(nd-1)
     dimY = dims(nd-2)
     dimZ = dims(nd-3)

     if ( nd .eq. 4 ) then
       if ( z_in(0,dimZ-1,0,0) .lt. z_in(0,dimZ-2,0,0) ) then
         ; We must be interpolating to pressure
           ; This routine needs input field and level in hPa - lets make sure of this
         if ( max(z_in) .gt. 2000. ) then
           ; looks like we have Pa as input - make this hPa
           z_in = z_in * 0.01
         end if
       end if
       z = z_in(0,:,:,:)
     else
       if ( z_in(dimZ-1,0,0) .lt. z_in(dimZ-2,0,0) ) then
         ; We must be interpolating to pressure
           ; This routine needs input field and level in hPa - lets make sure of this
         if ( z_in(0,0,0) .gt. 2000. ) then
           ; looks like we have Pa as input - make this hPa
           z_in = z_in * 0.01
         end if
       end if
       z = z_in
     end if

       
; set vertical cross section
     if (opts) then
       xy = wrf_user_set_xy( z, loc_param(0)-1, loc_param(1)-1, \    ; the -1 is for NCL dimensions
                                loc_param(2)-1, loc_param(3)-1, \
                                angle, opts )
     else
       xy = wrf_user_set_xy( z, loc_param(0), loc_param(1), \
                                0.0, 0.0, angle, opts )
     end if
     xp = dimsizes(xy)
  
  
; first we interp z
     var2dz   = wrf_interp_2d_xy( z, xy)

;  interp to constant z grid
     if(var2dz(0,0) .gt. var2dz(1,0) ) then  ; monotonically decreasing coordinate
        z_max = floor(max(z)/10)*10     ; bottom value
        z_min = ceil(min(z)/10)*10      ; top value
        dz = 10
        nlevels = tointeger( (z_max-z_min)/dz)
        z_var2d = new( (/nlevels/), typeof(z))
        z_var2d(0) = z_max
        dz = -dz
     else
        z_max = max(z)
        z_min = 0.
        dz = 0.01 * z_max
        nlevels = tointeger( z_max/dz )
        z_var2d = new( (/nlevels/), typeof(z))
        z_var2d(0) = z_min
     end if
  
     do i=1, nlevels-1
        z_var2d(i) = z_var2d(0)+i*dz
     end do
  

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
       var2d!1 = "Vertical"
       var2d!2 = "Horizontal"
     else
       var2d = new( (/nlevels, xp(0)/), typeof(var2dz))
       var2dtmp = wrf_interp_2d_xy( var3d, xy)
       do i=0,xp(0)-1
          var2d(:,i) = wrf_interp_1d( var2dtmp(:,i), var2dz(:,i), z_var2d)
       end do
       var2d!0 = "Vertical"
       var2d!1 = "Horizontal"
     end if



     st_x = tointeger(xy(0,0)) + 1
     st_y = tointeger(xy(0,1)) + 1
     ed_x = tointeger(xy(xp(0)-1,0)) + 1
     ed_y = tointeger(xy(xp(0)-1,1)) + 1
     if (opts) then
       var2d@Orientation = "Cross-Section: (" + \
                            st_x + "," + st_y + ") to (" + \
                            ed_x + "," + ed_y + ")"
     else
       var2d@Orientation = "Cross-Section: (" + \
                            st_x + "," + st_y + ") to (" + \
                            ed_x + "," + ed_y + ") ; center=(" + \
                            loc_param(0) + "," + loc_param(1) + \
                            ") ; angle=" + angle
     end if

     return(var2d)
  end if



end

;--------------------------------------------------------------------------------
 
undef("wrf_user_intrp2d")
function wrf_user_intrp2d( var2d:numeric, \
                           loc_param:numeric, angle:numeric, opts:logical )

; var2d      - 2d field to interpolate 
; loc_param  - plane for vertical plots (2 values representing an xy point 
;              on the model domain through which the vertical plane will pass
;              OR 4 values specifying start and end values
; angle      - 0.0 for horizontal plots, and 
;              an angle for vertical plots - 90 represent a WE cross section
; opts         Used IF opts is TRUE, else use loc_param and angle to determine crosssection

begin

     dims = dimsizes(var2d)
     nd = dimsizes(dims)
     dimX = dims(nd-1)
     dimY = dims(nd-2)
     dimT = 1
     if ( nd .eq. 3 ) then
       dimT = dims(nd-3)
     end if
     var2dtmp = new( (/ 1, dimY, dimX /), typeof(var2d) )
     

; set vertical cross section
     if ( nd .eq. 3 ) then
       var2dtmp(0,:,:) = var2d(0,:,:)
     else
       var2dtmp(0,:,:) = var2d(:,:)
     end if
     if (opts) then
       xy = wrf_user_set_xy( var2dtmp, \
                             loc_param(0)-1, loc_param(1)-1, \    ; the -1 is for NCL dimensions
                             loc_param(2)-1, loc_param(3)-1, \
                             angle, opts )
     else
       xy = wrf_user_set_xy( var2dtmp, \
                             loc_param(0), loc_param(1), \
                             0.0, 0.0, angle, opts )
     end if
     xp = dimsizes(xy)

     var2dout = new( (/ dimT, xp(0) /), typeof(var2d) )

     var1dtmp = wrf_interp_2d_xy( var2dtmp, xy )
     var2dout(0,:) = var1dtmp(0,:)
     if ( dimT .eq. 1 ) then
       var2dout!1 = "Horizontal"
       return ( var2dout(0,:) )
     end if

     do it = 1,dimT-1
       var2dtmp(0,:,:) = var2d(it,:,:)
       var1dtmp = wrf_interp_2d_xy( var2dtmp, xy )
       var2dout(it,:) = var1dtmp(0,:)
     end do
     var2dout!0 = "Time"
     var2dout!1 = "Horizontal"
     return(var2dout)

end

;--------------------------------------------------------------------------------

undef("wrf_user_unstagger")
function wrf_user_unstagger( varin:numeric, unstagDim:string )

begin

  dims = dimsizes(varin)
  nd = dimsizes(dims)

  if ( unstagDim .eq. "X" .or. unstagDim .eq. "U" ) then
    dimU = dims(nd-1)
    if ( nd .eq. 5 ) then
      varout = 0.5*(varin(:,:,:,:,:dimU-2) + varin(:,:,:,:,1:dimU-1))
    end if
    if ( nd .eq. 4 ) then
      varout = 0.5*(varin(:,:,:,:dimU-2) + varin(:,:,:,1:dimU-1))
    end if
    if ( nd .eq. 3 ) then
      varout = 0.5*(varin(:,:,:dimU-2) + varin(:,:,1:dimU-1))
    end if
    if ( nd .eq. 2 ) then
      varout = 0.5*(varin(:,:dimU-2) + varin(:,1:dimU-1))
    end if
    do i = 0,nd-2
      varout!i = varin!i
    end do
    i = nd-1
    varout!i = "west_east"
    copy_VarAtts(varin,varout)
    varout@coordinates = "XLONG XLAT"
    varout@stagger = " "
  end if

  if ( unstagDim .eq. "Y" .or. unstagDim .eq. "V" ) then
    dimV = dims(nd-2)
    if ( nd .eq. 5 ) then
      varout = 0.5*(varin(:,:,:,:dimV-2,:)+varin(:,:,:,1:dimV-1,:))
    end if
    if ( nd .eq. 4 ) then
      varout = 0.5*(varin(:,:,:dimV-2,:)+varin(:,:,1:dimV-1,:))
    end if
    if ( nd .eq. 3 ) then
      varout = 0.5*(varin(:,:dimV-2,:)+varin(:,1:dimV-1,:))
    end if
    if ( nd .eq. 2 ) then
      varout = 0.5*(varin(:dimV-2,:)+varin(1:dimV-1,:))
    end if
    do i = 0,nd-1
      varout!i = varin!i
    end do
    i = nd-2
    varout!i = "south_north"
    copy_VarAtts(varin,varout)
    varout@coordinates = "XLONG XLAT"
    varout@stagger = " "
  end if

  if ( unstagDim .eq. "Z" ) then
    dimW = dims(nd-3)
    if ( nd .eq. 5 ) then
      varout = 0.5*(varin(:,:,0:dimW-2,:,:)+varin(:,:,1:dimW-1,:,:))
    end if
    if ( nd .eq. 4 ) then
      varout = 0.5*(varin(:,0:dimW-2,:,:)+varin(:,1:dimW-1,:,:))
    end if
    if ( nd .eq. 3 ) then
      varout = 0.5*(varin(0:dimW-2,:,:)+varin(1:dimW-1,:,:))
    end if
    do i = 0,nd-1
      varout!i = varin!i
    end do
    i = nd-3
    varout!i = "bottom_top" 
    copy_VarAtts(varin,varout)
    varout@coordinates = "XLONG XLAT"
    varout@stagger = " "
  end if

  if( any( unstagDim .eq. (/"X","U","Y","V","Z"/) ) ) then
    return(varout) 
  else
    print("NOTE: No unstaggering required, as the input field is already on mass points.")
    return(varin)
  end if

end

;--------------------------------------------------------------------------------
; This function was modified in May 2011 to allow a list of files.
;
undef("wrf_user_getvar")
function wrf_user_getvar( file_handle, varin[*]:string, time_in:integer )
local variable, time, var, u, v, u_in, v_in, pii, radians_per_degree, \
dims, nd, latitude, longitude, rank
begin

;---As of NCL V6.0.0, wrf_user_getvar can now handle a file or a list of files.
  if(typeof(file_handle).eq."file") then
    ISFILE = True
    nc_file = file_handle
  else if(typeof(file_handle).eq."list") then
    ISFILE = False
    nc_file = file_handle[0]
  else
    print("wrf_user_getvar: error: the first argument must be a file or a list of files opened with addfile or addfiles")
    return
  end if
  end if

  variable = varin(0)
  time = time_in(0)

  if( (variable .eq. "uvmet") .or. (variable .eq. "uvmet10") ) then
       ;; Calculate winds rotated to earth coord. 

       pii = 3.14159265
       radians_per_degree = pii/180.

       if( (variable .eq. "uvmet") )  then
         getU = "U" 
         getV = "V" 
         if(.not. isfilevar(nc_file,"U")) then
           if(isfilevar(nc_file,"UU")) then
             getU = "UU"
             getV = "VV"
           end if
         end if
         if ( time .eq. -1 ) then
           if(ISFILE) then
             u_in = nc_file->$getU$
             v_in = nc_file->$getV$
           else
             u_in = file_handle[:]->$getU$
             v_in = file_handle[:]->$getV$
           end if
           u = wrf_user_unstagger(u_in,u_in@stagger)
           v = wrf_user_unstagger(v_in,v_in@stagger)
       else
           if(ISFILE) then
             u_in   = nc_file->$getU$(time_in,:,:,:)
             v_in   = nc_file->$getV$(time_in,:,:,:)
           else
             u_in   = file_handle[:]->$getU$(time_in,:,:,:)
             v_in   = file_handle[:]->$getV$(time_in,:,:,:)
           end if
           u = wrf_user_unstagger(u_in,u_in@stagger)
           v = wrf_user_unstagger(v_in,v_in@stagger)
         end if
       end if  

       if( (variable .eq. "uvmet10") )  then
         if(isfilevar(nc_file,"U10")) then
           if ( time .eq. -1 ) then
             if(ISFILE) then
               u_in   = nc_file->U10
               v_in   = nc_file->V10
             else
               u_in   = file_handle[:]->U10
               v_in   = file_handle[:]->V10
             end if
             u = wrf_user_unstagger(u_in,u_in@stagger)
             v = wrf_user_unstagger(v_in,v_in@stagger)
           else
             if(ISFILE) then
               u_in  = nc_file->U10(time_in,:,:)
               v_in  = nc_file->V10(time_in,:,:)
             else
               u_in  = file_handle[:]->U10(time_in,:,:)
               v_in  = file_handle[:]->V10(time_in,:,:)
             end if
             u = wrf_user_unstagger(u_in,u_in@stagger)
             v = wrf_user_unstagger(v_in,v_in@stagger)
           end if
         else  ; may be a met file, so get lowest level of UU and VV
           if(isfilevar(nc_file,"UU")) then
             print("wrf_user_getvar: Assume this is a met_em file - getting lowest level from UU and VV fields")
             if ( time .eq. -1 ) then
               if(ISFILE) then
                 u_in  = nc_file->UU(:,0,:,:)
                 v_in  = nc_file->VV(:,0,:,:)
               else
                 u_in  = file_handle[:]->UU(:,0,:,:)
                 v_in  = file_handle[:]->VV(:,0,:,:)
               end if
               u = wrf_user_unstagger(u_in,u_in@stagger)
               v = wrf_user_unstagger(v_in,v_in@stagger)
             else
               if(ISFILE) then
                 u_in  = nc_file->UU(time_in,0,:,:)
                 v_in  = nc_file->VV(time_in,0,:,:)
               else
                 u_in  = file_handle[:]->UU(time_in,0,:,:)
                 v_in  = file_handle[:]->VV(time_in,0,:,:)
               end if
               u = wrf_user_unstagger(u_in,u_in@stagger)
               v = wrf_user_unstagger(v_in,v_in@stagger)
             end if
           end if
         end if  
       end if  


       map_projection = nc_file@MAP_PROJ


       if( any(map_projection.eq.(/0,3,6/)) ) then  ; no rotation needed
         dims = dimsizes(u)
         nd = dimsizes(dims)
         if ( nd .eq. 5 ) then
           uvmet = new( (/ 2, dims(0), dims(1), dims(2), dims(3), dims(4) /), typeof(u))
           uvmet(0,:,:,:,:,:) = u(:,:,:,:,:) 
           uvmet(1,:,:,:,:,:) = v(:,:,:,:,:)
         end if
         if ( nd .eq. 4 ) then
           uvmet = new( (/ 2, dims(0), dims(1), dims(2), dims(3) /), typeof(u))
           uvmet(0,:,:,:,:) = u(:,:,:,:) 
           uvmet(1,:,:,:,:) = v(:,:,:,:)
         end if
         if ( nd .eq. 3 ) then
           uvmet = new( (/ 2, dims(0), dims(1), dims(2) /), typeof(u))
           uvmet(0,:,:,:) = u(:,:,:) 
           uvmet(1,:,:,:) = v(:,:,:)
         end if
         if ( nd .eq. 2 ) then
           uvmet = new( (/ 2, dims(0), dims(1) /), typeof(u))
           uvmet(0,:,:) = u(:,:) 
           uvmet(1,:,:) = v(:,:)
         end if
         delete_VarAtts(u,(/"description","units"/))
         copy_VarAtts(u,uvmet)
         uvmet@description = " u,v met velocity"
         uvmet!0 = "u_v"
       end if


       if( any(map_projection.eq.(/1,2/)) ) then  ; no rotation needed
         cen_lat  = nc_file@CEN_LAT
         if(isatt(nc_file,"STAND_LON")) then
             cen_long = nc_file@STAND_LON
         else
             cen_long = nc_file@CEN_LON
         end if
         true_lat1 = nc_file@TRUELAT1
         true_lat2 = nc_file@TRUELAT2
  
         getLAT = "XLAT" 
         getLON = "XLONG" 
         if(.not. isfilevar(nc_file,"XLAT")) then
           if(isfilevar(nc_file,"XLAT_M")) then
             getLAT = "XLAT_M"
             getLON = "XLONG_M"
           end if
         end if
         if ( time .eq. -1 ) then
           if(ISFILE) then
             latitude  = nc_file->$getLAT$
             longitude = nc_file->$getLON$
           else
             latitude  = file_handle[:]->$getLAT$
             longitude = file_handle[:]->$getLON$
           end if
         else
           if(ISFILE) then
             latitude  = nc_file->$getLAT$(time_in,:,:)
             longitude = nc_file->$getLON$(time_in,:,:)
           else
             latitude  = file_handle[:]->$getLAT$(time_in,:,:)
             longitude = file_handle[:]->$getLON$(time_in,:,:)
           end if
         end if

         cone = 1.
         if( map_projection .eq. 1) then    ; Lambert Conformal mapping
           if( (fabs(true_lat1 - true_lat2) .gt. 0.1) .and.  \
               (fabs(true_lat2 - 90. )      .gt. 0.1)       ) then
               cone = log(cos(true_lat1*radians_per_degree)) \ 
                     -log(cos(true_lat2*radians_per_degree)) 
               cone = cone/( log(tan(( 45. -fabs(true_lat1/2.))*radians_per_degree)) - \ 
                             log(tan((45. -fabs(true_lat2/2.))*radians_per_degree)) ) 
           else
               cone = sin(fabs(true_lat1)*radians_per_degree)
           end if
         end if
         if(map_projection .eq. 2) then      ; polar stereographic
           cone = 1.
         end if
         if(map_projection .eq. 3) then      ; Mercator
           cone = 0.
         end if

         uvmet = wrf_uvmet( u, v, latitude, longitude, cen_long, cone )
         delete_VarAtts(u,(/"description","units"/))
         copy_VarAtts(u,uvmet)

       end if  
       
       if( (variable .eq. "uvmet10") )  then
         uvmet@description = " u10,v10 met velocity"
       end if
       return(uvmet)

  end if



  if( variable .eq. "ua" ) then
       ;  U interpolated to mass points
       getTHIS = "U" 
       if(.not. isfilevar(nc_file,"U")) then
         if(isfilevar(nc_file,"UU")) then
           getTHIS = "UU"
         end if
       end if
       if ( time .eq. -1 ) then
         if(ISFILE) then
           var = nc_file->$getTHIS$
         else
           var = file_handle[:]->$getTHIS$
         end if
       else
         if(ISFILE) then
           var = nc_file->$getTHIS$(time_in,:,:,:)
         else
           var = file_handle[:]->$getTHIS$(time_in,:,:,:)
         end if
       end if

       ua = wrf_user_unstagger(var,var@stagger)

       return(ua)
  end if



  if( variable .eq. "va" ) then
       ; V interpolated to mass points
       getTHIS = "V" 
       if(.not. isfilevar(nc_file,"V")) then
         if(isfilevar(nc_file,"VV")) then
           getTHIS = "VV"
         end if
       end if
       if ( time .eq. -1 ) then
         if(ISFILE) then
           var = nc_file->$getTHIS$
         else
           var = file_handle[:]->$getTHIS$
         end if
       else
         if(ISFILE) then
           var = nc_file->$getTHIS$(time_in,:,:,:)
         else
           var = file_handle[:]->$getTHIS$(time_in,:,:,:)
         end if
       end if

       va = wrf_user_unstagger(var,var@stagger)

       return(va)
  end if



  if( variable .eq. "wa" ) then
       ; W interpolated to mass points
       if ( time .eq. -1 ) then
         if(ISFILE) then
           var = nc_file->W
         else
           var = file_handle[:]->W
         end if
       else
         if(ISFILE) then
           var = nc_file->W(time_in,:,:,:)
         else
           var = file_handle[:]->W(time_in,:,:,:)
         end if
       end if

       wa = wrf_user_unstagger(var,var@stagger)

       return(wa)
  end if



  if( any( variable .eq. (/"p","pres","pressure"/) ) ) then
       ; Full model pressure [=base pressure (PB) + pertubation pressure (P)]
       if(isfilevar(nc_file,"P")) then
         if ( time .eq. -1 ) then
           if(ISFILE) then
             var = nc_file->P
             PB  = nc_file->PB
           else
             var = file_handle[:]->P
             PB  = file_handle[:]->PB
           end if
         else
           if(ISFILE) then
             var = nc_file->P(time_in,:,:,:)
             PB  = nc_file->PB(time_in,:,:,:)
           else
             var = file_handle[:]->P(time_in,:,:,:)
             PB  = file_handle[:]->PB(time_in,:,:,:)
           end if
         end if
         var = var + PB
       else
         ;; may be a met_em file - see if we can get PRES
         if(isfilevar(nc_file,"PRES")) then
           if ( time .eq. -1 ) then
             if(ISFILE) then
               var = nc_file->PRES
             else
               var = file_handle[:]->PRES
             end if
           else
             if(ISFILE) then
               var = nc_file->PRES(time_in,:,:,:)
             else
               var = file_handle[:]->PRES(time_in,:,:,:)
             end if
           end if
         end if
       end if
       var@description = "Pressure"
       if( variable .eq. "pressure" ) then
         var = var * 0.01
         var@units = "hPa"
       end if
       return(var)
  end if


  if( any( variable .eq. (/"geopt","geopotential","z","height"/) ) ) then
       ; Height [=full geopotentail height / 9.81]
       if(isfilevar(nc_file,"PH")) then
         if ( time .eq. -1 ) then
           if(ISFILE) then
             var = nc_file->PH
             PHB = nc_file->PHB
           else
             var = file_handle[:]->PH
             PHB = file_handle[:]->PHB
           end if
         else
           if(ISFILE) then
             var = nc_file->PH(time,:,:,:)
             PHB = nc_file->PHB(time,:,:,:)
           else
             var = file_handle[:]->PH(time,:,:,:)
             PHB = file_handle[:]->PHB(time,:,:,:)
           end if
         end if

         var = var + PHB
         z = wrf_user_unstagger(var,var@stagger)
         z@description = "Geopotential"

       else
         ;; may be a met_em file - see if we can get GHT - data in met_em file is Height in M
         if(isfilevar(nc_file,"GHT")) then
           if ( time .eq. -1 ) then
             if(ISFILE) then
               z = nc_file->GHT
             else
               z = file_handle[:]->GHT
             end if
           else
             if(ISFILE) then
               z = nc_file->GHT(time,:,:,:)
             else
               z = file_handle[:]->GHT(time,:,:,:)
             end if
           end if
           z = z * 9.81
           z@description = "Geopotential"
           z@units       = "m2 s-2"
         end if
       end if

       if( any( variable .eq. (/"z","height"/) ) ) then
         z = z / 9.81
         z@description = "Height"
         z@units       = "m"
       end if
       return(z)
  end if



  if( any( variable .eq. (/"th","theta"/) ) ) then
       ; Potentail Temperature is model output T + 300K
       if ( time .eq. -1 ) then
         if(ISFILE) then
           var = nc_file->T
         else
           var = file_handle[:]->T
         end if
       else
         if(ISFILE) then
           var = nc_file->T(time_in,:,:,:)
         else
           var = file_handle[:]->T(time_in,:,:,:)
         end if
       end if  
       var = var + 300.
       var@description = "Potential Temperature (theta) "
       return(var)
  end if



  if( any( variable .eq. (/"tk","tc"/) ) ) then
       ;; function wrf_tk needs theta and pressure (Pa) on input and returns temperature in K on return
       if(isfilevar(nc_file,"T")) then
         if ( time .eq. -1 ) then
           if(ISFILE) then
             T  = nc_file->T
             P  = nc_file->P
             PB = nc_file->PB
           else
             T  = file_handle[:]->T
             P  = file_handle[:]->P
             PB = file_handle[:]->PB
           end if
         else
           if(ISFILE) then
             T  = nc_file->T(time_in,:,:,:)
             P  = nc_file->P(time_in,:,:,:)
             PB = nc_file->PB(time_in,:,:,:)
           else
             T  = file_handle[:]->T(time_in,:,:,:)
             P  = file_handle[:]->P(time_in,:,:,:)
             PB = file_handle[:]->PB(time_in,:,:,:)
           end if
         end if  
         T = T + 300.
         P = P + PB
         t = wrf_tk( P , T )
         delete_VarAtts(T,(/"description"/))
         copy_VarAtts(T,t)
       else
         ;; may be a met_em file - see if we can get TT
         if(isfilevar(nc_file,"TT")) then
           if ( time .eq. -1 ) then
             if(ISFILE) then
               t = nc_file->TT
             else
               t = file_handle[:]->TT
             end if
           else
             if(ISFILE) then
               t = nc_file->TT(time_in,:,:,:)
             else
               t = file_handle[:]->TT(time_in,:,:,:)
             end if
           end if  
         end if  
       end if  
       if( variable .eq. "tc" ) then
         t = t - 273.16
         t@units = "C"     ; Overwrite return units
       end if
       return(t)
  end if



   if( variable .eq. "eth" ) then
   ;Equivalent Potential Temperature in degrees K
       if(isfilevar(nc_file,"T")) then
         if ( time .eq. -1 ) then
           if(ISFILE) then
             T  = nc_file->T
             P  = nc_file->P
             PB = nc_file->PB
             QV = nc_file->QVAPOR
           else
             T  = file_handle[:]->T
             P  = file_handle[:]->P
             PB = file_handle[:]->PB
             QV = file_handle[:]->QVAPOR
           end if
         else
           if(ISFILE) then
             T  = nc_file->T(time_in,:,:,:)
             P  = nc_file->P(time_in,:,:,:)
             PB = nc_file->PB(time_in,:,:,:)
             QV = nc_file->QVAPOR(time_in,:,:,:)
           else
             T  = file_handle[:]->T(time_in,:,:,:)
             P  = file_handle[:]->P(time_in,:,:,:)
             PB = file_handle[:]->PB(time_in,:,:,:)
             QV = file_handle[:]->QVAPOR(time_in,:,:,:)
           end if
         end if  
         T  = T + 300.        ; potential temperature in K.
         P  = P + PB          ; full pressure in Pa.   
         tk = wrf_tk( P , T ) ; temperature in K.
         eth = wrf_eth ( QV, tk, P )
         delete_VarAtts(T,(/"description"/))
         copy_VarAtts(T,eth)
         return(eth)
       else
         print("This diagnostic only for with WRF data file and not for WPS files")
         exit
       end if  

   end if



  if( variable .eq. "td" ) then
       ;; function wrf_td needs qv and pressure (Pa) on input and returns dewpoint temperature on return
       if ( time .eq. -1 ) then
         if(ISFILE) then
           P      = nc_file->P
           PB     = nc_file->PB
           QVAPOR = nc_file->QVAPOR
         else
           P      = file_handle[:]->P
           PB     = file_handle[:]->PB
           QVAPOR = file_handle[:]->QVAPOR
         end if
       else
         if(ISFILE) then
           P      = nc_file->P(time_in,:,:,:)
           PB     = nc_file->PB(time_in,:,:,:)
           QVAPOR = nc_file->QVAPOR(time_in,:,:,:)
         else
           P      = file_handle[:]->P(time_in,:,:,:)
           PB     = file_handle[:]->PB(time_in,:,:,:)
           QVAPOR = file_handle[:]->QVAPOR(time_in,:,:,:)
         end if
       end if
       P = P + PB
       td  = wrf_td( P , QVAPOR )
       delete_VarAtts(QVAPOR,(/"description","units"/))
       copy_VarAtts(QVAPOR,td)
       return(td)
  end if



  if( variable .eq. "td2" ) then
       ;; function wrf_td needs qv and pressure (Pa) on input and returns dewpoint temperature on return
       if ( time .eq. -1 ) then
         if(ISFILE) then
           PSFC = nc_file->PSFC
           Q2   = nc_file->Q2
         else
           PSFC = file_handle[:]->PSFC
           Q2   = file_handle[:]->Q2
         end if
       else
         if(ISFILE) then
           PSFC = nc_file->PSFC(time_in,:,:) 
           Q2   = nc_file->Q2(time_in,:,:)
         else
           PSFC = file_handle[:]->PSFC(time_in,:,:) 
           Q2   = file_handle[:]->Q2(time_in,:,:)
         end if
       end if
       td = wrf_td( PSFC , Q2 )
       delete_VarAtts(Q2,(/"description","units"/))
       copy_VarAtts(Q2,td)
       td@description = "2m Dewpoint Temperature"   ; Overwrite return description
       return(td)
  end if



  if( variable .eq. "slp" ) then
       if(isfilevar(nc_file,"T")) then
         ;; first compute theta - function wrf_tk needs theta and pressure (Pa) on input
         ;; THEN compute sea level pressure, from qv, p (Pa), tk, z
         if ( time .eq. -1 ) then
           if(ISFILE) then
             T      = nc_file->T
             P      = nc_file->P
             PB     = nc_file->PB
             QVAPOR = nc_file->QVAPOR
             PH     = nc_file->PH
             PHB    = nc_file->PHB
           else
             T      = file_handle[:]->T
             P      = file_handle[:]->P
             PB     = file_handle[:]->PB
             QVAPOR = file_handle[:]->QVAPOR
             PH     = file_handle[:]->PH
             PHB    = file_handle[:]->PHB
           end if
         else
           if(ISFILE) then
             T      = nc_file->T(time_in,:,:,:)
             P      = nc_file->P(time_in,:,:,:)
             PB     = nc_file->PB(time_in,:,:,:)
             QVAPOR = nc_file->QVAPOR(time_in,:,:,:)
             PH     = nc_file->PH(time_in,:,:,:)
             PHB    = nc_file->PHB(time_in,:,:,:)
           else
             T      = file_handle[:]->T(time_in,:,:,:)
             P      = file_handle[:]->P(time_in,:,:,:)
             PB     = file_handle[:]->PB(time_in,:,:,:)
             QVAPOR = file_handle[:]->QVAPOR(time_in,:,:,:)
             PH     = file_handle[:]->PH(time_in,:,:,:)
             PHB    = file_handle[:]->PHB(time_in,:,:,:)
           end if
         end if
         T = T + 300.
         P = P + PB
         QVAPOR = QVAPOR > 0.000
         PH    = ( PH + PHB ) / 9.81
         z = wrf_user_unstagger(PH,PH@stagger)

         tk = wrf_tk( P , T )    ; calculate TK
         slp   = wrf_slp( z, tk, P, QVAPOR )  ; calculate slp
         delete_VarAtts(T,(/"description","units"/))
         copy_VarAtts(T,slp)
       else
         ;; may be a met_em file - see if we can get PMSL
         if(isfilevar(nc_file,"PMSL")) then
           if ( time .eq. -1 ) then
             if(ISFILE) then
               slp = nc_file->PMSL
             else
               slp = file_handle[:]->PMSL
             end if
           else
             if(ISFILE) then
               slp = nc_file->PMSL(time_in,:,:)
             else
               slp = file_handle[:]->PMSL(time_in,:,:)
             end if
           end if
         end if
       end if  

       return(slp)
  end if



  if( variable .eq. "rh" ) then
       if(isfilevar(nc_file,"T")) then
         ;; first compute theta - function wrf_tk needs theta and pressure (Pa) on input
         ;; THEN compute rh, using tk, p (Pa), QVAPOR
         if ( time .eq. -1 ) then
           if(ISFILE) then
             T      = nc_file->T
             P      = nc_file->P
             PB     = nc_file->PB
             QVAPOR = nc_file->QVAPOR
           else
             T      = file_handle[:]->T
             P      = file_handle[:]->P
             PB     = file_handle[:]->PB
             QVAPOR = file_handle[:]->QVAPOR
           end if
         else
           if(ISFILE) then
             T      = nc_file->T(time_in,:,:,:)
             P      = nc_file->P(time_in,:,:,:)
             PB     = nc_file->PB(time_in,:,:,:)
             QVAPOR = nc_file->QVAPOR(time_in,:,:,:)
           else
             T      = file_handle[:]->T(time_in,:,:,:)
             P      = file_handle[:]->P(time_in,:,:,:)
             PB     = file_handle[:]->PB(time_in,:,:,:)
             QVAPOR = file_handle[:]->QVAPOR(time_in,:,:,:)
           end if
         end if  
         T = T + 300.
         P  = P + PB
         Pdims = dimsizes(P)
         QVAPOR = QVAPOR > 0.000
         tk = wrf_tk( P , T )
         rh = wrf_rh( QVAPOR, P, tk )
         delete_VarAtts(T,(/"description","units"/))
         copy_VarAtts(T,rh)
       else
         ;; may be a met_em file - see if we can get RH
         if(isfilevar(nc_file,"RH")) then
           if ( time .eq. -1 ) then
             if(ISFILE) then
               rh = nc_file->RH
             else
               rh = file_handle[:]->RH
             end if
           else
             if(ISFILE) then
               rh = nc_file->RH(time_in,:,:,:)
             else
               rh = file_handle[:]->RH(time_in,:,:,:)
             end if
           end if
         end if
       end if  
       return(rh)
  end if



  if( variable .eq. "rh2" ) then
       if(isfilevar(nc_file,"T2")) then
         ;; Compute rh2, using T2, PSFC, Q2
         if ( time .eq. -1 ) then
           if(ISFILE) then
             T2   = nc_file->T2
             PSFC = nc_file->PSFC
             Q2   = nc_file->Q2
           else
             T2   = file_handle[:]->T2
             PSFC = file_handle[:]->PSFC
             Q2   = file_handle[:]->Q2
           end if
         else
           if(ISFILE) then
             T2   = nc_file->T2(time_in,:,:)
             PSFC = nc_file->PSFC(time_in,:,:)
             Q2   = nc_file->Q2(time_in,:,:)
           else
             T2   = file_handle[:]->T2(time_in,:,:)
             PSFC = file_handle[:]->PSFC(time_in,:,:)
             Q2   = file_handle[:]->Q2(time_in,:,:)
           end if
         end if  
         Q2 = Q2 > 0.000
         rh = wrf_rh( Q2, PSFC, T2 )
         delete_VarAtts(T2,(/"description","units"/))
         copy_VarAtts(T2,rh)
         rh@description = "2m Relative Humidity"           
       else
         ;; may be a met_em file - see if we can get RH
         if(isfilevar(nc_file,"RH")) then
           print("Probably a met_em file - get lowerst level from RH field")
           if ( time .eq. -1 ) then
             if(ISFILE) then
               rh2 = nc_file->RH(:,0,:,:)
             else
               rh2 = file_handle[:]->RH(:,0,:,:)
             end if
           else
             if(ISFILE) then
               rh2 = nc_file->RH(time_in,0,:,:)
             else
               rh2 = file_handle[:]->RH(time_in,0,:,:)
             end if
           end if
         end if
       end if  
       return(rh)
  end if



  if( variable .eq. "pvo" ) then
       if ( time .eq. -1 ) then
         if(ISFILE) then
           U    = nc_file->U
           V    = nc_file->V
           T    = nc_file->T
           P    = nc_file->P
           PB   = nc_file->PB
           MSFU = nc_file->MAPFAC_U
           MSFV = nc_file->MAPFAC_V
           MSFM = nc_file->MAPFAC_M
           COR  = nc_file->F
         else
           U    = file_handle[:]->U
           V    = file_handle[:]->V
           T    = file_handle[:]->T
           P    = file_handle[:]->P
           PB   = file_handle[:]->PB
           MSFU = file_handle[:]->MAPFAC_U
           MSFV = file_handle[:]->MAPFAC_V
           MSFM = file_handle[:]->MAPFAC_M
           COR  = file_handle[:]->F
         end if
       else
         if(ISFILE) then
           U    = nc_file->U(time_in,:,:,:)
           V    = nc_file->V(time_in,:,:,:)
           T    = nc_file->T(time_in,:,:,:)
           P    = nc_file->P(time_in,:,:,:)
           PB   = nc_file->PB(time_in,:,:,:)
           MSFU = nc_file->MAPFAC_U(time_in,:,:)
           MSFV = nc_file->MAPFAC_V(time_in,:,:)
           MSFM = nc_file->MAPFAC_M(time_in,:,:)
           COR  = nc_file->F(time_in,:,:)
         else
           U    = file_handle[:]->U(time_in,:,:,:)
           V    = file_handle[:]->V(time_in,:,:,:)
           T    = file_handle[:]->T(time_in,:,:,:)
           P    = file_handle[:]->P(time_in,:,:,:)
           PB   = file_handle[:]->PB(time_in,:,:,:)
           MSFU = file_handle[:]->MAPFAC_U(time_in,:,:)
           MSFV = file_handle[:]->MAPFAC_V(time_in,:,:)
           MSFM = file_handle[:]->MAPFAC_M(time_in,:,:)
           COR  = file_handle[:]->F(time_in,:,:)
         end if
       end if  
       T = T + 300.
       P = P + PB
       DX = nc_file@DX
       DY = nc_file@DY

       pvo = wrf_pvo( U, V, T, P, MSFU, MSFV, MSFM, COR, DX, DY, 0)

       delete_VarAtts(T,(/"description","units"/))
       copy_VarAtts(T,pvo)
       return(pvo)
  end if



  if( variable .eq. "avo" ) then
       if ( time .eq. -1 ) then
         if(ISFILE) then
           U    = nc_file->U
           V    = nc_file->V
           MSFU = nc_file->MAPFAC_U
           MSFV = nc_file->MAPFAC_V
           MSFM = nc_file->MAPFAC_M
           COR  = nc_file->F
         else
           U    = file_handle[:]->U
           V    = file_handle[:]->V
           MSFU = file_handle[:]->MAPFAC_U
           MSFV = file_handle[:]->MAPFAC_V
           MSFM = file_handle[:]->MAPFAC_M
           COR  = file_handle[:]->F
         end if
       else
         if(ISFILE) then
           U    = nc_file->U(time_in,:,:,:)
           V    = nc_file->V(time_in,:,:,:)
           MSFU = nc_file->MAPFAC_U(time_in,:,:)
           MSFV = nc_file->MAPFAC_V(time_in,:,:)
           MSFM = nc_file->MAPFAC_M(time_in,:,:)
           COR  = nc_file->F(time_in,:,:)
         else
           U    = file_handle[:]->U(time_in,:,:,:)
           V    = file_handle[:]->V(time_in,:,:,:)
           MSFU = file_handle[:]->MAPFAC_U(time_in,:,:)
           MSFV = file_handle[:]->MAPFAC_V(time_in,:,:)
           MSFM = file_handle[:]->MAPFAC_M(time_in,:,:)
           COR  = file_handle[:]->F(time_in,:,:)
         end if
       end if  
       DX = nc_file@DX
       DY = nc_file@DY

       avo = wrf_avo( U, V, MSFU, MSFV, MSFM, COR, DX, DY, 0)

       delete_VarAtts(COR,(/"description","units"/))
       copy_VarAtts(COR,avo)
       return(avo)
  end if



  if( variable .eq. "dbz" .or. variable .eq. "mdbz" ) then
       ; calculate dbz
       ivarint  = 0 
       iliqskin = 0
       dim_vars = dimsizes(varin)
       do idims = 1,dim_vars-1
         if ( idims .eq. 1 ) then
           if ( varin(idims) .eq. "1" ) then
             ivarint = 1
           end if
         end if
         if ( idims .eq. 2 ) then
           if ( varin(idims) .eq. "1" ) then
             iliqskin = 1
           end if
         end if
       end do

       if ( time .eq. -1 ) then
         if(ISFILE) then
           T  = nc_file->T
           P  = nc_file->P
           PB = nc_file->PB
           qv = nc_file->QVAPOR
           qr = nc_file->QRAIN
         else
           T  = file_handle[:]->T
           P  = file_handle[:]->P
           PB = file_handle[:]->PB
           qv = file_handle[:]->QVAPOR
           qr = file_handle[:]->QRAIN
         end if
         if(isfilevar(nc_file,"QSNOW"))
           if(ISFILE) then
             qs = nc_file->QSNOW
           else
             qs = file_handle[:]->QSNOW
           end if
         end if
         if(isfilevar(nc_file,"QGRAUP"))
           if(ISFILE) then
             qg = nc_file->QGRAUP
           else
             qg = file_handle[:]->QGRAUP
           end if
         end if
       else
         if(ISFILE) then
           T  = nc_file->T(time_in,:,:,:)
           P  = nc_file->P(time_in,:,:,:)
           PB = nc_file->PB(time_in,:,:,:)
           qv = nc_file->QVAPOR(time_in,:,:,:)
           qr = nc_file->QRAIN(time_in,:,:,:)
         else
           T  = file_handle[:]->T(time_in,:,:,:)
           P  = file_handle[:]->P(time_in,:,:,:)
           PB = file_handle[:]->PB(time_in,:,:,:)
           qv = file_handle[:]->QVAPOR(time_in,:,:,:)
           qr = file_handle[:]->QRAIN(time_in,:,:,:)
         end if
         if(isfilevar(nc_file,"QSNOW"))
           if(ISFILE) then
             qs = nc_file->QSNOW(time_in,:,:,:)
           else
             qs = file_handle[:]->QSNOW(time_in,:,:,:)
           end if
         end if
         if(isfilevar(nc_file,"QGRAUP"))
           if(ISFILE) then
             qg = nc_file->QGRAUP(time_in,:,:,:)
           else
             qg = file_handle[:]->QGRAUP(time_in,:,:,:)
           end if
         end if
       end if
       T  = T + 300.
       P  = P + PB
       tk = wrf_tk( P , T )

       if ( .not. isvar("qs") ) then
         qs = qv
         qs = 0.0
       end if
       if ( .not. isvar("qg") ) then
         qg = qv
         qg = 0.0
       end if

       dbz = wrf_dbz ( P, tk, qv, qr, qs, qg, ivarint, iliqskin)
       delete(qs)
       delete(qg)

       delete_VarAtts(T,(/"description","units"/))
       copy_VarAtts(T,dbz)


       if ( variable .eq. "mdbz") then
         dims = getvardims(dbz)
         rank = dimsizes(dims)
         if ( rank .eq. 5 ) then
            mdbz = dim_max ( dbz($dims(0)$|:,$dims(1)$|:,$dims(3)$|:,$dims(4)$|:,$dims(2)$|:) )
            mdbz!0 = dbz!0
            mdbz!1 = dbz!1
         end if
         if ( rank .eq. 4 ) then
            mdbz = dim_max ( dbz($dims(0)$|:,$dims(2)$|:,$dims(3)$|:,$dims(1)$|:) )
            mdbz!0 = dbz!0
         end if
         if ( rank .eq. 3 ) then
            mdbz = dim_max ( dbz($dims(1)$|:,$dims(2)$|:,$dims(0)$|:) )
         end if
         nn = rank-1
         nm = rank-2
         mdbz!nm = dbz!nn
         nn = rank-2
         nm = rank-3
         mdbz!nm = dbz!nn
         copy_VarAtts(dbz,mdbz)
         mdbz@description = "Max Reflectivity"
         return(mdbz)
       else
         return(dbz)
       end if

  end if



  if( any( variable .eq. (/"cape_3d","cape_2d"/) ) ) then
       if ( time .eq. -1 ) then
         if(ISFILE) then
           T  = nc_file->T
           P  = nc_file->P
           PB = nc_file->PB
           QV = nc_file->QVAPOR
           PH  = nc_file->PH
           PHB = nc_file->PHB
           HGT = nc_file->HGT
           PSFC = nc_file->PSFC
         else
           T  = file_handle[:]->T
           P  = file_handle[:]->P
           PB = file_handle[:]->PB
           QV = file_handle[:]->QVAPOR
           PH  = file_handle[:]->PH
           PHB = file_handle[:]->PHB
           HGT = file_handle[:]->HGT
           PSFC = file_handle[:]->PSFC
         end if
       else
         if(ISFILE) then
           T  = nc_file->T(time_in,:,:,:)
           P  = nc_file->P(time_in,:,:,:)
           PB = nc_file->PB(time_in,:,:,:)
           QV = nc_file->QVAPOR(time_in,:,:,:)
           PH  = nc_file->PH(time_in,:,:,:)
           PHB = nc_file->PHB(time_in,:,:,:)
           HGT = nc_file->HGT(time_in,:,:)
           PSFC = nc_file->PSFC(time_in,:,:)
         else
           T  = file_handle[:]->T(time_in,:,:,:)
           P  = file_handle[:]->P(time_in,:,:,:)
           PB = file_handle[:]->PB(time_in,:,:,:)
           QV = file_handle[:]->QVAPOR(time_in,:,:,:)
           PH  = file_handle[:]->PH(time_in,:,:,:)
           PHB = file_handle[:]->PHB(time_in,:,:,:)
           HGT = file_handle[:]->HGT(time_in,:,:)
           PSFC = file_handle[:]->PSFC(time_in,:,:)
         end if
       end if
       T = T + 300.
       P  = P + PB 
       tk = wrf_tk( P , T )
       PH =  PH + PHB 
       z = wrf_user_unstagger(PH,PH@stagger)
       z = z/9.81

       if( variable .eq. "cape_3d" ) then
         cape = wrf_cape_3d( P, tk, QV, z, HGT, PSFC, True )
         cape@description = "cape ; cin"
       end if
       if( variable .eq. "cape_2d" ) then
         cape = wrf_cape_2d( P, tk, QV, z, HGT, PSFC, True )
         delete_VarAtts(T,(/"MemoryOrder"/))
         cape@MemoryOrder = "XY"
         cape@description = "mcape ; mcin ; lcl ; lfc"
       end if
       delete_VarAtts(T,(/"description","units"/))
       copy_VarAtts(T,cape)

       return(cape)
  end if



  if( any( variable .eq. (/"pw"/) ) ) then
       ;Precipitable Water 
        print("calculating precipitable water")
        gas_const = 287.  ; J/K/kg
        Cp        = 1004. ; J/K/kg     
       if ( time .eq. -1 ) then
        if(ISFILE) then
          T   = nc_file->T
          P   = nc_file->P
          PB  = nc_file->PB
          PH  = nc_file->PH
          PHB = nc_file->PHB
          QV  = nc_file->QVAPOR
        else
          T   = file_handle[:]->T
          P   = file_handle[:]->P
          PB  = file_handle[:]->PB
          PH  = file_handle[:]->PH
          PHB = file_handle[:]->PHB
          QV  = file_handle[:]->QVAPOR
        end if
       else
        if(ISFILE) then
          T   = nc_file->T(time_in,:,:,:)
          P   = nc_file->P(time_in,:,:,:)
          PB  = nc_file->PB(time_in,:,:,:)
          PH  = nc_file->PH(time_in,:,:,:)
          PHB = nc_file->PHB(time_in,:,:,:)
          QV  = nc_file->QVAPOR(time_in,:,:,:)
        else
          T   = file_handle[:]->T(time_in,:,:,:)
          P   = file_handle[:]->P(time_in,:,:,:)
          PB  = file_handle[:]->PB(time_in,:,:,:)
          PH  = file_handle[:]->PH(time_in,:,:,:)
          PHB = file_handle[:]->PHB(time_in,:,:,:)
          QV  = file_handle[:]->QVAPOR(time_in,:,:,:)
        end if
       end if  
       
       pres   =  P + PB
       height = (PH + PHB)/9.8           ; height at full levels 
       theta  =  T + 300.
       temp   = theta * (pres/100000) ^ (gas_const/Cp)
       vtemp  = (1 + 0.61*QV) * temp  ; virtual temp

       dims = dimsizes(T)
       nd = dimsizes(dims)

       if ( nd .eq. 4 ) then
         zdiff  = height(:,0,:,:)
         zdiff  = 0.
         pw_sfc_ptop = height(:,0,:,:)
         pw_sfc_ptop = 0.
         do k = 0,dims(1)-1
           zdiff(:,:,:) = (height(:,k+1,:,:) - height(:,k,:,:))
           pw_sfc_ptop(:,:,:) =  pw_sfc_ptop(:,:,:) + ((pres(:,k,:,:)/(gas_const * vtemp(:,k,:,:))) * QV(:,k,:,:) * zdiff(:,:,:))
         end do
       end if
       if ( nd .eq. 3 ) then
         zdiff  = height(0,:,:)
         zdiff  = 0.
         pw_sfc_ptop = height(0,:,:)
         pw_sfc_ptop = 0.
         do k = 0,dims(0)-1
           zdiff(:,:) = (height(k+1,:,:) - height(k,:,:))
           pw_sfc_ptop(:,:) =  pw_sfc_ptop(:,:) + ((pres(k,:,:)/(gas_const * vtemp(k,:,:))) * QV(k,:,:) * zdiff(:,:))
         end do
       end if

       pw_sfc_ptop@description = "Precipitable Water"
       return(pw_sfc_ptop)
  end if



  if( any( variable .eq. (/"helicity"/) ) ) then
       getU = "U" 
       getV = "V" 
       if(.not. isfilevar(nc_file,"U")) then
         if(isfilevar(nc_file,"UU")) then
           getU = "UU"
           getV = "VV"
        end if
       end if
       if ( time .eq. -1 ) then
        if(ISFILE) then
          u_in  = nc_file->$getU$
          v_in  = nc_file->$getV$
          PH    = nc_file->PH
          geopt = nc_file->PHB 
          ter   = nc_file->HGT
        else
          u_in  = file_handle[:]->$getU$
          v_in  = file_handle[:]->$getV$
          PH    = file_handle[:]->PH
          geopt = file_handle[:]->PHB 
          ter   = file_handle[:]->HGT
        end if
       else
        if(ISFILE) then
          u_in  = nc_file->$getU$(time_in,:,:,:)
          v_in  = nc_file->$getV$(time_in,:,:,:)
          PH    = nc_file->PH(time_in,:,:,:)
          geopt = nc_file->PHB(time_in,:,:,:)
          ter   = nc_file->HGT(time_in,:,:)
        else
          u_in  = file_handle[:]->$getU$(time_in,:,:,:)
          v_in  = file_handle[:]->$getV$(time_in,:,:,:)
          PH    = file_handle[:]->PH(time_in,:,:,:)
          geopt = file_handle[:]->PHB(time_in,:,:,:)
          ter   = file_handle[:]->HGT(time_in,:,:)
        end if
       end if

       ua    = wrf_user_unstagger(u_in,u_in@stagger)
       va    = wrf_user_unstagger(v_in,v_in@stagger)
       geopt  = geopt + PH
       za     = wrf_user_unstagger(geopt,geopt@stagger)
       za     = za / 9.81  ; change to height

       ua1    = ua(::-1,:,:)
       va1    = va(::-1,:,:)
       za1    = za(::-1,:,:)

       top_ok = 0
       top = 3000.

       dim_vars = dimsizes(varin)
       if(dim_vars .eq. 2) then
        if( varin(1) .eq. "3000" ) then
           top = 3000.
           top_ok = 1
        end if

        if( varin(1) .eq. "1000" ) then
           top = 1000.
           top_ok = 1
        end if

        if(top_ok .eq. 0) then
          print("Top values of 1000 or 3000 are accepted.")
          top = 3000.
        end if
       end if

       print("Calculations are done with a top of " + top)
       sreh = wrf_helicity(ua1, va1, za1, ter, top)
       return(sreh)
  end if



  if( any(variable .eq. (/"updraft_helicity"/) ) )then
       getU = "U"
       getV = "V"
       if(.not. isfilevar(nc_file,"U")) then
         if(isfilevar(nc_file,"UU")) then
           getU = "UU"
           getV = "VV"
         end if
       end if
       if ( time .eq. -1 ) then
         if(ISFILE) then
            u_in = nc_file->$getU$
            v_in = nc_file->$getV$
            w    = nc_file->W
            ph   = nc_file->PH
            phb  = nc_file->PHB
         else
            u_in = file_handle[:]->$getU$
            v_in = file_handle[:]->$getV$
            w    = file_handle[:]->W
            ph   = file_handle[:]->PH
            phb  = file_handle[:]->PHB
         end if
       else
         if(ISFILE) then
            u_in = nc_file->$getU$(time_in,:,:,:)
            v_in = nc_file->$getV$(time_in,:,:,:)
            w    = nc_file->W(time_in,:,:,:)
            ph   = nc_file->PH(time_in,:,:,:)
            phb  = nc_file->PHB(time_in,:,:,:)
         else
            u_in = file_handle[:]->$getU$(time_in,:,:,:)
            v_in = file_handle[:]->$getV$(time_in,:,:,:)
            w    = file_handle[:]->W(time_in,:,:,:)
            ph   = file_handle[:]->PH(time_in,:,:,:)
            phb  = file_handle[:]->PHB(time_in,:,:,:)
         end if
       end if
       ua   = wrf_user_unstagger(u_in,u_in@stagger)
       va   = wrf_user_unstagger(v_in,v_in@stagger)
       mapfct = nc_file->MAPFAC_M(0,:,:)
       zp = ph + phb
       dx = nc_file@DX
       dy = nc_file@DY

       uh_opt = True
       uh_opt@uhmnhgt = 2000.
       uh_opt@uhmxhgt = 5000.
       dim_vars = dimsizes(varin)

       print("Calculating updraft helicity")
       if(dim_vars .eq. 1) then
         print("  Using defaults for the integration limits")
       end if
       
       if(dim_vars .eq. 2) then
         print("  Please enter both the minimum and maximum integration limits")
         print("  Going to use defaults for the integration limits")
       end if
       
       if(dim_vars .eq. 3) then
         if ( stringtofloat(varin(1)) .lt. 1000. ) then
           print("  Integration limits needs to be greater than 1000 meter")
           print("  Going to use defaults for the integration limits")
         else
           uh_opt@uhmnhgt = stringtofloat(varin(1))
           uh_opt@uhmxhgt = stringtofloat(varin(2))
           print("  Setting custom integration limits")
         end if
       end if
       
       print("  min = "+uh_opt@uhmnhgt+" max = "+uh_opt@uhmxhgt)
       uh = wrf_updraft_helicity(zp, mapfct, ua, va, w, dx, dy, uh_opt)
       delete(uh_opt)
       
       return(uh)
  end if


  if( variable .eq. "twb" ) then
       ; Wet bulb temperature
       if ( time .eq. -1 ) then
         if(ISFILE) then
           T  = nc_file->T
           P  = nc_file->P
           PB = nc_file->PB
           QV = nc_file->QVAPOR
         else
           T  = file_handle[:]->T
           P  = file_handle[:]->P
           PB = file_handle[:]->PB
           QV = file_handle[:]->QVAPOR
         end if
       else
         if(ISFILE) then
           T  = nc_file->T(time_in,:,:,:)
           P  = nc_file->P(time_in,:,:,:)
           PB = nc_file->PB(time_in,:,:,:)
           QV = nc_file->QVAPOR(time_in,:,:,:)
         else
           T  = file_handle[:]->T(time_in,:,:,:)
           P  = file_handle[:]->P(time_in,:,:,:)
           PB = file_handle[:]->PB(time_in,:,:,:)
           QV = file_handle[:]->QVAPOR(time_in,:,:,:)
         end if
       end if
       T = T + 300.
       P = P + PB
       t = wrf_tk(P,T)

       twb = wrf_wetbulb(P,t,QV)

       delete_VarAtts(T,(/"description"/))
       copy_VarAtts(T,twb)
       return(twb)
  end if


  if( variable .eq. "omg" ) then
       ; Omega
       if ( time .eq. -1 ) then
         if(ISFILE) then
           T  = nc_file->T
           P  = nc_file->P
           W  = nc_file->W
           PB = nc_file->PB
           QV = nc_file->QVAPOR
         else
           T  = file_handle[:]->T
           P  = file_handle[:]->P
           W  = file_handle[:]->W
           PB = file_handle[:]->PB
           QV = file_handle[:]->QVAPOR
         end if
       else
         if(ISFILE) then
           T  = nc_file->T(time_in,:,:,:)
           P  = nc_file->P(time_in,:,:,:)
           W  = nc_file->W(time_in,:,:,:)
           PB = nc_file->PB(time_in,:,:,:)
           QV = nc_file->QVAPOR(time_in,:,:,:)
         else
           T  = file_handle[:]->T(time_in,:,:,:)
           P  = file_handle[:]->P(time_in,:,:,:)
           W  = file_handle[:]->W(time_in,:,:,:)
           PB = file_handle[:]->PB(time_in,:,:,:)
           QV = file_handle[:]->QVAPOR(time_in,:,:,:)
         end if
       end if
       T = T + 300.
       P = P + PB
       t = wrf_tk(P,T)
       wa = wrf_user_unstagger(W,W@stagger)

       omg = wrf_omega(QV,t,wa,P)

       delete_VarAtts(T,(/"description","units"/))
       copy_VarAtts(T,omg)
       return(omg)
  end if


  if( variable .eq. "tv" ) then
       ; Virtual temperature
       if ( time .eq. -1 ) then
         if(ISFILE) then
           T  = nc_file->T
           P  = nc_file->P
           PB = nc_file->PB
           QV = nc_file->QVAPOR
         else
           T  = file_handle[:]->T
           P  = file_handle[:]->P
           PB = file_handle[:]->PB
           QV = file_handle[:]->QVAPOR
         end if
       else
         if(ISFILE) then
           T  = nc_file->T(time_in,:,:,:)
           P  = nc_file->P(time_in,:,:,:)
           PB = nc_file->PB(time_in,:,:,:)
           QV = nc_file->QVAPOR(time_in,:,:,:)
         else
           T  = file_handle[:]->T(time_in,:,:,:)
           P  = file_handle[:]->P(time_in,:,:,:)
           PB = file_handle[:]->PB(time_in,:,:,:)
           QV = file_handle[:]->QVAPOR(time_in,:,:,:)
         end if
       end if
       T = T + 300.
       P = P + PB
       t = wrf_tk(P,T)

       tv = wrf_virtual_temp(t,QV)

       delete_VarAtts(T,(/"description"/))
       copy_VarAtts(T,tv)
       return(tv)
  end if

  if( variable .eq. "ctt") then
      if ( time .eq. -1 ) then
         explicit_time = 0
         if(ISFILE) then
            P     = nc_file->P
            PB    = nc_file->PB
            T     = nc_file->T

;geopotential height
            PHB   = nc_file->PHB
            PH    = nc_file->PH

            if(isfilevar(nc_file,"QVAPOR")) then
               qvp =  nc_file->QVAPOR * 1000.  ; kg/kg -> g/kg
            else
               print("wrf_user_getvar: QVAPOR is needed to calculate the cloud top temperature.")
               print("It has not been found in the data set.")
               exit
            end if


            haveqci = 1
            if(isfilevar(nc_file,"QICE")) then
               qci = nc_file->QICE * 1000.  ; kg/kg -> g/kg
            else
               qci = new( (/Pdims(0),Pdims(1),Pdims(2),Pdims(3)/),float)
               haveqci = 0
            end if


            if(isfilevar(nc_file,"QCLOUD")) then
               qcw = nc_file->QCLOUD * 1000.  ; kg/kg -> g/kg
            else
               print("wrf_user_getvar: QCLOUD is needed to calculate the cloud top temperature.")
               print("It has not been found in the data set.")
               exit
            end if

            ter = nc_file->HGT(0,:,:)

         else ; the else for ISFILE
            P      = file_handle[:]->P
            PB     = file_handle[:]->PB
            T      = file_handle[:]->T
            PHB    = file_handle[:]->PHB
            PH     = file_handle[:]->PH
            nc_file =  file_handle[0]
            ter    = nc_file->HGT(0,:,:)
            if(isfilevar(nc_file,"QVAPOR")) then
               qvp = file_handle[:]->QVAPOR * 1000.  ; kg/kg -> g/kg
            else
               print("wrf_user_getvar: QVAPOR is needed to calculate the cloud top temperature.")
               print("It has not been found in the data set.")
               exit
            end if 
 
            haveqci = 1
            if(isfilevar(nc_file,"QICE")) then
               qci = file_handle[:]->QICE * 1000.  ; kg/kg -> g/kg
            else
               qci = new( (/Pdims(0),Pdims(1),Pdims(2),Pdims(3)/),float)
               haveqci = 0
            end if


            if(isfilevar(nc_file,"QCLOUD")) then
               qcw = file_handle[:]->QCLOUD * 1000.  ; kg/kg -> g/kg
            else
               print("wrf_user_getvar: QCLOUD is needed to calculate the cloud top temperature.")
               print("It has not been found in the data set.")
               exit
            end if

         end if ;if ISFILE
       else     ;the else for time = -1
          if(ISFILE) then
            explicit_time = 1
            P     = nc_file->P(time_in,:,:,:)
            PB    = nc_file->PB(time_in,:,:,:)
            T     = nc_file->T(time_in,:,:,:)

;geopotential height
            PHB   = nc_file->PHB(time_in,:,:,:)
            PH    = nc_file->PH(time_in,:,:,:)

            if(isfilevar(nc_file,"QVAPOR")) then
               qvp =  nc_file->QVAPOR(time_in,:,:,:) * 1000.  ; kg/kg -> g/kg
            else
               print("wrf_user_getvar: QVAPOR is needed to calculate the cloud top temperature.")
               print("It has not been found in the data set")
               exit
            end if


            haveqci = 1
            if(isfilevar(nc_file,"QICE")) then
               qci = nc_file->QICE(time_in,:,:,:) * 1000.  ; kg/kg -> g/kg
            else
               qci = new( (/Pdims(0),Pdims(1),Pdims(2),Pdims(3)/),float)
               haveqci = 0
            end if


            if(isfilevar(nc_file,"QCLOUD")) then
               qcw = nc_file->QCLOUD(time_in,:,:,:) * 1000.  ; kg/kg -> g/kg
            else
               print("wrf_user_getvar: QCLOUD is needed to calculate the cloud top temperature.")
               print("It has not been found in the data set.")
               exit
            end if

            ter = nc_file->HGT(0,:,:)
          end if ;end if for ISFILE
       end if   ;time = -1


;Get total pressure
         Pdims = dimsizes(P)
         pres  = P + PB
         pres  = pres * 0.01

;Get temperature in degrees K
         T     = T +300.
         tk    = wrf_tk( pres , T )

;Get geopotential height on mass points
         stag_ght   = PH
         stag_ght   = (PHB + PH)/9.81
         ght        =  wrf_user_unstagger(stag_ght,PHB@stagger)

         fctt = wrf_ctt(pres,tk,qci,qcw,qvp,ght,ter,haveqci)
; fctt will have description, units, and dimension names attached
       delete_VarAtts(T,(/"description","units"/))
       copy_VarAtts(T,fctt)
       return(fctt)
  end if         ;variable is ctt
 

  if( any( variable .eq. (/"ter","HGT","HGT_M"/) ) ) then
    variable = "HGT" 
    if(.not. isfilevar(nc_file,"HGT")) then
      variable = "HGT_M"
    end if
  end if


  if( any( variable .eq. (/"lat","XLAT","XLAT_M"/) ) ) then
    variable = "XLAT" 
    if(.not. isfilevar(nc_file,"XLAT")) then
      variable = "XLAT_M"
    end if
  end if


  if( any( variable .eq. (/"lon","long","XLONG","XLONG_M"/) ) ) then
    variable = "XLONG" 
    if(.not. isfilevar(nc_file,"XLONG")) then
      variable = "XLONG_M"
    end if
  end if


  if( any( variable .eq. (/"times"/) ) ) then
    if ( time .eq. -1 ) then
      if(ISFILE) then
        var = nc_file->Times
      else
        var = file_handle[:]->Times
      end if
    else
      if(ISFILE) then
        var = nc_file->Times(time_in,:)
      else
        var = file_handle[:]->Times(time_in,:)
      end if
    end if
    dims = dimsizes(var)
    times = new(dims(0),string)
    do i=0,dims(0)-1
      times(i) = chartostring(var(i,:))
    end do
    times@description = "times in file"
    return (times)
  end if



;  end of diagnostic variable list - we must want a variable already in the file.  

  if ( time .eq. -1 ) then
    if(ISFILE) then
      var = nc_file->$variable$
    else
      var = file_handle[:]->$variable$
    end if
  else
    ;  check variable dimensionality and pull proper time  out of file
    ndims = dimsizes(filevardimsizes(nc_file,variable))
    if( ndims .eq. 4) then
      if(ISFILE) then
        var = nc_file->$variable$(time_in,:,:,:)
      else
        var = file_handle[:]->$variable$(time_in,:,:,:)
      end if
    end if
    if( ndims .eq. 3) then
      if(ISFILE) then
        var = nc_file->$variable$(time_in,:,:)
      else
        var = file_handle[:]->$variable$(time_in,:,:)
      end if
    end if
    if( ndims .eq. 2) then
      if(ISFILE) then
        var = nc_file->$variable$(time_in,:)
      else
        var = file_handle[:]->$variable$(time_in,:)
      end if
    end if
    if( ndims .eq. 1) then
      if(ISFILE) then
      var = nc_file->$variable$(time_in)
      else
      var = file_handle[:]->$variable$(time_in)
      end if
    end if
  end if

  return(var)

end


;--------------------------------------------------------------------------------
undef("wrf_user_getvar_from_files")
function wrf_user_getvar_from_files( nc_files:string, varin[*]:string, time[*]:integer, opts_args:logical )

begin

  numFiles = dimsizes(nc_files)
  f = addfiles(nc_files,"r")

  var_tmp = wrf_user_getvar(f[0],varin,-1)
  varDims = dimsizes(var_tmp)
  varRank = dimsizes(varDims)
  time_req = 0
  time_inc = 1
  time_end = varDims(0)*numFiles
  if ( dimsizes(time) .eq. 1 ) then
    if ( time(0) .ge. 0 ) then
      numTimes = 1
      time_req = time(0)
      time_end = time(0)
      time_inc = 1
    else
      numTimes = varDims(0)*numFiles
    end if
  end if
  if ( dimsizes(time) .eq. 2 ) then
    if ( time(0) .ge. 0 ) then
      numTimes = (time(1)-time(0))+1
      time_req = time(0)
      time_end = time(1)
      time_inc = 1
    else
      numTimes = varDims(0)*numFiles
    end if
  end if
  if ( dimsizes(time) .eq. 3 ) then
    if ( time(0) .ge. 0 ) then
      numTimes = ((time(1)-time(0))+1)/time(2)
      time_req = time(0)
      time_end = time(1)
      time_inc = time(2)
    else
      numTimes = varDims(0)*numFiles
    end if
  end if
  outtime = 0
  outtime_avail = 0


  if ( varRank .eq. 4 )
    varout = new ( (/numTimes,varDims(1),varDims(2),varDims(3)/), typeof(var_tmp) )
    varout = 0
    do it = 0,varDims(0)-1
      if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
        varout(outtime,:,:,:) = var_tmp(it,:,:,:)
        outtime = outtime + 1
        time_req = time_req + time_inc
      end if
      outtime_avail = outtime_avail + 1
    end do
    delete(var_tmp)
    do ifil = 1,numFiles-1
      var_tmp = wrf_user_getvar(f[ifil],varin,-1)
      dimLoop = dimsizes(var_tmp)
      do it = 0,dimLoop(0)-1
        if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
          varout(outtime,:,:,:) = var_tmp(it,:,:,:)
          outtime = outtime + 1
          time_req = time_req + time_inc
        end if
        outtime_avail = outtime_avail + 1
      end do
      delete(var_tmp)
    end do
  end if
  if ( varRank .eq. 3 )
    varout = new ( (/numTimes,varDims(1),varDims(2)/), typeof(var_tmp) )
    varout = 0
    do it = 0,varDims(0)-1
      if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
        varout(outtime,:,:) = var_tmp(it,:,:)
        outtime = outtime + 1
        time_req = time_req + time_inc
      end if
      outtime_avail = outtime_avail + 1
    end do
    delete(var_tmp)
    do ifil = 1,numFiles-1
      var_tmp = wrf_user_getvar(f[ifil],varin,-1)
      dimLoop = dimsizes(var_tmp)
      do it = 0,dimLoop(0)-1
        if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
          varout(outtime,:,:) = var_tmp(it,:,:)
          outtime = outtime + 1
          time_req = time_req + time_inc
        end if
        outtime_avail = outtime_avail + 1
      end do
      delete(var_tmp)
    end do
  end if
  if ( varRank .eq. 2 )
    varout = new ( (/numTimes,varDims(1)/), typeof(var_tmp) )
    if ( typeof(var_tmp) .eq. "float" .or. typeof(var_tmp) .eq. "integer" ) then
      varout = 0
    end if
    do it = 0,varDims(0)-1
      if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
        varout(outtime,:) = var_tmp(it,:)
        outtime = outtime + 1
        time_req = time_req + time_inc
      end if
      outtime_avail = outtime_avail + 1
    end do
    delete(var_tmp)
    do ifil = 1,numFiles-1
      var_tmp = wrf_user_getvar(f[ifil],varin,-1)
      dimLoop = dimsizes(var_tmp)
      do it = 0,dimLoop(0)-1
        if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
          varout(outtime,:) = var_tmp(it,:)
          outtime = outtime + 1
          time_req = time_req + time_inc
        end if
        outtime_avail = outtime_avail + 1
      end do
      delete(var_tmp)
    end do
  end if
  if ( varRank .eq. 1 )
    varout = new ( (/numTimes/), typeof(var_tmp) )
    if ( typeof(var_tmp) .eq. "float" .or. typeof(var_tmp) .eq. "integer" ) then
      varout = 0
    end if
    do it = 0,varDims(0)-1
      if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
        varout(outtime) = var_tmp(it)
        outtime = outtime + 1
        time_req = time_req + time_inc
      end if
      outtime_avail = outtime_avail + 1
    end do
    delete(var_tmp)
    do ifil = 1,numFiles-1
      var_tmp = wrf_user_getvar(f[ifil],varin,-1)
      dimLoop = dimsizes(var_tmp)
      do it = 0,dimLoop(0)-1
        if ( (outtime_avail.eq.time_req) .and. (time_req.le.time_end) ) then
          varout(outtime) = var_tmp(it)
          outtime = outtime + 1
          time_req = time_req + time_inc
        end if
        outtime_avail = outtime_avail + 1
      end do
      delete(var_tmp)
    end do
  end if

  return(varout)
end




;--------------------------------------------------------------------------------
; This function was modified in May 2011 to allow a list of files.
;
undef("wrf_user_list_times")
function wrf_user_list_times( nc_file )

local times, times_in_file, dims, i
begin

;---As of NCL V6.0.0, wrf_user_getvar can now handle a file or a list of files.
  if(all(typeof(nc_file).ne.(/"file","list"/))) then
    print("wrf_user_list_times: error: the input argument must be a file or a list of files opened with addfile or addfiles")
    return
  end if

  if(typeof(nc_file).eq."file") then
    times_in_file = nc_file->Times
  else
    times_in_file = nc_file[:]->Times
  end if
  dims = dimsizes(times_in_file)
  times = new(dims(0),string)
  do i=0,dims(0)-1
    times(i) = chartostring(times_in_file(i,:))
  end do
  times@description = "times in file"
  print(times)
  return(times)

end

;--------------------------------------------------------------------------------

undef("wrf_user_latlon_to_ij")
function wrf_user_latlon_to_ij( nc_file:file, latitude:numeric, \
                                longitude:numeric )

begin
  WE       = "WEST-EAST_GRID_DIMENSION"
  SN       = "SOUTH-NORTH_GRID_DIMENSION"
  wedim    = nc_file@$WE$
  sndim    = nc_file@$SN$

  if(isfilevar(nc_file,"XLAT"))
    XLAT  = nc_file->XLAT(0,:,:)
    XLONG = nc_file->XLONG(0,:,:)
  else
    XLAT  = nc_file->XLAT_M(0,:,:)
    XLONG = nc_file->XLONG_M(0,:,:)
  end if

  loc = wrf_latlon_to_ij( XLAT, XLONG, latitude, longitude )

  loc!0 = "j & i locations"
  return(loc)

end

;--------------------------------------------------------------------------------

undef("wrf_user_ll_to_ij")
function wrf_user_ll_to_ij( file_handle, longitude:numeric, latitude:numeric, \
                            opts_args:logical )

begin
;
; As of NCL V6.0.0, wrf_user_ll_to_ij can now handle a file 
; or a list of files.
;
  if(typeof(file_handle).eq."file") then
    ISFILE = True
    nc_file = file_handle
  else if(typeof(file_handle).eq."list") then
    ISFILE = False
    nc_file = file_handle[0]
  else
    print("wrf_user_ll_to_ij: error: the first argument must be a file or a list of files opened with addfile or addfiles")
    return
  end if
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
      XLAT  = file_handle[useT]->XLAT
      XLONG = file_handle[useT]->XLONG
    end if
  else
    if(ISFILE) then
      XLAT  = nc_file->XLAT_M(useT,:,:)
      XLONG = nc_file->XLONG_M(useT,:,:)
    else
      XLAT  = file_handle[useT]->XLAT_M
      XLONG = file_handle[useT]->XLONG_M
    end if
  end if


  if(dimsizes(dimsizes(XLAT)).eq.2) then
; Rank 2
    res@REF_LAT = XLAT(0,0)
    res@REF_LON = XLONG(0,0)
  else
; Rank 3
    res@REF_LAT = XLAT(0,0,0)
    res@REF_LON = XLONG(0,0,0)
  end if
  res@KNOWNI  = 1.0
  res@KNOWNJ  = 1.0

  loc = wrf_ll_to_ij (longitude, latitude, res)

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

undef("wrf_user_ij_to_ll")
function wrf_user_ij_to_ll( file_handle, i:numeric, j:numeric, \
                            opts_args:logical )

begin
;
; As of NCL V6.0.0, wrf_user_ll_to_ij can now handle a file 
; or a list of files.
;
  if(typeof(file_handle).eq."file") then
    ISFILE = True
    nc_file = file_handle
  else if(typeof(file_handle).eq."list") then
    ISFILE = False
    nc_file = file_handle[0]
  else
    print("wrf_user_ij_to_ll: error: the first argument must be a file or a list of files opened with addfile or addfiles")
    return
  end if
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
      XLAT  = file_handle[useT]->XLAT
      XLONG = file_handle[useT]->XLONG
    end if
  else
    if(ISFILE) then
      XLAT  = nc_file->XLAT_M(useT,:,:)
      XLONG = nc_file->XLONG_M(useT,:,:)
    else
      XLAT  = file_handle[useT]->XLAT_M
      XLONG = file_handle[useT]->XLONG_M
    end if
  end if

  if(dimsizes(dimsizes(XLAT)).eq.2) then
; Rank 2
    res@REF_LAT = XLAT(0,0)
    res@REF_LON = XLONG(0,0)
  else
; Rank 3
    res@REF_LAT = XLAT(0,0,0)
    res@REF_LON = XLONG(0,0,0)
  end if
  res@KNOWNI  = 1.0
  res@KNOWNJ  = 1.0

  loc = wrf_ij_to_ll (i,j,res)

  return(loc)


end

;--------------------------------------------------------------------------------
undef("wrf_user_vert_interp")
function wrf_user_vert_interp(file_handle,field:numeric, \
                              vert_coordinate[1]:string, \
                              interp_levels[*]:numeric,opts[1]:logical)
local valid_vert_coords, nc_file, valid_field_types
begin

  valid_vert_coords = (/"pressure","pres","ght_msl","ght_agl","theta","theta-e"/)
  if(.not.any(vert_coordinate.eq.valid_vert_coords)) then
    print("wrf_user_vert_interp: Unrecognized vertical coordinate.")
    print("   Accepted vertical coordinates are:")
    print(     "pressure, pres   hPa")
    print(     "ght_msl          km")
    print(     "ght_agl          km")
    print(     "theta            K")
    print(     "theta-e          K")
    exit
  end if

  if(typeof(file_handle).eq."file") then
    ISFILE = True
    nc_file = file_handle
  else if(typeof(file_handle).eq."list") then
    ISFILE = False
    nc_file = file_handle[0]
  else
    print("wrf_user_vert_interp: error: the first argument must be a file or a list of files opened with addfile or addfiles")
    return
  end if
  end if

  rgas    = 287.04     ;J/K/kg
  ussalr  = .0065      ; deg C per m
  sclht   = rgas*256./9.81
  
;read grid sizes from the nc_file.  Check to make
;sure that we don't have a staggered field.
  dNames      = getvardims(nc_file)
  dSizes      = getfiledimsizes(nc_file)
  thedims     = dimsizes(dSizes)
  number_dims = thedims(0)

  ew = dSizes(ind(dNames .eq. "west_east"))
  ns = dSizes(ind(dNames .eq. "south_north"))
  nz = dSizes(ind(dNames .eq. "bottom_top"))

  field_dims     = dimsizes(field)
  num_field_dims = dimsizes(field_dims)

 if(num_field_dims .lt. 3 .or. num_field_dims .gt. 4) then
    print("wrf_user_vert_interp: We can only interpolate a 3D or 4D field")
    exit
 end if

; Check to see if the field comes in unstaggered.
  if(field_dims(num_field_dims-3) .ne. nz) then
     print("wrf_user_vert_interp: Please unstagger the field in the vertical")
     exit
  end if

  if(field_dims(num_field_dims-2) .ne. ns .or. \
     field_dims(num_field_dims-1) .ne. ew) then
     print("wrf_user_vert_interp: Please unstagger the field")
     exit
  end if

; See what options we have
  extrap_field    = get_res_value_keep(opts,"extrapolate",False)
  field_type      = str_lower(get_res_value_keep(opts,"field_type","none"))
  log_of_Pressure = get_res_value_keep(opts,"logP",False)
  debug           = get_res_value_keep(opts,"debug",False)
  
  valid_field_types = (/"none","pressure","pres","p","z","t","ght"/)
  if(.not.any(field_type.eq.valid_field_types)) then
    print("wrf_user_vert_interp: Unrecognized field type.")
    print("Valid field types are: " + str_join(valid_field_types,", "))
    exit
  end if

  if(log_of_Pressure) then
     logP = 1      ; this is for passing into Fortran
  else
     logP = 0
  end if

  icase  = 0
  extrap = 0   ; This is for passing into Fortran
  if(extrap_field) then
     extrap = 1
     if(any(field_type .eq. (/"p","pres","pressure"/))) then
        icase = 1
     end if

     if(field_type .eq. "z") then
        icase = 2
     end if

     if(field_type .eq. "ght") then
        icase = 2
     end if
;
; For temperature we may have just temperature, potential temperature
; or equivalent potential temperature.  Check the field description attribute 
; to see which one we have.
;
     if(field_type .eq. "t") then
        if(isatt(field,"description")) then

           if(field@description .eq. "Temperature") then
              if(field@units .eq. "C") then
                 icase = 3
              else
                 icase = 4
              end if
           end if

           if(field@description .eq. "Potential Temperature (theta) ") then
              icase = 5
           end if

           if(field@description .eq. "Equivalent Potential Temperature") then
              icase = 6
           end if
        end if ;the endif for checking for the field description attribute

     end if;end if for the field_type being T or t


  end if; the endif for extrap_field .eq. True


  numlevels = dimsizes(interp_levels)

;We will need some basic fields for the interpolation
;regardless of the field requested. Get all time periods
;of the fields.
  if(ISFILE) then
      P     = nc_file->P + nc_file->PB
      Pdims = dimsizes(P)
      ght   = wrf_user_getvar(nc_file,"height",-1)
      tk    = wrf_user_getvar(nc_file,"tk",-1)
      qvp   = nc_file->QVAPOR
      terht = nc_file->HGT
      sfp   = nc_file->PSFC * 0.01
    else
      P     = file_handle[:]->P + file_handle[:]->PB
      Pdims = dimsizes(P)
      tmpz  = file_handle[:]->PH
      PHB   = file_handle[:]->PHB
      tmpz  = (tmpz + PHB)/9.81
      ght   = wrf_user_unstagger(tmpz,"Z")
      T     = file_handle[:]->T
      T     = T + 300.
      tk    = wrf_tk( P , T )
      qvp   = file_handle[:]->QVAPOR
      terht = file_handle[:] ->HGT
      sfp   = file_handle[:] ->PSFC * 0.01
    end if
    smsfp = sfp
    wrf_smooth_2d(smsfp,3)  

      
;Initialize an array for the vertical coordinate
 ntimes  = Pdims(0)

;Get the vertical coordinate type
  vcor  = 0
  logp  = 0
  if(any(vert_coordinate .eq. (/"pressure","pres"/))) then
     vcor = 1
     vcord_array = P * 0.01
  end if

  if(vert_coordinate .eq. "ght_msl") then
    vcor = 2
    vcord_array = exp(-ght/sclht)
  end if

  if(vert_coordinate .eq. "ght_agl") then
    vcor = 3
    rtemp       = new( (/nz,ns,ew/),float)
    vcord_array = new((/ntimes,nz,ns,ew/),float)
    do it = 0, ntimes - 1
       do ilev = 0,nz-1
          rtemp(ilev,:,:) = ght(it,ilev,:,:) - terht(0,:,:)
       end do
       vcord_array(it,:,:,:) = exp(-rtemp/sclht)
    end do
    delete(rtemp)
  end if

  if(vert_coordinate .eq. "theta")   then
     vcor = 4
     idir   = 1       
     icorsw = 0       
     delta  = 0.01    
     if(ISFILE) then        
        coriolis   = nc_file->F(0,:,:)   
        theta      = wrf_user_getvar(nc_file,"theta",-1)
     else
        coriolis   =  file_handle[0]->F(0,:,:)   
        theta      = T
     end if        
     preshPa    = P * 0.01
     vcord_array = wrf_monotonic(theta,preshPa,coriolis,idir,delta,icorsw)      
;
;  We only extrapolate temperature fields below ground if we are interpolating
;  to pressure or height vertical surfaces.
;
     icase = 0       
  end if

  if(vert_coordinate .eq. "theta-e") then
     vcor   = 5
     icorsw = 0
     idir   = 1 
     delta  = 0.01 
     if(ISFILE) then  
        coriolis = nc_file->F(0,:,:) 
        eqpot    = wrf_user_getvar(nc_file,"eth",-1)
     else 
        coriolis =  file_handle[0]->F(0,:,:)
        eqpot    = wrf_eth ( qvp, tk, P )
     end if
     preshPa    = P * 0.01

     vcord_array = wrf_monotonic(eqpot,preshPa,coriolis,idir,delta,icorsw)
;  We only extrapolate temperature fields below ground if we are interpolating
;  to pressure or height vertical surfaces.
     icase = 0
  end if

  if(debug) then
    print("icase = " + icase + " extrap = " + extrap + " vcor = " + vcor + " logP = " + logP)
  end if

  field_out = wrf_vintrp(field,P,tk,qvp,ght,terht(0,:,:),sfp,smsfp,\
                         vcord_array,interp_levels,icase,extrap,vcor,logP)

; Add metadata to return array
  copy_VarMeta(field,field_out)

; Add new levels as a coordinate array
  lev_field           = num_field_dims-3
  field_out!lev_field = "interp_levels"
  field_out&$field_out!lev_field$ = interp_levels(::-1)   
  field_out@vert_interp_type = vert_coordinate

  return(field_out)
end

;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
undef("write_wrf_debug_info")
procedure write_wrf_debug_info(wks[1]:graphic,data1,data2, \
                               debug_file[1]:string, \
                               wrf_opts,wrf_func_name[1]:string)
;
; This procedure writes resources and data variables used to
; create a WRF plot to a NetCDF file. This file can be read
; in later to recreate the NCL script that created that plot,
; using gsn_xxx functions. This is for debug purposes. You can't
; tell what resources a WRF script is using, so this is a way
; of seeing them all in a NetCDF file.
;
begin
  if(.not.isatt(wks,"WRFDebug")) then
    first_time = True
    wks@WRFDebug = True
  else
    first_time = False
  end if
;
; The default will be "wrfdebug.ncl" and "wrfdebug.nc" 
; unless otherwise specified.
;
  cdf_debug_file = debug_file +  ".nc"
  ncl_debug_file = debug_file +  ".ncl"
  res_debug_file = debug_file +  ".res"

;
; If this is the first time writing debug information to the file,
; then create the file and add the first set of information.
;
; Make sure the files don't already exist.
;
  if(first_time) then
    if(fileexists(cdf_debug_file).or.fileexists(ncl_debug_file).or.\
       fileexists(res_debug_file)) then
      print("write_wrf_debug_info: error: debug files '" + cdf_debug_file + "',")
      print("'" + ncl_debug_file + "' and/or " + res_debug_file + " exist.")
      print("Please remove file(s) and start script again.")
      exit
    end if
    dbgfile = addfile(cdf_debug_file,"c")
  else
;
; If this is not the first time, open the file as read/write.
;
    dbgfile = addfile(cdf_debug_file,"w")
  end if

;
; If this is not the first time, then we need to append the information.
;  
  if(.not.first_time) then
;
; The variables should be wrf_var_1, wrf_var_2, etc. We need to get the
; highest one already in use, so we can create the next highest one.
;
    wrf_debug_vars = getfilevarnames(dbgfile)
    max_num = max(stringtointeger(str_get_field(wrf_debug_vars,3,"_"))) + 1
  else
    max_num = 1
  end if

; This will be name of the logical variable to hold all the resources
; for this dataset.
  wrf_res_name  = "wrf_res_" + max_num

;
; Write the plot data to the netCDF file.  If the data contains the
; special 2D lat2d/lon2d arrays, we have to write these as 1D arrays
; and reconstruct them as 2D later.
;
  if(typeof(data1).ne."logical".and.typeof(data2).eq."logical") then
; For non u,v data
    wrf_data_name = "wrf_data_" + max_num
    add_latlon2d_debug_info(data1)
    dbgfile->$wrf_data_name$ = (/data1/)           ; Write the data
  end if
  if(typeof(data1).ne."logical".and.typeof(data2).ne."logical") then
; For u,v data
    add_latlon2d_debug_info(data1)
    add_latlon2d_debug_info(data2)
    wrf_data_name = "wrf_udata_" + max_num
    dbgfile->$wrf_data_name$ = (/data1/)           ; Write the U data
    wrf_data_name = "wrf_vdata_" + max_num
    dbgfile->$wrf_data_name$ = (/data2/)           ; Write the V data
  end if

; Retain the name of the wrf function that called this routine.
  tmp = "wrf_func_name_" + max_num
  dbgfile@$tmp$ = wrf_func_name

;
; Get plot resources, if any.
;
  wattnames = getvaratts(wrf_opts)
  if(.not.any(ismissing(wattnames))) then
    natt = dimsizes(wattnames)
  else
    natt = 0
  end if

;
; Check if any of the plot attributes are ones that can contain
; big data arrays, like sfXArray, vfYArray, etc.
;
; If so, then write these to the netCDF file. Otherwise, write them
; to the file as attributes.
;
  array_resources = (/"sfXArray","sfYArray","vfXArray","vfYArray"/)

  if(natt.gt.0) then
    tmp_opts = 1
    do i=0,natt-1
      if(any(wattnames(i).eq.array_resources)) then
        tmp = "wrf_data_coord_name_" + wattnames(i) + "_" + max_num
        dbgfile->$tmp$ = wrf_opts@$wattnames(i)$
      else
; Can't write "logical" to a NetCDF file.
        if(typeof(wrf_opts@$wattnames(i)$).eq."logical") then
          if(wrf_opts@$wattnames(i)$) then
            tmp_opts@$wattnames(i)$ = 1
          else
            tmp_opts@$wattnames(i)$ = 0
          end if
        else
; Just write the resource.
          tmp_opts@$wattnames(i)$ = wrf_opts@$wattnames(i)$
        end if
      end if
    end do
    dbgfile->$wrf_res_name$ = tmp_opts
  end if

; Close the NetCDF file
  delete(dbgfile)
end

undef("write_wrf_debug_script")
procedure write_wrf_debug_script(wks,debug_file,wrf_func_name)
begin
  dbgfile = addfile(debug_file+".nc","r")
  print(getvaratts(dbgfile))
end

undef("delete_attrs")
procedure delete_attrs(opts:logical)

; This procedure does some cleanup by removing unneeded attributes
; so they don't get passed to other routines by accident.

begin
  list_attrs = (/"MainTitle","MainTitlePos","MainTitlePosF", \
                 "InitTime","ValidTime","TimePos","TimePosF", \
                 "NoHeaderFooter","TimeLabel","LevelLabel", \
                 "FieldTitle","UnitLabel","NumVectors","AspectRatio", \
                 "SubFieldTitle","PlotOrientation","PlotLevelID", \
                 "mpNestTime","ContourParameters","FontHeightF","Footer", \
                 "start_lat","start_lon","end_lat","end_lon", \
                 "proj","map_proj","stand_lon","truelat1","truelat2","cenlat", \
                 "pole_lat","pole_lon","ref_lat","ref_lon","ref_x","ref_y", \
                 "e_we","e_sn","parent_id","parent_grid_ratio", \
                 "i_parent_start","j_parent_start", \
                 "dx","dy","max_dom" \
                 /)

  do i=0,dimsizes(list_attrs)-1
    if(isatt(opts,list_attrs(i))) then
      delete(opts@$list_attrs(i)$)
    end if
  end do
end

;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
undef("set_cn_resources")
function set_cn_resources (data[*][*]:numeric, res:logical)

begin

  opts = res

; The ContourParameters resource can either be a scalar that
; represents the contour level spacing, or it can be an array
; of three elements that represent the minimum level, the maximum
; level, and the level spacing.
;
  mx = max(data)
  mn = min(data)

  if(mn.ne.mx.and.opts.and.isatt(opts,"ContourParameters")) then
    if(dimsizes(opts@ContourParameters) .eq. 1) then

; Only the contour interval is specified.
      nlev   = tointeger((mx-mn)/opts@ContourParameters)+1
      levels = nice_mnmxintvl(mn,mx,nlev,True)
      if(levels(0) .lt. 0.) then
        ; Set a zero contour.
        nlev = tointeger(levels(0)/opts@ContourParameters) - 1
        levels(0) = nlev*opts@ContourParameters
      end if
      nlev = tointeger((levels(1)-levels(0))/opts@ContourParameters)+1
      levels(1) = levels(0) + nlev*opts@ContourParameters
      levels(2) = opts@ContourParameters

; Min level, max level, and level spacing are specified by user.
    else
      if(dimsizes(opts@ContourParameters) .eq. 3) then
        levels = opts@ContourParameters
      else
        print("wrf_contour: Warning: illegal setting for ContourParameters attribute")
      end if
    end if

  end if

; Contour levels
  if(isvar("levels")) then
    opts@cnLevelSelectionMode = get_res_value_keep(opts, "cnLevelSelectionMode", "ManualLevels")
    opts@cnMinLevelValF       = get_res_value_keep(opts, "cnMinLevelValF", levels(0))
    opts@cnMaxLevelValF       = get_res_value_keep(opts, "cnMaxLevelValF", levels(1))
    opts@cnLevelSpacingF      = get_res_value_keep(opts, "cnLevelSpacingF",levels(2))
    delete(levels)
  end if


; Set the default zero line thickness to 2, and the negative contour
; line dash pattern to 1 (0 is solid).
  opts@gsnContourZeroLineThicknessF = get_res_value_keep(opts, "gsnContourZeroLineThicknessF",2.0)
  opts@gsnContourNegLineDashPattern = get_res_value_keep(opts, "gsnContourNegLineDashPattern",1)


; Set resources that apply for both filled and line contour plots.
  opts@cnFillDrawOrder            = get_res_value_keep(opts,"cnFillDrawOrder", "PreDraw")


  opts@cnLineLabelAngleF          = get_res_value_keep(opts,"cnLineLabelAngleF", 0.0)
  opts@cnLineLabelFontHeightF     = get_res_value_keep(opts,"cnLineLabelFontHeightF", 0.015)
  opts@cnInfoLabelFontHeightF     = get_res_value_keep(opts,"cnInfoLabelFontHeightF", 0.015)
  opts@cnLineLabelPerimOn         = get_res_value_keep(opts,"cnLineLabelPerimOn", True)
  opts@cnInfoLabelPerimOn         = get_res_value_keep(opts,"cnInfoLabelPerimOn", False)
  opts@cnLineLabelBackgroundColor = get_res_value_keep(opts,"cnLineLabelBackgroundColor", -1)
  opts@cnHighLabelBackgroundColor = get_res_value_keep(opts,"cnHighLabelBackgroundColor", -1)
  opts@cnLowLabelBackgroundColor  = get_res_value_keep(opts,"cnLowLabelBackgroundColor", -1)
  opts@cnLineColor                = get_res_value_keep(opts,"cnLineColor", "Black")
  opts@cnLineLabelFontColor       = opts@cnLineColor
  opts@cnLineLabelPerimColor      = opts@cnLineColor
  opts@cnInfoLabelFontColor       = opts@cnLineColor
  opts@cnHighLabelFontColor       = opts@cnLineColor
  opts@cnLowLabelFontColor        = opts@cnLineColor


; Set field Title and levels if available
  if(.not.isatt(opts,"cnInfoLabelString")) then
    info_string = " Contours: $CMN$ to $CMX$ by $CIU$"
    if(isatt(opts,"FieldTitle")) then
        opts@cnInfoLabelString = opts@FieldTitle + info_string
    else if(isatt(data,"description")) then
        opts@cnInfoLabelString = data@description + info_string
    else
        opts@cnInfoLabelString = info_string
    end if
    end if
  end if


  return(opts)
end
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
undef("set_lb_resources")
function set_lb_resources (data[*][*]:numeric, res:logical)

begin

  opts = res


; Somewhat convoluted way to see if a labelbar is not desired.
  if(check_attr(opts,"pmTickMarkDisplayMode","Never",True).or.\
     check_attr(opts,"pmTickMarkDisplayMode",-1,False).or.\
     check_attr(opts,"pmTickMarkDisplayMode",0,False).or. \
     check_attr(opts,"lbLabelBarOn",False,False).or.\
     check_attr(opts,"lbLabelBarOn",0,False)) then
    lbar_on = False
  else
    lbar_on = True
  end if
  atmp = get_res_value(opts,"lbLabelBarOn",True) ; Remove this resource
  delete(atmp)                                   ; just in case.


; Possible title for the labelbar
  if(isatt(opts,"FieldTitle")) then
    lb_desc = opts@FieldTitle 
  else
    if(isatt(data,"description")) then
      lb_desc = data@description 
    else
      lb_desc = ""  
    end if
  end if

  if(isatt(opts,"UnitLabel") ) then
    lb_desc = lb_desc + "  (" + opts@UnitLabel + ")"
  else
    if(isatt(data,"units") .and. .not.(data@units.eq."")) then
      lb_desc = lb_desc + "  (" + data@units + ")"
    end if
  end if


  if(.not.isatt(opts,"cnFillColors")) then
    opts@gsnSpreadColors = get_res_value_keep(opts, "gsnSpreadColors", True)
  end if
  opts@cnInfoLabelOn   = get_res_value_keep(opts,"cnInfoLabelOn", False)
  opts@cnLinesOn       = get_res_value_keep(opts,"cnLinesOn", False)
  opts@cnLineLabelsOn  = get_res_value_keep(opts,"cnLineLabelsOn", False)

; Labelbar resources
  if(lbar_on) then
    opts@pmLabelBarDisplayMode    = get_res_value_keep(opts,"pmLabelBarDisplayMode", "Always")
    opts@pmLabelBarSide           = get_res_value_keep(opts,"pmLabelBarSide", "Bottom")
    opts@lbAutoManage             = get_res_value_keep(opts,"lbAutoManage",False)
    opts@lbOrientation            = get_res_value_keep(opts,"lbOrientation", "Horizontal")
    opts@lbPerimOn                = get_res_value_keep(opts,"lbPerimOn", False)
    opts@lbLabelJust              = get_res_value_keep(opts,"lbLabelJust", "BottomCenter")
    opts@lbLabelAutoStride        = get_res_value_keep(opts,"lbLabelAutoStride",True)
    opts@lbBoxMinorExtentF        = get_res_value_keep(opts,"lbBoxMinorExtentF", 0.13)
    opts@lbTitleFontHeightF       = get_res_value_keep(opts,"lbTitleFontHeightF", 0.015)
    opts@lbLabelFontHeightF       = get_res_value_keep(opts,"lbLabelFontHeightF", 0.015)
    opts@pmLabelBarOrthogonalPosF = get_res_value_keep(opts,"pmLabelBarOrthogonalPosF", -0.1)

    opts@lbTitleOn                = get_res_value_keep(opts,"lbTitleOn", True)
    if(lb_desc.ne."" .and. opts@lbTitleOn) then
      opts@lbTitleOn              = get_res_value_keep(opts,"lbTitleOn", True)
      opts@lbTitleString          = get_res_value_keep(opts,"lbTitleString", lb_desc)
      opts@lbTitleJust            = get_res_value_keep(opts,"lbTitleJust", "BottomCenter")
      opts@lbTitleOffsetF         = get_res_value_keep(opts,"lbTitleOffsetF", -0.5)
    else
      opts@lbTitleOn = False
    end if
  end if


  return(opts)
end
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
undef("set_title_resources")
function set_title_resources (data[*][*]:numeric, res:logical)

begin

  opts = res

; Set field Title and levels if available
  if(isatt(opts,"FieldTitle")) then
    SubTitles = opts@FieldTitle
  else
    if(isatt(data,"description")) then
      SubTitles = data@description
    else
      SubTitles = "UnKnown"
    end if
  end if

  if(isatt(opts,"SubFieldTitle")) then
    SubTitles = SubTitles + " " +  opts@SubFieldTitle
  end if
  if(isatt(opts,"UnitLabel")) then
    SubTitles = SubTitles + "   (" +  opts@UnitLabel + ")"
  else
     if(isatt(data,"units") .and. .not.(data@units.eq."")) then
      SubTitles = SubTitles + "   (" +  data@units + ")"
    end if
  end if
  if (isatt(opts,"PlotLevelID")) then
    SubTitles = SubTitles + "     at   " +  opts@PlotLevelID 
  else
    if (isatt(data,"PlotLevelID")) then
       SubTitles = SubTitles + "     at   " +  data@PlotLevelID 
    end if
  end if
  opts@tiMainString = SubTitles      

  return(opts)
end
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
undef("set_vc_resources")
function set_vc_resources (res:logical)

begin

  opts = res

  if ( isatt(opts,"vpWidthF") ) then
  ; num_vectors is used for vcMinDistanceF and vcRefLengthF
    width = opts@vpWidthF
    num_vectors = get_res_value(opts,"NumVectors",25.0)
    opts@vcMinDistanceF   = get_res_value_keep(opts,"vcMinDistanceF", width/num_vectors)
    opts@vcRefLengthF     = get_res_value_keep(opts,"vcRefLengthF", width/num_vectors)
  else
    opts@vcMinDistanceF   = get_res_value_keep(opts,"vcMinDistanceF", 0.02)
    opts@vcRefLengthF     = get_res_value_keep(opts,"vcRefLengthF", 0.02)
  end if


  opts@vcGlyphStyle     = get_res_value_keep(opts,"vcGlyphStyle", "WindBarb")
  opts@vcWindBarbColor  = get_res_value_keep(opts,"vcWindBarbColor", "Black")
  opts@vcRefAnnoOn      = get_res_value_keep(opts,"vcRefAnnoOn", False)
  opts@vcMinFracLengthF = get_res_value_keep(opts,"vcMinFracLengthF", .2)
 
  return(opts)
end
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
undef("set_mp_resources")
function set_mp_resources (res:logical)

begin

  opts = res

; "LowRes" is the default that NCL uses, so you don't need to
; set it here. However, if you want a higher resolution, use 
; "MediumRes". If you want higher resolution for the coastlines,
; then set it to "HighRes", but then you also need to download
; the RANGS-GSHHS database.  Higher resolutions take longer to
; draw.

      opts@mpDataBaseVersion     = get_res_value_keep(opts, "mpDataBaseVersion","MediumRes")
     ;opts@mpOutlineBoundarySets = get_res_value_keep(opts, "mpOutlineBoundarySets", "AllBoundaries")
      opts@mpOutlineBoundarySets = get_res_value_keep(opts, "mpOutlineBoundarySets", "GeophysicalAndUSStates")
      opts@mpPerimLineThicknessF = get_res_value_keep(opts, "mpPerimLineThicknessF", 1.0)
      opts@tmXBLabelFontHeightF  = get_res_value_keep(opts, "tmXBLabelFontHeightF", 0.01)
      opts@tmYLLabelFontHeightF  = get_res_value_keep(opts, "tmYLLabelFontHeightF", 0.01)

; Select portion of the map to view.
      opts@mpLimitMode           = get_res_value_keep(opts, "mpLimitMode","Corners")
      opts@mpLeftCornerLatF      = get_res_value_keep(opts, "mpLeftCornerLatF", opts@start_lat)
      opts@mpLeftCornerLonF      = get_res_value_keep(opts, "mpLeftCornerLonF", opts@start_lon)
      opts@mpRightCornerLatF     = get_res_value_keep(opts, "mpRightCornerLatF",opts@end_lat)
      opts@mpRightCornerLonF     = get_res_value_keep(opts, "mpRightCornerLonF",opts@end_lon)

      if ( opts@mpRightCornerLonF .lt. 0.0 ) then
        opts@mpRightCornerLonF  = opts@mpRightCornerLonF + 360.0
      end if 

; Set some other resources for line colors and grid spacing.

      opts@mpGeophysicalLineColor      = get_res_value_keep(opts, "mpGeophysicalLineColor","Gray")
      opts@mpGeophysicalLineThicknessF = get_res_value_keep(opts, "mpGeophysicalLineThicknessF",0.5)
      opts@mpGridLineColor             = get_res_value_keep(opts, "mpGridLineColor","Gray")
      opts@mpGridLineThicknessF        = get_res_value_keep(opts, "mpGridLineThicknessF",0.5)
     ;opts@mpGridMaskMode              = get_res_value_keep(opts, "mpGridMaskMode",3)
      opts@mpGridSpacingF              = get_res_value_keep(opts, "mpGridSpacingF",5)
      opts@mpLimbLineColor             = get_res_value_keep(opts, "mpLimbLineColor","Gray")
      opts@mpLimbLineThicknessF        = get_res_value_keep(opts, "mpLimbLineThicknessF",0.5)
      opts@mpNationalLineColor         = get_res_value_keep(opts, "mpNationalLineColor","Gray")
      opts@mpNationalLineThicknessF    = get_res_value_keep(opts, "mpNationalLineThicknessF",0.5)
      opts@mpPerimLineColor            = get_res_value_keep(opts, "mpPerimLineColor","Gray")
      opts@mpPerimOn                   = get_res_value_keep(opts, "mpPerimOn",True)
      opts@mpUSStateLineColor          = get_res_value_keep(opts, "mpUSStateLineColor","Gray")
      opts@mpUSStateLineThicknessF     = get_res_value_keep(opts, "mpUSStateLineThicknessF",0.5)
      opts@pmTickMarkDisplayMode       = get_res_value_keep(opts, "pmTickMarkDisplayMode","Always")

; Tick mark resources

     ;opts@tmXBMajorLengthF      = get_res_value(opts, "tmXBMajorLengthF",-0.03)
     ;opts@tmYLMajorLengthF      = get_res_value(opts, "tmYLMajorLengthF",-0.03)
      opts@tmXTOn                = get_res_value(opts,"tmXTOn",False)
      opts@tmYROn                = get_res_value(opts,"tmYROn",False)
      opts@tmYRLabelsOn          = get_res_value(opts,"tmYRLabelsOn",True)
      opts@tmXBBorderOn          = get_res_value(opts,"tmXBBorderOn",True)
      opts@tmXTBorderOn          = get_res_value(opts,"tmXTBorderOn",True)
      opts@tmYLBorderOn          = get_res_value(opts,"tmYLBorderOn",True)
      opts@tmYRBorderOn          = get_res_value(opts,"tmYRBorderOn",True)

  return(opts)
end
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------

undef("_SetMainTitle")
procedure _SetMainTitle(nc_file:file,wks[1]:graphic,cn[1]:graphic,opts)
 
; This procedure checks the input data for certain attributes, and
; based on those, sets MainTitle, InitTime and ValidTime
;
;  Attributes recognized by this procedure:
;     MainTitle                 (main title - top left)
;                               (with Init time top right)
;     TimeLabel                 (valid time  - right  under init time)
;     NoHeaderFooter            (switch all headers and footers off - mainly for panels)
;
; If the "NoHeaderFooter" attribute exists and is set True, then
; don't create any titles.

begin
;
  if(opts.and.isatt(opts,"NoHeaderFooter").and.opts@NoHeaderFooter) then
    return
  end if

;
; Set basic plot font
;
  font_height = get_res_value_keep(opts,"FontHeightF",0.01)
;
;
; If a MainTitle attribute hasn't been set, then set to "WRF"
; Also set an Initial time
;
; MAIN Header of plot
  opts@MainTitle    = get_res_value_keep(opts,"MainTitle", " ") 
  opts@MainTitlePos = get_res_value_keep(opts,"MainTitlePos", "Left") 
  opts@InitTime     = get_res_value_keep(opts,"InitTime", True) 
  opts@ValidTime    = get_res_value_keep(opts,"ValidTime", True) 
  opts@TimePos      = get_res_value_keep(opts,"TimePos", "Right") 
  opts@Footer       = get_res_value_keep(opts,"Footer", True) 


  if (opts@MainTitlePos .eq. "Left")
    opts@MainTitlePos  = "CenterLeft"
    opts@MainTitlePosF = 0.0
  end if
  if (opts@MainTitlePos .eq. "Center")
    opts@MainTitlePos  = "CenterCenter"
    opts@MainTitlePosF = 0.5
  end if
  if (opts@MainTitlePos .eq. "Right")
    opts@MainTitlePos  = "CenterRight"
    opts@MainTitlePosF = 1.0
  end if

  if (opts@TimePos .eq. "Left")
    MTOPosF = 0.30
  else
    MTOPosF = 0.20
  end if

  txt0 = create "MainPlotTitle" textItemClass wks
    "txString"      : opts@MainTitle
    "txFontHeightF" : font_height*1.5
  end create
  anno = NhlAddAnnotation(cn,txt0)
  setvalues anno
    "amZone"           : 3
    "amSide"           : "Top"
    "amJust"           : opts@MainTitlePos
    "amParallelPosF"   : opts@MainTitlePosF
    "amOrthogonalPosF" : MTOPosF
    "amResizeNotify"   : False
  end setvalues

; Time information on plot
  if (opts@TimePos .eq. "Left")
    opts@TimePos  = "CenterLeft"
    opts@TimePosF = 0.0
    if (opts@MainTitlePos .eq. "CenterLeft")
     MTOPosF = MTOPosF - 0.05
    end if
  end if
  if (opts@TimePos .eq. "Right")
    opts@TimePos  = "CenterRight"
    opts@TimePosF = 1.0
    if (opts@MainTitlePos .eq. "CenterRight")
     MTOPosF = MTOPosF - 0.05
    end if
  end if

  if( isatt(nc_file,"START_DATE") ) then
     model_start_time = nc_file@START_DATE
  else
     if( isatt(nc_file,"SIMULATION_START_DATE") ) then
       model_start_time = nc_file@SIMULATION_START_DATE
     else
       opts@InitTime = False
     end if
  end if   
  if( opts@InitTime ) then
    InitTime = "Init: " + model_start_time
    txt1 = create "InitTime" textItemClass wks
      "txString"      : InitTime
      "txFontHeightF" : font_height
    end create
    anno = NhlAddAnnotation(cn,txt1)
    setvalues anno
      "amZone"           : 3
      "amSide"           : "Top"
      "amJust"           : opts@TimePos
      "amParallelPosF"   : opts@TimePosF        
      "amOrthogonalPosF" : MTOPosF
      "amResizeNotify"   : False
    end setvalues
  end if

  plot_narrow = False
  if((opts).and.(isatt(opts,"vpWidthF")).and.(isatt(opts,"vpHeightF"))) then
     ph = opts@vpHeightF  
     pw = opts@vpWidthF  
     phw = ph/pw
     if ( phw .gt. 1.8 ) then
       plot_narrow = True
     end if
  end if

  if( opts@ValidTime .and. isatt(opts,"TimeLabel") ) then

    ValidTime = "Valid: " + opts@TimeLabel

    MTOPosF = MTOPosF - 0.03
    txt2 = create "ValidTime" textItemClass wks
      "txString"      : ValidTime
      "txFontHeightF" : font_height
    end create
    anno = NhlAddAnnotation(cn,txt2)
    setvalues anno
      "amZone"           : 3
      "amSide"           : "Top"
      "amJust"           : opts@TimePos
      "amParallelPosF"   : opts@TimePosF        
      "amOrthogonalPosF" : MTOPosF           
      "amResizeNotify"   : False
    end setvalues
  end if


; Add Footer if called for
  if( opts@Footer ) then
    footer1 = nc_file@TITLE
    dis =  nc_file@DX / 1000.0
    WE = "WEST-EAST_GRID_DIMENSION"
    SN = "SOUTH-NORTH_GRID_DIMENSION"
    BT = "BOTTOM-TOP_GRID_DIMENSION"
    footer2 =  " WE = " + nc_file@$WE$ + \
               " ; SN = " + nc_file@$SN$ + \
               " ; Levels = " + nc_file@$BT$ + \
               " ; Dis = " + dis + "km"
    if ( isatt(nc_file,"MP_PHYSICS")) then
      footer2 = footer2 +  " ; Phys Opt = " + nc_file@MP_PHYSICS 
    end if
    if ( isatt(nc_file,"BL_PBL_PHYSICS")) then
      footer2 = footer2 +  " ; PBL Opt = " + nc_file@BL_PBL_PHYSICS 
    end if
    if ( isatt(nc_file,"CU_PHYSICS")) then
      footer2 = footer2 +  " ; Cu Opt = " + nc_file@CU_PHYSICS 
    end if
    Footer = footer1 + "~C~" + footer2
   else
     Footer = " "
   end if
    txt3 = create "Footer" textItemClass wks
      "txString"      : Footer 
      "txFontHeightF" : font_height*.9
    end create
    anno = NhlAddAnnotation(cn,txt3)
    setvalues anno
      "amZone"           : 1
;      "amZone"           : 7
      "amJust"           : "TopLeft"
      "amSide"           : "Bottom"
      "amParallelPosF"   : 0.0
      "amOrthogonalPosF" : -0.55
      "amResizeNotify"   : False
    end setvalues


; Add X-setion information if needed
  if(opts.and.isatt(opts,"PlotOrientation")) then
    ;Xsection = "Cross-Section Orientation : " + opts@PlotOrientation
    Xsection =  opts@PlotOrientation
    txt4 = create "Xsection" textItemClass wks
      "txString"      : Xsection 
      "txFontHeightF" : font_height*.9
    end create
    anno = NhlAddAnnotation(cn,txt4)
    setvalues anno
      "amZone"           : 3
      "amSide"           : "Top"
      "amJust"           : "CenterRight"
      "amParallelPosF"   : 1.0
      "amOrthogonalPosF" : 0.005
      "amResizeNotify"   : False
    end setvalues
  end if

end

;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
;--------------------------------------------------------------------------------
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
        opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",in_file@CEN_LON)
        if( isatt(in_file,"POLE_LAT") ) then
          opts@mpCenterRotF = get_res_value_keep(opts, "mpCenterRotF", 90.0 - in_file@POLE_LAT) 
          delete(opts@mpCenterLonF)
          calcen = -190.
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF", calcen )
        end if
      end if
  
    end if

  return(opts)                                     ; Return.

end
;--------------------------------------------------------------------------------

undef("wrf_contour_ps")
function wrf_contour_ps(nc_file:file,wks[1]: graphic, data[*][*]:numeric, \
                     opt_args[1]:logical)

begin

  callFrame = True
  if ( isatt(opt_args,"FrameIT") ) then
    if ( .not.opt_args@FrameIT ) then
      callFrame = False
    end if
    delete (opt_args@FrameIT)
  end if

  lat2U = nc_file->XLAT_U(0,:,:)
  lon2U = nc_file->XLONG_U(0,:,:)

  opts = opt_args
  opts@sfXArray = lon2U
  opts@sfYArray = lat2U
  opts@sfDataArray = data
  opts@mpProjection         = "Stereographic"
  opts@mpEllipticalBoundary = True
  opts@mpFillOn = False
 
  opts@mpGeophysicalLineColor      = get_res_value_keep(opts, "mpGeophysicalLineColor","Gray")
  opts@mpGeophysicalLineThicknessF = get_res_value_keep(opts, "mpGeophysicalLineThicknessF",2.0)


  ; Set the contour resources
  opts = set_cn_resources(data,opts)
  opts@cnInfoLabelFontHeightF = 0.012
  opts@cnLineLabelPerimOn = False
  opts@cnInfoLabelPerimOn = True


  ; Find out if we are working with a contour or a shaded plot
  ; fill_on = False : line contour plot
  ; fill_on = True  : filled contour plot
  fill_on = get_res_value_keep(opts,"cnFillOn",False)
  if(fill_on) then      ; set lb resources if needed
    opts@pmLabelBarDisplayMode    = get_res_value_keep(opts,"pmLabelBarDisplayMode", "Never")
    opts = set_lb_resources(data,opts)
    opts@pmLabelBarOrthogonalPosF = 0.0
    opts@lbTitleJust = "BottomLeft"
  end if
  


  opts@gsnDraw = False   
  opts@gsnFrame = False
  opts@gsnMaximize = False
  delete_attrs(opts)
  cn = gsn_csm_contour_map_polar(wks,data,opts)          ; Create the plot.

  draw(cn)


  if ( callFrame ) then
    frame(wks)
  end if

return (cn)
end
;--------------------------------------------------------------------------------
undef("wrf_vector_ps")
function wrf_vector_ps(nc_file:file,wks[1]: graphic, \
                       data_u[*][*]:numeric, data_v[*][*]:numeric, \
                       opt_args[1]:logical)

begin

  callFrame = True
  if ( isatt(opt_args,"FrameIT") ) then
    if ( .not.opt_args@FrameIT ) then
      callFrame = False
    end if
    delete (opt_args@FrameIT)
  end if

  if(isfilevar(nc_file,"XLAT"))
    lat2T = nc_file->XLAT(0,:,:)
    lon2T = nc_file->XLONG(0,:,:)
  else
    lat2T = nc_file->XLAT_M(0,:,:)
    lon2T = nc_file->XLONG_M(0,:,:)
  end if

  opts = opt_args
  opts@vfXArray = lon2T
  opts@vfYArray = lat2T
  opts@vfUDataArray = data_u
  opts@vfVDataArray = data_v
  opts@mpProjection         = "Stereographic"
  opts@mpEllipticalBoundary = True
  opts@mpFillOn = False

  opts@mpGeophysicalLineColor      = get_res_value_keep(opts, "mpGeophysicalLineColor","Gray")
  opts@mpGeophysicalLineThicknessF = get_res_value_keep(opts, "mpGeophysicalLineThicknessF",2.0)


; Set vector resources
  opts = set_vc_resources(opts) 


 opts@gsnDraw = False   
 opts@gsnFrame = False
 opts@gsnMaximize = False
 delete_attrs(opts)
 cn = gsn_csm_vector_map_polar(wks,data_u,data_v,opts)          ; Create the plot.

  draw(cn)


  if ( callFrame ) then
    frame(wks)
  end if

return (cn)
end
;--------------------------------------------------------------------------------

undef("wrf_contour")
function wrf_contour(nc_file:file,wks[1]: graphic, data[*][*]:numeric, \
                     opt_args[1]:logical)
 
; This function creates a contour plot and adds some titles to it.
;
;   1. Determine width to height ratio of plot.
;
;   2. First determine if this is to be a filled or line
;      contour plot (fill_on)
;
;   3. If the ContourParameters attribute is set, then calculate
;      the contour levels.
;
;   4. Set two resources for setting the zero contour line to 
;      a larger thickness, and for changing the negative contour
;      lines to a dashed pattern.
;
;   5. If doing a filled contour plot, set a title for the labelbar
;      based on whether a units attribute is set.
;
;   6. Make a copy of the resource list, and set some additional
;      resources for filled contour plots.
;
;   7. Create the contour plot, attach the titles, and draw
;      and advance the frame (if requested).
 
local dims
begin
  opts = opt_args      ; Make a copy of the resource list.

  if(opts.and.isatt(opts,"gsnDebugWriteFileName")) then
    wrf_debug_file = get_res_value(opts,"gsnDebugWriteFileName",  "")    
  end if
    
  if(opts.and.isatt(opts,"mpOutlineBoundarySets")) then
    delete(opts@mpOutlineBoundarySets)
  end if


; Calculate ratio of plot width and height. Note that this doesn't
; affect the setting of gsnMaximize to True, because gsnMaximize will
; retain the aspect ratio of the plot.

  if(opts.and.isatt(opts,"AspectRatio")) then
    ratio = opts@AspectRatio
  else
    dims  = dimsizes(data)
    ratio = 1.*dims(0)/dims(1)
    if(ratio .gt. 1.2) then
      ratio = 1.2
    end if
    if(ratio .lt. 0.6667) then
      ratio = 0.6667
    end if
  end if

  if(ratio .gt. 1) 
    width  = 0.65 * 1.0/ratio
    height = 0.65
  else
    width  = 0.85 
    height = 0.85 * ratio
  end if

  opts@vpWidthF = get_res_value_keep(opts,"vpWidthF",  width)
  opts@vpHeightF = get_res_value_keep(opts,"vpHeightF", height)


; Set some basic contour resources 
  opts = set_cn_resources(data,opts)


; Find out if we are working with a contour or a shaded plot
; fill_on = False : line contour plot
; fill_on = True  : filled contour plot
  fill_on = get_res_value_keep(opts,"cnFillOn",False)
  if(fill_on) then      ; set lb resources if needed
    opts = set_lb_resources(data,opts)
    atmp = get_res_value(opts,"lbLabelBarOn",True) ; Remove this resource
    delete(atmp)                                   ; just in case.
  end if


; Set Title resources
  opts = set_title_resources(data,opts)


; Setting gsnScale to True ensures that the tickmark lengths and labels
; will be the same size on both axes.
  opts@gsnScale    = get_res_value_keep(opts,"gsnScale",    True)


; The default is not to draw the plot or advance the frame, and
; to maximize the plot in the frame.
  opts@gsnDraw          = False    ; Make sure don't draw or frame or,
  opts@gsnFrame         = False    ; maximize, b/c we'll do this later.
  opts@gsnMaximize      = False


  opts2 = opts
  delete_attrs(opts2)                       ; Clean up.
  cn = gsn_contour(wks,data,opts2)          ; Create the plot.
  _SetMainTitle(nc_file,wks,cn,opts)        ; Set some titles

  if(isvar("wrf_debug_file")) then
    write_wrf_debug_info(wks,data,False,wrf_debug_file,opts2,"wrf_contour")
  end if

  opts2@gsnDraw     = get_res_value_keep(opts2,"gsnDraw",     False)
  opts2@gsnFrame    = get_res_value_keep(opts2,"gsnFrame",    False)
  opts2@gsnMaximize = get_res_value_keep(opts2,"gsnMaximize", True)
  draw_and_frame(wks,cn,opts2@gsnDraw,opts2@gsnFrame,False,opts2@gsnMaximize)

  return(cn)                                ; Return

end

;--------------------------------------------------------------------------------

undef("wrf_vector")
function wrf_vector(nc_file:file,wks[1]: graphic, data_u[*][*]:numeric, \
                    data_v[*][*]:numeric, opt_args[1]:logical)
;
; This function creates a vector plot and adds some titles to it.
;
;   1. Determine width to height ratio of plot. Will also be use
;      to calculate values for vector resources later.
;
;   2. Make a copy of the resource list, and set some additional
;      resources.
;
;   3. Create the vector plot, attach the titles, and draw
;      and advance the frame (if requested).
 
local dims
begin
  opts = opt_args      ; Make a copy of the resource list.

  if(opts.and.isatt(opts,"gsnDebugWriteFileName")) then
    wrf_debug_file = get_res_value(opts,"gsnDebugWriteFileName",  "")    
  end if

  if(opts.and.isatt(opts,"mpOutlineBoundarySets")) then
    delete(opts@mpOutlineBoundarySets)
  end if
;
; The ratio is used to determine the width and height of the
; plot, and also to determine the value for the vcMinDistanceF
; resource.
;
  if(opts.and.isatt(opts,"AspectRatio")) then
    ratio = get_res_value(opts,"AspectRatio",0.)
  else
    dims  = dimsizes(data_u)
    ratio = 1.*dims(0)/dims(1)
    if(ratio .gt. 1.2) then
      ratio = 1.2
    end if
    if(ratio .lt. 0.6667) then
      ratio = 0.6667
    end if
  end if

  if(ratio .gt. 1) 
    width  = 0.65/ratio
    height = 0.65
  else
    width  = 0.95 
    height = 0.95 * ratio
  end if

  opts@vpWidthF  = get_res_value_keep(opts,"vpWidthF",  width)
  opts@vpHeightF = get_res_value_keep(opts,"vpHeightF", height)


; Set Title resources
  opts = set_title_resources(data_u,opts)


; Set vector resources
  opts = set_vc_resources(opts)


; Setting gsnScale to True ensures that the tickmark lengths and labels
; will be the same size on both axes.
  opts@gsnScale = get_res_value_keep(opts,"gsnScale",    True)


; The default is not to draw the plot or advance the frame, and
; to maximize the plot in the frame.
  opts@gsnDraw          = False    ; Make sure don't draw or frame or,
  opts@gsnFrame         = False    ; maximize, b/c we'll do this later.
  opts@gsnMaximize      = False

  opts2 = opts
  delete_attrs(opts2)                             ; Clean up.
  vct = gsn_vector(wks,data_u,data_v,opts2)       ; Create vector plot.
  _SetMainTitle(nc_file,wks,vct,opts) 

  if(isvar("wrf_debug_file")) then
    write_wrf_debug_info(wks,data_u,data_v,wrf_debug_file,opts2,"wrf_vector")
  end if

  opts2@gsnDraw     = get_res_value_keep(opts2,"gsnDraw",     False)
  opts2@gsnFrame    = get_res_value_keep(opts2,"gsnFrame",    False)
  opts2@gsnMaximize = get_res_value_keep(opts2,"gsnMaximize", True)
  draw_and_frame(wks,vct,opts2@gsnDraw,opts2@gsnFrame,False, \
                 opts2@gsnMaximize)

  return(vct)                                     ; Return.
end

;--------------------------------------------------------------------------------
undef("wrf_wps_map")
function wrf_wps_map(wks[1]:graphic,opt_args[1]:logical)   

begin
;
;   1. Make a copy of the resource list, and set some resources
;      common to all map projections.
;
;   2. Determine the projection being used, and set resources based
;      on that projection.
;
;   3. Create the map plot, and draw and advance the frame
;      (if requested).

    opts = opt_args      ; Make a copy of the resource list
    opts  = True

; Set some resources depending on what kind of map projection is 
; chosen.
;
;   MAP_PROJ = 0 : "CylindricalEquidistant"
;   MAP_PROJ = 1 : "LambertConformal"
;   MAP_PROJ = 2 : "Stereographic"
;   MAP_PROJ = 3 : "Mercator"
;   MAP_PROJ = 6 : "Lat/Lon"

;   CylindricalEquidistant
      if(opts@map_proj .eq. 0)
        projection          = "CylindricalEquidistant"
        opts@mpGridSpacingF = 45
        opts@mpCenterLatF   = get_res_value_keep(opts, "mpCenterLatF", 0.0)
        opts@mpCenterLonF   = get_res_value_keep(opts, "mpCenterLonF",opts@stand_lon)
      end if

;   LambertConformal projection
      if(opts@map_proj .eq. 1)
        projection               = "LambertConformal"
        opts@mpLambertParallel1F = get_res_value_keep(opts, "mpLambertParallel1F",opts@truelat1)
        opts@mpLambertParallel2F = get_res_value_keep(opts, "mpLambertParallel2F",opts@truelat2)
        opts@mpLambertMeridianF  = get_res_value_keep(opts, "mpLambertMeridianF",opts@stand_lon)
      end if
  
;   Stereographic projection
      if(opts@map_proj .eq. 2)
        projection          = "Stereographic"
        if( isatt(opts,"cenlat") ) then
          opts@mpCenterLatF  = get_res_value_keep(opts, "mpCenterLatF",opts@cenlat)
        else
          opts@mpCenterLatF  = get_res_value_keep(opts, "mpCenterLatF",opts@ref_lat)
        end if
        opts@mpCenterLonF   = get_res_value_keep(opts, "mpCenterLonF",opts@stand_lon)
      end if
  
;   Mercator projection
      if(opts@map_proj .eq. 3)
        projection          = "Mercator"
        opts@mpCenterLatF   = get_res_value_keep(opts, "mpCenterLatF", 0.0)
        opts@mpCenterLonF   = get_res_value_keep(opts, "mpCenterLonF",opts@stand_lon)
      end if
  
;    global WRF CylindricalEquidistant
      if(opts@map_proj .eq. 6)
        projection          = "CylindricalEquidistant"
        opts@mpGridSpacingF = 45
        if( isatt(opts,"cenlon") ) then
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",opts@cenlon)
        else
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF",opts@ref_lon)
        end if
        if( isatt(opts,"pole_lat") ) then
          delete(opts@mpCenterLonF)
          opts@mpCenterLonF  = get_res_value_keep(opts, "mpCenterLonF", - 190. )
          opts@mpCenterRotF = get_res_value_keep(opts, "mpCenterRotF", 90.0 - opts@pole_lat) 
        else
          opts@mpCenterRotF = get_res_value_keep(opts, "mpCenterRotF", 0.0) 
        end if
      end if

; Set some resources common to all map projections.
  opts = set_mp_resources(opts)

; The default is not to draw the plot or advance the frame, and
; to maximize the plot in the frame.

  opts@gsnDraw       = get_res_value_keep(opts,"gsnDraw",     False)
  opts@gsnFrame      = get_res_value_keep(opts,"gsnFrame",    False)
  opts@gsnMaximize   = get_res_value_keep(opts,"gsnMaximize", True)

  delete_attrs(opts)                             ; Clean up.
  mp = gsn_map(wks,projection,opts)              ; Create map plot.
  return(mp)                                     ; Return.

end

;--------------------------------------------------------------------------------
undef("wrf_wps_dom")
function wrf_wps_dom(wks[1]:graphic,opt_args[1]:logical,lnres[1]:logical,txres[1]:logical) 

begin

  mpres = opt_args
  ;BPR BEGIN
  ;grid_to_plot = 0 => plot using the corner mass grid points
  ;grid_to_plot = 1 => plot using the edges of the corner grid cells
  ;  This uses the locations of the corner mass grid points +/- 0.5 grid cells
  ;  to get from the mass grid point to the edge of the grid cell
  grid_to_plot = 1
  ;grid_to_plot = mpres@grid_to_plot
  if(grid_to_plot.eq.0) then
  ; print("Plotting using corner mass grid points")
  else if(grid_to_plot.eq.1) then
  ; print("Plotting using edges of the corner grid cells")
  else
   print("ERROR: Invalid value for grid_to_plot = "+grid_to_plot)
  end if
  end if
  ;BPR END

  res = True
  res@DX = mpres@dx
  res@DY = mpres@dy
  res@LATINC = 0.0
  res@LONINC = 0.0
  if ( mpres@map_proj .eq. "lambert") then
    mpres@map_proj = 1
    res@MAP_PROJ = 1
  end if
  if ( mpres@map_proj .eq. "polar") then
    mpres@map_proj = 2
    res@MAP_PROJ = 2
  end if
  if ( mpres@map_proj .eq. "mercator") then
    mpres@map_proj = 3
    res@MAP_PROJ = 3
  end if
  if ( mpres@map_proj .eq. "lat-lon") then
    mpres@map_proj = 6
    res@MAP_PROJ = 6
    res@LATINC = mpres@dy
    res@LONINC = mpres@dx
  end if
  res@TRUELAT1 = mpres@truelat1
  res@TRUELAT2 = mpres@truelat2
  res@STAND_LON = mpres@stand_lon

  res@REF_LAT = mpres@ref_lat
  res@REF_LON = mpres@ref_lon
  if ( isatt(mpres,"ref_x") ) then
    res@KNOWNI = mpres@ref_x
  else
    res@KNOWNI = int2flt(mpres@e_we(0))/2.
  end if
  if ( isatt(mpres,"ref_y") ) then
    res@KNOWNJ = mpres@ref_y
  else
    res@KNOWNJ = int2flt(mpres@e_sn(0))/2.
  end if

  if ( isatt(mpres,"pole_lat") ) then
    res@POLE_LAT = mpres@pole_lat
  else
    res@POLE_LAT = 90.0
  end if
  if ( isatt(mpres,"pole_lon") ) then
    res@POLE_LON = mpres@pole_lon
  else
    res@POLE_LON = 0.0
  end if

  ;BPR BEGIN
  ;Determine adjustment needed to convert from mass grid to chosen grid
  if(grid_to_plot.eq.0) then
   adjust_grid = 0.0 
  else if(grid_to_plot.eq.1) then
   adjust_grid = 0.5 
  else 
   print("ERROR: Invalid value for grid_to_plot = "+grid_to_plot)
   adjust_grid = 0.0 
  end if
  end if

  xx = 1.0 - adjust_grid
  yy = 1.0 - adjust_grid
  ;xx = 1.0
  ;yy = 1.0
  ;BPR END
  loc = wrf_ij_to_ll (xx,yy,res)
  start_lon = loc(0)
  start_lat = loc(1)
  ;BPR BEGIN
  ;e_we is the largest U grid point and e_sn the largest V gridpoint
  ;xx = int2flt(mpres@e_we(0))
  ;yy = int2flt(mpres@e_sn(0))
  ;Change it so it is in terms of mass grid points since wrf_ij_to_ll is
  ;in terms of mass grid points 
  xx = int2flt(mpres@e_we(0)-1) + adjust_grid
  yy = int2flt(mpres@e_sn(0)-1) + adjust_grid
  ;BPR END
  loc = wrf_ij_to_ll (xx,yy,res)
  end_lon = loc(0)
  end_lat = loc(1)

  mpres@start_lat   = start_lat
  mpres@start_lon   = start_lon
  mpres@end_lat     = end_lat
  mpres@end_lon     = end_lon


  mp = wrf_wps_map(wks,mpres)
  draw(mp)


  if ( mpres@max_dom .gt. 1 ) then

    numLineColors = 0
    if ( isatt(lnres,"domLineColors") ) then
       numLineColors = dimsizes(lnres@domLineColors)
    end if

    do idom = 1,mpres@max_dom-1

      if ( numLineColors .gt. 0 ) then
        if ( numLineColors .ge. idom ) then
          lnres@gsLineColor = lnres@domLineColors(idom-1)
          txres@txFontColor = lnres@domLineColors(idom-1)
        else
          lnres@gsLineColor = lnres@domLineColors(numLineColors-1)
          txres@txFontColor = lnres@domLineColors(numLineColors-1)
        end if
      end if


      ; nest start and end points in large domain space
      if ( mpres@parent_id(idom) .eq. 1) then
        ; corner value
        ;BPR BEGIN
        ;Due to the alignment of nests we need goffset in order to
        ;find the location of (1,1) in the fine domain in coarse domain
        ;coordinates
        ;i_start = mpres@i_parent_start(idom)
        ;j_start = mpres@j_parent_start(idom)
        goffset = 0.5*(1-(1.0/mpres@parent_grid_ratio(idom)))
        i_start = mpres@i_parent_start(idom)-goffset
        j_start = mpres@j_parent_start(idom)-goffset
        ; end point
        ;Change to mass point
        ;i_end = (mpres@e_we(idom)-1)/mpres@parent_grid_ratio(idom) + i_start
        ;j_end = (mpres@e_sn(idom)-1)/mpres@parent_grid_ratio(idom) + j_start
        i_end = (mpres@e_we(idom)-2)/(1.0*mpres@parent_grid_ratio(idom)) + i_start
        j_end = (mpres@e_sn(idom)-2)/(1.0*mpres@parent_grid_ratio(idom)) + j_start
        
        if(grid_to_plot.eq.0) then
         adjust_grid = 0.0 
        else if(grid_to_plot.eq.1) then
         adjust_grid = 0.5/(1.0*mpres@parent_grid_ratio(idom))
        else 
         print("ERROR: Invalid value for grid_to_plot = "+grid_to_plot)
         adjust_grid = 0.0 
        end if
        end if
        
        ;BPR END
      end if
      if ( mpres@parent_id(idom) .ge. 2) then
        ; corner value
        nd = mpres@parent_id(idom)
        ;BPR BEGIN
        ;i_points = ((mpres@e_we(idom)-1)/mpres@parent_grid_ratio(idom))
        ;j_points = ((mpres@e_sn(idom)-1)/mpres@parent_grid_ratio(idom))
        i_points = ((mpres@e_we(idom)-2)/(1.0*mpres@parent_grid_ratio(idom)))
        j_points = ((mpres@e_sn(idom)-2)/(1.0*mpres@parent_grid_ratio(idom)))
        goffset = 0.5*(1-(1.0/(1.0*mpres@parent_grid_ratio(idom))))
        ai_start = mpres@i_parent_start(idom)*1.0-goffset
        aj_start = mpres@j_parent_start(idom)*1.0-goffset
        ;ai_start = mpres@i_parent_start(idom)*1.0
        ;aj_start = mpres@j_parent_start(idom)*1.0

        if(grid_to_plot.eq.0) then
         adjust_grid = 0.0 
        else if(grid_to_plot.eq.1) then
         adjust_grid = 0.5/(1.0*mpres@parent_grid_ratio(idom)) 
        else 
         print("ERROR: Invalid value for grid_to_plot = "+grid_to_plot)
         adjust_grid = 0.0 
        end if
        end if

        do while ( nd .gt. 1)
          ;Note that nd-1 is used in the following because the WPS namelist is
          ;one-based but arrays in NCL are zero-based
          goffset = 0.5*(1-(1.0/(1.0*mpres@parent_grid_ratio(nd-1))))
          ;ai_start = ai_start/mpres@parent_grid_ratio(nd-1) + mpres@i_parent_start(nd-1)
          ;aj_start = aj_start/mpres@parent_grid_ratio(nd-1) + mpres@j_parent_start(nd-1)
          ai_start = (ai_start-1)/(1.0*mpres@parent_grid_ratio(nd-1)) + mpres@i_parent_start(nd-1)-goffset
          aj_start = (aj_start-1)/(1.0*mpres@parent_grid_ratio(nd-1)) + mpres@j_parent_start(nd-1)-goffset
          ;i_points = (i_points/mpres@parent_grid_ratio(nd-1))
          ;j_points = (j_points/mpres@parent_grid_ratio(nd-1))
          i_points = (i_points/(1.0*mpres@parent_grid_ratio(nd-1)))
          j_points = (j_points/(1.0*mpres@parent_grid_ratio(nd-1)))
          if(grid_to_plot.eq.0) then
           adjust_grid = 0.0 
          else if(grid_to_plot.eq.1) then
           adjust_grid = adjust_grid/(1.0*mpres@parent_grid_ratio(nd-1))
          else 
           print("ERROR: Invalid value for grid_to_plot = "+grid_to_plot)
           adjust_grid = 0.0 
          end if
          end if

          ;nd = nd - 1
          nd = mpres@parent_id(nd-1)
        end do
        ;i_start = tointeger(ai_start + .5 )
        ;j_start = tointeger(aj_start + .5 )
        i_start = ai_start
        j_start = aj_start
        ; end point
        ;i_end = i_points + i_start + 1
        ;j_end = j_points + j_start + 1
        i_end = i_points + i_start
        j_end = j_points + j_start
        ;BPR END
      end if

      ; get the four corners

      xx = i_start - adjust_grid
      yy = j_start - adjust_grid
      ;xx = int2flt(i_start)
      ;yy = int2flt(j_start)
      loc = wrf_ij_to_ll (xx,yy,res)
      lon_SW = loc(0)
      lat_SW = loc(1)

      xx = i_end + adjust_grid
      yy = j_start - adjust_grid
      ;xx = int2flt(i_end)
      ;yy = int2flt(j_start)
      loc = wrf_ij_to_ll (xx,yy,res)
      lon_SE = loc(0)
      lat_SE = loc(1)

      xx = i_start - adjust_grid
      yy = j_end + adjust_grid
      ;xx = int2flt(i_start)
      ;yy = int2flt(j_end)
      loc = wrf_ij_to_ll (xx,yy,res)
      lon_NW = loc(0)
      lat_NW = loc(1)

      ;xx = int2flt(i_end)
      ;yy = int2flt(j_end)
      xx = i_end + adjust_grid
      yy = j_end + adjust_grid
      ;BPR END
      loc = wrf_ij_to_ll (xx,yy,res)
      lon_NE = loc(0)
      lat_NE = loc(1)

      xbox = (/lon_SW, lon_SE, lon_NE, lon_NW, lon_SW /)
      ybox = (/lat_SW, lat_SE, lat_NE, lat_NW, lat_SW /)
      x_out = new(dimsizes(xbox),typeof(xbox))
      y_out = new(dimsizes(ybox),typeof(ybox))
      datatondc(mp, xbox, ybox, x_out, y_out)
      gsn_polyline_ndc(wks, x_out, y_out, lnres)

      idd = idom + 1
      dom_text = "d0"+idd
      if ( txres@txJust .eq. "BottomLeft" ) then
        gsn_text(wks,mp,dom_text,lon_NW,lat_NW,txres)
      else
        gsn_text_ndc(wks,dom_text,x_out(3)+0.01,y_out(3)-0.01,txres)
      end if

    end do

  end if

  return(mp)

end

;--------------------------------------------------------------------------------
undef("wrf_map_resources")
function wrf_map_resources(in_file[1]:file,map_args[1]:logical)   
local lat, lon, x1, x2, y1, y2, dims, ii, jj
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
        map_args@mpProjection          = "CylindricalEquidistant"
        map_args@mpGridSpacingF = 45
        map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF",in_file@CEN_LON)
        if( isatt(in_file,"POLE_LAT") ) then
          map_args@mpCenterRotF = get_res_value_keep(map_args, "mpCenterRotF", 90.0 - in_file@POLE_LAT) 
          delete(map_args@mpCenterLonF)
          calcen = -190.
          map_args@mpCenterLonF  = get_res_value_keep(map_args, "mpCenterLonF", calcen )
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

        delete(map_args@ZoomIn)
      end if

      return(map_args)
end

;--------------------------------------------------------------------------------
undef("wrf_map")
function wrf_map(wks[1]:graphic,in_file[1]:file,opt_args[1]:logical)   

begin
;
; This function creates a map plot, and bases the projection on
; the MAP_PROJ attribute in the given file.
;
;   1. Make a copy of the resource list, and call a function to set
;      some resources common to all map projections.
;
;   2. Determine the projection being used, and set resources based
;      on that projection.
;
;   3. Create the map plot, and draw and advance the frame
;      (if requested).

    if(opt_args.and.isatt(opt_args,"gsnDebugWriteFileName")) then
      wrf_debug_file = get_res_value(opt_args,"gsnDebugWriteFileName",  "")    
    end if

    opts = opt_args      ; Make a copy of the resource list
    opts  = True

;---Set some map resources based on parameters and variables in input file.
    opts = wrf_map_resources(in_file,opts)

    if(.not.isatt(opts,"mpProjection")) then
      print("wrf_map: Error: no MAP_PROJ attribute in input file")
      return(new(1,graphic))
    else
      projection = opts@mpProjection
    end if

; The default is not to draw the plot or advance the frame, and
; to maximize the plot in the frame.

  opts@gsnDraw       = get_res_value_keep(opts,"gsnDraw",     False)
  opts@gsnFrame      = get_res_value_keep(opts,"gsnFrame",    False)
  opts@gsnMaximize   = get_res_value_keep(opts,"gsnMaximize", True)

  delete_attrs(opts)                             ; Clean up.
  mp = gsn_map(wks,projection,opts)              ; Create map plot.

  if(isvar("wrf_debug_file")) then
    opts@mpProjection = projection
    write_wrf_debug_info(wks,False,False,wrf_debug_file,opts,"wrf_map")
  end if

  return(mp)                                     ; Return.

end

;--------------------------------------------------------------------------------

undef("wrf_map_overlays")
function wrf_map_overlays(in_file[1]:file, \
                          wks:graphic, \
                          plots[*]:graphic, \
                          opt_arg[1]:logical, \
                          opt_mp[1]:logical) 
 
; This procedure takes an array of plots and overlays them on a
; base plot - map background.
;
; It will advance the plot and cleanup, unless you set the
; PanelPlot resource to True.
;
;  Attributes recognized by this procedure:
;     FramePlot
;     PanelPlot
;     NoTitles                  (don't do any titles) 
;     CommonTitle & PlotTitle is used to overwrite field titles
;        CommonTitle will supercede NoTitles
;     LatLonOverlay
;
; If FramePlot False, then Draw the plot but do not Frame.
; In this case a user want to add to the drawing, and will
; have to advance the Frame manually in the script.
;
; If the "NoTitles" attribute exists and is set True, then
; don't create the top-left titles, and leave the main titles alone.
; This resource can be useful if you are planning to panel
; the plots.
;
; If PanelPlot is set to True, then this flags to wrf_map_overlays
; that these plots are going to be eventually paneled (likely
; by gsn_panel), and hence 1) draw and frame should not be called
; (unless gsnDraw and/or gsnFrame are explicitly set to True),
; and 2) the overlays and titles should not be removed with
; NhlRemoveOverlay and NhlRemoveAnnotation.
;
; If LatLonOverlay is set to True, then this means the user is
; using the 2D lat/lon coordinates to do the overlay, and hence
; tfDoNDCOverlay should not be set to True. (The default is False.)
;
begin

  opts     = opt_arg         ; Make a copy of the resource lists
  opt_mp_2 = opt_mp      

  if(opts.and.isatt(opts,"gsnDebugWriteFileName")) then
    opt_mp_2@gsnDebugWriteFileName = get_res_value(opts, \
                                            "gsnDebugWriteFileName",  "")    
    wrf_debug_file = opt_mp_2@gsnDebugWriteFileName
  end if

  ; Let's make the map first
  base = wrf_map(wks,in_file,opt_mp_2)

  no_titles  = get_res_value(opts,"NoTitles",False)     ; Do we want field titles?
  com_title  = get_res_value(opts,"CommonTitle",False)     ; Do we have a common title?
  if ( com_title ) then
    plot_title = get_res_value(opts,"PlotTitle","  ")
    no_titles = True
  end if
  
  call_draw  = True
  call_frame = get_res_value(opts,"FramePlot",True)     ; Do we want to frame the plot?
  panel_plot = get_res_value(opts,"PanelPlot",False)    ; Are we paneling?
  latlon_overlay = get_res_value(opts,"LatLonOverlay",False) ; Lat/lon Overlay?
  opts@gsnMaximize = get_res_value_keep(opts,"gsnMaximize", True)

  nplots = dimsizes(plots)
;  font_color = "Black"

  do i=0,nplots-1
    if(.not.ismissing(plots(i))) then
;      class_name = NhlClassName(plots(i))
;      print(class_name)
;      if(class_name.eq."contourPlotClass") then
;        getvalues plots(i)
;          "cnFillOn"    : fill_on
;          "cnLineColor" : line_color
;        end getvalues
;        if (.not.fill_on) then
;          font_color = line_color      
;        end if       
;      end if
      if(.not.no_titles) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        if(i.eq.0) then
          SubTitles = SubTitle
        else
          SubTitles = SubTitles + "~C~" + SubTitle
        end if
      end if
      if(com_title .and. i .eq. nplots-1) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        SubTitles = plot_title
      end if
      if(.not.latlon_overlay) then
        setvalues plots(i)
          "tfDoNDCOverlay" : True
          "tiMainOn"       : False
        end setvalues
      else
        setvalues plots(i)
          "tiMainOn"       : False
        end setvalues
      end if
      overlay(base,plots(i))
    else
      print("wrf_map_overlays: Warning: overlay plot #" + i + " is not valid.")
    end if
  end do

  if(.not.no_titles .or. com_title) then
    font_height = get_res_value_keep(opts,"FontHeightF",0.01)
    txt = create "map_titles" textItemClass wks
      "txString"      : SubTitles
      "txFontHeightF" : font_height
     ;"txFontColor"   : font_color
    end create
    anno = NhlAddAnnotation(base,txt)
    setvalues anno
      "amZone"           : 3
      "amJust"           : "BottomLeft"
      "amSide"           : "Top"
      "amParallelPosF"   : 0.005
      "amOrthogonalPosF" : 0.03
      "amResizeNotify"   : False
    end setvalues
    base@map_titles = anno
  end if
;
; gsnDraw and gsnFrame default to False if panel plot.
;
  if(panel_plot) then
    call_draw = False
    call_frame= False
  end if


  opts@gsnDraw     = get_res_value_keep(opts,"gsnDraw",     call_draw)
  opts@gsnFrame    = get_res_value_keep(opts,"gsnFrame",    call_frame)

  draw_and_frame(wks,base,opts@gsnDraw,opts@gsnFrame,False, \
                 opts@gsnMaximize)

  if(.not.panel_plot) then
    do i=0,nplots-1
      if(.not.ismissing(plots(i))) then
        NhlRemoveOverlay(base,plots(i),False)
      else
        print("wrf_remove_map_overlays: Warning: overlay plot #" + i + " is not valid.")
        print("                        Nothing to remove.")
      end if
    end do
  end if

  if(isvar("wrf_debug_file")) then
    write_wrf_debug_script(wks,wrf_debug_file,"wrf_map_overlays")
  end if

  if(.not.no_titles.and..not.panel_plot) then
    if(isatt(base,"map_titles")) then
      NhlRemoveAnnotation(base,base@map_titles)
      delete(base@map_titles)
    end if
  end if

return(base)
end

;--------------------------------------------------------------------------------

undef("wrf_overlays")
function wrf_overlays(in_file[1]:file, \
                     wks:graphic, plots[*]:graphic, \
                     opt_arg[1]:logical)
 
; This procedure takes an array of plots and overlays them.
;
; It will advance the plot and cleanup, unless you set the
; PanelPlot resource to True.
; 
;  Attributes recognized by this procedure:
;     FramePlot
;     PanelPlot
;     NoTitles                  (don't do any titles) 
;     CommonTitle & PlotTile is used to overwrite field titles
;        CommonTitle will super-seed NoTitles
;
; If FramePlot False, then Draw the plot but do not Frame.
; In this case a user want to add to the drawing, and will
; have to advance the Frame manually in the script.
;
; If the "NoTitles" attribute exists and is set True, then
; don't create the top-left titles, and leave the main titles alone.
; This resource can be useful if you are planning to panel
; the plots.
;
; If PanelPlot is set to True, then this flags to wrf_overlays
; that these plots are going to be eventually paneled (likely
; by gsn_panel), and hence 1) draw and frame should not be called
; (unless gsnDraw and/or gsnFrame are explicitly set to True),
; and 2) the overlays and titles should not be removed with
; NhlRemoveOverlay and NhlRemoveAnnotation.
;
; If LatLonOverlay is set to True, then this means the user is
; using the 2D lat/lon coordinates to do the overlay, and hence
; tfDoNDCOverlay should not be set to True. (The default is False.)
;
begin
  opts = opt_arg      ; Make a copy of the resource list.

  no_titles  = get_res_value(opts,"NoTitles",False)     ; Do we want field titles?
  com_title  = get_res_value(opts,"CommonTitle",False)     ; Do we have a common title?
  latlon_overlay = get_res_value(opts,"LatLonOverlay",False) ; Lat/lon Overlay?
  if ( com_title ) then
    plot_title = get_res_value(opts,"PlotTitle","  ")
    no_titles = True
  end if

  call_draw  = True
  call_frame = get_res_value(opts,"FramePlot",True)     ; Do we want to frame the plot?
  panel_plot = get_res_value(opts,"PanelPlot",False)    ; Are we paneling?
  opts@gsnMaximize = get_res_value_keep(opts,"gsnMaximize", True)

  nplots = dimsizes(plots)
    
  base = plots(0)
  if(.not.no_titles) then
    getvalues plots(0)
      "tiMainString" : SubTitle
    end getvalues
    SubTitles = SubTitle
    if(.not.latlon_overlay) then
      setvalues plots(0)
        "tfDoNDCOverlay" : True
        "tiMainOn"       : False
      end setvalues
    else
      setvalues plots(0)
        "tiMainOn"       : False
      end setvalues
    end if
  else
    if(.not.latlon_overlay) then
      setvalues plots(0)
        "tfDoNDCOverlay" : True
      end setvalues
    end if
  end if

  if (nplots.eq.1) then 
    blank = create "BlankPlot" logLinPlotClass wks
      ;"cnConstFLabelOn" : False
    end create
    overlay(base,blank)
  end if

  do i=1,nplots-1
    if(.not.ismissing(plots(i))) then
      if(.not.no_titles) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        if(i.eq.0) then
          SubTitles = SubTitle
        else
          SubTitles = SubTitles + "~C~" + SubTitle
        end if
      end if
      if(com_title .and. i .eq. nplots-1) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        SubTitles = plot_title
      end if
      if(.not.latlon_overlay) then
        setvalues plots(i)
          "tfDoNDCOverlay" : True
          "tiMainOn"       : False
        end setvalues
      else
        setvalues plots(i)
          "tiMainOn"       : False
        end setvalues
      end if
      overlay(base,plots(i))
    else
      print("wrf_overlays: Warning: overlay plot #" + i + " is not valid.")
    end if
  end do

  if(.not.no_titles .or. com_title) then
    font_height = get_res_value_keep(opts,"FontHeightF",0.01)

    txt = create "map_titles" textItemClass wks
      "txString"      : SubTitles
      "txFontHeightF" : font_height
    end create
    anno = NhlAddAnnotation(base,txt)
    setvalues anno
      "amZone"           : 3
      "amJust"           : "BottomLeft"
      "amSide"           : "Top"
      "amParallelPosF"   : 0.005
      "amOrthogonalPosF" : 0.03
      "amResizeNotify"   : False
    end setvalues
    base@map_titles = anno
  end if

;
; gsnDraw and gsnFrame should default to True if not a panel plot.
;  
  if(panel_plot) then
    call_draw = False
    call_frame= False
  end if

  opts@gsnDraw     = get_res_value_keep(opts,"gsnDraw",     call_draw)
  opts@gsnFrame    = get_res_value_keep(opts,"gsnFrame",    call_frame)
  opts@gsnMaximize = get_res_value_keep(opts,"gsnMaximize", True)

  draw_and_frame(wks,base,opts@gsnDraw,opts@gsnFrame,False, \
                 opts@gsnMaximize)

  if(.not.no_titles.and..not.panel_plot) then
    NhlRemoveAnnotation(base,base@map_titles)
    delete(base@map_titles)
  end if

  if(.not.panel_plot) then
    if ( nplots .ge. 2 ) then
      do i=1,nplots-1
        if(.not.ismissing(plots(i))) then
          NhlRemoveOverlay(base,plots(i),False)
        else
          print("wrf_remove_overlays: Warning: overlay plot #" + i + " is not valid.")
          print("                        Nothing to remove.")
        end if
      end do
    end if
  end if

return(base)
end

;--------------------------------------------------------------------------------

undef("wrf_map_zoom")
function wrf_map_zoom(wks[1]:graphic,in_file[1]:file,opt_args[1]:logical, \
                      y1:integer,y2:integer,x1:integer,x2:integer)

; As of version 5.0.1, this routine is redundant. Use the special "ZoomIn"
; resource in wrf_map to accomplish the same thing. This function is
; being kept for backwards capability. There should be no need for it
; except to run old WRF-NCL codes.  Do not make any changes to it except
; possibly to fix bugs.
;
begin
      print("wrf_map_zoom: Warning: This function is obsolete. Consider using")
      print("                       the 'ZoomIn' resource in wrf_map instead.")
      if(isfilevar(in_file,"XLAT"))
        lat = in_file->XLAT(0,:,:)
        lon = in_file->XLONG(0,:,:)
      else
        lat = in_file->XLAT_M(0,:,:)
        lon = in_file->XLONG_M(0,:,:)
      end if
      opts  = opt_args      ; Make a copy of the resource list
      opts  = True
      opts@mpLeftCornerLatF      = lat(y1,x1)
      opts@mpLeftCornerLonF      = lon(y1,x1)
      opts@mpRightCornerLatF     = lat(y2,x2)
      opts@mpRightCornerLonF     = lon(y2,x2)
      mz = wrf_map(wks,in_file,opts)
      return(mz)
end

;--------------------------------------------------------------------------------

undef("wrf_map_overlay")
procedure wrf_map_overlay(wks:graphic,base[1]:graphic, \
                          plots[*]:graphic, \
                          opt_arg[1]:logical)
 
; As of version 5.0.1, this procedure is obsolete. Use wrf_map_overlays
; instead. It is being kept for backwards capability. Do not make any
; changes to it except possibly to fix bugs.
;
; This procedure takes an array of plots and overlays them on a
; base plot - map background.
;
; It will advance the plot and cleanup, unless you set the
; PanelPlot resource to True.
;
;  Attributes recognized by this procedure:
;     NoTitles                  (don't do any titles) 
;     PanelPlot
;
; If the "NoTitles" attribute exists and is set True, then
; don't create the top-left titles, and leave the main titles alone.
; This resource can be useful if you are planning to panel
; the plots.
;
; If PanelPlot is set to True, then this flags to wrf_map_overlay
; that these plots are going to be eventually paneled (likely
; by gsn_panel), and hence 1) draw and frame should not be called
; (unless gsnDraw and/or gsnFrame are explicitly set to True),
; and 2) the overlays and titles should not be removed with
; NhlRemoveOverlay and NhlRemoveAnnotation.
;
begin
  print("wrf_map_overlay: Warning: This procedure is obsolete. Consider" )
  print("                          using wrf_map_overlays instead.")

  opts = opt_arg      ; Make a copy of the resource list
  opts = True

  if(opts.and.isatt(opts,"NoTitles").and.opts@NoTitles) then
    no_titles = True
  else
    no_titles = False
  end if

  panel_plot = get_res_value(opts,"PanelPlot",False)    ; Are we paneling?
  latlon_overlay = get_res_value(opts,"LatLonOverlay",False) ; Lat/lon Overlay?

  nplots = dimsizes(plots)
;  font_color = "Black"

  do i=0,nplots-1
    if(.not.ismissing(plots(i))) then
;      class_name = NhlClassName(plots(i))
;      print(class_name)
;      if(class_name.eq."contourPlotClass") then
;        getvalues plots(i)
;          "cnFillOn"    : fill_on
;          "cnLineColor" : line_color
;        end getvalues
;        if (.not.fill_on) then
;          font_color = line_color      
;        end if       
;      end if
      if(.not.no_titles) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        if(i.eq.0) then
          SubTitles = SubTitle
        else
          SubTitles = SubTitles + "~C~" + SubTitle
        end if
        if(.not.latlon_overlay) then
    	  setvalues plots(i)
            "tfDoNDCOverlay" : True
            "tiMainOn"       : False
          end setvalues
        else
    	  setvalues plots(i)
            "tiMainOn"       : False
          end setvalues
        end if
      else
        if(.not.latlon_overlay) then
  	  setvalues plots(i)
            "tfDoNDCOverlay" : True
          end setvalues
        end if
      end if
      overlay(base,plots(i))
    else
      print("wrf_map_overlay: Warning: overlay plot #" + i + " is not valid.")
    end if
  end do

  if(.not.no_titles) then
    font_height = get_res_value_keep(opts,"FontHeightF",0.01)
    txt = create "map_titles" textItemClass wks
      "txString"      : SubTitles
      "txFontHeightF" : font_height
     ;"txFontColor"   : font_color
    end create
    anno = NhlAddAnnotation(base,txt)
    setvalues anno
      "amZone"           : 3
      "amJust"           : "BottomLeft"
      "amSide"           : "Top"
      "amParallelPosF"   : 0.005
      "amOrthogonalPosF" : 0.03
      "amResizeNotify"   : False
    end setvalues
    base@map_titles = anno
  end if
;
; gsnDraw and gsnFrame should default to True if not a panel plot.
; gsnFrame will default to False if opt_arg is False.
;
  if(panel_plot.or..not.opt_arg) then
    call_frame = False
  else
    call_frame = True
  end if
  if(panel_plot) then
    call_draw = False
  else
    call_draw = True
  end if
  opts@gsnDraw     = get_res_value_keep(opts,"gsnDraw",     call_draw)
  opts@gsnFrame    = get_res_value_keep(opts,"gsnFrame",    call_frame)
  opts@gsnMaximize = get_res_value_keep(opts,"gsnMaximize", True)

  draw_and_frame(wks,base,opts@gsnDraw,opts@gsnFrame,False, \
                 opts@gsnMaximize)

  if(.not.panel_plot) then
    do i=0,nplots-1
      if(.not.ismissing(plots(i))) then
        NhlRemoveOverlay(base,plots(i),False)
      else
        print("wrf_remove_map_overlay: Warning: overlay plot #" + i + " is not valid.")
        print("                        Nothing to remove.")
      end if
    end do
  end if

  if(.not.no_titles.and..not.panel_plot) then
    if(isatt(base,"map_titles")) then
      NhlRemoveAnnotation(base,base@map_titles)
      delete(base@map_titles)
    end if
  end if
end

;--------------------------------------------------------------------------------

undef("wrf_overlay")
procedure wrf_overlay(wks:graphic, plots[*]:graphic, \
                          opt_arg[1]:logical)
 
; As of version 5.0.1, this procedure is obsolete. Use wrf_overlays
; instead. It is being kept for backwards capability. Do not make any
; changes to it except possibly to fix bugs.
;
; This procedure takes an array of plots and overlays them.
;
; It will advance the plot and cleanup, unless you set the
; PanelPlot resource to True.
;
;  Attributes recognized by this procedure:
;     NoTitles                  (don't do any titles) 
;     PanelPlot
;
; If the "NoTitles" attribute exists and is set True, then
; don't create the top-left titles, and leave the main titles alone.
; This resource can be useful if you are planning to panel
; the plots.
;
; If PanelPlot is set to True, then this flags to wrf_overlay
; that these plots are going to be eventually paneled (likely
; by gsn_panel), and hence 1) draw and frame should not be called
; (unless gsnDraw and/or gsnFrame are explicitly set to True),
; and 2) the overlays and titles should not be removed with
; NhlRemoveOverlay and NhlRemoveAnnotation.
;
; If LatLonOverlay is set to True, then this means the user is
; using the 2D lat/lon coordinates to do the overlay, and hence
; tfDoNDCOverlay should not be set to True. (The default is False.)
;
begin
  print("wrf_overlay: Warning: This procedure is obsolete. Consider using")
  print("                      wrf_overlays instead.")

  opts = opt_arg      ; Make a copy of the resource list.
  opts = True

  if(opts.and.isatt(opts,"NoTitles").and.opts@NoTitles) then
    no_titles = True
  else
    no_titles = False
  end if

  panel_plot = get_res_value(opts,"PanelPlot",False)    ; Are we paneling?
  latlon_overlay = get_res_value(opts,"LatLonOverlay",False)  ; lat/lon overlay?

  nplots = dimsizes(plots)
    
  base = plots(0)
  if(.not.no_titles) then
    getvalues plots(0)
      "tiMainString" : SubTitle
    end getvalues
    SubTitles = SubTitle
    if(.not.latlon_overlay) then
      setvalues plots(0)
        "tfDoNDCOverlay" : True
        "tiMainOn"       : False
      end setvalues
    else
      setvalues plots(0)
        "tiMainOn"       : False
      end setvalues
    end if
  else
    setvalues plots(0)
      "tfDoNDCOverlay" : True
    end setvalues
  end if

  if (nplots.eq.1) then 
    blank = create "BlankPlot" logLinPlotClass wks
      ;"cnConstFLabelOn" : False
    end create
    overlay(base,blank)
  end if

  do i=1,nplots-1
    if(.not.ismissing(plots(i))) then
      if(.not.no_titles) then
        getvalues plots(i)
          "tiMainString" : SubTitle
        end getvalues
        SubTitles = SubTitles + "~C~" + SubTitle
        if(.not.latlon_overlay) then
          setvalues plots(i)
            "tfDoNDCOverlay" : True
            "tiMainOn"       : False
          end setvalues
        else
          setvalues plots(i)
            "tiMainOn"       : False
          end setvalues
        end if
      else
        if(.not.latlon_overlay) then
          setvalues plots(i)
            "tfDoNDCOverlay" : True
          end setvalues
        end if
      end if
      overlay(base,plots(i))
    else
      print("wrf_overlay: Warning: overlay plot #" + i + " is not valid.")
    end if
  end do

  if(.not.no_titles) then
    font_height = get_res_value_keep(opts,"FontHeightF",0.01)

    txt = create "map_titles" textItemClass wks
      "txString"      : SubTitles
      "txFontHeightF" : font_height
    end create
    anno = NhlAddAnnotation(base,txt)
    setvalues anno
      "amZone"           : 3
      "amJust"           : "BottomLeft"
      "amSide"           : "Top"
      "amParallelPosF"   : 0.005
      "amOrthogonalPosF" : 0.03
      "amResizeNotify"   : False
    end setvalues
    base@map_titles = anno
  end if

;
; gsnDraw and gsnFrame should default to True if not a panel plot.
; gsnFrame will default to False if opt_arg is False.
;
  if(panel_plot.or..not.opt_arg) then
    call_frame = False
  else
    call_frame = True
  end if
  if(panel_plot) then
    call_draw = False
  else
    call_draw = True
  end if
  opts@gsnDraw     = get_res_value_keep(opts,"gsnDraw",     call_draw)
  opts@gsnFrame    = get_res_value_keep(opts,"gsnFrame",    call_frame)
  opts@gsnMaximize = get_res_value_keep(opts,"gsnMaximize", True)

  draw_and_frame(wks,base,opts@gsnDraw,opts@gsnFrame,False, \
                 opts@gsnMaximize)

  if(.not.no_titles.and..not.panel_plot) then
    NhlRemoveAnnotation(base,base@map_titles)
    delete(base@map_titles)
  end if

  if(.not.panel_plot) then
    if ( nplots .ge. 2 ) then
      do i=1,nplots-1
        if(.not.ismissing(plots(i))) then
          NhlRemoveOverlay(base,plots(i),False)
        else
          print("wrf_remove_overlay: Warning: overlay plot #" + i + " is not valid.")
          print("                        Nothing to remove.")
        end if
      end do
    end if
  end if
end

;--------------------------------------------------------------------------------

undef("add_white_space")
function add_white_space(str:string,maxlen:integer)

begin
  cstr = stringtochar(str)
  len  = dimsizes(cstr)-1
  ws   = ""
  if(len.lt.maxlen) then
    do i=1,maxlen-len
      ws = ws + " "
    end do
  end if
  return(ws)

end

;--------------------------------------------------------------------------------

undef("print_opts")
procedure print_opts(opts_name,opts,debug)

begin
  if(.not.debug) then
    return
  end if
  varatts = getvaratts(opts)
;
; Sort attributes alphabetically/
;
  sqsort(varatts)
;
; Get length of longest attribute name.
;
  cvaratts = stringtochar(varatts)
  cmaxlen  = dimsizes(cvaratts(0,:))-1

  print("------------------------------------------------------------")
  print(opts_name + "...")    ; Print name of option variable.
;
; Loop through each attribute in the list. If its value is an array,
; then print out the array with '(/' and '/)' around it.
;
; If the value is a string, then put ticks (') around each string.
;
  do i=0,dimsizes(varatts)-1
    x = opts@$varatts(i)$
;
; Use add_white_space to make sure all the equal signs line up.
;
    tmp_str = "    " + varatts(i) + \
              add_white_space(varatts(i),cmaxlen) + " = "
;
; Check if attribute is an array.
;    
    if(dimsizes(x).gt.1) then
      tmp_str = tmp_str + "(/"
      do j=0,dimsizes(x)-1
        if(typeof(x(j)).eq."string") then
          tmp_str = tmp_str + "'" + x(j) + "'"
        else
          tmp_str = tmp_str + x(j)
        end if
        if(j.lt.dimsizes(x)-1) then
          tmp_str = tmp_str + ","
        else
          tmp_str = tmp_str + "/)"
        end if
      end do
    else      if(typeof(x).eq."string") then
        tmp_str = tmp_str + "'" + x + "'"
      else
        tmp_str = tmp_str + x
      end if
    end if
    print("" + tmp_str)
    delete(x)
  end do

end

;--------------------------------------------------------------------------------

undef("print_header")
procedure print_header(icount:integer,debug)
begin
  icount = icount + 1
  if(.not.debug) then
    return
  end if
  print("END plot #" + icount)
  print("------------------------------------------------------------")

end

;--------------------------------------------------------------------------------
undef("wrf_cloud_fraction")
function wrf_cloud_fraction(file_handle,time:integer)
begin

   print("in wrf cloud fraction")
   if(typeof(file_handle).eq."file") then
      ISFILE = True
      nc_file = file_handle
   else if(typeof(file_handle).eq."list") then
      ISFILE = False
      nc_file = file_handle[0]
   else
      print("wrf_user_getvar: error: the first argument must be a file or a list of files opened with addfile or addfiles")
      return
   end if
   end if

   time_in = time
   if ( time .eq. -1 ) then
      if(ISFILE) then
         P    = nc_file->P
         PB   = nc_file->PB
         rh   = wrf_user_getvar(nc_file,"rh",-1)
      else
         P    = file_handle[:]->P
         PB   = file_handle[:]->PB
         print("here now before rh")
         rh   = wrf_user_getvar(file_handle,"rh",-1)
      end if
   else
      if(ISFILE) then
         P    = nc_file->P(time_in,:,:,:)
         PB   = nc_file->PB(time_in,:,:,:)
         rh   = wrf_user_getvar(file_handle,"rh",time_in)
      else
         P    = file_handle[:]->P(time_in,:,:,:)
         PB   = file_handle[:]->PB(time_in,:,:,:)
         rh   = wrf_user_getvar(file_handle,"rh",time_in)
      end if
  end if


  pres = (P + PB)
  thedims = dimsizes(pres)


  if(time .eq. -1) then
    nt = thedims(0)
    nz = thedims(1)
    ns = thedims(2)
    ew = thedims(3)
    out_array = new( (/nt,3,ns,ew/),float)
  else
    nt = 1
    nz = thedims(0)
    ns = thedims(1)
    ew = thedims(2)
    out_array = new( (/3,ns,ew/),float)
  end if

   lowc  = new((/ns,ew/),float)
   midc  = new((/ns,ew/),float)
   highc = new((/ns,ew/),float)
   lowc  = 0.0
   midc  = 0.0
   highc = 0.0

   if(nt .gt. 1 ) then
      do it = 0, nt-1
         print("Working on Time "+it)
         fortran_clouds::cloud_frac(pres(it,:,:,:),rh(it,:,:,:),lowc(:,:),midc(:,:),highc(:,:),nz,ns,ew)
         out_array(it,0,:,:) = lowc(:,:)
         out_array(it,1,:,:) = midc(:,:)
         out_array(it,2,:,:) = highc(:,:)
      end do
   else
      fortran_clouds::cloud_frac(pres(:,:,:),rh(:,:,:),lowc(:,:),midc(:,:),highc(:,:),nz,ns,ew)
      out_array(0,:,:) = lowc(:,:)
      out_array(1,:,:) = midc(:,:)
      out_array(2,:,:) = highc(:,:)      
   end if

   out_array@units = "%"
   out_array@description = "Low, Med High Clouds"
   print("Done with clouds")
   return(out_array)
end