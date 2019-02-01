  a = addfile("wrfout_d03_2017-04-03_06:00:00_ctrl","r")
   
  time = 0
  refl_10cm = wrf_user_getvar(a,"REFL_10CM",time)
  z = wrf_user_getvar(a, "z", time)
  lat = wrf_user_getvar(a, "lat", time)
  lon = wrf_user_getvar(a, "lon", time)
  
  ; convert the lat/lon to x,y
  start_lat = 20.9
  start_lon = 92.5
  end_lat = 29.2
  end_lon = 92.5
  
  opt = True
  
  start_ij = wrf_user_ll_to_ij(a, start_lon, start_lat, opt)
  start_ij = start_ij - 1
  
  end_ij = wrf_user_ll_to_ij(a, end_lon, end_lat, opt)
  end_ij = end_ij - 1
  
  start_end =  (/start_ij(0), start_ij(1), end_ij(0), end_ij(1)/)
  
  lat_line = wrf_user_intrp2d(lat,start_end,0.0,True)
  nlat = dimsizes(lat_line)
  
  lon_line = wrf_user_intrp2d(lon,start_end,0.0,True)

  refl_cross = wrf_user_intrp3d(refl_10cm,z,"v",start_end,0.,True)
  
  ; Need to make a vertical coordinate by using the same code as the 
  ; cross section
  
  ; Currently, the vertical coordinate is not set, so let's do it 
  ; manually. This will be fixed in the next version of NCL.
  ; If you want to set these levels yourself, you'll need to copy the 
  ; code I sent before and manually set the levels in the cross section 
  ; routine, then do it again here.
  
   z_max = max(z)
   z_min = 0.
   dz = 0.01 * z_max
   nlevels = tointeger( z_max/dz )
   z_var2d = new( (/nlevels/), typeof(z))
   z_var2d(0) = z_min
  
   do i=1, nlevels-1
      z_var2d(i) = z_var2d(0)+i*dz
   end do
   
   refl_cross&Vertical = z_var2d
   
   wks = gsn_open_wks("png","cross")
   cmap     := read_colormap_file("BlAqGrYeOrReVi200")
   cmap(0,:) = (/0,0,0,0/)    ; make first color fully transparent

   resx = True
   resx@gsnMaximize          = True
   resx@lbLabelAutoStride    = True            ; default v6.1.0

   resx@cnFillOn             = True            ; turn on color fill
   resx@cnLinesOn            = False           ; turn lines on/off         ; True is default
   resx@cnLineLabelsOn       = False          ; turn line labels on/off   ; True is default
   resx@cnFillPalette        = cmap
   nLabels                   = 8              ; arbitrary
   resx@tmXBLabels           = new(nLabels,"string")
   resx@tmXBMode             = "Explicit"

   resx@tmXBValues           := toint(fspan(0,nlat-1,nLabels))
   do i=0,nLabels-1
      x = lon_line(i)
      y = lat_line(i)
      resx@tmXBLabels(i) = sprintf("%5.1f", y)+"~C~"+sprintf("%5.1f", x)
   end do

   resx@tiMainString = "Full South-North Grid Line X-Section"


   plot1  = gsn_csm_contour(wks, refl_cross, resx )
   