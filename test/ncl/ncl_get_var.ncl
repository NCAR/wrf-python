    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
    
    ;system("printenv")
    
    if (.not. isvar("dir")) then
        dir = "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/moving_nest"
        ;in_file = "/Users/ladwig/Documents/wrf_files/wrfout_d01_2010-06-13_21:00:00.nc"
        ;in_file = "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/moving_nest/wrfout_d02_2005-08-28_00:00:00.nc"
    end if 
    
    if (.not. isvar("pattern")) then
        pattern = "wrfout_d02_*"
    end if
    
    if (.not. isvar("out_file")) then
        out_file = "/tmp/wrftest.nc"
    end if
    
    pat = dir + "/" + pattern
    cmd = "ls " + pat
    fils = systemfunc (cmd) ; file paths
    input_file = addfiles(fils,"r")
    
    system("/bin/rm -f " + out_file) ; remove if exists
    fout = addfile(out_file, "c")
    
    time = -1
    
    wrf_vars = [/"avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", \
                "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", \
                "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", \
                "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", \
                "wa", "uvmet10", "uvmet", "z", "cfrac", "height_agl", \
                "wspd_wdir", "wspd_wdir10", "uvmet_wspd_wdir", \
                "uvmet10_wspd_wdir" /]
    
    unique_dimname_list = NewList("fifo")
    unique_dimsize_list = NewList("fifo")
    full_vardimname_list = NewList("fifo") ; Workaround for issue where NCL
                                           ; is dropping the dim names from
                                           ; the array stored in a list
    vardata_list = NewList("fifo")
    
    ; NCL lists need unique variable names to be inserted, so using these
    ; variables to create unique named attributes
    vardata = True
    vardimnamedata = True
    
    ; Note:  The list type seems to only work correctly when inserting
    ; variables with unique names.  This is the reason for all of the 
    ; name attribute stuff below.     
    do i = 0, ListCount(wrf_vars) - 1
       print("working on: " + wrf_vars[i])
       v := wrf_user_getvar(input_file, wrf_vars[i], time)
       fout->$wrf_vars[i]$ = v
    end do 
    
    
    ; Do the interpolation routines manually
    ;;;;;;;;;;;;;;;;;;; 3D vertical cross section
    time = -1
    
    z := wrf_user_getvar(input_file, "z", time)        ; grid point height
    p := wrf_user_getvar(input_file, "pressure", time) ; total pressure
    
    
    dimsz = dimsizes(z)
    pivot = (/ dimsz(3)/2, dimsz(2)/2 /)    ; pivot point is center of domain
    
    ht_cross = wrf_user_intrp3d(z,p,"v",pivot,90.0,False)
    
    p_cross = wrf_user_intrp3d(p,z,"v",pivot,90.0,False)
    p_cross!0 = "Vertical_p"
    p_cross!1 = "Horizontal_p"
    
    fout->ht_cross = ht_cross
    fout->p_cross = p_cross
 

    time = 0
    
    ; For the new cross section routine
    xopt = True
    xopt@use_pivot = True
    xopt@angle = 90.0
    xopt@file_handle = input_file
    xopt@timeidx = time
    xopt@linecoords = True
    
    ht_vertcross1 = wrf_user_vert_cross(z, p, pivot, xopt)
    
    fout->ht_vertcross1 = ht_vertcross1
    
    ; For the new cross section routine
    xopt := True
    xopt@use_pivot = True
    xopt@angle = 90.0
    xopt@levels = (/1000., 850., 700., 500., 250./)
    xopt@file_handle = input_file
    xopt@timeidx = time
    xopt@linecoords = True
    
    ht_vertcross2 = wrf_user_vert_cross(z, p, pivot, xopt)
    ht_vertcross2!1 = "vertical2"
    ht_vertcross2!2 = "cross_line_idx2"
    
    fout->ht_vertcross2 = ht_vertcross2
    
    ; Can only use a single time for lat/lon version at this time
    
    vertdims = dimsizes(ht_vertcross2)
    htdims = dimsizes(z)
    
    lats = wrf_user_getvar(input_file, "lat", 0)
    lons = wrf_user_getvar(input_file, "lon", 0) 
        
    start_lat = min(lats) + .25d*(max(lats) - min(lats))
    end_lat = min(lats) + .65d*(max(lats) - min(lats))
                
    start_lon = min(lons) + .25d*(max(lons) - min(lons))
    end_lon = min(lons) + .65d*(max(lons) - min(lons))
        
    start_end = (/ start_lon, start_lat, end_lon, end_lat /)
    
    ; For the new cross section routine
    xopt := True
    xopt@use_pivot = False
    xopt@latlon = True
    xopt@file_handle = input_file
    xopt@timeidx = 0
    xopt@linecoords = True
    xopt@autolevels = 1000
        
    ht_vertcross3 = wrf_user_vert_cross(z, p, start_end, xopt)
    
    ht_vertcross3!0 = "Time"
    ht_vertcross3!1 = "vertical3"
    ht_vertcross3!2 = "cross_line_idx3"
    
    fout->ht_vertcross3 = ht_vertcross3
    
    ; Test the moving nest with lat/lon over time
    
    times = wrf_user_getvar(input_file, "times", -1)
    ntimes = dimsizes(times)
    
    do i=0,ntimes-1
      xopt@timeidx = i
      name = sprinti("ht_vertcross_t%i", i)
      p_var := p(i,:,:,:)
      z_var := z(i,:,:,:)

      ht_vertcross := wrf_user_vert_cross(z_var, p_var, start_end, xopt)
      
      dim0name = sprinti("vertical_t%i",i)
      dim1name = sprinti("cross_line_idx_t%i",i)
      ht_vertcross!0 = dim0name
      ht_vertcross!1 = dim1name
      
      fout->$name$ = ht_vertcross
    end do
    
    ;;;;;;;;;;;;;;;;;;;;;;;; 3D horizontal interpolation
    
    time = -1
    
    z := wrf_user_getvar(input_file, "z", time)        ; grid point height
    p := wrf_user_getvar(input_file, "pressure", time) ; total pressure
    
    ; First, check backwards compat
    plev = 500.   ; 500 MB
    hlev = 5000; ; 5000 m
    
    z_500 = wrf_user_intrp3d(z,p,"h",plev,0.,False)
    p_5000 = wrf_user_intrp3d(p,z,"h",hlev,0.,False)
    
    fout->z_500 = z_500
    fout->p_5000 = p_5000
    
    plev := (/1000., 850., 500., 250./) 
    hlev := (/500., 2500., 5000., 10000. /)
    z_multi = wrf_user_intrp3d(z,p,"h",plev,0.,False)
    p_multi = wrf_user_intrp3d(p,z,"h",hlev,0.,False)
    
    fout->z_multi = z_multi
    fout->p_multi = p_multi
    
    ; Now check the new routine
    
    plev := 500.   ; 500 MB
    hlev := 5000   ; 5000 m
    
    z2_500 = wrf_user_interp_level(z,p,plev,False)
    p2_5000 = wrf_user_interp_level(p,z,hlev,False)
    
    fout->z2_500 = z2_500
    fout->p2_5000 = p2_5000
    
    
    plev := (/1000., 850., 500., 250./) 
    hlev := (/500., 2500., 5000., 10000. /)
    z2_multi = wrf_user_interp_level(z,p,plev,False)
    p2_multi = wrf_user_interp_level(p,z,hlev,False)
    
    fout->z2_multi = z2_multi
    fout->p2_multi = p2_multi
    
    pblh = wrf_user_getvar(input_file, "PBLH", time)
    opts := False
    opts@inc2dlevs = True
    p_lev2d = wrf_user_interp_level(p, z, pblh, opts)
    
    fout->p_lev2d = p_lev2d
    
    
    ;;;;;;;;;;;;;;;;;;;;;;;; 2D interpolation along line
    
    time = -1
    
    t2 = wrf_user_getvar(input_file, "T2", time)
    dimst2 = dimsizes(t2)
    pivot = (/ dimst2(2)/2, dimst2(1)/2 /)
    
    t2_line = wrf_user_intrp2d(t2, pivot, 90.0, False)
    
    fout->t2_line = t2_line
    
    ; For the new interplevel routine
    xopt := True
    xopt@use_pivot = True
    xopt@angle = 90.0
    xopt@latlon = False
    xopt@file_handle = input_file
    xopt@timeidx = 0
    xopt@linecoords = True
    
    t2_line2 = wrf_user_interp_line(t2, pivot, xopt)
    
    fout->t2_line2 = t2_line2
    
    lats = wrf_user_getvar(input_file, "lat", 0)
    lons = wrf_user_getvar(input_file, "lon", 0) 
        
    start_lat = min(lats) + .25d*(max(lats) - min(lats))
    end_lat = min(lats) + .65d*(max(lats) - min(lats))
                
    start_lon = min(lons) + .25d*(max(lons) - min(lons))
    end_lon = min(lons) + .65d*(max(lons) - min(lons))
        
    start_end = (/ start_lon, start_lat, end_lon, end_lat /)
    
    ; For the new line routine
    xopt := True
    xopt@use_pivot = False
    xopt@latlon = True
    xopt@file_handle = input_file
    xopt@timeidx = 0
    xopt@linecoords = True
    
    t2_line3 = wrf_user_interp_line(t2, start_end, xopt)
    t2_line3!1 = "line_idx_t2_line3"
    
    fout->t2_line3 = t2_line3
    
    times = wrf_user_getvar(input_file, "times", -1)
    ntimes = dimsizes(times)
    
    do i=0,ntimes-1
      xopt@timeidx = i
      name = sprinti("t2_line_t%i", i)
      dim0name = sprinti("lineidx_t%i",i)
      var := t2(i,:,:)
      t2_line := wrf_user_interp_line(var, start_end, xopt)
      t2_line!0 = dim0name
      fout->$name$ = t2_line
    end do
    
    ; Make sure the 1 time case still works
    t2 := wrf_user_getvar(input_file, "T2", 0)
    
    ; For the new line routine
    xopt := True
    xopt@use_pivot = False
    xopt@latlon = True
    xopt@file_handle = input_file
    xopt@timeidx = 0
    xopt@linecoords = True
    
    t2_line4 = wrf_user_interp_line(t2, start_end, xopt)
    t2_line4!0 = "t2_line4_idx"
    
    fout->t2_line4 = t2_line4
    
    
    ;;;;;;;;;;;;;;;;;;;;;;; 3D interpolation to new vertical coordinates
    time = -1
    
    ; interp t to theta
    fld1 = wrf_user_getvar(input_file, "tk", time)
    vert_coord       = "theta"
    interp_levels    = (/200,300,500,1000/)
 
    opts             = True
    opts@extrapolate = True 
    opts@field_type  = "T"
    opts@logP        = True 
    
    fld1_intrp = wrf_user_vert_interp(input_file,fld1,vert_coord,interp_levels,opts)
    fld1_intrp!1 = "interp_levels1"
    
    fout->fld_tk_theta = fld1_intrp
    
    ; interp t to theta-e
    fld2 = fld1
    vert_coord := "theta-e"
    fld2_intrp = wrf_user_vert_interp(input_file,fld2,vert_coord,interp_levels,opts)
    fld2_intrp!1 = "interp_levels2"
    
    fout->fld_tk_theta_e = fld2_intrp
    
    
    ; interp t to pressure
    fld3 = fld1
    vert_coord := "pressure"
    interp_levels    := (/850,500/)
    fld3_intrp = wrf_user_vert_interp(input_file,fld3,vert_coord,interp_levels,opts)
    fld3_intrp!1 = "interp_levels3"
    
    fout->fld_tk_pres = fld3_intrp
    
    
    ; interp t to ght_msl
    fld4 = fld1
    vert_coord := "ght_msl"
    interp_levels    := (/1,2/)
    fld4_intrp = wrf_user_vert_interp(input_file,fld4,vert_coord,interp_levels,opts)
    fld4_intrp!1 = "interp_levels4"
    
    fout->fld_tk_ght_msl = fld4_intrp
    
    
    ; interp t to ght_agl
    fld5 = fld1
    vert_coord := "ght_agl"
    interp_levels    := (/1,2/)
    fld5_intrp = wrf_user_vert_interp(input_file,fld1,vert_coord,interp_levels,opts)
    fld5_intrp!1 = "interp_levels5"
    
    fout->fld_tk_ght_agl = fld5_intrp
    
    ; interp ht to pres
    fld6 = wrf_user_getvar(input_file, "height", time)
    vert_coord := "pressure"
    opts@field_type  = "ght"
    interp_levels    := (/500,50/)
    fld6_intrp = wrf_user_vert_interp(input_file,fld6,vert_coord,interp_levels,opts)
    fld6_intrp!1 = "interp_levels6"
    
    fout->fld_ht_pres = fld6_intrp
    
    
    ; interp pres to theta
    fld7 = wrf_user_getvar(input_file, "pressure", time)
    vert_coord := "theta"
    opts@field_type  = "pressure"
    interp_levels    := (/200,300,500,1000/)
    fld7_intrp = wrf_user_vert_interp(input_file,fld7,vert_coord,interp_levels,opts)
    fld7_intrp!1 = "interp_levels7"
    
    fout->fld_pres_theta = (/fld7_intrp/)
    
    
    ; interp theta-e to pressure
    fld8 = wrf_user_getvar(input_file, "eth", time)
    vert_coord := "pressure"
    opts@field_type  = "T"
    interp_levels    := (/850,500,5/)
    fld8_intrp = wrf_user_vert_interp(input_file,fld8,vert_coord,interp_levels,opts)
    fld8_intrp!1 = "interp_levels8"
    
    fout->fld_thetae_pres = fld8_intrp
    
    
    ;;;;;;;;;;;;;;;;;;; lat/lon to x/y and x/y to lat/lon routines
    
    lats := (/22.0, 25.0, 27.0 /)
    lons := (/-90.0, -87.5, -83.75 /)
    x_s = (/10, 50, 90 /)
    y_s = (/10, 50, 90 /)
    
    opt = True
    opt@useTime = -1
    opt@returnInt = False
    
    xy1 = wrf_user_ll_to_xy(input_file, lons, lats, opt) 
    ll1 = wrf_user_xy_to_ll(input_file, x_s, y_s, opt)
    
    fout->xy1 = xy1
    fout->ll1 = ll1
    
    opt = True
    opt@useTime = 8
    opt@returnInt = True
    
    xy2 = wrf_user_ll_to_xy(input_file, lons, lats, opt) 
    ll2 = wrf_user_xy_to_ll(input_file, x_s, y_s, opt)
    
    fout->xy2 = xy2
    fout->ll2 = ll2
    
    delete(fout)







