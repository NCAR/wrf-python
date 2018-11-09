    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
    
    ;system("printenv")
    
    if (.not. isvar("in_file")) then
        in_file = "/Users/ladwig/Documents/wrf_files/wrfout_d01_2010-06-13_21:00:00.nc"
        ;in_file = "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-28_00:00:00.nc"
    end if 
    
    if (.not. isvar("out_file")) then
        out_file = "/tmp/wrftest.nc"
    end if
    input_file = addfile(in_file,"r")
    
    system("/bin/rm -f " + out_file) ; remove if exists
    fout = addfile(out_file, "c")
    
    time = 0
    
    wrf_vars = [/"avo", "eth", "cape_2d", "cape_3d", "ctt", "dbz", "mdbz", \
                "geopt", "helicity", "lat", "lon", "omg", "p", "pressure", \
                "pvo", "pw", "rh2", "rh", "slp", "ter", "td2", "td", "tc", \
                "theta", "tk", "tv", "twb", "updraft_helicity", "ua", "va", \
                "wa", "uvmet10", "uvmet", "z", "cfrac", "height_agl" /]
    
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
    ; 3D vertical cross section
    z := wrf_user_getvar(input_file, "z", 0)        ; grid point height
    p := wrf_user_getvar(input_file, "pressure", 0) ; total pressure
    
    dimsz = dimsizes(z)
    pivot = (/ dimsz(2)/2, dimsz(1)/2 /)    ; pivot point is center of domain
    
    ht_cross = wrf_user_intrp3d(z,p,"v",pivot,90.0,False)
    
    p_cross = wrf_user_intrp3d(p,z,"v",pivot,90.0,False)
    p_cross!0 = "Vertical_p"
    p_cross!1 = "Horizontal_p"
    
    fout->ht_cross = ht_cross
    fout->p_cross = p_cross
;    
    ; For the new cross section routine
    xopt = True
    xopt@use_pivot = True
    xopt@angle = 90.0
    xopt@file_handle = input_file
    xopt@timeidx = 0
    xopt@linecoords = True
    
    ht_vertcross1 = wrf_user_vertcross(z, p, pivot, xopt)
    
    fout->ht_vertcross1 = ht_vertcross1
    
    ; For the new cross section routine
    xopt := True
    xopt@use_pivot = True
    xopt@angle = 90.0
    xopt@levels = (/1000., 850., 700., 500., 250./)
    xopt@file_handle = input_file
    xopt@timeidx = 0
    xopt@linecoords = True
    
    ht_vertcross2 = wrf_user_vertcross(z, p, pivot, xopt)
    ht_vertcross2!0 = "vertical2"
    ht_vertcross2!1 = "cross_line_idx2"
    
    
    fout->ht_vertcross2 = ht_vertcross2
    
    lats = wrf_user_getvar(input_file, "lat", 0)
    lons = wrf_user_getvar(input_file, "lon", 0) 
    
    start_lat = min(lats) + .25d*(max(lats) - min(lats))
    end_lat = min(lats) + .75d*(max(lats) - min(lats))
            
    start_lon = min(lons) + .25d*(max(lons) - min(lons))
    end_lon = min(lons) + .75d*(max(lons) - min(lons))
    
    start_end = (/ start_lon, start_lat, end_lon, end_lat /)
    
    ; For the new cross section routine
    xopt := True
    xopt@use_pivot = False
    xopt@latlon = True
    xopt@file_handle = input_file
    xopt@timeidx = 0
    xopt@linecoords = True
    xopt@autolevels = 1000
    
    ht_vertcross3 = wrf_user_vertcross(z, p, start_end, xopt)
    ht_vertcross3!0 = "vertical3"
    ht_vertcross3!1 = "cross_line_idx3"
    
    fout->ht_vertcross3 = ht_vertcross3
    
    ; 3D horizontal interpolation
    
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
    
    z2_500 = wrf_user_interplevel(z,p,plev,False)
    p2_5000 = wrf_user_interplevel(p,z,hlev,False)
    
    fout->z2_500 = z2_500
    fout->p2_5000 = p2_5000
    
    
    plev := (/1000., 850., 500., 250./) 
    hlev := (/500., 2500., 5000., 10000. /)
    z2_multi = wrf_user_interplevel(z,p,plev,False)
    p2_multi = wrf_user_interplevel(p,z,hlev,False)
    
    fout->z2_multi = z2_multi
    fout->p2_multi = p2_multi
    
    pblh = wrf_user_getvar(input_file, "PBLH", 0)
    opts := False
    opts@inc2dlevs = True
    p_lev2d = wrf_user_interplevel(p, z, pblh, opts)
    
    fout->p_lev2d = p_lev2d
    
    
    ; 2D interpolation along line
    t2 = wrf_user_getvar(input_file, "T2", 0)
    dimst2 = dimsizes(t2)
    pivot = (/ dimst2(1)/2, dimst2(0)/2 /)
    
    t2_line = wrf_user_intrp2d(t2, pivot, 90.0, False)
    
    fout->t2_line = (/t2_line/)
    
    
    ; 3D interpolation to new vertical coordinates
    ; interp t to theta
    fld1 = wrf_user_getvar(input_file, "tk", -1)
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
    fld6 = wrf_user_getvar(input_file, "height", -1)
    vert_coord := "pressure"
    opts@field_type  = "ght"
    interp_levels    := (/500,50/)
    fld6_intrp = wrf_user_vert_interp(input_file,fld6,vert_coord,interp_levels,opts)
    fld6_intrp!1 = "interp_levels6"
    
    fout->fld_ht_pres = fld6_intrp
    
    
    ; interp pres to theta
    fld7 = wrf_user_getvar(input_file, "pressure", -1)
    vert_coord := "theta"
    opts@field_type  = "pressure"
    interp_levels    := (/200,300,500,1000/)
    fld7_intrp = wrf_user_vert_interp(input_file,fld7,vert_coord,interp_levels,opts)
    fld7_intrp!1 = "interp_levels7"
    
    fout->fld_pres_theta = (/fld7_intrp/)
    
    
    ; interp theta-e to pressure
    fld8 = wrf_user_getvar(input_file, "eth", -1)
    vert_coord := "pressure"
    opts@field_type  = "T"
    interp_levels    := (/850,500,5/)
    fld8_intrp = wrf_user_vert_interp(input_file,fld8,vert_coord,interp_levels,opts)
    fld8_intrp!1 = "interp_levels8"
    
    fout->fld_thetae_pres = fld8_intrp
    
    
    ; lat/lon to x/y and x/y to lat/lon routines
    lats := (/-55, -60, -65 /)
    lons := (/25, 30, 35 /)
    i_s = (/10, 100, 150 /)
    j_s = (/10, 100, 150 /)
    ij = wrf_user_ll_to_ij(input_file, lons, lats, True) 
    ll = wrf_user_ij_to_ll(input_file, i_s, j_s, True)
    
    fout->ij = ij
    fout->ll = ll
    
    delete(fout)







