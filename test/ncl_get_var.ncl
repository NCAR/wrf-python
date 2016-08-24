    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
    
    ;system("printenv")
    
    if (.not. isvar("in_file")) then
        in_file = "/Users/ladwig/Documents/wrf_files/wrfout_d02_2010-06-13_21:00:00.nc"
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
                "wa", "uvmet10", "uvmet", "z", "cfrac"/]
    
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
       
       ;if (wrf_vars[i] .eq. "avo") then
       ;    print(v)
       ;end if
       
       ; pw is written in pure NCL and does not contain dimension names
       ; so manually creating the dimension names here
       if (wrf_vars[i] .eq. "pw") then
           dim_names := (/"south_north", "west_east"/)
           dim_sizes := dimsizes(v)
       else
           dim_names := getvardims(v)
           dim_sizes := dimsizes(v)
       end if 
       
       vardata@$wrf_vars[i]$ := v
       vardimnamedata@$wrf_vars[i]$ := dim_names
       ListAppend(vardata_list,vardata@$wrf_vars[i]$) 
       ListAppend(full_vardimname_list, vardimnamedata@$wrf_vars[i]$)
       ;print(vardata_list)
       
       dimname=True
       dimsize=True
       
       ; Determine the unique dimensions names, which will be used when 
       ; creating the output NetCDF file
       do j=0, dimsizes(dim_sizes)-1
           ;print(dim_names)
           ;print(dim_names(j))
           
           name_id = sprintf("dimname_%i",i*j)
           size_id = sprintf("dimsize_%i",i*j)
           
           dimname@$name_id$ = dim_names(j) 
           dimsize@$size_id$ = dim_sizes(j)
           
           has_name = False
           do k=0, ListCount(unique_dimname_list)-1
               if ((/unique_dimname_list[k]/) .eq. (/dimname@$name_id$/)) then
                   has_name = True
               end if
           end do

           if (.not. has_name) then 
               ;print("inserting: " + dimname@$name_id$)
               ListAppend(unique_dimname_list, dimname@$name_id$)
               ListAppend(unique_dimsize_list, dimsize@$size_id$)
           end if
           
       end do
    end do 
        
    setfileoption(fout,"DefineMode",True)
    
    ; Set global attributes
    f_att               = True            ; assign file attributes
    f_att@title         = "NCL generated netCDF file"
    f_att@Conventions   = "None"
    fileattdef(fout, f_att)      ; copy file attributes
    
    ; Set up the NetCDF dimensions
    d_names = new(ListCount(unique_dimname_list), string)
    d_sizes = new(ListCount(unique_dimname_list), integer)
    d_unlim  = new(ListCount(unique_dimname_list), logical)
    
    ; Note: Need to do this copy since NCL can't coerce the list data to 
    ; array data
    do i=0, ListCount(unique_dimname_list) - 1
        d_names(i) = unique_dimname_list[i]
        d_sizes(i) = unique_dimsize_list[i]
        d_unlim(i) = False
    end do
    
    filedimdef(fout, d_names, d_sizes, d_unlim)
    
    ; Save the variables to the NetCDF file
    do i=0, ListCount(vardata_list)-1       
       d := vardata_list[i]
       filevardef(fout, wrf_vars[i], typeof(d), full_vardimname_list[i])
       filevarattdef(fout,wrf_vars[i], d)
       fout->$wrf_vars[i]$ = (/d/)
    end do
    
    ; Do the interpolation routines manually
    ; 3D vertical cross section
    z  = wrf_user_getvar(input_file, "z", 0)        ; grid point height
    p  = wrf_user_getvar(input_file, "pressure", 0) ; total pressure
    
    dimsz = dimsizes(z)
    pivot = (/ dimsz(2)/2, dimsz(1)/2 /)    ; pivot point is center of domain
    
    ht_cross = wrf_user_intrp3d(z,p,"v",pivot,90.0,False)
    p_cross = wrf_user_intrp3d(p,z,"v",pivot,90.0,False)
    
    ht_cross_dims = dimsizes(ht_cross)
    p_cross_dims = dimsizes(p_cross)
    
    ; 3D horizontal interpolation
    plev = 500.   ; 500 MB
    z_500 = wrf_user_intrp3d(z,p,"h",plev,0.,False)
    z_500_dims = dimsizes(z_500)
    
    
    ; 2D interpolation along line
    t2 = wrf_user_getvar(input_file, "T2", 0)
    dimst2 = dimsizes(t2)
    pivot = (/ dimst2(1)/2, dimst2(0)/2 /)
    
    t2_line = wrf_user_intrp2d(t2, pivot, 90.0, False)
    t2_line_dims = dimsizes(t2_line) 
    
    filedimdef(fout, (/"ht_cross_vert", "ht_cross_horiz", "p_cross_vert", "p_cross_horiz"/), \
                     (/ht_cross_dims(0), ht_cross_dims(1), p_cross_dims(0), p_cross_dims(1)/), \ 
                     (/False,False,False,False/))
    filedimdef(fout, (/"z500_vert", "z500_horiz"/), \
                     (/z_500_dims(0), z_500_dims(1) /), \
                     (/False, False/))                 
    filedimdef(fout, (/"t2_line_horiz"/), \
                     (/t2_line_dims(0) /), \
                     (/False/))
    
    filevardef(fout, "ht_cross", typeof(ht_cross), (/"ht_cross_vert", "ht_cross_horiz"/))
    filevardef(fout, "p_cross", typeof(p_cross), (/"p_cross_vert", "p_cross_horiz"/))
    filevardef(fout, "z_500", typeof(z_500), (/"z500_vert", "z500_horiz"/))
    filevardef(fout, "t2_line", typeof(t2_line), (/"t2_line_horiz"/))
    
    filevarattdef(fout,"ht_cross", ht_cross)
    filevarattdef(fout,"p_cross", p_cross)
    filevarattdef(fout, "z_500", z_500)
    filevarattdef(fout, "t2_line", t2_line)
    fout->ht_cross = (/ht_cross/)
    fout->p_cross = (/p_cross/)
    fout->z_500 = (/z_500/)
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
    
    fld1_dims = dimsizes(fld1_intrp)
    
    filedimdef(fout, (/"fld1_time", "fld1_levs", "fld1_sn", "fld1_we"/), \
                     (/fld1_dims(0), fld1_dims(1), fld1_dims(2), fld1_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_tk_theta", typeof(fld1_intrp), (/"fld1_time", "fld1_levs", "fld1_sn", "fld1_we"/))
    filevarattdef(fout,"fld_tk_theta", fld1_intrp)
    fout->fld_tk_theta = (/fld1_intrp/)
    
    ; interp t to theta-e
    fld2 = fld1
    vert_coord := "theta-e"
    fld2_intrp = wrf_user_vert_interp(input_file,fld2,vert_coord,interp_levels,opts)
    
    
    fld2_dims = dimsizes(fld2_intrp)
    
    filedimdef(fout, (/"fld2_time", "fld2_levs", "fld2_sn", "fld2_we"/), \
                     (/fld2_dims(0), fld2_dims(1), fld2_dims(2), fld2_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_tk_theta_e", typeof(fld2_intrp), (/"fld2_time", "fld2_levs", "fld2_sn", "fld2_we"/))
    filevarattdef(fout,"fld_tk_theta_e", fld2_intrp)
    fout->fld_tk_theta_e = (/fld2_intrp/)
    
    
    ; interp t to pressure
    fld3 = fld1
    vert_coord := "pressure"
    interp_levels    := (/850,500/)
    fld3_intrp = wrf_user_vert_interp(input_file,fld3,vert_coord,interp_levels,opts)
    
    
    fld3_dims = dimsizes(fld3_intrp)
    
    filedimdef(fout, (/"fld3_time", "fld3_levs", "fld3_sn", "fld3_we"/), \
                     (/fld3_dims(0), fld3_dims(1), fld3_dims(2), fld3_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_tk_pres", typeof(fld3_intrp), (/"fld3_time", "fld3_levs", "fld3_sn", "fld3_we"/))
    filevarattdef(fout,"fld_tk_pres", fld3_intrp)
    fout->fld_tk_pres = (/fld3_intrp/)
    
    ; interp t to ght_msl
    fld4 = fld1
    vert_coord := "ght_msl"
    interp_levels    := (/1,2/)
    fld4_intrp = wrf_user_vert_interp(input_file,fld4,vert_coord,interp_levels,opts)
    
    
    fld4_dims = dimsizes(fld4_intrp)
    
    filedimdef(fout, (/"fld4_time", "fld4_levs", "fld4_sn", "fld4_we"/), \
                     (/fld4_dims(0), fld4_dims(1), fld4_dims(2), fld4_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_tk_ght_msl", typeof(fld4_intrp), (/"fld4_time", "fld4_levs", "fld4_sn", "fld4_we"/))
    filevarattdef(fout,"fld_tk_ght_msl", fld4_intrp)
    fout->fld_tk_ght_msl = (/fld4_intrp/)
    
    
    ; interp t to ght_agl
    fld5 = fld1
    vert_coord := "ght_agl"
    interp_levels    := (/1,2/)
    fld5_intrp = wrf_user_vert_interp(input_file,fld1,vert_coord,interp_levels,opts)
    
    
    fld5_dims = dimsizes(fld5_intrp)
    
    filedimdef(fout, (/"fld5_time", "fld5_levs", "fld5_sn", "fld5_we"/), \
                     (/fld5_dims(0), fld5_dims(1), fld5_dims(2), fld5_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_tk_ght_agl", typeof(fld5_intrp), (/"fld5_time", "fld5_levs", "fld5_sn", "fld5_we"/))
    filevarattdef(fout,"fld_tk_ght_agl", fld5_intrp)
    fout->fld_tk_ght_agl = (/fld5_intrp/)
    
    ; interp ht to pres
    fld6 = wrf_user_getvar(input_file, "height", -1)
    vert_coord := "pressure"
    opts@field_type  = "ght"
    interp_levels    := (/500,50/)
    fld6_intrp = wrf_user_vert_interp(input_file,fld6,vert_coord,interp_levels,opts)
    
    
    fld6_dims = dimsizes(fld6_intrp)
    
    filedimdef(fout, (/"fld6_time", "fld6_levs", "fld6_sn", "fld6_we"/), \
                     (/fld6_dims(0), fld6_dims(1), fld6_dims(2), fld6_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_ht_pres", typeof(fld6_intrp), (/"fld6_time", "fld6_levs", "fld6_sn", "fld6_we"/))
    filevarattdef(fout,"fld_ht_pres", fld6_intrp)
    fout->fld_ht_pres = (/fld6_intrp/)
    
    
    ; interp pres to theta
    fld7 = wrf_user_getvar(input_file, "pressure", -1)
    vert_coord := "theta"
    opts@field_type  = "pressure"
    interp_levels    := (/200,300,500,1000/)
    fld7_intrp = wrf_user_vert_interp(input_file,fld7,vert_coord,interp_levels,opts)
    
    
    fld7_dims = dimsizes(fld7_intrp)
    
    filedimdef(fout, (/"fld7_time", "fld7_levs", "fld7_sn", "fld7_we"/), \
                     (/fld7_dims(0), fld7_dims(1), fld7_dims(2), fld7_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_pres_theta", typeof(fld7_intrp), (/"fld7_time", "fld7_levs", "fld7_sn", "fld7_we"/))
    filevarattdef(fout,"fld_pres_theta", fld7_intrp)
    fout->fld_pres_theta = (/fld7_intrp/)
    
    ; interp theta-e to pressure
    fld8 = wrf_user_getvar(input_file, "eth", -1)
    vert_coord := "pressure"
    opts@field_type  = "T"
    interp_levels    := (/850,500,5/)
    fld8_intrp = wrf_user_vert_interp(input_file,fld8,vert_coord,interp_levels,opts)
    
    
    fld8_dims = dimsizes(fld8_intrp)
    
    filedimdef(fout, (/"fld8_time", "fld8_levs", "fld8_sn", "fld8_we"/), \
                     (/fld8_dims(0), fld8_dims(1), fld8_dims(2), fld8_dims(3)/), \ 
                     (/False,False,False,False/))
                     
    filevardef(fout, "fld_thetae_pres", typeof(fld8_intrp), (/"fld8_time", "fld8_levs", "fld8_sn", "fld8_we"/))
    filevarattdef(fout,"fld_thetae_pres", fld8_intrp)
    fout->fld_thetae_pres = (/fld8_intrp/)
    
    
    delete(fout)







