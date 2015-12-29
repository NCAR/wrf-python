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
                "wa", "uvmet10", "uvmet", "z"/]
    
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
    
    delete(fout)







