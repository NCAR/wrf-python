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
    
    delete(fout)







