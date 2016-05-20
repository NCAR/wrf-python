    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
    
    in_files = (/"/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-28_00:00:00", \
                 "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-28_12:00:00", \
                 "/Users/ladwig/Documents/wrf_files/wrf_vortex_multi/wrfout_d02_2005-08-29_00:00:00"/)
    
    input_file = addfiles(in_files, "r")
    
    timeidx = 6
    
    ; 3D interpolation to new vertical coordinates
    ; interp t to theta
    fld1 = wrf_user_getvar(input_file, "tk", timeidx)
    vert_coord       = "theta"
    ;interp_levels    = (/200,300,500,1000/)
     interp_levels = ispan(200,1000,50)
    opts             = True
    opts@extrapolate = True 
    opts@field_type  = "T"
    opts@logP        = True
    opts@time        = timeidx  
    
    fld1_intrp = wrf_user_vert_interp(input_file,fld1,vert_coord,interp_levels,opts)
    
    printVarSummary(fld1_intrp)
    
    ; interp t to theta-e
    fld2 = fld1
    vert_coord := "theta-e"
    fld2_intrp = wrf_user_vert_interp(input_file,fld2,vert_coord,interp_levels,opts)
    
    printVarSummary(fld2_intrp)
    
    ; interp t to pressure
    fld3 = fld1
    vert_coord := "pressure"
    interp_levels    := (/850,500/)
    fld3_intrp = wrf_user_vert_interp(input_file,fld3,vert_coord,interp_levels,opts)
    
    printVarSummary(fld3_intrp)
    
    ; interp t to ght_msl
    fld4 = fld1
    vert_coord := "ght_msl"
    interp_levels    := (/1,2/)
    fld4_intrp = wrf_user_vert_interp(input_file,fld4,vert_coord,interp_levels,opts)
    
    printVarSummary(fld4_intrp)
    
    ; interp t to ght_agl
    fld5 = fld1
    vert_coord := "ght_agl"
    interp_levels    := (/1,2/)
    fld5_intrp = wrf_user_vert_interp(input_file,fld1,vert_coord,interp_levels,opts)
    
    printVarSummary(fld5_intrp)
    
    ; interp ht to pres
    fld6 = wrf_user_getvar(input_file, "height", timeidx)
    vert_coord := "pressure"
    opts@field_type  = "ght"
    interp_levels    := (/500,50/)
    fld6_intrp = wrf_user_vert_interp(input_file,fld6,vert_coord,interp_levels,opts)
    
    printVarSummary(fld6_intrp)
    
    ; interp pres to theta
    fld7 = wrf_user_getvar(input_file, "pressure", timeidx)
    vert_coord := "theta"
    opts@field_type  = "pressure"
    interp_levels    := (/200,300,500,1000/)
    fld7_intrp = wrf_user_vert_interp(input_file,fld7,vert_coord,interp_levels,opts)
    
    printVarSummary(fld7_intrp)
    
    ; interp theta-e to pressure
    fld8 = wrf_user_getvar(input_file, "eth", timeidx)
    vert_coord := "pressure"
    opts@field_type  = "T"
    interp_levels    := (/850,500,5/)
    fld8_intrp = wrf_user_vert_interp(input_file,fld8,vert_coord,interp_levels,opts)
    
    printVarSummary(fld8_intrp)
    
    