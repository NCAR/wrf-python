input_file = addfile("/Users/ladwig/Documents/wrf_files/wrfout_d02_2010-06-13_21:00:00.nc", "r")

z  = wrf_user_getvar(input_file, "z", 0)        ; grid point height
p  = wrf_user_getvar(input_file, "pressure", 0) ; total pressure

dimsz = dimsizes(z)
pivot = (/ (dimsz(2)-1)/2, (dimsz(1)-1)/2 /)    ; pivot point is center of domain

; For the new cross section routine
xopt = True
xopt@use_pivot = True
xopt@angle = 45.0
;xopt@levels = 
;xopt@latlon = 
xopt@file_handle = input_file
;xopt@timeidx = 
xopt@linecoords = True
    
ht_vertcross = wrf_user_vertcross(z, p, pivot, xopt)

printVarSummary(ht_vertcross)
print(min(ht_vertcross@lats))
print(min(ht_vertcross@lons))
print(max(ht_vertcross@lats))
print(max(ht_vertcross@lons))


xopt@use_pivot = False
xopt@angle = 0.0
;xopt@levels = 
xopt@latlon = True
xopt@file_handle = input_file
xopt@timeidx = 0 
xopt@linecoords = True

loc_param = (/-104.3632, 32.8562, -95.15308, 40.06575 /)    ; pivot point is center of domain
ht_vertcross2 = wrf_user_vertcross(z, p, loc_param, xopt)

printVarSummary(ht_vertcross2)
print(min(ht_vertcross2@lats))
print(min(ht_vertcross2@lons))
print(max(ht_vertcross2@lats))
print(max(ht_vertcross2@lons))

print(ht_vertcross2@lats(190))
print(ht_vertcross2@lons(190))

xopt@use_pivot = True
xopt@angle = 45.0
;xopt@levels = 
xopt@latlon = True
xopt@file_handle = input_file
xopt@timeidx = 0 
xopt@linecoords = True

loc_param := (/-99.98572, 36.54949 /)    ; pivot point is center of domain
ht_vertcross3 = wrf_user_vertcross(z, p, loc_param, xopt)

printVarSummary(ht_vertcross3)
print(min(ht_vertcross3@lats))
print(min(ht_vertcross3@lons))
print(max(ht_vertcross3@lats))
print(max(ht_vertcross3@lons))


xopt@use_pivot = True
xopt@angle = 45.0
xopt@levels = (/1000., 850., 700., 500., 250. /)
xopt@latlon = True
xopt@file_handle = input_file
xopt@timeidx = 0 
xopt@linecoords = True

loc_param := (/-99.98572, 36.54949 /)    ; pivot point is center of domain
ht_vertcross4 = wrf_user_vertcross(z, p, loc_param, xopt)

printVarSummary(ht_vertcross4)
print(min(ht_vertcross4@lats))
print(min(ht_vertcross4@lons))
print(max(ht_vertcross4@lats))
print(max(ht_vertcross4@lons))

o = True
o@returnInt = False
o@useTime = 0
l = wrf_user_ll_to_xy(input_file, -99.98572, 36.54949, o)
print(l)


l1 = wrf_user_xy_to_ll(input_file, l(1), l(0), o)
print(l1)

