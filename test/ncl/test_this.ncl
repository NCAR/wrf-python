ifils = systemfunc ("ls /Users/ladwig/Documents/wrf_files/wrf_vortex_single/wrfout_d02_*")
print(ifils)
a = addfiles(ifils, "r")
;a = addfile("/Users/ladwig/Documents/wrf_files/wrfout_d01_2010-06-13_21:00:00.nc", "r")

lats := (/22.0, 25.0, 27.0 /)
lons := (/-90.0, -87.5, -83.75 /)

opt = True
opt@useTime = -1 
opt@returnI = False
xy = wrf_user_ll_to_xy(a, lons, lats, opt)

print(xy)

x_s = (/10, 50, 90 /)
y_s = (/10, 50, 90 /)

ll = wrf_user_xy_to_ll(a, x_s, y_s, opt)
print(ll)

