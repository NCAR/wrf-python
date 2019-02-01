
;---Open file and calculate slp.
a = addfile("/Users/ladwig/Documents/wrf_files/rotated_pole_test.nc","r")

t2 = wrf_user_getvar(a,"T2",0)

;---Start the graphics
wks = gsn_open_wks("x11","wrf")

;---Set some graphical resources
res             = True
res@gsnMaximize = True
res@cnFillOn    = True
res@tfDoNDCOverlay = True     ; This is necessary if you don't
                              ; set sfXArray/sfYArray

;---Add additional map resources
res = wrf_map_resources(a,res)

;---Change some of the resources that were set (these were set to "gray")
res@mpGeophysicalLineColor  = "black"
res@mpNationalLineColor     = "black"
res@mpUSStateLineColor      = "black"

plot = gsn_csm_contour_map(wks,t2,res)

