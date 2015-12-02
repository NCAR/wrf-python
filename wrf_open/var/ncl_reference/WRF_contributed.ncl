;*************************************************************************
; Note: several of the functions/procedures are used
;       to invoke old [ugly!] function names. 
;*************************************************************************
; D. Shea
;
;      convert WRF character variable "Times" to udunits
;      2001-06-11_12:00:00
;
;      convert WRF character variable "Times" to a coordinate variable  "Time"
;      opt can be "integer" or "string"
;      .   integer:  opt = 0 : hours since initial time: Times(0,:)
;      .             opt = 1 : hours since 1901-01-01 00:00:00
;      .   string:   opt = 'any udunits compatible string'
;
undef ("WRF_Times2Udunits_c")
function WRF_Times2Udunits_c(Times:character, opt)
local dimT, rank, year, month, day, hour, minute, sec, units, time
begin
   
  dimT   = dimsizes(Times)
  rank   = dimsizes(dimT)
  if (rank.ne.2) then
      print("===> WRF_contributed.ncl: WRF_Times2Udunits_c expects 2D array: rank="+rank)
      exit
  end if

  if (.not.(typeof(opt).eq."integer" .or. typeof(opt).eq."string")) then
      print("===> WRF_contributed.ncl: opt must be integer or string: type="+typeof(opt))
      exit
  end if

  year   = stringtointeger((/Times(:, 0:3) /))
  month  = stringtointeger((/Times(:, 5:6) /))
  day    = stringtointeger((/Times(:, 8:9) /))
  hour   = stringtointeger((/Times(:,11:12)/))
  minute = stringtointeger((/Times(:,14:15)/))
  sec    = stringtointeger((/Times(:,17:18)/))

  if (typeof(opt).eq."integer") then
      if (opt.eq.0) then
          units  = "hours since "+year(0)+"-" \
                                 +sprinti("%0.2i",month(0)) +"-" \
                                 +sprinti("%0.2i",day(0))   +" " \
                                 +sprinti("%0.2i",hour(0))  +":" \
                                 +sprinti("%0.2i",minute(0))+":" \
                                 +sprinti("%0.2i",sec(0))       
      else  
          units  = "hours since 1901-01-01 00:00:00"
      end if
  else
          units  = opt   ; opt is udunits compatible string
  end if

  Time   = ut_inv_calendar(year,month,day,hour,minute,sec, units, 0)

  Time!0          = "Time"
  Time@long_name  = "Time"
  Time@description= "Time"
  Time@units      =  units
  Time&Time       =  Time        ; make coordinate variable
  return (Time)
end 
;*************************************************************************
; D. Shea
;      interface to WRF_Times2Udunits_c
undef ("WRF_Times_to_udunits")
function WRF_Times_to_udunits(Times:character, opt)
begin
  return( WRF_Times2Udunits_c(Times, 0) )    ; old name 
end

;*************************************************************************
; D. Shea
;      convert WRF character variable "Times" to 
;      a coordinate variable of type double 
;      time(double) =              yyyymmddhhmnss
;      2001-06-11_12:00:00     ==> 20010611120000
;
;      opt: currently not used [dummy]
;
undef ("WRF_Times2double_c")
function WRF_Times2double_c(Times:character, opt)
local dimT, rank, N, time, i, tmp_c
begin
  dimT   = dimsizes(Times)
  rank   = dimsizes(dimT)
  if (rank.ne.2) then
      print("===> WRF_contributed.ncl: WRF_Times2Udunits_c expects 2D array: rank="+rank)
      exit
  end if

  N      = dimT(0)

  Time   = new( N ,"double") ; preset to "double"
  delete(Time@_FillValue)    ; coord variables should not have a _FillValue

  Time   = stringtointeger((/Times(:,0:3)/))  *1d10 + \  ; yyyy
           stringtointeger((/Times(:,5:6)/))  *1d8  + \  ; mm
           stringtointeger((/Times(:,8:9)/))  *1d6  + \  ; dd
           stringtointeger((/Times(:,11:12)/))*1d4  + \  ; hh
           stringtointeger((/Times(:,14:15)/))*1d2  + \  ; mn
           stringtointeger((/Times(:,17:18)/))*1d0       ; ss

  Time!0          = "Time"
  Time@long_name  = "Time"
  Time@description= "Time"
  Time@units      = "yyyymmddhhmnss"
  Time&Time       =  Time        ; make coordinate variable
  return (Time)
end 
;*************************************************************************
; D. Shea
;      interface to WRF_Times2double_c
;      more explicit function name
undef ("WRF_Times_to_ymdhms")
function WRF_Times_to_ymdhms(Times:character, opt)
begin
  return( WRF_Times2double_c(Times, 0) )    ; old name 
end

;*************************************************************************
; D. Shea
;      convert WRF character variable "Times" to 
;      a coordinate variable of type integer 
;      time(integer)=              yyyymmddhh   [->ymdh]
;      2001-06-11_12:00:00     ==> 2001061112
;
;      Note: mminute and second are not part of the returned time
;
;      opt: currently not used [dummy]
;
undef ("WRF_Times_to_ymdh")
function WRF_Times_to_ymdh(Times:character, opt)
local dimT, rank, N, time, i, tmp_c
begin
  dimT   = dimsizes(Times)
  rank   = dimsizes(dimT)
  if (rank.ne.2) then
      print("===> WRF_contributed.ncl: WRF_Times2yyyymmddhh expects 2D array: rank="+rank)
      exit
  end if

  N      = dimT(0)

  Time   = new( N ,"integer") 
  delete(Time@_FillValue)    ; coord variables should not have a _FillValue

  Time   = stringtointeger((/Times(:,0:3)/))  *1000000 + \  ; yyyy
           stringtointeger((/Times(:,5:6)/))  *10000   + \  ; mm
           stringtointeger((/Times(:,8:9)/))  *100     + \  ; dd
           stringtointeger((/Times(:,11:12)/))              ; hh

  Time!0          = "Time"
  Time@long_name  = "Time"
  Time@description= "Time"
  Time@units      = "yyyymmddhh"
  Time&Time       =  Time        ; make coordinate variable
  return (Time)
end 

;*************************************************************************
; D. Shea
; This is a driver that selects the appropriate 
; mapping function based upon the file attribute: MAP_PROJ
; MAP_PROJ=1 [Lambert Conformal]; =2 [Stereographic]; =3 [Mercator]
;
; opt: currently not used [potentail use: time counter for XLAT/XLONG]
;
; Sample usage:
;             ncdf = addfile("...", r")
;             res  = True
;             WRF_map_c (ncdf, res, 0)
;             res  = ...
;
undef("WRF_map_c")
procedure WRF_map_c (f:file, res:logical, opt)
local rank, dimll, nlat, mlon, lat2d, lon2d
begin
  if (isatt(f,"MAP_PROJ")) then
      if (f@MAP_PROJ.eq.1) then
          res@mpProjection = "LambertConformal"
      end if
      if (f@MAP_PROJ.eq.2) then
          res@mpProjection = "Stereographic"
       end if
      if (f@MAP_PROJ.eq.3) then
          res@mpProjection = "Mercator"
      end if
  else
      print ("WRF_mapProj: no MAP_PROJ attribute")
  end if

  rank  = dimsizes(filevardimsizes(f,"XLAT")) ; # of dimensions
  if (rank.eq.3) then
      lat2d = f->XLAT(0,:,:)       ; opt could bt "nt"  f->XLAT(opt,:,:)
      lon2d = f->XLONG(0,:,:)
  else
    if (rank.eq.2) then
        lat2d = f->XLAT
        lon2d = f->XLONG
    else
        print ("WRF_resLamCon_c: unexpected lat/lon rank: rank="+rank)
        exit
    end if
  end if

  lat2d@units = "degrees_north"    ; not needed
  lon2d@units = "degrees_east"

  dimll = dimsizes(lat2d)
  nlat  = dimll(0)
  mlon  = dimll(1)

  res@mpLimitMode            = "Corners"  
  res@mpLeftCornerLatF       = lat2d(0,0)
  res@mpLeftCornerLonF       = lon2d(0,0)
  res@mpRightCornerLatF      = lat2d(nlat-1,mlon-1)
  res@mpRightCornerLonF      = lon2d(nlat-1,mlon-1)

  res@mpCenterLonF           = f@CEN_LON
  res@mpCenterLatF           = f@CEN_LAT        ; default

  if (res@mpProjection.eq."Mercator") then
      res@mpCenterLatF = 0.0    ; Cindy Bruyere MMM/WRF 24 Mar 2006
  end if

  if (res@mpProjection.eq."LambertConformal") then
      res@mpLambertParallel1F    = f@TRUELAT1
      res@mpLambertParallel2F    = f@TRUELAT2
      if (isatt(f, "STAND_LON") ) then
          res@mpLambertMeridianF = f@STAND_LON  ; use if present
                                                ; CB MMM/WRF 4 Aug 2006
      else
          if (isatt(f, "CEN_LON") ) then
              res@mpLambertMeridianF    = f@CEN_LON
          else
              print("WRF_map_c: STAND_LON and CEN_LON missing")
          end if
      end if
  end if

  res@mpFillOn              = False            ; turn off map fill
  res@mpOutlineDrawOrder    = "PostDraw"       ; draw continental outline last
  res@mpOutlineBoundarySets = "GeophysicalAndUSStates" ; state boundaries
  res@mpPerimDrawOrder      = "PostDraw"       ; force map perim
                                               ; commented 5/17/2007 
;;res@tfDoNDCOverlay        = True             ; True for 'native' grid
                                               ; some WRF are not native
  res@gsnAddCyclic          = False            ; data are not cyclic
end 

;*************************************************************************
; D. Shea
; interface for backward compatibility
undef("WRF_resLamCon_c")
procedure WRF_resLamCon_c (f:file, res:logical, opt)
begin
  WRF_map_c (f, res, opt)
end

;*************************************************************************
; D. Shea
; interface for newly named procedure
undef("wrf_mapres_c")
procedure wrf_mapres_c(f:file, res:logical, opt)
begin
  WRF_map_c (f, res, opt)
end
;************************************************************************* 
; D. Shea                                                    
; single interface to convert WRF character variable "Times"
; to user specified numeric values               
;
; M. Haley
; At some point we decided to rename this from WRF_Times to wrf_times_c
; Also added error check for opt.
;
undef ("wrf_times_c")                             
function wrf_times_c(Times:character, opt:integer)
begin                                          
   if (opt.ge.0 .and. opt.le.1) then 
       return(WRF_Times2Udunits_c(Times, opt) ) 
   end if

   if (opt.eq.2) then
       return(WRF_Times2double_c(Times, opt) )
   end if                                    

   if (opt.eq.3) then 
       return(WRF_Times_to_ymdh(Times, opt) ) 
   end if  
end       
