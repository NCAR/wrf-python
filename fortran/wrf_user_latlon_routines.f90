!NCLFORTSTART
SUBROUTINE ROTATECOORDS(ilat, ilon, olat, olon, lat_np, lon_np, lon_0, direction)
    USE wrf_constants, ONLY : PI, RAD_PER_DEG, DEG_PER_RAD

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: olat, olon

    REAL(KIND=8), INTENT(IN) :: ilat, ilon
    REAL(KIND=8), INTENT(OUT) :: olat, olon
    REAL(KIND=8), INTENT(IN) :: lat_np, lon_np, lon_0
    INTEGER, INTENT(IN) :: direction

! NCLFORTEND

    !  >=0, default : computational -> geographical
    !  < 0          : geographical  -> computational

    REAL(KIND=8) :: rlat, rlon
    REAL(KIND=8) :: phi_np, lam_np, lam_0, dlam
    REAL(KIND=8) :: sinphi, cosphi, coslam, sinlam
    !REAL(KIND=8), PARAMETER :: PI=3.141592653589793D0
    !REAL(KIND=8), PARAMETER :: RAD_PER_DEG=PI/180.D0
    !REAL(KIND=8), PARAMETER :: DEG_PER_RAD=180.D0/PI

    !convert all angles to radians
    phi_np = lat_np*RAD_PER_DEG
    lam_np = lon_np*RAD_PER_DEG
    lam_0 = lon_0*RAD_PER_DEG
    rlat = ilat*RAD_PER_DEG
    rlon = ilon*RAD_PER_DEG

    IF (direction .LT. 0) THEN
    ! the equations are exactly the same except for one
    ! small difference with respect to longitude ...
        dlam = pi - lam_0
    ELSE
        dlam = lam_np
    END IF

    sinphi = COS(phi_np)*COS(rlat)*COS(rlon - dlam) + SIN(phi_np)*SIN(rlat)
    cosphi = SQRT(1.D0 - sinphi*sinphi)
    coslam = SIN(phi_np)*COS(rlat)*COS(rlon - dlam) - COS(phi_np)*SIN(rlat)
    sinlam = COS(rlat)*SIN(rlon - dlam)

    IF (cosphi.NE.0.D0) THEN
        coslam = coslam/cosphi
        sinlam = sinlam/cosphi
    END IF

    olat = DEG_PER_RAD*ASIN(sinphi)
    olon = DEG_PER_RAD*(ATAN2(sinlam,coslam) - dlam - lam_0 + lam_np)

    RETURN

END SUBROUTINE ROTATECOORDS


!NCLFORTSTART
SUBROUTINE DLLTOIJ(map_proj, truelat1, truelat2, stdlon, lat1, lon1,&
                   pole_lat, pole_lon, knowni, knownj, dx, dy, latinc,&
                   loninc, lat, lon, loc, errstat, errmsg)
    USE wrf_constants, ONLY : ALGERR, PI, RAD_PER_DEG, DEG_PER_RAD, WRF_EARTH_RADIUS

    ! Converts input lat/lon values to the cartesian (i,j) value
    ! for the given projection.

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: loc


    INTEGER, INTENT(IN) :: map_proj
    REAL(KIND=8), INTENT(IN) :: stdlon
    REAL(KIND=8), INTENT(IN) ::lat1, lon1, pole_lat, pole_lon, knowni, knownj
    REAL(KIND=8), INTENT(IN) ::dx, dy, latinc, loninc
    REAL(KIND=8), INTENT(INOUT) :: lat, lon, truelat1, truelat2 ! these might get modified
    REAL(KIND=8), DIMENSION(2), INTENT(OUT) :: loc
    INTEGER, INTENT(INOUT) :: errstat
    CHARACTER(LEN=*), INTENT(INOUT) :: errmsg

!NCLEND

    REAL(KIND=8) :: deltalon1
    REAL(KIND=8) :: tl1r
    REAL(KIND=8) :: clain, dlon, rsw, deltalon, deltalat
    REAL(KIND=8) :: reflon, scale_top, ala1, alo1, ala, alo, rm, polei, polej
    ! earth radius divided by dx
    REAL(KIND=8) :: rebydx
    REAL(KIND=8) :: ctl1r, arg, cone, hemi
    REAL(KIND=8) :: i, j
    REAL(KIND=8) :: lat1n, lon1n, olat, olon

    ! Contants
    !REAL(KIND=8),PARAMETER :: PI=3.141592653589793D0
    !REAL(KIND=8),PARAMETER :: RAD_PER_DEG=PI/180.D0
    !REAL(KIND=8),PARAMETER :: DEG_PER_RAD=180.D0/PI
    !REAL(KIND=8),PARAMETER :: RE_M=6370000.D0 ! radius of earth in meters

    !      lat1     ! sw latitude (1,1) in degrees (-90->90n)
    !      lon1     ! sw longitude (1,1) in degrees (-180->180e)
    !      dx       ! grid spacing in meters at truelats
    !      dlat     ! lat increment for lat/lon grids
    !      dlon     ! lon increment for lat/lon grids
    !      stdlon   ! longitude parallel to y-axis (-180->180e)
    !      truelat1 ! first true latitude (all projections)
    !      truelat2 ! second true lat (lc only)
    !      hemi     ! 1 for nh, -1 for sh
    !      cone     ! cone factor for lc projections
    !      polei    ! computed i-location of pole point
    !      polej    ! computed j-location of pole point
    !      rsw      ! computed radius to sw corner
    !      knowni   ! x-location of known lat/lon
    !      knownj   ! y-location of known lat/lon
    !      re_m     ! radius of spherical earth, meters
    !      rebydx   ! earth radius divided by dx

    errstat = INT(dy) ! remove compiler warning since dy not used
    errstat = 0

    rebydx = WRF_EARTH_RADIUS/dx

    ! Get rid of compiler warnings
    i=0
    j=0

    hemi = 1.0D0
    IF (truelat1 .LT. 0.0D0) THEN
        hemi = -1.0D0
    END IF

    ! mercator
    IF (map_proj.EQ.3) THEN

        ! preliminary variables
        clain = COS(RAD_PER_DEG*truelat1)
        dlon = dx/(WRF_EARTH_RADIUS*clain)

        ! compute distance from equator to origin, and store in
        ! the rsw tag.
        rsw = 0.D0
        IF (lat1 .NE. 0.D0) THEN
            rsw = (LOG(TAN(0.5D0*((lat1 + 90.D0)*RAD_PER_DEG))))/dlon
        END IF

        deltalon = lon - lon1
        IF (deltalon .LT. -180.D0) deltalon = deltalon + 360.D0
        IF (deltalon .GT. 180.D0) deltalon = deltalon - 360.D0
        i = knowni + (deltalon/(dlon*DEG_PER_RAD))
        j = knownj + (LOG(TAN(0.5D0*((lat + 90.D0)*RAD_PER_DEG))))/dlon - rsw

    ! ps
    ELSE IF (map_proj .EQ. 2) THEN

        reflon = stdlon + 90.D0

        ! compute numerator term of map scale factor
        scale_top = 1.D0 + hemi*SIN(truelat1*RAD_PER_DEG)

        ! compute radius to lower-left (sw) corner
        ala1 = lat1*RAD_PER_DEG
        rsw = rebydx*COS(ala1)*scale_top/(1.D0 + hemi*SIN(ala1))

        ! find the pole point
        alo1 = (lon1 - reflon)*RAD_PER_DEG
        polei = knowni - rsw*COS(alo1)
        polej = knownj - hemi*rsw*SIN(alo1)

        ! find radius to desired point
        ala = lat*RAD_PER_DEG
        rm = rebydx*COS(ala)*scale_top/(1.D0 + hemi*SIN(ala))
        alo = (lon - reflon)*RAD_PER_DEG
        i = polei + rm*COS(alo)
        j = polej + hemi*rm*SIN(alo)

    ! lambert
    ELSE IF (map_proj .EQ. 1) THEN

        IF (ABS(truelat2) .GT. 90.D0) THEN
            truelat2 = truelat1
        END IF

        IF (ABS(truelat1 - truelat2) .GT. 0.1D0) THEN
            cone = (LOG(COS(truelat1*RAD_PER_DEG))-LOG(COS(truelat2*RAD_PER_DEG)))/&
                 (LOG(TAN((90.D0 - ABS(truelat1))*RAD_PER_DEG*0.5D0))-&
                 LOG(TAN((90.D0 - ABS(truelat2))*RAD_PER_DEG*0.5D0)))
        ELSE
            cone = SIN(ABS(truelat1)*RAD_PER_DEG)
        END IF

        ! compute longitude differences and ensure we stay
        ! out of the forbidden "cut zone"
        deltalon1 = lon1 - stdlon
        IF (deltalon1 .GT. +180.D0) deltalon1 = deltalon1 - 360.D0
        IF (deltalon1 .LT. -180.D0) deltalon1 = deltalon1 + 360.D0

        ! convert truelat1 to radian and compute cos for later use
        tl1r = truelat1*RAD_PER_DEG
        ctl1r = COS(tl1r)

        ! compute the radius to our known lower-left (sw) corner
        rsw = rebydx*ctl1r/cone*(TAN((90.D0*hemi - lat1)*RAD_PER_DEG/2.D0)/&
            TAN((90.D0*hemi - truelat1)*RAD_PER_DEG/2.D0))**cone

        ! find pole point
        arg = cone*(deltalon1*RAD_PER_DEG)
        polei = hemi*knowni - hemi*rsw*SIN(arg)
        polej = hemi*knownj + rsw*COS(arg)

        ! compute deltalon between known longitude and standard
        ! lon and ensure it is not in the cut zone
        deltalon = lon - stdlon
        IF (deltalon .GT. +180.D0) deltalon = deltalon - 360.D0
        IF (deltalon .LT. -180.D0) deltalon = deltalon + 360.D0

        ! radius to desired point
        rm = rebydx*ctl1r/cone*(TAN((90.D0*hemi - lat)*RAD_PER_DEG/2.D0)/&
                  TAN((90.D0*hemi - truelat1)*RAD_PER_DEG/2.D0))**cone

        arg = cone*(deltalon*RAD_PER_DEG)
        i = polei + hemi*rm*SIN(arg)
        j = polej - rm*COS(arg)

        ! finally, if we are in the southern hemisphere, flip the
        ! i/j values to a coordinate system where (1,1) is the sw
        ! corner (what we assume) which is different than the
        ! original ncep algorithms which used the ne corner as
        ! the origin in the southern hemisphere (left-hand vs.
        ! right-hand coordinate?)
        i = hemi*i
        j = hemi*j

    ! lat-lon
    ELSE IF (map_proj .EQ. 6) THEN

        IF (pole_lat .NE. 90.D0) THEN
            CALL ROTATECOORDS(lat, lon, olat, olon, pole_lat, pole_lon, stdlon, -1)
            lat = olat
            lon = olon + stdlon
        END IF

        ! make sure center lat/lon is good
        IF (pole_lat .NE. 90.D0) THEN
            CALL ROTATECOORDS(lat1, lon1, olat, olon, pole_lat, pole_lon, stdlon, -1)
            lat1n = olat
            lon1n = olon + stdlon
            deltalat = lat - lat1n
            deltalon = lon - lon1n
        ELSE
            deltalat = lat - lat1
            deltalon = lon - lon1
        END IF

        ! compute i/j
        i = deltalon/loninc
        j = deltalat/latinc

        i = i + knowni
        j = j + knownj

    ELSE
        errstat = ALGERR
        WRITE(errmsg, *) "Do not know map projection ", map_proj
        RETURN
    END IF

    loc(1) = j
    loc(2) = i

    RETURN

END SUBROUTINE DLLTOIJ


!NCLFORTSTART
SUBROUTINE DIJTOLL(map_proj, truelat1, truelat2, stdlon, lat1, lon1,&
                   pole_lat, pole_lon, knowni, knownj, dx, dy, latinc,&
                   loninc, ai, aj, loc, errstat, errmsg)
    USE wrf_constants, ONLY : ALGERR, PI, RAD_PER_DEG, DEG_PER_RAD, WRF_EARTH_RADIUS

    ! converts input lat/lon values to the cartesian (i,j) value
    ! for the given projection.
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: loc

    INTEGER, INTENT(IN) :: map_proj
    REAL(KIND=8), INTENT(IN) :: stdlon
    REAL(KIND=8), INTENT(IN) :: lat1, lon1, pole_lat, pole_lon, knowni, knownj
    REAL(KIND=8), INTENT(IN) :: dx, dy, latinc, loninc, ai, aj
    REAL(KIND=8), INTENT(INOUT) :: truelat1, truelat2
    REAL(KIND=8), DIMENSION(2), INTENT(OUT) :: loc
    INTEGER, INTENT(INOUT) :: errstat
    CHARACTER(LEN=*), INTENT(INOUT) :: errmsg

!NCLEND

    REAL(KIND=8) :: gi2
    REAL(KIND=8) :: arccos
    REAL(KIND=8) :: deltalon1
    REAL(KIND=8) :: tl1r
    REAL(KIND=8) :: clain, dlon, rsw, deltalon, deltalat
    REAL(KIND=8) :: reflon, scale_top, ala1, alo1, polei, polej
    ! earth radius divided by dx
    REAL(KIND=8) :: rebydx
    REAL(KIND=8) :: ctl1r, cone, hemi

    !REAL(KIND=8),PARAMETER :: PI = 3.141592653589793D0
    !REAL(KIND=8),PARAMETER :: RAD_PER_DEG = PI/180.D0
    !REAL(KIND=8),PARAMETER :: DEG_PER_RAD = 180.D0/PI
    !REAL(KIND=8),PARAMETER :: RE_M = 6370000.D0 ! radius of sperical earth

    REAL(KIND=8) :: inew, jnew, r, r2
    REAL(KIND=8) :: chi, chi1, chi2
    REAL(KIND=8) :: xx, yy, lat, lon

    REAL(KIND=8) :: olat, olon, lat1n, lon1n


    !     lat1     ! sw latitude (1,1) in degrees (-90->90n)
    !     lon1     ! sw longitude (1,1) in degrees (-180->180e)
    !     dx       ! grid spacing in meters at truelats
    !     dlat     ! lat increment for lat/lon grids
    !     dlon     ! lon increment for lat/lon grids
    !     stdlon   ! longitude parallel to y-axis (-180->180e)
    !     truelat1 ! first true latitude (all projections)
    !     truelat2 ! second true lat (lc only)
    !     hemi     ! 1 for nh, -1 for sh
    !     cone     ! cone factor for lc projections
    !     polei    ! computed i-location of pole point
    !     polej    ! computed j-location of pole point
    !     rsw      ! computed radius to sw corner
    !     knowni   ! x-location of known lat/lon
    !     knownj   ! y-location of known lat/lon
    !     re_m     ! radius of spherical earth, meters
    !     rebydx   ! earth radius divided by dx

    errstat = INT(dy) ! Remove compiler warning since dy not used
    errstat = 0

    rebydx = WRF_EARTH_RADIUS/dx

    hemi = 1.0D0
    IF (truelat1 .LT. 0.0D0) THEN
        hemi = -1.0D0
    END IF

    ! mercator
    IF (map_proj .EQ. 3) THEN
        ! preliminary variables
        clain = COS(RAD_PER_DEG*truelat1)
        dlon = dx/(WRF_EARTH_RADIUS*clain)

        ! compute distance from equator to origin, and store in
        ! the rsw tag.
        rsw = 0.D0
        IF (lat1 .NE. 0.D0) THEN
            rsw = (LOG(TAN(0.5D0*((lat1 + 90.D0)*RAD_PER_DEG))))/dlon
        END IF

        lat = 2.0D0*ATAN(EXP(dlon*(rsw + aj - knownj)))*DEG_PER_RAD - 90.D0
        lon = (ai - knowni)*dlon*DEG_PER_RAD + lon1
        IF (lon .GT. 180.D0) lon = lon - 360.D0
        IF (lon .LT. -180.D0) lon = lon + 360.D0

        ! ps
    ELSE IF (map_proj .EQ. 2) THEN

        ! compute the reference longitude by rotating 90 degrees to
        ! the east to find the longitude line parallel to the
        ! positive x-axis.
        reflon = stdlon + 90.D0

        ! compute numerator term of map scale factor
        scale_top = 1.D0 + hemi*SIN(truelat1*RAD_PER_DEG)

        ! compute radius to known point
        ala1 = lat1*RAD_PER_DEG
        rsw = rebydx*COS(ala1)*scale_top/(1.D0 + hemi*SIN(ala1))

        ! find the pole point
        alo1 = (lon1 - reflon)*RAD_PER_DEG
        polei = knowni - rsw*COS(alo1)
        polej = knownj - hemi*rsw*SIN(alo1)

        ! compute radius to point of interest
        xx = ai - polei
        yy = (aj - polej)*hemi
        r2 = xx**2 + yy**2

        ! now the magic code
        IF (r2 .EQ. 0.D0) THEN
            lat = hemi*90.D0
            lon = reflon
        ELSE
            gi2 = (rebydx*scale_top)**2.D0
            lat = DEG_PER_RAD*hemi*ASIN((gi2 - r2)/(gi2 + r2))
            arccos = ACOS(xx/SQRT(r2))
            IF (yy .GT. 0) THEN
                lon = reflon + DEG_PER_RAD*arccos
            ELSE
                lon = reflon - DEG_PER_RAD*arccos
            END IF
      END IF

        ! convert to a -180 -> 180 east convention
        IF (lon .GT. 180.D0) lon = lon - 360.D0
        IF (lon .LT. -180.D0) lon = lon + 360.D0

    !     !lambert
    ELSE IF (map_proj .EQ. 1) THEN

        IF (ABS(truelat2) .GT. 90.D0) THEN
            truelat2 = truelat1
        END IF

        IF (ABS(truelat1 - truelat2) .GT. 0.1D0) THEN
            cone = (LOG(COS(truelat1*RAD_PER_DEG)) - LOG(COS(truelat2*RAD_PER_DEG)))/&
                 (LOG(TAN((90.D0 - ABS(truelat1))*RAD_PER_DEG*0.5D0)) - &
                  LOG(TAN((90.D0 - ABS(truelat2))*RAD_PER_DEG*0.5D0)))
        ELSE
            cone = SIN(ABS(truelat1)*RAD_PER_DEG)
        END IF

        ! compute longitude differences and ensure we stay out of the
        ! forbidden "cut zone"
        deltalon1 = lon1 - stdlon
        IF (deltalon1 .GT. +180.D0) deltalon1 = deltalon1 - 360.D0
        IF (deltalon1 .LT. -180.D0) deltalon1 = deltalon1 + 360.D0

        ! convert truelat1 to radian and compute cos for later use
        tl1r = truelat1*RAD_PER_DEG
        ctl1r = COS(tl1r)

        ! compute the radius to our known point
        rsw = rebydx*ctl1r/cone*(TAN((90.D0*hemi - lat1)*RAD_PER_DEG/2.D0)/&
                         TAN((90.D0*hemi - truelat1)*RAD_PER_DEG/2.D0))**cone

        ! find pole point
        alo1 = cone*(deltalon1*RAD_PER_DEG)
        polei = hemi*knowni - hemi*rsw*SIN(alo1)
        polej = hemi*knownj + rsw*COS(alo1)

        chi1 = (90.D0 - hemi*truelat1)*RAD_PER_DEG
        chi2 = (90.D0 - hemi*truelat2)*RAD_PER_DEG

        ! see if we are in the southern hemispere and flip the
        ! indices if we are.
        inew = hemi*ai
        jnew = hemi*aj

        ! compute radius**2 to i/j location
        reflon = stdlon + 90.D0
        xx = inew - polei
        yy = polej - jnew
        r2 = (xx*xx + yy*yy)
        r = sqrt(r2)/rebydx

        ! convert to lat/lon
        IF (r2 .EQ. 0.D0) THEN
            lat = hemi*90.D0
            lon = stdlon
        ELSE
            lon = stdlon + DEG_PER_RAD*ATAN2(hemi*xx,yy)/cone
            lon = dmod(lon + 360.D0, 360.D0)
            IF (chi1 .EQ. chi2) THEN
                chi = 2.0D0*ATAN((r/TAN(chi1))**(1.D0/cone)*TAN(chi1*0.5D0))
            ELSE
                chi = 2.0D0*ATAN((r*cone/SIN(chi1))**(1.D0/cone)*TAN(chi1*0.5D0))
            END IF
            lat = (90.0D0 - chi*DEG_PER_RAD)*hemi
        END IF

        IF (lon .GT. +180.D0) lon = lon - 360.D0
        IF (lon .LT. -180.D0) lon = lon + 360.D0

    !     !lat-lon
    ELSE IF (map_proj .EQ. 6) THEN
        inew = ai - knowni
        jnew = aj - knownj

        IF (inew .LT. 0.D0) inew = inew + 360.D0/loninc
        IF (inew .GE. 360.D0/dx) inew = inew - 360.D0/loninc

        ! compute deltalat and deltalon
        deltalat = jnew*latinc
        deltalon = inew*loninc

        IF (pole_lat .NE. 90.D0) THEN
            CALL ROTATECOORDS(lat1, lon1, olat, olon, pole_lat, pole_lon, stdlon, -1)
            lat1n = olat
            lon1n = olon + stdlon
            lat = deltalat + lat1n
            lon = deltalon + lon1n
        ELSE
            lat = deltalat + lat1
            lon = deltalon + lon1
        END IF

        IF (pole_lat .NE. 90.D0) THEN
            lon = lon - stdlon
            CALL ROTATECOORDS(lat, lon, olat, olon, pole_lat, pole_lon, stdlon, 1)
            lat = olat
            lon = olon
        END IF

        IF (lon .LT. -180.D0) lon = lon + 360.D0
        IF (lon .GT. 180.D0) lon = lon - 360.D0

    ELSE
        errstat = ALGERR
        WRITE(errmsg, *) "Do not know map projection ", map_proj
        RETURN
    END IF

    loc(1) = lat
    loc(2) = lon

    RETURN

END SUBROUTINE DIJTOLL
