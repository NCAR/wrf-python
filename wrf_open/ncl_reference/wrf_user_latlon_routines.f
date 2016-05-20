C NCLFORTSTART
      SUBROUTINE DLLTOIJ(MAP_PROJ,TRUELAT1,TRUELAT2,STDLON,LAT1,LON1,
     +                   POLE_LAT,POLE_LON,KNOWNI,KNOWNJ,DX,DY,LATINC,
     +                   LONINC,LAT,LON,LOC)
      DOUBLE PRECISION DELTALON1
      DOUBLE PRECISION TL1R


ccc       Converts input lat/lon values to the cartesian (i,j) value
ccc       for the given projection.

      INTEGER MAP_PROJ
      DOUBLE PRECISION TRUELAT1,TRUELAT2,STDLON
      DOUBLE PRECISION LAT1,LON1,POLE_LAT,POLE_LON,KNOWNI,KNOWNJ
      DOUBLE PRECISION DX,DY,LATINC,LONINC,LAT,LON,LOC(2)
C NCLEND

      DOUBLE PRECISION CLAIN,DLON,RSW,DELTALON,DELTALAT
      DOUBLE PRECISION REFLON,SCALE_TOP,ALA1,ALO1,ALA,ALO,RM,POLEI,POLEJ
c Earth radius divided by dx
      DOUBLE PRECISION REBYDX
      DOUBLE PRECISION DELTALON1TL1R,CTL1R,ARG,CONE,HEMI
      DOUBLE PRECISION I,J
      DOUBLE PRECISION LAT1N,LON1N,OLAT,OLON

      DOUBLE PRECISION PI,RAD_PER_DEG,DEG_PER_RAD,RE_M

ccc      lat1     ! SW latitude (1,1) in degrees (-90->90N)
ccc      lon1     ! SW longitude (1,1) in degrees (-180->180E)
ccc      dx       ! Grid spacing in meters at truelats
ccc      dlat     ! Lat increment for lat/lon grids
ccc      dlon     ! Lon increment for lat/lon grids
ccc      stdlon   ! Longitude parallel to y-axis (-180->180E)
ccc      truelat1 ! First true latitude (all projections)
ccc      truelat2 ! Second true lat (LC only)
ccc      hemi     ! 1 for NH, -1 for SH
ccc      cone     ! Cone factor for LC projections
ccc      polei    ! Computed i-location of pole point
ccc      polej    ! Computed j-location of pole point
ccc      rsw      ! Computed radius to SW corner
ccc      knowni   ! X-location of known lat/lon
ccc      knownj   ! Y-location of known lat/lon
ccc      RE_M     ! Radius of spherical earth, meters
ccc      REbydx   ! Earth radius divided by dx

      PI = 3.141592653589793D0
      RAD_PER_DEG = PI/180.D0
      DEG_PER_RAD = 180.D0/PI
c Radius of spherical earth, meters
      RE_M = 6370000.D0
      REBYDX = RE_M/DX

      HEMI = 1.0D0
      IF (TRUELAT1.LT.0.0D0) THEN
          HEMI = -1.0D0
      END IF


ccc      !MERCATOR
      IF (MAP_PROJ.EQ.3) THEN

ccc         !  Preliminary variables
          CLAIN = COS(RAD_PER_DEG*TRUELAT1)
          DLON = DX/ (RE_M*CLAIN)

ccc         ! Compute distance from equator to origin, and store in 
ccc         ! the rsw tag.
          RSW = 0.D0
          IF (LAT1.NE.0.D0) THEN
              RSW = (DLOG(TAN(0.5D0* ((LAT1+90.D0)*RAD_PER_DEG))))/DLON
          END IF

          DELTALON = LON - LON1
          IF (DELTALON.LT.-180.D0) DELTALON = DELTALON + 360.D0
          IF (DELTALON.GT.180.D0) DELTALON = DELTALON - 360.D0
          I = KNOWNI + (DELTALON/ (DLON*DEG_PER_RAD))
          J = KNOWNJ + (DLOG(TAN(0.5D0* ((LAT+90.D0)*RAD_PER_DEG))))/
     +        DLON - RSW

ccc      !PS
      ELSE IF (MAP_PROJ.EQ.2) THEN

          REFLON = STDLON + 90.D0

ccc         ! Compute numerator term of map scale factor
          SCALE_TOP = 1.D0 + HEMI*SIN(TRUELAT1*RAD_PER_DEG)

ccc         ! Compute radius to lower-left (SW) corner
          ALA1 = LAT1*RAD_PER_DEG
          RSW = REBYDX*COS(ALA1)*SCALE_TOP/ (1.D0+HEMI*SIN(ALA1))

ccc         ! Find the pole point
          ALO1 = (LON1-REFLON)*RAD_PER_DEG
          POLEI = KNOWNI - RSW*COS(ALO1)
          POLEJ = KNOWNJ - HEMI*RSW*SIN(ALO1)

ccc         ! Find radius to desired point
          ALA = LAT*RAD_PER_DEG
          RM = REBYDX*COS(ALA)*SCALE_TOP/ (1.D0+HEMI*SIN(ALA))
          ALO = (LON-REFLON)*RAD_PER_DEG
          I = POLEI + RM*COS(ALO)
          J = POLEJ + HEMI*RM*SIN(ALO)

ccc      !LAMBERT
      ELSE IF (MAP_PROJ.EQ.1) THEN

          IF (ABS(TRUELAT2).GT.90.D0) THEN
              TRUELAT2 = TRUELAT1
          END IF

          IF (ABS(TRUELAT1-TRUELAT2).GT.0.1D0) THEN
              CONE = (DLOG(COS(TRUELAT1*RAD_PER_DEG))-
     +               DLOG(COS(TRUELAT2*RAD_PER_DEG)))/
     +               (DLOG(TAN((90.D0-ABS(TRUELAT1))*RAD_PER_DEG*
     +               0.5D0))-DLOG(TAN((90.D0-ABS(TRUELAT2))*RAD_PER_DEG*
     +               0.5D0)))
          ELSE
              CONE = SIN(ABS(TRUELAT1)*RAD_PER_DEG)
          END IF

ccc         ! Compute longitude differences and ensure we stay
ccc         ! out of the forbidden "cut zone"
          DELTALON1 = LON1 - STDLON
          IF (DELTALON1.GT.+180.D0) DELTALON1 = DELTALON1 - 360.D0
          IF (DELTALON1.LT.-180.D0) DELTALON1 = DELTALON1 + 360.D0

ccc         ! Convert truelat1 to radian and compute COS for later use
          TL1R = TRUELAT1*RAD_PER_DEG
          CTL1R = COS(TL1R)

ccc         ! Compute the radius to our known lower-left (SW) corner
          RSW = REBYDX*CTL1R/CONE* (TAN((90.D0*HEMI-
     +          LAT1)*RAD_PER_DEG/2.D0)/TAN((90.D0*HEMI-
     +          TRUELAT1)*RAD_PER_DEG/2.D0))**CONE

ccc         ! Find pole point
          ARG = CONE* (DELTALON1*RAD_PER_DEG)
          POLEI = HEMI*KNOWNI - HEMI*RSW*SIN(ARG)
          POLEJ = HEMI*KNOWNJ + RSW*COS(ARG)

ccc         ! Compute deltalon between known longitude and standard
ccc         ! lon and ensure it is not in the cut zone
          DELTALON = LON - STDLON
          IF (DELTALON.GT.+180.D0) DELTALON = DELTALON - 360.D0
          IF (DELTALON.LT.-180.D0) DELTALON = DELTALON + 360.D0

ccc         ! Radius to desired point
          RM = REBYDX*CTL1R/CONE* (TAN((90.D0*HEMI-
     +         LAT)*RAD_PER_DEG/2.D0)/TAN((90.D0*HEMI-
     +         TRUELAT1)*RAD_PER_DEG/2.D0))**CONE

          ARG = CONE* (DELTALON*RAD_PER_DEG)
          I = POLEI + HEMI*RM*SIN(ARG)
          J = POLEJ - RM*COS(ARG)

ccc         ! Finally, if we are in the southern hemisphere, flip the
ccc         ! i/j values to a coordinate system where (1,1) is the SW
ccc         ! corner (what we assume) which is different than the
ccc         ! original NCEP algorithms which used the NE corner as
ccc         ! the origin in the southern hemisphere (left-hand vs.
ccc         ! right-hand coordinate?)
          I = HEMI*I
          J = HEMI*J


ccc     !lat-lon
      ELSE IF (MAP_PROJ.EQ.6) THEN

          IF (POLE_LAT.NE.90.D0) THEN
              CALL ROTATECOORDS(LAT,LON,OLAT,OLON,POLE_LAT,POLE_LON,
     +                          STDLON,-1)
              LAT = OLAT
              LON = OLON + STDLON
          END IF

c         ! make sure center lat/lon is good
          IF (POLE_LAT.NE.90.D0) THEN
              CALL ROTATECOORDS(LAT1,LON1,OLAT,OLON,POLE_LAT,POLE_LON,
     +                          STDLON,-1)
              LAT1N = OLAT
              LON1N = OLON + STDLON
              DELTALAT = LAT - LAT1N
              DELTALON = LON - LON1N
          ELSE
              DELTALAT = LAT - LAT1
              DELTALON = LON - LON1
          END IF

c         ! Compute i/j
          I = DELTALON/LONINC
          J = DELTALAT/LATINC

          I = I + KNOWNI
          J = J + KNOWNJ

      ELSE

          PRINT *,'ERROR: Do not know map projection ',MAP_PROJ

      END IF

      LOC(1) = J
      LOC(2) = I

      RETURN
      END


C NCLFORTSTART
      SUBROUTINE DIJTOLL(MAP_PROJ,TRUELAT1,TRUELAT2,STDLON,LAT1,LON1,
     +                   POLE_LAT,POLE_LON,KNOWNI,KNOWNJ,DX,DY,LATINC,
     +                   LONINC,AI,AJ,LOC)
      DOUBLE PRECISION GI2
      DOUBLE PRECISION ARCCOS
      DOUBLE PRECISION DELTALON1
      DOUBLE PRECISION TL1R

ccc     ! Converts input lat/lon values to the cartesian (i,j) value
ccc     ! for the given projection.

      INTEGER MAP_PROJ
      DOUBLE PRECISION TRUELAT1,TRUELAT2,STDLON
      DOUBLE PRECISION LAT1,LON1,POLE_LAT,POLE_LON,KNOWNI,KNOWNJ
      DOUBLE PRECISION DX,DY,LATINC,LONINC,AI,AJ,LOC(2)
C NCLEND

      DOUBLE PRECISION CLAIN,DLON,RSW,DELTALON,DELTALAT
      DOUBLE PRECISION REFLON,SCALE_TOP,ALA1,ALO1,ALA,ALO,RM,POLEI,POLEJ
c Earth radius divided by dx
      DOUBLE PRECISION REBYDX
      DOUBLE PRECISION DELTALON1TL1R,CTL1R,ARG,CONE,HEMI

      DOUBLE PRECISION PI,RAD_PER_DEG,DEG_PER_RAD,RE_M

      DOUBLE PRECISION INEW,JNEW,R,R2
      DOUBLE PRECISION CHI,CHI1,CHI2
      DOUBLE PRECISION XX,YY,LAT,LON

      DOUBLE PRECISION RLAT,RLON,OLAT,OLON,LAT1N,LON1N
      DOUBLE PRECISION PHI_NP,LAM_NP,LAM_0,DLAM
      DOUBLE PRECISION SINPHI,COSPHI,COSLAM,SINLAM


ccc     lat1     ! SW latitude (1,1) in degrees (-90->90N)
ccc     lon1     ! SW longitude (1,1) in degrees (-180->180E)
ccc     dx       ! Grid spacing in meters at truelats
ccc     dlat     ! Lat increment for lat/lon grids
ccc     dlon     ! Lon increment for lat/lon grids
ccc     stdlon   ! Longitude parallel to y-axis (-180->180E)
ccc     truelat1 ! First true latitude (all projections)
ccc     truelat2 ! Second true lat (LC only)
ccc     hemi     ! 1 for NH, -1 for SH
ccc     cone     ! Cone factor for LC projections
ccc     polei    ! Computed i-location of pole point
ccc     polej    ! Computed j-location of pole point
ccc     rsw      ! Computed radius to SW corner
ccc     knowni   ! X-location of known lat/lon
ccc     knownj   ! Y-location of known lat/lon
ccc     RE_M     ! Radius of spherical earth, meters
ccc     REbydx   ! Earth radius divided by dx

      PI = 3.141592653589793D0
      RAD_PER_DEG = PI/180.D0
      DEG_PER_RAD = 180.D0/PI
c Radius of spherical earth, meters
      RE_M = 6370000.D0
      REBYDX = RE_M/DX

      HEMI = 1.0D0
      IF (TRUELAT1.LT.0.0D0) THEN
          HEMI = -1.0D0
      END IF


ccc     !MERCATOR
      IF (MAP_PROJ.EQ.3) THEN

ccc       !  Preliminary variables
          CLAIN = COS(RAD_PER_DEG*TRUELAT1)
          DLON = DX/ (RE_M*CLAIN)

ccc       ! Compute distance from equator to origin, and store in 
ccc       ! the rsw tag.
          RSW = 0.D0
          IF (LAT1.NE.0.D0) THEN
              RSW = (DLOG(TAN(0.5D0* ((LAT1+90.D0)*RAD_PER_DEG))))/DLON
          END IF

          LAT = 2.0D0*ATAN(EXP(DLON* (RSW+AJ-KNOWNJ)))*DEG_PER_RAD -
     +          90.D0
          LON = (AI-KNOWNI)*DLON*DEG_PER_RAD + LON1
          IF (LON.GT.180.D0) LON = LON - 360.D0
          IF (LON.LT.-180.D0) LON = LON + 360.D0


ccc     !PS
      ELSE IF (MAP_PROJ.EQ.2) THEN

ccc       ! Compute the reference longitude by rotating 90 degrees to
ccc       ! the east to find the longitude line parallel to the 
ccc       ! positive x-axis.
          REFLON = STDLON + 90.D0

ccc       ! Compute numerator term of map scale factor
          SCALE_TOP = 1.D0 + HEMI*SIN(TRUELAT1*RAD_PER_DEG)

ccc       ! Compute radius to known point
          ALA1 = LAT1*RAD_PER_DEG
          RSW = REBYDX*COS(ALA1)*SCALE_TOP/ (1.D0+HEMI*SIN(ALA1))

ccc       ! Find the pole point
          ALO1 = (LON1-REFLON)*RAD_PER_DEG
          POLEI = KNOWNI - RSW*COS(ALO1)
          POLEJ = KNOWNJ - HEMI*RSW*SIN(ALO1)

ccc       ! Compute radius to point of interest
          XX = AI - POLEI
          YY = (AJ-POLEJ)*HEMI
          R2 = XX**2 + YY**2

ccc       ! Now the magic code
          IF (R2.EQ.0.D0) THEN
              LAT = HEMI*90.D0
              LON = REFLON
          ELSE
              GI2 = (REBYDX*SCALE_TOP)**2.D0
              LAT = DEG_PER_RAD*HEMI*ASIN((GI2-R2)/ (GI2+R2))
              ARCCOS = ACOS(XX/SQRT(R2))
              IF (YY.GT.0) THEN
                  LON = REFLON + DEG_PER_RAD*ARCCOS
              ELSE
                  LON = REFLON - DEG_PER_RAD*ARCCOS
              END IF
          END IF

ccc       ! Convert to a -180 -> 180 East convention
          IF (LON.GT.180.D0) LON = LON - 360.D0
          IF (LON.LT.-180.D0) LON = LON + 360.D0

ccc     !LAMBERT
      ELSE IF (MAP_PROJ.EQ.1) THEN

          IF (ABS(TRUELAT2).GT.90.D0) THEN
              TRUELAT2 = TRUELAT1
          END IF

          IF (ABS(TRUELAT1-TRUELAT2).GT.0.1D0) THEN
              CONE = (DLOG(COS(TRUELAT1*RAD_PER_DEG))-
     +               DLOG(COS(TRUELAT2*RAD_PER_DEG)))/
     +               (DLOG(TAN((90.D0-ABS(TRUELAT1))*RAD_PER_DEG*
     +               0.5D0))-DLOG(TAN((90.D0-ABS(TRUELAT2))*RAD_PER_DEG*
     +               0.5D0)))
          ELSE
              CONE = SIN(ABS(TRUELAT1)*RAD_PER_DEG)
          END IF

ccc       ! Compute longitude differences and ensure we stay out of the
ccc       ! forbidden "cut zone"
          DELTALON1 = LON1 - STDLON
          IF (DELTALON1.GT.+180.D0) DELTALON1 = DELTALON1 - 360.D0
          IF (DELTALON1.LT.-180.D0) DELTALON1 = DELTALON1 + 360.D0

ccc       ! Convert truelat1 to radian and compute COS for later use
          TL1R = TRUELAT1*RAD_PER_DEG
          CTL1R = COS(TL1R)

ccc       ! Compute the radius to our known point
          RSW = REBYDX*CTL1R/CONE* (TAN((90.D0*HEMI-
     +          LAT1)*RAD_PER_DEG/2.D0)/TAN((90.D0*HEMI-
     +          TRUELAT1)*RAD_PER_DEG/2.D0))**CONE

ccc       ! Find pole point
          ALO1 = CONE* (DELTALON1*RAD_PER_DEG)
          POLEI = HEMI*KNOWNI - HEMI*RSW*SIN(ALO1)
          POLEJ = HEMI*KNOWNJ + RSW*COS(ALO1)

          CHI1 = (90.D0-HEMI*TRUELAT1)*RAD_PER_DEG
          CHI2 = (90.D0-HEMI*TRUELAT2)*RAD_PER_DEG

ccc       ! See if we are in the southern hemispere and flip the 
ccc       ! indices if we are.
          INEW = HEMI*AI
          JNEW = HEMI*AJ

ccc       ! Compute radius**2 to i/j location
          REFLON = STDLON + 90.D0
          XX = INEW - POLEI
          YY = POLEJ - JNEW
          R2 = (XX*XX+YY*YY)
          R = SQRT(R2)/REBYDX

ccc       ! Convert to lat/lon
          IF (R2.EQ.0.D0) THEN
              LAT = HEMI*90.D0
              LON = STDLON
          ELSE
              LON = STDLON + DEG_PER_RAD*ATAN2(HEMI*XX,YY)/CONE
              LON = DMOD(LON+360.D0,360.D0)
              IF (CHI1.EQ.CHI2) THEN
                  CHI = 2.0D0*ATAN((R/TAN(CHI1))** (1.D0/CONE)*
     +                  TAN(CHI1*0.5D0))
              ELSE
                  CHI = 2.0D0*ATAN((R*CONE/SIN(CHI1))** (1.D0/CONE)*
     +                  TAN(CHI1*0.5D0))
              END IF
              LAT = (90.0D0-CHI*DEG_PER_RAD)*HEMI
          END IF

          IF (LON.GT.+180.D0) LON = LON - 360.D0
          IF (LON.LT.-180.D0) LON = LON + 360.D0


ccc     !lat-lon
      ELSE IF (MAP_PROJ.EQ.6) THEN

          INEW = AI - KNOWNI
          JNEW = AJ - KNOWNJ

          IF (INEW.LT.0.D0) INEW = INEW + 360.D0/LONINC
          IF (INEW.GE.360.D0/DX) INEW = INEW - 360.D0/LONINC
c
ccc       ! Compute deltalat and deltalon
          DELTALAT = JNEW*LATINC
          DELTALON = INEW*LONINC

          IF (POLE_LAT.NE.90.D0) THEN
              CALL ROTATECOORDS(LAT1,LON1,OLAT,OLON,POLE_LAT,POLE_LON,
     +                          STDLON,-1)
              LAT1N = OLAT
              LON1N = OLON + STDLON
              LAT = DELTALAT + LAT1N
              LON = DELTALON + LON1N
          ELSE
              LAT = DELTALAT + LAT1
              LON = DELTALON + LON1
          END IF


          IF (POLE_LAT.NE.90.D0) THEN
              LON = LON - STDLON
              CALL ROTATECOORDS(LAT,LON,OLAT,OLON,POLE_LAT,POLE_LON,
     +                          STDLON,1)
              LAT = OLAT
              LON = OLON
          END IF

          IF (LON.LT.-180.D0) LON = LON + 360.D0
          IF (LON.GT.180.D0) LON = LON - 360.D0

      ELSE

          PRINT *,'ERROR: Do not know map projection ',MAP_PROJ

      END IF

      LOC(1) = LAT
      LOC(2) = LON
      RETURN

      END


C NCLFORTSTART
      SUBROUTINE ROTATECOORDS(ILAT,ILON,OLAT,OLON,LAT_NP,LON_NP,LON_0,
     +                        DIRECTION)
      DOUBLE PRECISION ILAT,ILON
      DOUBLE PRECISION OLAT,OLON
      DOUBLE PRECISION LAT_NP,LON_NP,LON_0
      INTEGER DIRECTION
C NCLEND

c       ! >=0, default : computational -> geographical
c       ! < 0          : geographical  -> computational

      DOUBLE PRECISION RLAT,RLON
      DOUBLE PRECISION PHI_NP,LAM_NP,LAM_0,DLAM
      DOUBLE PRECISION SINPHI,COSPHI,COSLAM,SINLAM
      DOUBLE PRECISION PI,RAD_PER_DEG,DEG_PER_RAD

      PI = 3.141592653589793D0
      RAD_PER_DEG = PI/180.D0
      DEG_PER_RAD = 180.D0/PI

c       ! Convert all angles to radians
      PHI_NP = LAT_NP*RAD_PER_DEG
      LAM_NP = LON_NP*RAD_PER_DEG
      LAM_0 = LON_0*RAD_PER_DEG
      RLAT = ILAT*RAD_PER_DEG
      RLON = ILON*RAD_PER_DEG

      IF (DIRECTION.LT.0) THEN
c          ! The equations are exactly the same except for one
c          ! small difference with respect to longitude ...
          DLAM = PI - LAM_0
      ELSE
          DLAM = LAM_NP
      END IF
      SINPHI = COS(PHI_NP)*COS(RLAT)*COS(RLON-DLAM) +
     +         SIN(PHI_NP)*SIN(RLAT)
      COSPHI = SQRT(1.D0-SINPHI*SINPHI)
      COSLAM = SIN(PHI_NP)*COS(RLAT)*COS(RLON-DLAM) -
     +         COS(PHI_NP)*SIN(RLAT)
      SINLAM = COS(RLAT)*SIN(RLON-DLAM)
      IF (COSPHI.NE.0.D0) THEN
          COSLAM = COSLAM/COSPHI
          SINLAM = SINLAM/COSPHI
      END IF
      OLAT = DEG_PER_RAD*ASIN(SINPHI)
      OLON = DEG_PER_RAD* (ATAN2(SINLAM,COSLAM)-DLAM-LAM_0+LAM_NP)

      END
