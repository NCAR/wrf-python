C NCLFORTSTART
      SUBROUTINE DCOMPUTEPI(PI,PRESSURE,NX,NY,NZ)
      IMPLICIT NONE
      INTEGER NX,NY,NZ
      DOUBLE PRECISION PI(NX,NY,NZ)
      DOUBLE PRECISION PRESSURE(NX,NY,NZ)
C NCLEND

      INTEGER I,J,K
      DOUBLE PRECISION P1000MB,R_D,CP
      PARAMETER (P1000MB=100000.D0,R_D=287.D0,CP=7.D0*R_D/2.D0)

      DO K = 1,NZ
          DO J = 1,NY
              DO I = 1,NX
                  PI(I,J,K) = (PRESSURE(I,J,K)/P1000MB)** (R_D/CP)
              END DO
          END DO
      END DO

      END

C NCLFORTSTART
      SUBROUTINE DCOMPUTETK(TK,PRESSURE,THETA,NX)
      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION PI
      DOUBLE PRECISION PRESSURE(NX)
      DOUBLE PRECISION THETA(NX)
      DOUBLE PRECISION TK(NX)
C NCLEND

      INTEGER I
      DOUBLE PRECISION P1000MB,R_D,CP
      PARAMETER (P1000MB=100000.D0,R_D=287.D0,CP=7.D0*R_D/2.D0)

      DO I = 1,NX
         PI = (PRESSURE(I)/P1000MB)** (R_D/CP)
         TK(I) = PI*THETA(I)
      END DO

      END

C NCLFORTSTART
      SUBROUTINE DINTERP3DZ(V3D,V2D,Z,LOC,NX,NY,NZ,VMSG)
      IMPLICIT NONE
      INTEGER NX,NY,NZ
      DOUBLE PRECISION V3D(NX,NY,NZ),V2D(NX,NY)
      DOUBLE PRECISION Z(NX,NY,NZ)
      DOUBLE PRECISION LOC
      DOUBLE PRECISION VMSG
C NCLEND

      INTEGER I,J,KP,IP,IM
      LOGICAL INTERP
      DOUBLE PRECISION HEIGHT,W1,W2

      HEIGHT = LOC

c does vertical coordinate increase or decrease with increasing k?
c set offset appropriately

      IP = 0
      IM = 1
      IF (Z(1,1,1).GT.Z(1,1,NZ)) THEN
          IP = 1
          IM = 0
      END IF

      DO I = 1,NX
          DO J = 1,NY
C Initialize to missing.  Was initially hard-coded to -999999.
              V2D(I,J) = VMSG
              INTERP = .false.
              KP = NZ

              DO WHILE ((.NOT.INTERP) .AND. (KP.GE.2))

                  IF (((Z(I,J,KP-IM).LE.HEIGHT).AND. (Z(I,J,
     +                KP-IP).GT.HEIGHT))) THEN
                      W2 = (HEIGHT-Z(I,J,KP-IM))/
     +                     (Z(I,J,KP-IP)-Z(I,J,KP-IM))
                      W1 = 1.D0 - W2
                      V2D(I,J) = W1*V3D(I,J,KP-IM) + W2*V3D(I,J,KP-IP)
                      INTERP = .true.
                  END IF
                  KP = KP - 1

              END DO

          END DO
      END DO

      RETURN
      END

C NCLFORTSTART
      SUBROUTINE DZSTAG(ZNEW,NX,NY,NZ,Z,NXZ,NYZ,NZZ,TERRAIN)
      IMPLICIT NONE
      INTEGER NX,NY,NZ,NXZ,NYZ,NZZ
      DOUBLE PRECISION ZNEW(NX,NY,NZ),Z(NXZ,NYZ,NZZ)
      DOUBLE PRECISION TERRAIN(NXZ,NYZ)
C NCLEND

      INTEGER I,J,K,II,IM1,JJ,JM1

c check for u, v, or w (x,y,or z) staggering
c
c for x and y stag, avg z to x, y, point
c
      IF (NX.GT.NXZ) THEN

          DO K = 1,NZ
              DO J = 1,NY
                  DO I = 1,NX
                      II = MIN0(I,NXZ)
                      IM1 = MAX0(I-1,1)
                      ZNEW(I,J,K) = 0.5D0* (Z(II,J,K)+Z(IM1,J,K))
                  END DO
              END DO
          END DO

      ELSE IF (NY.GT.NYZ) THEN

          DO K = 1,NZ
              DO J = 1,NY
                  JJ = MIN0(J,NYZ)
                  JM1 = MAX0(J-1,1)
                  DO I = 1,NX
                      ZNEW(I,J,K) = 0.5D0* (Z(I,JJ,K)+Z(I,JM1,K))
                  END DO
              END DO
          END DO
c
c w (z) staggering
c
      ELSE IF (NZ.GT.NZZ) THEN

          DO J = 1,NY
              DO I = 1,NX
                  ZNEW(I,J,1) = TERRAIN(I,J)
              END DO
          END DO

          DO K = 2,NZ
              DO J = 1,NY
                  DO I = 1,NX
                      ZNEW(I,J,K) = ZNEW(I,J,K-1) +
     +                              2.D0* (Z(I,J,K-1)-ZNEW(I,J,K-1))
                  END DO
              END DO
          END DO

      END IF

      RETURN
      END

C NCLFORTSTART
      SUBROUTINE DINTERP2DXY(V3D,V2D,XY,NX,NY,NZ,NXY)
      IMPLICIT NONE
      INTEGER NX,NY,NZ,NXY
      DOUBLE PRECISION V3D(NX,NY,NZ),V2D(NXY,NZ)
      DOUBLE PRECISION XY(2,NXY)
C NCLEND

      INTEGER I,J,K,IJ
      DOUBLE PRECISION W11,W12,W21,W22,WX,WY

      DO IJ = 1,NXY

          I = MAX0(1,MIN0(NX-1,INT(XY(1,IJ)+1)))
          J = MAX0(1,MIN0(NY-1,INT(XY(2,IJ)+1)))
          WX = DBLE(I+1) - (XY(1,IJ)+1)
          WY = DBLE(J+1) - (XY(2,IJ)+1)
          W11 = WX*WY
          W21 = (1.D0-WX)*WY
          W12 = WX* (1.D0-WY)
          W22 = (1.D0-WX)* (1.D0-WY)
          DO K = 1,NZ
              V2D(IJ,K) = W11*V3D(I,J,K) + W21*V3D(I+1,J,K) +
     +                    W12*V3D(I,J+1,K) + W22*V3D(I+1,J+1,K)
          END DO
      END DO

      RETURN
      END

C NCLFORTSTART
      SUBROUTINE DINTERP1D(V_IN,V_OUT,Z_IN,Z_OUT,NZ_IN,NZ_OUT,VMSG)
      IMPLICIT NONE
      INTEGER NZ_IN,NZ_OUT
      DOUBLE PRECISION V_IN(NZ_IN),Z_IN(NZ_IN)
      DOUBLE PRECISION V_OUT(NZ_OUT),Z_OUT(NZ_OUT)
      DOUBLE PRECISION VMSG
C NCLEND

      INTEGER KP,K,IM,IP
      LOGICAL INTERP
      DOUBLE PRECISION HEIGHT,W1,W2

c does vertical coordinate increase of decrease with increasing k?
c set offset appropriately

      IP = 0
      IM = 1
      IF (Z_IN(1).GT.Z_IN(NZ_IN)) THEN
          IP = 1
          IM = 0
      END IF

      DO K = 1,NZ_OUT
          V_OUT(K) = VMSG

          INTERP = .false.
          KP = NZ_IN
          HEIGHT = Z_OUT(K)

          DO WHILE ((.NOT.INTERP) .AND. (KP.GE.2))

              IF (((Z_IN(KP-IM).LE.HEIGHT).AND.
     +            (Z_IN(KP-IP).GT.HEIGHT))) THEN
                  W2 = (HEIGHT-Z_IN(KP-IM))/ (Z_IN(KP-IP)-Z_IN(KP-IM))
                  W1 = 1.D0 - W2
                  V_OUT(K) = W1*V_IN(KP-IM) + W2*V_IN(KP-IP)
                  INTERP = .true.
              END IF
              KP = KP - 1

          END DO

      END DO

      RETURN
      END

c---------------------------------------------

c Bill,
c This routine assumes
c    index order is (i,j,k)
c    wrf staggering
C
c    units: pressure (Pa), temperature(K), height (m), mixing ratio
c     (kg kg{-1}) availability of 3d p, t, and qv; 2d terrain; 1d 
c half-level zeta string
c    output units of SLP are Pa, but you should divide that by 100 for the
c          weather weenies.
c    virtual effects are included
c

C NCLFORTSTART
      SUBROUTINE DCOMPUTESEAPRS(NX,NY,NZ,Z,T,P,Q,SEA_LEVEL_PRESSURE,
     +                          T_SEA_LEVEL,T_SURF,LEVEL)
      IMPLICIT NONE
c     Estimate sea level pressure.
      INTEGER NX,NY,NZ
      DOUBLE PRECISION Z(NX,NY,NZ)
      DOUBLE PRECISION T(NX,NY,NZ),P(NX,NY,NZ),Q(NX,NY,NZ)
c     The output is the 2d sea level pressure.
      DOUBLE PRECISION SEA_LEVEL_PRESSURE(NX,NY)
      INTEGER LEVEL(NX,NY)
      DOUBLE PRECISION T_SURF(NX,NY),T_SEA_LEVEL(NX,NY)
C NCLEND
c     Some required physical constants:

      DOUBLE PRECISION R,G,GAMMA
      PARAMETER (R=287.04D0,G=9.81D0,GAMMA=0.0065D0)

c     Specific constants for assumptions made in this routine:

      DOUBLE PRECISION TC,PCONST
      PARAMETER (TC=273.16D0+17.5D0,PCONST=10000)
      LOGICAL RIDICULOUS_MM5_TEST
      PARAMETER (RIDICULOUS_MM5_TEST=.TRUE.)
c      PARAMETER (ridiculous_mm5_test = .false.)

c     Local variables:

      INTEGER I,J,K
      INTEGER KLO,KHI


      DOUBLE PRECISION PLO,PHI,TLO,THI,ZLO,ZHI
      DOUBLE PRECISION P_AT_PCONST,T_AT_PCONST,Z_AT_PCONST
      DOUBLE PRECISION Z_HALF_LOWEST

      LOGICAL L1,L2,L3,FOUND

C
c  Find least zeta level that is PCONST Pa above the surface.  We
c  later use this level to extrapolate a surface pressure and 
c  temperature, which is supposed to reduce the effect of the diurnal
c  heating cycle in the pressure field.

      DO J = 1,NY
          DO I = 1,NX
              LEVEL(I,J) = -1

              K = 1
              FOUND = .false.
              DO WHILE ((.NOT.FOUND) .AND. (K.LE.NZ))
                  IF (P(I,J,K).LT.P(I,J,1)-PCONST) THEN
                      LEVEL(I,J) = K
                      FOUND = .true.
                  END IF
                  K = K + 1
              END DO

              IF (LEVEL(I,J).EQ.-1) THEN
                  PRINT '(A,I4,A)','Troubles finding level ',
     +              NINT(PCONST)/100,' above ground.'
                  PRINT '(A,I4,A,I4,A)','Problems first occur at (',I,
     +              ',',J,')'
                  PRINT '(A,F6.1,A)','Surface pressure = ',P(I,J,1)/100,
     +              ' hPa.'
                  STOP 'Error_in_finding_100_hPa_up'
              END IF


          END DO
      END DO

c     Get temperature PCONST Pa above surface.  Use this to extrapolate
c     the temperature at the surface and down to sea level.

      DO J = 1,NY
          DO I = 1,NX

              KLO = MAX(LEVEL(I,J)-1,1)
              KHI = MIN(KLO+1,NZ-1)

              IF (KLO.EQ.KHI) THEN
                  PRINT '(A)','Trapping levels are weird.'
                  PRINT '(A,I3,A,I3,A)','klo = ',KLO,', khi = ',KHI,
     +              ': and they should not be equal.'
                  STOP 'Error_trapping_levels'
              END IF

              PLO = P(I,J,KLO)
              PHI = P(I,J,KHI)
              TLO = T(I,J,KLO)* (1.D0+0.608D0*Q(I,J,KLO))
              THI = T(I,J,KHI)* (1.D0+0.608D0*Q(I,J,KHI))
c         zlo = zetahalf(klo)/ztop*(ztop-terrain(i,j))+terrain(i,j)
c         zhi = zetahalf(khi)/ztop*(ztop-terrain(i,j))+terrain(i,j)
              ZLO = Z(I,J,KLO)
              ZHI = Z(I,J,KHI)

              P_AT_PCONST = P(I,J,1) - PCONST
              T_AT_PCONST = THI - (THI-TLO)*LOG(P_AT_PCONST/PHI)*
     +                      LOG(PLO/PHI)
              Z_AT_PCONST = ZHI - (ZHI-ZLO)*LOG(P_AT_PCONST/PHI)*
     +                      LOG(PLO/PHI)

              T_SURF(I,J) = T_AT_PCONST* (P(I,J,1)/P_AT_PCONST)**
     +                      (GAMMA*R/G)
              T_SEA_LEVEL(I,J) = T_AT_PCONST + GAMMA*Z_AT_PCONST

          END DO
      END DO

C
c If we follow a traditional computation, there is a correction to the
c sea level temperature if both the surface and sea level 
c temperatures are *too* hot.

      IF (RIDICULOUS_MM5_TEST) THEN
          DO J = 1,NY
              DO I = 1,NX
                  L1 = T_SEA_LEVEL(I,J) .LT. TC
                  L2 = T_SURF(I,J) .LE. TC
                  L3 = .NOT. L1
                  IF (L2 .AND. L3) THEN
                      T_SEA_LEVEL(I,J) = TC
                  ELSE
                      T_SEA_LEVEL(I,J) = TC -
     +                                   0.005D0* (T_SURF(I,J)-TC)**2
                  END IF
              END DO
          END DO
      END IF

c     The grand finale: ta da!

      DO J = 1,NY
          DO I = 1,NX
c   z_half_lowest=zetahalf(1)/ztop*(ztop-terrain(i,j))+terrain(i,j)
              Z_HALF_LOWEST = Z(I,J,1)

C Convert to hPa in this step, by multiplying by 0.01. The original
C Fortran routine didn't do this, but the NCL script that called it
C did, so we moved it here.
              SEA_LEVEL_PRESSURE(I,J) = 0.01 * (P(I,J,1)*
     +                                  EXP((2.D0*G*Z_HALF_LOWEST)/
     +                                  (R* (T_SEA_LEVEL(I,J)+T_SURF(I,
     +                                  J)))))
          END DO
      END DO

c     print *,'sea pres input at weird location i=20,j=1,k=1'
c     print *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
c     print *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
c     print *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
c     print *,'slp=',sea_level_pressure(20,1),
c    *         sea_level_pressure(20,2),sea_level_pressure(20,3)

      END


c---------------------------------------------------

C
C Double precision version. If you make a change here, you
C must make the same change below to filter2d.
C
C NCLFORTSTART
      SUBROUTINE DFILTER2D(A,B,NX,NY,IT,MISSING)
      IMPLICIT NONE
c     Estimate sea level pressure.
      INTEGER NX,NY,IT
      DOUBLE PRECISION A(NX,NY),B(NX,NY),MISSING
C NCLEND

      DOUBLE PRECISION COEF
      PARAMETER (COEF=0.25D0)
      INTEGER I,J,ITER

      DO ITER = 1,IT
          DO J = 1,NY
              DO I = 1,NX
                  B(I,J) = A(I,J)
              END DO
          END DO
          DO J = 2,NY - 1
             DO I = 1,NX
                IF ( B(I,J-1).EQ.MISSING .OR. B(I,J).EQ.MISSING .OR.
     +               B(I,J+1).EQ.MISSING ) THEN
                   A(I,J) = A(I,J)
                ELSE
                   A(I,J) = A(I,J) + COEF* (B(I,J-1)-2*B(I,J)+B(I,J+1))
                END IF
             END DO
          END DO
          DO J = 1,NY
             DO I = 2,NX - 1
                IF ( B(I-1,J).EQ.MISSING .OR. B(I,J).EQ.MISSING .OR.
     +               B(I+1,J).EQ.MISSING ) THEN
                   A(I,J) = A(I,J)
                ELSE
                   A(I,J) = A(I,J) + COEF* (B(I-1,J)-2*B(I,J)+B(I+1,J))
                END IF
             END DO
          END DO
c        do j=1,ny
c        do i=1,nx
c          b(i,j) = a(i,j)
c        enddo
c        enddo
c        do j=2,ny-1
c        do i=1,nx
c          a(i,j) = a(i,j) - .99*coef*(b(i,j-1)-2*b(i,j)+b(i,j+1))
c        enddo
c        enddo
c        do j=1,ny
c        do i=2,nx-1
c          a(i,j) = a(i,j) - .99*coef*(b(i-1,j)-2*b(i,j)+b(i+1,j))
c        enddo
c        enddo
      END DO
      RETURN
      END

C
C Single precision version. If you make a change here, you
C must make the same change above to dfilter2d.
C
C NCLFORTSTART
      SUBROUTINE filter2d( a, b, nx , ny , it, missing)
      IMPLICIT NONE
c     Estimate sea level pressure.
      INTEGER nx , ny, it
      REAL    a(nx,ny),b(nx,ny), missing
C NCLEND

      REAL coef
      parameter( coef = 0.25)
      INTEGER i,j,iter

      do iter=1, it
        do j=1,ny
        do i=1,nx
          b(i,j) = a(i,j)
        enddo
        enddo
        do j=2,ny-1
        do i=1,nx
          if ( b(i,j-1).eq.missing .or. b(i,j).eq.missing .or.
     +         b(i,j+1).eq.missing ) then
             a(i,j) = a(i,j)
          else
             a(i,j) = a(i,j) + coef*(b(i,j-1)-2*b(i,j)+b(i,j+1))
          end if
        enddo
        enddo
        do j=1,ny
        do i=2,nx-1
           if ( b(i-1,j).eq.missing .or. b(i,j).eq.missing .or.
     +          b(i+1,j).eq.missing ) then
              a(i,j) = a(i,j)
           else
              a(i,j) = a(i,j) + coef*(b(i-1,j)-2*b(i,j)+b(i+1,j))
           end if
        enddo
        enddo
c        do j=1,ny
c        do i=1,nx
c          b(i,j) = a(i,j)
c        enddo
c        enddo
c        do j=2,ny-1
c        do i=1,nx
c          a(i,j) = a(i,j) - .99*coef*(b(i,j-1)-2*b(i,j)+b(i,j+1))
c        enddo
c        enddo
c        do j=1,ny
c        do i=2,nx-1
c          a(i,j) = a(i,j) - .99*coef*(b(i-1,j)-2*b(i,j)+b(i+1,j))
c        enddo
c        enddo
      enddo
      return
      end
c---------------------------------------------------------

C NCLFORTSTART
      SUBROUTINE DCOMPUTERH(QV,P,T,RH,NX)

      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION QV(NX),P(NX),T(NX),RH(NX)
C NCLEND
      DOUBLE PRECISION SVP1,SVP2,SVP3,SVPT0
      PARAMETER (SVP1=0.6112D0,SVP2=17.67D0,SVP3=29.65D0,SVPT0=273.15D0)
      INTEGER I
      DOUBLE PRECISION QVS,ES,PRESSURE,TEMPERATURE
      DOUBLE PRECISION EP_2,R_D,R_V
      PARAMETER (R_D=287.D0,R_V=461.6D0,EP_2=R_D/R_V)
      DOUBLE PRECISION EP_3
      PARAMETER (EP_3=0.622D0)

      DO I = 1,NX
         PRESSURE = P(I)
         TEMPERATURE = T(I)
c       es  = 1000.*svp1*
         ES = 10.D0*SVP1*EXP(SVP2* (TEMPERATURE-SVPT0)/
     +        (TEMPERATURE-SVP3))
c       qvs = ep_2*es/(pressure-es)
         QVS = EP_3*ES/ (0.01D0*PRESSURE- (1.D0-EP_3)*ES)
c        rh = 100*amax1(1., qv(i)/qvs)
c       rh(i) = 100.*qv(i)/qvs
         RH(I) = 100.D0*DMAX1(DMIN1(QV(I)/QVS,1.0D0),0.0D0)
      END DO

      RETURN
      END

c----------------------------------------------

C NCLFORTSTART
      SUBROUTINE DGETIJLATLONG(LAT_ARRAY,LONG_ARRAY,LAT,LONGITUDE,
     +                         II,JJ,NX,NY,IMSG)
      IMPLICIT NONE
      INTEGER NX,NY,II,JJ,IMSG
      DOUBLE PRECISION LAT_ARRAY(NX,NY),LONG_ARRAY(NX,NY)
      DOUBLE PRECISION LAT,LONGITUDE
C NCLEND
      DOUBLE PRECISION LONGD,LATD
      INTEGER I,J
      DOUBLE PRECISION IR,JR
      DOUBLE PRECISION DIST_MIN,DIST

C Init to missing. Was hard-coded to -999 initially.
      IR = IMSG
      JR = IMSG

      DIST_MIN = 1.D+20
      DO J = 1,NY
          DO I = 1,NX
              LATD = (LAT_ARRAY(I,J)-LAT)**2
              LONGD = (LONG_ARRAY(I,J)-LONGITUDE)**2
C             LONGD = DMIN1((LONG_ARRAY(I,J)-LONGITUDE)**2,
C    +                (LONG_ARRAY(I,J)+LONGITUDE)**2)
              DIST = SQRT(LATD+LONGD)
              IF (DIST_MIN.GT.DIST) THEN
                  DIST_MIN = DIST
                  IR = DBLE(I)
                  JR = DBLE(J)
              END IF
          END DO
      END DO
C
C The original version of this routine returned IR and JR. But, then
C the NCL script that called this routine was converting IR and JR
C to integer, so why not just return II and JJ?
C
C Also, I'm subtracing 1 here, because it will be returned to NCL
C script which has 0-based indexing.
C 
      IF(IR.ne.IMSG.and.JR.ne.IMSG) then
        II = NINT(IR)-1
        JJ = NINT(JR)-1
      ELSE
        II = IMSG
        JJ = IMSG
      END IF

c we will just return the nearest point at present

      RETURN
      END

C NCLFORTSTART
      SUBROUTINE DCOMPUTEUVMET(U,V,UVMET,LONGCA,LONGCB,FLONG,FLAT,
     +                         CEN_LONG,CONE,RPD,NX,NY,NXP1,NYP1,
     +                         ISTAG,IS_MSG_VAL,UMSG,VMSG,UVMETMSG)
      IMPLICIT NONE

C ISTAG should be 0 if the U,V grids are not staggered.
C That is, NY = NYP1 and NX = NXP1.

      INTEGER NX,NY,NXP1,NYP1,ISTAG
      LOGICAL IS_MSG_VAL
      DOUBLE PRECISION U(NXP1,NY),V(NX,NYP1)
      DOUBLE PRECISION UVMET(NX,NY,2)
      DOUBLE PRECISION FLONG(NX,NY),FLAT(NX,NY)
      DOUBLE PRECISION LONGCB(NX,NY),LONGCA(NX,NY)
      DOUBLE PRECISION CEN_LONG,CONE,RPD
      DOUBLE PRECISION UMSG,VMSG,UVMETMSG
C NCLEND

      INTEGER I,J
      DOUBLE PRECISION UK,VK


c      WRITE (6,FMT=*) ' in compute_uvmet ',NX,NY,NXP1,NYP1,ISTAG

      DO J = 1,NY
          DO I = 1,NX

              LONGCA(I,J) = FLONG(I,J) - CEN_LONG
              IF (LONGCA(I,J).GT.180.D0) THEN
                  LONGCA(I,J) = LONGCA(I,J) - 360.D0
              END IF
              IF (LONGCA(I,J).LT.-180.D0) THEN
                  LONGCA(I,J) = LONGCA(I,J) + 360.D0
              END IF
              IF (FLAT(I,J).LT.0.D0) THEN
                  LONGCB(I,J) = -LONGCA(I,J)*CONE*RPD
              ELSE
                  LONGCB(I,J) = LONGCA(I,J)*CONE*RPD
              END IF

              LONGCA(I,J) = COS(LONGCB(I,J))
              LONGCB(I,J) = SIN(LONGCB(I,J))

          END DO
      END DO

c      WRITE (6,FMT=*) ' computing velocities '

      DO J = 1,NY
         DO I = 1,NX
            IF (ISTAG.EQ.1) THEN
               IF (IS_MSG_VAL.AND.(U(I,J).EQ.UMSG.OR.
     +                             V(I,J).EQ.VMSG.OR.
     +                             U(I+1,J).EQ.UMSG.OR.
     +                             V(I,J+1).EQ.VMSG)) THEN
                  UVMET(I,J,1) = UVMETMSG
                  UVMET(I,J,2) = UVMETMSG
               ELSE
                  UK = 0.5D0* (U(I,J)+U(I+1,J))
                  VK = 0.5D0* (V(I,J)+V(I,J+1))
                  UVMET(I,J,1) = VK*LONGCB(I,J) + UK*LONGCA(I,J)
                  UVMET(I,J,2) = VK*LONGCA(I,J) - UK*LONGCB(I,J)
               END IF
            ELSE
               IF (IS_MSG_VAL.AND.(U(I,J).EQ.UMSG.OR.
     +                             V(I,J).EQ.VMSG)) THEN
                  UVMET(I,J,1) = UVMETMSG
                  UVMET(I,J,2) = UVMETMSG
               ELSE
                  UK = U(I,J)
                  VK = V(I,J)
                  UVMET(I,J,1) = VK*LONGCB(I,J) + UK*LONGCA(I,J)
                  UVMET(I,J,2) = VK*LONGCA(I,J) - UK*LONGCB(I,J)
               END IF
            END IF
         END DO
      END DO

      RETURN
      END

C NCLFORTSTART
C
C This was originally a routine that took 2D input arrays. Since
C the NCL C wrapper routine can handle multiple dimensions, it's
C not necessary to have anything bigger than 1D here.
C
      SUBROUTINE DCOMPUTETD(TD,PRESSURE,QV_IN,NX)
      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION PRESSURE(NX)
      DOUBLE PRECISION QV_IN(NX)
      DOUBLE PRECISION TD(NX)
C NCLEND
      DOUBLE PRECISION QV,TDC

      INTEGER I

      DO I = 1,NX
          QV = DMAX1(QV_IN(I),0.D0)
c vapor pressure
          TDC = QV*PRESSURE(I)/ (.622D0+QV)

c avoid problems near zero
          TDC = DMAX1(TDC,0.001D0)
          TD(I) = (243.5D0*LOG(TDC)-440.8D0)/ (19.48D0-LOG(TDC))
      END DO

      RETURN
      END

C NCLFORTSTART
      SUBROUTINE DCOMPUTEICLW(ICLW,PRESSURE,QC_IN,NX,NY,NZ)
      IMPLICIT NONE
      INTEGER NX,NY,NZ
      DOUBLE PRECISION PRESSURE(NX,NY,NZ)
      DOUBLE PRECISION QC_IN(NX,NY,NZ)
      DOUBLE PRECISION ICLW(NX,NY)
      DOUBLE PRECISION QCLW,DP,GG
C NCLEND

      INTEGER I,J,K

      GG = 1000.D0/9.8D0

      DO J = 1,NY
          DO I = 1,NX
              ICLW(I,J) = 0.D0
          END DO
      END DO

      DO J = 3,NY - 2
          DO I = 3,NX - 2
              DO K = 1,NZ
                  QCLW = DMAX1(QC_IN(I,J,K),0.D0)
                  IF (K.EQ.1) THEN
                      DP = (PRESSURE(I,J,K-1)-PRESSURE(I,J,K))
                  ELSE IF (K.EQ.NZ) THEN
                      DP = (PRESSURE(I,J,K)-PRESSURE(I,J,K+1))
                  ELSE
                      DP = (PRESSURE(I,J,K-1)-PRESSURE(I,J,K+1))/2.D0
                  END IF
                  ICLW(I,J) = ICLW(I,J) + QCLW*DP*GG
              END DO
          END DO
      END DO

      RETURN
      END
