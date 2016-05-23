C
C premaptform.f and maptform.f copied from RIP/src
C By So-Young Ha on Sep 29, 2005.
C
C
C NCLFORTSTART
      SUBROUTINE DMAPTFORM(DSKMC,MIYCORS,MJXCORS,NPROJ,XLATC,XLONC,
     +                     TRUE1,TRUE2,RIY,RJX,RLAT,RLON,IDIR)
C
C Input vars:        DSKMC, MIYCORS, MJXCORS, NPROJ, XLATC, XLONC,
C                    NPROJ, IDIR
C Input/output vars: RIY, RIX, RLAT
C Output vars:       TRUE1, TRUE2, RLON 
C 
C
C Possible NCL interface:
C
C   wrf_maptform(dskmc, miycors, mjxcors, nproj, xlatc, xlonc, riy, rjx,
C                idir, rlat, rlon, opts)
C
C where opts could contain the TRUE1 and TRUE2 information in some fashion.
C
      DOUBLE PRECISION PI_MPTF
      DOUBLE PRECISION RPD_MPTF
      DOUBLE PRECISION REARTH_MPTF
      DOUBLE PRECISION DSKMC_MPTF
      DOUBLE PRECISION XLONC_MPTF
      DOUBLE PRECISION CIY_MPTF
      DOUBLE PRECISION CJX_MPTF
      DOUBLE PRECISION CONE_MPTF
      DOUBLE PRECISION CONEI_MPTF
      DOUBLE PRECISION C1_MPTF
      DOUBLE PRECISION C2_MPTF
      DOUBLE PRECISION YC_MPTF
      DOUBLE PRECISION COTRUE1
      DOUBLE PRECISION YPOINT
      DOUBLE PRECISION XPOINT
      DOUBLE PRECISION DLON
C
c   This routine converts a coarse domain dot grid point, <riy,rjx>,
c   into a lat/lon point <rlat,rlon> if idir=1, or vice versa if
c   idir=-1. It works for Lambert Conformal (LC,1),
c   Polar Stereographic (ST,2), or Mercator (ME,3) projections,
c   with any true latitide(s).
c   It is assumed that premaptform has been called prior to this so
c   that the proper constants have been placed in the common block
c   called mptf, which should be declared in (and only in) the
c   main program and routines maptform (this routine) and premaptform.
c

C Input, Output Args
      INTEGER MIYCORS,MJXCORS,NPROJ
      DOUBLE PRECISION DSKMC,XLATC,XLONC,TRUE1,TRUE2
      INTEGER IDIR
C Latitude (-90->90 deg N)
      DOUBLE PRECISION RLAT
C Longitude (-180->180 E)
      DOUBLE PRECISION RLON
C Cartesian X coordinate
      DOUBLE PRECISION RIY
C Cartesian Y coordinate
      DOUBLE PRECISION RJX
C NCLEND


c ===========
c premaptform
c ===========
C 3.1415...
      PI_MPTF = 4.D0*ATAN(1.D0)
C radians per degree
      RPD_MPTF = PI_MPTF/180.D0
C radius of planet, in km
      REARTH_MPTF = 6370.949D0
      DSKMC_MPTF = DSKMC
      XLONC_MPTF = XLONC
      NPROJ_MPTF = NPROJ
      CIY_MPTF = .5D0* (1.D0+MIYCORS)
      CJX_MPTF = .5D0* (1.D0+MJXCORS)
c
C Mercator
      IF (NPROJ_MPTF.EQ.3) THEN
c
          TRUE1 = 0.D0
          TRUE2 = 0.D0
          IHM_MPTF = 1
          CONE_MPTF = 1.D0
          CONEI_MPTF = 1.D0
          C1_MPTF = 1.D0
          C2_MPTF = 1.D0
          YC_MPTF = REARTH_MPTF*LOG((1.D0+SIN(RPD_MPTF*XLATC))/
     +              COS(RPD_MPTF*XLATC))
c
C Lambert Comformal or Polar Stereographic
      ELSE
c
c   Make sure xlatc, true1, and true2 are all in same hemisphere,
c      and calculate ihm_mptf.
c
          IF (XLATC.GT.0.D0 .AND. TRUE1.GT.0.D0 .AND.
     +        TRUE2.GT.0.D0) THEN
              IHM_MPTF = 1
          ELSE IF (XLATC.LT.0.D0 .AND. TRUE1.LT.0.D0 .AND.
     +             TRUE2.LT.0.D0) THEN
              IHM_MPTF = -1
          ELSE
              WRITE (*,FMT=*) 'Invalid latitude parameters for map.'
              STOP
          END IF
c
c   Calculate cone factor
c
          IF (NPROJ_MPTF.EQ.1) THEN
              IF (TRUE1.NE.TRUE2) THEN
                  CONE_MPTF = LOG10(COS(RPD_MPTF*TRUE1)/
     +                        COS(RPD_MPTF*TRUE2))/
     +                        LOG10(TAN(.25D0*PI_MPTF-
     +                        IHM_MPTF*.5D0*RPD_MPTF*TRUE1)/
     +                        TAN(.25D0*PI_MPTF-IHM_MPTF*.5D0*RPD_MPTF*
     +                        TRUE2))
              ELSE
                  CONE_MPTF = COS(RPD_MPTF* (90.D0-IHM_MPTF*TRUE1))
              END IF
          ELSE IF (NPROJ_MPTF.EQ.2) THEN
              CONE_MPTF = 1.D0
          END IF
c
c   Calculate other constants
c
          CONEI_MPTF = 1.D0/CONE_MPTF
          COTRUE1 = IHM_MPTF*90.D0 - TRUE1
          IF (NPROJ_MPTF.EQ.1) THEN
              C1_MPTF = REARTH_MPTF*SIN(RPD_MPTF*COTRUE1)/
     +                  (CONE_MPTF* (IHM_MPTF*TAN(.5D0*RPD_MPTF*
     +                  COTRUE1))**CONE_MPTF)
              C2_MPTF = TAN(.5D0*RPD_MPTF*COTRUE1)*
     +                  (CONE_MPTF/ (IHM_MPTF*REARTH_MPTF*SIN(RPD_MPTF*
     +                  COTRUE1)))**CONEI_MPTF
              YC_MPTF = -C1_MPTF* (IHM_MPTF*
     +                  TAN(.25D0* (IHM_MPTF*PI_MPTF-
     +                  2.D0*RPD_MPTF*XLATC)))**CONE_MPTF
          ELSE IF (NPROJ_MPTF.EQ.2) THEN
              C1_MPTF = 1.D0 + COS(RPD_MPTF*COTRUE1)
              C2_MPTF = 1.D0
              YC_MPTF = -REARTH_MPTF*SIN(.5D0*IHM_MPTF*PI_MPTF-
     +                  RPD_MPTF*XLATC)*C1_MPTF/
     +                  (1.D0+COS(.5D0*IHM_MPTF*PI_MPTF-RPD_MPTF*XLATC))
          END IF
c
      END IF

c ========
c maptform
c ========

      IF (RLAT.EQ.-90.D0) PRINT *,'maptform:',RIY,RJX,RLAT,RLON,IDIR

C First, deal with idir=1
      IF (IDIR.EQ.1) THEN
c
          YPOINT = (RIY-CIY_MPTF)*DSKMC_MPTF + YC_MPTF
          XPOINT = (RJX-CJX_MPTF)*DSKMC_MPTF
c
          IF (NPROJ_MPTF.EQ.3) THEN
              RLAT = (2.D0*ATAN(EXP(YPOINT/REARTH_MPTF))-.5D0*PI_MPTF)/
     +               RPD_MPTF
              RLON = XLONC_MPTF + (XPOINT/REARTH_MPTF)/RPD_MPTF
          ELSE IF (NPROJ_MPTF.EQ.1) THEN
              RLAT = (.5D0*IHM_MPTF*PI_MPTF-
     +               2.D0*ATAN(C2_MPTF* (SQRT(XPOINT**2+
     +               YPOINT**2))**CONEI_MPTF))/RPD_MPTF
              RLON = XLONC_MPTF + (CONEI_MPTF*
     +               ATAN2(XPOINT,-IHM_MPTF*YPOINT))/RPD_MPTF
          ELSE IF (NPROJ_MPTF.EQ.2) THEN
              RLAT = (.5D0*IHM_MPTF*PI_MPTF-
     +               IHM_MPTF*2.D0*ATAN(SQRT(XPOINT**2+
     +               YPOINT**2)/ (REARTH_MPTF*C1_MPTF)))/RPD_MPTF
              IF (XPOINT.EQ.0.D0 .AND. YPOINT.EQ.0.D0) THEN
                  RLON = XLONC_MPTF
              ELSE
                  RLON = XLONC_MPTF + (ATAN2(XPOINT,-IHM_MPTF*YPOINT))/
     +                   RPD_MPTF
              END IF
          END IF
          RLON = MOD(RLON+900.D0,360.D0) - 180.D0
c
C Otherwise, deal with idir=-1
      ELSE
c
          DLON = RLON - XLONC_MPTF
          IF (DLON.LT.-180.D0) DLON = DLON + 360
          IF (DLON.GT.180.D0) DLON = DLON - 360
          IF (NPROJ_MPTF.EQ.3) THEN
              YPOINT = REARTH_MPTF*LOG((1.D0+SIN(RPD_MPTF*RLAT))/
     +                 COS(RPD_MPTF*RLAT))
              XPOINT = DLON*RPD_MPTF*REARTH_MPTF
          ELSE IF (NPROJ_MPTF.EQ.1) THEN
              YPOINT = -C1_MPTF* (IHM_MPTF*
     +                 TAN(.25D0* (IHM_MPTF*PI_MPTF-2.D0*RPD_MPTF*
     +                 RLAT)))**CONE_MPTF*COS(CONE_MPTF*RPD_MPTF*DLON)
              XPOINT = IHM_MPTF*C1_MPTF* (IHM_MPTF*
     +                 TAN(.25D0* (IHM_MPTF*PI_MPTF-
     +                 2.D0*RPD_MPTF*RLAT)))**CONE_MPTF*
     +                 SIN(CONE_MPTF*RPD_MPTF*DLON)
          ELSE IF (NPROJ_MPTF.EQ.2) THEN
              YPOINT = -REARTH_MPTF*SIN(.5D0*IHM_MPTF*PI_MPTF-
     +                 RPD_MPTF*RLAT)*C1_MPTF/ (1.D0+
     +                 COS(.5D0*IHM_MPTF*PI_MPTF-RPD_MPTF*RLAT))*
     +                 COS(RPD_MPTF*DLON)
              XPOINT = IHM_MPTF*REARTH_MPTF*
     +                 SIN(.5D0*IHM_MPTF*PI_MPTF-RPD_MPTF*RLAT)*C1_MPTF/
     +                 (1.D0+COS(.5D0*IHM_MPTF*PI_MPTF-RPD_MPTF*RLAT))*
     +                 SIN(RPD_MPTF*DLON)
          END IF
          RIY = (YPOINT-YC_MPTF)/DSKMC_MPTF + CIY_MPTF
          RJX = XPOINT/DSKMC_MPTF + CJX_MPTF
c
      END IF

      RETURN
      END

C********************************************************
C NCLFORTSTART
      SUBROUTINE DBINT3D(DATA_OUT,OBSII,OBSJJ,DATA_IN,NX,NY,NZ,NOBSICRS,
     +                   NOBSJCRS,ICRS,JCRS)
C
C Possible NCL interface:
C
C    data_out = wrf_bint3d(data_in,obsii,obsjj,icrs,jcrs)
C
C     !!! 1_based_array (cols x rows) in fortran <=> 0_based_array
C      (rows x cols) in NCL !!!
C     !!! Include K-index to make a 3-D array !!!
C
C     INPUT VARIABLES
C     ---------------
      INTEGER ICRS,JCRS,NX,NY,NZ
      INTEGER NOBSJCRS,NOBSICRS
      DOUBLE PRECISION OBSII(NOBSICRS,NOBSJCRS)
      DOUBLE PRECISION OBSJJ(NOBSICRS,NOBSJCRS)
      DOUBLE PRECISION DATA_IN(NX,NY,NZ)

C     OUTPUT
C     ---------------
      DOUBLE PRECISION DATA_OUT(NOBSICRS,NOBSJCRS,NZ)
C NCLEND

C     LOCAL
      DOUBLE PRECISION OBSI,OBSJ
      DOUBLE PRECISION DATA_OBS
C

      DO K = 1,NZ
          DO J = 1,NOBSJCRS
              DO I = 1,NOBSICRS
C grid index in lon
                  OBSI = OBSII(I,J)
C grid index in lat
                  OBSJ = OBSJJ(I,J)
                  DATA_OBS = 0.0D0
                  CALL DBINT(DATA_OBS,OBSI,OBSJ,DATA_IN(1,1,K),NX,NY,
     +                       ICRS,JCRS)
                  DATA_OUT(I,J,K) = DATA_OBS
              END DO
          END DO
      END DO

      RETURN
      END


      SUBROUTINE DBINT(PP,XX,YY,LIST,III,JJJ,ICRS,JCRS)
      DOUBLE PRECISION PP
      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION C
      DOUBLE PRECISION D
      DOUBLE PRECISION E
      DOUBLE PRECISION F
      DOUBLE PRECISION G
      DOUBLE PRECISION H
      DOUBLE PRECISION QQ
C
C --- BI-LINEAR INTERPOLATION AMONG FOUR GRID VALUES
C
C     INPUT : LIST, XX, YY
C     OUTPUT: PP
C
      INTEGER ICRS,JCRS,III,JJJ
      DOUBLE PRECISION XX,YY
      DOUBLE PRECISION LIST(III,JJJ),STL(4,4)

C MASS GRID IN WRF (I-> west-east, J-> south-north)
C
      IB = III - ICRS
      JB = JJJ - JCRS
      PP = 0.0D0
      N = 0
      I = INT(XX+0.00001D0)
      J = INT(YY+0.00001D0)
      X = XX - I
      Y = YY - J
      IF ((ABS(X).GT.0.00001D0) .OR. (ABS(Y).GT.0.00001D0)) THEN
C
          DO 2 K = 1,4
              KK = I + K
              DO 2 L = 1,4
                  STL(K,L) = 0.D0
                  LL = J + L
                  IF ((KK.GE.1) .AND. (KK.LE.IB) .AND. (LL.LE.JB) .AND.
     +                (LL.GE.1)) THEN
                      STL(K,L) = LIST(KK,LL)
                      N = N + 1
C .. a zero value inside the domain being set to 1.E-12:
                      IF (STL(K,L).EQ.0.D0) STL(K,L) = 1.D-12
                  END IF
    2     CONTINUE
C
          CALL DONED(A,X,STL(1,1),STL(2,1),STL(3,1),STL(4,1))
          CALL DONED(B,X,STL(1,2),STL(2,2),STL(3,2),STL(4,2))
          CALL DONED(C,X,STL(1,3),STL(2,3),STL(3,3),STL(4,3))
          CALL DONED(D,X,STL(1,4),STL(2,4),STL(3,4),STL(4,4))
C
C .. CHECK TANGENT LINEAR OF ONED, SAVE BASIC STATE:
C      WRITE(20) XX,YY,Y,A,B,C,D
C
          CALL DONED(PP,Y,A,B,C,D)
          IF (N.NE.16) THEN
              CALL DONED(E,Y,STL(1,1),STL(1,2),STL(1,3),STL(1,4))
              CALL DONED(F,Y,STL(2,1),STL(2,2),STL(2,3),STL(2,4))
              CALL DONED(G,Y,STL(3,1),STL(3,2),STL(3,3),STL(3,4))
              CALL DONED(H,Y,STL(4,1),STL(4,2),STL(4,3),STL(4,4))
C .. CHECK TANGENT LINEAR OF ONED, SAVE BASIC STATE:
C      WRITE(20) XX,YY,X,E,F,G,H
C
              CALL DONED(QQ,X,E,F,G,H)
              PP = (PP+QQ)*0.5D0
          END IF
C
      ELSE
C
          PP = LIST(I,J)
      END IF
C
      RETURN
      END



      SUBROUTINE DONED(Y,X,A,B,C,D)
      DOUBLE PRECISION Y
      DOUBLE PRECISION X
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION C
      DOUBLE PRECISION D
      DOUBLE PRECISION ONE
C
C ..  Input : X, A, B, C, D
C     Output: Y
C       1, 2, 3, and 4 points interpolation:
C       In this subroutine, the zero value of A, B, C, D means that
C       point outside the domain.
C
C .. 1-point:
C .. take the value at the second point:
      IF (X.EQ.0.D0) THEN
          ONE = B
C .. take the value at the third point:
      ELSE IF (X.EQ.1.D0) THEN
          ONE = C
C .. the point X outside the range:
      ELSE IF (B*C.EQ.0.D0) THEN
          ONE = 0.D0
      ELSE
          IF (A*D.EQ.0.D0) THEN
C .. 3-point interpolation:
              IF (A.NE.0.D0) THEN
                  ONE = B + X* (0.5D0* (C-A)+X* (0.5D0* (C+A)-B))
              ELSE IF (D.NE.0.D0) THEN
                  ONE = C + (1.0D0-X)* (0.5D0* (B-D)+
     +                  (1.0D0-X)* (0.5D0* (B+D)-C))
              ELSE
C .. 2-point interpolation:
                  ONE = B* (1.0D0-X) + C*X
              END IF
          ELSE
C .. 4-point interpolation:
              ONE = (1.0D0-X)* (B+X* (0.5D0* (C-A)+X* (0.5D0* (C+A)-B)))
     +              + X* (C+ (1.0D0-X)* (0.5D0* (B-D)+ (1.0D0-
     +              X)* (0.5D0* (B+D)-C)))
          END IF
      END IF
C
      Y = ONE
C
      RETURN

      END
