c======================================================================
c
c !IROUTINE: capecalc3d -- Calculate CAPE and CIN
c
c !DESCRIPTION:
c
c   If i3dflag=1, this routine calculates CAPE and CIN (in m**2/s**2,
c   or J/kg) for every grid point in the entire 3D domain (treating
c   each grid point as a parcel).  If i3dflag=0, then it
c   calculates CAPE and CIN only for the parcel with max theta-e in
c   the column, (i.e. something akin to Colman's MCAPE).  By "parcel",
c   we mean a 500-m deep parcel, with actual temperature and moisture
c   averaged over that depth.
c
c   In the case of i3dflag=0,
c   CAPE and CIN are 2D fields that are placed in the k=mkzh slabs of
c   the cape and cin arrays.  Also, if i3dflag=0, LCL and LFC heights
c   are put in the k=mkzh-1 and k=mkzh-2 slabs of the cin array.
c
c ASSUMPTIONS:
c
c !REVISION HISTORY:
c     2005-May-15 - Mark T. Stoelinga - oringinal version from RIP4
c     2005-Nov-28 - J. Schramm - modified to run outside of RIP4 with
c     2012-Jul-18 - M. Haley - modified to change/add missing value.
c                                NCL
c
c !INTERFACE:
c ------------------------------------------------------------------
C NCLFORTSTART
      SUBROUTINE DCAPECALC3D(PRS,TMK,QVP,GHT,TER,SFP,CAPE,CIN,CMSG,
     +                       MIY,MJX,MKZH,I3DFLAG,TER_FOLLOW,PSAFILE)
c
      IMPLICIT NONE
      INTEGER MIY,MJX,MKZH,I3DFLAG,TER_FOLLOW
      DOUBLE PRECISION PRS(MIY,MJX,MKZH)
      DOUBLE PRECISION TMK(MIY,MJX,MKZH)
      DOUBLE PRECISION QVP(MIY,MJX,MKZH)
      DOUBLE PRECISION GHT(MIY,MJX,MKZH)
      DOUBLE PRECISION TER(MIY,MJX)
      DOUBLE PRECISION SFP(MIY,MJX)
      DOUBLE PRECISION CAPE(MIY,MJX,MKZH)
      DOUBLE PRECISION CIN(MIY,MJX,MKZH)
      DOUBLE PRECISION CMSG
      CHARACTER*(*) PSAFILE
C NCLEND
c Local variables
      INTEGER I,J,K,ILCL,IUP,KEL,KK,KLCL,KLEV,KLFC,KMAX,KPAR,KPAR1,KPAR2
      DOUBLE PRECISION DAVG,ETHMAX,Q,T,P,E,ETH,TLCL,ZLCL
      DOUBLE PRECISION CP,EPS,GAMMA,GAMMAMD,RGAS,RGASMD,TLCLC1,TLCLC2,
     +                 TLCLC3,TLCLC4
      DOUBLE PRECISION CPMD,THTECON1,THTECON2,THTECON3
      DOUBLE PRECISION CELKEL,EZERO,ESLCON1,ESLCON2,GRAV
      DOUBLE PRECISION PAVG,VIRTUAL,P1,P2,PP1,PP2,TH,TOTTHE,TOTQVP,
     +                 TOTPRS
      DOUBLE PRECISION CPM,DELTAP,ETHPARI,GAMMAM,GHTPARI,QVPPARI,
     +                 PRSPARI,TMKPARI
      DOUBLE PRECISION FACDEN,FAC1,FAC2,QVPLIFT,TMKLIFT,TVENV,TVLIFT,
     +                 GHTLIFT
      DOUBLE PRECISION ESLIFT,TMKENV,QVPENV,TONPSADIABAT
      DOUBLE PRECISION BENAMIN,DZ,PUP,PDN
      DOUBLE PRECISION BUOY(150),ZREL(150),BENACCUM(150),
     +                 PRSF(MIY,MJX,MKZH)
      DOUBLE PRECISION PSADITHTE(150),PSADIPRS(150),PSADITMK(150,150)
c
C The comments were taken from a Mark Stoelinga email, 23 Apr 2007,
C in response to a user getting the "Outside of lookup table bounds"
C error message. 
C
C TMKPARI  - Initial temperature of parcel, K
C    Values of 300 okay. (Not sure how much from this you can stray.)
C
C PRSPARI - Initial pressure of parcel, hPa
C    Values of 980 okay. (Not sure how much from this you can stray.)
C
C THTECON1, THTECON2, THTECON3
C     These are all constants, the first in K and the other two have
C     no units.  Values of 3376, 2.54, and 0.81 were stated as being
C     okay.
C
C TLCL - The temperature at the parcel's lifted condensation level, K
C        should be a reasonable atmospheric temperature around 250-300 K
C        (398 is "way too high")
C
C QVPPARI - The initial water vapor mixing ratio of the parcel,
C           kg/kg (should range from 0.000 to 0.025)
C

c Constants
      IUP = 6
      CELKEL = 273.15D0
      GRAV = 9.81D0
C hPa
      EZERO = 6.112D0
      ESLCON1 = 17.67D0
      ESLCON2 = 29.65D0
      EPS = 0.622D0
C J/K/kg
      RGAS = 287.04D0
C  J/K/kg  Note: not using Bolton's value of 1005.7
      CP = 1004.D0
      GAMMA = RGAS/CP
C  cp_moist=cp*(1.+cpmd*qvp)
      CPMD = .887D0
C  rgas_moist=rgas*(1.+rgasmd*qvp)
      RGASMD = .608D0
C  gamma_moist=gamma*(1.+gammamd*qvp)
      GAMMAMD = RGASMD - CPMD
      TLCLC1 = 2840.D0
      TLCLC2 = 3.5D0
      TLCLC3 = 4.805D0
      TLCLC4 = 55.D0
C  K
      THTECON1 = 3376.D0
      THTECON2 = 2.54D0
      THTECON3 = .81D0
c
c  Calculated the pressure at full sigma levels (a set of pressure
c  levels that bound the layers represented by the vertical grid points)

      CALL DPFCALC(PRS,SFP,PRSF,MIY,MJX,MKZH,TER_FOLLOW)
c
c  Before looping, set lookup table for getting temperature on
c  a pseudoadiabat.
c
      CALL DLOOKUP_TABLE(PSADITHTE,PSADIPRS,PSADITMK,PSAFILE)
c
C   do j=1,mjx-1
      DO J = 1,MJX
C   do i=1,miy-1
          DO I = 1,MIY
              CAPE(I,J,1) = 0.D0
              CIN(I,J,1) = 0.D0
c
              IF (I3DFLAG.EQ.1) THEN
                  KPAR1 = 2
                  KPAR2 = MKZH
              ELSE
c
c      Find parcel with max theta-e in lowest 3 km AGL.
c
                  ETHMAX = -1.D0
                  DO K = MKZH,1,-1
                      IF (GHT(I,J,K)-TER(I,J).LT.3000.D0) THEN
                          Q = MAX(QVP(I,J,K),1.D-15)
                          T = TMK(I,J,K)
                          P = PRS(I,J,K)
                          E = Q*P/ (EPS+Q)
                          TLCL = TLCLC1/ (LOG(T**TLCLC2/E)-TLCLC3) +
     +                           TLCLC4
                          ETH = T* (1000.D0/P)**
     +                          (GAMMA* (1.D0+GAMMAMD*Q))*
     +                          EXP((THTECON1/TLCL-THTECON2)*Q*
     +                          (1.D0+THTECON3*Q))
                          IF (ETH.GT.ETHMAX) THEN
                              KLEV = K
                              ETHMAX = ETH
                          END IF
                      END IF
                  END DO
                  KPAR1 = KLEV
                  KPAR2 = KLEV
c
c      Establish average properties of that parcel
c         (over depth of approximately davg meters)
c
c         davg=.1
                  DAVG = 500.D0
                  PAVG = DAVG*PRS(I,J,KPAR1)*GRAV/
     +                   (RGAS*VIRTUAL(TMK(I,J,KPAR1),QVP(I,J,KPAR1)))
                  P2 = MIN(PRS(I,J,KPAR1)+.5D0*PAVG,PRSF(I,J,MKZH))
                  P1 = P2 - PAVG
                  TOTTHE = 0.D0
                  TOTQVP = 0.D0
                  TOTPRS = 0.D0
                  DO K = MKZH,2,-1
                      IF (PRSF(I,J,K).LE.P1) GO TO 35
                      IF (PRSF(I,J,K-1).GE.P2) GO TO 34
                      P = PRS(I,J,K)
                      PUP = PRSF(I,J,K)
                      PDN = PRSF(I,J,K-1)
                      Q = MAX(QVP(I,J,K),1.D-15)
                      TH = TMK(I,J,K)* (1000.D0/P)**
     +                     (GAMMA* (1.D0+GAMMAMD*Q))
                      PP1 = MAX(P1,PDN)
                      PP2 = MIN(P2,PUP)
                      IF (PP2.GT.PP1) THEN
                          DELTAP = PP2 - PP1
                          TOTQVP = TOTQVP + Q*DELTAP
                          TOTTHE = TOTTHE + TH*DELTAP
                          TOTPRS = TOTPRS + DELTAP
                      END IF
   34                 CONTINUE
                  END DO
   35             CONTINUE
                  QVPPARI = TOTQVP/TOTPRS
                  TMKPARI = (TOTTHE/TOTPRS)*
     +                      (PRS(I,J,KPAR1)/1000.D0)** (GAMMA*
     +                      (1.D0+GAMMAMD*QVP(I,J,KPAR1)))
              END IF
c
              DO KPAR = KPAR1,KPAR2
c
c   Calculate temperature and moisture properties of parcel
c     (Note, qvppari and tmkpari already calculated above for 2D case.)
c
                  IF (I3DFLAG.EQ.1) THEN
                      QVPPARI = QVP(I,J,KPAR)
                      TMKPARI = TMK(I,J,KPAR)
                  END IF
                  PRSPARI = PRS(I,J,KPAR)
                  GHTPARI = GHT(I,J,KPAR)
                  GAMMAM = GAMMA* (1.D0+GAMMAMD*QVPPARI)
                  CPM = CP* (1.D0+CPMD*QVPPARI)
c
                  E = MAX(1.D-20,QVPPARI*PRSPARI/ (EPS+QVPPARI))
                  TLCL = TLCLC1/ (LOG(TMKPARI**TLCLC2/E)-TLCLC3) +
     +                   TLCLC4
                  ETHPARI = TMKPARI* (1000.D0/PRSPARI)**
     +                      (GAMMA* (1.D0+GAMMAMD*QVPPARI))*
     +                      EXP((THTECON1/TLCL-THTECON2)*QVPPARI*
     +                      (1.D0+THTECON3*QVPPARI))
                  ZLCL = GHTPARI + (TMKPARI-TLCL)/ (GRAV/CPM)
c
c   Calculate buoyancy and relative height of lifted parcel at
c   all levels, and store in bottom up arrays.  Add a level at the LCL,
c   and at all points where buoyancy is zero.
c
C  for arrays that go bottom to top
                  KK = 0
                  ILCL = 0
                  IF (GHTPARI.GE.ZLCL) THEN
c
c      initial parcel already saturated or supersaturated.
c
                      ILCL = 2
                      KLCL = 1
                  END IF
                  DO K = KPAR,1,-1
C  for arrays that go bottom to top
   33                 KK = KK + 1
C  model level is below LCL
                      IF (GHT(I,J,K).LT.ZLCL) THEN
                          QVPLIFT = QVPPARI
                          TMKLIFT = TMKPARI - GRAV/CPM*
     +                              (GHT(I,J,K)-GHTPARI)
                          TVENV = VIRTUAL(TMK(I,J,K),QVP(I,J,K))
                          TVLIFT = VIRTUAL(TMKLIFT,QVPLIFT)
                          GHTLIFT = GHT(I,J,K)
                      ELSE IF (GHT(I,J,K).GE.ZLCL .AND. ILCL.EQ.0) THEN
c
c     This model level and previous model level straddle the LCL,
c     so first create a new level in the bottom-up array, at the LCL.
c
                          TMKLIFT = TLCL
                          QVPLIFT = QVPPARI
                          FACDEN = GHT(I,J,K) - GHT(I,J,K+1)
                          FAC1 = (ZLCL-GHT(I,J,K+1))/FACDEN
                          FAC2 = (GHT(I,J,K)-ZLCL)/FACDEN
                          TMKENV = TMK(I,J,K+1)*FAC2 + TMK(I,J,K)*FAC1
                          QVPENV = QVP(I,J,K+1)*FAC2 + QVP(I,J,K)*FAC1
                          TVENV = VIRTUAL(TMKENV,QVPENV)
                          TVLIFT = VIRTUAL(TMKLIFT,QVPLIFT)
                          GHTLIFT = ZLCL
                          ILCL = 1
                      ELSE
                          TMKLIFT = TONPSADIABAT(ETHPARI,PRS(I,J,K),
     +                              PSADITHTE,PSADIPRS,PSADITMK,GAMMA)
                          ESLIFT = EZERO*EXP(ESLCON1* (TMKLIFT-CELKEL)/
     +                             (TMKLIFT-ESLCON2))
                          QVPLIFT = EPS*ESLIFT/ (PRS(I,J,K)-ESLIFT)
                          TVENV = VIRTUAL(TMK(I,J,K),QVP(I,J,K))
                          TVLIFT = VIRTUAL(TMKLIFT,QVPLIFT)
                          GHTLIFT = GHT(I,J,K)
                      END IF
C  buoyancy
                      BUOY(KK) = GRAV* (TVLIFT-TVENV)/TVENV
                      ZREL(KK) = GHTLIFT - GHTPARI
                      IF ((KK.GT.1).AND.
     +                    (BUOY(KK)*BUOY(KK-1).LT.0.0D0)) THEN
c
c   Parcel ascent curve crosses sounding curve, so create a new level
c   in the bottom-up array at the crossing.
c
                          KK = KK + 1
                          BUOY(KK) = BUOY(KK-1)
                          ZREL(KK) = ZREL(KK-1)
                          BUOY(KK-1) = 0.D0
                          ZREL(KK-1) = ZREL(KK-2) +
     +                                 BUOY(KK-2)/ (BUOY(KK-2)-
     +                                 BUOY(KK))* (ZREL(KK)-ZREL(KK-2))
                      END IF
                      IF (ILCL.EQ.1) THEN
                          KLCL = KK
                          ILCL = 2
                          GO TO 33
                      END IF
                  END DO
                  KMAX = KK
                  IF (KMAX.GT.150) THEN
                      print *,
     +                  'capecalc3d: kmax got too big. kmax=',KMAX
                      STOP
                  END IF
c
c   If no LCL was found, set klcl to kmax.  It is probably not really
c   at kmax, but this will make the rest of the routine behave
c   properly.
c
                  IF (ILCL.EQ.0) KLCL=KMAX
c
c   Get the accumulated buoyant energy from the parcel's starting
c   point, at all levels up to the top level.
c
                  BENACCUM(1) = 0.0D0
                  BENAMIN = 9D9
                  DO K = 2,KMAX
                      DZ = ZREL(K) - ZREL(K-1)
                      BENACCUM(K) = BENACCUM(K-1) +
     +                              .5D0*DZ* (BUOY(K-1)+BUOY(K))
                      IF (BENACCUM(K).LT.BENAMIN) THEN
                          BENAMIN = BENACCUM(K)
                      END IF
                  END DO
c
c     Determine equilibrium level (EL), which we define as the highest
c     level of non-negative buoyancy above the LCL. Note, this may be
c     the top level if the parcel is still buoyant there.
c
                  DO K = KMAX,KLCL,-1
                      IF (BUOY(K).GE.0.D0) THEN
C  k of equilibrium level
                          KEL = K
                          GO TO 50
                      END IF
                  END DO
c
c   If we got through that loop, then there is no non-negative
c   buoyancy above the LCL in the sounding.  In these situations,
c   both CAPE and CIN will be set to -0.1 J/kg. (See below about
c   missing values in V6.1.0). Also, where CAPE is
c   non-zero, CAPE and CIN will be set to a minimum of +0.1 J/kg, so
c   that the zero contour in either the CIN or CAPE fields will
c   circumscribe regions of non-zero CAPE.
c
c   In V6.1.0 of NCL, we added a _FillValue attribute to the return
c   value of this function. At that time we decided to change -0.1 
c   to a more appropriate missing value, which is passed into this 
c   routine as CMSG.
c
c                 CAPE(I,J,KPAR) = -0.1D0
c                 CIN(I,J,KPAR) = -0.1D0
                  CAPE(I,J,KPAR) = CMSG
                  CIN(I,J,KPAR)  = CMSG
                  KLFC = KMAX
c
                  GO TO 102
c
   50             CONTINUE
c
c   If there is an equilibrium level, then CAPE is positive.  We'll
c   define the level of free convection (LFC) as the point below the
c   EL, but at or above the LCL, where accumulated buoyant energy is a
c   minimum.  The net positive area (accumulated buoyant energy) from
c   the LFC up to the EL will be defined as the CAPE, and the net
c   negative area (negative of accumulated buoyant energy) from the
c   parcel starting point to the LFC will be defined as the convective
c   inhibition (CIN).
c
c   First get the LFC according to the above definition.
c
                  BENAMIN = 9D9
                  KLFC = KMAX
                  DO K = KLCL,KEL
                      IF (BENACCUM(K).LT.BENAMIN) THEN
                          BENAMIN = BENACCUM(K)
                          KLFC = K
                      END IF
                  END DO
c
c   Now we can assign values to cape and cin
c
                  CAPE(I,J,KPAR) = MAX(BENACCUM(KEL)-BENAMIN,0.1D0)
                  CIN(I,J,KPAR) = MAX(-BENAMIN,0.1D0)
c
c   CIN is uninteresting when CAPE is small (< 100 J/kg), so set
c   CIN to -0.1 (see note about missing values in V6.1.0) in 
c   that case.
c
c   In V6.1.0 of NCL, we added a _FillValue attribute to the return
c   value of this function. At that time we decided to change -0.1 
c   to a more appropriate missing value, which is passed into this 
c   routine as CMSG.
c
C                 IF (CAPE(I,J,KPAR).LT.100.D0) CIN(I,J,KPAR) = -0.1D0
                  IF (CAPE(I,J,KPAR).LT.100.D0) CIN(I,J,KPAR) = CMSG
  102             CONTINUE
c
              END DO
c
              IF (I3DFLAG.EQ.0) THEN
                  CAPE(I,J,MKZH) = CAPE(I,J,KPAR1)
                  CIN(I,J,MKZH) = CIN(I,J,KPAR1)
C  meters AGL
                  CIN(I,J,MKZH-1) = ZREL(KLCL) + GHTPARI - TER(I,J)
C  meters AGL
                  CIN(I,J,MKZH-2) = ZREL(KLFC) + GHTPARI - TER(I,J)
              END IF
c
          END DO
      END DO
c
      RETURN
      END
c                                                                     c
c*********************************************************************c
c                                                                     c
C NCLFORTSTART
      DOUBLE PRECISION FUNCTION TONPSADIABAT(THTE,PRS,PSADITHTE,
     &                                       PSADIPRS,PSADITMK,GAMMA)
      IMPLICIT NONE
      DOUBLE PRECISION THTE
      DOUBLE PRECISION PRS
      DOUBLE PRECISION PSADITHTE
      DOUBLE PRECISION PSADIPRS
      DOUBLE PRECISION PSADITMK
      DOUBLE PRECISION GAMMA
C NCLEND
      DOUBLE PRECISION FRACJT
      DOUBLE PRECISION FRACJT2
      DOUBLE PRECISION FRACIP
      DOUBLE PRECISION FRACIP2
      DIMENSION PSADITHTE(150),PSADIPRS(150),PSADITMK(150,150)
      INTEGER IP, IPCH, JT, JTCH
c                                                                     c
c   This function gives the temperature (in K) on a moist adiabat
c   (specified by thte in K) given pressure in hPa.  It uses a
c   lookup table, with data that was generated by the Bolton (1980)
c   formula for theta_e.
c
c     First check if pressure is less than min pressure in lookup table.
c     If it is, assume parcel is so dry that the given theta-e value can
c     be interpretted as theta, and get temperature from the simple dry
c     theta formula.
c
      IF (PRS.LE.PSADIPRS(150)) THEN
          TONPSADIABAT = THTE* (PRS/1000.D0)**GAMMA
          RETURN
      END IF
c
c   Otherwise, look for the given thte/prs point in the lookup table.
c
      DO JTCH = 1,150 - 1
          IF (THTE.GE.PSADITHTE(JTCH) .AND.
     +        THTE.LT.PSADITHTE(JTCH+1)) THEN
              JT = JTCH
              GO TO 213
          END IF
      END DO
      JT = -1
  213 CONTINUE
      DO IPCH = 1,150 - 1
          IF (PRS.LE.PSADIPRS(IPCH) .AND. PRS.GT.PSADIPRS(IPCH+1)) THEN
              IP = IPCH
              GO TO 215
          END IF
      END DO
      IP = -1
  215 CONTINUE
      IF (JT.EQ.-1 .OR. IP.EQ.-1) THEN
         print *,'capecalc3d: ',
     +           'Outside of lookup table bounds. prs,thte=',
     +      PRS,THTE
          STOP
      END IF
      FRACJT = (THTE-PSADITHTE(JT))/ (PSADITHTE(JT+1)-PSADITHTE(JT))
      FRACJT2 = 1.D0 - FRACJT
      FRACIP = (PSADIPRS(IP)-PRS)/ (PSADIPRS(IP)-PSADIPRS(IP+1))
      FRACIP2 = 1.D0 - FRACIP
      IF (PSADITMK(IP,JT).GT.1D9 .OR. PSADITMK(IP+1,JT).GT.1D9 .OR.
     +    PSADITMK(IP,JT+1).GT.1D9 .OR. PSADITMK(IP+1,JT+1).GT.1D9) THEN
          print *,'capecalc3d: ',
     +      'Tried to access missing temperature in lookup table.',
     +      'Prs and Thte probably unreasonable. prs,thte=',PRS,THTE
          STOP
      END IF
      TONPSADIABAT = FRACIP2*FRACJT2*PSADITMK(IP,JT) +
     +               FRACIP*FRACJT2*PSADITMK(IP+1,JT) +
     +               FRACIP2*FRACJT*PSADITMK(IP,JT+1) +
     +               FRACIP*FRACJT*PSADITMK(IP+1,JT+1)
c
      RETURN
      END
c                                                                     c
c*********************************************************************c
      SUBROUTINE DLOOKUP_TABLE(PSADITHTE,PSADIPRS,PSADITMK,FNAME)
      DOUBLE PRECISION PSADITHTE
      DOUBLE PRECISION PSADIPRS
      DOUBLE PRECISION PSADITMK
c   Set up lookup table for getting temperature on a pseudoadiabat.
c   (Borrow the unit number for the stationlist, just for the moment.)
c
C      CHARACTER*15 FNAME
      CHARACTER*(*) FNAME
      DIMENSION PSADITHTE(150),PSADIPRS(150),PSADITMK(150,150)

C      FNAME = 'psadilookup.dat'
      IUSTNLIST = 33
      OPEN (UNIT=IUSTNLIST,FILE=FNAME,FORM='formatted',STATUS='old')
      DO I = 1,14
          READ (IUSTNLIST,FMT=*)
      END DO
      READ (IUSTNLIST,FMT=*) NTHTE,NPRS
      IF (NTHTE.NE.150 .OR. NPRS.NE.150) THEN
          WRITE (IUP,FMT=*)
     +      'Number of pressure or theta_e levels in lookup table'
          WRITE (IUP,FMT=*) 'file not = 150.  Check lookup table file.'
          STOP
      END IF
      READ (IUSTNLIST,FMT=173) (PSADITHTE(JT),JT=1,NTHTE)
      READ (IUSTNLIST,FMT=173) (PSADIPRS(IP),IP=1,NPRS)
      READ (IUSTNLIST,FMT=173) ((PSADITMK(IP,JT),IP=1,NPRS),JT=1,NTHTE)
  173 FORMAT (5D15.7)
      CLOSE (IUSTNLIST)

      RETURN
      END
c                                                                     c
c*********************************************************************c
c                                                                     c
      SUBROUTINE DPFCALC(PRS,SFP,PF,MIY,MJX,MKZH,TER_FOLLOW)
      DOUBLE PRECISION PRS
      DOUBLE PRECISION SFP
      DOUBLE PRECISION PF
c
c     Historically, this routine calculated the pressure at full sigma
c     levels when RIP was specifically designed for MM4/MM5 output.
c     With the new generalized RIP (Feb '02), this routine is still
c     intended to calculate a set of pressure levels that bound the
c     layers represented by the vertical grid points, although no such
c     layer boundaries are assumed to be defined.  The routine simply
c     uses the midpoint between the pressures of the vertical grid
c     points as the bounding levels.  The array only contains mkzh
c     levels, so the pressure of the top of the uppermost layer is
c     actually excluded.  The kth value of pf is the lower bounding
c     pressure for the layer represented by kth data level.  At the
c     lower bounding level of the lowest model layer, it uses the
c     surface pressure, unless the data set is pressure-level data, in
c     which case it assumes the lower bounding pressure level is as far
c     below the lowest vertical level as the upper bounding pressure
c     level is above.
c
      DIMENSION PRS(MIY,MJX,MKZH),SFP(MIY,MJX),PF(MIY,MJX,MKZH)
      INTEGER TER_FOLLOW
c
C  do j=1,mjx-1  Artifact of MM5
      DO J = 1,MJX
C  do i=1,miy-1  staggered grid
          DO I = 1,MIY
              DO K = 1,MKZH
                  IF (K.EQ.MKZH) THEN
C  terrain-following data
                      IF (TER_FOLLOW.EQ.1) THEN
                          PF(I,J,K) = SFP(I,J)
C  pressure-level data
                      ELSE
                          PF(I,J,K) = .5D0* (3.D0*PRS(I,J,K)-
     +                                PRS(I,J,K-1))
                      END IF
                  ELSE
                      PF(I,J,K) = .5D0* (PRS(I,J,K+1)+PRS(I,J,K))
                  END IF
              END DO
          END DO
      END DO
c
      RETURN
      END
c======================================================================
c
c !IROUTINE: VIRTUAL -- Calculate virtual temperature (K)
c
c !DESCRIPTION:
c
c   This function returns a single value of virtual temperature in
c   K, given temperature in K and mixing ratio in kg/kg.  For an
c   array of virtual temperatures, use subroutine VIRTUAL_TEMP.
c
c !INPUT:
c    RATMIX - water vapor mixing ratio (kg/kg)
c    TEMP   - temperature (K)
c
c !OUTPUT:
c    TV     - Virtual temperature (K)
c
c !ASSUMPTIONS:
c
c !REVISION HISTORY:
c     2009-March  - Mark T. Stoelinga - from RIP4.5
c     2010-August - J. Schramm - modified to run with NCL and ARW wrf output
c
c ------------------------------------------------------------------
C NCLFORTSTART
      DOUBLE PRECISION FUNCTION VIRTUAL(TEMP,RATMIX)
      IMPLICIT NONE
      DOUBLE PRECISION TEMP,RATMIX
C NCLEND
      DOUBLE PRECISION EPS
      EPS = 0.622D0
      VIRTUAL = TEMP* (EPS+RATMIX)/ (EPS* (1.D0+RATMIX))
      RETURN
      END
