C NCLFORTSTART
      SUBROUTINE CALCDBZ(DBZ,PRS,TMK,QVP,QRA,QSN,QGR,WEDIM,SNDIM,BTDIM,
     +                   SN0,IVARINT,ILIQSKIN)
c
c     This routine computes equivalent reflectivity factor (in dBZ) at
c     each model grid point.  In calculating Ze, the RIP algorithm makes
c     assumptions consistent with those made in an early version
c     (ca. 1996) of the bulk mixed-phase microphysical scheme in the MM5
c     model (i.e., the scheme known as "Resiner-2").  For each species:
c
c     1. Particles are assumed to be spheres of constant density.  The
c     densities of rain drops, snow particles, and graupel particles are
c     taken to be rho_r = rho_l = 1000 kg m^-3, rho_s = 100 kg m^-3, and
c     rho_g = 400 kg m^-3, respectively. (l refers to the density of
c     liquid water.)
c
c     2. The size distribution (in terms of the actual diameter of the
c     particles, rather than the melted diameter or the equivalent solid
c     ice sphere diameter) is assumed to follow an exponential
c     distribution of the form N(D) = N_0 * exp( lambda*D ).
c
c     3. If ivarint=0, the intercept parameters are assumed constant
c     (as in early Reisner-2), with values of 8x10^6, 2x10^7, 
c    and 4x10^6 m^-4, for rain, snow, and graupel, respectively.
c    If ivarint=1, variable intercept parameters are used, as 
c    calculated in Thompson, Rasmussen, and Manning (2004, Monthly
c    Weather Review, Vol. 132, No. 2, pp. 519-542.)
c
c     4. If iliqskin=1, frozen particles that are at a temperature above
c     freezing are assumed to scatter as a liquid particle.
c
c     More information on the derivation of simulated reflectivity in
c     RIP can be found in Stoelinga (2005, unpublished write-up).
c     Contact Mark Stoelinga (stoeling@atmos.washington.edu) for a copy.
c

c   Arguments
      INTEGER WEDIM,SNDIM,BTDIM
      INTEGER SN0,IVARINT,ILIQSKIN
      DOUBLE PRECISION DBZ(WEDIM,SNDIM,BTDIM)
      DOUBLE PRECISION PRS(WEDIM,SNDIM,BTDIM)
      DOUBLE PRECISION TMK(WEDIM,SNDIM,BTDIM)
      DOUBLE PRECISION QVP(WEDIM,SNDIM,BTDIM)
      DOUBLE PRECISION QRA(WEDIM,SNDIM,BTDIM)
      DOUBLE PRECISION QSN(WEDIM,SNDIM,BTDIM)
      DOUBLE PRECISION QGR(WEDIM,SNDIM,BTDIM)

C NCLEND

c   Local Variables
      INTEGER I,J,K
      DOUBLE PRECISION TEMP_C,VIRTUAL_T
      DOUBLE PRECISION GONV,RONV,SONV
      DOUBLE PRECISION FACTOR_G,FACTOR_R,FACTOR_S
      DOUBLE PRECISION FACTORB_G,FACTORB_R,FACTORB_S
      DOUBLE PRECISION RHOAIR,Z_E

c   Constants used to calculate variable intercepts
      DOUBLE PRECISION R1,RON,RON2,SON,GON
      DOUBLE PRECISION RON_MIN,RON_QR0,RON_DELQR0
      DOUBLE PRECISION RON_CONST1R,RON_CONST2R
c   Constant intercepts
      DOUBLE PRECISION RN0_R,RN0_S,RN0_G
c   Other constants
      DOUBLE PRECISION RHO_R,RHO_S,RHO_G
      DOUBLE PRECISION GAMMA_SEVEN,ALPHA
      DOUBLE PRECISION RHOWAT,CELKEL,PI,RD


c   Constants used to calculate variable intercepts
      R1 = 1.D-15
      RON = 8.D6
      RON2 = 1.D10
      SON = 2.D7
      GON = 5.D7
      RON_MIN = 8.D6
      RON_QR0 = 0.00010D0
      RON_DELQR0 = 0.25D0*RON_QR0
      RON_CONST1R = (RON2-RON_MIN)*0.5D0
      RON_CONST2R = (RON2+RON_MIN)*0.5D0

c   Constant intercepts
      RN0_R = 8.D6
      RN0_S = 2.D7
      RN0_G = 4.D6

c   Other constants
      GAMMA_SEVEN = 720.D0
      RHOWAT = 1000.D0
      RHO_R = RHOWAT
      RHO_S = 100.D0
      RHO_G = 400.D0
      ALPHA = 0.224D0
      CELKEL = 273.15D0
      PI = 3.141592653589793D0
      RD = 287.04D0

c   Force all Q arrays to be 0.0 or greater.
      DO K = 1,BTDIM
         DO J = 1,SNDIM
            DO I = 1,WEDIM
               IF (QVP(I,J,K).LT.0.0) THEN
                  QVP(I,J,K) = 0.0
               END IF
               IF (QRA(I,J,K).LT.0.0) THEN
                  QRA(I,J,K) = 0.0
               END IF
               IF (QSN(I,J,K).LT.0.0) THEN
                  QSN(I,J,K) = 0.0
               END IF
               IF (QGR(I,J,K).LT.0.0) THEN
                  QGR(I,J,K) = 0.0
               END IF
            END DO
         END DO
      END DO

c   Input pressure is Pa, but we need hPa in calculations

      IF (SN0.EQ.0) THEN
          DO K = 1,BTDIM
              DO J = 1,SNDIM
                  DO I = 1,WEDIM
                      IF (TMK(I,J,K).LT.CELKEL) THEN
                          QSN(I,J,K) = QRA(I,J,K)
                          QRA(I,J,K) = 0.D0
                      END IF
                  END DO
              END DO
          END DO
      END IF


      FACTOR_R = GAMMA_SEVEN*1.D18* (1.D0/ (PI*RHO_R))**1.75D0
      FACTOR_S = GAMMA_SEVEN*1.D18* (1.D0/ (PI*RHO_S))**1.75D0*
     +           (RHO_S/RHOWAT)**2*ALPHA
      FACTOR_G = GAMMA_SEVEN*1.D18* (1.D0/ (PI*RHO_G))**1.75D0*
     +           (RHO_G/RHOWAT)**2*ALPHA


      DO K = 1,BTDIM
          DO J = 1,SNDIM
              DO I = 1,WEDIM

                  VIRTUAL_T = TMK(I,J,K)* (0.622D0+QVP(I,J,K))/
     +                        (0.622D0* (1.D0+QVP(I,J,K)))
                  RHOAIR = PRS(I,J,K) / (RD*VIRTUAL_T)

c      Adjust factor for brightband, where snow or graupel particle
c      scatters like liquid water (alpha=1.0) because it is assumed to
c      have a liquid skin.

                  IF (ILIQSKIN.EQ.1 .AND. TMK(I,J,K).GT.CELKEL) THEN
                      FACTORB_S = FACTOR_S/ALPHA
                      FACTORB_G = FACTOR_G/ALPHA
                  ELSE
                      FACTORB_S = FACTOR_S
                      FACTORB_G = FACTOR_G
                  END IF

c      Calculate variable intercept parameters

                  IF (IVARINT.EQ.1) THEN

                      TEMP_C = DMIN1(-0.001D0,TMK(I,J,K)-CELKEL)
                      SONV = DMIN1(2.0D8,2.0D6*EXP(-0.12D0*TEMP_C))

                      GONV = GON
                      IF (QGR(I,J,K).GT.R1) THEN
                          GONV = 2.38D0* (PI*RHO_G/
     +                           (RHOAIR*QGR(I,J,K)))**0.92D0
                          GONV = MAX(1.D4,MIN(GONV,GON))
                      END IF

                      RONV = RON2
                      IF (QRA(I,J,K).GT.R1) THEN
                          RONV = RON_CONST1R*TANH((RON_QR0-QRA(I,J,K))/
     +                           RON_DELQR0) + RON_CONST2R
                      END IF

                  ELSE

                      RONV = RN0_R
                      SONV = RN0_S
                      GONV = RN0_G

                  END IF

c      Total equivalent reflectivity factor (z_e, in mm^6 m^-3) is
c      the sum of z_e for each hydrometeor species:

                  Z_E = FACTOR_R* (RHOAIR*QRA(I,J,K))**1.75D0/
     +                  RONV**.75D0 + FACTORB_S*
     +                  (RHOAIR*QSN(I,J,K))**1.75D0/SONV**.75D0 +
     +                  FACTORB_G* (RHOAIR*QGR(I,J,K))**1.75D0/
     +                  GONV**.75D0

c      Adjust small values of Z_e so that dBZ is no lower than -30
                  Z_E = MAX(Z_E,.001D0)

c      Convert to dBZ
                  DBZ(I,J,K) = 10.D0*LOG10(Z_E)

              END DO
          END DO
      END DO

      RETURN
      END
