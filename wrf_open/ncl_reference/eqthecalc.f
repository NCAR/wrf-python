      SUBROUTINE DEQTHECALC(QVP,TMK,PRS,ETH,MIY,MJX,MKZH)
      DOUBLE PRECISION EPS
      DOUBLE PRECISION RGAS
      DOUBLE PRECISION RGASMD
      DOUBLE PRECISION CP
      DOUBLE PRECISION CPMD
      DOUBLE PRECISION GAMMA
      DOUBLE PRECISION GAMMAMD
      DOUBLE PRECISION TLCLC1
      DOUBLE PRECISION TLCLC2
      DOUBLE PRECISION TLCLC3
      DOUBLE PRECISION TLCLC4
      DOUBLE PRECISION THTECON1
      DOUBLE PRECISION THTECON2
      DOUBLE PRECISION THTECON3
      DOUBLE PRECISION Q
      DOUBLE PRECISION T
      DOUBLE PRECISION P
      DOUBLE PRECISION E
      DOUBLE PRECISION TLCL
c
c Input variables
c Qvapor [g/kg]
      DOUBLE PRECISION QVP(MIY,MJX,MKZH)
c Temperature [K]
      DOUBLE PRECISION TMK(MIY,MJX,MKZH)
c full pressure (=P+PB) [hPa]
      DOUBLE PRECISION PRS(MIY,MJX,MKZH)
c
c Output variable
c equivalent potential temperature [K]
      DOUBLE PRECISION ETH(MIY,MJX,MKZH)
c
c parameters
      PARAMETER (EPS=0.622D0)

c J/K/kg
      RGAS = 287.04D0
c rgas_moist=rgas*(1.+rgasmd*qvp)
      RGASMD = .608D0
c J/K/kg  Note: not using Bolton's value of 1005.7
      CP = 1004.D0
c cp_moist=cp*(1.+cpmd*qvp)
      CPMD = .887D0
      GAMMA = RGAS/CP
c gamma_moist=gamma*(1.+gammamd*qvp)
      GAMMAMD = RGASMD - CPMD

      TLCLC1 = 2840.D0
      TLCLC2 = 3.5D0
      TLCLC3 = 4.805D0
      TLCLC4 = 55.D0
c K
      THTECON1 = 3376.D0
      THTECON2 = 2.54D0
      THTECON3 = .81D0
c
      DO 1000 K = 1,MKZH
          DO 1000 J = 1,MJX
              DO 1000 I = 1,MIY
                  Q = MAX(QVP(I,J,K),1.D-15)
                  T = TMK(I,J,K)
                  P = PRS(I,J,K)/100.
                  E = Q*P/ (EPS+Q)
                  TLCL = TLCLC1/ (LOG(T**TLCLC2/E)-TLCLC3) + TLCLC4
                  ETH(I,J,K) = T* (1000.D0/P)**
     +                         (GAMMA* (1.D0+GAMMAMD*Q))*
     +                         EXP((THTECON1/TLCL-THTECON2)*Q*
     +                         (1.D0+THTECON3*Q))
 1000 CONTINUE
      RETURN
      END
