c ----------------------------------------------------------- 
C NCLFORTSTART
      SUBROUTINE DRCM2POINTS(NGRD,NYI,NXI,YI,XI,FI,NXYO,YO,XO,FO
     +                      ,XMSG,OPT,NCRIT,KVAL,IER)
      IMPLICIT NONE
      INTEGER NGRD,NXI,NYI,NXYO,OPT,NCRIT,KVAL,IER
      DOUBLE PRECISION XI(NXI,NYI),YI(NXI,NYI),FI(NXI,NYI,NGRD)
      DOUBLE PRECISION XO(NXYO),YO(NXYO),FO(NXYO,NGRD),XMSG
C NCLEND

C This is written  with GNU f77 acceptable extensions
c .   This could be improved considerably with f90

c nomenclature:
c .   nxi,nyi - lengths of xi,yi and dimensions of fi (must be >= 2)
c .   xi      - coordinates of fi (eg, lon [2D] )
c .   yi      - coordinates of fi (eg, lat [2D] )
c .   fi      - functional input values [2D]
c .   nxyo    - number of output points
c .   xo      - lon coordinates of fo (eg, lon [1D])
c .   yo      - lat coordinates of fo (eg, lat [1D])
c .   fo      - functional output values [interpolated]
c .   xmsg    - missing code
c .   opt     - 0/1 = inv distance, 2 = bilinear
c .   ier     - error code
c .             =0;   no error
c .             =1;   not enough points in input/output array
c .             =2/3; xi or yi are not monotonically increasing
c .             =4/5; xo or yo are not monotonically increasing
c
c                              local
      INTEGER NG,NX,NY,NXY,NEXACT,IX,IY,M,N,NW,NER,K
      DOUBLE PRECISION FW(2,2),W(2,2),SUMF,SUMW,CHKLAT(NYI),CHKLON(NXI)
      DOUBLE PRECISION DGCDIST, WX, WY
      DOUBLE PRECISION REARTH, DLAT, PI, RAD, DKM, DIST 
c                              error checking
      IER = 0
      IF (NXI.LE.1 .OR. NYI.LE.1 .OR. NXYO.LE.0) THEN
          IER = 1
          RETURN
      END IF
      IF (IER.NE.0) RETURN

      DO NY = 1,NYI
          CHKLAT(NY) = YI(1,NY)
c c c    print *,"chklat: ny=",ny,"  chklat=",chklat(ny)
      END DO
      CALL DMONOINC(CHKLAT,NYI,IER,NER)
      IF (IER.NE.0) RETURN

      DO NX = 1,NXI
          CHKLON(NX) = XI(NX,1)
c c c    print *,"chklon: nx=",nx,"  chklon=",chklon(nx)
      END DO
      CALL DMONOINC(CHKLAT,NYI,IER,NER)
      IF (IER.NE.0) RETURN

C ORIGINAL  (k = op, never implemented)
      IF (KVAL.LE.0) THEN
         K = 1
      ELSE
         K = KVAL
      END IF
      DO NG = 1,NGRD
        DO NXY = 1,NXYO
           FO(NXY,NG) = XMSG
        END DO
      END DO
c                              main loop [exact matches]
      NEXACT = 0
      DO NXY = 1,NXYO

          DO IY = 1,NYI
              DO IX = 1,NXI
                  IF (XO(NXY).EQ.XI(IX,IY) .AND.
     +                YO(NXY).EQ.YI(IX,IY)) THEN
                      DO NG = 1,NGRD
                         FO(NXY,NG) = FI(IX,IY,NG)
                         NEXACT     = NEXACT + 1
                      END DO
                      GO TO 10
                  END IF
              END DO
          END DO

   10     CONTINUE
      END DO

c c c print *, "nexact=",nexact
c                              main loop [interpolation]
      DO NXY = 1,NXYO

              DO IY = 1,NYI - K
                DO IX = 1,NXI - K
                   IF (XO(NXY).GE.XI(IX,IY) .AND.
     +                 XO(NXY).LE.XI(IX+K,IY) .AND.
     +                 YO(NXY).GE.YI(IX,IY) .AND.
     +                 YO(NXY).LE.YI(IX,IY+K)) THEN

                   IF (ABS(OPT).EQ.2) THEN
                       WX = (XO(NXY)-XI(IX,IY))/
     +                      (XI(IX+K,IY)-XI(IX,IY))
                       WY = (YO(NXY)-YI(IX,IY))/
     +                      (YI(IX,IY+K)-YI(IX,IY))
                       W(1,1) = (1.D0-WX)*(1.D0-WY)
                       W(2,1) = WX*(1.D0-WY)
                       W(1,2) = (1.D0-WX)*WY
                       W(2,2) = WX*WY
                   ELSE
                       W(1,1) = (1.D0/DGCDIST(YO(NXY),XO(NXY),
     +                           YI(IX,IY),XI(IX,IY),2))**2
                       W(2,1) = (1.D0/DGCDIST(YO(NXY),XO(NXY),
     +                           YI(IX+K,IY),XI(IX+K,IY),2))**2
                       W(1,2) = (1.D0/DGCDIST(YO(NXY),XO(NXY),
     +                           YI(IX,IY+K),XI(IX,IY+K),2))**2
                       W(2,2) = (1.D0/DGCDIST(YO(NXY),XO(NXY),
     +                           YI(IX+K,IY+K),XI(IX+K,IY+K),2))**2
                   END IF

                   DO NG = 1,NGRD
                      IF (FO(NXY,NG).EQ.XMSG) THEN
                          
                          FW(1,1) = FI(IX,IY,NG)
                          FW(2,1) = FI(IX+K,IY,NG)
                          FW(1,2) = FI(IX,IY+K,NG)
                          FW(2,2) = FI(IX+K,IY+K,NG)

                          NW = 0
                          SUMF = 0.0D0
                          SUMW = 0.0D0
                          DO N = 1,2
                              DO M = 1,2
                                  IF (FW(M,N).NE.XMSG) THEN
                                      SUMF = SUMF + FW(M,N)*W(M,N)
                                      SUMW = SUMW + W(M,N)
                                      NW = NW + 1
                                  END IF
                              END DO
                          END DO
c                                             nw >=3 arbitrary
                          IF (NW.GE.NCRIT .AND. SUMW.GT.0.D0) THEN
                              FO(NXY,NG) = SUMF/SUMW
                          END IF
                      END IF
                    END DO
                    GO TO 20
                  END IF
                END DO
              END DO

   20         CONTINUE
      END DO

C Are all the output points filled in? Check the 1st grid
C If so, return

      DO NG = 1,NGRD   
        DO NXY = 1,NXYO
           IF (FO(NXY,NG).EQ.XMSG) GO TO 30
        END DO
      END DO
      RETURN

C only enter if some points are not interpolated to
C DLAT is arbitrary.  It ould be made an option.
C DLAT is expressed in terms of degrees of latitude.
C DKM  is DLAT in KILOMETERS

   30 REARTH= 6371D0
      DLAT  = 5  
      PI    = 4D0*ATAN(1.0D0)
      RAD   = PI/180D0
      DKM   = DLAT*(2D0*PI*REARTH)/360D0

C LOOP OVER EACH GRID ... INEFFICIENT 
C THE RUB IS THAT SOME LEVELS COULD HAVE XMSG.

      DO NG = 1,NGRD   

        DO NXY = 1,NXYO
           IF(FO(NXY,NG).EQ.XMSG) THEN

C FIND ALL GRID POINTS WITHIN 'DKM' KILOMETERS OF PT 

              NW   = 0
              SUMF = 0.0D0
              SUMW = 0.0D0

              DO IY = 1,NYI
                DO IX = 1,NXI
                   IF ((YI(IX,IY).GE.YO(NXY)-DLAT)  .AND.
     +                 (YI(IX,IY).LE.YO(NXY)+DLAT)) THEN      
                        DIST = DGCDIST(YO(NXY),XO(NXY) 
     +                                ,YI(IX,IY),XI(IX,IY),2)
                        IF (DIST.LE.DKM .AND. DIST.GT.0.0D0 .AND.
     +                      FI(IX,IY,NG).NE.XMSG) THEN
                            DIST = 1.0D0/DIST**2
                            SUMF = SUMF + FI(IX,IY,NG)*DIST
                            SUMW = SUMW + DIST
                            NW   = NW + 1
                        END IF
                   END IF
                END DO
              END DO

C C C         IF (NW.GE.NCRIT .AND. SUMW.GT. 0.0D0) THEN
              IF (SUMW.GT.0.0D0) THEN
                  FO(NXY,NG) = SUMF/SUMW
              END IF
           END IF
        END DO
      END DO

      RETURN
      END
