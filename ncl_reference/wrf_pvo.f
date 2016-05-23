c--------------------------------------------------------
C NCLFORTSTART
      SUBROUTINE DCOMPUTEPV(PV,U,V,THETA,PRS,MSFU,MSFV,MSFT,COR,DX,DY,
     +                      NX,NY,NZ,NXP1,NYP1)

      IMPLICIT NONE
      INTEGER NX,NY,NZ,NXP1,NYP1
      DOUBLE PRECISION U(NXP1,NY,NZ),V(NX,NYP1,NZ),PRS(NX,NY,NZ)
      DOUBLE PRECISION THETA(NX,NY,NZ),PV(NX,NY,NZ)
      DOUBLE PRECISION MSFU(NXP1,NY),MSFV(NX,NYP1),MSFT(NX,NY)
      DOUBLE PRECISION COR(NX,NY)
      DOUBLE PRECISION DX,DY
C NCLEND
      INTEGER KP1,KM1,JP1,JM1,IP1,IM1,I,J,K
      DOUBLE PRECISION DSY,DSX,DP,DUDY,DVDX,DUDP,DVDP,DTHDP,AVORT
      DOUBLE PRECISION DTHDX,DTHDY,MM

c          print*,'nx,ny,nz,nxp1,nyp1'
c          print*,nx,ny,nz,nxp1,nyp1
      DO K = 1,NZ
          KP1 = MIN(K+1,NZ)
          KM1 = MAX(K-1,1)
          DO J = 1,NY
              JP1 = MIN(J+1,NY)
              JM1 = MAX(J-1,1)
              DO I = 1,NX
                  IP1 = MIN(I+1,NX)
                  IM1 = MAX(I-1,1)
c         print *,jp1,jm1,ip1,im1
                  DSX = (IP1-IM1)*DX
                  DSY = (JP1-JM1)*DY
                  MM = MSFT(I,J)*MSFT(I,J)
c         print *,j,i,u(i,jp1,k),msfu(i,jp1),u(i,jp1,k)/msfu(i,jp1)
                  DUDY = 0.5D0* (U(I,JP1,K)/MSFU(I,JP1)+
     +                   U(I+1,JP1,K)/MSFU(I+1,JP1)-
     +                   U(I,JM1,K)/MSFU(I,JM1)-
     +                   U(I+1,JM1,K)/MSFU(I+1,JM1))/DSY*MM
                  DVDX = 0.5D0* (V(IP1,J,K)/MSFV(IP1,J)+
     +                   V(IP1,J+1,K)/MSFV(IP1,J+1)-
     +                   V(IM1,J,K)/MSFV(IM1,J)-
     +                   V(IM1,J+1,K)/MSFV(IM1,J+1))/DSX*MM
                  AVORT = DVDX - DUDY + COR(I,J)
                  DP = PRS(I,J,KP1) - PRS(I,J,KM1)
                  DUDP = 0.5D0* (U(I,J,KP1)+U(I+1,J,KP1)-U(I,J,KM1)-
     +                   U(I+1,J,KM1))/DP
                  DVDP = 0.5D0* (V(I,J,KP1)+V(I,J+1,KP1)-V(I,J,KM1)-
     +                   V(I,J+1,KM1))/DP
                  DTHDP = (THETA(I,J,KP1)-THETA(I,J,KM1))/DP
                  DTHDX = (THETA(IP1,J,K)-THETA(IM1,J,K))/DSX*MSFT(I,J)
                  DTHDY = (THETA(I,JP1,K)-THETA(I,JM1,K))/DSY*MSFT(I,J)
                  PV(I,J,K) = -9.81D0* (DTHDP*AVORT-DVDP*DTHDX+
     +                        DUDP*DTHDY)*10000.D0
c               if(i.eq.300 .and. j.eq.300) then
c                 print*,'avort,dudp,dvdp,dthdp,dthdx,dthdy,pv'
c                 print*,avort,dudp,dvdp,dthdp,dthdx,dthdy,pv(i,j,k)
c               endif
                  PV(I,J,K) = PV(I,J,K)*1.D2
              END DO
          END DO
      END DO
      RETURN
      END

c--------------------------------------------------------
C NCLFORTSTART
      SUBROUTINE DCOMPUTEABSVORT(AV,U,V,MSFU,MSFV,MSFT,COR,DX,DY,NX,NY,
     +                           NZ,NXP1,NYP1)

      IMPLICIT NONE
      INTEGER NX,NY,NZ,NXP1,NYP1
      DOUBLE PRECISION U(NXP1,NY,NZ),V(NX,NYP1,NZ)
      DOUBLE PRECISION AV(NX,NY,NZ)
      DOUBLE PRECISION MSFU(NXP1,NY),MSFV(NX,NYP1),MSFT(NX,NY)
      DOUBLE PRECISION COR(NX,NY)
      DOUBLE PRECISION DX,DY
C NCLEND
      INTEGER KP1,KM1,JP1,JM1,IP1,IM1,I,J,K
      DOUBLE PRECISION DSY,DSX,DP,DUDY,DVDX,DUDP,DVDP,DTHDP,AVORT
      DOUBLE PRECISION DTHDX,DTHDY,MM

c          print*,'nx,ny,nz,nxp1,nyp1'
c          print*,nx,ny,nz,nxp1,nyp1
      DO K = 1,NZ
          DO J = 1,NY
              JP1 = MIN(J+1,NY)
              JM1 = MAX(J-1,1)
              DO I = 1,NX
                  IP1 = MIN(I+1,NX)
                  IM1 = MAX(I-1,1)
c         print *,jp1,jm1,ip1,im1
                  DSX = (IP1-IM1)*DX
                  DSY = (JP1-JM1)*DY
                  MM = MSFT(I,J)*MSFT(I,J)
c         print *,j,i,u(i,jp1,k),msfu(i,jp1),u(i,jp1,k)/msfu(i,jp1)
                  DUDY = 0.5D0* (U(I,JP1,K)/MSFU(I,JP1)+
     +                   U(I+1,JP1,K)/MSFU(I+1,JP1)-
     +                   U(I,JM1,K)/MSFU(I,JM1)-
     +                   U(I+1,JM1,K)/MSFU(I+1,JM1))/DSY*MM
                  DVDX = 0.5D0* (V(IP1,J,K)/MSFV(IP1,J)+
     +                   V(IP1,J+1,K)/MSFV(IP1,J+1)-
     +                   V(IM1,J,K)/MSFV(IM1,J)-
     +                   V(IM1,J+1,K)/MSFV(IM1,J+1))/DSX*MM
                  AVORT = DVDX - DUDY + COR(I,J)
                  AV(I,J,K) = AVORT*1.D5
              END DO
          END DO
      END DO
      RETURN
      END
