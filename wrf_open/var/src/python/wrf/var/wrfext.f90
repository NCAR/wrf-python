
SUBROUTINE f_interpz3d(data3d,zdata,desiredloc,missingval,out2d,nx,ny,nz)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx,ny,nz
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) ::  data3d
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(OUT) :: out2d
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: zdata
    REAL(KIND=8), INTENT(IN) :: desiredloc
    REAL(KIND=8), INTENT(IN) :: missingval

    INTEGER :: i,j,kp,ip,im
    LOGICAL :: dointerp
    REAL(KIND=8) :: height,w1,w2

    height = desiredloc

    ! does vertical coordinate increase or decrease with increasing k?
    ! set offset appropriately

    ip = 0
    im = 1
    IF (zdata(1,1,1).GT.zdata(1,1,nz)) THEN
        ip = 1
        im = 0
    END IF

    DO i = 1,nx
        DO j = 1,ny
            ! Initialize to missing.  Was initially hard-coded to -999999.
            out2d(i,j) = missingval
            dointerp = .FALSE.
            kp = nz

            DO WHILE ((.NOT. dointerp) .AND. (kp >= 2))
                IF (((zdata(i,j,kp-im) < height) .AND. (zdata(i,j,kp-ip) > height))) THEN
                    w2 = (height-zdata(i,j,kp-im))/(zdata(i,j,kp-ip)-zdata(i,j,kp-im))
                    w1 = 1.D0 - w2
                    out2d(i,j) = w1*data3d(i,j,kp-im) + w2*data3d(i,j,kp-ip)
                    dointerp = .TRUE.
                END IF
                kp = kp - 1
            END DO

        END DO
    END DO

    RETURN

END SUBROUTINE f_interpz3d

SUBROUTINE f_interp2dxy(v3d,xy,v2d,nx,ny,nz,nxy)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nx,ny,nz,nxy
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: v3d
    REAL(KIND=8),DIMENSION(nxy,nz),INTENT(OUT) :: v2d
    REAL(KIND=8),DIMENSION(2,nxy),INTENT(IN) :: xy

    INTEGER :: i,j,k,ij
    REAL(KIND=8) :: w11,w12,w21,w22,wx,wy

    DO ij = 1,nxy
      i = MAX0(1,MIN0(nx-1,INT(xy(1,ij)+1)))
      j = MAX0(1,MIN0(ny-1,INT(xy(2,ij)+1)))
      wx = DBLE(i+1) - (xy(1,ij)+1)
      wy = DBLE(j+1) - (xy(2,ij)+1)
      w11 = wx*wy
      w21 = (1.D0-wx)*wy
      w12 = wx* (1.D0-wy)
      w22 = (1.D0-wx)* (1.D0-wy)
      DO k = 1,nz
          v2d(ij,k) = w11*v3d(i,j,k) + w21*v3d(i+1,j,k) + &
              w12*v3d(i,j+1,k) + w22*v3d(i+1,j+1,k)
      END DO
    END DO

    RETURN

END SUBROUTINE f_interp2dxy

SUBROUTINE f_interp1d(v_in,z_in,z_out,vmsg,v_out,nz_in,nz_out)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nz_in, nz_out
    REAL(KIND=8),DIMENSION(nz_in),INTENT(IN) :: v_in,z_in
    REAL(KIND=8),DIMENSION(nz_out),INTENT(IN) :: z_out
    REAL(KIND=8),DIMENSION(nz_out),INTENT(OUT) :: v_out
    REAL(KIND=8),INTENT(IN) :: vmsg

    INTEGER :: kp,k,im,ip
    LOGICAL :: interp
    REAL(KIND=8) :: height,w1,w2

    ! does vertical coordinate increase of decrease with increasing k?
    ! set offset appropriately

    ip = 0
    im = 1
    IF (z_in(1).GT.z_in(nz_in)) THEN
      ip = 1
      im = 0
    END IF

    DO k = 1,nz_out
      v_out(k) = vmsg

      interp = .FALSE.
      kp = nz_in
      height = z_out(k)

      DO WHILE ((.NOT.interp) .AND. (kp.GE.2))
          IF (((z_in(kp-im).LE.height).AND.(z_in(kp-ip).GT.height))) THEN
              w2 = (height-z_in(kp-im))/(z_in(kp-ip)-z_in(kp-im))
              w1 = 1.d0 - w2
              v_out(k) = w1*v_in(kp-im) + w2*v_in(kp-ip)
              interp = .TRUE.
          END IF
          kp = kp - 1
      END DO
    END DO

    RETURN

END SUBROUTINE f_interp1d

! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!
!    units: pressure (Pa), temperature(K), height (m), mixing ratio
!     (kg kg{-1}) availability of 3d p, t, and qv; 2d terrain; 1d
! half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!

SUBROUTINE f_computeslp(z,t,p,q,t_sea_level,t_surf,level,throw_exception,&
                        sea_level_pressure,nx,ny,nz)

    IMPLICIT NONE

    EXTERNAL throw_exception
    !     Estimate sea level pressure.
    INTEGER, INTENT(IN) :: nx,ny,nz
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: z
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: t,p,q
    !     The output is the 2d sea level pressure.
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(OUT) :: sea_level_pressure
    INTEGER,DIMENSION(nx,ny), INTENT(INOUT) ::  level
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(INOUT) :: t_surf,t_sea_level

    !     Some required physical constants:

    REAL(KIND=8), PARAMETER :: R=287.04D0, G=9.81D0, GAMMA=0.0065D0

    !     Specific constants for assumptions made in this routine:

    REAL(KIND=8), PARAMETER :: TC=273.16D0+17.5D0, PCONST=10000

    LOGICAL, PARAMETER :: ridiculous_mm5_test=.TRUE.
    !      PARAMETER (ridiculous_mm5_test = .FALSE.)

    !     Local variables:

    INTEGER :: i,j,k
    INTEGER :: klo,khi

    REAL(KIND=8) :: plo,phi,tlo,thi,zlo,zhi
    REAL(KIND=8) :: p_at_pconst,t_at_pconst,z_at_pconst
    REAL(KIND=8) :: z_half_lowest

    LOGICAL :: l1,l2,l3,found

    !  Find least zeta level that is PCONST Pa above the surface.  We
    !  later use this level to extrapolate a surface pressure and
    !  temperature, which is supposed to reduce the effect of the diurnal
    !  heating cycle in the pressure field.

    DO j = 1,ny
      DO i = 1,nx
          level(i,j) = -1

          k = 1
          found = .FALSE.
          DO WHILE ((.NOT. found) .AND. (k <= nz))
              IF (p(i,j,k) < p(i,j,1)-PCONST) THEN
                  level(i,j) = k
                  found = .TRUE.
              END IF
              k = k + 1
          END DO

          IF (level(i,j) == -1) THEN
              !PRINT '(A,I4,A)','Troubles finding level ', NINT(PCONST)/100,' above ground.'
              !PRINT '(A,I4,A,I4,A)','Problems first occur at (',I,',',J,')'
              !PRINT '(A,F6.1,A)','Surface pressure = ',p(i,j,1)/100,' hPa.'
              CALL throw_exception('Error in finding 100 hPa up')
          END IF
      END DO
    END DO

    !     Get temperature PCONST Pa above surface.  Use this to extrapolate
    !     the temperature at the surface and down to sea level.

    DO J = 1,ny
      DO I = 1,nx

          klo = MAX(level(i,j)-1,1)
          khi = MIN(klo+1,nz-1)

          IF (klo == khi) THEN
              !PRINT '(A)','Trapping levels are weird.'
              !PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi,': and they should not be equal.'
              CALL throw_exception('Error trapping levels')
          END IF

          plo = p(i,j,klo)
          phi = p(i,j,khi)
          tlo = t(i,j,klo)* (1.D0+0.608D0*q(i,j,klo))
          thi = t(i,j,khi)* (1.D0+0.608D0*q(i,j,khi))
    !         zlo = zetahalf(klo)/ztop*(ztop-terrain(i,j))+terrain(i,j)
    !         zhi = zetahalf(khi)/ztop*(ztop-terrain(i,j))+terrain(i,j)
          zlo = z(i,j,klo)
          zhi = z(i,j,khi)

          p_at_pconst = p(i,j,1) - PCONST
          t_at_pconst = thi - (thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
          z_at_pconst = zhi - (zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

          t_surf(i,j) = t_at_pconst * (p(i,j,1)/p_at_pconst)**(GAMMA*R/G)
          t_sea_level(i,j) = t_at_pconst + GAMMA*z_at_pconst

      END DO
    END DO

    ! If we follow a traditional computation, there is a correction to the
    ! sea level temperature if both the surface and sea level
    ! temperatures are *too* hot.

    IF (ridiculous_mm5_test) THEN
      DO J = 1,ny
          DO I = 1,nx
              l1 = t_sea_level(i,j) < TC
              l2 = t_surf(i,j) <= TC
              l3 = .NOT. l1
              IF (l2 .AND. l3) THEN
                  t_sea_level(i,j) = TC
              ELSE
                  t_sea_level(i,j) = TC - 0.005D0* (t_surf(i,j)-TC)**2
              END IF
          END DO
      END DO
    END IF

    !     The grand finale: ta da!

    DO J = 1,ny
      DO I = 1,nx
    !   z_half_lowest=zetahalf(1)/ztop*(ztop-terrain(i,j))+terrain(i,j)
          z_half_lowest = z(i,j,1)

    ! Convert to hPa in this step, by multiplying by 0.01. The original
    ! Fortran routine didn't do this, but the NCL script that called it
    ! did, so we moved it here.
          sea_level_pressure(i,j) = 0.01 * (p(i,j,1)*EXP((2.D0*G*z_half_lowest)/&
                (R*(t_sea_level(i,j)+t_surf(i,j)))))
      END DO
    END DO

    !     PRINT *,'sea pres input at weird location i=20,j=1,k=1'
    !     PRINT *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
    !     PRINT *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
    !     PRINT *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
    !     PRINT *,'slp=',sea_level_pressure(20,1),
    !    *         sea_level_pressure(20,2),sea_level_pressure(20,3)

    RETURN

END SUBROUTINE f_computeslp

! Temperature from potential temperature in kelvin.
SUBROUTINE f_computetk(pressure,theta,tk,nx)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx
    REAL(KIND=8) :: pi
    REAL(KIND=8), DIMENSION(nx), INTENT(IN) :: pressure
    REAL(KIND=8), DIMENSION(nx), INTENT(IN) :: theta
    REAL(KIND=8), DIMENSION(nx), INTENT(OUT) :: tk

    INTEGER :: i
    REAL(KIND=8), PARAMETER :: P1000MB=100000.D0, R_D=287.D0, CP=7.D0*R_D/2.D0

    DO i = 1,nx
        pi = (pressure(i)/P1000MB)**(R_D/CP)
        tk(i) = pi*theta(i)
    END DO

    RETURN

END SUBROUTINE f_computetk

! Dewpoint.  Note:  1D array arguments.
SUBROUTINE f_computetd(pressure,qv_in,td,nx)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx
    REAL(KIND=8),DIMENSION(nx),INTENT(IN) :: pressure
    REAL(KIND=8),DIMENSION(nx),INTENT(IN) :: qv_in
    REAL(KIND=8),DIMENSION(nx),INTENT(OUT) :: td

    REAL(KIND=8) :: qv,tdc

    INTEGER :: i

    DO i = 1,nx
      qv = DMAX1(qv_in(i),0.D0)
      ! vapor pressure
      tdc = qv*pressure(i)/ (.622D0 + qv)

      ! avoid problems near zero
      tdc = DMAX1(tdc,0.001D0)
      td(i) = (243.5D0*LOG(tdc)-440.8D0)/ (19.48D0-LOG(tdc))
    END DO

    RETURN

END SUBROUTINE f_computetd

! Relative Humidity.  Note:  1D array arguments
SUBROUTINE f_computerh(qv,p,t,rh,nx)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx
    REAL(KIND=8),DIMENSION(nx),INTENT(IN) :: qv,p,t
    REAL(KIND=8),DIMENSION(nx),INTENT(OUT) :: rh

    REAL(KIND=8), PARAMETER :: SVP1=0.6112D0,SVP2=17.67D0,SVP3=29.65D0,SVPT0=273.15D0

    INTEGER :: i
    REAL(KIND=8) :: qvs,es,pressure,temperature
    REAL(KIND=8), PARAMETER :: R_D=287.D0,R_V=461.6D0,EP_2=R_D/R_V
    REAL(KIND=8), PARAMETER :: EP_3=0.622D0

    DO i = 1,nx
        pressure = p(i)
        temperature = t(i)
        ! es  = 1000.*svp1*
        es = 10.D0*SVP1*EXP(SVP2* (temperature-SVPT0)/(temperature-SVP3))
        ! qvs = ep_2*es/(pressure-es)
        qvs = EP_3*es/ (0.01D0*pressure- (1.D0-EP_3)*es)
        ! rh = 100*amax1(1., qv(i)/qvs)
        ! rh(i) = 100.*qv(i)/qvs
        rh(i) = 100.D0*DMAX1(DMIN1(qv(i)/qvs,1.0D0),0.0D0)
    END DO

    RETURN

END SUBROUTINE f_computerh

SUBROUTINE f_computeabsvort(u,v,msfu,msfv,msft,cor,dx,dy,av,nx,ny,nz,nxp1,nyp1)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx,ny,nz,nxp1,nyp1
    REAL(KIND=8),DIMENSION(nxp1,ny,nz),INTENT(IN) :: u
    REAL(KIND=8),DIMENSION(nx,nyp1,nz),INTENT(IN) :: v
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(OUT) :: av
    REAL(KIND=8),DIMENSION(nxp1,ny),INTENT(IN):: msfu
    REAL(KIND=8),DIMENSION(nx,nyp1),INTENT(IN) :: msfv
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: msft
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: cor
    REAL(KIND=8) :: dx,dy

    INTEGER :: jp1,jm1,ip1,im1,i,j,k
    REAL(KIND=8) :: dsy,dsx,dudy,dvdx,avort
    REAL(KIND=8) :: mm

    !          PRINT*,'nx,ny,nz,nxp1,nyp1'
    !          PRINT*,nx,ny,nz,nxp1,nyp1
    DO k = 1,nz
      DO j = 1,ny
          jp1 = MIN(j+1,ny)
          jm1 = MAX(j-1,1)
          DO i = 1,nx
              ip1 = MIN(i+1,nx)
              im1 = MAX(i-1,1)
    !         PRINT *,jp1,jm1,ip1,im1
              dsx = (ip1-im1)*dx
              dsy = (jp1-jm1)*dy
              mm = msft(i,j)*msft(i,j)
    !         PRINT *,j,i,u(i,jp1,k),msfu(i,jp1),u(i,jp1,k)/msfu(i,jp1)
              dudy = 0.5D0* (u(i,jp1,k)/msfu(i,jp1)+u(i+1,jp1,k)/&
                            msfu(i+1,jp1)-u(i,jm1,k)/&
                            msfu(i,jm1)-u(i+1,jm1,k)/&
                            msfu(i+1,jm1))/dsy*mm
              dvdx = 0.5D0* (v(ip1,j,k)/msfv(ip1,j)+v(ip1,j+1,k)/&
                            msfv(ip1,j+1)-v(im1,j,k)/&
                            msfv(im1,j)-v(im1,j+1,k)/&
                            msfv(im1,j+1))/dsx*mm
              avort = dvdx - dudy + cor(i,j)
              av(I,J,K) = avort*1.D5
          END DO
      END DO
    END DO

    RETURN

END SUBROUTINE f_computeabsvort


SUBROUTINE f_computepvo(u,v,theta,prs,msfu,msfv,msft,cor,dx,dy,pv,nx,ny,nz,nxp1,nyp1)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nx,ny,nz,nxp1,nyp1
    REAL(KIND=8),DIMENSION(nxp1,ny,nz),INTENT(IN) :: u
    REAL(KIND=8),DIMENSION(nx,nyp1,nz),INTENT(IN) :: v
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: prs
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: theta
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(OUT) :: pv
    REAL(KIND=8),DIMENSION(nxp1,ny),INTENT(IN) ::  msfu
    REAL(KIND=8),DIMENSION(nx,nyp1),INTENT(IN) :: msfv
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: msft
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: cor
    REAL(KIND=8) :: dx,dy


    INTEGER :: kp1,km1,jp1,jm1,ip1,im1,i,j,k
    REAL(KIND=8) :: dsy,dsx,dp,dudy,dvdx,dudp,dvdp,dthdp,avort
    REAL(KIND=8) :: dthdx,dthdy,mm

    !          PRINT*,'nx,ny,nz,nxp1,nyp1'
    !          PRINT*,nx,ny,nz,nxp1,nyp1
    DO k = 1,nz
      kp1 = MIN(k+1,nz)
      km1 = MAX(k-1,1)
      DO J = 1,ny
          jp1 = MIN(j+1,ny)
          jm1 = MAX(j-1,1)
          DO i = 1,nx
              ip1 = MIN(i+1,nx)
              im1 = MAX(i-1,1)
    !         PRINT *,jp1,jm1,ip1,im1
              dsx = (ip1-im1)*dx
              dsy = (jp1-jm1)*dy
              mm = msft(i,j)*msft(i,j)
    !         PRINT *,j,i,u(i,jp1,k),msfu(i,jp1),u(i,jp1,k)/msfu(i,jp1)
              dudy = 0.5D0* (u(i,jp1,k)/msfu(i,jp1)+u(i+1,jp1,k)/&
                            msfu(i+1,jp1)-u(i,jm1,k)/&
                            msfu(i,jm1)-u(i+1,jm1,k)/&
                            msfu(i+1,jm1))/dsy*mm
              dvdx = 0.5D0* (v(ip1,j,k)/msfv(ip1,j)+v(ip1,j+1,k)/&
                            msfv(ip1,j+1)-v(im1,j,k)/&
                            msfv(im1,j)-v(im1,j+1,k)/&
                            msfv(im1,j+1))/dsx*mm
              avort = dvdx - dudy + cor(i,j)
              dp = prs(i,j,kp1) - prs(i,j,km1)
              dudp = 0.5D0* (u(i,j,kp1)+u(i+1,j,kp1)-u(i,j,km1)-u(i+1,j,km1))/dp
              dvdp = 0.5D0* (v(i,j,kp1)+v(i,j+1,kp1)-v(i,j,km1)-v(i,J+1,km1))/dp
              dthdp = (theta(i,j,kp1)-theta(i,j,km1))/dp
              dthdx = (theta(ip1,j,k)-theta(im1,j,k))/dsx*msft(i,j)
              dthdy = (theta(i,jp1,k)-theta(i,jm1,k))/dsy*msft(i,j)
              pv(i,j,k) = -9.81D0* (dthdp*avort-dvdp*dthdx+dudp*dthdy)*10000.D0
    !               if(i.eq.300 .and. j.eq.300) then
    !                 PRINT*,'avort,dudp,dvdp,dthdp,dthdx,dthdy,pv'
    !                 PRINT*,avort,dudp,dvdp,dthdp,dthdx,dthdy,pv(i,j,k)
    !               endif
              pv(i,j,k) = pv(i,j,k)*1.D2
          END DO
      END DO
    END DO

    RETURN

END SUBROUTINE f_computepvo

! Theta-e
SUBROUTINE f_computeeth(qvp,tmk,prs,eth,miy,mjx,mkzh)

    IMPLICIT NONE

    ! Input variables
    ! Qvapor [g/kg]
    REAL(KIND=8),DIMENSION(miy,mjx,mkzh),INTENT(IN) :: qvp
    ! Temperature [K]
    REAL(KIND=8),DIMENSION(miy,mjx,mkzh),INTENT(IN) :: tmk
    ! full pressure (=P+PB) [hPa]
    REAL(KIND=8),DIMENSION(miy,mjx,mkzh),INTENT(IN) :: prs
    !
    ! Output variable
    ! equivalent potential temperature [K]
    REAL(KIND=8),DIMENSION(miy,mjx,mkzh),INTENT(OUT) :: eth
    ! Sizes
    INTEGER,INTENT(IN) :: miy,mjx,mkzh

    ! local variables
    REAL(KIND=8) :: q
    REAL(KIND=8) :: t
    REAL(KIND=8) :: p
    REAL(KIND=8) :: e
    REAL(KIND=8) :: tlcl
    INTEGER :: i,j,k

    ! parameters
    REAL(KIND=8),PARAMETER :: EPS = 0.622D0
    REAL(KIND=8),PARAMETER :: RGAS = 287.04D0
    REAL(KIND=8),PARAMETER :: RGASMD = .608D0
    REAL(KIND=8),PARAMETER :: CP = 1004.D0
    REAL(KIND=8),PARAMETER :: CPMD = .887D0
    REAL(KIND=8),PARAMETER :: GAMMA = RGAS/CP
    REAL(KIND=8),PARAMETER :: GAMMAMD = RGASMD - CPMD
    REAL(KIND=8),PARAMETER :: TLCLC1 = 2840.D0
    REAL(KIND=8),PARAMETER :: TLCLC2 = 3.5D0
    REAL(KIND=8),PARAMETER :: TLCLC3 = 4.805D0
    REAL(KIND=8),PARAMETER :: TLCLC4 = 55.D0
    REAL(KIND=8),PARAMETER :: THTECON1 = 3376.D0
    REAL(KIND=8),PARAMETER :: THTECON2 = 2.54D0
    REAL(KIND=8),PARAMETER :: THTECON3 = .81D0

    DO k = 1,mkzh
      DO j = 1,mjx
          DO i = 1,miy
              q = MAX(qvp(i,j,k),1.D-15)
              t = tmk(i,j,k)
              p = prs(i,j,k)/100.
              e = q*p / (EPS+q)
              tlcl = TLCLC1 / (LOG(t**TLCLC2/e)-TLCLC3) + TLCLC4
              eth(i,j,k) = t * (1000.D0/p)**(GAMMA * (1.D0+GAMMAMD*q))*&
                        EXP((THTECON1/tlcl-THTECON2)*q*(1.D0+THTECON3*q))
          END DO
      END DO
    END DO

    RETURN

END SUBROUTINE f_computeeth

SUBROUTINE f_computeuvmet(u,v,longca,longcb,flong,flat,&
                        cen_long,cone,rpd,istag,is_msg_val,umsg,vmsg,uvmetmsg,&
                        uvmet,nx,ny,nxp1,nyp1,nz)

    IMPLICIT NONE

    ! ISTAG should be 0 if the U,V grids are not staggered.
    ! That is, NY = NYP1 and NX = NXP1.

    INTEGER,INTENT(IN) :: nx,ny,nz,nxp1,nyp1,istag
    LOGICAL,INTENT(IN) :: is_msg_val
    REAL(KIND=8),DIMENSION(nxp1,ny,nz),INTENT(IN):: u
    REAL(KIND=8),DIMENSION(nx,nyp1,nz),INTENT(IN) :: v
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: flong
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: flat
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(INOUT) :: longca
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(INOUT) :: longcb
    REAL(KIND=8),INTENT(IN) :: cen_long,cone,rpd
    REAL(KIND=8),INTENT(IN) :: umsg,vmsg,uvmetmsg
    REAL(KIND=8),DIMENSION(nx,ny,nz,2),INTENT(OUT) :: uvmet

    INTEGER :: i,j,k
    REAL(KIND=8) :: uk,vk

    ! msg stands for missing value in this code
    !      WRITE (6,FMT=*) ' in compute_uvmet ',NX,NY,NXP1,NYP1,ISTAG

    DO J = 1,ny
      DO I = 1,nx

          longca(i,j) = flong(i,j) - cen_long
          IF (longca(i,j).GT.180.D0) THEN
              longca(i,j) = longca(i,j) - 360.D0
          END IF
          IF (longca(i,j).LT.-180.D0) THEN
              longca(i,j) = longca(i,j) + 360.D0
          END IF
          IF (flat(i,j).LT.0.D0) THEN
              longcb(i,j) = -longca(i,j)*cone*rpd
          ELSE
              longcb(i,j) = longca(i,j)*cone*rpd
          END IF

          longca(i,j) = COS(longcb(i,j))
          longcb(i,j) = SIN(longcb(i,j))

      END DO
    END DO

    !      WRITE (6,FMT=*) ' computing velocities '
    DO k = 1,nz
        DO j = 1,ny
         DO i = 1,nx
            IF (istag.EQ.1) THEN
               IF (is_msg_val .AND. (u(i,j,k) .EQ. umsg .OR. v(i,j,k) .EQ. vmsg &
                   .OR. u(i+1,j,k) .EQ. umsg .OR. v(i,j+1,k) .EQ. vmsg)) THEN
                  uvmet(i,j,k,1) = uvmetmsg
                  uvmet(i,j,k,2) = uvmetmsg
               ELSE
                  uk = 0.5D0* (u(i,j,k)+u(i+1,j,k))
                  vk = 0.5D0* (v(i,j,k)+v(i,j+1,k))
                  uvmet(i,j,k,1) = vk*longcb(i,j) + uk*longca(i,j)
                  uvmet(i,j,k,2) = vk*longca(i,j) - uk*longcb(i,j)
               END IF
            ELSE
               IF (is_msg_val .AND. (u(i,j,k) .EQ. umsg .OR. v(i,j,k) .EQ. vmsg)) THEN
                  uvmet(i,j,k,1) = uvmetmsg
                  uvmet(i,j,k,2) = uvmetmsg
               ELSE
                  uk = u(i,j,k)
                  vk = v(i,j,k)
                  uvmet(i,j,k,1) = vk*longcb(i,j) + uk*longca(i,j)
                  uvmet(i,j,k,2) = vk*longca(i,j) - uk*longcb(i,j)
               END IF
            END IF
         END DO
        END DO
    END DO

    RETURN

END SUBROUTINE f_computeuvmet


SUBROUTINE f_computeomega(qvp,tmk,www,prs,omg,mx,my,mz)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: mx, my, mz
    REAL(KIND=8),INTENT(IN),DIMENSION(mx,my,mz) :: qvp
    REAL(KIND=8),INTENT(IN),DIMENSION(mx,my,mz) :: tmk
    REAL(KIND=8),INTENT(IN),DIMENSION(mx,my,mz) :: www
    REAL(KIND=8),INTENT(IN),DIMENSION(mx,my,mz) :: prs
    REAL(KIND=8),INTENT(OUT),DIMENSION(mx,my,mz) :: omg

    ! Local variables
    INTEGER :: i, j, k
    REAL(KIND=8),PARAMETER :: GRAV=9.81, RGAS=287.04, EPS=0.622

    DO k=1,mz
        DO j=1,my
          DO i=1,mx
              omg(i,j,k)=-GRAV*prs(i,j,k) / &
              (RGAS*((tmk(i,j,k)*(EPS+qvp(i,j,k))) / &
              (EPS*(1.+qvp(i,j,k)))))*www(i,j,k)
          END DO
        END DO
    END DO
    !
    RETURN

END SUBROUTINE f_computeomega

SUBROUTINE f_computetv(temp,qv,tv,nx,ny,nz)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx,ny,nz
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: temp
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: qv
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(OUT) :: tv

    INTEGER :: i,j,k
    REAL(KIND=8),PARAMETER :: EPS = 0.622D0

    DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tv(i,j,k) = temp(i,j,k) * (EPS+qv(i,j,k)) / (EPS * (1.D0+qv(i,j,k)))
          END DO
        END DO
    END DO

    RETURN

END SUBROUTINE f_computetv

! Need to deal with the fortran stop statements
SUBROUTINE f_computewetbulb(prs,tmk,qvp,PSADITHTE,PSADIPRS,PSADITMK,throw_exception,twb,nx,ny,nz)

    IMPLICIT NONE

    EXTERNAL throw_exception
    INTEGER,INTENT(IN) :: nx, ny, nz
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: prs
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: tmk
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: qvp
    REAL(KIND=8),DIMENSION(150),INTENT(IN) :: PSADITHTE
    REAL(KIND=8),DIMENSION(150),INTENT(IN) ::PSADIPRS
    REAL(KIND=8),DIMENSION(150,150),INTENT(IN) :: PSADITMK
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(OUT) :: twb
    !EXTERNAL throw_exception

    INTEGER :: i,j,k
    INTEGER :: jtch,jt,ipch,ip
    REAL(KIND=8) :: q, t, p, e, tlcl, eth
    REAL(KIND=8) :: fracip,fracip2,fracjt,fracjt2
    REAL(KIND=8) :: tonpsadiabat

    REAL(KIND=8),PARAMETER :: EPS=0.622
    REAL(KIND=8),PARAMETER :: RGAS=287.04
    REAL(KIND=8),PARAMETER :: RGASMD=.608
    REAL(KIND=8),PARAMETER :: CP=1004.
    REAL(KIND=8),PARAMETER :: CPMD=.887
    REAL(KIND=8),PARAMETER :: GAMMA=RGAS/CP
    REAL(KIND=8),PARAMETER :: GAMMAMD=RGASMD-CPMD
    REAL(KIND=8),PARAMETER :: CELKEL=273.15
    REAL(KIND=8),PARAMETER :: TLCLC1=2840.
    REAL(KIND=8),PARAMETER :: TLCLC2=3.5
    REAL(KIND=8),PARAMETER :: TLCLC3=4.805
    REAL(KIND=8),PARAMETER :: TLCLC4=55.
    REAL(KIND=8),PARAMETER :: THTECON1=3376.
    REAL(KIND=8),PARAMETER :: THTECON2=2.54
    REAL(KIND=8),PARAMETER :: THTECON3=.81

    DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            q=dmax1(qvp(i,j,k),1.d-15)
            t=tmk(i,j,k)
            p=prs(i,j,k)/100.
            e=q*p/(EPS+q)
            tlcl=TLCLC1/(DLOG(t**TLCLC2/e)-TLCLC3)+TLCLC4
            eth=t*(1000./p)**(GAMMA*(1.+GAMMAMD*q))*&
                EXP((THTECON1/tlcl-THTECON2)*q*(1.+THTECON3*q))

        !   Now we need to find the temperature (in K) on a moist adiabat
        !   (specified by eth in K) given pressure in hPa.  It uses a
        !   lookup table, with data that was generated by the Bolton (1980)
        !   formula for theta_e.

        !     First check if pressure is less than min pressure in lookup table.
        !     If it is, assume parcel is so dry that the given theta-e value can
        !     be interpretted as theta, and get temperature from the simple dry
        !     theta formula.
        !

            IF (p.LE.psadiprs(150)) THEN
              tonpsadiabat=eth*(p/1000.)**GAMMA
            ELSE

        !   Otherwise, look for the given thte/prs point in the lookup table.

            jt=-1
            DO jtch=1,150-1
              IF (eth.GE.psadithte(jtch).AND.eth.LT.psadithte(jtch+1)) THEN
                 jt=jtch
                 EXIT
              ENDIF
            END DO
    !        jt=-1
    !    213        CONTINUE
            ip=-1
            DO ipch=1,150-1
              IF (p.LE.psadiprs(ipch).AND.p.GT.psadiprs(ipch+1)) THEN
                 ip=ipch
                 EXIT
              ENDIF
            END DO
    !        ip=-1
    !    215        CONTINUE
            IF (jt.EQ.-1.OR.ip.EQ.-1) THEN
                CALL throw_exception('Outside of lookup table bounds. prs,thte=',p,eth)
              !STOP ! TODO: Need to make python throw an exception here
            ENDIF
            fracjt=(eth-psadithte(jt))/(psadithte(jt+1)-psadithte(jt))
            fracjt2=1.-fracjt
            fracip=(psadiprs(ip)-p)/(psadiprs(ip)-psadiprs(ip+1))
            fracip2=1.-fracip
            IF (psaditmk(ip,jt).GT.1e9.OR.psaditmk(ip+1,jt).GT.1e9.OR.psaditmk(ip,jt+1).GT.1e9.OR.psaditmk(ip+1,jt+1).GT.1e9) THEN
                !PRINT*,'Tried to access missing tmperature in lookup table.'
                CALL throw_exception('Prs and Thte probably unreasonable. prs,thte=',p,eth)
               !STOP ! TODO: Need to make python throw an exception here
            ENDIF
            tonpsadiabat=fracip2*fracjt2*psaditmk(ip,jt)+fracip*&
                    fracjt2*psaditmk(ip+1,jt)+fracip2*fracjt*psaditmk(ip,jt+1)+fracip*&
                    fracjt *psaditmk(ip+1,jt+1)
            ENDIF

            twb(i,j,k)=tonpsadiabat

          END DO
        END DO
    END DO

    RETURN

END SUBROUTINE f_computewetbulb

SUBROUTINE f_computesrh(u, v, ght, ter, top, sreh, miy, mjx, mkzh)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: miy, mjx, mkzh
    REAL(KIND=8),DIMENSION(miy,mjx,mkzh),INTENT(IN) :: u, v, ght
    REAL(KIND=8),INTENT(IN) :: top
    REAL(KIND=8),DIMENSION(miy,mjx),INTENT(IN) :: ter
    REAL(KIND=8),DIMENSION(miy,mjx),INTENT(OUT) :: sreh

    ! This helicity code was provided by Dr. Craig Mattocks, and
    ! verified by Cindy Bruyere to produce results equivalent to
    ! those generated by RIP4. (The code came from RIP4?)

    REAL(KIND=8) :: dh, sdh, su, sv, ua, va, asp, adr, bsp, bdr
    REAL(KIND=8) :: cu, cv, x, sum
    INTEGER :: i, j, k, k10, k3, ktop
    REAL(KIND=8),PARAMETER :: PI=3.14159265d0, DTR=PI/180.d0, DPR=180.d0/PI

    DO j = 1, mjx-1
        DO i = 1, miy-1
            sdh = 0.d0
            su = 0.d0
            sv = 0.d0
            k3 = 0
            k10 = 0
            ktop = 0
            DO k = mkzh, 2, -1
                IF (((ght(i,j,k) - ter(i,j)) .GT. 10000.d0) .AND.(k10 .EQ. 0)) THEN
                    k10 = k
                    EXIT
                ENDIF
                IF (((ght(i,j,k) - ter(i,j)) .GT. top) .AND.(ktop .EQ. 0)) THEN
                    ktop = k
                ENDIF
                IF (((ght(i,j,k) - ter(i,j)) .GT. 3000.d0) .AND.(k3 .EQ. 0)) THEN
                    k3 = k
                ENDIF
            END DO

            IF (k10 .EQ. 0) THEN
                k10 = 2
            ENDIF
            DO k = k3, k10, -1
                dh = ght(i,j,k-1) - ght(i,j,k)
                sdh = sdh + dh
                su = su + 0.5d0*dh*(u(i,j,k-1)+u(i,j,k))
                sv = sv + 0.5d0*dh*(v(i,j,k-1)+v(i,j,k))
            END DO
            ua = su / sdh
            va = sv / sdh
            asp = SQRT(ua*ua + va*va)
            IF (ua .EQ. 0.d0 .AND. va .EQ. 0.d0) THEN
                adr = 0.d0
            ELSE
                adr = DPR * (PI + ATAN2(ua,va))
            ENDIF
            bsp = 0.75d0 * asp
            bdr = adr + 30.d0
            IF (bdr .GT. 360.d0) THEN
                bdr = bdr-360.d0
            ENDIF
            cu = -bsp * SIN(bdr*dtr)
            cv = -bsp * COS(bdr*dtr)
            sum = 0.d0
            DO k = mkzh-1, ktop, -1
                x = ((u(i,j,k)-cu) * (v(i,j,k)-v(i,j,k+1))) - ((v(i,j,k)-cv) * (u(i,j,k)-u(i,j,k+1)))
                sum = sum + x
            END DO
            sreh(i,j) = -sum
        END DO
    END DO

    RETURN

END SUBROUTINE f_computesrh

SUBROUTINE f_computeuh(zp,mapfct,dx,dy,uhmnhgt,uhmxhgt,us,vs,w,tem1,tem2,uh,nx,ny,nz,nzp1)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz,nzp1
  REAL(KIND=8),DIMENSION(nx,ny,nzp1),INTENT(IN)  :: zp
  REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: mapfct
  REAL(KIND=8),INTENT(IN)  :: dx,dy
  REAL(KIND=8),INTENT(IN)  :: uhmnhgt,uhmxhgt
  REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN)  :: us
  REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN)  :: vs
  REAL(KIND=8),DIMENSION(nx,ny,nzp1),INTENT(IN)  :: w
  REAL(KIND=8),DIMENSION(nx,ny),INTENT(OUT) :: uh
  REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: tem1
  REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: tem2
!
! Misc local variables
!
  INTEGER :: i,j,k,kbot,ktop
  REAL(KIND=8) :: twodx,twody,wgtlw,sum,wmean,wsum,wavg
  REAL(KIND=8) :: helbot,heltop,wbot,wtop
  REAL(KIND=8) :: zbot,ztop
!
! Initialize arrays
!
  uh=0.0
  tem1=0.0
!
! Calculate vertical component of helicity at scalar points
!   us: u at scalar points
!   vs: v at scalar points
!
  twodx=2.0*dx
  twody=2.0*dy
  DO k=2,nz-2
    DO j=2,ny-1
      DO i=2,nx-1
        wavg=0.5*(w(i,j,k)+w(i,j,k+1))
        tem1(i,j,k)=wavg *                                      &
            ((vs(i+1,j,k)-vs(i-1,j,k))/(twodx*mapfct(i,j))  -   &
             (us(i,j+1,k)-us(i,j-1,k))/(twody*mapfct(i,j)))
        tem2(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO
!
! Integrate over depth uhminhgt to uhmxhgt AGL
!
!  WRITE(6,'(a,f12.1,a,f12.1,a)') &
!        'Calculating UH from ',uhmnhgt,' to ',uhmxhgt,' m AGL'
  DO j=2,ny-2
    DO i=2,nx-2
      zbot=zp(i,j,2)+uhmnhgt
      ztop=zp(i,j,2)+uhmxhgt
!
! Find wbar, weighted-mean vertical velocity in column
! Find w at uhmnhgt AGL (bottom)
!
      DO k=2,nz-3
        IF(zp(i,j,k) > zbot) EXIT
      END DO
      kbot=k
      wgtlw=(zp(i,j,kbot)-zbot)/(zp(i,j,kbot)-zp(i,j,kbot-1))
      wbot=(wgtlw*w(i,j,kbot-1))+((1.-wgtlw)*w(i,j,kbot))
!
! Find w at uhmxhgt AGL (top)
!
      DO k=2,nz-3
        IF(zp(i,j,k) > ztop) EXIT
      END DO
      ktop=k
      wgtlw=(zp(i,j,ktop)-ztop)/(zp(i,j,ktop)-zp(i,j,ktop-1))
      wtop=(wgtlw*w(i,j,ktop-1))+((1.-wgtlw)*w(i,j,ktop))
!
! First part, uhmnhgt to kbot
!
      wsum=0.5*(w(i,j,kbot)+wbot)*(zp(i,j,kbot)-zbot)
!
! Integrate up through column
!
      DO k=(kbot+1),(ktop-1)
        wsum=wsum+0.5*(w(i,j,k)+w(i,j,k-1))*(zp(i,j,k)-zp(i,j,k-1))
      END DO
!
! Last part, ktop-1 to uhmxhgt
!
      wsum=wsum+0.5*(wtop+w(i,j,ktop-1))*(ztop-zp(i,j,ktop-1))
      wmean=wsum/(uhmxhgt-uhmnhgt)

      IF(wmean > 0.) THEN    ! column updraft, not downdraft
!
! Find helicity at uhmnhgt AGL (bottom)
!
        DO k=2,nz-3
          IF(tem2(i,j,k) > zbot) EXIT
        END DO
        kbot=k
        wgtlw=(tem2(i,j,kbot)-zbot)/(tem2(i,j,kbot)-tem2(i,j,kbot-1))
        helbot=(wgtlw*tem1(i,j,kbot-1))+((1.-wgtlw)*tem1(i,j,kbot))
!
! Find helicity at uhmxhgt AGL (top)
!
        DO k=2,nz-3
          IF(tem2(i,j,k) > ztop) EXIT
        END DO
        ktop=k
        wgtlw=(tem2(i,j,ktop)-ztop)/(tem2(i,j,ktop)-tem2(i,j,ktop-1))
        heltop=(wgtlw*tem1(i,j,ktop-1))+((1.-wgtlw)*tem1(i,j,ktop))
!
! First part, uhmnhgt to kbot
!
        sum=0.5*(tem1(i,j,kbot)+helbot)*(tem2(i,j,kbot)-zbot)
!
! Integrate up through column
!
        DO k=(kbot+1),(ktop-1)
          sum=sum+0.5*(tem1(i,j,k)+tem1(i,j,k-1))*(tem2(i,j,k)-tem2(i,j,k-1))
        END DO
!
! Last part, ktop-1 to uhmxhgt
!
        uh(i,j)=sum+0.5*(heltop+tem1(i,j,ktop-1))*(ztop-tem2(i,j,ktop-1))
      END IF
    END DO
  END DO

  uh = uh * 1000.   ! Scale according to Kain et al. (2008)

  RETURN

END SUBROUTINE f_computeuh

SUBROUTINE f_computepw(p,tv,qv,ht,zdiff,pw,nx,ny,nz,nzh)

    IMPLICIT NONE

    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: p,tv,qv
    REAL(KIND=8),DIMENSION(nx,ny,nzh),INTENT(IN) :: ht
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(OUT) :: pw
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(INOUT) :: zdiff
    INTEGER,INTENT(IN) :: nx,ny,nz,nzh

    INTEGER :: i,j,k
    REAL(KIND=8),PARAMETER :: R=287.06

    pw = 0
    DO k=1,nz
        DO j=1,ny
            DO i=1,nx
                zdiff(i,j) = ht(i,j,k+1) - ht(i,j,k)
                pw(i,j) = pw(i,j) + ((p(i,j,k)/(R * tv(i,j,k))) * qv(i,j,k) * zdiff(i,j))
            END DO
        END DO
    END DO


    RETURN

END SUBROUTINE f_computepw


SUBROUTINE f_computedbz(prs,tmk,qvp,qra,qsn,qgr,sn0,ivarint,iliqskin,dbz,nx,ny,nz)

    IMPLICIT NONE

    !   Arguments
    INTEGER, INTENT(IN) :: nx,ny,nz
    INTEGER, INTENT(IN) :: sn0,ivarint,iliqskin
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(OUT) :: dbz
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: prs
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: tmk
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: qvp
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: qra
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: qsn
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: qgr

    !   Local Variables
    INTEGER :: i,j,k
    REAL(KIND=8) :: temp_c,virtual_t
    REAL(KIND=8) :: gonv,ronv,sonv
    REAL(KIND=8) :: factor_g,factor_r,factor_s
    REAL(KIND=8) :: factorb_g,factorb_s
    REAL(KIND=8) :: rhoair,z_e

    !   Constants used to calculate variable intercepts
    REAL(KIND=8),PARAMETER :: R1 = 1.D-15
    REAL(KIND=8),PARAMETER :: RON = 8.D6
    REAL(KIND=8),PARAMETER :: RON2 = 1.D10
    REAL(KIND=8),PARAMETER :: SON = 2.D7
    REAL(KIND=8),PARAMETER :: GON = 5.D7
    REAL(KIND=8),PARAMETER :: RON_MIN = 8.D6
    REAL(KIND=8),PARAMETER :: RON_QR0 = 0.00010D0
    REAL(KIND=8),PARAMETER :: RON_DELQR0 = 0.25D0*RON_QR0
    REAL(KIND=8),PARAMETER :: RON_CONST1R = (RON2-RON_MIN)*0.5D0
    REAL(KIND=8),PARAMETER :: RON_CONST2R = (RON2+RON_MIN)*0.5D0

    !   Constant intercepts
    REAL(KIND=8),PARAMETER :: RN0_R = 8.D6
    REAL(KIND=8),PARAMETER :: RN0_S = 2.D7
    REAL(KIND=8),PARAMETER :: RN0_G = 4.D6

    !   Other constants
    REAL(KIND=8),PARAMETER :: GAMMA_SEVEN = 720.D0
    REAL(KIND=8),PARAMETER :: RHOWAT = 1000.D0
    REAL(KIND=8),PARAMETER :: RHO_R = RHOWAT
    REAL(KIND=8),PARAMETER :: RHO_S = 100.D0
    REAL(KIND=8),PARAMETER :: RHO_G = 400.D0
    REAL(KIND=8),PARAMETER :: ALPHA = 0.224D0
    REAL(KIND=8),PARAMETER :: CELKEL = 273.15D0
    REAL(KIND=8),PARAMETER :: PI = 3.141592653589793D0
    REAL(KIND=8),PARAMETER :: RD = 287.04D0

    !   Force all Q arrays to be 0.0 or greater.
    DO k = 1,nz
     DO j = 1,ny
        DO i = 1,nx
           IF (qvp(i,j,k).LT.0.0) THEN
              qvp(i,j,k) = 0.0
           END IF
           IF (qra(i,j,k).LT.0.0) THEN
              qra(i,j,k) = 0.0
           END IF
           IF (qsn(i,j,k).LT.0.0) THEN
              qsn(i,j,k) = 0.0
           END IF
           IF (qgr(i,j,k).LT.0.0) THEN
              qgr(i,j,k) = 0.0
           END IF
        END DO
     END DO
    END DO

    !   Input pressure is Pa, but we need hPa in calculations

    IF (sn0.EQ.0) THEN
      DO k = 1,nz
          DO j = 1,ny
              DO i = 1,nx
                  IF (tmk(i,j,k).LT.CELKEL) THEN
                      qsn(i,j,k) = qra(i,j,k)
                      qra(i,j,k) = 0.D0
                  END IF
              END DO
          END DO
      END DO
    END IF


    factor_r = GAMMA_SEVEN*1.D18* (1.D0/ (PI*RHO_R))**1.75D0
    factor_s = GAMMA_SEVEN*1.D18* (1.D0/ (PI*RHO_S))**1.75D0*(RHO_S/RHOWAT)**2*ALPHA
    factor_g = GAMMA_SEVEN*1.D18* (1.D0/ (PI*RHO_G))**1.75D0*(RHO_G/RHOWAT)**2*ALPHA


    DO k = 1,nz
      DO j = 1,ny
          DO i = 1,nx

              virtual_t = tmk(i,j,k)* (0.622D0+qvp(i,j,k))/(0.622D0* (1.D0+qvp(i,j,k)))
              rhoair = prs(i,j,k) / (RD*virtual_t)

    !      Adjust factor for brightband, where snow or graupel particle
    !      scatters like liquid water (alpha=1.0) because it is assumed to
    !      have a liquid skin.

              IF (iliqskin.EQ.1 .AND. tmk(i,j,k).GT.CELKEL) THEN
                  factorb_s = factor_s/ALPHA
                  factorb_g = factor_g/ALPHA
              ELSE
                  factorb_s = factor_s
                  factorb_g = factor_g
              END IF

    !      Calculate variable intercept parameters

              IF (ivarint.EQ.1) THEN

                  temp_c = DMIN1(-0.001D0,tmk(i,j,k)-CELKEL)
                  sonv = DMIN1(2.0D8,2.0D6*EXP(-0.12D0*temp_c))

                  gonv = gon
                  IF (qgr(i,j,k).GT.R1) THEN
                      gonv = 2.38D0* (PI*RHO_G/(rhoair*qgr(i,j,k)))**0.92D0
                      gonv = MAX(1.D4,MIN(gonv,GON))
                  END IF

                  ronv = RON2
                  IF (qra(i,j,k).GT.R1) THEN
                      ronv = RON_CONST1R*TANH((RON_QR0-qra(i,j,k))/RON_DELQR0) + RON_CONST2R
                  END IF

              ELSE

                  ronv = RN0_R
                  sonv = RN0_S
                  gonv = RN0_G

              END IF

    !      Total equivalent reflectivity factor (z_e, in mm^6 m^-3) is
    !      the sum of z_e for each hydrometeor species:

              z_e = factor_r*(rhoair*qra(i,j,k))**1.75D0/ronv**.75D0 + &
                factorb_s*(rhoair*qsn(i,j,k))**1.75D0/sonv**.75D0 + &
                factorb_g* (rhoair*qgr(i,j,k))**1.75D0/gonv**.75D0

    !      Adjust small values of Z_e so that dBZ is no lower than -30
              z_e = MAX(z_e,.001D0)

    !      Convert to dBZ
              dbz(i,j,k) = 10.D0*LOG10(z_e)

          END DO
      END DO
    END DO

    RETURN

END SUBROUTINE f_computedbz

SUBROUTINE rotatecoords(ilat,ilon,olat,olon,lat_np,lon_np,lon_0,direction)

    IMPLICIT NONE

    REAL(KIND=8),INTENT(IN) :: ilat,ilon
    REAL(KIND=8),INTENT(OUT) :: olat,olon
    REAL(KIND=8),INTENT(IN) :: lat_np,lon_np,lon_0
    INTEGER,INTENT(IN) :: direction

    !  >=0, default : computational -> geographical
    !  < 0          : geographical  -> computational

    REAL(KIND=8) :: rlat,rlon
    REAL(KIND=8) :: phi_np,lam_np,lam_0,dlam
    REAL(KIND=8) :: sinphi,cosphi,coslam,sinlam
    REAL(KIND=8),PARAMETER :: PI=3.141592653589793d0
    REAL(KIND=8),PARAMETER :: RAD_PER_DEG=PI/180.d0
    REAL(KIND=8),PARAMETER :: DEG_PER_RAD=180.d0/PI

    !convert all angles to radians
    phi_np = lat_np*RAD_PER_DEG
    lam_np = lon_np*RAD_PER_DEG
    lam_0 = lon_0*RAD_PER_DEG
    rlat = ilat*RAD_PER_DEG
    rlon = ilon*RAD_PER_DEG

    IF (direction.LT.0) THEN
    ! the equations are exactly the same except for one
    ! small difference with respect to longitude ...
      dlam = pi - lam_0
    ELSE
      dlam = lam_np
    END IF

    sinphi = COS(phi_np)*COS(rlat)*COS(rlon-dlam) + SIN(phi_np)*SIN(rlat)
    cosphi = SQRT(1.d0-sinphi*sinphi)
    coslam = SIN(phi_np)*COS(rlat)*COS(rlon-dlam) - COS(phi_np)*SIN(rlat)
    sinlam = COS(rlat)*SIN(rlon-dlam)

    IF (cosphi.NE.0.d0) THEN
      coslam = coslam/cosphi
      sinlam = sinlam/cosphi
    END IF

    olat = DEG_PER_RAD*ASIN(sinphi)
    olon = DEG_PER_RAD*(ATAN2(sinlam,coslam)-dlam-lam_0+lam_np)

    RETURN

END SUBROUTINE rotatecoords

SUBROUTINE f_lltoij(map_proj,truelat1,truelat2,stdlon,lat1,lon1,&
                   pole_lat,pole_lon,knowni,knownj,dx,latinc,&
                   loninc,lat,lon,throw_exception,loc)

    IMPLICIT NONE

    EXTERNAL throw_exception
    ! Converts input lat/lon values to the cartesian (i,j) value
    ! for the given projection.

    INTEGER,INTENT(IN) :: map_proj
    REAL(KIND=8),INTENT(IN) :: stdlon
    REAL(KIND=8),INTENT(IN) ::lat1,lon1,pole_lat,pole_lon,knowni,knownj
    REAL(KIND=8),INTENT(IN) ::dx,latinc,loninc
    REAL(KIND=8),INTENT(INOUT) :: lat,lon,truelat1,truelat2 ! these might get modified
    REAL(KIND=8),DIMENSION(2),INTENT(OUT) :: loc

    REAL(KIND=8) :: deltalon1
    REAL(KIND=8) :: tl1r
    REAL(KIND=8) :: clain,dlon,rsw,deltalon,deltalat
    REAL(KIND=8) :: reflon,scale_top,ala1,alo1,ala,alo,rm,polei,polej
    ! earth radius divided by dx
    REAL(KIND=8) :: rebydx
    REAL(KIND=8) :: ctl1r,arg,cone,hemi
    REAL(KIND=8) :: i,j
    REAL(KIND=8) :: lat1n,lon1n,olat,olon

    ! Contants
    REAL(KIND=8),PARAMETER :: PI=3.141592653589793d0
    REAL(KIND=8),PARAMETER :: RAD_PER_DEG=PI/180.d0
    REAL(KIND=8),PARAMETER :: DEG_PER_RAD=180.d0/PI
    REAL(KIND=8),PARAMETER :: RE_M=6370000.d0 ! radius of earth in meters

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

    rebydx = RE_M/dx

    ! Get rid of compiler warnings
    i=0
    j=0

    hemi = 1.0d0
    IF (truelat1.LT.0.0d0) THEN
      hemi = -1.0d0
    END IF

    ! mercator
    IF (map_proj.EQ.3) THEN

      ! preliminary variables
      clain = COS(RAD_PER_DEG*truelat1)
      dlon = dx/(RE_M*clain)

      ! compute distance from equator to origin, and store in
      ! the rsw tag.
      rsw = 0.d0
      IF (lat1.NE.0.d0) THEN
          rsw = (DLOG(TAN(0.5d0*((lat1+90.d0)*RAD_PER_DEG))))/dlon
      END IF

      deltalon = lon - lon1
      IF (deltalon.LT.-180.d0) deltalon = deltalon + 360.d0
      IF (deltalon.GT.180.d0) deltalon = deltalon - 360.d0
      i = knowni + (deltalon/(dlon*DEG_PER_RAD))
      j = knownj + (DLOG(TAN(0.5d0*((lat+90.d0)*RAD_PER_DEG))))/dlon - rsw

    ! ps
    ELSE IF (map_proj.EQ.2) THEN

      reflon = stdlon + 90.d0

      ! compute numerator term of map scale factor
      scale_top = 1.d0 + hemi*SIN(truelat1*RAD_PER_DEG)

      ! compute radius to lower-left (sw) corner
      ala1 = lat1*RAD_PER_DEG
      rsw = rebydx*COS(ala1)*scale_top/(1.d0+hemi*SIN(ala1))

      ! find the pole point
      alo1 = (lon1-reflon)*RAD_PER_DEG
      polei = knowni - rsw*COS(alo1)
      polej = knownj - hemi*rsw*SIN(alo1)

      ! find radius to desired point
      ala = lat*RAD_PER_DEG
      rm = rebydx*COS(ala)*scale_top/ (1.d0+hemi*SIN(ala))
      alo = (lon-reflon)*RAD_PER_DEG
      i = polei + rm*COS(alo)
      j = polej + hemi*rm*SIN(alo)

    ! lambert
    ELSE IF (map_proj.EQ.1) THEN

      IF (ABS(truelat2).GT.90.d0) THEN
          truelat2 = truelat1
      END IF

      IF (ABS(truelat1-truelat2).GT.0.1d0) THEN
          cone = (DLOG(COS(truelat1*RAD_PER_DEG))-DLOG(COS(truelat2*RAD_PER_DEG)))/&
                 (DLOG(TAN((90.d0-ABS(truelat1))*RAD_PER_DEG*0.5d0))-&
                 DLOG(TAN((90.d0-ABS(truelat2))*RAD_PER_DEG*0.5d0)))
      ELSE
          cone = SIN(ABS(truelat1)*RAD_PER_DEG)
      END IF

      ! compute longitude differences and ensure we stay
      ! out of the forbidden "cut zone"
      deltalon1 = lon1 - stdlon
      IF (deltalon1.GT.+180.d0) deltalon1 = deltalon1 - 360.d0
      IF (deltalon1.LT.-180.d0) deltalon1 = deltalon1 + 360.d0

      ! convert truelat1 to radian and compute cos for later use
      tl1r = truelat1*RAD_PER_DEG
      ctl1r = COS(tl1r)

      ! compute the radius to our known lower-left (sw) corner
      rsw = rebydx*ctl1r/cone*(TAN((90.d0*hemi-lat1)*RAD_PER_DEG/2.d0)/&
            TAN((90.d0*hemi-truelat1)*RAD_PER_DEG/2.d0))**cone

      ! find pole point
      arg = cone* (deltalon1*RAD_PER_DEG)
      polei = hemi*knowni - hemi*rsw*SIN(arg)
      polej = hemi*knownj + rsw*COS(arg)

      ! compute deltalon between known longitude and standard
      ! lon and ensure it is not in the cut zone
      deltalon = lon - stdlon
      IF (deltalon.GT.+180.d0) deltalon = deltalon - 360.d0
      IF (deltalon.LT.-180.d0) deltalon = deltalon + 360.d0

      ! radius to desired point
      rm = rebydx*ctl1r/cone*(TAN((90.d0*hemi-lat)*RAD_PER_DEG/2.d0)/&
                  TAN((90.d0*hemi-truelat1)*RAD_PER_DEG/2.d0))**cone

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
    ELSE IF (map_proj.EQ.6) THEN

      IF (pole_lat.NE.90.d0) THEN
          CALL rotatecoords(lat,lon,olat,olon,pole_lat,pole_lon,stdlon,-1)
          lat = olat
          lon = olon + stdlon
      END IF

      ! make sure center lat/lon is good
      IF (pole_lat.NE.90.d0) THEN
          CALL rotatecoords(lat1,lon1,olat,olon,pole_lat,pole_lon,stdlon,-1)
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

      CALL throw_exception('Do not know map projection ', map_proj)
      ! TODO throw exception

    END IF

    loc(1) = j
    loc(2) = i

    RETURN

END SUBROUTINE f_lltoij


SUBROUTINE f_ijtoll(map_proj,truelat1,truelat2,stdlon,lat1,lon1,&
                   pole_lat,pole_lon,knowni,knownj,dx,latinc,&
                   loninc,ai,aj,throw_exception,loc)

    ! converts input lat/lon values to the cartesian (i,j) value
    ! for the given projection.
    IMPLICIT NONE

    EXTERNAL throw_exception
    INTEGER,INTENT(IN) :: map_proj
    REAL(KIND=8),INTENT(IN) :: stdlon
    REAL(KIND=8),INTENT(IN) :: lat1,lon1,pole_lat,pole_lon,knowni,knownj
    REAL(KIND=8),INTENT(IN) :: dx,latinc,loninc,ai,aj
    REAL(KIND=8),INTENT(INOUT) :: truelat1,truelat2
    REAL(KIND=8),DIMENSION(2),INTENT(OUT) :: loc

    REAL(KIND=8) :: gi2
    REAL(KIND=8) :: arccos
    REAL(KIND=8) :: deltalon1
    REAL(KIND=8) :: tl1r
    REAL(KIND=8) :: clain,dlon,rsw,deltalon,deltalat
    REAL(KIND=8) :: reflon,scale_top,ala1,alo1,polei,polej
    ! earth radius divided by dx
    REAL(KIND=8) :: rebydx
    REAL(KIND=8) :: ctl1r,cone,hemi

    REAL(KIND=8),PARAMETER :: PI = 3.141592653589793d0
    REAL(KIND=8),PARAMETER :: RAD_PER_DEG = PI/180.d0
    REAL(KIND=8),PARAMETER :: DEG_PER_RAD = 180.d0/PI
    REAL(KIND=8),PARAMETER :: RE_M = 6370000.d0 ! radius of sperical earth

    REAL(KIND=8) :: inew,jnew,r,r2
    REAL(KIND=8) :: chi,chi1,chi2
    REAL(KIND=8) :: xx,yy,lat,lon

    REAL(KIND=8) :: olat,olon,lat1n,lon1n


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

    rebydx = RE_M/dx

    hemi = 1.0d0
    IF (truelat1.LT.0.0d0) THEN
      hemi = -1.0d0
    END IF

    ! mercator
    IF (map_proj.EQ.3) THEN

      ! preliminary variables
      clain = COS(RAD_PER_DEG*truelat1)
      dlon = dx/(RE_M*clain)

      ! compute distance from equator to origin, and store in
      ! the rsw tag.
      rsw = 0.d0
      IF (lat1.NE.0.d0) THEN
          rsw = (DLOG(TAN(0.5d0*((lat1+90.d0)*RAD_PER_DEG))))/dlon
      END IF

      lat = 2.0d0*ATAN(EXP(dlon*(rsw+aj-knownj)))*DEG_PER_RAD - 90.d0
      lon = (ai-knowni)*dlon*DEG_PER_RAD + lon1
      IF (lon.GT.180.d0) lon = lon - 360.d0
      IF (lon.LT.-180.d0) lon = lon + 360.d0

    ! ps
    ELSE IF (map_proj.EQ.2) THEN

      ! compute the reference longitude by rotating 90 degrees to
      ! the east to find the longitude line parallel to the
      ! positive x-axis.
      reflon = stdlon + 90.d0

      ! compute numerator term of map scale factor
      scale_top = 1.d0 + hemi*SIN(truelat1*RAD_PER_DEG)

      ! compute radius to known point
      ala1 = lat1*RAD_PER_DEG
      rsw = rebydx*COS(ala1)*scale_top/(1.d0+hemi*SIN(ala1))

      ! find the pole point
      alo1 = (lon1-reflon)*RAD_PER_DEG
      polei = knowni - rsw*COS(alo1)
      polej = knownj - hemi*rsw*SIN(alo1)

      ! compute radius to point of interest
      xx = ai - polei
      yy = (aj-polej)*hemi
      r2 = xx**2 + yy**2

      ! now the magic code
      IF (r2.EQ.0.d0) THEN
          lat = hemi*90.d0
          lon = reflon
      ELSE
          gi2 = (rebydx*scale_top)**2.d0
          lat = DEG_PER_RAD*hemi*ASIN((gi2-r2)/ (gi2+r2))
          arccos = ACOS(xx/SQRT(r2))
          IF (yy.GT.0) THEN
              lon = reflon + DEG_PER_RAD*arccos
          ELSE
              lon = reflon - DEG_PER_RAD*arccos
          END IF
      END IF

      ! convert to a -180 -> 180 east convention
      IF (lon.GT.180.d0) lon = lon - 360.d0
      IF (lon.LT.-180.d0) lon = lon + 360.d0

    !     !lambert
    ELSE IF (map_proj.EQ.1) THEN

      IF (ABS(truelat2).GT.90.d0) THEN
          truelat2 = truelat1
      END IF

      IF (ABS(truelat1-truelat2).GT.0.1d0) THEN
          cone = (DLOG(COS(truelat1*RAD_PER_DEG))-DLOG(COS(truelat2*RAD_PER_DEG)))/&
                 (DLOG(TAN((90.d0-ABS(truelat1))*RAD_PER_DEG*0.5d0))-&
                  DLOG(TAN((90.d0-ABS(truelat2))*RAD_PER_DEG*0.5d0)))
      ELSE
          cone = SIN(ABS(truelat1)*RAD_PER_DEG)
      END IF

      ! compute longitude differences and ensure we stay out of the
      ! forbidden "cut zone"
      deltalon1 = lon1 - stdlon
      IF (deltalon1.GT.+180.d0) deltalon1 = deltalon1 - 360.d0
      IF (deltalon1.LT.-180.d0) deltalon1 = deltalon1 + 360.d0

      ! convert truelat1 to radian and compute cos for later use
      tl1r = truelat1*RAD_PER_DEG
      ctl1r = COS(tl1r)

      ! compute the radius to our known point
      rsw = rebydx*ctl1r/cone*(TAN((90.d0*hemi-lat1)*RAD_PER_DEG/2.d0)/&
                         TAN((90.d0*hemi-truelat1)*RAD_PER_DEG/2.d0))**cone

      ! find pole point
      alo1 = cone* (deltalon1*RAD_PER_DEG)
      polei = hemi*knowni - hemi*rsw*SIN(alo1)
      polej = hemi*knownj + rsw*COS(alo1)

      chi1 = (90.d0-hemi*truelat1)*RAD_PER_DEG
      chi2 = (90.d0-hemi*truelat2)*RAD_PER_DEG

      ! see if we are in the southern hemispere and flip the
      ! indices if we are.
      inew = hemi*ai
      jnew = hemi*aj

      ! compute radius**2 to i/j location
      reflon = stdlon + 90.d0
      xx = inew - polei
      yy = polej - jnew
      r2 = (xx*xx+yy*yy)
      r = sqrt(r2)/rebydx

      ! convert to lat/lon
      IF (r2.EQ.0.d0) THEN
          lat = hemi*90.d0
          lon = stdlon
      ELSE
          lon = stdlon + DEG_PER_RAD*ATAN2(hemi*xx,yy)/cone
          lon = dmod(lon+360.d0,360.d0)
          IF (chi1.EQ.chi2) THEN
              chi = 2.0d0*ATAN((r/TAN(chi1))** (1.d0/cone)*TAN(chi1*0.5d0))
          ELSE
              chi = 2.0d0*ATAN((r*cone/SIN(chi1))** (1.d0/cone)*TAN(chi1*0.5d0))
          END IF
          lat = (90.0d0-chi*DEG_PER_RAD)*hemi
      END IF

      IF (lon.GT.+180.d0) lon = lon - 360.d0
      IF (lon.LT.-180.d0) lon = lon + 360.d0

    !     !lat-lon
    ELSE IF (map_proj.EQ.6) THEN

      inew = ai - knowni
      jnew = aj - knownj

      IF (inew.LT.0.d0) inew = inew + 360.d0/loninc
      IF (inew.GE.360.d0/dx) inew = inew - 360.d0/loninc

      ! compute deltalat and deltalon
      deltalat = jnew*latinc
      deltalon = inew*loninc

      IF (pole_lat.NE.90.d0) THEN
          CALL rotatecoords(lat1,lon1,olat,olon,pole_lat,pole_lon,stdlon,-1)
          lat1n = olat
          lon1n = olon + stdlon
          lat = deltalat + lat1n
          lon = deltalon + lon1n
      ELSE
          lat = deltalat + lat1
          lon = deltalon + lon1
      END IF


      IF (pole_lat.NE.90.d0) THEN
          lon = lon - stdlon
          CALL rotatecoords(lat,lon,olat,olon,pole_lat,pole_lon,stdlon,1)
          lat = olat
          lon = olon
      END IF

      IF (lon.LT.-180.d0) lon = lon + 360.d0
      IF (lon.GT.180.d0) lon = lon - 360.d0

    ELSE

      CALL throw_exception('Do not know map projection ',map_proj)

    END IF

    loc(1) = lat
    loc(2) = lon

    RETURN

END SUBROUTINE f_ijtoll

! Eta = (p - ptop) / (psfc - ptop)
! Potential temperature:
! theta = T * (Po/P)^(R/CP)
! Hypsometric equation:
! h = z2-z1 = R*Tbar/G * ln(p1/p2)
! where z1, p1 are the surface
SUBROUTINE f_converteta(full_t, znu, psfc, ptop, pcalc, mean_t, temp_t,&
                        z, nx,ny,nz)

    IMPLICIT NONE

    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(IN) :: full_t
    REAL(KIND=8),DIMENSION(nz),INTENT(IN) :: znu
    REAL(KIND=8),DIMENSION(nx,ny),INTENT(IN) :: psfc
    REAL(KIND=8),INTENT(IN) :: ptop
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: pcalc
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: mean_t
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(INOUT) :: temp_t
    REAL(KIND=8),DIMENSION(nx,ny,nz),INTENT(OUT) :: z
    INTEGER,INTENT(IN) :: nx,ny,nz

    REAL(KIND=8) :: s, cnt, avg
    REAL(KIND=8),PARAMETER :: R=287.06, G=9.81, CP=1005.0, P0=100000.0
    INTEGER i,j,k,kk

    DO k=1,nz
        DO j=1,ny
            DO i=1,nx
                pcalc(i,j,k) = (znu(k) * (psfc(i,j) - ptop)) + ptop
            END DO
        END DO
    END DO

    DO k=1,nz
        DO j=1,ny
            DO i=1,nx
                temp_t(i,j,k) = (full_t(i,j,k)) / ((P0 / (pcalc(i,j,k)))**(R/CP))
            END DO
        END DO
    END DO

    DO k=1,nz
        DO j = 1, ny
            DO i = 1, nx
                s = 0
                cnt = 0
                DO kk=1,k
                    s = s + temp_t(i,j,kk)
                    cnt = cnt + 1
                END DO
                avg = s / cnt
                mean_t(i,j,k) = avg
            END DO
        END DO
    END DO

    DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            z(i,j,k) = ((R*mean_t(i,j,k))/G) * LOG(psfc(i,j)/pcalc(i,j,k))
          END DO
        END DO
    END DO

    RETURN

END SUBROUTINE f_converteta


SUBROUTINE f_computectt(prs,tk,qci,qcw,qvp,ght,ter,haveqci,ctt,ew,ns,nz)

    IMPLICIT NONE

    REAL(KIND=8),DIMENSION(ew,ns,nz), INTENT(IN) :: ght, prs, tk, qci, qcw, qvp
    REAL(KIND=8),DIMENSION(ew,ns), INTENT(IN) :: ter
    REAL(KIND=8),DIMENSION(ew,ns), INTENT(OUT) :: ctt
    INTEGER, INTENT(IN) :: nz,ns,ew,haveqci
    !     REAL(KIND=8) ::     znfac(nz)

    ! LOCAL VARIABLES
    INTEGER i,j,k,ripk
    !INTEGER :: mjx,miy,mkzh
    REAL(KIND=8) :: vt,opdepthu,opdepthd,dp
    REAL(KIND=8) :: ratmix,arg1,arg2,agl_hgt
    REAL(KIND=8) :: fac,prsctt
    !REAL(KIND=8) :: eps,ussalr,rgas,grav,abscoefi,abscoef,celkel,wrfout
    !REAL(KIND=8) ::    ght(ew,ns,nz),stuff(ew,ns)
    !REAL(KIND=8), DIMENSION(ew,ns,nz) ::     pf(ns,ew,nz),p1,p2
    REAL(KIND=8), DIMENSION(ew,ns,nz) :: pf
    REAL(KIND=8) :: p1, p2

    REAL(KIND=8), PARAMETER :: EPS = 0.622d0
    REAL(KIND=8), PARAMETER :: USSALR = .0065d0      ! deg C per m
    REAL(KIND=8), PARAMETER :: RGAS = 287.04d0     !J/K/kg
    REAL(KIND=8), PARAMETER :: GRAV = 9.81d0
    REAL(KIND=8), PARAMETER :: ABSCOEFI = .272d0  ! cloud ice absorption coefficient in m^2/g
    REAL(KIND=8), PARAMETER :: ABSCOEF =.145d0   ! cloud water absorption coefficient in m^2/g
    REAL(KIND=8), PARAMETER :: CELKEL = 273.15d0
    REAL(KIND=8), PARAMETER :: WRFOUT = 1

    !mjx = ew
    !miy = ns
    !mkzh = nz

    prsctt = 0 ! removes the warning

! Calculate the surface pressure
    DO j=1,ns
        DO i=1,ew
           ratmix     = .001d0*qvp(i,j,1)
           arg1       = EPS + ratmix
           arg2       = EPS*(1.+ratmix)
           vt         = tk(i,j,1) * arg1/arg2 !Virtual temperature
           agl_hgt    = ght(i,j,nz) - ter(i,j)
           arg1       = -GRAV/(RGAS*USSALR)
           pf(i,j,nz) = prs(i,j,1)*(vt/(vt+USSALR*(agl_hgt)))**(arg1)
        END DO
    END DO

    DO k=1,nz-1
        DO j=1,ns
          DO i=1,ew
            ripk = nz-k+1
            pf(i,j,k)=.5d0*(prs(i,j,ripk)+prs(i,j,ripk-1))
          END DO
        END DO
    END DO

    DO j=1,ns
        DO i=1,ew
            opdepthd=0.d0
            k=0

!      Integrate downward from model top, calculating path at full
!      model vertical levels.

!20          opdepthu=opdepthd

            DO k=1, nz
                opdepthu=opdepthd
                !k=k+1
                ripk = nz-k+1

                IF (k.EQ.1) THEN
                    dp=200.d0*(pf(i,j,1)-prs(i,j,nz))  ! should be in Pa
                ELSE
                    dp=100.d0*(pf(i,j,k)-pf(i,j,k-1))  ! should be in Pa
                END IF

                IF (haveqci .EQ. 0) then
                    IF (tk(i,j,k) .LT. CELKEL) then
                        ! Note: abscoefi is m**2/g, qcw is g/kg, so no convrsion needed
                        opdepthd=opdepthu+ABSCOEFI*qcw(i,j,k)*dp/GRAV
                    ELSE
                        opdepthd=opdepthu+ABSCOEF*qcw(i,j,k)*dp/GRAV
                    END IF
                ELSE
                    opdepthd=opdepthd+(ABSCOEF*qcw(i,j,ripk)+ABSCOEFI*qci(i,j,ripk))*dp/GRAV
                END IF

                IF (opdepthd.LT.1. .AND. k.LT.nz) THEN
                    !GOTO 20
                    CYCLE

                ELSE IF (opdepthd.LT.1. .AND. k.EQ.nz) THEN
                    prsctt=prs(i,j,1)
                    EXIT
                ELSE
                    fac=(1.-opdepthu)/(opdepthd-opdepthu)
                    prsctt=pf(i,j,k-1)+fac*(pf(i,j,k)-pf(i,j,k-1))
                    prsctt=min(prs(i,j,1),max(prs(i,j,nz),prsctt))
                    EXIT
                END IF
            END DO

            DO k=2,nz
                ripk = nz-k+1
                p1 = prs(i,j,ripk+1)
                p2 = prs(i,j,ripk)
                IF (prsctt .GE. p1 .AND. prsctt .LE. p2) THEN
                    fac=(prsctt-p1)/(p2-p1)
                    arg1 = fac*(tk(i,j,ripk) - tk(i,j,ripk+1)) - CELKEL
                    ctt(i,j) = tk(i,j,ripk+1) + arg1
                    !GOTO 40
                    EXIT
                END IF
            END DO
        END DO
    END DO
!   30    CONTINUE
!   40    CONTINUE
! 190  CONTINUE
    RETURN

END SUBROUTINE f_computectt

SUBROUTINE f_filter2d(a, b, missing, it, nx, ny)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx, ny, it
    REAL(KIND=8), DIMENSION(nx, ny), INTENT(INOUT) :: a
    REAL(KIND=8), INTENT(IN) :: missing
    REAL(KIND=8), DIMENSION(nx, ny), INTENT(INOUT) :: b

    REAL(KIND=8), PARAMETER :: COEF=0.25D0

    INTEGER :: i, j, iter

    DO iter=1,it
        DO j=1,ny
            DO i = 1,nx
                b(i,j) = a(i,j)
            END DO
        END DO

        DO j=2,ny-1
            DO i=1,nx
                IF (b(i,j-1) .EQ. missing .OR. b(i,j) .EQ. missing .OR. &
                    b(i,j+1) .EQ. missing) THEN
                    a(i,j) = a(i,j)
                ELSE
                    a(i,j) = a(i,j) + COEF*(b(i,j-1)-2*b(i,j) + b(i,j+1))
                END IF
            END DO
        END DO

        DO j=1,ny
            DO i=2,nx-1
                IF (b(i-1,j) .EQ. missing .OR. b(i,j) .EQ. missing .OR. &
                    b(i+1,j) .EQ. missing) THEN
                    a(i,j) = a(i,j)
                ELSE
                    a(i,j) = a(i,j) + COEF*(b(i-1,j)-2*b(i,j)+b(i+1,j))
                END IF
            END DO
        END DO
    !        do j=1,ny
    !        do i=1,nx
    !          b(i,j) = a(i,j)
    !        enddo
    !        enddo
    !        do j=2,ny-1
    !        do i=1,nx
    !          a(i,j) = a(i,j) - .99*coef*(b(i,j-1)-2*b(i,j)+b(i,j+1))
    !        enddo
    !        enddo
    !        do j=1,ny
    !        do i=2,nx-1
    !          a(i,j) = a(i,j) - .99*coef*(b(i-1,j)-2*b(i,j)+b(i+1,j))
    !        enddo
    !        enddo
    END DO

    RETURN

END SUBROUTINE f_filter2d

SUBROUTINE f_monotonic(out,in,lvprs,cor,idir,delta,ew,ns,nz,icorsw)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: idir,ew,ns,nz,icorsw
    REAL(KIND=8), INTENT(IN) :: delta
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(INOUT) :: in
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(OUT) :: out
    REAL(KIND=8), DIMENSION(ew,ns,nz) :: lvprs
    REAL(KIND=8), DIMENSION(ew,ns) :: cor

    INTEGER :: i,j,k,ripk,k300

    k300 = 0 ! removes the warning

    DO j=1,ns
        DO i=1,ew
            IF (icorsw .EQ. 1 .AND. cor(i,j) .LT. 0.) THEN
                DO k=1,nz
                    in(i,j,k) = -in(i,j,k)
                END DO
            END IF

            ! First find k index that is at or below (height-wise)
            ! the 300 hPa level.
            DO k = 1,nz
                ripk = nz-k+1
                IF (lvprs(i,j,k) .LE. 300.d0) THEN
                    k300 = k
                    EXIT
                END IF
            END DO

            DO k = k300, 1,-1
                IF (idir .EQ. 1) THEN
                    out(i,j,k) = MIN(in(i,j,k), in(i,j,k+1)+delta)
                ELSE IF (idir .EQ. -1) THEN
                    out(i,j,k) = MAX(in(i,j,k), in(i,j,k+1)-delta)
                END IF
            END DO

            DO k = k300+1, nz
                IF (idir .EQ. 1) THEN
                    out(i,j,k) = MAX(in(i,j,k), in(i,j,k-1)-delta)
                ELSE IF (idir .EQ. -1) THEN
                    out(i,j,k) = MIN(in(i,j,k), in(i,j,k-1)+delta)
                END IF
            END DO
        END DO
    END DO

    RETURN

END SUBROUTINE f_monotonic


FUNCTION f_intrpvalue (wvalp0,wvalp1,vlev,vcp0,vcp1,icase)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: icase
    REAL(KIND=8), INTENT(IN) :: wvalp0,wvalp1,vlev,vcp0,vcp1
    REAL(KIND=8) :: f_intrpvalue

    REAL(KIND=8) :: valp0,valp1,rvalue
    REAL(KIND=8) :: chkdiff

    REAL(KIND=8), PARAMETER :: RGAS=287.04d0
    REAL(KIND=8), PARAMETER :: USSALR=0.0065d0
    REAL(KIND=8), PARAMETER :: SCLHT=RGAS*256.d0/9.81d0

    valp0 = wvalp0
    valp1 = wvalp1
    IF ( icase .EQ. 2) THEN  !GHT
       valp0=EXP(-wvalp0/SCLHT)
       valp1=EXP(-wvalp1/SCLHT)
    END IF

    ! TODO:  Remove this STOP
    chkdiff = vcp1 - vcp0
    IF(chkdiff .EQ. 0) THEN
        PRINT *,"bad difference in vcp's"
        STOP
    END IF

    rvalue = (vlev-vcp0)*(valp1-valp0)/(vcp1-vcp0)+valp0
    IF (icase .EQ. 2) THEN  !GHT
        f_intrpvalue = -SCLHT*LOG(rvalue)
    ELSE
        f_intrpvalue = rvalue
    END IF

    RETURN

END FUNCTION f_intrpvalue

! NOTES:
! vcarray is the array holding the values for the vertical coordinate.
! It will always come in with the dimensions of the staggered U and V grid.

SUBROUTINE  f_vintrp(datain,dataout,pres,tk,qvp,ght,terrain,&
                       sfp,smsfp,vcarray,interp_levels,numlevels,&
                       icase,ew,ns,nz,extrap,vcor,logp,rmsg)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ew,ns,nz,icase,extrap
    INTEGER, INTENT(IN) :: vcor,numlevels,logp
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: datain,pres,tk,qvp
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: ght
    REAL(KIND=8), DIMENSION(ew,ns), INTENT(IN) :: terrain,sfp,smsfp
    REAL(KIND=8), DIMENSION(ew,ns,numlevels), INTENT(OUT) :: dataout
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: vcarray
    REAL(KIND=8), DIMENSION(numlevels), INTENT(IN) :: interp_levels
    REAL(KIND=8), INTENT(IN) :: rmsg

    INTEGER :: nreqlvs,ripk !njx,niy
    INTEGER :: i,j,k,kupper !itriv
    INTEGER :: ifound,isign !miy,mjx
    REAL(KIND=8), DIMENSION(ew,ns) :: tempout
    REAL(KIND=8) :: rlevel,vlev,diff
    REAL(KIND=8) :: tmpvlev
    REAL(KIND=8) :: vcp1,vcp0,valp0,valp1
!    REAL(KIND=8) :: cvc
    REAL(KIND=8) :: vclhsl,vctophsl !qvlhsl,ttlhsl
    REAL(KIND=8) :: f_intrpvalue
    REAL(KIND=8) :: plhsl,zlhsl,ezlhsl,tlhsl,psurf,pratio,tlev
    REAL(KIND=8) :: ezsurf,psurfsm,zsurf,qvapor,vt
    REAL(KIND=8) :: ezlev,plev,zlev,ptarget,dpmin,dp
    REAL(KIND=8) :: pbot,zbot,tbotextrap,e
    REAL(KIND=8) :: tlcl,gammam
    CHARACTER :: cvcord*1

    REAL(KIND=8), PARAMETER :: RGAS    = 287.04d0     !J/K/kg
    REAL(KIND=8), PARAMETER :: RGASMD  = .608d0
    REAL(KIND=8), PARAMETER :: USSALR  = .0065d0      ! deg C per m
    REAL(KIND=8), PARAMETER :: SCLHT   = RGAS*256.d0/9.81d0
    REAL(KIND=8), PARAMETER :: EPS     = 0.622d0
    REAL(KIND=8), PARAMETER :: RCONST  = -9.81d0/(RGAS * USSALR)
    REAL(KIND=8), PARAMETER :: EXPON   =  RGAS*USSALR/9.81d0
    REAL(KIND=8), PARAMETER :: EXPONI  =  1./EXPON
    REAL(KIND=8), PARAMETER :: TLCLC1   = 2840.d0
    REAL(KIND=8), PARAMETER :: TLCLC2   = 3.5d0
    REAL(KIND=8), PARAMETER :: TLCLC3   = 4.805d0
    REAL(KIND=8), PARAMETER :: TLCLC4   = 55.d0
    REAL(KIND=8), PARAMETER :: THTECON1 = 3376.d0 ! K
    REAL(KIND=8), PARAMETER :: THTECON2 = 2.54d0
    REAL(KIND=8), PARAMETER :: THTECON3 = 0.81d0
    REAL(KIND=8), PARAMETER :: CP       = 1004.d0
    REAL(KIND=8), PARAMETER :: CPMD     = 0.887d0
    REAL(KIND=8), PARAMETER :: GAMMA    = RGAS/CP
    REAL(KIND=8), PARAMETER :: GAMMAMD  = RGASMD-CPMD
    REAL(KIND=8), PARAMETER :: CELKEL   = 273.16d0

    ! Removes the warnings for uninitialized variables
    cvcord = ''
    plev = 0
    zlev = 0
    vlev = 0

    IF (vcor .EQ. 1) THEN
        cvcord = 'p'
    ELSE IF ((vcor .EQ. 2) .OR. (vcor .EQ. 3)) THEN
        cvcord = 'z'
    ELSE IF ((vcor .EQ. 4) .OR. (vcor .EQ. 5)) THEN
        cvcord = 't'
    END IF

    !miy = ns
    !mjx = ew
    !njx = ew
    !niy = ns

    DO j = 1,ns
        DO i = 1,ew
            tempout(i,j) = rmsg
        END DO
    END DO

    DO nreqlvs = 1,numlevels
        IF (cvcord .EQ. 'z') THEN
            ! Convert rlevel to meters from km
            rlevel = interp_levels(nreqlvs) * 1000.d0
            vlev = EXP(-rlevel/SCLHT)
        ELSE IF (cvcord .EQ. 'p') THEN
            vlev = interp_levels(nreqlvs)
        ELSE IF (cvcord .EQ. 't') THEN
            vlev = interp_levels(nreqlvs)
        END IF

        DO j=1,ns
            DO i=1,ew
                ! Get the interpolated value that is within the model domain
                ifound = 0
                DO k = 1,nz-1
                    ripk = nz-k+1
                    vcp1 = vcarray(i,j,ripk-1)
                    vcp0 = vcarray(i,j,ripk)
                    valp0 = datain(i,j,ripk)
                    valp1 = datain(i,j,ripk-1)

                    IF ((vlev .GE. vcp0 .AND. vlev .LE. vcp1) .OR. &
                        (vlev .LE. vcp0 .AND. vlev .GE. vcp1)) THEN
                        ! print *,i,j,valp0,valp1
                        IF ((valp0 .EQ. rmsg) .OR. (valp1 .EQ. rmsg)) THEN
                            tempout(i,j) = rmsg
                            ifound = 1
                        ELSE
                            IF (logp .EQ. 1) THEN
                                vcp1 = LOG(vcp1)
                                vcp0 = LOG(vcp0)
                                IF (vlev .EQ. 0.0d0) THEN
                                    PRINT *,"Pressure value = 0"
                                    PRINT *,"Unable to take log of 0"
                                    STOP
                                END IF
                                tmpvlev = LOG(vlev)
                            ELSE
                                tmpvlev = vlev
                            END IF
                            tempout(i,j) = f_intrpvalue(valp0,valp1,tmpvlev,vcp0,vcp1,icase)

                            ! print *,"one ",i,j,tempout(i,j)
                            ifound = 1
                        END IF
                        !GOTO 115 ! EXIT
                        EXIT
                    END IF
                END DO !end for the k loop
 !115  CONTINUE

                IF (ifound .EQ. 1) THEN !Grid point is in the model domain
                    !GOTO 333 ! CYCLE
                    CYCLE
                END IF

                !If the user has requested no extrapolatin then just assign
                !all values above or below the model level to rmsg.
                IF (extrap .EQ. 0) THEN
                    tempout(i,j) = rmsg
                    !GOTO 333 ! CYCLE
                    CYCLE
                END IF

                ! The grid point is either above or below the model domain
                ! First we will check to see if the grid point is above the
                ! model domain.
                vclhsl = vcarray(i,j,1)  !lowest model level
                vctophsl = vcarray(i,j,nz) !highest model level
                diff = vctophsl - vclhsl
                isign = NINT(diff/ABS(diff))

                IF (isign*vlev .GE. isign*vctophsl) THEN
                    ! Assign the highest model level to the out array
                    tempout(i,j) = datain(i,j,nz)
                    ! print *,"at warn",i,j,tempout(i,j)
                    !GOTO 333 ! CYCLE
                    CYCLE
                END IF

                ! Only remaining possibility is that the specified level is below
                ! lowest model level.  If lowest model level value is missing,
                ! set interpolated value to missing.

                IF (datain(i,j,1) .EQ. rmsg) THEN
                    tempout(i,j) = rmsg
                    !GOTO 333 ! CYCLE
                    CYCLE
                END IF

                ! If the field comming in is not a pressure,temperature or height
                ! field we can set the output array to the value at the lowest
                ! model level.

                tempout(i,j) = datain(i,j,1)

                ! For the special cases of pressure on height levels or height on
                ! pressure levels, or temperature-related variables on pressure or
                ! height levels, perform a special extrapolation based on
                ! US Standard Atmosphere.  Here we calcualate the surface pressure
                ! with the altimeter equation.  This is how RIP calculates the
                ! surface pressure.
                IF (icase .GT. 0) THEN
                    plhsl = pres(i,j,1) * 0.01d0  !pressure at lowest model level
                    zlhsl = ght(i,j,1)            !grid point height a lowest model level
                    ezlhsl = EXP(-zlhsl/SCLHT)
                    tlhsl = tk(i,j,1)             !temperature in K at lowest model level
                    zsurf = terrain(i,j)
                    qvapor = MAX((qvp(i,j,1)*.001d0),1.e-15)
                    ! virtual temperature
                    ! vt     = tlhsl * (eps + qvapor)/(eps*(1.0 + qvapor))
                    ! psurf  = plhsl * (vt/(vt+USSALR * (zlhsl-zsurf)))**rconst
                    psurf = sfp(i,j)
                    psurfsm = smsfp(i,j)
                    ezsurf = EXP(-zsurf/SCLHT)

                    ! The if for checking above ground
                    IF ((cvcord .EQ. 'z' .AND. vlev .LT. ezsurf) .OR. &
                        (cvcord .EQ. 'p' .AND. vlev .LT. psurf)) THEN

                    ! We are below the lowest data level but above the ground.
                    ! Use linear interpolation (linear in prs and exp-height).

                        IF (cvcord .EQ. 'p') THEN
                            plev = vlev
                            ezlev = ((plev-plhsl)*&
                                    ezsurf+(psurf-plev)*ezlhsl)/(psurf-plhsl)
                            zlev = -sclht*LOG(ezlev)
                            IF (icase .EQ. 2) THEN
                                tempout(i,j) = zlev
                                !GOTO 333 ! CYCLE
                                CYCLE
                            END IF

                        ELSE IF (cvcord .EQ. 'z') THEN
                            ezlev = vlev
                            zlev = -sclht*LOG(ezlev)
                            plev = ((ezlev-ezlhsl)*&
                                   psurf+(ezsurf-ezlev)*plhsl)/(ezsurf-ezlhsl)
                            IF (icase .EQ. 1) THEN
                                tempout(i,j) = plev
                                !GOTO 333 ! CYCLE
                                CYCLE
                            END IF
                        END IF

                    ELSE   !else for checking above ground
                        ptarget = psurfsm - 150.d0
                        dpmin = 1.e4
                        DO k=1,nz
                            ripk = nz-k+1
                            dp = ABS((pres(i,j,ripk) * 0.01d0) - ptarget)
                            IF (dp .GT. dpmin) THEN
                                !GOTO 334 ! EXIT
                                EXIT
                            END IF
                            dpmin = MIN(dpmin,dp)
                        END DO
         !334
                        kupper = k-1

                        ripk = nz - kupper + 1
                        pbot = MAX(plhsl,psurf)
                        zbot = MIN(zlhsl,zsurf)
                        pratio = pbot/(pres(i,j,ripk) * 0.01d0)
                        tbotextrap = tk(i,j,ripk)*(pratio)**EXPON
                        ! virtual temperature
                        vt = tbotextrap * (EPS + qvapor)/(EPS*(1.0d0 + qvapor))
                        IF (cvcord .EQ. 'p') THEN
                            plev = vlev
                            zlev = zbot + vt/USSALR*(1. - (vlev/pbot)**EXPON)
                            IF (icase .EQ. 2) THEN
                                tempout(i,j) = zlev
                                !GOTO 333 ! CYCLE
                                CYCLE
                            END IF
                        ELSE IF (cvcord .EQ. 'z') THEN
                            zlev = -sclht*LOG(vlev)
                            plev = pbot*(1. + USSALR/vt*(zbot-zlev))**EXPONI
                            IF (icase .EQ. 1) THEN
                                tempout(i,j) = plev
                                !GOTO 333 ! CYCLE
                                CYCLE
                            END IF
                        END IF
                    END IF !end if for checking above ground
                END IF !for icase gt 0

                IF (icase .GT. 2) THEN !extrapolation for temperature
                    tlev = tlhsl + (zlhsl - zlev)*USSALR
                    qvapor = MAX(qvp(i,j,1),1.e-15)
                    gammam = GAMMA*(1. + GAMMAMD*qvapor)
                    IF (icase .EQ. 3) THEN
                        tempout(i,j) = tlev - CELKEL
                    ELSE IF (icase .EQ. 4) THEN
                        tempout(i,j) = tlev
                    ! Potential temperature - theta
                    ELSE IF (icase .EQ. 5) THEN
                        tempout(i,j) = tlev*(1000.d0/plev)**gammam
                    ! extraolation for equivalent potential temperature
                    ELSE IF (icase .EQ. 6) THEN
                        e = qvapor*plev/(EPS + qvapor)
                        tlcl = TLCLC1/(LOG(tlev**TLCLC2/e) - TLCLC3) + TLCLC4
                        tempout(i,j)=tlev*(1000.d0/plev)**(gammam)*&
                                     EXP((THTECON1/tlcl - THTECON2)*&
                                     qvapor*(1. + THTECON3*qvapor))
                    END IF
                END IF

 !333  CONTINUE

            END DO
        END DO

        ! print *,"----done----",interp_levels(nreqlvs)
        DO j = 1,ns
            DO i = 1,ew
                dataout(i,j,nreqlvs) = tempout(i,j)
            END DO
        END DO

    END DO !end for the nreqlvs

    RETURN

END SUBROUTINE f_vintrp




