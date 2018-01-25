!NCLFORTSTART
SUBROUTINE DCOMPUTEABSVORT(av, u, v, msfu, msfv, msft, cor, dx, dy, nx, ny, nz,&
                          nxp1, nyp1)

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: av

    INTEGER, INTENT(IN) :: nx, ny, nz, nxp1, nyp1
    REAL(KIND=8), DIMENSION(nxp1,ny,nz), INTENT(IN) :: u
    REAL(KIND=8), DIMENSION(nx,nyp1,nz), INTENT(IN) :: v
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: av
    REAL(KIND=8), DIMENSION(nxp1,ny), INTENT(IN):: msfu
    REAL(KIND=8), DIMENSION(nx,nyp1), INTENT(IN) :: msfv
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: msft
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: cor
    REAL(KIND=8) :: dx, dy

!NCLEND

    INTEGER :: jp1, jm1, ip1, im1, i, j, k
    REAL(KIND=8) :: dsy, dsx, dudy, dvdx, avort
    REAL(KIND=8) :: mm

    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i, j, k, jp1, jm1, ip1, im1, &
    !$OMP dsx, dsy, mm, dudy, dvdx, avort) SCHEDULE(runtime)
    DO k = 1,nz
        DO j = 1,ny
            DO i = 1,nx
                jp1 = MIN(j+1, ny)
                jm1 = MAX(j-1, 1)
                ip1 = MIN(i+1, nx)
                im1 = MAX(i-1, 1)
                dsx = (ip1 - im1) * dx
                dsy = (jp1 - jm1) * dy
                mm = msft(i,j)*msft(i,j)
                dudy = 0.5D0*(u(i,jp1,k)/msfu(i,jp1) + u(i+1,jp1,k)/msfu(i+1,jp1) - &
                     u(i,jm1,k)/msfu(i,jm1) - u(i+1,jm1,k)/msfu(i+1,jm1))/dsy*mm
                dvdx = 0.5D0*(v(ip1,j,k)/msfv(ip1,j) + v(ip1,j+1,k)/msfv(ip1,j+1) - &
                     v(im1,j,k)/msfv(im1,j) - v(im1,j+1,k)/msfv(im1,j+1))/dsx*mm
                avort = dvdx - dudy + cor(i,j)
                av(i,j,k) = avort*1.D5

            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DCOMPUTEABSVORT


!NCLFORTSTART
SUBROUTINE DCOMPUTEPV(pv, u, v, theta, prs, msfu, msfv, msft, cor, dx, dy, nx, &
                      ny, nz, nxp1, nyp1)
    USE wrf_constants, ONLY : G

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: pv

    INTEGER,INTENT(IN) :: nx, ny, nz, nxp1, nyp1
    REAL(KIND=8), DIMENSION(nxp1,ny,nz), INTENT(IN) :: u
    REAL(KIND=8), DIMENSION(nx,nyp1,nz), INTENT(IN) :: v
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: prs
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: theta
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: pv
    REAL(KIND=8), DIMENSION(nxp1,ny), INTENT(IN) ::  msfu
    REAL(KIND=8), DIMENSION(nx,nyp1), INTENT(IN) :: msfv
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: msft
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: cor
    REAL(KIND=8) :: dx,dy

!NCLEND

    INTEGER :: kp1, km1, jp1, jm1, ip1, im1, i, j, k
    REAL(KIND=8) :: dsy, dsx, dp, dudy, dvdx, dudp, dvdp, dthdp, avort
    REAL(KIND=8) :: dthdx, dthdy, mm

    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i, j, k, kp1, km1, jp1, jm1, ip1, &
    !$OMP im1, dsx, dsy, mm, dudy, dvdx, avort, &
    !$OMP dp, dudp, dvdp, dthdp, dthdx, dthdy) SCHEDULE(runtime)
    DO k = 1,nz
        DO J = 1,ny
            DO i = 1,nx
                kp1 = MIN(k+1, nz)
                km1 = MAX(k-1, 1)
                jp1 = MIN(j+1, ny)
                jm1 = MAX(j-1, 1)
                ip1 = MIN(i+1, nx)
                im1 = MAX(i-1, 1)

                dsx = (ip1 - im1)*dx
                dsy = (jp1 - jm1)*dy
                mm = msft(i,j)*msft(i,j)

                dudy = 0.5D0*(u(i,jp1,k)/msfu(i,jp1) + u(i+1,jp1,k)/msfu(i+1,jp1) - &
                       u(i,jm1,k)/msfu(i,jm1) - u(i+1,jm1,k)/msfu(i+1,jm1))/dsy*mm
                dvdx = 0.5D0*(v(ip1,j,k)/msfv(ip1,j) + v(ip1,j+1,k)/msfv(ip1,j+1) - &
                       v(im1,j,k)/msfv(im1,j) - v(im1,j+1,k)/msfv(im1,j+1))/dsx*mm
                avort = dvdx - dudy + cor(i,j)
                dp = prs(i,j,kp1) - prs(i,j,km1)
                dudp = 0.5D0*(u(i,j,kp1) + u(i+1,j,kp1) - u(i,j,km1) - u(i+1,j,km1))/dp
                dvdp = 0.5D0*(v(i,j,kp1) + v(i,j+1,kp1) - v(i,j,km1) - v(i,J+1,km1))/dp
                dthdp = (theta(i,j,kp1) - theta(i,j,km1))/dp
                dthdx = (theta(ip1,j,k) - theta(im1,j,k))/dsx*msft(i,j)
                dthdy = (theta(i,jp1,k) - theta(i,jm1,k))/dsy*msft(i,j)
                pv(i,j,k) = -G*(dthdp*avort - dvdp*dthdx + dudp*dthdy)*10000.D0
                pv(i,j,k) = pv(i,j,k)*1.D2
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DCOMPUTEPV
