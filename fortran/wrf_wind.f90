! NCLFORTSTART
SUBROUTINE DCOMPUTEWSPD(wspd, u, v, nx, ny)

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: wspd

    INTEGER, INTENT(IN) :: nx, ny
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(OUT) :: wspd
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: u, v
! NCLEND

    INTEGER i, j

    !$OMP PARALLEL DO COLLAPSE(2)
    DO j = 1,ny
        DO i = 1,nx
            wspd(i,j) = SQRT(u(i,j)*u(i,j) + v(i,j)*v(i,j))
        END DO
    END DO
    !$OMP END PARALLEL DO

END SUBROUTINE DCOMPUTEWSPD


! NCLFORTSTART
SUBROUTINE DCOMPUTEWDIR(wdir, u, v, nx, ny)
    USE wrf_constants, ONLY : DEG_PER_RAD

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: wdir

    INTEGER, INTENT(IN) :: nx, ny
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(OUT) :: wdir
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: u, v
! NCLEND

    INTEGER i, j

    !$OMP PARALLEL DO COLLAPSE(2)
    DO j = 1,ny
        DO i = 1,nx
            wdir(i,j) = MOD(270.0 - ATAN2(v(i,j), u(i,j)) * DEG_PER_RAD, 360.)
        END DO
    END DO
    !$OMP END PARALLEL DO

END SUBROUTINE DCOMPUTEWDIR

