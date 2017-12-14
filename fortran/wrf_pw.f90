!NCLFORTSTART
SUBROUTINE DCOMPUTEPW(p, tv, qv, ht, pw, nx, ny, nz, nzh)
    USE wrf_constants, ONLY : RD

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: pw

    INTEGER, INTENT(IN) :: nx, ny, nz, nzh
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: p, tv, qv
    REAL(KIND=8), DIMENSION(nx,ny,nzh), INTENT(IN) :: ht
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(OUT) :: pw

!NCLEND

    INTEGER :: i, j, k
    !REAL(KIND=8),PARAMETER :: R=287.06

    pw = 0

    !$OMP PARALLEL

    DO k=1,nz
        !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
        DO j=1,ny
            DO i=1,nx
                pw(i,j) = pw(i,j) + ((p(i,j,k)/(RD*tv(i,j,k)))*qv(i,j,k)*(ht(i,j,k+1) - ht(i,j,k)))
            END DO
        END DO
        !$OMP END DO
    END DO

    !$OMP END PARALLEL

    RETURN

END SUBROUTINE DCOMPUTEPW
