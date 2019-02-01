! NCLFORTSTART
SUBROUTINE DCOMPUTEWSPD(wspd, u, v, n)

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: wspd

    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), DIMENSION(n), INTENT(OUT) :: wspd
    REAL(KIND=8), DIMENSION(n), INTENT(IN) :: u, v
! NCLEND

    INTEGER i

    !$OMP PARALLEL DO SCHEDULE(runtime)
    DO i = 1,n
        wspd(i) = SQRT(u(i)*u(i) + v(i)*v(i))
    END DO
    !$OMP END PARALLEL DO

END SUBROUTINE DCOMPUTEWSPD


! NCLFORTSTART
SUBROUTINE DCOMPUTEWDIR(wdir, u, v, n)
    USE wrf_constants, ONLY : DEG_PER_RAD

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: wdir

    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), DIMENSION(n), INTENT(OUT) :: wdir
    REAL(KIND=8), DIMENSION(n), INTENT(IN) :: u, v
! NCLEND

    INTEGER i

    !$OMP PARALLEL DO SCHEDULE(runtime)
    DO i = 1,n
        wdir(i) = MOD(270.0 - ATAN2(v(i), u(i)) * DEG_PER_RAD, 360.)
    END DO
    !$OMP END PARALLEL DO

END SUBROUTINE DCOMPUTEWDIR

