SUBROUTINE testfunc(a, b, c, nx, ny, nz, errstat, errmsg)
    USE wrf_constants, ONLY : ALGERR
    IMPLICIT NONE

    !threadsafe
    !f2py intent(in,out) :: b

    INTEGER, INTENT(IN) :: nx, ny, nz

    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: a
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: b
    CHARACTER(LEN=*), INTENT(IN) :: c
    INTEGER, INTENT(INOUT) :: errstat
    REAL(KIND=8), PARAMETER :: blah=123.45
    CHARACTER(LEN=*), INTENT(INOUT) :: errmsg

    INTEGER :: i,j,k

    DO k=1,nz
        DO j=1,ny
            DO i=1,nx
                b(i,j,k) = a(i,j,k)
            END DO
        END DO
    END DO

    errstat = ALGERR
    WRITE(errmsg, *) c(1:20), blah

    RETURN

END SUBROUTINE testfunc
