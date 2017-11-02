!NCLFORTSTART
SUBROUTINE wrfcttcalc(prs, tk, qci, qcw, qvp, ght, ter, ctt, haveqci, nz, ns, ew)
    USE wrf_constants, ONLY : EPS, USSALR, RD, G, ABSCOEFI, ABSCOEF, CELKEL

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: ctt

    INTEGER, INTENT(IN) :: nz, ns, ew, haveqci
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: ght, prs, tk, qci, qcw, qvp
    REAL(KIND=8), DIMENSION(ew,ns), INTENT(IN) :: ter
    REAL(KIND=8), DIMENSION(ew,ns), INTENT(OUT) :: ctt

!NCLEND

    !     REAL(KIND=8) ::     znfac(nz)

    ! LOCAL VARIABLES
    INTEGER i,j,k,ripk
    REAL(KIND=8) :: opdepthu, opdepthd, dp, arg1, fac, prsctt, ratmix
    REAL(KIND=8) :: arg2, agl_hgt, vt
    REAL(KIND=8), DIMENSION(ew,ns,nz) :: pf
    REAL(KIND=8) :: p1, p2

    !$OMP PARALLEL

    ! Calculate the surface pressure
    !$OMP DO COLLAPSE(2)
    DO j=1,ns
        DO i=1,ew
           ratmix = .001D0*qvp(i,j,1)
           arg1 = EPS + ratmix
           arg2 = EPS*(1. + ratmix)
           vt = tk(i,j,1)*arg1/arg2 !Virtual temperature
           agl_hgt = ght(i,j,nz) - ter(i,j)
           arg1 = -G/(RD*USSALR)
           pf(i,j,nz) = prs(i,j,1)*(vt/(vt + USSALR*(agl_hgt)))**(arg1)
        END DO
    END DO
    !$OMP END DO

    !$OMP DO COLLAPSE(3)
    DO k=1,nz-1
        DO j=1,ns
            DO i=1,ew
                ripk = nz-k+1
                pf(i,j,k) = .5D0*(prs(i,j,ripk) + prs(i,j,ripk-1))
            END DO
        END DO
    END DO
    !$OMP END DO

    !$OMP DO COLLAPSE(2) PRIVATE(i, j, k, ripk, opdepthd, opdepthu, &
    !$OMP prsctt, dp, p1, p2, fac, arg1)
    DO j=1,ns
        DO i=1,ew
            opdepthd = 0.D0
            k = 0
            prsctt = 0

            ! Integrate downward from model top, calculating path at full
            ! model vertical levels.

            DO k=1, nz
                opdepthu = opdepthd
                ripk = nz - k + 1

                IF (k .NE. 1) THEN
                    dp = 100.D0*(pf(i,j,k) - pf(i,j,k-1))  ! should be in Pa
                ELSE
                    dp = 200.D0*(pf(i,j,1) - prs(i,j,nz))  ! should be in Pa
                END IF

                IF (haveqci .EQ. 0) then
                    IF (tk(i,j,k) .LT. CELKEL) then
                        ! Note: abscoefi is m**2/g, qcw is g/kg, so no convrsion needed
                        opdepthd = opdepthu + ABSCOEFI*qcw(i,j,k) * dp/G
                    ELSE
                        opdepthd = opdepthu + ABSCOEF*qcw(i,j,k) * dp/G
                    END IF
                ELSE
                    opdepthd = opdepthd + (ABSCOEF*qcw(i,j,ripk) + ABSCOEFI*qci(i,j,ripk))*dp/G
                END IF

                IF (opdepthd .LT. 1. .AND. k .LT. nz) THEN
                    CYCLE

                ELSE IF (opdepthd .LT. 1. .AND. k .EQ. nz) THEN
                    prsctt = prs(i,j,1)
                    EXIT
                ELSE
                    fac = (1. - opdepthu)/(opdepthd - opdepthu)
                    prsctt = pf(i,j,k-1) + fac*(pf(i,j,k) - pf(i,j,k-1))
                    prsctt = MIN(prs(i,j,1), MAX(prs(i,j,nz), prsctt))
                    EXIT
                END IF
            END DO

            DO k=2,nz
                ripk = nz - k + 1
                p1 = prs(i,j,ripk+1)
                p2 = prs(i,j,ripk)
                IF (prsctt .GE. p1 .AND. prsctt .LE. p2) THEN
                    fac = (prsctt - p1)/(p2 - p1)
                    arg1 = fac*(tk(i,j,ripk) - tk(i,j,ripk+1)) - CELKEL
                    ctt(i,j) = tk(i,j,ripk+1) + arg1
                    EXIT
                END IF
            END DO
        END DO
    END DO
    !$OMP END DO

    !$OMP END PARALLEL
    RETURN

END SUBROUTINE wrfcttcalc
