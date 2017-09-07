! NCLFORTSTART
SUBROUTINE DCLOUDFRAC(pres, rh, lowc, midc, highc, nz, ns, ew)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: lowc, midc, highc

    INTEGER  nz, ns, ew
    REAL(KIND=8), DIMENSION(ew, ns, nz), INTENT(IN) :: pres, rh
    REAL(KIND=8), DIMENSION(ew, ns), INTENT(OUT) :: lowc, midc, highc

! NCLEND

    INTEGER i, j, k
    INTEGER kchi, kcmi, kclo

    ! Remove compiler warnings
    kchi = 0
    kcmi = 0
    kclo = 0

    DO j = 1,ns
        DO i = 1,ew
            DO k = 1,nz-1
                IF ( pres(i,j,k) .GT. 97000. ) kclo=k
                IF ( pres(i,j,k) .GT. 80000. ) kcmi=k
                IF ( pres(i,j,k) .GT. 45000. ) kchi=k
            END DO

        DO k = 1,nz-1
            IF (k .GE. kclo .AND. k .LT. kcmi) THEN
                lowc(i,j) = MAX(rh(i,j,k), lowc(i,j))
            ELSE IF (k .GE. kcmi .AND. k .LT. kchi) THEN ! mid cloud
                midc(i,j) = MAX(rh(i,j,k), midc(i,j))
            ELSE if (k .GE. kchi) THEN                  ! high cloud
                highc(i,j) = MAX(rh(i,j,k), highc(i,j))
            END IF
        END DO

        lowc(i,j)  = 4.0*lowc(i,j)/100. - 3.0
        midc(i,j)  = 4.0*midc(i,j)/100. - 3.0
        highc(i,j) = 2.5*highc(i,j)/100. - 1.5

        lowc(i,j)  = MIN(lowc(i,j), 1.0)
        lowc(i,j)  = MAX(lowc(i,j), 0.0)
        midc(i,j)  = MIN(midc(i,j), 1.0)
        midc(i,j)  = MAX(midc(i,j), 0.0)
        highc(i,j) = MIN(highc(i,j), 1.0)
        highc(i,j) = MAX(highc(i,j), 0.0)

       END DO
    END DO

    RETURN

END SUBROUTINE DCLOUDFRAC
