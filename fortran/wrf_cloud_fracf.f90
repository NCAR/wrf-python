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
    lowc = 0
    midc = 0
    highc = 0

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, k, kchi, kcmi, kclo) &
    !$OMP SCHEDULE(runtime)
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
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DCLOUDFRAC


! NCLFORTSTART
SUBROUTINE DCLOUDFRAC2(vert, rh, vert_inc_w_height, low_thresh, mid_thresh, &
                       high_thresh, msg, lowc, midc, highc, nz, ns, ew)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: lowc, midc, highc

    INTEGER  nz, ns, ew
    REAL(KIND=8), DIMENSION(ew, ns, nz), INTENT(IN) :: rh, vert
    REAL(KIND=8), INTENT(IN) :: low_thresh, mid_thresh, high_thresh, msg
    INTEGER, INTENT(IN) :: vert_inc_w_height
    REAL(KIND=8), DIMENSION(ew, ns), INTENT(OUT) :: lowc, midc, highc

! NCLEND

    INTEGER i, j, k
    INTEGER kchi, kcmi, kclo

    ! Initialize the output
    lowc = 0
    midc = 0
    highc = 0

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, k, kchi, kcmi, kclo) &
    !$OMP SCHEDULE(runtime)
    DO j = 1,ns
        DO i = 1,ew
            ! A value of -1 means 'not found'.  This is needed to handle
            ! the mountains, where the level thresholds are below the lowest
            ! model level.
            kchi = -1
            kcmi = -1
            kclo = -1

            IF (vert_inc_w_height .NE. 0) THEN ! Vert coord increase with height
                DO k = 1,nz
                    IF (vert(i,j,k) .LT. low_thresh) kclo=k
                    IF (vert(i,j,k) .LT. mid_thresh) kcmi=k
                    IF (vert(i,j,k) .LT. high_thresh) kchi=k
                END DO
            ELSE ! Vert coord decrease with height
                DO k = 1,nz
                    IF (vert(i,j,k) .GT. low_thresh) kclo=k
                    IF (vert(i,j,k) .GT. mid_thresh) kcmi=k
                    IF (vert(i,j,k) .GT. high_thresh) kchi=k
                END DO
            ENDIF

            DO k = 1,nz
                IF (k .GE. kclo .AND. k .LT. kcmi) THEN
                    lowc(i,j) = MAX(rh(i,j,k), lowc(i,j))
                ELSE IF (k .GE. kcmi .AND. k .LT. kchi) THEN ! mid cloud
                    midc(i,j) = MAX(rh(i,j,k), midc(i,j))
                ELSE if (k .GE. kchi) THEN                  ! high cloud
                    highc(i,j) = MAX(rh(i,j,k), highc(i,j))
                END IF
            END DO

            ! Only do this when a cloud threshold is in the model vertical
            ! domain, otherwise it will be set to missing
            IF (kclo .GE. 1) THEN
                lowc(i,j)  = 4.0*lowc(i,j)/100. - 3.0
                lowc(i,j)  = MIN(lowc(i,j), 1.0)
                lowc(i,j)  = MAX(lowc(i,j), 0.0)
            ELSE
                lowc(i,j) = msg
            END IF

            IF (kcmi .GE. 1) THEN
                midc(i,j)  = 4.0*midc(i,j)/100. - 3.0
                midc(i,j)  = MIN(midc(i,j), 1.0)
                midc(i,j)  = MAX(midc(i,j), 0.0)
            ELSE
                midc(i,j) = msg
            END IF

            IF (kchi .GE. 1) THEN
                highc(i,j) = 2.5*highc(i,j)/100. - 1.5
                highc(i,j) = MIN(highc(i,j), 1.0)
                highc(i,j) = MAX(highc(i,j), 0.0)
            ELSE
                highc(i,j) = msg
            END IF

       END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DCLOUDFRAC2
