!NCLFORTSTART
SUBROUTINE CALCDBZ(prs, tmk, qvp, qra, qsn, qgr, sn0, ivarint, iliqskin, dbz, nx, ny, nz)
    USE constants, ONLY : GAMMA_SEVEN, RHOWAT, RHO_R, RHO_S, RHO_G, ALPHA, &
                          CELKEL, PI, RD

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: dbz

    !   Arguments
    INTEGER, INTENT(IN) :: nx, ny, nz
    INTEGER, INTENT(IN) :: sn0, ivarint, iliqskin
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: dbz
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: prs
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: tmk
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(INOUT) :: qvp
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(INOUT) :: qra
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(INOUT) :: qsn
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(INOUT) :: qgr

!NCLEND

    !   Local Variables
    INTEGER :: i, j, k
    REAL(KIND=8) :: temp_c, virtual_t
    REAL(KIND=8) :: gonv, ronv, sonv
    REAL(KIND=8) :: factor_g, factor_r, factor_s
    REAL(KIND=8) :: factorb_g, factorb_s
    REAL(KIND=8) :: rhoair, z_e

    !   Constants used to calculate variable intercepts
    REAL(KIND=8), PARAMETER :: R1 = 1.D-15
    REAL(KIND=8), PARAMETER :: RON = 8.D6
    REAL(KIND=8), PARAMETER :: RON2 = 1.D10
    REAL(KIND=8), PARAMETER :: SON = 2.D7
    REAL(KIND=8), PARAMETER :: GON = 5.D7
    REAL(KIND=8), PARAMETER :: RON_MIN = 8.D6
    REAL(KIND=8), PARAMETER :: RON_QR0 = 0.00010D0
    REAL(KIND=8), PARAMETER :: RON_DELQR0 = 0.25D0*RON_QR0
    REAL(KIND=8), PARAMETER :: RON_CONST1R = (RON2-RON_MIN)*0.5D0
    REAL(KIND=8), PARAMETER :: RON_CONST2R = (RON2+RON_MIN)*0.5D0

    !   Constant intercepts
    REAL(KIND=8), PARAMETER :: RN0_R = 8.D6
    REAL(KIND=8), PARAMETER :: RN0_S = 2.D7
    REAL(KIND=8), PARAMETER :: RN0_G = 4.D6

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

    IF (sn0 .EQ. 0) THEN
        DO k = 1,nz
            DO j = 1,ny
                DO i = 1,nx
                    IF (tmk(i,j,k) .LT. CELKEL) THEN
                        qsn(i,j,k) = qra(i,j,k)
                        qra(i,j,k) = 0.D0
                    END IF
                END DO
            END DO
        END DO
    END IF

    factor_r = GAMMA_SEVEN*1.D18*(1.D0/(PI*RHO_R))**1.75D0
    factor_s = GAMMA_SEVEN*1.D18*(1.D0/(PI*RHO_S))**1.75D0*(RHO_S/RHOWAT)**2*ALPHA
    factor_g = GAMMA_SEVEN*1.D18*(1.D0/(PI*RHO_G))**1.75D0*(RHO_G/RHOWAT)**2*ALPHA

    DO k = 1,nz
        DO j = 1,ny
            DO i = 1,nx
                virtual_t = tmk(i,j,k)*(0.622D0 + qvp(i,j,k))/(0.622D0*(1.D0 + qvp(i,j,k)))
                rhoair = prs(i,j,k) / (RD*virtual_t)

                ! Adjust factor for brightband, where snow or graupel particle
                ! scatters like liquid water (alpha=1.0) because it is assumed to
                ! have a liquid skin.
                IF (iliqskin .EQ. 1 .AND. tmk(i,j,k) .GT. CELKEL) THEN
                    factorb_s = factor_s/ALPHA
                    factorb_g = factor_g/ALPHA
                ELSE
                    factorb_s = factor_s
                    factorb_g = factor_g
                END IF

                ! Calculate variable intercept parameters
                IF (ivarint .EQ. 1) THEN

                    temp_c = DMIN1(-0.001D0, tmk(i,j,k)-CELKEL)
                    sonv = DMIN1(2.0D8, 2.0D6*EXP(-0.12D0*temp_c))

                    gonv = gon
                    IF (qgr(i,j,k) .GT. R1) THEN
                        gonv = 2.38D0 * (PI*RHO_G/(rhoair*qgr(i,j,k)))**0.92D0
                        gonv = MAX(1.D4, MIN(gonv,GON))
                    END IF

                    ronv = RON2
                    IF (qra(i,j,k) .GT. R1) THEN
                        ronv = RON_CONST1R*TANH((RON_QR0-qra(i,j,k))/RON_DELQR0) + RON_CONST2R
                    END IF

                ELSE
                    ronv = RN0_R
                    sonv = RN0_S
                    gonv = RN0_G
                END IF

                ! Total equivalent reflectivity factor (z_e, in mm^6 m^-3) is
                ! the sum of z_e for each hydrometeor species:

                z_e = factor_r*(rhoair*qra(i,j,k))**1.75D0/ronv**.75D0 + &
                    factorb_s*(rhoair*qsn(i,j,k))**1.75D0/sonv**.75D0 + &
                    factorb_g* (rhoair*qgr(i,j,k))**1.75D0/gonv**.75D0

                ! Adjust small values of Z_e so that dBZ is no lower than -30
                z_e = MAX(z_e, .001D0)

                ! Convert to dBZ
                dbz(i,j,k) = 10.D0*LOG10(z_e)
            END DO
        END DO
    END DO

    RETURN

END SUBROUTINE CALCDBZ

