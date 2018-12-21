! The subroutines in this file were taken directly from RIP code written
! by Dr. Mark Stoelinga.  They were modified by Sherrie
! Fredrick(NCAR/MMM) to work with NCL February 2015.

!NCLFORTSTART
SUBROUTINE wrf_monotonic(out, in, lvprs, cor, idir, delta, ew, ns, nz, icorsw)

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: out

    INTEGER, INTENT(IN) :: idir, ew, ns, nz, icorsw
    REAL(KIND=8), INTENT(IN) :: delta
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(INOUT) :: in
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(OUT) :: out
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: lvprs
    REAL(KIND=8), DIMENSION(ew,ns), INTENT(IN) :: cor


!NCLEND

    INTEGER :: i, j, k, k300

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, k, k300) SCHEDULE(runtime)
    DO j=1,ns
        DO i=1,ew
            k300 = -1

            IF (icorsw .EQ. 1 .AND. cor(i,j) .LT. 0.) THEN
                DO k=1,nz
                    in(i,j,k) = -in(i,j,k)
                END DO
            END IF

            ! First find k index that is at or below (height-wise)
            ! the 300 hPa level.
            DO k = 1,nz-1
                IF (lvprs(i,j,k) .LE. 300.D0) THEN
                    k300 = k
                    EXIT
                END IF
            END DO

            ! If the search fails for some reason, use the second to last
            ! k index
            IF (k300 .EQ. -1) THEN
                k300 = nz-1
            END IF

            DO k = k300,1,-1
                IF (idir .EQ. 1) THEN
                    out(i,j,k) = MIN(in(i,j,k), in(i,j,k+1) + delta)
                ELSE IF (idir .EQ. -1) THEN
                    out(i,j,k) = MAX(in(i,j,k), in(i,j,k+1) - delta)
                END IF
            END DO

            DO k = k300+1, nz
                IF (idir .EQ. 1) THEN
                    out(i,j,k) = MAX(in(i,j,k), in(i,j,k-1) - delta)
                ELSE IF (idir .EQ. -1) THEN
                    out(i,j,k) = MIN(in(i,j,k), in(i,j,k-1) + delta)
                END IF
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE wrf_monotonic


!NCLFORTSTART
FUNCTION wrf_intrp_value(wvalp0, wvalp1, vlev, vcp0, vcp1, icase, errstat)
    USE wrf_constants, ONLY : ALGERR, SCLHT

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: icase
    REAL(KIND=8), INTENT(IN) :: wvalp0, wvalp1, vlev, vcp0, vcp1
    INTEGER, INTENT(INOUT) :: errstat
    REAL(KIND=8) :: wrf_intrp_value

!NCLEND

    REAL(KIND=8) :: valp0, valp1, rvalue
    REAL(KIND=8) :: chkdiff

    errstat = 0

    valp0 = wvalp0
    valp1 = wvalp1
    IF ( icase .EQ. 2) THEN  !GHT
        valp0=EXP(-wvalp0/SCLHT)
        valp1=EXP(-wvalp1/SCLHT)
    END IF

    chkdiff = vcp1 - vcp0
    IF(chkdiff .EQ. 0) THEN
        errstat = ALGERR
        !errmsg = "bad difference in vcp's"
        wrf_intrp_value = 0
        RETURN
        !PRINT *,"bad difference in vcp's"
        !STOP
    END IF

    rvalue = (vlev - vcp0)*(valp1 - valp0)/(vcp1 - vcp0) + valp0
    IF (icase .EQ. 2) THEN  !GHT
        wrf_intrp_value = -SCLHT*LOG(rvalue)
    ELSE
        wrf_intrp_value = rvalue
    END IF

    RETURN

END FUNCTION wrf_intrp_value

! NOTES:
! vcarray is the array holding the values for the vertical coordinate.
! It will always come in with the dimensions of the staggered U and V grid.

!NCLFORTSTART
SUBROUTINE wrf_vintrp(datain, dataout, pres, tk, qvp, ght, terrain,&
                      sfp, smsfp, vcarray, interp_levels, numlevels,&
                      icase, ew, ns, nz, extrap, vcor, logp, tempout, rmsg,&
                      errstat, errmsg)
    USE wrf_constants, ONLY : ALGERR, SCLHT, EXPON, EXPONI, GAMMA, GAMMAMD, TLCLC1, &
                          TLCLC2, TLCLC3, TLCLC4, THTECON1, THTECON2, THTECON3, &
                          CELKEL, EPS, USSALR

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: dataout

    INTEGER, INTENT(IN) :: ew, ns, nz, icase, extrap
    INTEGER, INTENT(IN) :: vcor, numlevels, logp
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: datain, pres, tk, qvp
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: ght
    REAL(KIND=8), DIMENSION(ew,ns), INTENT(IN) :: terrain, sfp, smsfp
    REAL(KIND=8), DIMENSION(ew,ns,numlevels), INTENT(OUT) :: dataout
    REAL(KIND=8), DIMENSION(ew,ns,nz), INTENT(IN) :: vcarray
    REAL(KIND=8), DIMENSION(numlevels), INTENT(IN) :: interp_levels

    REAL(KIND=8), DIMENSION(ew,ns), INTENT(INOUT) :: tempout

    REAL(KIND=8), INTENT(IN) :: rmsg
    INTEGER, INTENT(INOUT) :: errstat
    CHARACTER(LEN=*), INTENT(INOUT) :: errmsg

!NCLEND

    INTEGER :: nreqlvs, ripk !njx,niy
    INTEGER :: i, j, k, kupper !itriv
    INTEGER :: ifound, isign !miy,mjx
    INTEGER :: log_errcnt, interp_errcnt, interp_errstat
    REAL(KIND=8) :: rlevel, vlev, diff
    REAL(KIND=8) :: tmpvlev
    REAL(KIND=8) :: vcp1, vcp0, valp0, valp1
!    REAL(KIND=8) :: cvc
    REAL(KIND=8) :: vclhsl, vctophsl !qvlhsl,ttlhsl
    REAL(KIND=8) :: wrf_intrp_value
    REAL(KIND=8) :: plhsl, zlhsl, ezlhsl, tlhsl, psurf, pratio, tlev
    REAL(KIND=8) :: ezsurf, psurfsm, zsurf, qvapor, vt
    REAL(KIND=8) :: ezlev, plev, zlev, ptarget, dpmin, dp
    REAL(KIND=8) :: pbot, zbot, tbotextrap, e
    REAL(KIND=8) :: tlcl, gammam
    CHARACTER(LEN=1) :: cvcord

    ! Removes the warnings for uninitialized variables
    cvcord = ''
    plev = 0
    zlev = 0
    vlev = 0
    errstat = 0
    interp_errcnt = 0
    interp_errstat = 0
    log_errcnt = 0

    IF (vcor .EQ. 1) THEN
        cvcord = 'p'
    ELSE IF ((vcor .EQ. 2) .OR. (vcor .EQ. 3)) THEN
        cvcord = 'z'
    ELSE IF ((vcor .EQ. 4) .OR. (vcor .EQ. 5)) THEN
        cvcord = 't'
    END IF

    DO nreqlvs = 1,numlevels
        IF (cvcord .EQ. 'z') THEN
            ! Convert rlevel to meters from km
            rlevel = interp_levels(nreqlvs) * 1000.D0
            vlev = EXP(-rlevel/SCLHT)
        ELSE IF (cvcord .EQ. 'p') THEN
            vlev = interp_levels(nreqlvs)
        ELSE IF (cvcord .EQ. 't') THEN
            vlev = interp_levels(nreqlvs)
        END IF

        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i, j, k, ifound, &
        !$OMP ripk, vcp1, vcp0, valp0, valp1, tmpvlev, interp_errstat, &
        !$OMP vclhsl, vctophsl, diff, isign, plhsl, zlhsl, ezlhsl, tlhsl, &
        !$OMP zsurf, qvapor, psurf, psurfsm, ezsurf, plev, ezlev, zlev, &
        !$OMP ptarget, dpmin, kupper, pbot, zbot, pratio, tbotextrap, &
        !$OMP vt, tlev, gammam, e, tlcl) REDUCTION (+:log_errcnt, interp_errcnt) &
        !$OMP SCHEDULE(runtime)
        DO j=1,ns
            DO i=1,ew
                tempout(i,j) = rmsg
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
                                IF (vlev .NE. 0.0D0) THEN
                                    tmpvlev = LOG(vlev)
                                ELSE
                                    log_errcnt = log_errcnt + 1
                                    tmpvlev = rmsg
                                END IF
                            ELSE
                                tmpvlev = vlev
                            END IF

                            IF (tmpvlev .NE. rmsg) THEN
                                tempout(i,j) = wrf_intrp_value(valp0, valp1, tmpvlev, vcp0, &
                                                               vcp1, icase, interp_errstat)

                                IF (interp_errstat .NE. 0) THEN
                                    tempout(i,j) = rmsg
                                    interp_errcnt = interp_errcnt + 1
                                END IF

                                ifound = 1
                            END IF
                        END IF

                        EXIT
                    END IF
                END DO !end for the k loop

                IF (ifound .EQ. 1) THEN !Grid point is in the model domain
                    CYCLE
                END IF

                !If the user has requested no extrapolatin then just assign
                !all values above or below the model level to rmsg.
                IF (extrap .EQ. 0) THEN
                    tempout(i,j) = rmsg
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
                    CYCLE
                END IF

                ! Only remaining possibility is that the specified level is below
                ! lowest model level.  If lowest model level value is missing,
                ! set interpolated value to missing.

                IF (datain(i,j,1) .EQ. rmsg) THEN
                    tempout(i,j) = rmsg
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
                    plhsl = pres(i,j,1) * 0.01D0  !pressure at lowest model level
                    zlhsl = ght(i,j,1)            !grid point height a lowest model level
                    ezlhsl = EXP(-zlhsl/SCLHT)
                    tlhsl = tk(i,j,1)             !temperature in K at lowest model level
                    zsurf = terrain(i,j)
                    qvapor = MAX((qvp(i,j,1)*.001D0),1.e-15)
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
                            ezlev = ((plev - plhsl)*&
                                    ezsurf + (psurf - plev)*ezlhsl)/(psurf - plhsl)
                            zlev = -SCLHT*LOG(ezlev)
                            IF (icase .EQ. 2) THEN
                                tempout(i,j) = zlev
                                CYCLE
                            END IF

                        ELSE IF (cvcord .EQ. 'z') THEN
                            ezlev = vlev
                            zlev = -SCLHT*LOG(ezlev)
                            plev = ((ezlev - ezlhsl)*&
                                   psurf + (ezsurf - ezlev)*plhsl)/(ezsurf - ezlhsl)
                            IF (icase .EQ. 1) THEN
                                tempout(i,j) = plev * 100.
                                CYCLE
                            END IF
                        END IF

                    ELSE   !else for checking above ground
                        ptarget = psurfsm - 150.D0
                        dpmin = 1.E4
                        DO k=1,nz
                            ripk = nz-k+1
                            dp = ABS((pres(i,j,ripk) * 0.01D0) - ptarget)
                            IF (dp .GT. dpmin) THEN
                                EXIT
                            END IF
                            dpmin = MIN(dpmin, dp)
                        END DO

                        kupper = k-1

                        ripk = nz - kupper + 1
                        pbot = MAX(plhsl,psurf)
                        zbot = MIN(zlhsl,zsurf)
                        pratio = pbot/(pres(i,j,ripk) * 0.01D0)
                        tbotextrap = tk(i,j,ripk)*(pratio)**EXPON
                        ! virtual temperature
                        vt = tbotextrap * (EPS + qvapor)/(EPS*(1.0D0 + qvapor))
                        IF (cvcord .EQ. 'p') THEN
                            plev = vlev
                            zlev = zbot + vt/USSALR*(1. - (vlev/pbot)**EXPON)
                            IF (icase .EQ. 2) THEN
                                tempout(i,j) = zlev
                                CYCLE
                            END IF
                        ELSE IF (cvcord .EQ. 'z') THEN
                            zlev = -sclht*LOG(vlev)
                            plev = pbot*(1. + USSALR/vt*(zbot - zlev))**EXPONI
                            IF (icase .EQ. 1) THEN
                                tempout(i,j) = plev * 100.
                                CYCLE
                            END IF
                        END IF
                    END IF !end if for checking above ground
                END IF !for icase gt 0

                IF (icase .GT. 2) THEN !extrapolation for temperature
                    tlev = tlhsl + (zlhsl - zlev)*USSALR
                    qvapor = MAX(qvp(i,j,1), 1.e-15)
                    gammam = GAMMA*(1. + GAMMAMD*qvapor)
                    IF (icase .EQ. 3) THEN
                        tempout(i,j) = tlev - CELKEL
                    ELSE IF (icase .EQ. 4) THEN
                        tempout(i,j) = tlev
                    ! Potential temperature - theta
                    ELSE IF (icase .EQ. 5) THEN
                        tempout(i,j) = tlev*(1000.D0/plev)**gammam
                    ! extrapolation for equivalent potential temperature
                    ELSE IF (icase .EQ. 6) THEN
                        e = qvapor*plev/(EPS + qvapor)
                        tlcl = TLCLC1/(LOG(tlev**TLCLC2/e) - TLCLC3) + TLCLC4
                        tempout(i,j)=tlev*(1000.D0/plev)**(gammam)*&
                                     EXP((THTECON1/tlcl - THTECON2)*&
                                     qvapor*(1. + THTECON3*qvapor))
                    END IF
                END IF

 !333  CONTINUE

            END DO
        END DO
        !$OMP END PARALLEL DO

        IF (log_errcnt > 0) THEN
            errstat = ALGERR
            WRITE(errmsg, *) "Pres=0.  Unable to take log of 0."
            RETURN
        END IF

        IF (interp_errcnt > 0) THEN
            errstat = ALGERR
            WRITE(errmsg, *) "bad difference in vcp's"
            RETURN
        END IF

        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(runtime)
        DO j = 1,ns
            DO i = 1,ew
                dataout(i,j,nreqlvs) = tempout(i,j)
            END DO
        END DO
        !$OMP END PARALLEL DO

    END DO !end for the nreqlvs

    RETURN

END SUBROUTINE wrf_vintrp
