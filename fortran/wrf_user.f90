! NCLFORTSTART
SUBROUTINE DCOMPUTEPI(pi, pressure, nx, ny, nz)
    USE wrf_constants, ONLY : P1000MB, RD, CP

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: pi

    INTEGER, INTENT(IN) :: nx, ny, nz
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: pi
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: pressure
! NCLEND

    INTEGER i, j, k

    !REAL(KIND=8), PARAMETER :: P1000MB=100000.D0, R_D=287.D0, CP=7.D0*R_D/2.D0

    !$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(runtime)
    DO k = 1,nz
        DO j = 1,ny
            DO i = 1,nx
                pi(i,j,k) = (pressure(i,j,k)/P1000MB)**(RD/CP)
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO

END SUBROUTINE DCOMPUTEPI

! Temperature from potential temperature in kelvin.
!NCLFORTSTART
SUBROUTINE DCOMPUTETK(tk, pressure, theta, nx)
    USE wrf_constants, ONLY : P1000MB, RD, CP

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: tk

    INTEGER, INTENT(IN) :: nx
    !REAL(KIND=8) :: pi
    REAL(KIND=8), DIMENSION(nx), INTENT(IN) :: pressure
    REAL(KIND=8), DIMENSION(nx), INTENT(IN) :: theta
    REAL(KIND=8), DIMENSION(nx), INTENT(OUT) :: tk

! NCLEND

    INTEGER :: i

    !REAL(KIND=8), PARAMETER :: P1000MB=100000.D0, RD=287.D0, CP=7.D0*RD/2.D0

    !$OMP PARALLEL DO SCHEDULE(runtime)
    DO i = 1,nx
        !pi = (pressure(i)/P1000MB)**(RD/CP)
        !tk(i) = pi * theta(i)
        tk(i) = (pressure(i)/P1000MB)**(RD/CP) * theta(i)
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DCOMPUTETK


! NCLFORTSTART
SUBROUTINE DINTERP3DZ(data3d, out2d, zdata, levels, nx, ny, nz, nlev, missingval)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: out2d

    INTEGER, INTENT(IN) :: nx, ny, nz, nlev
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) ::  data3d
    REAL(KIND=8), DIMENSION(nx,ny,nlev), INTENT(OUT) :: out2d
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: zdata
    REAL(KIND=8), DIMENSION(nlev), INTENT(IN) :: levels
    REAL(KIND=8), INTENT(IN) :: missingval

! NCLEND

    INTEGER :: i,j,kp,ip,im,lev
    LOGICAL :: dointerp
    REAL(KIND=8) :: w1,w2,desiredloc

    ! does vertical coordinate increase or decrease with increasing k?
    ! set offset appropriately

    ip = 0
    im = 1
    IF (zdata(1,1,1) .GT. zdata(1,1,nz)) THEN
        ip = 1
        im = 0
    END IF

    !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,lev,kp,dointerp,w1,w2,desiredloc) &
    !$OMP FIRSTPRIVATE(ip,im) SCHEDULE(runtime)
    DO lev = 1,nlev
        DO i = 1,nx
            DO j = 1,ny
                ! Initialize to missing.  Was initially hard-coded to -999999.
                out2d(i,j,lev) = missingval
                dointerp = .FALSE.
                kp = nz
                desiredloc = levels(lev)

                DO WHILE ((.NOT. dointerp) .AND. (kp >= 2))
                    IF (((zdata(i,j,kp-im) < desiredloc) .AND. (zdata(i,j,kp-ip) > desiredloc))) THEN
                        w2 = (desiredloc - zdata(i,j,kp-im))/(zdata(i,j,kp-ip) - zdata(i,j,kp-im))
                        w1 = 1.D0 - w2
                        out2d(i,j,lev) = w1*data3d(i,j,kp-im) + w2*data3d(i,j,kp-ip)
                        dointerp = .TRUE.
                    END IF
                    kp = kp - 1
                END DO
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DINTERP3DZ


! NCLFORTSTART
SUBROUTINE DINTERP3DZ_2DLEV(data3d, out2d, zdata, levs2d, nx, ny, nz, missingval)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: out2d

    INTEGER, INTENT(IN) :: nx, ny, nz
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) ::  data3d
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(OUT) :: out2d
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: zdata
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: levs2d
    REAL(KIND=8), INTENT(IN) :: missingval

! NCLEND

    INTEGER :: i,j,kp,ip,im
    LOGICAL :: dointerp
    REAL(KIND=8) :: w1,w2,desiredloc

    ! does vertical coordinate increase or decrease with increasing k?
    ! set offset appropriately

    ip = 0
    im = 1
    IF (zdata(1,1,1) .GT. zdata(1,1,nz)) THEN
        ip = 1
        im = 0
    END IF

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,kp,dointerp,w1,w2,desiredloc) &
    !$OMP FIRSTPRIVATE(ip,im) SCHEDULE(runtime)
    DO i = 1,nx
        DO j = 1,ny
            ! Initialize to missing.  Was initially hard-coded to -999999.
            out2d(i,j) = missingval
            dointerp = .FALSE.
            kp = nz
            desiredloc = levs2d(i,j)

            DO WHILE ((.NOT. dointerp) .AND. (kp >= 2))
                IF (((zdata(i,j,kp-im) < desiredloc) .AND. (zdata(i,j,kp-ip) > desiredloc))) THEN
                    w2 = (desiredloc - zdata(i,j,kp-im))/(zdata(i,j,kp-ip) - zdata(i,j,kp-im))
                    w1 = 1.D0 - w2
                    out2d(i,j) = w1*data3d(i,j,kp-im) + w2*data3d(i,j,kp-ip)
                    dointerp = .TRUE.
                END IF
                kp = kp - 1
            END DO
        END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DINTERP3DZ_2DLEV

! PORT DZSTAG HERE

! NCLFORTSTART
SUBROUTINE DZSTAG(znew, nx, ny, nz, z, nxz, nyz ,nzz, terrain)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: znew

    INTEGER, INTENT(IN) :: nx, ny, nz, nxz, nyz, nzz
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(OUT) :: znew
    REAL(KIND=8), DIMENSION(nxz,nyz,nzz), INTENT(IN) :: z
    REAL(KIND=8), DIMENSION(nxz,nyz), INTENT(IN) :: terrain
! NCLEND

    INTEGER :: i,j,k,ii,im1,jj,jm1

    ! check for u, v, or w (x,y,or z) staggering
    ! for x and y stag, avg z to x, y, point
    IF (nx .GT. nxz) THEN
        DO k = 1,nz
            DO j = 1,ny
                DO i = 1,nx
                    ii = MIN(i,nxz)
                    im1 = MAX(i-1,1)
                    znew(i,j,k) = 0.5D0*(z(ii,j,k) + z(im1,j,k))
                END DO
            END DO
        END DO

    ELSE IF (ny .GT. nyz) THEN
        DO k = 1,nz
            DO j = 1,NY
                jj = MIN(j,nyz)
                jm1 = MAX(j-1,1)
                DO i = 1,nx
                    znew(i,j,k) = 0.5D0*(z(i,jj,k) + z(i,jm1,k))
                END DO
            END DO
        END DO

    ! w (z) staggering
    ELSE IF (nz .GT. nzz) THEN
        DO j = 1,ny
            DO i = 1,nx
                znew(i,j,1) = terrain(i,j)
            END DO
        END DO

        DO k = 2,nz
            DO j = 1,ny
                DO i = 1,nx
                    znew(i,j,k) = znew(i,j,k-1) + 2.D0*(z(i,j,k-1) - znew(i,j,k-1))
                END DO
            END DO
        END DO

    END IF

    RETURN

    END SUBROUTINE DZSTAG

! NCLFORTSTART
SUBROUTINE DINTERP2DXY(v3d, v2d, xy, nx, ny, nz, nxy)

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: v2d

    INTEGER, INTENT(IN) :: nx, ny, nz, nxy
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: v3d
    REAL(KIND=8), DIMENSION(nxy,nz), INTENT(OUT) :: v2d
    REAL(KIND=8), DIMENSION(2,nxy), INTENT(IN) :: xy

! NCLEND

    INTEGER :: i, j, k, ij
    REAL(KIND=8) :: w11, w12, w21, w22, wx, wy

    !$OMP PARALLEL DO PRIVATE(i,j,k,ij,w11,w12,w21,w22,wx,wy) &
    !$OMP SCHEDULE(runtime)
    DO ij = 1,nxy
        i = MAX(1,MIN(nx-1,INT(xy(1,ij)+1)))
        j = MAX(1,MIN(ny-1,INT(xy(2,ij)+1)))
        wx = DBLE(i+1) - (xy(1,ij)+1)
        wy = DBLE(j+1) - (xy(2,ij)+1)
        w11 = wx*wy
        w21 = (1.D0-wx)*wy
        w12 = wx*(1.D0-wy)
        w22 = (1.D0-wx)* (1.D0-wy)
        DO k = 1,nz
            v2d(ij,k) = w11*v3d(i,j,k) + w21*v3d(i+1,j,k) + &
            w12*v3d(i,j+1,k) + w22*v3d(i+1,j+1,k)
        END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DINTERP2DXY


! NCLFORTSTART
SUBROUTINE DINTERP1D(v_in, v_out, z_in, z_out, vmsg, nz_in, nz_out)

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: v_out

    INTEGER, INTENT(IN) :: nz_in, nz_out
    REAL(KIND=8), DIMENSION(nz_in), INTENT(IN) :: v_in, z_in
    REAL(KIND=8), DIMENSION(nz_out), INTENT(IN) :: z_out
    REAL(KIND=8), DIMENSION(nz_out), INTENT(OUT) :: v_out
    REAL(KIND=8), INTENT(IN) :: vmsg

! NCLEND

    INTEGER :: kp,k,im,ip
    LOGICAL :: interp
    REAL(KIND=8) :: height,w1,w2

    ! does vertical coordinate increase of decrease with increasing k?
    ! set offset appropriately

    ip = 0
    im = 1
    IF (z_in(1) .GT. z_in(nz_in)) THEN
        ip = 1
        im = 0
    END IF

    !$OMP PARALLEL DO PRIVATE(kp, k, interp, height, w1, w2) &
    !$OMP FIRSTPRIVATE(ip, im) SCHEDULE(runtime)
    DO k = 1,nz_out
        v_out(k) = vmsg

        interp = .FALSE.
        kp = nz_in
        height = z_out(k)

        DO WHILE ((.NOT. interp) .AND. (kp .GE. 2))
            IF (((z_in(kp-im) .LE. height) .AND. (z_in(kp-ip) .GT. height))) THEN
                w2 = (height - z_in(kp-im))/(z_in(kp-ip) - z_in(kp-im))
                w1 = 1.D0 - w2
                v_out(k) = w1*v_in(kp-im) + w2*v_in(kp-ip)
                interp = .TRUE.
            END IF
            kp = kp - 1
        END DO
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DINTERP1D

! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!
!    units: pressure (Pa), temperature(K), height (m), mixing ratio
!     (kg kg{-1}) availability of 3d p, t, and qv; 2d terrain; 1d
! half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!

! NCLFORTSTART
SUBROUTINE DCOMPUTESEAPRS(nx, ny, nz, z, t, p, q, sea_level_pressure, &
                          t_sea_level, t_surf, level, errstat, errmsg)
    USE wrf_constants, ONLY : ALGERR, RD, G, USSALR

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: sea_level_pressure

    !     Estimate sea level pressure.
    INTEGER, INTENT(IN) :: nx, ny, nz
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: z
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: t, p, q
    !     The output is the 2d sea level pressure.
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(OUT) :: sea_level_pressure
    INTEGER, DIMENSION(nx,ny), INTENT(INOUT) ::  level
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(INOUT) :: t_surf, t_sea_level
    INTEGER, INTENT(INOUT) :: errstat
    CHARACTER(LEN=*), INTENT(INOUT) :: errmsg

! NCLFORTEND

    !     Some required physical constants:

    !REAL(KIND=8), PARAMETER :: R=287.04D0, G=9.81D0, GAMMA=0.0065D0

    !     Specific constants for assumptions made in this routine:
    REAL(KIND=8), PARAMETER :: TC=273.16D0+17.5D0, PCONST=10000

    LOGICAL, PARAMETER :: ridiculous_mm5_test=.TRUE.
    !      PARAMETER (ridiculous_mm5_test = .FALSE.)

    !     Local variables:

    INTEGER :: i, j, k
    INTEGER :: klo, khi
    INTEGER :: errcnt, bad_i, bad_j

    REAL(KIND=8) :: bad_sfp
    REAL(KIND=8) :: plo, phi, tlo, thi, zlo, zhi
    REAL(KIND=8) :: p_at_pconst, t_at_pconst, z_at_pconst
    !REAL(KIND=8) :: z_half_lowest

    LOGICAL :: l1, l2, l3, found

    !  Find least zeta level that is PCONST Pa above the surface.  We
    !  later use this level to extrapolate a surface pressure and
    !  temperature, which is supposed to reduce the effect of the diurnal
    !  heating cycle in the pressure field.

    errstat = 0
    errcnt = 0
    bad_i = -1
    bad_j = -1
    bad_sfp = -1

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k,found) REDUCTION(+:errcnt) &
    !$OMP SCHEDULE(runtime)
    DO j = 1,ny
        DO i = 1,nx
            level(i,j) = -1

            k = 1
            found = .FALSE.
            DO WHILE ((.NOT. found) .AND. (k <= nz))
                IF (p(i,j,k) < p(i,j,1)-PCONST) THEN
                    level(i,j) = k
                    found = .TRUE.
                END IF
                k = k + 1
            END DO

            IF (level(i,j) == -1) THEN
                errcnt = errcnt + 1
                !$OMP CRITICAL
                ! Only do this the first time
                IF (bad_i .EQ. -1) THEN
                    bad_i = i
                    bad_j = j
                    bad_sfp = p(i,j,1) / 100.
                END IF
                !$OMP END CRITICAL
            END IF
        END DO
    END DO
    !$OMP END PARALLEL DO

    IF (errcnt > 0) THEN
        errstat = ALGERR
        WRITE(errmsg,*) "Error in finding 100 hPa up.  i=", bad_i, "j=", bad_j, "sfc_p=", bad_sfp
        RETURN
    END IF

    !     Get temperature PCONST Pa above surface.  Use this to extrapolate
    !     the temperature at the surface and down to sea level.
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,klo,khi,plo, &
    !$OMP phi,tlo,thi,zlo,zhi,p_at_pconst,t_at_pconst,z_at_pconst) &
    !$OMP REDUCTION(+:errcnt) SCHEDULE(runtime)
    DO j = 1,ny
        DO i = 1,nx

            klo = MAX(level(i,j)-1, 1)
            khi = MIN(klo+1, nz-1)

            IF (klo == khi) THEN
                errcnt = errcnt + 1
                !$OMP CRITICAL
                IF (bad_i .EQ. -1) THEN
                    bad_i = i
                    bad_j = j
                END IF
                !$OMP END CRITICAL
            END IF

            plo = p(i,j,klo)
            phi = p(i,j,khi)
            tlo = t(i,j,klo)*(1.D0 + 0.608D0*q(i,j,klo))
            thi = t(i,j,khi)*(1.D0 + 0.608D0*q(i,j,khi))
            zlo = z(i,j,klo)
            zhi = z(i,j,khi)
            p_at_pconst = p(i,j,1) - PCONST
            t_at_pconst = thi - (thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
            z_at_pconst = zhi - (zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
            !
            t_surf(i,j) = t_at_pconst * (p(i,j,1)/p_at_pconst)**(USSALR*RD/G)
            t_sea_level(i,j) = t_at_pconst + USSALR*z_at_pconst
        END DO
    END DO
    !$OMP END PARALLEL DO

    IF (errcnt > 0) THEN
        errstat = ALGERR
        WRITE(errmsg,*) "Error trapping levels at i=", bad_i, "j=", bad_j
        RETURN
    END IF

    ! If we follow a traditional computation, there is a correction to the
    ! sea level temperature if both the surface and sea level
    ! temperatures are *too* hot.

    IF (ridiculous_mm5_test) THEN
        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(l1,l2,l3) SCHEDULE(runtime)
        DO j = 1,ny
            DO i = 1,nx
                l1 = t_sea_level(i,j) < TC
                l2 = t_surf(i,j) <= TC
                l3 = .NOT. l1
                IF (l2 .AND. l3) THEN
                    t_sea_level(i,j) = TC
                ELSE
                    t_sea_level(i,j) = TC - 0.005D0*(t_surf(i,j)-TC)**2
                END IF
            END DO
        END DO
        !$OMP END PARALLEL DO
    END IF

    !     The grand finale: ta da!
    !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(runtime)
    DO j = 1,ny
        DO i = 1,nx
            !z_half_lowest = z(i,j,1)

    ! Convert to hPa in this step, by multiplying by 0.01. The original
    ! Fortran routine didn't do this, but the NCL script that called it
    ! did, so we moved it here.
            sea_level_pressure(i,j) = 0.01*(p(i,j,1)*EXP((2.D0*G*z(i,j,1))/&
                (RD*(t_sea_level(i,j) + t_surf(i,j)))))
        END DO
    END DO
    !$OMP END PARALLEL DO

    !     PRINT *,'sea pres input at weird location i=20,j=1,k=1'
    !     PRINT *,'t=',t(20,1,1),t(20,2,1),t(20,3,1)
    !     PRINT *,'z=',z(20,1,1),z(20,2,1),z(20,3,1)
    !     PRINT *,'p=',p(20,1,1),p(20,2,1),p(20,3,1)
    !     PRINT *,'slp=',sea_level_pressure(20,1),
    !    *         sea_level_pressure(20,2),sea_level_pressure(20,3)

    RETURN

END SUBROUTINE DCOMPUTESEAPRS


! Double precision version. If you make a change here, you
! must make the same change below to filter2d.

! NCLFORTSTART
SUBROUTINE DFILTER2D(a, b, nx, ny, it, missing, cenweight)

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: a

    INTEGER, INTENT(IN) :: nx, ny, it
    REAL(KIND=8), DIMENSION(nx, ny), INTENT(INOUT) :: a
    REAL(KIND=8), INTENT(IN) :: missing, cenweight
    REAL(KIND=8), DIMENSION(nx, ny), INTENT(INOUT) :: b

! NCLEND

    INTEGER :: i, j, iter
    REAL(KIND=8) :: cenmult, coef

    cenmult = (cenweight) / 2.
    coef = 1.0 / (4. + cenweight)

    DO iter=1,it
        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(runtime)
        DO j=1,ny
            DO i = 1,nx
                b(i,j) = a(i,j)
            END DO
        END DO
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(runtime)
        DO j=2,ny-1
            DO i=2,nx-1
                IF (b(i,j-1) .EQ. missing .OR. b(i,j) .EQ. missing .OR. &
                    b(i,j+1) .EQ. missing) THEN
                    a(i,j) = a(i,j)
                ELSE
                    a(i,j) = coef*(b(i,j-1) + cenmult*b(i,j) + b(i,j+1))
                END IF
            END DO
        END DO
        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(runtime)
        DO j=2,ny-1
            DO i=2,nx-1
                IF (b(i-1,j) .EQ. missing .OR. b(i,j) .EQ. missing .OR. &
                    b(i+1,j) .EQ. missing) THEN
                    a(i,j) = a(i,j)
                ELSE
                    a(i,j) = a(i,j) + coef*(b(i-1,j) + cenmult*b(i,j) + b(i+1,j))
                END IF
            END DO
        END DO
        !$OMP END PARALLEL DO
    !        do j=1,ny
    !        do i=1,nx
    !          b(i,j) = a(i,j)
    !        enddo
    !        enddo
    !        do j=2,ny-1
    !        do i=1,nx
    !          a(i,j) = a(i,j) - .99*coef*(b(i,j-1)-2*b(i,j)+b(i,j+1))
    !        enddo
    !        enddo
    !        do j=1,ny
    !        do i=2,nx-1
    !          a(i,j) = a(i,j) - .99*coef*(b(i-1,j)-2*b(i,j)+b(i+1,j))
    !        enddo
    !        enddo
    END DO

    RETURN

END SUBROUTINE DFILTER2D

! Single precision version. If you make a change here, you
! must make the same change below to dfilter2d.

! NCLFORTSTART
SUBROUTINE FILTER2D(a, b, nx, ny, it, missing, cenweight)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: a

    INTEGER, INTENT(IN) :: nx, ny, it
    REAL(KIND=4), DIMENSION(nx, ny), INTENT(INOUT) :: a
    REAL(KIND=4), INTENT(IN) :: missing, cenweight
    REAL(KIND=4), DIMENSION(nx, ny), INTENT(INOUT) :: b

! NCLEND

    !REAL(KIND=8), PARAMETER :: COEF = .125
    INTEGER :: i, j, iter
    REAL(KIND=8) :: cenmult, coef

    cenmult = (cenweight) / 2.
    coef = 1.0 / (4. + cenweight)

    !$OMP PARALLEL

    DO iter=1,it
        !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
        DO j=1,ny
            DO i = 1,nx
                b(i,j) = a(i,j)
            END DO
        END DO
        !$OMP END DO

        !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
        DO j=2,ny-1
            DO i=2,nx-1
                IF (b(i,j-1) .EQ. missing .OR. b(i,j) .EQ. missing .OR. &
                    b(i,j+1) .EQ. missing) THEN
                    a(i,j) = a(i,j)
                ELSE
                    a(i,j) = coef*(b(i,j-1) + cenmult*b(i,j) + b(i,j+1))
                END IF
            END DO
        END DO
        !$OMP END DO

        !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
        DO j=2,ny-1
            DO i=2,nx-1
                IF (b(i-1,j) .EQ. missing .OR. b(i,j) .EQ. missing .OR. &
                    b(i+1,j) .EQ. missing) THEN
                    a(i,j) = a(i,j)
                ELSE
                    a(i,j) = a(i,j) + coef*(b(i-1,j) + cenmult*b(i,j) + b(i+1,j))
                END IF
            END DO
        END DO
        !$OMP END DO
    !        do j=1,ny
    !        do i=1,nx
    !          b(i,j) = a(i,j)
    !        enddo
    !        enddo
    !        do j=2,ny-1
    !        do i=1,nx
    !          a(i,j) = a(i,j) - .99*coef*(b(i,j-1)-2*b(i,j)+b(i,j+1))
    !        enddo
    !        enddo
    !        do j=1,ny
    !        do i=2,nx-1
    !          a(i,j) = a(i,j) - .99*coef*(b(i-1,j)-2*b(i,j)+b(i+1,j))
    !        enddo
    !        enddo
    END DO

    !$OMP END PARALLEL

    RETURN

END SUBROUTINE FILTER2D


! NCLFORTSTART
SUBROUTINE DCOMPUTERH(qv, p, t, rh, nx)
    USE wrf_constants, ONLY : EZERO, ESLCON1, ESLCON2, CELKEL, RD, RV, EPS

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: rh

    INTEGER, INTENT(IN) :: nx
    REAL(KIND=8), DIMENSION(nx), INTENT(IN) :: qv, p, t
    REAL(KIND=8), DIMENSION(nx), INTENT(OUT) :: rh

! NCLEND

    INTEGER :: i
    REAL(KIND=8) :: qvs,es,pressure,temperature

    !$OMP PARALLEL DO PRIVATE(qvs, es, pressure, temperature) SCHEDULE(runtime)
    DO i = 1,nx
        pressure = p(i)
        temperature = t(i)
        ! es  = 1000.*svp1*
        es = EZERO*EXP(ESLCON1*(temperature - CELKEL)/(temperature - ESLCON2))
        ! qvs = ep_2*es/(pressure-es)
        qvs = EPS*es/(0.01D0*pressure - (1.D0 - EPS)*es)
        ! rh = 100*MAX(1., qv(i)/qvs)
        ! rh(i) = 100.*qv(i)/qvs
        rh(i) = 100.D0*MAX(MIN(qv(i)/qvs, 1.0D0), 0.0D0)
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DCOMPUTERH


! NCLFORTSTART
SUBROUTINE DGETIJLATLONG(lat_array, long_array, lat, longitude, ii, jj, nx, ny, imsg)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: ii,jj

    INTEGER, INTENT(IN) :: nx, ny, imsg
    INTEGER, INTENT(OUT) :: ii, jj
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: lat_array, long_array
    REAL(KIND=8) :: lat, longitude
! NCLEND

    REAL(KIND=8) :: longd,latd
    INTEGER :: i,j
    REAL(KIND=8) :: ir,jr
    REAL(KIND=8) :: dist_min,dist

    ! init to missing. was hard-coded to -999 initially.
    ir = imsg
    jr = imsg

    dist_min = 1.D+20
    DO j = 1,ny
        DO i = 1,nx
            latd = (lat_array(i,j) - lat)**2
            longd = (long_array(i,j) - longitude)**2
            ! longd = MIN((long_array(i,j)-longitude)**2, &
            !               (long_array(i,j)+longitude)**2)
            dist = SQRT(latd + longd)
            IF (dist_min .GT. dist) THEN
                dist_min = dist
                ir = DBLE(i)
                jr = DBLE(j)
            END IF
        END DO
    END DO

    ! The original version of this routine returned IR and JR. But, then
    ! the NCL script that called this routine was converting IR and JR
    ! to integer, so why not just return II and JJ?

    ! Also, I'm subtracing 1 here, because it will be returned to NCL
    ! script which has 0-based indexing.

    IF (ir .NE. imsg .AND. jr .NE. imsg) THEN
        ii = NINT(ir) - 1
        jj = NINT(jr) - 1
    ELSE
        ii = imsg
        jj = imsg
    END IF

    ! we will just return the nearest point at present

    RETURN

END SUBROUTINE DGETIJLATLONG

! You need to modify the C-WRAPPER in NCL for this to work

! NCLFORTSTART
SUBROUTINE DCOMPUTEUVMET(u, v, uvmet, longca,longcb,flong,flat, &
                        cen_long, cone, rpd, nx, ny, nxp1, nyp1, &
                        istag, is_msg_val, umsg, vmsg, uvmetmsg)
    IMPLICIT NONE

    ! ISTAG should be 0 if the U,V grids are not staggered.
    ! That is, NY = NYP1 and NX = NXP1.

    !f2py threadsafe
    !f2py intent(in,out) :: uvmet

    INTEGER,INTENT(IN) :: nx, ny, nxp1, nyp1, istag
    LOGICAL,INTENT(IN) :: is_msg_val
    REAL(KIND=8), DIMENSION(nxp1,ny), INTENT(IN):: u
    REAL(KIND=8), DIMENSION(nx,nyp1), INTENT(IN) :: v
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: flong
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(IN) :: flat
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(INOUT) :: longca
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(INOUT) :: longcb
    REAL(KIND=8), INTENT(IN) :: cen_long, cone, rpd
    REAL(KIND=8), INTENT(IN) :: umsg, vmsg, uvmetmsg
    REAL(KIND=8), DIMENSION(nx,ny,2), INTENT(OUT) :: uvmet


! NCLEND

    INTEGER :: i,j
    !REAL(KIND=8) :: uk, vk

    ! msg stands for missing value in this code
    ! WRITE (6,FMT=*) ' in compute_uvmet ',NX,NY,NXP1,NYP1,ISTAG

    !$OMP PARALLEL

    !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
    DO j = 1,ny
        DO i = 1,nx

            longca(i,j) = flong(i,j) - cen_long
            IF (longca(i,j).GT.180.D0) THEN
                longca(i,j) = longca(i,j) - 360.D0
            END IF
            IF (longca(i,j).LT.-180.D0) THEN
                longca(i,j) = longca(i,j) + 360.D0
            END IF
            IF (flat(i,j).LT.0.D0) THEN
                longcb(i,j) = -longca(i,j)*cone*rpd
            ELSE
                longcb(i,j) = longca(i,j)*cone*rpd
            END IF

            longca(i,j) = COS(longcb(i,j))
            longcb(i,j) = SIN(longcb(i,j))

        END DO
    END DO
    !$OMP END DO

    ! Note:  Intentionally removed as many IF statements as possible from loops
    ! to improve vectorization.

    IF (istag .EQ. 0) THEN ! Not staggered
        IF (.NOT. is_msg_val) THEN ! No missing values used
            !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
            DO j = 1,ny
                DO i = 1,nx
                    uvmet(i,j,1) = v(i,j)*longcb(i,j) + u(i,j)*longca(i,j)
                    uvmet(i,j,2) = v(i,j)*longca(i,j) - u(i,j)*longcb(i,j)
                END DO
            END DO
            !$OMP END DO
        ELSE ! Missing values used
            !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
            DO j = 1,ny
                DO i = 1,nx
                    IF ((u(i,j) .NE. umsg .AND. v(i,j) .NE. vmsg)) THEN
                        uvmet(i,j,1) = v(i,j)*longcb(i,j) + u(i,j)*longca(i,j)
                        uvmet(i,j,2) = v(i,j)*longca(i,j) - u(i,j)*longcb(i,j)
                    ELSE
                        uvmet(i,j,1) = uvmetmsg
                        uvmet(i,j,2) = uvmetmsg
                    END IF
                END DO
            END DO
            !$OMP END DO
        END IF
    ELSE ! Staggered
        IF (.NOT. is_msg_val) THEN ! No missing values used
            !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
            DO j = 1,ny
                DO i = 1,nx
                    ! This is the more readable version.
                    !uk = 0.5D0*(u(i,j) + u(i+1,j))
                    !vk = 0.5D0*(v(i,j) + v(i,j+1))
                    !uvmet(i,j,1) = vk*longcb(i,j) + uk*longca(i,j)
                    !uvmet(i,j,2) = vk*longca(i,j) - uk*longcb(i,j)

                    uvmet(i,j,1) = (0.5D0*(v(i,j) + v(i,j+1)))*longcb(i,j) + &
                                   (0.5D0*(u(i,j) + u(i+1,j)))*longca(i,j)
                    uvmet(i,j,2) = (0.5D0*(v(i,j) + v(i,j+1)))*longca(i,j) - &
                                   (0.5D0*(u(i,j) + u(i+1,j)))*longcb(i,j)

                END DO
            END DO
            !$OMP END DO
        ELSE ! Missing values used
            !$OMP DO COLLAPSE(2) SCHEDULE(runtime)
            DO j = 1,ny
                DO i = 1,nx
                    IF (u(i,j) .NE. umsg .AND. v(i,j) .NE. vmsg .AND. u(i+1,j) .NE. umsg .AND. v(i,j+1) .NE. vmsg) THEN
                        ! This is the more readable version.
                        !uk = 0.5D0*(u(i,j) + u(i+1,j))
                        !vk = 0.5D0*(v(i,j) + v(i,j+1))
                        !uvmet(i,j,1) = vk*longcb(i,j) + uk*longca(i,j)
                        !uvmet(i,j,2) = vk*longca(i,j) - uk*longcb(i,j)

                        uvmet(i,j,1) = (0.5D0*(v(i,j) + v(i,j+1)))*longcb(i,j) + &
                                       (0.5D0*(u(i,j) + u(i+1,j)))*longca(i,j)
                        uvmet(i,j,2) = (0.5D0*(v(i,j) + v(i,j+1)))*longca(i,j) - &
                                       (0.5D0*(u(i,j) + u(i+1,j)))*longcb(i,j)
                    ELSE
                        uvmet(i,j,1) = uvmetmsg
                        uvmet(i,j,2) = uvmetmsg
                    END IF
                END DO
            END DO
            !$OMP END DO
        END IF
    END IF

    !$OMP END PARALLEL

    RETURN

END SUBROUTINE DCOMPUTEUVMET




! This was originally a routine that took 2D input arrays. Since
! the NCL C wrapper routine can handle multiple dimensions, it's
! not necessary to have anything bigger than 1D here.

! NCLFORTSTART
SUBROUTINE DCOMPUTETD(td, pressure, qv_in, nx)
    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: td

    INTEGER, INTENT(IN) :: nx
    REAL(KIND=8), DIMENSION(nx), INTENT(IN) :: pressure
    REAL(KIND=8), DIMENSION(nx), INTENT(IN) :: qv_in
    REAL(KIND=8), DIMENSION(nx), INTENT(OUT) :: td

! NCLEND

    REAL(KIND=8) :: qv,tdc

    INTEGER :: i

    !$OMP PARALLEL DO PRIVATE(i,qv,tdc) SCHEDULE(runtime)
    DO i = 1,nx
        qv = MAX(qv_in(i), 0.D0)
         ! vapor pressure
         tdc = qv*pressure(i)/(.622D0 + qv)

         ! avoid problems near zero
         tdc = MAX(tdc, 0.001D0)
         td(i) = (243.5D0*LOG(tdc) - 440.8D0)/(19.48D0 - LOG(tdc))
    END DO
    !$OMP END PARALLEL DO

    RETURN

END SUBROUTINE DCOMPUTETD

! NCLFORTSTART
SUBROUTINE DCOMPUTEICLW(iclw, pressure, qc_in, nx, ny, nz)
    USE wrf_constants, ONLY : G

    IMPLICIT NONE

    !f2py threadsafe
    !f2py intent(in,out) :: iclw

    INTEGER, INTENT(IN) :: nx,ny,nz
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: pressure
    REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(IN) :: qc_in
    REAL(KIND=8), DIMENSION(nx,ny), INTENT(OUT) :: iclw

! NCLEND

    REAL(KIND=8) :: qclw, dp
    REAL(KIND=8), PARAMETER :: GG = 1000.D0/G
    INTEGER i,j,k

    iclw = 0

    DO j = 3,ny - 2
        DO i = 3,nx - 2
            DO k = 1,nz
                qclw = MAX(qc_in(i,j,k), 0.D0)
                IF (k.EQ.1) THEN
                    dp = pressure(i,j,k-1) - pressure(i,j,k)
                ELSE IF (k.EQ.nz) then
                    dp = pressure(i,j,k) - pressure(i,j,k+1)
                ELSE
                    dp = (pressure(i,j,k-1) - pressure(i,j,k+1))/2.D0
                END IF
                iclw(i,j) = iclw(i,j) + qclw*dp*GG
            END DO
        END DO
    END DO

    RETURN

END SUBROUTINE DCOMPUTEICLW










