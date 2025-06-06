! These are chosen to match the wrf module_model_constants.F where
! applicable
MODULE wrf_constants
    INTEGER, PARAMETER :: ERRLEN=512
    INTEGER, PARAMETER :: ALGERR=64

    REAL(KIND=8), PARAMETER :: WRF_EARTH_RADIUS = 6370000.D0
    REAL(KIND=8), PARAMETER :: T_BASE = 300.0D0
    REAL(KIND=8), PARAMETER :: PI = 3.1415926535897932384626433D0
    REAL(KIND=8), PARAMETER :: RAD_PER_DEG = PI/180.D0
    REAL(KIND=8), PARAMETER :: DEG_PER_RAD = 180.D0/PI
    REAL(KIND=8), PARAMETER :: DEFAULT_FILL = 9.9692099683868690D36
    INTEGER(KIND=1), PARAMETER :: DEFAULT_FILL_INT8 = -127
    INTEGER(KIND=2), PARAMETER :: DEFAULT_FILL_INT16 = -32767
    INTEGER(KIND=4), PARAMETER :: DEFAULT_FILL_INT32 = -2147483647
    INTEGER(KIND=8), PARAMETER :: DEFAULT_FILL_INT64 = INT(-9223372036854775806D0, KIND=8)
    REAL(KIND=4), PARAMETER :: DEFAULT_FILL_FLOAT = 9.9692099683868690E36
    REAL(KIND=8), PARAMETER :: DEFAULT_FILL_DOUBLE = 9.9692099683868690D36
    CHARACTER(LEN=1), PARAMETER :: DEFAULT_FILL_CHAR = ACHAR(0)


    REAL(KIND=8), PARAMETER :: P1000MB = 100000.D0
    ! j/k/kg
    REAL(KIND=8), PARAMETER :: RD = 287.D0
    REAL(KIND=8), PARAMETER :: RV = 461.6D0
    !REAL(KIND=8), PARAMETER :: RV = 461.5D0
    !  j/k/kg  note: not using bolton's value of 1005.7
    REAL(KIND=8), PARAMETER :: CP = 1004.5D0
    !REAL(KIND=8), PARAMETER :: CP = 1004D0

    REAL(KIND=8), PARAMETER :: G = 9.81D0
    REAL(KIND=8), PARAMETER :: USSALR = 0.0065D0  ! deg C per m

    REAL(KIND=8), PARAMETER :: CELKEL = 273.15D0
    REAL(KIND=8), PARAMETER :: CELKEL_TRIPLE = 273.16D0
    !REAL(KIND=8), PARAMETER :: GRAV = 9.81D0
    ! hpa
    REAL(KIND=8), PARAMETER :: EZERO = 6.112D0
    REAL(KIND=8), PARAMETER :: ESLCON1 = 17.67D0
    REAL(KIND=8), PARAMETER :: ESLCON2 = 29.65D0
    REAL(KIND=8), PARAMETER :: EPS = 0.622D0
    REAL(KIND=8), PARAMETER :: GAMMA = RD/CP
    !  cp_moist=cp*(1.+cpmd*qvp)
    REAL(KIND=8), PARAMETER :: CPMD = .887D0
    !  rgas_moist=rgas*(1.+rgasmd*qvp)
    REAL(KIND=8), PARAMETER :: RGASMD = .608D0
    !  gamma_moist=gamma*(1.+gammamd*qvp)
    REAL(KIND=8), PARAMETER :: GAMMAMD = RGASMD - CPMD
    REAL(KIND=8), PARAMETER :: TLCLC1 = 2840.D0
    REAL(KIND=8), PARAMETER :: TLCLC2 = 3.5D0
    REAL(KIND=8), PARAMETER :: TLCLC3 = 4.805D0
    REAL(KIND=8), PARAMETER :: TLCLC4 = 55.D0
    !  k
    REAL(KIND=8), PARAMETER :: THTECON1 = 3376.D0
    REAL(KIND=8), PARAMETER :: THTECON2 = 2.54D0
    REAL(KIND=8), PARAMETER :: THTECON3 = .81D0

    REAL(KIND=8), PARAMETER :: ABSCOEFI = .272D0  ! cloud ice absorption coefficient in m^2/g
    REAL(KIND=8), PARAMETER :: ABSCOEF = .145D0   ! cloud water absorption coefficient in m^2/g

    REAL(KIND=8), PARAMETER :: GAMMA_SEVEN = 720.D0
    REAL(KIND=8), PARAMETER :: RHOWAT = 1000.D0
    REAL(KIND=8), PARAMETER :: RHO_R = RHOWAT
    REAL(KIND=8), PARAMETER :: RHO_S = 100.D0
    REAL(KIND=8), PARAMETER :: RHO_G = 400.D0
    REAL(KIND=8), PARAMETER :: ALPHA = 0.224D0

    REAL(KIND=8), PARAMETER :: SCLHT = RD*256.D0/G
    REAL(KIND=8), PARAMETER :: EXPON =  RD*USSALR/G
    REAL(KIND=8), PARAMETER :: EXPONI =  1./EXPON

  contains
    subroutine have_wrf_constants()
    end subroutine have_wrf_constants

END MODULE wrf_constants

