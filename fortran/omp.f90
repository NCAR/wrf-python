
MODULE omp_constants
    INTEGER, PARAMETER :: fomp_sched_kind = 4
    INTEGER, PARAMETER :: fomp_lock_kind = 4
    INTEGER, PARAMETER :: fomp_nest_lock_kind = 8
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_static = 1
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_dynamic = 2
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_guided = 3
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_auto = 4
contains
  subroutine have_omp_constants()
  end subroutine have_omp_constants
END MODULE omp_constants


FUNCTION fomp_enabled()

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_enabled

    fomp_enabled = .FALSE.

END FUNCTION fomp_enabled


SUBROUTINE fomp_set_num_threads(num_threads)

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: num_threads
    IF (.FALSE.) PRINT *, num_threads

END SUBROUTINE fomp_set_num_threads


FUNCTION fomp_get_num_threads()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_num_threads

    fomp_get_num_threads = -1

END FUNCTION fomp_get_num_threads


FUNCTION fomp_get_max_threads()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_max_threads

    fomp_get_max_threads = -1

END FUNCTION fomp_get_max_threads


FUNCTION fomp_get_thread_num()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_thread_num

    fomp_get_thread_num = -1

END FUNCTION fomp_get_thread_num


FUNCTION fomp_get_num_procs()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_num_procs

    fomp_get_num_procs = -1

END FUNCTION fomp_get_num_procs


FUNCTION fomp_in_parallel()

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_in_parallel

    fomp_in_parallel = .FALSE.

END FUNCTION fomp_in_parallel


SUBROUTINE fomp_set_dynamic(dynamic_threads)

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL, INTENT(IN) :: dynamic_threads
    IF (.FALSE.) PRINT *, dynamic_threads

END SUBROUTINE fomp_set_dynamic


FUNCTION fomp_get_dynamic()

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_get_dynamic

    fomp_get_dynamic = .FALSE.

END FUNCTION fomp_get_dynamic


SUBROUTINE fomp_set_nested(nested)

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL, INTENT(IN) :: nested
    IF (.FALSE.) PRINT *, nested

END SUBROUTINE fomp_set_nested


FUNCTION fomp_get_nested()

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_get_nested

    fomp_get_nested = .FALSE.

END FUNCTION fomp_get_nested


SUBROUTINE fomp_set_schedule(kind, modifier)

    USE omp_constants, ONLY : fomp_sched_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_sched_kind), INTENT(IN) :: kind
    INTEGER, INTENT(IN) :: modifier
    IF (.FALSE.) PRINT *, kind, modifier

END SUBROUTINE fomp_set_schedule


SUBROUTINE fomp_get_schedule(kind, modifier)

    USE omp_constants, ONLY : fomp_sched_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_sched_kind), INTENT(OUT) :: kind
    INTEGER, INTENT(OUT) :: modifier

    kind = -1
    modifier = -1

END SUBROUTINE fomp_get_schedule


FUNCTION fomp_get_thread_limit()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_thread_limit

    fomp_get_thread_limit = -1

END FUNCTION fomp_get_thread_limit


SUBROUTINE fomp_set_max_active_levels(max_levels)

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: max_levels
    IF (.FALSE.) PRINT *, max_levels

END SUBROUTINE fomp_set_max_active_levels


FUNCTION fomp_get_max_active_levels()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_max_active_levels

    fomp_get_max_active_levels = -1

END FUNCTION fomp_get_max_active_levels


FUNCTION fomp_get_level()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_level

    fomp_get_level = -1

END FUNCTION fomp_get_level


FUNCTION fomp_get_ancestor_thread_num(level)

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: level
    INTEGER :: fomp_get_ancestor_thread_num
    IF (.FALSE.) PRINT *, level

    fomp_get_ancestor_thread_num = -1

END FUNCTION fomp_get_ancestor_thread_num


FUNCTION fomp_get_team_size(level)

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: level
    INTEGER :: fomp_get_team_size
    IF (.FALSE.) PRINT *, level

    fomp_get_team_size = -1

END FUNCTION fomp_get_team_size


FUNCTION fomp_get_active_level()

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_active_level

    fomp_get_active_level = -1

END FUNCTION fomp_get_active_level


FUNCTION fomp_in_final()

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_in_final

    fomp_in_final = .FALSE.

END FUNCTION fomp_in_final


SUBROUTINE fomp_init_lock(svar)

    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(OUT) :: svar

    svar = -1

END SUBROUTINE fomp_init_lock


SUBROUTINE fomp_init_nest_lock(nvar)

    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(OUT) :: nvar

    nvar = -1

END SUBROUTINE fomp_init_nest_lock


SUBROUTINE fomp_destroy_lock(svar)

    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar
    IF (.FALSE.) PRINT *, svar

END SUBROUTINE fomp_destroy_lock


SUBROUTINE fomp_destroy_nest_lock(nvar)

    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar
    IF (.FALSE.) PRINT *, nvar

END SUBROUTINE fomp_destroy_nest_lock


SUBROUTINE fomp_set_lock(svar)

    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar
    IF (.FALSE.) PRINT *, svar

END SUBROUTINE fomp_set_lock


SUBROUTINE fomp_set_nest_lock(nvar)

    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar
    IF (.FALSE.) PRINT *, nvar

END SUBROUTINE fomp_set_nest_lock


SUBROUTINE fomp_unset_lock(svar)

    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar
    IF (.FALSE.) PRINT *, svar

END SUBROUTINE fomp_unset_lock


SUBROUTINE fomp_unset_nest_lock(nvar)

    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar
    IF (.FALSE.) PRINT *, nvar

END SUBROUTINE fomp_unset_nest_lock


FUNCTION fomp_test_lock(svar)

    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar
    LOGICAL :: fomp_test_lock
    IF (.FALSE.) PRINT *, svar

    fomp_test_lock = .FALSE.

END FUNCTION fomp_test_lock


FUNCTION fomp_test_nest_lock(nvar)

    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar
    INTEGER :: fomp_test_nest_lock
    IF (.FALSE.) PRINT *, nvar

    fomp_test_nest_lock = -1

END FUNCTION fomp_test_nest_lock


FUNCTION fomp_get_wtime()

    IMPLICIT NONE

    !f2py threadsafe

    REAL (KIND=8) :: fomp_get_wtime

    fomp_get_wtime = -1

END FUNCTION fomp_get_wtime


FUNCTION fomp_get_wtick()

    IMPLICIT NONE

    !f2py threadsafe

    REAL (KIND=8) :: fomp_get_wtick

    fomp_get_wtick = -1

END FUNCTION fomp_get_wtick


