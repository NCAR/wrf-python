MODULE omp_constants
#ifdef _OPENMP
    USE omp_lib
    INTEGER, PARAMETER :: fomp_sched_kind = omp_sched_kind
    INTEGER, PARAMETER :: fomp_nest_lock_kind = omp_nest_lock_kind
    INTEGER, PARAMETER :: fomp_lock_kind = omp_lock_kind
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_static = omp_sched_static
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_dynamic = omp_sched_dynamic
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_guided = omp_sched_guided
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_auto = omp_sched_auto
#else
    INTEGER, PARAMETER :: fomp_sched_kind = 4
    INTEGER, PARAMETER :: fomp_lock_kind = 4
    INTEGER, PARAMETER :: fomp_nest_lock_kind = 8
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_static = 1
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_dynamic = 2
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_guided = 3
    INTEGER(KIND=fomp_sched_kind), PARAMETER :: fomp_sched_auto = 4
#endif

END MODULE omp_constants


SUBROUTINE fomp_set_num_threads(num_threads)
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: num_threads

#ifdef _OPENMP
    CALL omp_set_num_threads(num_threads)
#endif

    RETURN

END SUBROUTINE fomp_set_num_threads


FUNCTION fomp_get_num_threads()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_num_threads

#ifdef _OPENMP
    fomp_get_num_threads = omp_get_num_threads()
#else
    fomp_get_num_threads = -1
#endif

END FUNCTION fomp_get_num_threads


FUNCTION fomp_get_max_threads()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_max_threads

#ifdef _OPENMP
    fomp_get_max_threads = omp_get_max_threads()
#else
    fomp_get_max_threads = -1
#endif

END FUNCTION fomp_get_max_threads


FUNCTION fomp_get_thread_num()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_thread_num

#ifdef _OPENMP
    fomp_get_thread_num = omp_get_thread_num()
#else
    fomp_get_thread_num = -1
#endif

END FUNCTION fomp_get_thread_num


FUNCTION fomp_get_num_procs()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_num_procs

#ifdef _OPENMP
    fomp_get_num_procs = omp_get_num_procs()
#else
    fomp_get_num_procs = -1
#endif

END FUNCTION fomp_get_num_procs


FUNCTION fomp_in_parallel()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_in_parallel

#ifdef _OPENMP
    fomp_in_parallel = omp_in_parallel()
#else
    fomp_in_parallel = .FALSE.
#endif

END FUNCTION fomp_in_parallel


SUBROUTINE fomp_set_dynamic(dynamic_threads)
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL, INTENT(IN) :: dynamic_threads

#ifdef _OPENMP
    CALL omp_set_dynamic(dynamic_threads)
#endif

    RETURN

END SUBROUTINE fomp_set_dynamic


FUNCTION fomp_get_dynamic()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_get_dynamic

#ifdef _OPENMP
    fomp_get_dynamic = omp_get_dynamic()
#else
    fomp_get_dynamic = .FALSE.
#endif

END FUNCTION fomp_get_dynamic


SUBROUTINE fomp_set_nested(nested)
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL, INTENT(IN) :: nested

#ifdef _OPENMP
    CALL omp_set_nested(nested)
#endif

    RETURN

END SUBROUTINE fomp_set_nested


FUNCTION fomp_get_nested()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_get_nested

#ifdef _OPENMP
    fomp_get_nested = omp_get_nested()
#else
    fomp_get_nested = .FALSE.
#endif

END FUNCTION fomp_get_nested


SUBROUTINE fomp_set_schedule(kind, modifier)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_sched_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_sched_kind), INTENT(IN) :: kind
    INTEGER, INTENT(IN) :: modifier

#ifdef _OPENMP
    CALL omp_set_schedule(kind, modifier)
#endif

    RETURN

END SUBROUTINE fomp_set_schedule


SUBROUTINE fomp_get_schedule(kind, modifier)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_sched_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_sched_kind), INTENT(OUT) :: kind
    INTEGER, INTENT(OUT) :: modifier

#ifdef _OPENMP
    CALL omp_get_schedule(kind, modifier)
#else
    kind = -1
    modifier = -1
#endif

    RETURN

END SUBROUTINE fomp_get_schedule


FUNCTION fomp_get_thread_limit()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_thread_limit

#ifdef _OPENMP
    fomp_get_thread_limit = omp_get_thread_limit()
#else
    fomp_get_thread_limit = -1
#endif

END FUNCTION fomp_get_thread_limit


SUBROUTINE fomp_set_max_active_levels(max_levels)
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

#ifdef _OPENMP
    INTEGER, INTENT(IN) :: max_levels
#else
    INTEGER, INTENT(IN) :: max_levels

    max_levels = -1
#endif

#ifdef _OPENMP
    CALL omp_set_max_active_levels(max_levels)
#endif

    RETURN

END SUBROUTINE fomp_set_max_active_levels


FUNCTION fomp_get_max_active_levels()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_max_active_levels

#ifdef _OPENMP
    fomp_get_max_active_levels = omp_get_max_active_levels()
#else
    fomp_get_max_active_levels = -1
#endif

END FUNCTION fomp_get_max_active_levels


FUNCTION fomp_get_level()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_level

#ifdef _OPENMP
    fomp_get_level = omp_get_level()
#else
    fomp_get_level = -1
#endif

END FUNCTION fomp_get_level


FUNCTION fomp_get_ancestor_thread_num(level)
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: level
    INTEGER :: fomp_get_ancestor_thread_num

#ifdef _OPENMP
    fomp_get_ancestor_thread_num = omp_get_ancestor_thread_num(level)
#else
    fomp_get_ancestor_thread_num = -1
#endif

END FUNCTION fomp_get_ancestor_thread_num


FUNCTION fomp_get_team_size(level)
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER, INTENT(IN) :: level
    INTEGER :: fomp_get_team_size

#ifdef _OPENMP
    fomp_get_team_size = omp_get_team_size(level)
#else
    fomp_get_team_size = -1
#endif

END FUNCTION fomp_get_team_size


FUNCTION fomp_get_active_level()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER :: fomp_get_active_level

#ifdef _OPENMP
    fomp_get_active_level = omp_get_active_level()
#else
    fomp_get_active_level = -1
#endif

END FUNCTION fomp_get_active_level


FUNCTION fomp_in_final()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    LOGICAL :: fomp_in_final

#ifdef _OPENMP
    fomp_in_final = omp_in_final()
#else
    fomp_in_final = .FALSE.
#endif

END FUNCTION fomp_in_final


SUBROUTINE fomp_init_lock(svar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(OUT) :: svar

#ifdef _OPENMP
    CALL omp_init_lock(svar)
#endif

    RETURN

END SUBROUTINE fomp_init_lock


SUBROUTINE fomp_init_nest_lock(nvar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(OUT) :: nvar

#ifdef _OPENMP
    CALL omp_init_nest_lock(nvar)
#endif

    RETURN

END SUBROUTINE fomp_init_nest_lock


SUBROUTINE fomp_destroy_lock(svar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar

#ifdef _OPENMP
    CALL omp_destroy_lock(svar)
#endif

    RETURN

END SUBROUTINE fomp_destroy_lock


SUBROUTINE fomp_destroy_nest_lock(nvar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar

#ifdef _OPENMP
    CALL omp_destroy_nest_lock(nvar)
#endif

    RETURN

END SUBROUTINE fomp_destroy_nest_lock


SUBROUTINE fomp_set_lock(svar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar

#ifdef _OPENMP
    CALL omp_set_lock(svar)
#endif

    RETURN

END SUBROUTINE fomp_set_lock


SUBROUTINE fomp_set_nest_lock(nvar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar

#ifdef _OPENMP
    CALL omp_set_nest_lock(nvar)
#endif

    RETURN

END SUBROUTINE fomp_set_nest_lock


SUBROUTINE fomp_unset_lock(svar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar

#ifdef _OPENMP
    CALL omp_unset_lock(svar)
#endif

    RETURN

END SUBROUTINE fomp_unset_lock


SUBROUTINE fomp_unset_nest_lock(nvar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar

#ifdef _OPENMP
    CALL omp_unset_nest_lock(nvar)
#endif

    RETURN

END SUBROUTINE fomp_unset_nest_lock


FUNCTION fomp_test_lock(svar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_lock_kind), INTENT(INOUT) :: svar
    LOGICAL :: fomp_test_lock

#ifdef _OPENMP
    fomp_test_lock =  omp_test_lock(svar)
#else
    fomp_test_lock = .FALSE.
#endif

    RETURN

END FUNCTION fomp_in_final


FUNCTION fomp_test_nest_lock(nvar)
#ifdef _OPENMP
    USE omp_lib
#endif
    USE omp_constants, ONLY : fomp_nest_lock_kind

    IMPLICIT NONE

    !f2py threadsafe

    INTEGER(KIND=fomp_nest_lock_kind), INTENT(INOUT) :: nvar
    INTEGER :: fomp_test_nest_lock

#ifdef _OPENMP
    fomp_test_lock =  omp_test_nest_lock(nvar)
#else
    fomp_test_lock = -1
#endif

    RETURN

END FUNCTION fomp_test_nest_lock


FUNCTION fomp_get_wtime()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    REAL (KIND=8) :: fomp_get_wtime

#ifdef _OPENMP
    fomp_get_wtime =  omp_get_wtime()
#else
    fomp_get_wtime = -1
#endif

    RETURN

END FUNCTION fomp_get_wtime


FUNCTION fomp_get_wtick()
#ifdef _OPENMP
    USE omp_lib
#endif

    IMPLICIT NONE

    !f2py threadsafe

    REAL (KIND=8) :: fomp_get_wtick

#ifdef _OPENMP
    fomp_get_wtick =  omp_get_wtick()
#else
    fomp_get_wtick = -1
#endif

    RETURN

END FUNCTION fomp_get_wtick


