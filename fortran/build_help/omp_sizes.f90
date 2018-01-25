PROGRAM sizes

    USE omp_lib

    PRINT *, "SIZES"
    PRINT *, INT(OMP_SCHED_KIND, kind=OMP_SCHED_KIND)
    PRINT *, INT(OMP_LOCK_KIND, kind=OMP_LOCK_KIND)
    PRINT *, INT(OMP_NEST_LOCK_KIND, kind=OMP_NEST_LOCK_KIND)
    PRINT *, INT(omp_sched_static, kind=OMP_SCHED_KIND)
    PRINT *, INT(omp_sched_dynamic, kind=OMP_SCHED_KIND)
    PRINT *, INT(omp_sched_guided, kind=OMP_SCHED_KIND)
    PRINT *, INT(omp_sched_auto, kind=OMP_SCHED_KIND)

END PROGRAM sizes
