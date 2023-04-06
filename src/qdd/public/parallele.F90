
!-----init_parallele-------------------------------------------------

SUBROUTINE init_parallele()

! Initializes parallele computing and determines actual node
! NUMBER 'myn'. The actual 'myn' is communicated via module 'params'.

  USE params
  IMPLICIT NONE
  INTEGER :: nprocs

!------------------------------------------------------------------

  myn = 0


! Set the number of OpenMP threads to use equal to the number of CPUs or,
! IF the OMP_NUM_THREADS env. var. is set, use that value.
! Also, the dynamic adjustment of the number of threads
! available for execution of parallel regions is enabled to prevent
! creating more threads than there are CPU's
#if(omp)
  IF (numthr > 0) setdyn = .FALSE.
!  CALL OMP_SET_DYNAMIC(setdyn)
  numthr = OMP_GET_MAX_THREADS()
!  CALL OMP_SET_NUM_THREADS(numthr) ! set the system wide number of OMP threads
  nthr = numthr - 1
#if(omp_debug)
  WRITE (*, *) 'setdyn = ', setdyn
  WRITE (*, *) 'numthr = ', numthr
  WRITE (*, *) 'OMP_GET_MAX_THREADS() = ', OMP_GET_MAX_THREADS()
  WRITE (*, *) 'OMP_GET_NUM_PROCS() = ', OMP_GET_NUM_PROCS()
  WRITE (*, *) 'OMP_GET_DYNAMIC() = ', OMP_GET_DYNAMIC()
  WRITE (*, *) 'OMP_GET_NESTED() = ', OMP_GET_NESTED()
!$OMP PARALLEL
  WRITE (*, *) 'OMP_GET_NUM_THREADS() = ', OMP_GET_NUM_THREADS()
  WRITE (*, *) 'OMP_GET_THREAD_NUM() = ', OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
#endif
  WRITE (*, *) 'INIT_PARALLELE(): numthr=', numthr
  WRITE (*, *) 'INIT_PARALLELE(): nthr=', nthr
#else
  nthr = 0
#endif

  RETURN
END SUBROUTINE init_parallele

