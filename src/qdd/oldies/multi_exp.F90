
SUBROUTINE cranknicolson_exp(q0, aloc, rho, it)

! Propagation with the Crank-Nicolson schem

! The problem is written as a linear system
! Ax=y
! WHERE y = [1-iH(t)dt/2-...]*psi(t); A= 1+iH(t+dt)dt/2+...(taylor series); x=psi(t+dt)
! at each time step:
! 1) compute y
! 2) predict H and psi at t+dt
! 3) compute the error Ax-y and verify that its norm tends TO 0.threshold criterion.
! other wise change x TO x-ERR.r; WHERE r is the chosen relaxation factor
! predict H(t+dt)
! DO this step until norm(Ax-y)> threshold THEN DO at the next time step.

!
! q0 = s.p. wavefunctions TO be propagated
! aloc = array for local mean field
! rho = array for local density
! it = nr. of time step
! qwork = work space for s.p. wavefunctions

  USE params
  USE util
  USE twost, ONLY: tnearest

  IMPLICIT NONE

  COMPLEX(DP), INTENT(INOUT) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(INOUT) :: aloc(2*kdfull2)
  REAL(DP), INTENT(INOUT) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it
! COMPLEX(DP), INTENT(OUT) :: qwork(kdfull2,kstate)

  INTEGER :: iter_cn, nb
  REAL(DP) :: ro, rnorm_err, relax, threshold1, threshold

  COMPLEX(DP) :: q1(kdfull2)
  COMPLEX(DP) :: cdtact
  COMPLEX(DP), DIMENSION(KDFULL2, KSTATE) :: y, errcn, x

! The PARAMETER 'tnorotate' de-activates the subtraction of the
! Lagrangian matrix IN the SIC step. The version of exponential
! evolution with subtraction of the Lagrangian matrix is found
! IN 'exp_evolp'. The strategy needs yet testing. It was not
! really beneficial so far.

  LOGICAL, PARAMETER :: tnorotate = .true.

#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
#endif

  cdtact = CMPLX(dt1/2D0, 0D0)

!----------------------------------------------------------------------

#if(mpi)
  STOP 'exponential evolution not yet MPI parallelized'
#else
  myn = 0
#endif

#if(raregas)
  IF (nc + NE + nk > 0) STOP 'TSTEP_EXP not appropriate for rare gas'
#endif

! one half time step TO define new mean field
! USE exponential evolution TO second order

  DO nb = 1, nstate
    CALL exp_evol(q0(1, nb), aloc, nb, 1, cdtact, q1)
  END DO
! SAVE the (1-iH(t+dt/2)dt/2-...)*psi(t) for the right term of Crank-Nicolson scheme

  y = q0

  CALL dyn_mfield(rho, aloc, y, dt1, it)

  DO nb = 1, nstate
    ! qwork(:,nb) = q0(:,nb)
    CALL exp_evol(q0(1, nb), aloc, nb, 1, cdtact, q1)
  END DO

! second half time step TO estimate x and H(t+dt)
!estimate mean field of full time step
!CALL dyn_mfield(rho,aloc,qwork,dt1,it)
  x = q0
! One negative half time step TO initialize Ax
!DO nb = 1,nstate
! CALL exp_evol(qwork(1,nb),aloc,nb,1,-cdtact,q1)
!END DO

!error
  rnorm_err = 1D0
  relax = 0.25D0
  threshold1 = 1D-2
  threshold = 1D-9
  iter_CN = 0

! begin of self consistent Crank Nicolson relaxation loop
  DO WHILE (rnorm_err .GT. threshold1)
    errcn = q0 - y
    ro = rnorm_err
    rnorm_err = 0D0
    DO nb = 1, 1
      rnorm_err = rnorm_err + wfnorm(errcn(:, nb))
    END DO
    IF (ro .LT. rnorm_err) THEN
      relax = relax*0.7D0
    ELSE
      relax = relax*1.2D0
    END IF
    x = x - relax*errcn
    q0 = x
    iter_CN = iter_CN + 1
    IF (iter_CN .GT. 20) THEN
      WRITE (*, *) ' warning Crank-Nicolson 1', iter_CN
      rnorm_err = 0D0
    ELSE
      IF (iter_CN .LT. 1) THEN
        rnorm_err = 1D0
      ELSE
        IF (MOD(iter_CN + 1, 30) .EQ. 0) CALL dyn_mfield(rho, aloc, x, dt1, it)

! One negative half time step TO compute Ax
        DO nb = 1, nstate
          CALL exp_evol(q0(1, nb), aloc, nb, 1, -cdtact, q1)
        END DO
      END IF
    END IF
  END DO
! END of first self consistent Crank Nicolson relaxation loop

  q0 = x

! compute mean field at new time
!CALL dyn_mfield(rho,aloc,qwork,dt1)
  RETURN
END SUBROUTINE cranknicolson_exp
