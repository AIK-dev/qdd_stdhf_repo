!-----density_criterion---------------------------------------------------------

SUBROUTINE density_criterion(iter, rho, aloc, denscrit)

!
! Change of density and potential as convergence criterion.
! This is the long version with 12 detailed criteria.
!
! Input:
! iter = iteration NUMBER at which routine is called
! rho = density
! aloc = local potentials for spin-up and spin-down
!
! Output:
! denscrit = array of criteria
! last index --> 1=rho, 2=aloc_up, 3=aloc_down.
! 4 = radius etc
! first index --> IN CASE of laste index 1-3
! 1=average |field-field_old|
! 2=average |field-field_old|**2
! 3= Max|field-field_old|
! IN CASE of last index 4
! 1= radius
! 2= quadrupole
! 3= dipole (rms)

  USE params

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iter
  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  REAL(DP), INTENT(OUT) :: denscrit(3, 4)

  REAL(DP) :: diffnorm, diffiter
  INTEGER, SAVE :: iter_old, nmin, nmax
  LOGICAL :: tfirst = .TRUE.
  REAL(DP), ALLOCATABLE, SAVE :: rho_old(:), aloc_old(:), qe_old(:)

  IF (tfirst) THEN
    ALLOCATE (rho_old(2*kdfull2))
    ALLOCATE (aloc_old(2*kdfull2))
    ALLOCATE (qe_old(2:7))
    rho_old = rho
    aloc_old = aloc
    iter_old = iter
    qe_old(2:7) = qe(2:7)
    denscrit = -1D0 ! TO indicate undefined staus
    tfirst = .FALSE.
    RETURN
  END IF

  diffiter = 1D0/(iter - iter_old)
  diffnorm = diffiter/kdfull2
  denscrit(1, 1) = SUM(ABS(rho(1:kdfull2) - rho_old(1:kdfull2)))*diffnorm
  denscrit(2, 1) = SQRT(SUM((rho(1:kdfull2) - rho_old(1:kdfull2))**2)*diffnorm)
  denscrit(3, 1) = MAXVAL(ABS(rho(1:kdfull2) - rho_old(1:kdfull2)))*diffiter

  nmin = 1
  nmax = kdfull2
  denscrit(1, 2) = SUM(ABS(aloc(nmin:nmax) - aloc_old(nmin:nmax)))*diffnorm
  denscrit(2, 2) = SQRT(SUM((aloc(nmin:nmax) - aloc_old(nmin:nmax))**2)*diffnorm)
  denscrit(3, 2) = MAXVAL(ABS(aloc(nmin:nmax) - aloc_old(nmin:nmax)))*diffiter

  nmin = kdfull2 + 1
  nmax = kdfull2*2
  denscrit(1, 3) = SUM(ABS(aloc(nmin:nmax) - aloc_old(nmin:nmax)))*diffnorm
  denscrit(2, 3) = SQRT(SUM((aloc(nmin:nmax) - aloc_old(nmin:nmax))**2)*diffnorm)
  denscrit(3, 3) = MAXVAL(ABS(aloc(nmin:nmax) - aloc_old(nmin:nmax)))*diffiter

  denscrit(1, 4) = (SQRT(SUM(qe(5:7))) - SQRT(SUM(qe_old(5:7))))*diffiter
  denscrit(2, 4) = (SQRT(SUM(qe(2:4)**2)) - SQRT(SUM(qe_old(2:4)**2)))*diffiter
  denscrit(3, 4) = ((2D0*qe(5) - qe(6) - qe(7)) - (2D0*qe_old(5) - qe_old(6) - qe_old(7))) &
                   *diffiter

  rho_old = rho
  aloc_old = aloc
  iter_old = iter
  qe_old(2:7) = qe(2:7)

  RETURN

END SUBROUTINE density_criterion
