
#include "falr.F90"
#include "coulex.F90"


MODULE coulsolv

! Driver for Coulomb solver. Switches to FALR or exact.

  USE coulsolv_f
  USE coulsolv_e
  LOGICAL, PUBLIC :: tcoulfalr = .FALSE.

CONTAINS

  SUBROUTINE init_coul(dxi, dyi, dzi, nxi, nyi, nzi)

! Switch to actual initialization routine

    REAL(DP) :: dxi, dyi, dzi
    INTEGER ::nxi, nyi, nzi

    IF (tcoulfalr) THEN
      CALL init_coul_f(dxi, dyi, dzi, nxi, nyi, nzi)
! solv_poisson => solv_poisson_f
    ELSE
      CALL init_coul_e(dxi, dyi, dzi, nxi, nyi, nzi)
! solv_poisson => solv_poisson_e
    END IF

  END SUBROUTINE init_coul

  SUBROUTINE solv_poisson(rhoinp, chpfalr, kdum)

! Switch to actual solver routine

    REAL(8), INTENT(IN) :: rhoinp(*)
    REAL(8), INTENT(OUT) :: chpfalr(*)
    INTEGER, INTENT(IN) :: kdum ! dummy variable

    IF (tcoulfalr) THEN
      CALL solv_poisson_f(rhoinp, chpfalr, kdum)
    ELSE
      CALL solv_poisson_e(rhoinp, chpfalr, kdum)
    END IF

  END SUBROUTINE solv_poisson

END MODULE coulsolv

