
#ifdef REALSWITCH

MODULE twost_util

  USE params
  USE kinetic
  USE orthmat
  USE util, ONLY: wfovlp, project

  INTEGER, SAVE, PRIVATE :: kdim

  INTEGER, SAVE, PRIVATE :: ndim(2)

  INTERFACE superpose_state
! automated choice between REAL, COMPLEX and complexsic versions of superpose_state
    MODULE PROCEDURE superpose_state_r, superpose_state_c, superpose_state_rc!, superpose_state_rc_dummy
  END INTERFACE superpose_state

CONTAINS
! ******************************

!----------------------------------------------------------------------
! REAL version
!----------------------------------------------------------------------
  SUBROUTINE superpose_state_r(wfsup, coeff, q0, is)

! Superposition TO new state:
! wfsup new single-particle state
! coeff superposition coefficients
! q0 set of s.p.states TO be combined
! is spin of states

    USE params
    USE kinetic
    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: wfsup(kdfull2)
    REAL(DP), INTENT(IN) :: coeff(kstate)
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    INTEGER, INTENT(IN) ::is

    INTEGER :: i, na, nas
!---------------------------------------------------------------------
    wfsup(1:nxyz) = 0.0_DP

    DO nas = 1, ndims(is)
      na = nas + (is - 1)*ndims(1)
      DO i = 1, nxyz
        wfsup(i) = wfsup(i) + coeff(nas)*q0(i, na)
      END DO
    END DO

    RETURN
  END SUBROUTINE superpose_state_r

!----------------------------------------------------------------------
! COMPLEX version
!----------------------------------------------------------------------
  SUBROUTINE superpose_state_c(wfsup, coeff, q0, is)

! Superposition TO new state:
! wfsup new single-particle state
! coeff superposition coefficients
! q0 set of s.p.states TO be combined
! is spin of states

    USE params
    USE kinetic
    IMPLICIT NONE

    COMPLEX(DP), INTENT(OUT) :: wfsup(kdfull2)
    COMPLEX(DP), INTENT(IN) :: coeff(kstate)
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
    INTEGER, INTENT(IN) ::is

    INTEGER :: i, na, nas
!---------------------------------------------------------------------

    wfsup(1:nxyz) = (0.0_DP, 0.0_DP)

    DO nas = 1, ndims(is)
      na = nas + (is - 1)*ndims(1)
      DO i = 1, nxyz
        wfsup(i) = wfsup(i) + coeff(nas)*q0(i, na)
      END DO
    END DO

    RETURN
  END SUBROUTINE superpose_state_c

!----------------------------------------------------------------------
! cmplxsic version
!----------------------------------------------------------------------
  SUBROUTINE superpose_state_rc(wfsup, coeff, q0, is)

! Superposition TO new state:
! wfsup new single-particle state
! coeff superposition coefficients
! q0 set of s.p.states TO be combined
! is spin of states

    USE params
    USE kinetic
    IMPLICIT NONE

    COMPLEX(DP), INTENT(OUT) :: wfsup(kdfull2)
    COMPLEX(DP), INTENT(IN) :: coeff(kstate)
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    INTEGER, INTENT(IN) ::is

    INTEGER :: i, na, nas
!---------------------------------------------------------------------

    wfsup(1:nxyz) = (0.0_DP, 0.0_DP)

    DO nas = 1, ndims(is)
      na = nas + (is - 1)*ndims(1)
      DO i = 1, nxyz
        wfsup(i) = wfsup(i) + coeff(nas)*q0(i, na)
      END DO
    END DO

    RETURN
  END SUBROUTINE superpose_state_rc

END MODULE twost_util

#endif

