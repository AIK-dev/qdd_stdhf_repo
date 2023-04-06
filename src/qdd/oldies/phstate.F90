
!-----phstate---------------------------------------------------------

SUBROUTINE phstate(psi)

! Picks one particular 1ph state from the given wavefunctions
! 'psi' and mixes this TO a rotated ph configuration. The state
! is rotated by 'phangle'.
! The particle state is 'npstate' and the hole state 'nhstate'.
! The old hole wavefunction is stored and returned on 'oldhole'.
! The parameters 'phangle', 'npstate' and 'nhstate' are communicated
! through MODULE PARAMS.
! 'oldhole' is also communicated through PARAMS and allocated here.

  USE params
  IMPLICIT REAL(DP) (A - H, O - Z)

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)

  COMPLEX(DP) :: csave
  INTEGER :: i

  IF (npstate > nstate) STOP ' PHSTATE: particle state OUT of range'
  IF (nhstate > nstate) STOP ' PHSTATE: hole state OUT of range'
  IF (occup(npstate) > 0.5D0) STOP 'PHSTATE: particle state already occupied'
  IF (occup(nhstate) < 0.5D0) STOP 'PHSTATE: hole state not occupied'

  ALLOCATE (oldhole(kdfull2))

  oldhole = psi(:, nhstate)
  ca = cos(phangle)
  sa = sin(phangle)
  psi(:, nhstate) = ca*psi(:, nhstate) + sa*psi(:, npstate)
  psi(:, npstate) = -sa*oldhole + ca*psi(:, npstate)

END SUBROUTINE phstate
