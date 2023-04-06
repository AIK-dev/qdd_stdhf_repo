SUBROUTINE V2ph(q0, rscreen, emin, emax)
  USE params
  IMPLICIT NONE

! Scans all 2ph states IN a given energy window, computes the
! matrix elements of the residual interaction, and prints them.

  COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate) ! set of s.p. wavefunctions
  REAL(DP), INTENT(IN) :: rscreen ! screening length
  REAL(DP), INTENT(IN) :: emin, emax ! energy window

  INTEGER :: n1, n2, n3, n4
  REAL(DP) :: edif, v1, v2, v2b
  COMPLEX(DP) :: V2phYuk, V2phdelta, cv2, cv2b

  WRITE (912, '(a,f7.4,i5/a)') '# 2ph matrix elements at time', time*0.0484, nstate, &
    ' n1 n2 n3 n4 edif VYukawa Vdelta'
  WRITE (913, '(a,f7.4/a)') '# 2ph matrix elements at time', time*0.0484, &
    ' n1 n2 n3 n4 edif VYukawa Vdelta Vdelta0'
  DO n1 = 1, nstate
    IF (occup(n1) < 0.9D0) CYCLE
    DO n2 = n1 + 1, nstate
      IF (occup(n2) < 0.9D0) CYCLE
      DO n3 = 1, nstate
        IF (occup(n3) > 0.1D0) CYCLE
        DO n4 = n3 + 1, nstate
          IF (occup(n4) > 0.1D0) CYCLE
          edif = amoy(n3) + amoy(n4) - amoy(n1) - amoy(n2)
          IF (edif < emin .OR. edif > emax) CYCLE
          v1 = ABS(V2phYuk(q0, n1, n2, n3, n4, rscreen))
          CALL delta(q0, n1, n2, n3, n4, rscreen, cv2, cv2b)
          v2 = ABS(cv2)
          v2b = ABS(cv2b)
          IF (v1 > 1D-10 .AND. v2 > 1D-10) &
            WRITE (912, '(8i4,3(1pg13.5),a,1pg13.5)') &
            n1, n2, n3, n4, ispin(n1), ispin(n2), &
            ispin(n3), ispin(n4), edif, v1, v2, ' . ', v2b
          IF (v1 > 1D-10 .AND. v2 == 0D0) &
            WRITE (912, '(8i4,1pg13.5,a,3(1pg13.5))') &
            n1, n2, n3, n4, ispin(n1), ispin(n2), &
            ispin(n3), ispin(n4), edif, ' 0 ', v2, v1, v2b
          IF (v1 > 1D-10) &
            WRITE (913, '(8i4,4(1pg13.5))') &
            n1, n2, n3, n4, ispin(n1), ispin(n2), &
            ispin(n3), ispin(n4), edif, v1, v2, v2b
          CALL flush (912)
          CALL flush (913)
        END DO
      END DO
    END DO
  END DO
  WRITE (912, '(1x/1x)')

END SUBROUTINE V2ph

!------------------------------------------------------------
COMPLEX(DP) FUNCTION V2phYuk(q0, n1, n2, n3, n4, rscreen)
!------------------------------------------------------------
  USE params
  USE kinetic
  IMPLICIT NONE

! Computes the 2ph matrix element between the four states n1...n4
! for the screened Coulomb interaction with screening length rscreen.
!
! Input:
! q0 = set of s.p. wavefunctions
! n1...n4 = the four s.p. states for which V2ph is computed
! rscreen = screening length

  COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate) ! set of s.p. wavefunctions
  REAL(DP), INTENT(IN) :: rscreen ! screening length
  INTEGER, INTENT(IN) :: n1, n2, n3, n4

  COMPLEX(DP), ALLOCATABLE :: crho(:), crhof(:)
  INTEGER :: ind
  REAL(DP) :: rsinv2
  COMPLEX(DP) :: vres1, vres2

  ALLOCATE (crho(kdfull2), crhof(kdfull2))

  IF (ispin(n1) == ispin(n3) .AND. ispin(n2) == ispin(n4)) THEN
    DO ind = 1, kdfull2
      crho(ind) = CONJG(q0(ind, n2))*q0(ind, n4)
    END DO

    rsinv2 = 1D0/rscreen**2
    CALL fftf(crho, crhof)
    crhof = 4D0*pi/(akv + rsinv2)*crhof
    CALL fftback(crhof, crho)

    vres1 = SUM(CONJG(q0(:, n1))*q0(:, n3)*crho)*dvol
  ELSE
    vres1 = (0D0, 0D0)
  END IF

  IF (ispin(n1) == ispin(n4) .AND. ispin(n2) == ispin(n3)) THEN
    DO ind = 1, kdfull2
      crho(ind) = CONJG(q0(ind, n2))*q0(ind, n3)
    END DO

    rsinv2 = 1D0/rscreen**2
    CALL fftf(crho, crhof)
    crhof = 4D0*pi/(akv + rsinv2)*crhof
    CALL fftback(crhof, crho)

    vres2 = -SUM(CONJG(q0(:, n1))*q0(:, n4)*crho)*dvol
  ELSE
    vres2 = (0D0, 0D0)
  END IF

  DEALLOCATE (crho, crhof)

  V2phYuk = vres1 + vres2

  RETURN
END FUNCTION V2phYuk

!------------------------------------------------------------
SUBROUTINE delta(q0, n1, n2, n3, n4, rscreen, V2phdelta, V2phbas)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! Computes the 2ph matrix element between the four states n1...n4
! for a zero-range interaction with strength equiavlent TO a
! screened Coulomb interaction
!
! Input:
! q0 = set of s.p. wavefunctions
! n1...n4 = the four s.p. states for which V2ph is computed
! rscreen = screening length

  COMPLEX(DP), INTENT(OUT) :: V2phdelta, V2phbas
  COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate) ! set of s.p. wavefunctions
  REAL(DP), INTENT(IN) :: rscreen ! screening length
  INTEGER, INTENT(IN) :: n1, n2, n3, n4

  REAL(DP) :: strength
  COMPLEX(DP) :: vres

  strength = 4D0*pi*rscreen**2

  IF (ispin(n1) == ispin(n3) .AND. ispin(n2) == ispin(n4) .AND. &
      ispin(n1) == ispin(n4)) THEN
    V2phdelta = 0D0
    V2phbas = SUM(CONJG(q0(:, n1)*q0(:, n2))*q0(:, n3)*q0(:, n4))*dvol*strength
  ELSE
    IF (ispin(n1) == ispin(n3) .AND. ispin(n2) == ispin(n4)) THEN
      vres = SUM(CONJG(q0(:, n1)*q0(:, n2))*q0(:, n3)*q0(:, n4))*dvol*strength
      V2phdelta = vres
      V2phbas = vres
    ELSE IF (ispin(n1) == ispin(n4) .AND. ispin(n2) == ispin(n3)) THEN
      V2phdelta = -SUM(CONJG(q0(:, n1)*q0(:, n2))*q0(:, n3)*q0(:, n4))*dvol*strength
      V2phbas = SUM(CONJG(q0(:, n1)*q0(:, n2))*q0(:, n3)*q0(:, n4))*dvol*strength
    ELSE
      V2phdelta = 0D0
      V2phbas = 0D0
    END IF
  END IF

  RETURN
END SUBROUTINE delta
