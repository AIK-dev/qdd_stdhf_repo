!-----afterburn---------------------------------------------------------

SUBROUTINE afterburn(psir, rho, aloc)

! Improving static solution by doing a couple of steps with
! imaginary-time propatation.
!
! Input/Output:
! psir = REAL wavefunctions
! rho = electron density
! aloc = local mean-field potential

  USE params
  USE twost, ONLY: tnearest, init_fsic, init_vecs, end_fsic, expdabvol_rotate_init

  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: psir(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)

  COMPLEX(DP), ALLOCATABLE :: psiw(:, :), psi(:, :)

  INTEGER :: it, ion

! ****************************************************
! imaginary-time iteration (TO improve statics)
! ****************************************************

  WRITE (*, *) ' start afterburn. isitmax=', isitmax
  CALL flush (6)
  ALLOCATE (psi(kdfull2, kstate))
  IF (ifsicp > 7 .AND. nelect > 0) CALL init_fsic()
  IF (ifsicp >= 8) CALL expdabvol_rotate_init
  CALL restart2(psi, outnam, .true.) ! READ static wf's
  IF (ifsicp > 7) CALL init_vecs()
  CALL calcpseudo()
  CALL calclocal(rho, aloc) ! ??
  IF (ifsicp > 0) CALL calc_sic(rho, aloc, psi) ! USE some TYPE of SIC
  IF (ipsptyp == 1) THEN ! full Goedecker pseudopotentials require projectors
    DO ion = 1, nion
      IF (iswitch_interpol) THEN
        CALL calc_projFine(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
        CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
      ELSE
        CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
      END IF
    END DO
    IF (iswitch_interpol) CALL mergetabs
  END IF
  CALL info(psi, rho, aloc, -1)

  IF (ifexpevol) ALLOCATE (psiw(kdfull2, kstate))

  ! Start imaginary time iterations
  DO it = 1, isitmax
    WRITE (*, *) ' afterburn. iteration=', it
    IF (ifexpevol) THEN
      CALL tstep_exp(psi, aloc, rho, it, psiw, .TRUE.)
    ELSE
      STOP 'imaginary-time step requires exponential evolution'
! CALL tstep(psi,aloc,rho,it)
    END IF
    CALL cschmidt(psi) ! schmidt normalization of results
    IF (MOD(it, istinf) == 0) CALL info(psi, rho, aloc, it)
  END DO
  DEALLOCATE (psiw)
  CALL info(psi, rho, aloc, -1)
  CALL SAVE(psi, -1, outnam)
  DEALLOCATE (psi)
  IF (nelect > 0 .AND. ifsicp > 7) CALL end_fsic()
  IF (itmax <= 0) STOP ' terminate with afterburn '

  RETURN

END SUBROUTINE afterburn
