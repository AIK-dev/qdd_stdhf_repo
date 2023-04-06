
! ***********************

SUBROUTINE calcpseudo()

! Driver routine for local part of pseudopotentials (PsP).
! Switches between soft Gaussian and Goedecker depending on
! 'ipsptyp' communicated via module 'params'.

  USE params
  USE kinetic
  IMPLICIT NONE

! Choice of ionic background :
  IF (nion2 == 2) THEN ! READ background from a FILE (potion.dat)
    CALL pseudo_external()
    RETURN
  ELSE IF (nion2 == 0) THEN ! Jellium (Homogeneous Electron Gas)
    RETURN
  END IF
! Or else, background from ionic PsP :

  SELECT CASE (ipsptyp)
  CASE (0) ! soft local PsP
    CALL pseudosoft() ! soft Gaussian Psp

  CASE (1:4) ! Goedecker PsP (soft or local)
    CALL pseudogoed()

  CASE DEFAULT
    STOP ' CALCPSEUDO: this TYPE of PsP not yet implemented'
  END SELECT


  RETURN
END SUBROUTINE calcpseudo

!------pseudo_external----------------------------------------------

SUBROUTINE pseudo_external()

! in this routine we read the PsP from FILE 'potion.dat'

  USE params
  IMPLICIT NONE

  OPEN (48, FORM='UNFORMATTED', FILE='potion.dat')
  READ (48) potion
  CLOSE (48)

  RETURN
END SUBROUTINE pseudo_external

!------pseudosoft----------------------------------------------

SUBROUTINE pseudosoft()

! Local pseudopotentials as sum of two Gaussians.
! I/O is handled via module 'params.

  USE params
  USE kinetic
  USE coulsolv, ONLY: solv_poisson_f, solv_poisson_e, tcoulfalr

  IMPLICIT NONE

! size of subgrid in units of mesh size

  INTEGER :: ind, ist, ix, iy, iz
  REAL(DP) :: cfac1, cfac2, exfac1, exfac2, rr, rx, ry, rz

  REAL(DP), DIMENSION(:), ALLOCATABLE :: pseudorho, potsave, potshort
  INTEGER, EXTERNAL :: conv3to1
  INTEGER, EXTERNAL :: getnearestgridpoint
  REAL(DP), EXTERNAL :: v_soft
!--------------------------------------------------------------

  IF (tfreezekspot .AND. tfs > 0) RETURN

  ALLOCATE (pseudorho(kdfull2))
  ALLOCATE (potsave(kdfull2))
  ALLOCATE (potshort(kdfull2))

  DO ind = 1, nxyz
    potion(ind) = 0D0
    potsave(ind) = 0D0
    potshort(ind) = 0D0
  END DO


! pseudo-potentials directly

    DO ist = 1, nion

      CALL addfunctofield1(potion, v_soft, cx(ist), cy(ist), cz(ist), &
                           chg1(np(ist))*e2, sgm1(np(ist))*sq2)
      CALL addfunctofield1(potion, v_soft, cx(ist), cy(ist), cz(ist), &
                           chg2(np(ist))*e2, sgm2(np(ist))*sq2)

    END DO




! now add contribution from fixed ions and out-of-box ions

  DO ind = 1, kdfull2
    potion(ind) = potion(ind) + potfixedion(ind) + potsave(ind) + potshort(ind)
  END DO

  DEALLOCATE (pseudorho)
  DEALLOCATE (potsave)
  DEALLOCATE (potshort)

  RETURN
END SUBROUTINE pseudosoft

!-----V_soft------------------------------------------------------------

REAL(DP) FUNCTION v_soft(r, sigma)
  USE params, ONLY: DP, PI
  IMPLICIT NONE

! Soft Coulomb potential from Gaussian density,
! Input:
! r = distance at which potential is computed
! sigma = width parameter of underlying Gaussian
  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: sigma

  REAL(DP) :: rabs
!------------------------------------------------------------------------

  rabs = ABS(r)
  v_soft = erf(rabs/sigma)/rabs

  RETURN
END FUNCTION v_soft

!-------------------------------------------------------------------------

REAL(DP) FUNCTION dvsdr(r, sigma)

! First derivative of V_soft

  USE params, ONLY: DP, pi
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: sigma

  REAL(DP), EXTERNAL :: v_soft
  REAL(DP):: fac, rabs

  rabs = ABS(r)

  fac = 2D0/(SQRT(pi)*sigma*rabs)

  dvsdr = fac*EXP(-rabs*rabs/(sigma*sigma)) - v_soft(rabs, sigma)/rabs

  RETURN
END FUNCTION dvsdr

REAL(DP) FUNCTION d2vsdr2(r, sigma)
  USE params, ONLY: DP, pi
  IMPLICIT NONE

! Second derivative of V_soft

  REAL(DP), INTENT(IN OUT) :: r
  REAL(DP), INTENT(IN) :: sigma

  REAL(DP) :: fac
  REAL(DP), EXTERNAL :: dvsdr
  REAL(DP), EXTERNAL :: v_soft
  r = ABS(r)

  fac = -2D0/(SQRT(pi)*sigma)

  d2vsdr2 = fac*EXP(-r*r/(sigma*sigma))*(r**(-2D0) + 2*sigma**(-2D0)) &
            + r**(-2D0)*v_soft(r, sigma) - dvsdr(r, sigma)/r

  RETURN
END FUNCTION d2vsdr2

!-----V_ion_ion------------------------------------------------------------

REAL(DP) FUNCTION v_ion_ion(dist, n1, n2)

  USE params
  USE kinetic
  IMPLICIT NONE

! effective ion-ion potential
! dist = distance at which potential is computed
! n1 = index of the first ion (1<=n1<=nion)
! n2 = index of the second ion (1<=n2<=nion)
! (see array np(*) in 'init.F')

  REAL(DP), INTENT(IN) :: dist
  INTEGER, INTENT(IN) :: n1
  INTEGER, INTENT(IN) :: n2

  INTEGER :: npmin, npmax
!------------------------------------------------------------------------

  npmin = MIN(np(n1), np(n2))
  npmax = MAX(np(n1), np(n2))
  v_ion_ion = e2*ch(npmin)*ch(npmax)/dist

  RETURN
END FUNCTION v_ion_ion

!------------------------------------------------------------

REAL(DP) FUNCTION dv_softdr(r, s)

! returns the derivative of erf(r/s)/r by finite differences

  USE params, ONLY: DP
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: s

  REAL(DP):: ftemp, rder
  REAL(DP), EXTERNAL::v_soft
  rder = 1.0D-5

  ftemp = v_soft(r + rder, s)
  ftemp = ftemp - v_soft(r - rder, s)
  ftemp = ftemp/(rder + rder)

  dv_softdr = ftemp

  RETURN
END FUNCTION dv_softdr
