
! ***********************

SUBROUTINE calcpseudo()

! Driver routine for local part of pseudopotentials (PsP).
! Switches between soft Gaussian and Goedecker depending on
! 'ipsptyp' communicated via module 'params'.

  USE params
  USE kinetic
  IMPLICIT NONE
#if(raregas)
  INTEGER :: ind
#endif

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
#if(raregas)
    IF (idielec == 0) THEN
      CALL pseudosoft() ! soft Gaussian Psp
    ELSE
      CALL pseudosoft_dielec() ! Gaussian PsP with dielectric support
    END IF
#else
    CALL pseudosoft() ! soft Gaussian Psp
#endif

  CASE (1:4) ! Goedecker PsP (soft or local)
    CALL pseudogoed()

  CASE DEFAULT
    STOP ' CALCPSEUDO: this TYPE of PsP not yet implemented'
  END SELECT

#if(raregas)
  IF (idielec == 1) THEN
    ! Dielectric layer
    DO ind = 1, kdfull2
      rfieldtmp(ind) = potion(ind) ! potion will be erased in pseudosoft2
    END DO
    IF (ipsptyp == 0) THEN
      CALL pseudosoft2() !!! Not clear what happens here : seems to be similar to pseudosoft_dielec, but outside the layer (and slightly different formula)
    ELSE IF (ipsptyp >= 1) THEN
      STOP ' Goedecker PsP and dielectric layer incompatible'
    END IF
    DO ind = 1, kdfull2
      CALL conv1to3(ind)
      IF (iindtmp(1) > nint(xdielec/dx) + nx) THEN
        potion(ind) = rfieldtmp(ind) ! Outside the layer, take PsP as it was first computed
      END IF
    END DO
  END IF
#endif

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
#if(extended)
  INTEGER, EXTERNAL :: isoutofbox
#endif
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

#if(extended)
  IF (ipseudo) THEN

! traditional PsP of the Na cores by pseudodensities:

    pseudorho(1:nxyz) = 0D0

    DO ist = 1, nion

      IF (isoutofbox(cx(ist), cy(ist), cz(ist)) == 0) THEN

        ind = getnearestgridpoint(cx(ist), cy(ist), cz(ist))

        CALL conv1to3(ind)

        cfac1 = chg1(np(ist))/(pi**1.5D0*2D0**1.5D0*sgm1(np(ist))**3.D0)
        cfac2 = chg2(np(ist))/(pi**1.5D0*2D0**1.5D0*sgm2(np(ist))**3.D0)
        exfac1 = 1D0/(2D0*sgm1(np(ist))**2D0)
        exfac2 = 1D0/(2D0*sgm2(np(ist))**2D0)
        DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
          rz = (iz - nzsh)*dz - cz(ist)
          DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
            ry = (iy - nysh)*dy - cy(ist)
            DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg
              rx = (ix - nxsh)*dx - cx(ist)
              rr = rx*rx + ry*ry + rz*rz

              ind = conv3to1(ix, iy, iz)

              pseudorho(ind) = pseudorho(ind) &
                               + cfac1*EXP(-rr*exfac1) + cfac2*EXP(-rr*exfac2)

            END DO
          END DO
        END DO

      ELSE ! isOutOfBox not equal 0

        CALL addfunctofield1(potsave, v_soft, cx(ist), cy(ist), &
                             cz(ist), chg1(np(ist))*e2, sgm1(np(ist))*sq2)
        CALL addfunctofield1(potsave, v_soft, cx(ist), cy(ist), &
                             cz(ist), chg2(np(ist))*e2, sgm2(np(ist))*sq2)

      END IF ! isOutOfBox

    END DO ! sodium core loop

#if(raregas)
    IF (isurf /= 0) CALL pseudosoft_substrate(pseudorho, potsave)
#endif

    IF (tcoulfalr) THEN
      WRITE (*, *) 'IN PSEUDOSOFT FALR switch'
      CALL solv_poisson_f(pseudorho, potion, kdfull2)
      WRITE (*, *) 'IN PSEUDOSOFT after FALR switch'
    ELSE
      WRITE (*, *) 'IN PSEUDOSOFT COULEX switch'
      CALL solv_poisson_e(pseudorho, potion, kdfull2)
      WRITE (*, *) 'IN PSEUDOSOFT after COULEX switch'
    END IF

#if(raregas)
    IF (isurf /= 0 .AND. ivdw /= 2) THEN
      CALL addshortrepulsivepotonsubgrid(potshort, 1)
    ELSE IF (isurf /= 0 .AND. ivdw == 2) THEN
      CALL addshortrepulsivepot(potshort, 1)
    END IF
#endif

  ELSE ! ipseudo=.FALSE.
#endif

! pseudo-potentials directly

    DO ist = 1, nion

      CALL addfunctofield1(potion, v_soft, cx(ist), cy(ist), cz(ist), &
                           chg1(np(ist))*e2, sgm1(np(ist))*sq2)
      CALL addfunctofield1(potion, v_soft, cx(ist), cy(ist), cz(ist), &
                           chg2(np(ist))*e2, sgm2(np(ist))*sq2)

    END DO

#if(raregas)
    IF (isurf /= 0) THEN
      CALL addgsmpot(potion, 1)
      CALL addshortrepulsivepot(potshort, 1)
    END IF
#endif

#if(extended)
  END IF ! ipseudo
#endif

#if(raregas)
  IF (nc > 0 .AND. ivdw == 1) CALL getvdwpot
#endif

! now add contribution from fixed ions and out-of-box ions

  DO ind = 1, kdfull2
    potion(ind) = potion(ind) + potfixedion(ind) + potsave(ind) + potshort(ind)
#if(raregas)
    rfieldtmp(ind) = potshort(ind)
#endif
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
#if(raregas)
  REAL(DP), EXTERNAL:: v_soft
  REAL(DP), EXTERNAL:: v_ar_ar
  REAL(DP), EXTERNAL:: v_ar_na
#endif
!------------------------------------------------------------------------

  npmin = MIN(np(n1), np(n2))
  npmax = MAX(np(n1), np(n2))
#if(raregas)
  IF (npmin == -18) THEN
    IF (npmax == -18) THEN
      v_ion_ion = ch(-18)**2*e2*v_soft(dist, 2D0*sgm1(-18))
    ELSE IF (npmax == 18) THEN
      IF (ABS(n1 - n2) == nrare) THEN
        v_ion_ion = 0.5D0*c_dipmod*dist*dist
      ELSE
        v_ion_ion = ch(18)*ch(-18)*e2*v_soft(dist, 2D0*sgm1(18))
      END IF
    ELSE
      v_ion_ion = e2*ch(-18)*ch(npmax)*v_soft(dist, sgm1(18)*sq2)
    END IF
  ELSE IF (npmin /= 18) THEN
    IF (npmax == 18) THEN
      v_ion_ion = e2*ch(18)*ch(npmin)*v_soft(dist, sgm1(18)*sq2) + v_ar_na(dist)
    ELSE
      v_ion_ion = e2*ch(npmin)*ch(npmax)/dist
    END IF
  ELSE
    v_ion_ion = v_ar_ar(dist)
  END IF
#else
  v_ion_ion = e2*ch(npmin)*ch(npmax)/dist
#endif

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
