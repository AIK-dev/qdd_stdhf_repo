
!------------------------------------------------------------

SUBROUTINE getforces(rho, psi, it, iflag)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! switchboard for calculating forces on all particles
! except DFT-electrons depending in 'ipsptype'

! rho = input electronic density
! psi = input s.p. wavefunctions
! it = iteration number
! iflag = flag parameters, used only for 'ipsptype=0'

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  INTEGER, INTENT(IN) :: it
  INTEGER, INTENT(IN) :: iflag

  CHARACTER(LEN=1) :: ext
  INTEGER :: i, ion

  WRITE (6, *) 'Entering getforces: ipsptyp=', ipsptyp, &
    ', tnonlocany=', tnonlocany

! in case of Goedecker PsP, switch to the corresponding routine
! in the old force.F FILE

  IF (ipsptyp == 1 .AND. tnonlocany) THEN
    CALL calcf_goenonl(rho, it, psi)
  ELSE IF ((ipsptyp == 1 .AND. .NOT. tnonlocany) .OR. ipsptyp == 2) THEN
    CALL calcf_goeloc(rho)
  ELSE
    CALL calcf_softlocal(rho, iflag)
  END IF

! force from colliding projectile
  CALL forceproject()

! force from laser
  CALL laserf()

! protocol

  IF (jforce /= 0 .AND. MOD(it, jforce) == 0 .AND. it >= 0) THEN
    IF (projcharge /= 0D0) THEN
      DO ion = 1, nion
        WRITE (ext, '(i1)') ion
        OPEN (8886, POSITION='append', FILE='projforce'//ext//'.'//outnam)
        WRITE (8886, '(4f13.5)') tfs, fprojx(ion), fprojy(ion), fprojz(ion)
      END DO
    END IF
    IF (jlaser > 0) THEN
      DO ion = 1, nion
        WRITE (ext, '(i1)') ion
        OPEN (255, POSITION='append', FILE='plforce.'//ext//'.'//outnam)
        WRITE (255, '(4f13.5)') tfs, flx(ion), fly(ion), flz(ion)
      END DO
    END IF
    DO ion = 1, nion
      WRITE (ext, '(i1)') ion
      OPEN (24, POSITION='append', FILE='pforce.'//ext//'.'//outnam)
      WRITE (24, '(4f13.5)') tfs, fx(ion), fy(ion), fz(ion)
    END DO
  END IF

  RETURN
END SUBROUTINE getforces
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE calcf_softlocal(rho, iflag)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! Calculates forces on all particles except DFT-electrons
! for the caseof soft-local PsP

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: iflag

  INTEGER :: i


  IF (iflag == 0) THEN
    fx(1:nion) = 0D0; fy(1:nion) = 0D0; fz(1:nion) = 0D0
  END IF


! pure cluster part

  IF (iflag == 0) THEN
    IF (nelect > 0) CALL getforceelna(rho)
    CALL getforcenana()
  END IF

! GSM relevant parts


  RETURN
END SUBROUTINE calcf_softlocal
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getforcenana()

! Force between ionic cores
! Input and output communicated via module 'params'

  USE params
  IMPLICIT NONE

  INTEGER :: ii, jj
  REAL(DP) :: dist, dist2, radfor, forcex, forcey, forcez
  REAL(DP) :: xi, yi, zi, xr, yr, zr

! Na(core)-Na(core) forces

  DO ii = 1, nion - 1
    xi = cx(ii)
    yi = cy(ii)
    zi = cz(ii)

    DO jj = ii + 1, nion

      xr = xi - cx(jj)
      yr = yi - cy(jj)
      zr = zi - cz(jj)

      dist2 = xr**2 + yr**2 + zr**2
      dist = SQRT(dist2)

      radfor = -e2*ch(np(ii))*ch(np(jj))/dist2

      forcex = radfor*xr/dist
      forcey = radfor*yr/dist
      forcez = radfor*zr/dist

      fx(ii) = fx(ii) - forcex
      fy(ii) = fy(ii) - forcey
      fz(ii) = fz(ii) - forcez

      fx(jj) = fx(jj) + forcex
      fy(jj) = fy(jj) + forcey
      fz(jj) = fz(jj) + forcez


    END DO

  END DO


  RETURN
END SUBROUTINE getforcenana
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getforceelna(rho)

! Force from electron cloud on ions.
!
! Input:
! rho = electronic local density
! Other I/O is communicated via module 'params'

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)

  INTEGER :: ii, ind
  REAL(DP) :: prefac

  EXTERNAL v_soft, gauss

  DO ii = 1, nion

      CALL foldgradfunc(rho, v_soft, cx(ii), cy(ii), cz(ii), sgm1(np(ii))*sq2)
! contribution from first gaussian
! the plus sign in the forces is really a plus because
! rho is positive but the electronic charge is -rho

      fx(ii) = fx(ii) + e2*chg1(np(ii))*rvectmp(1)
      fy(ii) = fy(ii) + e2*chg1(np(ii))*rvectmp(2)
      fz(ii) = fz(ii) + e2*chg1(np(ii))*rvectmp(3)

      CALL foldgradfunc(rho, v_soft, cx(ii), cy(ii), cz(ii), sgm2(np(ii))*sq2)
! contribution from second gaussian

      fx(ii) = fx(ii) + e2*chg2(np(ii))*rvectmp(1)
      fy(ii) = fy(ii) + e2*chg2(np(ii))*rvectmp(2)
      fz(ii) = fz(ii) + e2*chg2(np(ii))*rvectmp(3)




  END DO

  RETURN
END SUBROUTINE getforceelna
!------------------------------------------------------------

! ******************************

SUBROUTINE calcf_goeloc(rho)

! Force from local part of Goedecketr PsP on ions
! Input:
! rho = electronic local density
! it = time step nr. in calling routine
! Other I/O is communicated via module 'params'

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)

  CHARACTER(LEN=1) :: ext
  INTEGER :: ind, ion, ion1, is, ix, iy, iz
  REAL(DP) :: c1, c2, chpddr, dist2, dist3, pch
  REAL(DP) :: r2, rdn, rloc, zion
  REAL(DP) :: rr, rr3, rx, ry, rz, x1, y1, z1

  REAL(DP), EXTERNAL :: v_ion_el_lgoed

  REAL(DP), PARAMETER :: rder = 1D-1 ! radius step for finite difference

  DO ion = 1, nion
    fx(ion) = 0D0
    fy(ion) = 0D0
    fz(ion) = 0D0
  END DO

! derivatives of the n goedecker pseudos

  DO is = 1, nion
    c1 = cc1(np(is))
    c2 = cc2(np(is))
    rloc = crloc(np(is))*SQRT(2D0)
    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          IF ((ix /= nx2) .AND. (iy /= ny2) .AND. (iz /= nz2)) THEN
            rdn = 2D0/SQRT(pi)*e2
            zion = ch(np(is))*e2
            rx = x1 - cx(is)
            ry = y1 - cy(is)
            rz = z1 - cz(is)
            r2 = rx*rx + ry*ry + rz*rz + 1D-12
            rr = SQRT(r2)
            rr3 = rr*r2
            chpddr = -(v_ion_el_lgoed(rr + rder, rloc, c1, c2, zion) &
                       - v_ion_el_lgoed(rr - rder, rloc, c1, c2, zion))/(rder + rder)/rr
            fx(is) = fx(is) + (chpddr*rho(ind))*rx
            fy(is) = fy(is) + (chpddr*rho(ind))*ry
            fz(is) = fz(is) + (chpddr*rho(ind))*rz
          END IF
        END DO
      END DO
    END DO
  END DO
  DO ion = 1, nion
    fx(ion) = fx(ion)*dvol
    fy(ion) = fy(ion)*dvol
    fz(ion) = fz(ion)*dvol
  END DO

! force from ion-ion interactions

  DO ion = 1, nion
    DO ion1 = 1, nion
      IF (ion1 /= ion) THEN
        dist2 = (cx(ion1) - cx(ion))**2 + (cy(ion1) - cy(ion))**2 &
                + (cz(ion1) - cz(ion))**2
        dist3 = -SQRT(dist2)*dist2
        pch = e2*ch(np(ion))*ch(np(ion1))
        fx(ion) = fx(ion) - pch*(cx(ion) - cx(ion1))/dist3
        fy(ion) = fy(ion) - pch*(cy(ion) - cy(ion1))/dist3
        fz(ion) = fz(ion) - pch*(cz(ion) - cz(ion1))/dist3
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE calcf_goeloc

SUBROUTINE forceproject()

! Force from colliding projectile (IF any).
! Input:
! rho = electronic local density
! Other I/O is communicated via module 'params'

  USE params
  IMPLICIT NONE

  INTEGER :: ion
  REAL(DP) :: dist2, dist3, pch

  IF (ABS(projcharge) <= 1.0D-5) RETURN

  DO ion = 1, nion
    dist2 = (cx(ion) - (projvelx*tfs + projinix))**2 + &
            (cy(ion) - (projvely*tfs + projiniy))**2 + &
            (cz(ion) - (projvelz*tfs + projiniz))**2
    dist3 = -SQRT(dist2)*dist2
    pch = e2*ch(np(ion))*projcharge

    fprojx(ion) = -pch*(cx(ion) - (projvelx*tfs + projinix))/dist3
    fprojy(ion) = -pch*(cy(ion) - (projvely*tfs + projiniy))/dist3
    fprojz(ion) = -pch*(cz(ion) - (projvelz*tfs + projiniz))/dist3

    fx(ion) = fx(ion) + fprojx(ion)
    fy(ion) = fy(ion) + fprojy(ion)
    fz(ion) = fz(ion) + fprojz(ion)
  END DO

  RETURN

END SUBROUTINE forceproject

!----------------------------------------
SUBROUTINE laserf()
!----------------------------------------
! Force from EXTERNAL laser pulse on ionic cores.

! Since dyn_mfield is called previously in
! the routine dyn_propag (call of dyn_mfield then call of itstep
! that calls getforces), the array 'dataold' contains
! foft for 1st pulse (index 2), 2nd pulse (index 3) and XUV APT (index 4).
!
! Mind that foft is only the product of a cosine and a cosine square.
! Neither the field amplitude nor the polarization are included.

  USE params
  USE util, ONLY: laserp
  IMPLICIT NONE

  REAL(DP) :: tempx, tempy, tempz
  REAL(DP) :: foft1, foft2, foftXUV
  INTEGER :: i

  tempx = e0*e1x*foft1 + e0_2*e2x*foft2
  tempy = e0*e1y*foft1 + e0_2*e2y*foft2
  tempz = e0*e1z*foft1 + e0_2*e2z*foft2

  DO i = 1, nion
    flx(i) = -tempx*e2*ch(np(i))
    fly(i) = -tempy*e2*ch(np(i))
    flz(i) = -tempz*e2*ch(np(i))

    ! Adding forces from laser to total forces array
    fx(i) = fx(i) + flx(i)
    fy(i) = fy(i) + fly(i)
    fz(i) = fz(i) + flz(i)
  END DO


  RETURN
END SUBROUTINE laserf
! ********************************

SUBROUTINE calcf_goenonl(rho, it, psi)

! Force from non-local Goedecker PsP on ionic cores.
! Input:
! rho = electronic local density
! Other I/O is communicated via module 'params'

  USE params
  USE util, ONLY: realoverlap, realovsubgrid
  IMPLICIT NONE


  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

  INTEGER :: i, ion, ion1, nb
  REAL(DP) :: dist2, dist3, forceloc, forcenonl, pch
  REAL(DP) :: sumfor, sumfulp, sumfulm, sumslp, sumslm
  REAL(DP) :: xion, yion, zion
  REAL(DP) :: zshift = 0.001D0, yshift = 0.001D0, xshift = 0.001D0
  REAL(DP) :: zshinv = 1D3, yshinv = 1D3, xshinv = 1D3 ! must be inverse of zshift
  COMPLEX(DP), ALLOCATABLE :: q1(:)
  REAL(DP), ALLOCATABLE :: rhoslp(:), rhoslm(:)
  REAL(DP), ALLOCATABLE :: fxnl(:), fynl(:), fznl(:)

  CHARACTER(LEN=1) :: ext

  ALLOCATE (q1(kdfull2), rhoslp(kdfull2), rhoslm(kdfull2))
  ALLOCATE (fxnl(nion), fynl(nion), fznl(nion))

! force on the ions

  DO ion = 1, nion
    fx(ion) = 0D0
    fy(ion) = 0D0
    fz(ion) = 0D0
  END DO
  fznl = 0D0; fynl = 0D0; fxnl = 0D0

  DO ion = 1, nion

    xion = cx(ion)
    yion = cy(ion)
    zion = cz(ion)

    CALL rhopsg(xion, yion, zion + zshift*0.5D0, rhoslp, ion)
    CALL rhopsg(xion, yion, zion - zshift*0.5D0, rhoslm, ion)
    sumfor = 0D0
    DO i = 1, nxyz
      sumfor = sumfor + rho(i)*(rhoslp(i) - rhoslm(i))
    END DO
    fz(ion) = sumfor*dvol*zshinv

    IF (tnonloc(ion)) THEN
      CALL calc_proj(xion, yion, zion + zshift*0.5D0, xion, yion, zion, ion)
      sumslp = 0D0
      DO nb = 1, nstate
        CALL nonlocalc(psi(1, nb), q1, ion)
        sumslp = realovsubgrid(psi(1, nb), q1, ion) + sumslp
      END DO

      CALL calc_proj(xion, yion, zion - zshift*0.5D0, xion, yion, zion, ion)
      sumslm = 0D0
      DO nb = 1, nstate
        CALL nonlocalc(psi(1, nb), q1, ion)
        sumslm = realovsubgrid(psi(1, nb), q1, ion) + sumslm
      END DO
      fznl(ion) = (sumslm - sumslp)*zshinv
    END IF

    CALL rhopsg(xion, yion + yshift*0.5D0, zion, rhoslp, ion)
    CALL rhopsg(xion, yion - yshift*0.5D0, zion, rhoslm, ion)
    sumfor = 0D0
    DO i = 1, nxyz
      sumfor = sumfor + rho(i)*(rhoslp(i) - rhoslm(i))
    END DO
    fy(ion) = sumfor*dvol*yshinv

    IF (tnonloc(ion)) THEN
      CALL calc_proj(xion, yion + yshift*0.5D0, zion, xion, yion, zion, ion)
      sumslp = 0D0
      DO nb = 1, nstate
        CALL nonlocalc(psi(1, nb), q1, ion)
        sumslp = realovsubgrid(psi(1, nb), q1, ion) + sumslp
      END DO
      CALL calc_proj(xion, yion - yshift*0.5D0, zion, xion, yion, zion, ion)
      sumslm = 0D0
      DO nb = 1, nstate
        CALL nonlocalc(psi(1, nb), q1, ion)
        sumslm = realovsubgrid(psi(1, nb), q1, ion) + sumslm
      END DO
      fynl(ion) = (sumslm - sumslp)*yshinv
    END IF

    CALL rhopsg(xion + xshift*0.5D0, yion, zion, rhoslp, ion)
    CALL rhopsg(xion - xshift*0.5D0, yion, zion, rhoslm, ion)
    sumfor = 0D0
    DO i = 1, nxyz
      sumfor = sumfor + rho(i)*(rhoslp(i) - rhoslm(i))
    END DO
    fx(ion) = sumfor*dvol*xshinv
    forceloc = sumfor*xshinv ! for testing
    IF (tnonloc(ion)) THEN
      CALL calc_proj(xion + xshift*0.5D0, yion, zion, xion, yion, zion, ion)
      sumslp = 0D0
      sumfulp = 0D0
      DO nb = 1, nstate
        CALL nonlocalc(psi(1, nb), q1, ion)
        sumslp = realovsubgrid(psi(1, nb), q1, ion) + sumslp
        sumfulp = realoverlap(psi(1, nb), q1) + sumfulp
      END DO
      CALL calc_proj(xion - xshift*0.5D0, yion, zion, xion, yion, zion, ion)
      sumslm = 0D0
      sumfulm = 0D0
      DO nb = 1, nstate
        CALL nonlocalc(psi(1, nb), q1, ion)
        sumslm = realovsubgrid(psi(1, nb), q1, ion) + sumslm
        sumfulm = realoverlap(psi(1, nb), q1) + sumfulm
      END DO
      fxnl(ion) = (sumslm - sumslp)*xshinv
    END IF
    forcenonl = fxnl(ion)

    sumslm = sumslm*xshinv
    sumslp = sumslp*xshinv
    sumfulp = sumfulp*xshinv
    sumfulm = sumfulm*xshinv

! OPTIONAL prints

! restore projectors for original ionic positions
    CALL calc_proj(xion, yion, zion, xion, yion, zion, ion)

  END DO

! compose total force

  fx(1:nion) = fx(1:nion) + fxnl(1:nion)
  fy(1:nion) = fy(1:nion) + fynl(1:nion)
  fz(1:nion) = fz(1:nion) + fznl(1:nion)

  DEALLOCATE (fxnl, fynl, fznl)

  DO ion = 1, nion
    DO ion1 = 1, nion
      IF (ion1 /= ion) THEN
        dist2 = (cx(ion1) - cx(ion))**2 + (cy(ion1) - cy(ion))**2 &
                + (cz(ion1) - cz(ion))**2
        dist3 = -SQRT(dist2)*dist2
        pch = e2*ch(np(ion))*ch(np(ion1))
        fx(ion) = fx(ion) - pch*(cx(ion) - cx(ion1))/dist3
        fy(ion) = fy(ion) - pch*(cy(ion) - cy(ion1))/dist3
        fz(ion) = fz(ion) - pch*(cz(ion) - cz(ion1))/dist3
      END IF
    END DO
  END DO
  DEALLOCATE (q1, rhoslp, rhoslm)

  RETURN
END SUBROUTINE calcf_goenonl

!-----rhopsg-----------------------------------------------------------

! **********************************************

SUBROUTINE rhopsg(cxact, cyact, czact, rhopsp, is)

! Pseudo-density of ion 'is' related to a local PsP.
!
! Input:
! cxact,cyact,czact = coordinates of ion
! is = nr. of ion
! Output:
! rhopsp = emerging pseudo-density

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: is
  REAL(DP), INTENT(IN) :: cxact
  REAL(DP), INTENT(IN) :: cyact
  REAL(DP), INTENT(IN) :: czact
  REAL(DP), INTENT(OUT) :: rhopsp(kdfull2)

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: c1, c2, f1, f2, pt1, pt2, rloc, zion
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

  REAL(DP), EXTERNAL :: v_ion_el_lgoed
  REAL(DP), EXTERNAL :: v_soft

  DO ind = 1, nxyz
    rhopsp(ind) = 0D0
  END DO

! soft local PsP

  IF (ipsptyp == 0) THEN
    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          rx = x1 - cxact
          ry = y1 - cyact
          rz = z1 - czact
          rr = SQRT(rx*rx + ry*ry + rz*rz)
          rr = rr + 1D-10 ! avoid zero
          f1 = (v_soft(rr, (sgm1(np(is))*SQRT(2D0))))
          f2 = (v_soft(rr, (sgm2(np(is))*SQRT(2D0))))
          pt1 = e2*(f1)
          pt2 = e2*(f2)
          rhopsp(ind) = rhopsp(ind) + chg1(np(is))*pt1 + chg2(np(is))*pt2
        END DO
      END DO
    END DO

! local part of Goedecker

  ELSE IF (ipsptyp >= 1) THEN
    c1 = cc1(np(is))
    c2 = cc2(np(is))
    rloc = crloc(np(is))
    zion = ch(np(is))

    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          rx = x1 - cxact
          ry = y1 - cyact
          rz = z1 - czact
          rr = SQRT(rx*rx + ry*ry + rz*rz)
          rr = rr + 1D-6 ! avoid zero
          rhopsp(ind) = rhopsp(ind) + v_ion_el_lgoed(rr, rloc, c1, c2, zion)
        END DO
      END DO
    END DO

  ELSE
    STOP 'this type of PsP not yet provided'
  END IF

  RETURN
END SUBROUTINE rhopsg

REAL(DP) FUNCTION vgian(r)
  USE params, ONLY: DP
  IMPLICIT NONE

! returns the Gianozzi hydrogen PP in Rydberg units.
! Reference: F. Gygi, PRB 48, 11692 (1993).
! This PP was used for the calculations IN
! K"ummel, Kronik, Perdew, PRL 93, 213002 (2004)

  REAL(DP), INTENT(IN) :: r

  REAL(DP), PARAMETER :: rc1 = 0.25D0
  REAL(DP), PARAMETER :: rc2 = 0.284D0
  REAL(DP), PARAMETER :: a = -1.9287D0*2D0
  REAL(DP), PARAMETER :: b = 0.3374D0*2D0

  vgian = -2.d0*erf(r/rc1)/r + (a + b*r**2)*EXP(-(r/rc2)**2)

END FUNCTION vgian

!------------------------------------------------------------
!
REAL(DP) FUNCTION gauss(r, s)
!------------------------------------------------------------
  USE params, ONLY: DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: s
  REAL(DP)::rs


  rs = r/s
  gauss = EXP(-rs*rs)

  RETURN
END FUNCTION gauss
!------------------------------------------------------------

!------------------------------------------------------------
!
REAL(DP) FUNCTION dgaussdr(r, s)
!------------------------------------------------------------
! derivation of a Gaussian
  USE params, ONLY: DP
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: s
  REAL(DP)::ra

  ra = r/s
  dgaussdr = -2D0*ra*EXP(-(ra*ra))/s
  RETURN
END FUNCTION dgaussdr

