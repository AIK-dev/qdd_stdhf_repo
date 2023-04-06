
! ***********************

SUBROUTINE pseudogoed()

! Local part of Goedecker pseudo-potentials.
! I/O handled through module 'params'.

  USE params
  IMPLICIT NONE

  INTEGER :: ind, is, ix, iy, iz
  REAL(DP) :: c1, c2, rloc, zion
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1
  REAL(DP), EXTERNAL :: vgian
  REAL(DP), EXTERNAL :: v_ion_el_lgoed

! REAL space part of the pseudopotentials taken IN
! PRB 54(3)1703 (1996) by Goedecker et al.

! the local part of the pseudos

  DO ind = 1, nxyz
    potion(ind) = 0D0
  END DO

  DO is = 1, nion
    c1 = cc1(np(is))
    c2 = cc2(np(is))
    rloc = crloc(np(is))
    zion = ch(np(is))
    IF (ipsptyp == 2 .AND. np(is) == 1) THEN
      ind = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz
        DO iy = miny, maxy
          y1 = (iy - nysh)*dy
          DO ix = minx, maxx
            x1 = (ix - nxsh)*dx
            ind = ind + 1
            rx = x1 - cx(is)
            ry = y1 - cy(is)
            rz = z1 - cz(is)
            rr = SQRT(rx*rx + ry*ry + rz*rz)
            rr = rr + 1D-7
            potion(ind) = potion(ind) - vgian(rr)
          END DO
        END DO
      END DO
    ELSE
      ind = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz
        DO iy = miny, maxy
          y1 = (iy - nysh)*dy
          DO ix = minx, maxx
            x1 = (ix - nxsh)*dx
            ind = ind + 1
            rx = x1 - cx(is)
            ry = y1 - cy(is)
            rz = z1 - cz(is)
            rr = SQRT(rx*rx + ry*ry + rz*rz)
            rr = rr + 1D-7
            potion(ind) = potion(ind) + v_ion_el_lgoed(rr, rloc, c1, c2, zion)
          END DO
        END DO
      END DO
    END IF
  END DO
  RETURN
END SUBROUTINE pseudogoed

! ***********************

REAL(DP) FUNCTION v_ion_el_lgoed(rr, rloc, c1, c2, zion)

! function delivering the local Goedecker pseudo-potential.
!
! Input:
! rr = distance at which potential is evaluated.
! rloc = radius PARAMETER of PsP
! c1,c2 = parameters for short range parts
! zion = ionic charge (long range Coulomb)

  USE params
  IMPLICIT NONE
  REAL(DP), INTENT(IN):: rr
  REAL(DP), INTENT(IN):: rloc
  REAL(DP), INTENT(IN):: c1
  REAL(DP), INTENT(IN):: c2
  REAL(DP), INTENT(IN):: zion

  REAL(DP) :: f1, f3
  REAL(DP), EXTERNAL :: v_soft

  IF (rr <= 7.0D0*rloc) THEN
    f1 = -zion*v_soft(rr, sq2*rloc)
    f3 = (c1 + c2*((rr/rloc)**2D0))*EXP(-0.5D0*((rr/rloc)**2))
    v_ion_el_lgoed = -e2*(f1 + f3)
  ELSE
    v_ion_el_lgoed = e2*zion/rr
  END IF

  RETURN
END FUNCTION v_ion_el_lgoed

!-----calc_proj-------------------------------------------

SUBROUTINE calc_proj(cxa, cya, cza, cxg, cyg, czg, ion)

! Computes the projectors for the non-local part of
! the Goedecker PsP for ion 'ion' at positions
! 'cxa,'cya',cza'. The positions 'cxg', 'cyg', 'czg'
! define the reference point for the center of the
! subgrid. These two vector are usually the same,
! but differ when computing forces.

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN):: cxa, cya, cza, cxg, cyg, czg
  INTEGER, INTENT(IN) :: ion

  REAL(DP) :: h0_12
  REAL(DP), PARAMETER :: fac0_12 = -0.387298334621D0 ! -0.5D0*SQRT(3D0/5D0)

!---------------------------------------------------------

  IF (h0_12g(np(ion)) .gt. -1D10) THEN
    h0_12 = h0_12g(np(ion))
  ELSE
    h0_12 = fac0_12*h0_22g(np(ion))
  END IF

  IF (ABS(h2_11g(np(ion))) + ABS(h1_22g(np(ion))) &
      + ABS(h0_33g(np(ion))) > small) THEN
    CALL calpr4(cxa, cya, cza, cxg, cyg, czg, ion)
  ELSE IF (ABS(h1_11g(np(ion))) + ABS(h0_22g(np(ion))) + ABS(h0_12) > small) THEN
    CALL calpr3(cxa, cya, cza, cxg, cyg, czg, ion)
  ELSE IF (ABS(h0_11g(np(ion))) > small) THEN
    CALL calpr2(cxa, cya, cza, cxg, cyg, czg, ion)
  END IF

  RETURN
END SUBROUTINE calc_proj

!! ****************************************

SUBROUTINE calpr2(cxact, cyact, czact, cxg, cyg, czg, ion)

! Non-local part of Goedecker PsP, case of 2 projectors.
!
! Input:
! cxact,cyact,czact = positions of ionic core
! cxg,cyg,czg = reference for center of sub-grid
! ion = actial ion NUMBER
! Other I/O handled via module 'params'.

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN):: cxact, cyact, czact, cxg, cyg, czg
  INTEGER, INTENT(IN)::ion

  INTEGER :: i, icrsx, icrsy, icrsz, ii, il, IN, ind, inn, i1, i2, i3, i1l, i2l, i3l
  REAL(DP) :: r0, r1, radion, rr, x, y, z, xion, yion, zion
  REAL(DP) :: gamfac, proj, rfac, sum1

  r0 = r0g(np(ion))
  r1 = r1g(np(ion))
  radion = radiong(np(ion))

  i1l = nint((cxg - radion)/dx)
  i2l = nint((cyg - radion)/dy)
  i3l = nint((czg - radion)/dz)
  icrsx = 2*nint(radion/dx)
  icrsy = 2*nint(radion/dy)
  icrsz = 2*nint(radion/dz)

! compute projectors on these auxiliary grid:
  ind = 0
  DO i3 = 0, icrsz
    z = (i3l + i3)*dz
    zion = z - czact
    DO i2 = 0, icrsy
      y = (i2l + i2)*dy
      yion = y - cyact
      DO i1 = 0, icrsx
        x = (i1l + i1)*dx
        xion = x - cxact
        rr = SQRT(MAX(xion*xion + yion*yion + zion*zion, small))
        IF ((rr <= radion) &
        & .AND. (((i3l + i3 + nz - 1) < nz2) & ! avoid segfault bug if going
        & .AND. ((i2l + i2 + ny - 1) < ny2) &  ! out of box
        & .AND. ((i1l + i1 + nx - 1) < nx2))) THEN

          ind = ind + 1
          IF (ind > knl) STOP " CALCPR: subgrid exceeded. enhance KNL"
          p0_1(ind, ion) = 0D0

! counter of gridpoint in one-dimensional density array:

          ii = ((i3l + i3 + nz) - 1)*nxyf + ((i2l + i2 + ny) - 1)*nx2 &
              + (i1l + i1 + nx)
          icount(ind, ion) = ii

          DO in = 1, 2
            IF (IN == 1) THEN
! projectors p0_1:
              inn = IN
              il = 0
              rfac = r0
              gamfac = 0.5D0*SQRT(pi)
            END IF
            proj = SQRT(2D0)*rr**(il + 2*(inn - 1))* &
                   EXP(-(rr*rr)/(2D0*rfac**2D0))/(rfac**(il + (4*inn - 1)/2D0)* &
                   SQRT(gamfac))

! compose the product of the radial and the angular components:
            IF (IN == 1) p0_1(ind, ion) = proj
          END DO
        END IF

      END DO
    END DO
  END DO

! end of counter-array:

  ifin(ion) = ind

! re-normalize projectors on grid

  sum1 = 0D0
  DO i = 1, ifin(ion)
    sum1 = sum1 + p0_1(i, ion)*p0_1(i, ion)
  END DO
  sum1 = sum1*dvol/(4D0*pi)
  sum1 = 1D0/SQRT(sum1)
  DO i = 1, ifin(ion)
    p0_1(i, ion) = p0_1(i, ion)*sum1
  END DO

  RETURN
END SUBROUTINE calpr2

! ****************************************

SUBROUTINE calpr3(cxact, cyact, czact, cxg, cyg, czg, ion)

! Non-local part of Goedecker PsP, case of 3 projectors.
!
! Input:
! cxact,cyact,czact = POSITION of ionic core
! cxg,cyg,czg = reference for center of sub-grid
! ion = actial ion NUMBER
! Other I/O handled via module 'params'.

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN):: cxact, cyact, czact, cxg, cyg, czg
  INTEGER, INTENT(IN)::ion

  INTEGER :: i1, i2, i3, i1l, i2l, i3l, icrsx, icrsy, icrsz, ii, il, IN, ind, inn
  REAL(DP) :: gamfac, proj, r0, r1, radion, rfac, rr, xnorm, xion, yion, zion, x, y, z
  WRITE (*, *) ' in CALPR3'

  r0 = r0g(np(ion))
  r1 = r1g(np(ion))

  radion = radiong(np(ion))

! compute boundaries of auxiliary grid for each ion:
  i1l = nint((cxg - radion)/dx)
  i2l = nint((cyg - radion)/dy)
  i3l = nint((czg - radion)/dz)
  icrsx = 2*nint(radion/dx)
  icrsy = 2*nint(radion/dy)
  icrsz = 2*nint(radion/dz)

! compute projectors on these auxiliary grid:
  ind = 0
  DO i3 = 0, icrsz
    z = (i3l + i3)*dz
    zion = z - czact
    DO i2 = 0, icrsy
      y = (i2l + i2)*dy
      yion = y - cyact
      DO i1 = 0, icrsx
        x = (i1l + i1)*dx
        xion = x - cxact
        rr = SQRT(MAX(xion*xion + yion*yion + zion*zion, small))
        IF ((rr <= radion) &
        & .AND. (((i3l + i3 + nz - 1) < nz2) & ! avoid segfault bug if going
        & .AND. ((i2l + i2 + ny - 1) < ny2) & ! out of box
        & .AND. ((i1l + i1 + nx - 1) < nx2))) THEN

          ind = ind + 1
          p0_1(ind, ion) = 0D0
          p0_2(ind, ion) = 0D0
          p1_1(ind, ion) = 0D0
          p1_1x(ind, ion) = 0D0
          p1_1y(ind, ion) = 0D0
          p1_1z(ind, ion) = 0D0

! counter of gridpoint in one-dimensional density array:

          ii = ((i3l + i3 + nz) - 1)*nxyf + ((i2l + i2 + ny) - 1)*nx2 &
              + (i1l + i1 + nx)
          icount(ind, ion) = ii

          DO IN = 1, 3
            IF (IN == 1 .OR. IN == 2) THEN
! projectors p0_1,p0_2:
              inn = IN
              il = 0
              rfac = r0
              IF (IN == 1) THEN
                gamfac = 0.5D0*SQRT(pi)
              ELSE IF (IN == 2) THEN
                gamfac = 2.5D0*1.5D0*0.5D0*SQRT(pi)
              END IF
            ELSE IF (IN == 3) THEN
! projectors p1_1,p1_2:
              il = 1
              rfac = r1
              inn = 1
              gamfac = 1.5D0*0.5D0*SQRT(pi)
            END IF
            proj = SQRT(2D0)*rr**(il + 2*(inn - 1))* &
                   EXP(-(rr*rr)/(2D0*rfac**2D0))/(rfac**(il + (4*inn - 1)/2D0)* &
                                                  SQRT(gamfac))

! compose the product of the radial and the angular components:
! p_i -> p_i*Y_lm
            IF (IN == 1) p0_1(ind, ion) = proj
            IF (IN == 2) p0_2(ind, ion) = proj
            IF (IN == 3) THEN
              p1_1(ind, ion) = proj
              p1_1x(ind, ion) = proj*xion/rr ! *SQRT(3D0)
              p1_1y(ind, ion) = proj*yion/rr ! *SQRT(3D0)
              p1_1z(ind, ion) = proj*zion/rr ! *SQRT(3D0)
            END IF
          END DO
        END IF

      END DO
    END DO
  END DO

! end of counter-array:

! normalize

  xnorm = 1D0/SQRT(dvol*SUM(p0_1(:, ion)**2)/(4D0*pi))
  p0_1(:, ion) = xnorm*p0_1(:, ion)
  WRITE (6, '(a,1pg13.5)') ' norm of 1s projector:', xnorm**(-2)

  xnorm = 1D0/SQRT(dvol*SUM(p0_2(:, ion)**2)/(4D0*pi))
  p0_2(:, ion) = xnorm*p0_2(:, ion)
  WRITE (6, '(a,1pg13.5)') ' norm of 2s projector:', xnorm**(-2)

  xnorm = 1D0/SQRT(dvol*SUM(p1_1(:, ion)**2)/(4D0*pi))
  p1_1(:, ion) = xnorm*p1_1(:, ion)
  p1_1x(:, ion) = xnorm*p1_1x(:, ion)
  p1_1y(:, ion) = xnorm*p1_1y(:, ion)
  p1_1z(:, ion) = xnorm*p1_1z(:, ion)
  WRITE (6, '(a,1pg13.5)') ' norm of 1p projector:', xnorm**(-2)

  ifin(ion) = ind

  RETURN
END SUBROUTINE calpr3

! ****************************************

SUBROUTINE calpr4(cxact, cyact, czact, cxg, cyg, czg, ion)

! Non-local part of Goedecker PsP, case of 4 projectors.
!
! Input:
! cxact,cyact,czact = POSITION of ionic core
! cxg,cyg,czg = reference for center of sub-grid
! ion = actial ion NUMBER
! Other I/O handled via module 'params'.

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN):: cxact, cyact, czact, cxg, cyg, czg
  INTEGER, INTENT(IN)::ion

  INTEGER :: i1, i2, i3, i1l, i2l, i3l, icrsx, icrsy, icrsz, ii, il, IN, ind, inn
  REAL(DP) :: gamfac, proj, r0, r1, r2, radion, rfac, rr, xion, yion, zion, x, y, z

  r0 = r0g(np(ion))
  r1 = r1g(np(ion))
  r2 = r2g(np(ion))
  radion = 3D0

! compute boundaries of auxiliary grid for each ion:
  i1l = nint((cxg - radion)/dx)
  i2l = nint((cyg - radion)/dy)
  i3l = nint((czg - radion)/dz)
  icrsx = 2*nint(radion/dx)
  icrsy = 2*nint(radion/dy)
  icrsz = 2*nint(radion/dz)

! compute projectors on these auxiliary grid:
  ind = 0
  DO i3 = 0, icrsz
    z = (i3l + i3)*dz
    zion = z - czact
    DO i2 = 0, icrsy
      y = (i2l + i2)*dy
      yion = y - cyact
      DO i1 = 0, icrsx
        x = (i1l + i1)*dx
        xion = x - cxact
        rr = SQRT(MAX(xion*xion + yion*yion + zion*zion, small))
        IF ((rr <= radion) &
        & .AND. (((i3l + i3 + nz - 1) < nz2) & ! avoid segfault bug if going
        & .AND. ((i2l + i2 + ny - 1) < ny2) & ! out of box
        & .AND. ((i1l + i1 + nx - 1) < nx2))) THEN

          ind = ind + 1
          p0_1(ind, ion) = 0D0
          p0_2(ind, ion) = 0D0
          p0_3(ind, ion) = 0D0
          p1_1(ind, ion) = 0D0
          p1_1x(ind, ion) = 0D0
          p1_1y(ind, ion) = 0D0
          p1_1z(ind, ion) = 0D0
          p1_2(ind, ion) = 0D0
          p1_2x(ind, ion) = 0D0
          p1_2y(ind, ion) = 0D0
          p1_2z(ind, ion) = 0D0
          p2_1(ind, ion) = 0D0
          p2_xy(ind, ion) = 0D0
          p2_xz(ind, ion) = 0D0
          p2_yz(ind, ion) = 0D0
          p2_xy2(ind, ion) = 0D0
          p2_z2(ind, ion) = 0D0

! counter of gridpoint in one-dimensional density array:

          ii = ((i3l + i3 + nz) - 1)*nxyf + ((i2l + i2 + ny) - 1)*nx2 + (i1l + i1 + nx)
          icount(ind, ion) = ii

          DO IN = 1, 6
            IF (IN == 1 .OR. IN == 2 .OR. IN == 3) THEN
! projectors p0_1,p0_2,p0_3:
              inn = IN
              il = 0
              rfac = r0
              IF (IN == 1) THEN
                gamfac = 0.5D0*SQRT(pi)
              ELSE IF (IN == 2) THEN
                gamfac = 2.5D0*1.5D0*0.5D0*SQRT(pi)
              ELSE IF (IN == 3) THEN
                gamfac = 4.5D0*3.5D0*2.5D0*1.5D0*0.5D0*SQRT(pi)
              END IF
            ELSE IF (IN == 4 .OR. IN == 5) THEN
! projectors p1_1,p1_2:
              il = 1
              rfac = r1
              IF (IN == 4) THEN
                inn = 1
                gamfac = 1.5D0*0.5D0*SQRT(pi)
              ELSE IF (IN == 5) THEN
                inn = 2
                gamfac = 3.5D0*2.5D0*1.5D0*0.5D0*SQRT(pi)
              END IF
            ELSE IF (IN == 6) THEN
! projector p2_1:
              il = 2
              inn = 1
              rfac = r2
              gamfac = 2.5D0*1.5D0*0.5D0*SQRT(pi)
            END IF
            proj = SQRT(2D0)*rr**(il + 2*(inn - 1))* &
                   EXP(-(rr*rr)/(2D0*rfac**2D0))/(rfac**(il + (4*inn - 1)/2D0)* &
                                                  SQRT(gamfac))

! compose the product of the radial and the angular components:
! p_i -> p_i*Y_lm
            IF (IN == 1) p0_1(ind, ion) = proj
            IF (IN == 2) p0_2(ind, ion) = proj
            IF (IN == 3) p0_3(ind, ion) = proj
            IF (IN == 4) THEN
              p1_1(ind, ion) = proj
              p1_1x(ind, ion) = proj*xion/rr
              p1_1y(ind, ion) = proj*yion/rr
              p1_1z(ind, ion) = proj*zion/rr
            END IF
            IF (IN == 5) THEN
              p1_2(ind, ion) = proj
              p1_2x(ind, ion) = proj*xion/rr
              p1_2y(ind, ion) = proj*yion/rr
              p1_2z(ind, ion) = proj*zion/rr
            END IF
            IF (IN == 6) THEN
              p2_1(ind, ion) = proj
              p2_xy(ind, ion) = proj*xion*yion/(rr*rr)
              p2_xz(ind, ion) = proj*xion*zion/(rr*rr)
              p2_yz(ind, ion) = proj*yion*zion/(rr*rr)
              p2_xy2(ind, ion) = proj*SQRT(3D0)*(xion*xion - yion*yion)/(rr*rr)
              p2_z2(ind, ion) = proj*(2D0*zion*zion &
                                      - xion*xion - yion*yion)/(rr*rr)
            END IF
          END DO
        END IF

      END DO
    END DO
  END DO

! end of counter-array:

  ifin(ion) = ind
  RETURN
END SUBROUTINE calpr4

! ********************

SUBROUTINE checkproj(ion)

! Test projectors.

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ion

  INTEGER :: i, nrowact
  REAL(DP) :: h0_12, fac0_12, h0_11, h1_11, h0_22
  REAL(DP) :: erg1, erg2, erg3, erg4, erg5, erg6
  REAL(DP) :: sum1, sum2, sum3, sum4, sum5, sum6
  REAL(DP), PARAMETER :: plimit = 0.04D0

! determine case 'nrowact'

  IF (h0_12g(np(ion)) .gt. -1D10) THEN
    h0_12 = h0_12g(np(ion))
  ELSE
    h0_12 = fac0_12*h0_22
  END IF

  IF (ABS(h2_11g(np(ion))) + ABS(h1_22g(np(ion))) &
      + ABS(h0_33g(np(ion))) < small) THEN
    nrowact = 4
  ELSE IF (ABS(h1_11g(np(ion))) + ABS(h0_22g(np(ion))) + ABS(h0_12) < small) THEN
    nrowact = 3
  ELSE IF (ABS(h0_11g(np(ion))) > small) THEN
    nrowact = 2
  ELSE
    nrowact = 1
  END IF

! checks normalization of projectors

  IF (nrowact <= 2) THEN
    sum1 = 0D0
    sum4 = 0D0
    DO i = 1, ifin(ion)
      sum1 = sum1 + p0_1(i, ion)*p0_1(i, ion)
      sum4 = sum4 + p1_1(i, ion)*p1_1(i, ion)
    END DO
    sum1 = sum1*dvol/(4D0*pi)
    sum4 = sum4*dvol/(4D0*pi)
    erg1 = ABS(sum1 - 1D0)
    erg4 = ABS(sum4 - 1D0)
    WRITE (6, '(i2,2f12.7)') ion, sum1, sum4
    IF (erg1 >= plimit .OR. erg4 >= plimit) THEN
      WRITE (*, *) ' ion,sum1,sum4,r=', ion, sum1, sum4, cx(ion), cy(ion), cz(ion)
      STOP 'projectors are not normalized'
    END IF
  ELSE IF (nrowact == 3) THEN
    sum1 = 0D0
    sum2 = 0D0
    sum4 = 0D0
    DO i = 1, ifin(ion)
      sum1 = sum1 + p0_1(i, ion)*p0_1(i, ion)
      sum2 = sum2 + p0_2(i, ion)*p0_2(i, ion)
      sum4 = sum4 + p1_1(i, ion)*p1_1(i, ion)
    END DO
    sum1 = sum1*dvol/(4D0*pi)
    sum2 = sum2*dvol/(4D0*pi)
    sum4 = sum4*dvol/(4D0*pi)
    erg1 = ABS(sum1 - 1D0)
    erg2 = ABS(sum2 - 1D0)
    erg4 = ABS(sum4 - 1D0)
    IF (ABS(h0_11) < 1D-5) erg1 = 0D0
    IF (ABS(h0_22) < 1D-5) erg2 = 0D0
    IF (ABS(h1_11) < 1D-5) erg4 = 0D0
    WRITE (6, '(i2,3f12.7)') ion, sum1, sum2, sum4
    IF (erg1 >= plimit .OR. erg2 >= plimit .OR. erg4 >= plimit) &
      STOP 'projectors are not normalized'
  ELSE IF (nrowact == 4) THEN
    sum1 = 0D0
    sum2 = 0D0
    sum3 = 0D0
    sum4 = 0D0
    sum5 = 0D0
    sum6 = 0D0
    DO i = 1, ifin(ion)
      sum1 = sum1 + p0_1(i, ion)*p0_1(i, ion)
      sum2 = sum2 + p0_2(i, ion)*p0_2(i, ion)
      sum3 = sum3 + p0_3(i, ion)*p0_3(i, ion)
      sum4 = sum4 + p1_1(i, ion)*p1_1(i, ion)
      sum5 = sum5 + p1_2(i, ion)*p1_2(i, ion)
      sum6 = sum6 + p2_1(i, ion)*p2_1(i, ion)
    END DO
    sum1 = sum1*dvol/(4D0*pi)
    sum2 = sum2*dvol/(4D0*pi)
    sum3 = sum3*dvol/(4D0*pi)
    sum4 = sum4*dvol/(4D0*pi)
    sum5 = sum5*dvol/(4D0*pi)
    sum6 = sum6*dvol/(4D0*pi)
    erg1 = ABS(sum1 - 1D0)
    erg2 = ABS(sum2 - 1D0)
    erg3 = ABS(sum3 - 1D0)
    erg4 = ABS(sum4 - 1D0)
    erg5 = ABS(sum5 - 1D0)
    erg6 = ABS(sum6 - 1D0)
    WRITE (6, '(i2,6f12.7)') ion, sum1, sum2, sum3, sum4, sum5, sum6
    IF (erg1 >= plimit .OR. erg2 >= plimit &
        .OR. erg3 >= plimit .OR. erg4 >= plimit &
        .OR. erg5 >= plimit .OR. erg6 >= plimit) &
      STOP 'projectors are not normalized'
  END IF
  RETURN
END SUBROUTINE checkproj

