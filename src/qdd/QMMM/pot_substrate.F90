
#if(raregas)
!-----pseudosoft_substrate-----------------------------------------

SUBROUTINE pseudosoft_substrate(pseudorho, potsave)

! The routine computes the potentials from substrate and
! possibly its images IN a dielectric layer.
! Input and output are:
! pseudorho = accumulator for pseudo-density
! potsave = accumulator for pseudo-density

!--------------------------------------------------------------

  USE params
  IMPLICIT NONE
! DIMENSION rho(2*kdfull2)

  REAL(DP), INTENT(IN OUT) :: pseudorho(kdfull2)
  REAL(DP), INTENT(IN OUT) :: potsave(kdfull2)

! DIMENSION potshort(kdfull2)
! DIMENSION ri(3)
  INTEGER :: i
  INTEGER, EXTERNAL :: isoutofbox
  REAL(DP), EXTERNAL :: v_soft

!------------------------------------------------------------------

  IF (isurf == 0) RETURN
  DO i = 1, nc + NE + nk
    CALL getparas(i)

    IF (ioutofbox(i) == 0 .AND. imobtmp == 1) THEN

      IF (idielec == 0) THEN

        CALL addgsmdensity(pseudorho, rvectmp(1), rvectmp(2), &
                           rvectmp(3), sigtmp, chgtmp, 0)

      ELSE

        IF (rvectmp(1) < xdielec) THEN

          CALL addgsmdensity(pseudorho, rvectmp(1), rvectmp(2), &
                             rvectmp(3), sigtmp, 2D0/(epsdi + 1D0)*chgtmp, 0)

        ELSE
          CALL addgsmdensity(pseudorho, rvectmp(1), rvectmp(2), &
                             rvectmp(3), sigtmp, chgtmp, 0)

          IF (isoutofbox(-rvectmp(1) + 2D0*xdielec, &
                         rvectmp(2), rvectmp(3)) == 0) THEN
            CALL addgsmdensity(pseudorho, &
                               -rvectmp(1) + 2D0*xdielec, rvectmp(2), rvectmp(3), &
                               sigtmp, -chgtmp*(epsdi - 1)/(epsdi + 1), 0)
          ELSE ! image charge is outside
            CALL addfunctofield1(potsave, v_soft, -rvectmp(1) + 2D0*xdielec, &
                                 rvectmp(2), rvectmp(3), -e2*chgtmp*(epsdi - 1)/(epsdi + 1), sigtmp*sq2)
          END IF ! isOutOfBox
        END IF ! rVecTmp(1)

      END IF ! iDielec
    ELSE ! outside box or inside box and imobile

      IF (imobtmp == 1) THEN ! particle is outside but mobile

        IF (idielec == 0) THEN

          CALL addfunctofield1(potsave, v_soft, rvectmp(1), &
                               rvectmp(2), rvectmp(3), e2*chgtmp, sigtmp*sq2)
        ELSE

          IF (rvectmp(1) < xdielec) THEN
            CALL addfunctofield1(potsave, v_soft, &
                                 rvectmp(1), rvectmp(2), rvectmp(3), &
                                 2D0/(epsdi + 1D0)*e2*chgtmp, sigtmp*sq2)
          ELSE
            CALL addfunctofield1(potsave, v_soft, &
                                 rvectmp(1), rvectmp(2), rvectmp(3), e2*chgtmp, sigtmp*sq2)
            CALL addfunctofield1(potsave, v_soft, -rvectmp(1) + 2D0*xdielec, &
                                 rvectmp(2), rvectmp(3), -e2*chgtmp*(epsdi - 1)/(epsdi + 1), sigtmp*sq2)
          END IF

        END IF ! iDielec

      END IF ! imobTmp
    END IF ! iOutOfBox
  END DO

  RETURN
END SUBROUTINE pseudosoft_substrate

!************************************************************
SUBROUTINE getvdwpot
!************************************************************
  USE params
  IMPLICIT NONE

  INTEGER :: ind, is, ix, iy, iz
  REAL(DP) :: rr, x1, y1, z1
  REAL(DP) :: ri(3)

  REAL(DP), EXTERNAL :: v_vdw

  DO ind = 1, nxyz
    potvdw(ind) = 0D0
  END DO

  DO is = 1, nc
!test WRITE(6,'(a,2i4)') ' IN CALCPSEUDO: IS,NP=',is,np(is)
    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      ri(3) = zc(is) - z1
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        ri(2) = yc(is) - y1
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ri(1) = xc(is) - x1

          rr = SQRT(ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3))
          rr = MAX(rr, small) ! avoid zero
          ind = ind + 1

          potvdw(ind) = potvdw(ind) + v_vdw(ri, rr, is, 1D0)

        END DO
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE getvdwpot

!************************************************************

SUBROUTINE addgsmpot(field, iswitch)
!************************************************************
  USE params
  IMPLICIT NONE
! add electrostatic potentials from GSM particles TO potion(kdfull2)
! iswitch = 0 --> calculate potential of fixed ions
! 1 --> calculate potential of mobile ions

  REAL(DP), INTENT(OUT) :: field(kdfull2)
  INTEGER, INTENT(IN) :: iswitch

  INTEGER :: ind, is, ix, iy, iz
  REAL(DP) :: rr, rx, ry, rz, sigc, sigk, sigv, x1, y1, z1
  REAL(DP), EXTERNAL :: v_soft
! IF (nelect.eq.0) RETURN

  IF (idielec == 0) THEN

    DO is = 1, nc

      IF (imobc(is) == iswitch) THEN ! IF .eq. 0 THEN treat is as static potential

        sigc = sigmac

        chgtmp = chgc(is)*e2

        ind = 0
        DO iz = minz, maxz
          z1 = (iz - nzsh)*dz
          rz = z1 - zc(is)
          DO iy = miny, maxy
            y1 = (iy - nysh)*dy
            ry = y1 - yc(is)
            DO ix = minx, maxx
              x1 = (ix - nxsh)*dx
              rx = x1 - xc(is)
              rr = SQRT(rx*rx + ry*ry + rz*rz)
              rr = MAX(rr, small) ! avoid zero
              ind = ind + 1

              field(ind) = field(ind) + chgtmp*v_soft(rr, sigc*sq2)
! replaced by tabulated version of V_soft: !!! TOO SLOW!!!
! field(ind)=field(ind) +
! & chgtmp*e2*V_soft_tabc(rr)

            END DO
          END DO
        END DO

      END IF

    END DO

    DO is = 1, NE

      IF (imobe(is) == iswitch) THEN

        sigv = sigmav
        chgtmp = chge(is)*e2

        ind = 0
        DO iz = minz, maxz
          z1 = (iz - nzsh)*dz
          rz = z1 - ze(is)
          DO iy = miny, maxy
            y1 = (iy - nysh)*dy
            ry = y1 - ye(is)
            DO ix = minx, maxx
              x1 = (ix - nxsh)*dx
              rx = x1 - xe(is)
              rr = SQRT(rx*rx + ry*ry + rz*rz)
              rr = MAX(rr, small) ! avoid zero
              ind = ind + 1

              field(ind) = field(ind) + chgtmp*v_soft(rr, sigv*sq2)
! field(ind)=field(ind) +
! & chgtmp*e2*V_soft_tabe(rr)

            END DO
          END DO
        END DO

      END IF

    END DO

    DO is = 1, nk

      IF (imobk(is) == iswitch) THEN

        sigk = sigmak
        chgtmp = chgk(is)*e2

        ind = 0
        DO iz = minz, maxz
          z1 = (iz - nzsh)*dz
          rz = z1 - zk(is)
          DO iy = miny, maxy
            y1 = (iy - nysh)*dy
            ry = y1 - yk(is)
            DO ix = minx, maxx
              x1 = (ix - nxsh)*dx
              rx = x1 - xk(is)
              rr = SQRT(rx*rx + ry*ry + rz*rz)
              rr = MAX(rr, small) ! avoid zero
              ind = ind + 1

              field(ind) = field(ind) + chgtmp*v_soft(rr, sigk*sq2)
! field(ind)=field(ind) +
! & chgtmp*e2*V_soft_tabk(rr)

            END DO
          END DO
        END DO

      END IF

    END DO

  END IF ! iDielec.eq.0

  IF (idielec /= 0) THEN ! image potential

    DO is = 1, nc

      IF (imobc(is) == iswitch) THEN ! IF .eq. 0 THEN treat is as static potential
        sigc = sigmac

        chgtmp = chgc(is)

        IF (xc(is) > xdielec) THEN

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - zc(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - yc(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 - xc(is)
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) + chgtmp*e2*v_soft(rr, sigc*sq2)

              END DO
            END DO
          END DO

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - zc(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - yc(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 + xc(is) - 2D0*xdielec
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) - (epsdi - 1D0)/(epsdi + 1D0)* &
                             chgtmp*e2*v_soft(rr, sigc*sq2)

              END DO
            END DO
          END DO

        ELSE

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - zc(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - yc(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 - xc(is)
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) + 2D0/(epsdi + 1D0)* &
                             chgtmp*e2*v_soft(rr, sigc*sq2)

              END DO
            END DO
          END DO

        END IF
      END IF

    END DO

    DO is = 1, NE

      IF (imobe(is) == iswitch) THEN

        sigv = sigmav
        chgtmp = chge(is)

        IF (xe(is) > xdielec) THEN

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - ze(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - ye(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 - xe(is)
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) + chgtmp*e2*v_soft(rr, sigv*sq2)

              END DO
            END DO
          END DO

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - ze(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - ye(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 + xe(is) - 2D0*xdielec
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) - (epsdi - 1D0)/(epsdi + 1D0)* &
                             chgtmp*e2*v_soft(rr, sigv*sq2)

              END DO
            END DO
          END DO

        ELSE ! xe(is).lt.xDielec

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - ze(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - ye(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 - xe(is)
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) + 2D0/(epsdi + 1D0)* &
                             chgtmp*e2*v_soft(rr, sigv*sq2)

              END DO
            END DO
          END DO

        END IF ! xe(is).lt.xDielec
      END IF ! imobe(is).eq.iswitch

    END DO

    DO is = 1, nk

      IF (imobk(is) == iswitch) THEN

        sigk = sigmak
        chgtmp = chgk(is)

        IF (xk(is) > xdielec) THEN

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - zk(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - yk(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 - xk(is)
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) + chgtmp*e2*v_soft(rr, sigk*sq2)

              END DO
            END DO
          END DO

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - zk(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - yk(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 + xk(is) - 2D0*xdielec
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) - (epsdi - 1D0)/(epsdi + 1D0)* &
                             chgtmp*e2*v_soft(rr, sigk*sq2)

              END DO
            END DO
          END DO

        ELSE

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            rz = z1 - zk(is)
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              ry = y1 - yk(is)
              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                rx = x1 - xk(is)
                rr = SQRT(rx*rx + ry*ry + rz*rz)
                rr = MAX(rr, small) ! avoid zero
                ind = ind + 1

                field(ind) = field(ind) + 2D0/(epsdi + 1D0)* &
                             chgtmp*e2*v_soft(rr, sigk*sq2)

              END DO
            END DO
          END DO

        END IF ! xk(is).gt.xDielec
      END IF ! imobk(is).eq.iswitch

    END DO

  END IF

  RETURN
END SUBROUTINE addgsmpot

!************************************************************

SUBROUTINE addgsmdensity(field, xx, yy, zz, sigm, ccharge, iparit)
!************************************************************
! adds (iparit=even) or subtracts (iparit=odd) Gaussian
! density TO field on subgrid centered at xx,yy,zz

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz
  REAL(DP), INTENT(IN) :: sigm
  REAL(DP), INTENT(IN) :: ccharge
  INTEGER, INTENT(IN) :: iparit

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1
  INTEGER, EXTERNAL :: conv3to1, getnearestgridpoint

  ind = getnearestgridpoint(xx, yy, zz)

  CALL conv1to3(ind)

  DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
    z1 = (iz - nzsh)*dz
    rz = z1 - zz
    DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
      y1 = (iy - nysh)*dy
      ry = y1 - yy
      DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg
        x1 = (ix - nxsh)*dx
        rx = x1 - xx
        rr = rx*rx + ry*ry + rz*rz

        ind = conv3to1(ix, iy, iz)

        field(ind) = field(ind) &
                     + (-1D0)**iparit*ccharge*EXP(-rr/(2D0*sigm**2)) &
                     /(pi**1.5D0*2D0**1.5D0*sigm**3)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE addgsmdensity
!************************************************************

!************************************************************

SUBROUTINE addgsmdensities(field)
!************************************************************
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)

  INTEGER :: ind, is, ix, iy, iz
  REAL(DP) :: rr, rx, ry, rz, sigc, sigv, sigk, x1, y1, z1
  INTEGER, EXTERNAL :: conv3to1

  DO is = 1, nc

    IF (imobc(is) /= 0) THEN

      ind = isubgcenter(is)
      sigc = sigmac ! TO be replaced by array IN the future
      chgtmp = chgc(is)

      CALL conv1to3(ind)

      DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
        z1 = (iz - nzsh)*dz
        rz = z1 - zc(is)
        DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
          y1 = (iy - nysh)*dy
          ry = y1 - yc(is)
          DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg
            x1 = (ix - nxsh)*dx
            rx = x1 - xc(is)
            rr = rx*rx + ry*ry + rz*rz

            ind = conv3to1(ix, iy, iz)

            field(ind) = field(ind) + chgtmp*EXP(-rr/(2D0*sigc**2)) &
                         /(pi**1.5D0*2D0**1.5D0*sigc**3)
          END DO
        END DO
      END DO

    END IF

  END DO

  DO is = 1, NE

    IF (imobe(is) /= 0) THEN

      ind = isubgcenter(is + nc)
      sigv = sigmav ! TO be replaced by array IN the future
      chgtmp = chge(is)

      CALL conv1to3(ind)

      DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
        z1 = (iz - nzsh)*dz
        rz = z1 - ze(is)
        DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
          y1 = (iy - nysh)*dy
          ry = y1 - ye(is)
          DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg
            x1 = (ix - nxsh)*dx
            rx = x1 - xe(is)
            rr = rx*rx + ry*ry + rz*rz

            ind = conv3to1(ix, iy, iz)

            field(ind) = field(ind) + chgtmp*EXP(-rr/(2D0*sigv**2)) &
                         /(pi**1.5D0*2D0**1.5D0*sigv**3)
          END DO
        END DO
      END DO

    END IF

  END DO

  DO is = 1, nk

    IF (imobk(is) /= 0) THEN

      ind = isubgcenter(is + nc + NE)
      sigk = sigmak ! TO be replaced by array IN the future
      chgtmp = chgk(is)

      CALL conv1to3(ind)

      DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
        z1 = (iz - nzsh)*dz
        rz = z1 - zk(is)
        DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
          y1 = (iy - nysh)*dy
          ry = y1 - yk(is)
          DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg
            x1 = (ix - nxsh)*dx
            rx = x1 - xk(is)
            rr = rx*rx + ry*ry + rz*rz

            ind = conv3to1(ix, iy, iz)

            field(ind) = field(ind) + chgtmp*EXP(-rr/(2D0*sigk**2)) &
                         /(pi**1.5D0*2D0**1.5D0*sigk**3)
          END DO
        END DO
      END DO

    END IF

  END DO

  RETURN
END SUBROUTINE addgsmdensities

! ****************************

SUBROUTINE calc_frho(rho)

! ****************************

! frho(i,is) is the component i of the integral of (grad.Vsoft)/r over rho for
! the raregas atom is
! for computation of Van der Waals potential and forces

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)

  INTEGER :: i, ico, ind, is, ix, iy, iz
  REAL(DP) :: r, r2, rhodvsdr, sgm, x1, y1, z1
  REAL(DP) :: ri(3)

  REAL(DP), EXTERNAL:: dvsdr

  DO is = 1, nc
    DO i = 1, 3
      frho(is, i) = 0D0
    END DO
  END DO

  sgm = sigmac*sq2

  DO is = 1, nc
    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      ri(3) = zc(is) - z1
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        ri(2) = yc(is) - y1
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ri(1) = xc(is) - x1
          ind = ind + 1

          r2 = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
          r = SQRT(r2)
          r = MAX(r, small)

          rhodvsdr = rho(ind)*dvsdr(r, sgm)/r

          DO ico = 1, 3
            frho(is, ico) = frho(is, ico) + rhodvsdr*ri(ico)
          END DO
        END DO
      END DO
    END DO
  END DO
  DO is = 1, nc
    DO i = 1, 3
      frho(is, i) = frho(is, i)*dvol
    END DO
  END DO

  RETURN
END SUBROUTINE calc_frho

! ************************************

REAL(DP) FUNCTION v_vdw(ri, r, is, fac)

! ************************************

! computation of the Van der Waals potential

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: ri(3)
  REAL(DP), INTENT(IN) :: r
  INTEGER, INTENT(IN) :: is
  REAL(DP), INTENT(IN) :: fac

! frho is a vector of configuration space
  INTEGER :: ico
  REAL(DP) :: dvsdrtmp, elnum, sgm
  REAL(DP), PARAMETER :: alpha_ar = 10.6D0

  REAL(DP), EXTERNAL :: dvsdr

  elnum = nelect*1D0

  sgm = sigmac*sq2

  dvsdrtmp = dvsdr(r, sgm)

  v_vdw = -0.5D0*dvsdrtmp**2

  DO ico = 1, 3
    v_vdw = v_vdw + frho(is, ico)*dvsdrtmp*ri(ico)/(elnum*r)*fac
  END DO

  v_vdw = e2*alpha_ar*v_vdw

  RETURN
END FUNCTION v_vdw

!----V_Ar_el_core--------------------------------------------------------------

REAL(DP) FUNCTION v_ar_el_core(r)
! corepotential of Ar (sfort-range part)
!---------------------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), PARAMETER :: coreheight = 0.47D0
  REAL(DP), PARAMETER :: coreradius = 2.2D0
  REAL(DP), PARAMETER :: corewidth = 1.6941D0

  v_ar_el_core = -e2*coreheight*(1.0D0 + EXP(corewidth*(r - coreradius)))**(-1)

  RETURN
END FUNCTION v_ar_el_core

!-----V_Ar_el----------------------------------------------------------

REAL(DP) FUNCTION v_ar_el(r)
  USE params
  IMPLICIT NONE
! Ar electron potential
! from Xiao-Ping Li et al, J.Chem.Phy. 85 (1986) 3444,
! fitted dipole response subtracted.

  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: r2, el_ar_psp, v_dip
  REAL(DP) :: coef(0:14)
  ! the coefficients
  DATA coef(0)/7.26018748D-1/
  DATA coef(1)/-1.89605233D-2/
  DATA coef(2)/3.55049697D-4/
  DATA coef(3)/-1.18202645D-1/
  DATA coef(4)/-9.27313880D-2/
  DATA coef(5)/2.37970670D-1/
  DATA coef(6)/-1.15588294D-1/
  DATA coef(7)/-4.83199858D-2/
  DATA coef(8)/7.94876042D-2/
  DATA coef(9)/-4.07365689D-2/
  DATA coef(10)/1.16000480D-2/
  DATA coef(11)/-2.01816503D-3/
  DATA coef(12)/2.14061474D-4/
  DATA coef(13)/-1.27688158D-5/
  DATA coef(14)/3.29414225D-7/

!?? DATA Vdip/11.08/ ! Ar polarizability Ha
  REAL(DP), PARAMETER :: alpha_ar = 10.6D0 ! polarizability of Ar

! parameters of the fitted dipole response

  REAL(DP) :: rfit(0:2)
! DATA rfit/7.629608,15.47005,0.2018586/ ! fitted for a=11.08
  DATA rfit/8.161628D0, 15.60000D0, 0.4996526D0/ ! fitted for a=10.6

  REAL(DP), EXTERNAL :: v_soft
!-------------------------------------------------------------------

  r2 = r*r

! the full PsP

  IF (r < 4.275) THEN
    el_ar_psp = e2*(coef(0) + r*(coef(1) + &
                                 r*(coef(2) + r*(coef(3) + &
                                                 r*(coef(4) + r*(coef(5) + &
                                                                 r*(coef(6) + r*(coef(7) + &
                                                                                 r*(coef(8) + r*(coef(9) + &
                                                                                                 r*(coef(10) + r*(coef(11) + &
                                                                                                       r*(coef(12) + r*(coef(13) + &
                                                                                                            r*coef(14)))))))))))))))
  ELSE
    el_ar_psp = -e2*alpha_ar*0.5D0/(r2*r2)
  END IF

! the softened dipole part subtracted

  v_dip = -e2*0.5D0*alpha_ar*r2/ &
          (rfit(0) + r2*(rfit(1) + r2*(rfit(2) + r2)))

  v_ar_el = el_ar_psp - v_dip - &
            ch(18)*e2*v_soft(r, sgm1(18)*sq2)

! WRITE(6,'(1x,f8.2,3g13.5)') r,V_dip,el_ar_psp,V_Ar_el

  RETURN
END FUNCTION v_ar_el

!-----V_Ar_el_dyn-------------------------------------------------------

REAL(DP) FUNCTION v_ar_el_dyn(r)

  USE params
  IMPLICIT NONE

! computes the reduced effective Ar-electron potential.
! IN deviation from the notes, the 'core' potential is
! fully associated with the Ar and nothing TO the valence cloud.

  REAL(DP), INTENT(IN) :: r

  INTEGER :: iter
  REAL(DP) :: effc, effch
  REAL(DP) :: rstep, vdip, vrep, rcor, epsil, rsmall
  REAL(DP) :: rpmin, rpold, sigma_ar, v_core, v_dipa
!DATA c_dipmod/3.378/ ! spring constant of dipole model IN Ha
! DATA sigma_Ar/1.43/
  DATA rstep/1.0D-3/
  DATA vdip/11.08D0/ ! effective polar.pot. IN Ha
  DATA vrep/52.8D0/ ! repulsive core IN Ha
! DATA rcut/3.1D0/ ! cut radius for pol.pot.
  DATA rcor/1.35D0/ ! inverse radius**2 for core
! DIMENSION ch(18:18)
! DATA ch(18)/6.119/
! DATA e2/2.0/
  DATA epsil/1.0D-6/
  DATA rsmall/0.012D0/

! switch TO simple core potential

  LOGICAL :: tcore
  DATA tcore/.true./

  REAL(DP), EXTERNAL :: v_soft

!-------------------------------------------------------------------

  IF (r < rsmall) THEN
    v_ar_el_dyn = vrep*e2
    RETURN
  END IF
  sigma_ar = sgm1(18)*sq2
  effch = ch(18)*e2
  effc = c_dipmod*e2 ! IN Ry
  IF (tcore) THEN
    v_core = vrep*e2*EXP(-rcor*r*r)
  ELSE
    rpmin = -effch*(v_soft(r + rstep, sigma_ar) - v_soft(r - rstep, sigma_ar)) &
            /(2D0*rstep*effc)
! WRITE(6,'(1x,i4,2g14.5)') 0,rpmin,abs(rpold-rpmin)
    DO iter = 1, 100
      rpold = rpmin
      rpmin = -effch*(v_soft(r + rstep + rpmin, sigma_ar) &
                      - v_soft(r + rpmin - rstep, sigma_ar))/(2D0*rstep*effc)
! WRITE(6,'(1x,i4,2g14.5)') iter,rpmin,abs(rpold-rpmin)
      IF (ABS(rpold - rpmin) < epsil) EXIT
    END DO
    v_dipa = effch*(v_soft(r + rpmin, sigma_ar) - v_soft(r, sigma_ar)) &
             + 0.5D0*effc*rpmin*rpmin
    v_core = vrep*e2*EXP(-rcor*r*r) - vdip*e2/(1D0 + EXP((3.1D0/r)**8))/r**4 &
             - v_dipa
  END IF

  v_ar_el_dyn = v_core - effch*v_soft(r, sigma_ar)

!test WRITE(6,'(1x,f8.2,3g13.5)') r,V_dipa,V_Ar_el,V_Ar_el+V_dipa

  RETURN
END FUNCTION v_ar_el_dyn

!-----V_ion_el------------------------------------------------------------

REAL(DP) FUNCTION v_ion_el(rr, nptype)

  USE params
  IMPLICIT NONE

! effective ion-electron potential
! rr = distance at which potential is computed
! nptype = TYPE of ion (see array np(*) IN 'init.F')

  REAL(DP), INTENT(IN) :: rr
  INTEGER, INTENT(IN) :: nptype

  REAL(DP), EXTERNAL :: v_ar_el_core, v_soft
!------------------------------------------------------------------------

  IF (nptype == -18) THEN
    v_ion_el = chg1(-18)*e2*v_soft(rr, sgm1(-18)*sq2)
  ELSE IF (nptype < 18) THEN
    v_ion_el = chg1(nptype)*e2*v_soft(rr, sgm1(nptype)*sq2) &
               + chg2(nptype)*e2*v_soft(rr, sgm2(nptype)*sq2)
  ELSE
    v_ion_el = v_ar_el_core(rr) + chg1(18)*e2*v_soft(rr, sgm1(18)*sq2)
  END IF

  RETURN
END FUNCTION v_ion_el

!-----V_Ar_Ar-----------------------------------------------------

REAL(DP) FUNCTION v_ar_ar(r)

  USE params
  IMPLICIT NONE

! Ar-Ar potential

  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: alpha_ar, beta_ar, c6_ar, c8_ar, a_ar
  REAL(DP) :: rabs, core
  DATA alpha_ar, beta_ar, c6_ar, c8_ar, a_ar &
    /1.7301D0, 1.7966D0, 55.465D0, 3672.9D0, 794.21D0/

!----------------------------------------------------------------------

  rabs = MAX(ABS(r), small)
  core = a_ar*EXP(-alpha_ar*rabs)/rabs - 2D0/(1D0 + EXP((beta_ar/rabs)**2)) &
         *(c6_ar/rabs**6 + c8_ar/rabs**8)

  v_ar_ar = e2*(core)

! V_Ar_Ar = e2*(core + ch(18)**2*V_soft(r,2D0*sgm1(18)))

  RETURN
END FUNCTION v_ar_ar

!-----V_Ar_Na-----------------------------------------------------

REAL(DP) FUNCTION v_ar_na(r)

  USE params
  IMPLICIT NONE

! Ar-Na potential

! parameters for repulsive core (energies IN Ha, rescale TO Ry by e2 !)

! DATA vrep/214.14/ ! strength of repulsive core
! DATA rcor/1.6906/ ! inverse radius of repulsive core
! DATA vdip2/308.98/ ! 2*strength of higher dipole part
! DATA rdip2/2.20/ ! radius of cutoff IN higher dipole part

!----------------------------------------------------------------------

! rabs = abs(r)

! V_Ar_Na =
! & e2*(vrep*exp(-rcor*rabs)/rabs
! & - vdip2/((1.0d0+exp(rdip2/rabs)**2)*rabs**6)
! & + ch(18)*ch(11)*V_soft(r,sgm1(18))
! & )

  REAL(DP), INTENT(IN) :: r
  REAL(DP) a1, a2, a3, b1, b2
  DATA a1/334.85D0/
  DATA a2/52.5D0/
  DATA a3/1383.0D0/
  DATA b1/-1.7624D0/
  DATA b2/1.815D0/

  v_ar_na = a1*EXP(b1*r)/r - 2D0/(1D0 + EXP(b2/r))*(a2/r**6 + a3/r**8)
  v_ar_na = e2*v_ar_na

  RETURN
END FUNCTION v_ar_na
#endif
