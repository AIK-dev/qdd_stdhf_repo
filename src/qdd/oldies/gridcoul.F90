
!-----solv_FFT for finite difference package-----------------------------------
MODULE coulsolv
  USE params
  IMPLICIT REAL(DP) (A - H, O - Z)

  USE kinetic

CONTAINS
  SUBROUTINE solv_fft(rhoin, pot, deltax, deltay, deltaz)

    USE params

! For given density 'rhoin'
! solve Poisson equation with b.c. 0 at infinity
! using sine-FFT.

    INCLUDE 'work3D.inc'

    REAL, INTENT(IN OUT) :: rhoin(-maxx:maxx, -maxy:maxy, -
    REAL, INTENT(OUT) :: pot(-maxx:maxx, -maxy:maxy, -
    REAL, INTENT(IN OUT) :: deltax
    REAL, INTENT(IN OUT) :: deltay
    REAL, INTENT(IN OUT) :: deltaz

    REAL :: potbcxl(-maxy:maxy, -maxz:maxz)
    REAL :: potbcxu(-maxy:maxy, -maxz:maxz)
    REAL :: potbcyl(-maxx:maxx, -maxz:maxz)
    REAL :: potbcyu(-maxx:maxx, -maxz:maxz)
    REAL :: potbczl(-maxx:maxx, -maxy:maxy)
    REAL :: potbczu(-maxx:maxx, -maxy:maxy)

    REAL :: fac
    INTEGER :: ix, iy, iz

!-----------------------------------------------------------------------

! inhomogeneous equation with b.c. zero

    CALL d3sinfinverse(rhoin, pot, zero)
    fac = 4.0D0*pi
    DO iz = -maxmz, maxmz
      DO iy = -maxmy, maxmy
        DO ix = -maxmx, maxmx
          pot(ix, iy, iz) = fac*pot(ix, iy, iz)
        END DO
      END DO
    END DO

! Coulomb field at bounds, by multipole expansion

    CALL boundc(rhoin, deltax, deltay, deltaz, &
                potbcxl, potbcxu, potbcyl, potbcyu, potbczl, potbczu, pot, 0)

! add homogeneous solution from all six boundaries of the cube

    CALL homog_sol(pot, deltax, deltay, deltaz, &
                   potbcxl, potbcxu, potbcyl, potbcyu, potbczl, potbczu)

! scale by charge

    DO iz = -maxmz, maxmz
      DO iy = -maxmy, maxmy
        DO ix = -maxmx, maxmx
          pot(ix, iy, iz) = e2*pot(ix, iy, iz)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE solv_fft

!-----homog_sol-----------------------------------------------------------

  SUBROUTINE homog_sol(pot, deltax, deltay, deltaz, &
                       potbcxl, potbcxu, potbcyl, potbcyu, potbczl, potbczu)

    USE params

! Compute homogeneous solution for given potentials
! 'potbcxl',...,'potbczu' at the six boundary planes of the cube.
! Inhomogeneous solution is entered on 'pot',
! homogeneous solution is added, and RESULT is returned on 'pot'.
! The mesh spacing is entered through 'deltax,deltay,deltaz'
! and the mesh SIZE is prescribed IN 'all.inc'.

! Note: This code is optimized for readability not yet for speed.

    INCLUDE 'work3D.inc'

    REAL, INTENT(OUT) :: pot(-maxx:maxx, -maxy:maxy, -
    REAL, INTENT(IN) :: deltax
    REAL, INTENT(IN) :: deltay
    REAL, INTENT(IN) :: deltaz
    REAL, INTENT(IN OUT) :: potbcxl(-maxy:maxy, -maxz:maxz)
    REAL, INTENT(IN OUT) :: potbcxu(-maxy:maxy, -maxz:maxz)
    REAL, INTENT(IN OUT) :: potbcyl(-maxx:maxx, -maxz:maxz)
    REAL, INTENT(IN OUT) :: potbcyu(-maxx:maxx, -maxz:maxz)
    REAL, INTENT(IN OUT) :: potbczl(-maxx:maxx, -maxy:maxy)
    REAL, INTENT(IN OUT) :: potbczu(-maxx:maxx, -maxy:maxy)

! local variables

    LOGICAL :: testpr
    INTEGER :: ix, iy, iz
    REAL :: fac, xl, xu, yl, yu, zl, zu

    REAL :: gammax(-maxy:maxy, -maxz:maxz)
    REAL :: sinhlx(-maxy:maxy, -maxz:maxz)
    REAL :: gammay(-maxx:maxx, -maxz:maxz)
    REAL :: sinhly(-maxx:maxx, -maxz:maxz)
    REAL :: gammaz(-maxx:maxx, -maxy:maxy)
    REAL :: sinhlz(-maxx:maxx, -maxy:maxy)
    REAL :: potbcx(-maxy:maxy, -maxz:maxz)
    REAL :: potbcy(-maxx:maxx, -maxz:maxz)
    REAL :: potbcz(-maxx:maxx, -maxy:maxy)

    DATA testpr/.false./

! stabilised sinh as internal statement FUNCTION

    sins(x) = SINH(SIGN(MIN(ABS(x), 72.0), x))

!-----------------------------------------------------------------------

!test WRITE(6,'(a)') 'POTBCXU:'
!test WRITE(6,'(4g12.4)') potbcxu
!test WRITE(6,'(a)') 'POTBCXL:'
!test WRITE(6,'(4g12.4)') potbcxl
!test WRITE(6,'(a)') 'POTBCYU:'
!test WRITE(6,'(4g12.4)') potbcyu
!test WRITE(6,'(a)') 'POTBCYL:'
!test WRITE(6,'(4g12.4)') potbcyl
!test WRITE(6,'(a)') 'POTBCZU:'
!test WRITE(6,'(4g12.4)') potbczu
!test WRITE(6,'(a)') 'POTBCZL:'
!test WRITE(6,'(4g12.4)') potbczl
! WRITE(6,'(a)') 'KINX:'
! WRITE(6,'(4g12.4)') kinx
! WRITE(6,'(a)') 'KINY:'
! WRITE(6,'(4g12.4)') kiny
! WRITE(6,'(a)') 'KINZ:'
! WRITE(6,'(4g12.4)') kinz

    IF (testpr) WRITE (6, '(a/6g12.4)') ' bounds:', &
      potbczl(0, 0), potbczu(0, 0), potbcyl(0, 0), potbcyu(0, 0), potbcxl(0, 0), potbcxu(0, 0)

! contribution from z-boundaries

    CALL d2xysinft(potbczu, potbczu)
    CALL d2xysinft(potbczl, potbczl)
    DO iy = -maxmy, maxmy
      DO ix = -maxmx, maxmx
        gammaz(ix, iy) = SQRT(kinx(ix) + kiny(iy))
        sinhlz(ix, iy) = one/SINH(gammaz(ix, iy)*ndimpz*deltaz)
      END DO
    END DO
! WRITE(6,'(a)') 'SINHLZ:'
! WRITE(6,'(4g12.4)') sinhlz
! WRITE(6,'(a)') 'GAMMAZ:'
! WRITE(6,'(4g12.4)') gammaz

    fac = invnormx*invnormy
    DO iz = -maxmz, maxmz
      zu = deltaz*(iz + maxz)
      zl = deltaz*(maxz - iz)
      DO iy = -maxmy, maxmy
        DO ix = -maxmx, maxmx
          potbcz(ix, iy) = fac*sinhlz(ix, iy)*( &
                           potbczl(ix, iy)*sins(gammaz(ix, iy)*zl) &
                           + potbczu(ix, iy)*sins(gammaz(ix, iy)*zu))
! WRITE(6,'(3i4,5(1pg12.4))')
! & ix,iy,iz,potbcz(ix,iy),sinhlz(ix,iy),gammaz(ix,iy),
! & potbczl(ix,iy),potbczu(ix,iy)
        END DO
      END DO
      CALL d2xysinft(potbcz, potbcz)
      DO iy = -maxmy, maxmy
        DO ix = -maxmx, maxmx
          pot(ix, iy, iz) = potbcz(ix, iy) + pot(ix, iy, iz)
        END DO
      END DO
    END DO

! contribution from y-boundaries

    CALL d2xzsinft(potbcyu, potbcyu)
    CALL d2xzsinft(potbcyl, potbcyl)
    DO iz = -maxmz, maxmz
      DO ix = -maxmx, maxmx
        gammay(ix, iz) = SQRT(kinx(ix) + kiny(iz))
        sinhly(ix, iz) = one/SINH(gammay(ix, iz)*ndimpy*deltay)
      END DO
    END DO
! WRITE(6,'(a)') 'SINHLY:'
! WRITE(6,'(4g12.4)') sinhly
! WRITE(6,'(a)') 'GAMMAY:'
! WRITE(6,'(4g12.4)') gammay

    fac = invnormx*invnormz
    DO iy = -maxmy, maxmy
      yu = deltay*(iy + maxy)
      yl = deltay*(maxy - iy)
      DO iz = -maxmz, maxmz
        DO ix = -maxmx, maxmx
          potbcy(ix, iz) = fac*sinhly(ix, iz)*( &
                           potbcyl(ix, iz)*sins(gammay(ix, iz)*yl) &
                           + potbcyu(ix, iz)*sins(gammay(ix, iz)*yu))
        END DO
      END DO
      CALL d2xzsinft(potbcy, potbcy)
      DO iz = -maxmz, maxmz
        DO ix = -maxmx, maxmx
          pot(ix, iy, iz) = potbcy(ix, iz) + pot(ix, iy, iz)
        END DO
      END DO
    END DO

! contribution from x-boundaries

    CALL d2yzsinft(potbcxu, potbcxu)
    CALL d2yzsinft(potbcxl, potbcxl)
    DO iz = -maxmz, maxmz
      DO iy = -maxmy, maxmy
        gammax(iy, iz) = SQRT(kinx(iy) + kiny(iz))
        sinhlx(iy, iz) = one/SINH(gammax(iy, iz)*ndimpx*deltax)
      END DO
    END DO
! WRITE(6,'(a)') 'SINHLX:'
! WRITE(6,'(4g12.4)') sinhlx
! WRITE(6,'(a)') 'GAMMAX:'
! WRITE(6,'(4g12.4)') gammax

    fac = invnormx*invnormy
    DO ix = -maxmx, maxmx
      xu = deltax*(ix + maxx)
      xl = deltax*(maxx - ix)
      DO iz = -maxmz, maxmz
        DO iy = -maxmy, maxmy
          potbcx(iy, iz) = fac*sinhlx(iy, iz)*( &
                           potbcxl(iy, iz)*sins(gammax(iy, iz)*xl) &
                           + potbcxu(iy, iz)*sins(gammax(iy, iz)*xu))
        END DO
      END DO
      CALL d2yzsinft(potbcx, potbcx)
      DO iz = -maxmz, maxmz
        DO iy = -maxmy, maxmy
          pot(ix, iy, iz) = potbcx(iy, iz) + pot(ix, iy, iz)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE homog_sol

!-----boundc-----------------------------------------------------------

  SUBROUTINE boundc(rhoin, deltax, deltay, deltaz, &
                    potbcxl, potbcxu, potbcyl, potbcyu, potbczl, potbczu, &
                    potini, ifinit)

    USE params

! For given density 'rhoin'
! compute multipole moments 'mulmLM'
! and the Coulomb potential at the six boundary planes of the cube.
! The mesh spacing is entered through 'deltax,deltay,deltaz'
! and the mesh SIZE is prescribed IN 'all.inc'.

! Additionally, an initial guess for the potential is
! computed and returned on 'potini'.

! Note: This code is optimized for readability not yet for speed.

    INCLUDE 'work3D.inc'

    REAL, INTENT(IN) :: rhoin(-maxx:maxx, -maxy:maxy, -
    REAL, INTENT(IN) :: deltax
    REAL, INTENT(IN) :: deltay
    REAL, INTENT(IN) :: deltaz
    REAL, INTENT(OUT) :: potbcxl(-maxy:maxy, -maxz:maxz)
    REAL, INTENT(OUT) :: potbcxu(-maxy:maxy, -maxz:maxz)
    REAL, INTENT(OUT) :: potbcyl(-maxx:maxx, -maxz:maxz)
    REAL, INTENT(OUT) :: potbcyu(-maxx:maxx, -maxz:maxz)
    REAL, INTENT(OUT) :: potbczl(-maxx:maxx, -maxy:maxy)
    REAL, INTENT(OUT) :: potbczu(-maxx:maxx, -maxy:maxy)
    REAL, INTENT(OUT) :: potini(-maxx:maxx, -maxy:maxy, -
    INTEGER, INTENT(IN) :: ifinit

! local variables

    LOGICAL :: testpr
    REAL :: mulm00, mulm10, mulm11r, mulm11i, mulm20, mulm21r, mulm21i, mulm22r, mulm22i

    INTEGER :: ix, iy, iz
    REAL :: fac, x, y, z, x2, y2, z2, dr2, r, r2
    REAL :: y00, y10, y11r, y11i, y20, y21r, y21i, y22r, y22i

! damping radius for initialization

    REAL, PARAMETER :: rms2ini = 2.0D0

! forefactors of the spherical harmonics

    REAL, PARAMETER :: pfy00 = 0.282094791
    REAL, PARAMETER :: pfy10 = 0.488602511
    REAL, PARAMETER :: pfy11 = 0.488602511
    REAL, PARAMETER :: pfy20 = 0.315391565
    REAL, PARAMETER :: pfy21 = 1.092548431
    REAL, PARAMETER :: pfy22 = 0.546274215

    DATA testpr/.false./

! the (REAL) spherical harmonics as internal statement functions

    y00(x, y, z) = pfy00
    y10(x, y, z) = pfy10*z
    y11r(x, y, z) = -pfy11*x
    y11i(x, y, z) = -pfy11*y
    y20(x, y, z) = pfy20*(2.0D0*(z*z) - (x*x) - (y*y))
    y21r(x, y, z) = -pfy21*(x*z)
    y21i(x, y, z) = -pfy21*(y*z)
    y22r(x, y, z) = pfy22*((x*x) - (y*y))
    y22i(x, y, z) = pfy22*2.0*(x*y)

!-----------------------------------------------------------------------

! compute the multipole moments

    mulm00 = zero
    mulm10 = zero
    mulm11r = zero
    mulm11i = zero
    mulm20 = zero
    mulm21r = zero
    mulm21i = zero
    mulm22r = zero
    mulm22i = zero

    DO iz = -maxmz, maxmz
      z = deltaz*iz
      DO iy = -maxmy, maxmy
        y = deltay*iy
        DO ix = -maxmx, maxmx
          x = deltax*ix
          fac = rhoin(ix, iy, iz)
          mulm00 = fac*y00(x, y, z) + mulm00
          mulm10 = fac*y10(x, y, z) + mulm10
          mulm11r = fac*y11r(x, y, z) + mulm11r
          mulm11i = fac*y11i(x, y, z) + mulm11i
          mulm20 = fac*y20(x, y, z) + mulm20
          mulm21r = fac*y21r(x, y, z) + mulm21r
          mulm21i = fac*y21i(x, y, z) + mulm21i
          mulm22r = fac*y22r(x, y, z) + mulm22r
          mulm22i = fac*y22i(x, y, z) + mulm22i
        END DO
      END DO
    END DO

    fac = deltax*deltay*deltaz*4.0D0*pi
    mulm00 = fac*mulm00
    mulm10 = fac*mulm10
    mulm11r = fac*mulm11r
    mulm11i = fac*mulm11i
    mulm20 = fac*mulm20
    mulm21r = fac*mulm21r
    mulm21i = fac*mulm21i
    mulm22r = fac*mulm22r
    mulm22i = fac*mulm22i

    IF (testpr) THEN
      WRITE (6, '(a)') ' Multipole moments:'
      WRITE (6, *) mulm00, mulm10, mulm11r, mulm11i, &
        mulm20, mulm21r, mulm21i, mulm22r, mulm22i
    END IF

! z-boundaries

    z = deltaz*maxz
    z2 = z*z
    DO iy = -maxmy, maxmy
      y = deltay*iy
      dr2 = y*y + z2
      DO ix = -maxmx, maxmx
        x = deltax*ix
        r2 = one/(x*x + dr2)
        r = SQRT(r2)
        potbczu(ix, iy) = r*(mulm00*y00(x, y, z) + r2*(mulm10*y10(x, y, z) &
                                                       + mulm11r*y11r(x, y, z) + mulm11i*y11i(x, y, z) &
                                                       + r2*(mulm20*y20(x, y, z) + mulm21r*y21r(x, y, z) &
                                                             + mulm21i*y21i(x, y, z) + mulm22r*y22r(x, y, z) &
                                                             + mulm22i*y22i(x, y, z))))
        potbczl(ix, iy) = r*(mulm00*y00(x, y, -z) + r2*(mulm10*y10(x, y, -z) &
                                                        + mulm11r*y11r(x, y, -z) + mulm11i*y11i(x, y, -z) &
                                                        + r2*(mulm20*y20(x, y, -z) + mulm21r*y21r(x, y, -z) &
                                                              + mulm21i*y21i(x, y, -z) + mulm22r*y22r(x, y, -z) &
                                                              + mulm22i*y22i(x, y, -z))))
      END DO
    END DO

! y-boundaries

    y = deltay*maxy
    y2 = y*y
    DO iz = -maxmz, maxmz
      z = deltaz*iz
      dr2 = z*z + y2
      DO ix = -maxmx, maxmx
        x = deltax*ix
        r2 = one/(x*x + dr2)
        r = SQRT(r2)
        potbcyu(ix, iz) = r*(mulm00*y00(x, y, z) + r2*(mulm10*y10(x, y, z) &
                                                       + mulm11r*y11r(x, y, z) + mulm11i*y11i(x, y, z) &
                                                       + r2*(mulm20*y20(x, y, z) + mulm21r*y21r(x, y, z) &
                                                             + mulm21i*y21i(x, y, z) + mulm22r*y22r(x, y, z) &
                                                             + mulm22i*y22i(x, y, z))))
        potbcyl(ix, iz) = r*(mulm00*y00(x, -y, z) + r2*(mulm10*y10(x, -y, z) &
                                                        + mulm11r*y11r(x, -y, z) + mulm11i*y11i(x, -y, z) &
                                                        + r2*(mulm20*y20(x, -y, z) + mulm21r*y21r(x, -y, z) &
                                                              + mulm21i*y21i(x, -y, z) + mulm22r*y22r(x, -y, z) &
                                                              + mulm22i*y22i(x, -y, z))))
      END DO
    END DO

! x-boundaries

    x = deltax*maxx
    x2 = x*x
    DO iz = -maxmz, maxmz
      z = deltaz*iz
      dr2 = z*z + x2
      DO iy = -maxmy, maxmy
        y = deltay*iy
        r2 = one/(y*y + dr2)
        r = SQRT(r2)
        potbcxu(iy, iz) = r*(mulm00*y00(x, y, z) + r2*(mulm10*y10(x, y, z) &
                                                       + mulm11r*y11r(x, y, z) + mulm11i*y11i(x, y, z) &
                                                       + r2*(mulm20*y20(x, y, z) + mulm21r*y21r(x, y, z) &
                                                             + mulm21i*y21i(x, y, z) + mulm22r*y22r(x, y, z) &
                                                             + mulm22i*y22i(x, y, z))))
        potbcxl(iy, iz) = r*(mulm00*y00(-x, y, z) + r2*(mulm10*y10(-x, y, z) &
                                                        + mulm11r*y11r(-x, y, z) + mulm11i*y11i(-x, y, z) &
                                                        + r2*(mulm20*y20(-x, y, z) + mulm21r*y21r(-x, y, z) &
                                                              + mulm21i*y21i(-x, y, z) + mulm22r*y22r(-x, y, z) &
                                                              + mulm22i*y22i(-x, y, z))))
      END DO
    END DO

    IF (ifinit /= 1) RETURN

! initialization IN interior

    DO iz = -maxmz, maxmz
      z = deltaz*iz
      dr2 = z*z + rms2ini
      DO iy = -maxmy, maxmy
        y = deltay*iy
        DO ix = -maxmx, maxmx
          x = deltax*ix
          r2 = one/(x*x + y*y + dr2)
          r = SQRT(r2)
          potini(ix, iy, iz) = r*(mulm00*y00(x, y, z) + r2*(mulm10*y10(x, y, z) &
                                                            + mulm11r*y11r(x, y, z) + mulm11i*y11i(x, y, z) &
                                                            + r2*(mulm20*y20(x, y, z) + mulm21r*y21r(x, y, z) &
                                                                  + mulm21i*y21i(x, y, z) + mulm22r*y22r(x, y, z) &
                                                                  + mulm22i*y22i(x, y, z))))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE boundc

END MODULE coulsolv
