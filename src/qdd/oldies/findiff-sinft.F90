
MODULE coulsolv
  USE params, ONLY: maxx, maxy, maxz, DP, minx, miny, minz, zero
  SAVE
!!$PRIVATE
!!$!-----work3D------------------------------------------------------------
!!$REAL(DP) :: workx(2*maxx+1)
!!$REAL(DP) :: workxi(2*maxx+1)
!!$REAL(DP) :: wsavex (6*maxx+15)
!!$REAL(DP) :: kinx(-maxx:maxx)
!!$REAL(DP) :: pfacx, pinvfacx, invnormx
!!$INTEGER :: ndimx,ndimpx,maxmx
!!$
!!$REAL(DP) :: worky(2*maxy+1)
!!$REAL(DP) :: workyi(2*maxy+1)
!!$REAL(DP) :: wsavey (6*maxy+15)
!!$REAL(DP) :: kiny(-maxy:maxy)
!!$REAL(DP) :: pfacy, pinvfacy, invnormy
!!$INTEGER :: ndimy,ndimpy,maxmy
!!$
!!$REAL(DP) :: workz(2*maxz+1)
!!$REAL(DP) :: workzi(2*maxz+1)
!!$REAL(DP) :: wsavez (6*maxz+15)
!!$REAL(DP) :: kinz(-maxz:maxz)
!!$REAL(DP) :: pfacz, pinvfacz, invnormz
!!$INTEGER :: ndimz,ndimpz,maxmz

!COMMON /work3d/ workx,worky,workz,wsavex,wsavey,wsavez, &
! workxi,workyi,workzi, kinx,kiny,kinz, &
! pfacx, pinvfacx, invnormx, pfacy, pinvfacy, invnormy, &
! pfacz, pinvfacz, invnormz, ndimx,ndimpx,maxmx, &
! ndimy,ndimpy,maxmy, ndimz,ndimpz,maxmz

  REAL(DP), ALLOCATABLE, PRIVATE :: workx(:), worky(:), workz(:)
  REAL(DP), ALLOCATABLE, PRIVATE :: workxi(:), workyi(:), workzi(:)
  REAL(DP), ALLOCATABLE, PRIVATE :: wsavex(:), wsavey(:), wsavez(:)
  REAL(DP), ALLOCATABLE, PRIVATE :: kinx(:), kiny(:), kinz(:)
  REAL(DP), PRIVATE :: pfacx, pinvfacx, invnormx
  REAL(DP), PRIVATE :: pfacy, pinvfacy, invnormy
  REAL(DP), PRIVATE :: pfacz, pinvfacz, invnormz

  INTEGER, PRIVATE :: ndimx, ndimpx, maxmx
  INTEGER, PRIVATE :: ndimy, ndimpy, maxmy
  INTEGER, PRIVATE :: ndimz, ndimpz, maxmz

CONTAINS

  SUBROUTINE init_findiff_sinft(maxx, maxy, maxz)
    INTEGER, INTENT(IN)::maxx, maxy, maxz
    ALLOCATE (workx(2*maxx + 1), worky(2*maxx + 1), workz(2*maxx + 1))
    ALLOCATE (workxi(2*maxx + 1), workyi(2*maxx + 1), workzi(2*maxx + 1))
    ALLOCATE (wsavex(6*maxx + 15), wsavey(6*maxx + 15), wsavez(6*maxx + 15))
    ALLOCATE (kinx(-maxx:maxx), kiny(-maxy:maxx), kinz(-maxz:maxx))

  ENDSUBROUTINE init_findiff_sinft

  SUBROUTINE solv_fft(rhoin, pot, deltax, deltay, deltaz)

!USE params

! For given density 'rhoin'
! solve Poisson equation with b.c. 0 at infinity
! using sine-FFT.

!INCLUDE 'work3D.inc'

    REAL(DP), INTENT(IN OUT) :: rhoin(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(OUT) :: pot(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(IN OUT) :: deltax
    REAL(DP), INTENT(IN OUT) :: deltay
    REAL(DP), INTENT(IN OUT) :: deltaz

    REAL(DP):: potbcxl(-maxy:maxy, -maxz:maxz)
    REAL(DP):: potbcxu(-maxy:maxy, -maxz:maxz)
    REAL(DP):: potbcyl(-maxx:maxx, -maxz:maxz)
    REAL(DP):: potbcyu(-maxx:maxx, -maxz:maxz)
    REAL(DP):: potbczl(-maxx:maxx, -maxy:maxy)
    REAL(DP):: potbczu(-maxx:maxx, -maxy:maxy)

    REAL(DP):: fac
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

!USE params

! Compute homogeneous solution for given potentials
! 'potbcxl',...,'potbczu' at the six boundary planes of the cube.
! Inhomogeneous solution is entered on 'pot',
! homogeneous solution is added, and RESULT is returned on 'pot'.
! The mesh spacing is entered through 'deltax,deltay,deltaz'
! and the mesh SIZE is prescribed IN 'all.inc'.

! Note: This code is optimized for readability not yet for speed.

!INCLUDE 'work3D.inc'

    REAL(DP), INTENT(OUT) :: pot(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(IN) :: deltax
    REAL(DP), INTENT(IN) :: deltay
    REAL(DP), INTENT(IN) :: deltaz
    REAL(DP), INTENT(IN OUT) :: potbcxl(-maxy:maxy, -maxz:maxz)
    REAL(DP), INTENT(IN OUT) :: potbcxu(-maxy:maxy, -maxz:maxz)
    REAL(DP), INTENT(IN OUT) :: potbcyl(-maxx:maxx, -maxz:maxz)
    REAL(DP), INTENT(IN OUT) :: potbcyu(-maxx:maxx, -maxz:maxz)
    REAL(DP), INTENT(IN OUT) :: potbczl(-maxx:maxx, -maxy:maxy)
    REAL(DP), INTENT(IN OUT) :: potbczu(-maxx:maxx, -maxy:maxy)

! local variables

    LOGICAL :: testpr
    INTEGER :: ix, iy, iz
    REAL(DP):: fac, xl, xu, yl, yu, zl, zu

    REAL(DP):: gammax(-maxy:maxy, -maxz:maxz)
    REAL(DP):: sinhlx(-maxy:maxy, -maxz:maxz)
    REAL(DP):: gammay(-maxx:maxx, -maxz:maxz)
    REAL(DP):: sinhly(-maxx:maxx, -maxz:maxz)
    REAL(DP):: gammaz(-maxx:maxx, -maxy:maxy)
    REAL(DP):: sinhlz(-maxx:maxx, -maxy:maxy)
    REAL(DP):: potbcx(-maxy:maxy, -maxz:maxz)
    REAL(DP):: potbcy(-maxx:maxx, -maxz:maxz)
    REAL(DP):: potbcz(-maxx:maxx, -maxy:maxy)

    DATA testpr/.false./

! stabilised sinh as internal statement FUNCTION

!sins(x) = SINH(SIGN(MIN(ABS(x),72.0D0),x))

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

!USE params

! For given density 'rhoin'
! compute multipole moments 'mulmLM'
! and the Coulomb potential at the six boundary planes of the cube.
! The mesh spacing is entered through 'deltax,deltay,deltaz'
! and the mesh SIZE is prescribed IN 'all.inc'.

! Additionally, an initial guess for the potential is
! computed and returned on 'potini'.

! Note: This code is optimized for readability not yet for speed.

!INCLUDE 'work3D.inc'

    REAL(DP), INTENT(IN) :: rhoin(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(IN) :: deltax
    REAL(DP), INTENT(IN) :: deltay
    REAL(DP), INTENT(IN) :: deltaz
    REAL(DP), INTENT(OUT) :: potbcxl(-maxy:maxy, -maxz:maxz)
    REAL(DP), INTENT(OUT) :: potbcxu(-maxy:maxy, -maxz:maxz)
    REAL(DP), INTENT(OUT) :: potbcyl(-maxx:maxx, -maxz:maxz)
    REAL(DP), INTENT(OUT) :: potbcyu(-maxx:maxx, -maxz:maxz)
    REAL(DP), INTENT(OUT) :: potbczl(-maxx:maxx, -maxy:maxy)
    REAL(DP), INTENT(OUT) :: potbczu(-maxx:maxx, -maxy:maxy)
    REAL(DP), INTENT(OUT) :: potini(-maxx:maxx, -maxy:maxy, minz:maxz)
    INTEGER, INTENT(IN) :: ifinit

! local variables

    LOGICAL :: testpr
    REAL(DP):: mulm00, mulm10, mulm11r, mulm11i, mulm20, mulm21r, mulm21i, mulm22r, mulm22i

    INTEGER :: ix, iy, iz
    REAL(DP):: fac, x, y, z, x2, y2, z2, dr2, r, r2
    REAL(DP):: y00, y10, y11r, y11i, y20, y21r, y21i, y22r, y22i

! damping radius for initialization

    REAL(DP), PARAMETER :: rms2ini = 2.0D0

! forefactors of the spherical harmonics

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

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
    y22i(x, y, z) = pfy22*2.0D0*(x*y)

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

!-----d3sinfinit----------------------------------------------------------

  SUBROUTINE d3sinfinit(deltax, deltay, deltaz)

!USE params
! include 'param3D.inc'

! Initialization of Sine FFT IN 3D.
! 'deltax', 'deltay', and 'deltaz' are spatial mesh sizes.
! Actual grid sizes are given IN 'param3D.inc'.
! Resulting work spaces and deduced parameters are returned
! via 'work3D.inc'.

! the storage of the kinetic energy IN x-, y-, and z
! is shifted TO a negative left bound. This is TO accommodate
! TO loop parameters IN an actual all.
! Note that p-space counts from 1,...,ndimx
! WHILE x-space was from -maxmx,...,+maxmx

    REAL(DP), INTENT(IN OUT) :: deltax
    REAL(DP), INTENT(IN OUT) :: deltay
    REAL(DP), INTENT(IN OUT) :: deltaz

!INCLUDE 'work3D.inc'

    INTEGER :: n

!-----------------------------------------------------------------------

    ndimpx = 2*maxx
    ndimx = ndimpx - 1
    maxmx = maxx - 1
    pfacx = (pi/(deltax*DBLE(ndimpx)))**2 ! second derivative
    invnormx = half/DBLE(ndimpx) ! normalization
    CALL sinti(ndimx, wsavex)
    DO n = 1, ndimx
      kinx(n - maxx) = DBLE(n)**2*pfacx
!test WRITE(*,*) kinx(n-MAXX)
    END DO

    ndimpy = 2*maxy
    ndimy = ndimpy - 1
    maxmy = maxy - 1
    pfacy = (pi/(deltay*DBLE(ndimpy)))**2 ! second derivative
    invnormy = half/DBLE(ndimpy) ! normalization
    CALL sinti(ndimy, wsavey)
    DO n = 1, ndimy
      kiny(n - maxy) = DBLE(n)**2*pfacy
    END DO

    ndimpz = 2*maxz
    ndimz = ndimpz - 1
    maxmz = maxz - 1
    pfacz = (pi/(deltaz*DBLE(ndimpz)))**2 ! second derivative
    invnormz = half/DBLE(ndimpz) ! normalization
    CALL sinti(ndimz, wsavez)
    DO n = 1, ndimz
      kinz(n - maxz) = DBLE(n)**2*pfacz
    END DO

    WRITE (*, '(A,3I4)') &
      ' sine-FFT initialized: x-y-z dimensions =', ndimx, ndimy, ndimz

    RETURN
  END SUBROUTINE d3sinfinit

!-----d3sinfft----------------------------------------------------------

  SUBROUTINE d3sinft(fieldin, fieldout)

!USE params
! include 'param3D.inc'

! Sine FFT IN 3D by repeated calls TO 1D transformation
! 'fieldin' = input array
! 'fieldout' = output array
! Both fields can be identical IN the calling routine.
! Both array are zero at the bounds.
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    COMPLEX(DP), INTENT(IN OUT) :: fieldin(-maxx:maxx, -maxy:maxy, minz:maxz)
    COMPLEX(DP), INTENT(OUT) :: fieldout(-maxx:maxx, -maxy:maxy, minz:maxz)

!INCLUDE 'work3D.inc'

!old INTEGER klt
    INTEGER :: ix, iy, iz, n

!-----------------------------------------------------------------------

! x FFT

    DO iz = -maxmz, maxmz
      DO iy = -maxmy, maxmy
        n = 0
        DO ix = -maxmx, maxmx
          n = 1 + n
          workx(n) = REAL(fieldin(ix, iy, iz), DP)
          workxi(n) = AIMAG(fieldin(ix, iy, iz))
        END DO

! DO klt=1,ndimx
! WRITE(*,*) '1 ON',ndimx,workx(klt),Wsavex(klt)
! END DO

        CALL sint(ndimx, workx(1), wsavex)
        CALL sint(ndimx, workxi(1), wsavex)
! WRITE(*,*) '1 OFF'
        n = 0
        DO ix = -maxmx, maxmx
          n = 1 + n
          fieldout(ix, iy, iz) = CMPLX(workx(n), workxi(n), DP)
        END DO
      END DO
    END DO

! y FFT

    DO iz = -maxmz, maxmz
      DO ix = -maxmy, maxmy
        n = 0
        DO iy = -maxmy, maxmy
          n = 1 + n
          worky(n) = REAL(fieldout(ix, iy, iz), DP)
          workyi(n) = AIMAG(fieldout(ix, iy, iz))
        END DO
! WRITE(*,*) '2 ON'

! DO klt=1,ndimy
! WRITE(*,*) '2 ON',ndimy,worky(klt),Wsavey(klt)
! END DO

        CALL sint(ndimy, worky(1), wsavey)
        CALL sint(ndimy, workyi(1), wsavey)
! WRITE(*,*) '2 OFF'
        n = 0
        DO iy = -maxmy, maxmy
          n = 1 + n
          fieldout(ix, iy, iz) = CMPLX(worky(n), workyi(n), DP)
        END DO
      END DO
    END DO

! z FFT

    DO iy = -maxmy, maxmy
      DO ix = -maxmx, maxmx
        n = 0
        DO iz = -maxmz, maxmz
          n = 1 + n
          workz(n) = REAL(fieldout(ix, iy, iz), DP)
          workzi(n) = AIMAG(fieldout(ix, iy, iz))
        END DO
! WRITE(*,*) '3 ON'
        CALL sint(ndimz, workz(1), wsavez)
        CALL sint(ndimz, workzi(1), wsavez)
! WRITE(*,*) '3 OFF'
        n = 0
        DO iz = -maxmz, maxmz
          n = 1 + n
          fieldout(ix, iy, iz) = CMPLX(workz(n), workzi(n), DP)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE d3sinft

!-----d3sinflaplace----------------------------------------------------------

  SUBROUTINE kin3d_fr(field, kin_field)

!USE params
! include 'param3D.inc'

! ACTION of -Laplace OPERATOR by sine FFT.
! 'field' = input array
! 'kin_field' = resulting array
! Both array are zero at the bounds (as required by sine FFT)
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.
! Other grid parameters and momenta are transferred via 'work3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    COMPLEX(DP), INTENT(IN OUT) :: field(-maxx:maxx, -maxy:maxy, minz:maxz)
    COMPLEX(DP), INTENT(OUT) :: kin_field(-maxx:maxx, -maxy:maxy, minz:maxz)

!INCLUDE 'work3D.inc'

    INTEGER :: ix, iy, iz
    REAL(DP):: laplz, laplzy
    REAL(DP):: totalnorm

!-----------------------------------------------------------------------

! forward FFT

    CALL d3sinft(field, kin_field)

! ACTION of Laplacian

    totalnorm = invnormx*invnormy*invnormz
    DO iz = -maxmz, maxmz
      laplz = kinz(iz)
      DO iy = -maxmy, maxmy
        laplzy = laplz + kiny(iy)
        DO ix = -maxmx, maxmx
          kin_field(ix, iy, iz) = totalnorm*(laplzy + kinx(ix))*kin_field(ix, iy, iz)
        END DO
      END DO
    END DO

! backward FFT

    CALL d3sinft(kin_field, kin_field)

    RETURN
  END SUBROUTINE kin3d_fr

!-----d3sinfpropag----------------------------------------------------------

  SUBROUTINE d3sinfpropag(field, prop_field, del_time)

!USE params
! include 'param3D.inc'

! ACTION of kinetic propagator by sine FFT:

! prop_field = exp(-i*del_time*T) * field

! Lst parameters are:
! 'field' = input array , is destroyed
! 'inv_field' = resulting array , may be identical with input
! 'e0inv' = offset for inversion
! Both array are zero at the bounds (as required by sine FFT)
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.
! Other grid parameters and momenta are transferred via 'work3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    COMPLEX(DP), INTENT(IN OUT) :: field(-maxx:maxx, -maxy:maxy, minz:maxz)
    COMPLEX(DP), INTENT(OUT) :: prop_field(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(IN OUT) :: del_time

!INCLUDE 'work3D.inc'

    INTEGER :: ix, iy, iz
    REAL(DP):: kinactz, kinactzy
    REAL(DP):: totalnorm

!-----------------------------------------------------------------------

! WRITE(*,*) e0inv

! forward FFT

    CALL d3sinft(field, prop_field)

! ACTION of Laplacian

    totalnorm = invnormx*invnormy*invnormz
    DO iz = -maxmz, maxmz
      kinactz = kinz(iz)
      DO iy = -maxmy, maxmy
        kinactzy = kinactz + kiny(iy)
        DO ix = -maxmx, maxmx
          prop_field(ix, iy, iz) = totalnorm &
                                   *EXP(CMPLX(zero, -del_time*(kinactzy + kinx(ix)), DP))*prop_field(ix, iy, iz)
        END DO
      END DO
    END DO

! backward FFT

    CALL d3sinft(prop_field, prop_field)

    RETURN
  END SUBROUTINE d3sinfpropag

!-----d3sinfft----------------------------------------------------------

  SUBROUTINE d3sinftreal(fieldin, fieldout)

!USE params
! include 'param3D.inc'

! Sine FFT IN 3D by repeated calls TO 1D transformation
! 'fieldin' = input array
! 'fieldout' = output array
! Both fields can be identical IN the calling routine.
! Both array are zero at the bounds.
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    REAL(DP), INTENT(IN) :: fieldin(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(OUT) :: fieldout(-maxx:maxx, -maxy:maxy, minz:maxz)

!INCLUDE 'work3D.inc'

    INTEGER :: n
    INTEGER :: ix, iy, iz

!-----------------------------------------------------------------------

! x FFT

    DO iz = -maxmz, maxmz
      DO iy = -maxmy, maxmy
        n = 0
        DO ix = -maxmx, maxmx
          n = 1 + n
          workx(n) = fieldin(ix, iy, iz)
        END DO
        CALL sint(ndimx, workx(1), wsavex)
        n = 0
        DO ix = -maxmx, maxmx
          n = 1 + n
          fieldout(ix, iy, iz) = workx(n)
        END DO
      END DO
    END DO

! y FFT

    DO iz = -maxmz, maxmz
      DO ix = -maxmy, maxmy
        n = 0
        DO iy = -maxmy, maxmy
          n = 1 + n
          worky(n) = fieldout(ix, iy, iz)
        END DO
        CALL sint(ndimy, worky(1), wsavey)
        n = 0
        DO iy = -maxmy, maxmy
          n = 1 + n
          fieldout(ix, iy, iz) = worky(n)
        END DO
      END DO
    END DO

! z FFT

    DO iy = -maxmy, maxmy
      DO ix = -maxmx, maxmx
        n = 0
        DO iz = -maxmz, maxmz
          n = 1 + n
          workz(n) = fieldout(ix, iy, iz)
        END DO
        CALL sint(ndimz, workz(1), wsavez)
        n = 0
        DO iz = -maxmz, maxmz
          n = 1 + n
          fieldout(ix, iy, iz) = workz(n)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE d3sinftreal

!-----d3sinflaplace----------------------------------------------------------

  SUBROUTINE d3sinflaplace(field, lapl_field)

!USE params
! include 'param3D.inc'

! ACTION of -Laplace OPERATOR by sine FFT.
! 'field' = input array
! 'lapl_field' = resulting array
! Both array are zero at the bounds (as required by sine FFT)
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.
! Other grid parameters and momenta are transferred via 'work3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    REAL(DP), INTENT(IN OUT) :: field(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(OUT) :: lapl_field(-maxx:maxx, -maxy:maxy, minz:maxz)

!INCLUDE 'work3D.inc'

    INTEGER :: ix, iy, iz
    REAL(DP):: laplz, laplzy
    REAL(DP):: totalnorm

!-----------------------------------------------------------------------

! forward FFT

    CALL d3sinftreal(field, lapl_field)

! ACTION of Laplacian

    totalnorm = invnormx*invnormy*invnormz
    DO iz = -maxmz, maxmz
      laplz = kinz(iz)
      DO iy = -maxmy, maxmy
        laplzy = laplz + kiny(iy)
        DO ix = -maxmx, maxmx
          lapl_field(ix, iy, iz) = totalnorm*(laplzy + kinx(ix))*field(ix, iy, iz)
        END DO
      END DO
    END DO

! backward FFT

    CALL d3sinftreal(lapl_field, lapl_field)

    RETURN
  END SUBROUTINE d3sinflaplace

!-----d3sinfinverse----------------------------------------------------------

  SUBROUTINE d3sinfinverse(field, inv_field, e0inv)

!USE params
! include 'param3D.inc'

! ACTION of inverse Laplace OPERATOR by sine FFT:

! inv_field = ( -Laplace + e0inv ) * field

! Lst parameters are:
! 'field' = input array , is destroyed
! 'inv_field' = resulting array , may be identical with input
! 'e0inv' = offset for inversion
! Both array are zero at the bounds (as required by sine FFT)
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.
! Other grid parameters and momenta are transferred via 'work3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    REAL(DP), INTENT(IN OUT) :: field(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(OUT) :: inv_field(-maxx:maxx, -maxy:maxy, minz:maxz)
    REAL(DP), INTENT(IN) :: e0inv

!INCLUDE 'work3D.inc'

    INTEGER :: ix, iy, iz
    REAL(DP):: kinactz, kinactzy
    REAL(DP):: totalnorm

!-----------------------------------------------------------------------

! forward FFT

    CALL d3sinftreal(field, inv_field)

! ACTION of Laplacian

    totalnorm = invnormx*invnormy*invnormz
    DO iz = -maxmz, maxmz
      kinactz = kinz(iz) + e0inv
      DO iy = -maxmy, maxmy
        kinactzy = kinactz + kiny(iy)
        DO ix = -maxmx, maxmx
          inv_field(ix, iy, iz) = totalnorm/(kinactzy + kinx(ix))*inv_field(ix, iy, iz)
        END DO
      END DO
    END DO

! backward FFT

    CALL d3sinftreal(inv_field, inv_field)

    RETURN
  END SUBROUTINE d3sinfinverse

!-----d2xysinft----------------------------------------------------------

  SUBROUTINE d2xysinft(fieldin, fieldout)

!USE params
! include 'param3D.inc'

! Sine FFT IN 2D (x,y) by repeated calls TO 1D transformation
! 'fieldin' = input array
! 'fieldout' = output array
! Both fields can be identical IN the calling routine.
! Both array are zero at the bounds.
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    REAL(DP), INTENT(IN) :: fieldin(-maxx:maxx, -maxy:maxy)
    REAL(DP), INTENT(OUT) :: fieldout(-maxx:maxx, -maxy:maxy)

!INCLUDE 'work3D.inc'

    INTEGER :: n
    INTEGER :: ix, iy, iz

!-----------------------------------------------------------------------

! x FFT

    DO iy = -maxmy, maxmy
      n = 0
      DO ix = -maxmx, maxmx
        n = 1 + n
        workx(n) = fieldin(ix, iy)
      END DO
      CALL sint(ndimx, workx(1), wsavex)
      n = 0
      DO ix = -maxmx, maxmx
        n = 1 + n
        fieldout(ix, iy) = workx(n)
      END DO
    END DO

! y FFT

    DO ix = -maxmy, maxmy
      n = 0
      DO iy = -maxmy, maxmy
        n = 1 + n
        worky(n) = fieldout(ix, iy)
      END DO
      CALL sint(ndimy, worky(1), wsavey)
      n = 0
      DO iy = -maxmy, maxmy
        n = 1 + n
        fieldout(ix, iy) = worky(n)
      END DO
    END DO

    RETURN
  END SUBROUTINE d2xysinft

!-----d2xzsinft----------------------------------------------------------

  SUBROUTINE d2xzsinft(fieldin, fieldout)

!USE params
! include 'param3D.inc'

! Sine FFT IN 2D (x,z) by repeated calls TO 1D transformation
! 'fieldin' = input array
! 'fieldout' = output array
! Both fields can be identical IN the calling routine.
! Both array are zero at the bounds.
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    REAL(DP), INTENT(IN) :: fieldin(-maxx:maxx, -maxz:maxz)
    REAL(DP), INTENT(OUT) :: fieldout(-maxx:maxx, -maxz:maxz)

!INCLUDE 'work3D.inc'

    INTEGER :: n
    INTEGER :: ix, iy, iz

!-----------------------------------------------------------------------

! x FFT

    DO iz = -maxmz, maxmz
      n = 0
      DO ix = -maxmx, maxmx
        n = 1 + n
        workx(n) = fieldin(ix, iz)
      END DO
      CALL sint(ndimx, workx(1), wsavex)
      n = 0
      DO ix = -maxmx, maxmx
        n = 1 + n
        fieldout(ix, iz) = workx(n)
      END DO
    END DO

! z FFT

    DO ix = -maxmx, maxmx
      n = 0
      DO iz = -maxmz, maxmz
        n = 1 + n
        workz(n) = fieldout(ix, iz)
      END DO
      CALL sint(ndimz, workz(1), wsavez)
      n = 0
      DO iz = -maxmz, maxmz
        n = 1 + n
        fieldout(ix, iz) = workz(n)
      END DO
    END DO

    RETURN
  END SUBROUTINE d2xzsinft

!-----d2yzsinft----------------------------------------------------------

  SUBROUTINE d2yzsinft(fieldin, fieldout)

!USE params
! include 'param3D.inc'

! Sine FFT IN 2D (y,z) by repeated calls TO 1D transformation
! 'fieldin' = input array
! 'fieldout' = output array
! Both fields can be identical IN the calling routine.
! Both array are zero at the bounds.
! Actual bounds are identical with dimensions as fixed IN 'param3D.inc'.

! This routine requires that all work spaces have been initialized
! by a previous CALL TO 'd3sinfinit'.

    REAL(DP), INTENT(IN OUT) :: fieldin(-maxy:maxy, -maxz:maxz)
    REAL(DP), INTENT(IN OUT) :: fieldout(-maxy:maxy, -maxz:maxz)

!INCLUDE 'work3D.inc'

    INTEGER :: n
    INTEGER :: ix, iy, iz

!-----------------------------------------------------------------------

! y FFT

    DO iz = -maxmz, maxmz
      n = 0
      DO iy = -maxmy, maxmy
        n = 1 + n
        worky(n) = fieldout(iy, iz)
      END DO
      CALL sint(ndimy, worky(1), wsavey)
      n = 0
      DO iy = -maxmy, maxmy
        n = 1 + n
        fieldout(iy, iz) = worky(n)
      END DO
    END DO

! z FFT

    DO iy = -maxmy, maxmy
      n = 0
      DO iz = -maxmz, maxmz
        n = 1 + n
        workz(n) = fieldout(iy, iz)
      END DO
      CALL sint(ndimz, workz(1), wsavez)
      n = 0
      DO iz = -maxmz, maxmz
        n = 1 + n
        fieldout(iy, iz) = workz(n)
      END DO
    END DO

    RETURN
  END SUBROUTINE d2yzsinft

!---------------- sine FFT from Netlib ---------------------------------

  SUBROUTINE sinti(n, wsave)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(OUT) :: wsave(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    DATA pi/3.14159265358979D0/

    IF (n <= 1) RETURN
    ns2 = n/2
    np1 = n + 1
    dt = pi/DBLE(np1)
    DO k = 1, ns2
      wsave(k) = 2.0D0*SIN(k*dt)
    END DO
    CALL rffti(np1, wsave(ns2 + 1))
    RETURN
  END SUBROUTINE sinti

  SUBROUTINE rffti(n, wsave)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN OUT) :: wsave(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    IF (n == 1) RETURN
    CALL rffti1(n, wsave(n + 1), wsave(2*n + 1))
    RETURN
  END SUBROUTINE rffti

  SUBROUTINE rffti1(n, wa, ifac)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(OUT) :: wa(*)
    INTEGER, INTENT(OUT) :: ifac(3)
!IMPLICIT REAL(DP)(a-h,o-z)
    INTEGER :: ntryh(4)
    DATA ntryh(1), ntryh(2), ntryh(3), ntryh(4)/4, 2, 3, 5/

    nl = n
    nf = 0
    j = 0
101 j = j + 1
    IF (j - 4 > 0) THEN
      GO TO 103
    END IF
102 ntry = ntryh(j)
    GO TO 104
103 ntry = ntry + 2
104 nq = nl/ntry
    nr = nl - ntry*nq
    IF (nr == 0) THEN
      GO TO 105
    ELSE
      GO TO 101
    END IF
105 nf = nf + 1
    ifac(nf + 2) = ntry
    nl = nq
    IF (ntry /= 2) GO TO 107
    IF (nf == 1) GO TO 107
    DO i = 2, nf
      ib = nf - i + 2
      ifac(ib + 2) = ifac(ib + 1)
    END DO
    ifac(3) = 2
107 IF (nl /= 1) GO TO 104
    ifac(1) = n
    ifac(2) = nf
    tpi = 6.28318530717959D0
    argh = tpi/DBLE(n)
    is = 0
    nfm1 = nf - 1
    l1 = 1
    IF (nfm1 == 0) RETURN
    DO k1 = 1, nfm1
      ip = ifac(k1 + 2)
      ld = 0
      l2 = l1*ip
      ido = n/l2
      ipm = ip - 1
      DO j = 1, ipm
        ld = ld + l1
        i = is
        argld = DBLE(ld)*argh
        fi = 0.
        DO ii = 3, ido, 2
          i = i + 2
          fi = fi + 1.
          arg = fi*argld
          wa(i - 1) = COS(arg)
          wa(i) = SIN(arg)
        END DO
        is = is + ido
      END DO
      l1 = l2
    END DO
    RETURN
  END SUBROUTINE rffti1

  SUBROUTINE sint(n, x, wsave)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN OUT) :: x(1)
    REAL(DP), INTENT(IN OUT) :: wsave(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    np1 = n + 1
    iw1 = n/2 + 1
    iw2 = iw1 + np1
    iw3 = iw2 + np1
    CALL sint1(n, x, wsave, wsave(iw1), wsave(iw2), wsave(iw3))
    RETURN
  END SUBROUTINE sint

  SUBROUTINE radf2(ido, l1, cc, ch, wa1)

    INTEGER, INTENT(IN) :: ido
    INTEGER, INTENT(IN) :: l1
    REAL(DP), INTENT(IN) :: cc(ido, l1, 2)
    REAL(DP), INTENT(OUT) :: ch(ido, 2, l1)
    REAL(DP), INTENT(IN) :: wa1(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    DO k = 1, l1
      ch(1, 1, k) = cc(1, k, 1) + cc(1, k, 2)
      ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 2)
    END DO
    IF (ido - 2 < 0) THEN
      GO TO 107
    ELSE IF (ido - 2 == 0) THEN
      GO TO 105
    END IF
102 idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        tr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        ti2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        ch(i, 1, k) = cc(i, k, 1) + ti2
        ch(ic, 2, k) = ti2 - cc(i, k, 1)
        ch(i - 1, 1, k) = cc(i - 1, k, 1) + tr2
        ch(ic - 1, 2, k) = cc(i - 1, k, 1) - tr2
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
105 DO k = 1, l1
      ch(1, 2, k) = -cc(ido, k, 2)
      ch(ido, 1, k) = cc(ido, k, 1)
    END DO
107 RETURN
  END SUBROUTINE radf2

  SUBROUTINE radf3(ido, l1, cc, ch, wa1, wa2)

    INTEGER, INTENT(IN) :: ido
    INTEGER, INTENT(IN) :: l1
    REAL(DP), INTENT(IN) :: cc(ido, l1, 3)
    REAL(DP), INTENT(OUT) :: ch(ido, 3, l1)
    REAL(DP), INTENT(IN) :: wa1(1)
    REAL(DP), INTENT(IN) :: wa2(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    DATA taur, taui/-.5D0, .866025403784439D0/

    DO k = 1, l1
      cr2 = cc(1, k, 2) + cc(1, k, 3)
      ch(1, 1, k) = cc(1, k, 1) + cr2
      ch(1, 3, k) = taui*(cc(1, k, 3) - cc(1, k, 2))
      ch(ido, 2, k) = cc(1, k, 1) + taur*cr2
    END DO
    IF (ido == 1) RETURN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
        di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
        cr2 = dr2 + dr3
        ci2 = di2 + di3
        ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2
        ch(i, 1, k) = cc(i, k, 1) + ci2
        tr2 = cc(i - 1, k, 1) + taur*cr2
        ti2 = cc(i, k, 1) + taur*ci2
        tr3 = taui*(di2 - di3)
        ti3 = taui*(dr3 - dr2)
        ch(i - 1, 3, k) = tr2 + tr3
        ch(ic - 1, 2, k) = tr2 - tr3
        ch(i, 3, k) = ti2 + ti3
        ch(ic, 2, k) = ti3 - ti2
      END DO
    END DO
    RETURN
  END SUBROUTINE radf3

  SUBROUTINE radf4(ido, l1, cc, ch, wa1, wa2, wa3)

    INTEGER, INTENT(IN) :: ido
    INTEGER, INTENT(IN) :: l1
    REAL(DP), INTENT(IN) :: cc(ido, l1, 4)
    REAL(DP), INTENT(OUT) :: ch(ido, 4, l1)
    REAL(DP), INTENT(IN) :: wa1(1)
    REAL(DP), INTENT(IN) :: wa2(1)
    REAL(DP), INTENT(IN) :: wa3(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    DATA hsqt2/.7071067811865475D0/

    DO k = 1, l1
      tr1 = cc(1, k, 2) + cc(1, k, 4)
      tr2 = cc(1, k, 1) + cc(1, k, 3)
      ch(1, 1, k) = tr1 + tr2
      ch(ido, 4, k) = tr2 - tr1
      ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 3)
      ch(1, 3, k) = cc(1, k, 4) - cc(1, k, 2)
    END DO
    IF (ido - 2 < 0) THEN
      GO TO 107
    ELSE IF (ido - 2 == 0) THEN
      GO TO 105
    END IF
102 idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        cr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        ci2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        cr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
        ci3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
        cr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
        ci4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
        tr1 = cr2 + cr4
        tr4 = cr4 - cr2
        ti1 = ci2 + ci4
        ti4 = ci2 - ci4
        ti2 = cc(i, k, 1) + ci3
        ti3 = cc(i, k, 1) - ci3
        tr2 = cc(i - 1, k, 1) + cr3
        tr3 = cc(i - 1, k, 1) - cr3
        ch(i - 1, 1, k) = tr1 + tr2
        ch(ic - 1, 4, k) = tr2 - tr1
        ch(i, 1, k) = ti1 + ti2
        ch(ic, 4, k) = ti1 - ti2
        ch(i - 1, 3, k) = ti4 + tr3
        ch(ic - 1, 2, k) = tr3 - ti4
        ch(i, 3, k) = tr4 + ti3
        ch(ic, 2, k) = tr4 - ti3
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
105 CONTINUE
    DO k = 1, l1
      ti1 = -hsqt2*(cc(ido, k, 2) + cc(ido, k, 4))
      tr1 = hsqt2*(cc(ido, k, 2) - cc(ido, k, 4))
      ch(ido, 1, k) = tr1 + cc(ido, k, 1)
      ch(ido, 3, k) = cc(ido, k, 1) - tr1
      ch(1, 2, k) = ti1 - cc(ido, k, 3)
      ch(1, 4, k) = ti1 + cc(ido, k, 3)
    END DO
107 RETURN
  END SUBROUTINE radf4

  SUBROUTINE radf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)

    INTEGER, INTENT(IN) :: ido
    INTEGER, INTENT(IN) :: l1
    REAL(DP), INTENT(IN) :: cc(ido, l1, 5)
    REAL(DP), INTENT(OUT) :: ch(ido, 5, l1)
    REAL(DP), INTENT(IN) :: wa1(1)
    REAL(DP), INTENT(IN) :: wa2(1)
    REAL(DP), INTENT(IN) :: wa3(1)
    REAL(DP), INTENT(IN) :: wa4(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    DATA tr11, ti11, tr12, ti12/.309016994374947D0, .951056516295154D0, &
      -.809016994374947D0, .587785252292473D0/

    DO k = 1, l1
      cr2 = cc(1, k, 5) + cc(1, k, 2)
      ci5 = cc(1, k, 5) - cc(1, k, 2)
      cr3 = cc(1, k, 4) + cc(1, k, 3)
      ci4 = cc(1, k, 4) - cc(1, k, 3)
      ch(1, 1, k) = cc(1, k, 1) + cr2 + cr3
      ch(ido, 2, k) = cc(1, k, 1) + tr11*cr2 + tr12*cr3
      ch(1, 3, k) = ti11*ci5 + ti12*ci4
      ch(ido, 4, k) = cc(1, k, 1) + tr12*cr2 + tr11*cr3
      ch(1, 5, k) = ti12*ci5 - ti11*ci4
    END DO
    IF (ido == 1) RETURN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
        di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
        dr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
        di4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
        dr5 = wa4(i - 2)*cc(i - 1, k, 5) + wa4(i - 1)*cc(i, k, 5)
        di5 = wa4(i - 2)*cc(i, k, 5) - wa4(i - 1)*cc(i - 1, k, 5)
        cr2 = dr2 + dr5
        ci5 = dr5 - dr2
        cr5 = di2 - di5
        ci2 = di2 + di5
        cr3 = dr3 + dr4
        ci4 = dr4 - dr3
        cr4 = di3 - di4
        ci3 = di3 + di4
        ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2 + cr3
        ch(i, 1, k) = cc(i, k, 1) + ci2 + ci3
        tr2 = cc(i - 1, k, 1) + tr11*cr2 + tr12*cr3
        ti2 = cc(i, k, 1) + tr11*ci2 + tr12*ci3
        tr3 = cc(i - 1, k, 1) + tr12*cr2 + tr11*cr3
        ti3 = cc(i, k, 1) + tr12*ci2 + tr11*ci3
        tr5 = ti11*cr5 + ti12*cr4
        ti5 = ti11*ci5 + ti12*ci4
        tr4 = ti12*cr5 - ti11*cr4
        ti4 = ti12*ci5 - ti11*ci4
        ch(i - 1, 3, k) = tr2 + tr5
        ch(ic - 1, 2, k) = tr2 - tr5
        ch(i, 3, k) = ti2 + ti5
        ch(ic, 2, k) = ti5 - ti2
        ch(i - 1, 5, k) = tr3 + tr4
        ch(ic - 1, 4, k) = tr3 - tr4
        ch(i, 5, k) = ti3 + ti4
        ch(ic, 4, k) = ti4 - ti3
      END DO
    END DO
    RETURN
  END SUBROUTINE radf5

  SUBROUTINE radfg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)

    INTEGER, INTENT(IN) :: ido
    INTEGER, INTENT(IN) :: ip
    INTEGER, INTENT(IN) :: l1
    INTEGER, INTENT(IN) :: idl1
    REAL(DP), INTENT(OUT) :: cc(ido, ip, l1)
    REAL(DP), INTENT(OUT) :: c1(ido, l1, ip)
    REAL(DP), INTENT(OUT) :: c2(idl1, ip)
    REAL(DP), INTENT(OUT) :: ch(ido, l1, ip)
    REAL(DP), INTENT(OUT) :: ch2(idl1, ip)
    REAL(DP), INTENT(IN) :: wa(*)
!IMPLICIT REAL(DP)(a-h,o-z)

    DATA tpi/6.28318530717959D0/

    arg = tpi/DBLE(ip)
    dcp = COS(arg)
    dsp = SIN(arg)
    ipph = (ip + 1)/2
    ipp2 = ip + 2
    idp2 = ido + 2
    nbd = (ido - 1)/2
    IF (ido == 1) GO TO 119
    DO ik = 1, idl1
      ch2(ik, 1) = c2(ik, 1)
    END DO
    DO j = 2, ip
      DO k = 1, l1
        ch(1, k, j) = c1(1, k, j)
      END DO
    END DO
    IF (nbd > l1) GO TO 107
    is = -ido
    DO j = 2, ip
      is = is + ido
      idij = is
      DO i = 3, ido, 2
        idij = idij + 2
        DO k = 1, l1
          ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij)*c1(i, k, j)
          ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij)*c1(i - 1, k, j)
        END DO
      END DO
    END DO
    GO TO 111
107 is = -ido
    DO j = 2, ip
      is = is + ido
      DO k = 1, l1
        idij = is
        DO i = 3, ido, 2
          idij = idij + 2
          ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij)*c1(i, k, j)
          ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij)*c1(i - 1, k, j)
        END DO
      END DO
    END DO
111 IF (nbd < l1) GO TO 115
    DO j = 2, ipph
      jc = ipp2 - j
      DO k = 1, l1
        DO i = 3, ido, 2
          c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
          c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
          c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
        END DO
      END DO
    END DO
    GO TO 121
115 DO j = 2, ipph
      jc = ipp2 - j
      DO i = 3, ido, 2
        DO k = 1, l1
          c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
          c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
          c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
        END DO
      END DO
    END DO
    GO TO 121
119 DO ik = 1, idl1
      c2(ik, 1) = ch2(ik, 1)
    END DO
121 DO j = 2, ipph
      jc = ipp2 - j
      DO k = 1, l1
        c1(1, k, j) = ch(1, k, j) + ch(1, k, jc)
        c1(1, k, jc) = ch(1, k, jc) - ch(1, k, j)
      END DO
    END DO

    ar1 = 1D0
    ai1 = 0D0
    DO l = 2, ipph
      lc = ipp2 - l
      ar1h = dcp*ar1 - dsp*ai1
      ai1 = dcp*ai1 + dsp*ar1
      ar1 = ar1h
      DO ik = 1, idl1
        ch2(ik, l) = c2(ik, 1) + ar1*c2(ik, 2)
        ch2(ik, lc) = ai1*c2(ik, ip)
      END DO
      dc2 = ar1
      ds2 = ai1
      ar2 = ar1
      ai2 = ai1
      DO j = 3, ipph
        jc = ipp2 - j
        ar2h = dc2*ar2 - ds2*ai2
        ai2 = dc2*ai2 + ds2*ar2
        ar2 = ar2h
        DO ik = 1, idl1
          ch2(ik, l) = ch2(ik, l) + ar2*c2(ik, j)
          ch2(ik, lc) = ch2(ik, lc) + ai2*c2(ik, jc)
        END DO
      END DO
    END DO
    DO j = 2, ipph
      DO ik = 1, idl1
        ch2(ik, 1) = ch2(ik, 1) + c2(ik, j)
      END DO
    END DO

    IF (ido < l1) GO TO 132
    DO k = 1, l1
      DO i = 1, ido
        cc(i, 1, k) = ch(i, k, 1)
      END DO
    END DO
    GO TO 135
132 DO i = 1, ido
      DO k = 1, l1
        cc(i, 1, k) = ch(i, k, 1)
      END DO
    END DO
135 DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO k = 1, l1
        cc(ido, j2 - 2, k) = ch(1, k, j)
        cc(1, j2 - 1, k) = ch(1, k, jc)
      END DO
    END DO
    IF (ido == 1) RETURN
    IF (nbd < l1) GO TO 141
    DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO k = 1, l1
        DO i = 3, ido, 2
          ic = idp2 - i
          cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
          cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
          cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
        END DO
      END DO
    END DO
    RETURN
141 DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO i = 3, ido, 2
        ic = idp2 - i
        DO k = 1, l1
          cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
          cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
          cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE radfg

  SUBROUTINE rfftf1(n, c, ch, wa, ifac)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(OUT) :: c(1)
    REAL(DP), INTENT(IN) :: ch(1)
    REAL(DP), INTENT(IN OUT) :: wa(*)
    INTEGER, INTENT(IN) :: ifac(*)
!IMPLICIT REAL(DP)(a-h,o-z)

    nf = ifac(2)
    na = 1
    l2 = n
    iw = n
    DO k1 = 1, nf
      kh = nf - k1
      ip = ifac(kh + 3)
      l1 = l2/ip
      ido = n/l2
      idl1 = ido*l1
      iw = iw - (ip - 1)*ido
      na = 1 - na
      IF (ip /= 4) GO TO 102
      ix2 = iw + ido
      ix3 = ix2 + ido
      IF (na /= 0) GO TO 101
      CALL radf4(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
      GO TO 110
101   CALL radf4(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
      GO TO 110
102   IF (ip /= 2) GO TO 104
      IF (na /= 0) GO TO 103
      CALL radf2(ido, l1, c, ch, wa(iw))
      GO TO 110
103   CALL radf2(ido, l1, ch, c, wa(iw))
      GO TO 110
104   IF (ip /= 3) GO TO 106
      ix2 = iw + ido
      IF (na /= 0) GO TO 105
      CALL radf3(ido, l1, c, ch, wa(iw), wa(ix2))
      GO TO 110
105   CALL radf3(ido, l1, ch, c, wa(iw), wa(ix2))
      GO TO 110
106   IF (ip /= 5) GO TO 108
      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido
      IF (na /= 0) GO TO 107
      CALL radf5(ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      GO TO 110
107   CALL radf5(ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      GO TO 110
108   IF (ido == 1) na = 1 - na
      IF (na /= 0) GO TO 109
      CALL radfg(ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
      na = 1
      GO TO 110
109   CALL radfg(ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
      na = 0
110   l2 = l1
    END DO
    IF (na == 1) RETURN
    DO i = 1, n
      c(i) = ch(i)
    END DO
    RETURN
  END SUBROUTINE rfftf1

  SUBROUTINE sint1(n, war, was, xh, x, ifac)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN OUT) :: war(1)
    REAL(DP), INTENT(IN) :: was(1)
    REAL(DP), INTENT(OUT) :: xh(2)
    REAL(DP), INTENT(IN OUT) :: x(1)
    INTEGER, INTENT(IN OUT) :: ifac(1)
!IMPLICIT REAL(DP)(a-h,o-z)

    DATA sqrt3/1.73205080756888D0/

    DO i = 1, n
      xh(i) = war(i)
      war(i) = x(i)
    END DO
    IF (n - 2 < 0) THEN
      GO TO 101
    ELSE IF (n - 2 == 0) THEN
      GO TO 102
    ELSE
      GO TO 103
    END IF
101 xh(1) = xh(1) + xh(1)
    GO TO 106
102 xhold = sqrt3*(xh(1) + xh(2))
    xh(2) = sqrt3*(xh(1) - xh(2))
    xh(1) = xhold
    GO TO 106
103 np1 = n + 1
    ns2 = n/2
    x(1) = 0D0
    DO k = 1, ns2
      kc = np1 - k
      t1 = xh(k) - xh(kc)
      t2 = was(k)*(xh(k) + xh(kc))
      x(k + 1) = t1 + t2
      x(kc + 1) = t2 - t1
    END DO
    modn = MOD(n, 2)
    IF (modn /= 0) x(ns2 + 2) = 4.*xh(ns2 + 1)
    CALL rfftf1(np1, x, xh, war, ifac)
    xh(1) = .5D0*x(1)
    DO i = 3, n, 2
      xh(i - 1) = -x(i)
      xh(i) = xh(i - 2) + x(i - 1)
    END DO
    IF (modn /= 0) GO TO 106
    xh(n) = -x(n + 1)
106 DO i = 1, n
      x(i) = war(i)
      war(i) = xh(i)
    END DO
    RETURN
  END SUBROUTINE sint1

  FUNCTION sins(x)
    USE params, ONLY: DP

    REAL(DP)::sins
    REAL(DP), INTENT(IN)::x

! stabilised sinh as internal statement FUNCTION

    sins = SINH(SIGN(MIN(ABS(x), 72.0D0), x))

  ENDFUNCTION sins

!INCLUDE "findiff.F90"

END MODULE coulsolv
