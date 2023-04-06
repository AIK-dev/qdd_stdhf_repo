!CALL-----COUNET------------------------------------------------------
!
MODULE coulsolv
  USE params
  USE fishpack
  IMPLICIT NONE
CONTAINS
  SUBROUTINE solv_poisson(RHO, POTC, KDF)
!
! SOLVES THE POISSON EQUATION IN THE 3 DIMENSIONAL DIRECT SPACE
! USING THE FISHPACK HW3CRT ROUTINE
! INPUT :RHO NETTO CHARGE DENSITY
!
! OUTPUT:POTC COULOMB FIELD
!
    USE params
    IMPLICIT NONE
!---------------------------------------------------------------------
    REAL(DP), INTENT(IN) :: RHO(KDF)
    REAL(DP), INTENT(OUT) :: POTC(KDF)
    INTEGER, INTENT(IN) :: KDF
    REAL(DP), ALLOCATABLE :: W(:)
    REAL(DP), ALLOCATABLE :: POTC_MAT(:, :, :)
!---------------------------------------------------------------------

    INTEGER :: I1, I2, I3, IERROR, IND

    REAL(DP) :: FTPI, PER
    REAL(DP):: potbcxl(miny:maxy, minz:maxz)
    REAL(DP):: potbcxu(miny:maxy, minz:maxz)
    REAL(DP):: potbcyl(minx:maxx, minz:maxz)
    REAL(DP):: potbcyu(minx:maxx, minz:maxz)
    REAL(DP):: potbczl(minx:maxx, miny:maxy)
    REAL(DP):: potbczu(minx:maxx, miny:maxy)

    ALLOCATE (POTC_MAT(kxbox, kybox, kzbox))
    ALLOCATE (W(30 + kxbox + kybox + 5*kzbox + MAX(kxbox, kybox, kzbox) + 7*(INT((kxbox + 1)/2) + INT((kybox + 1)/2))))

!
! CALCULATION OF THE BOUNDARIES CONDITIONS
!

    CALL boundc(RHO, dx, dy, dz, potbcxl, potbcxu, potbcyl, potbcyu, potbczl, potbczu, POTC_MAT, 0)

!
! INCLUSION OF THE SOURCE TERM
!

    FTPI = -4.*3.1415927

    DO I3 = 2, (kzbox - 1)
      DO I2 = 2, (kybox - 1)
        DO I1 = 2, (kxbox - 1)

          IND = (I1) + (I2 - 1)*(kxbox) + (I3 - 1)*(kybox)*(kxbox)
          POTC_MAT(I1, I2, I3) = FTPI*RHO(IND)

        END DO
      END DO
    END DO

!WRITE(6,*) minx,maxx,kxbox
!STOP'test solvpoisson'

    CALL HW3CRT((minx - 1)*dx, (maxx - 1)*dx, kxbox - 1, 1, potbcxl, potbcxu, &
                (miny - 1)*dy, (maxy - 1)*dy, kybox - 1, 1, potbcyl, potbcyu, &
                (minz - 1)*dz, (maxz - 1)*dz, kzbox - 1, 1, potbczl, potbczu, &
                0.0D0, kxbox, kybox, POTC_MAT, PER, IERROR, W)

    IND = 0
    DO I3 = 1, kzbox
      DO I2 = 1, kybox
        DO I1 = 1, kxbox
          !IND= I1+(I2-1)*(kxbox)+(I3-1)*(kybox)*(kxbox)
          IND = IND + 1
          POTC(IND) = POTC_MAT(I1, I2, I3)
        END DO
      END DO
    END DO

    POTC = POTC*e2

    DEALLOCATE (POTC_MAT, W)
    RETURN
  END SUBROUTINE solv_poisson

!-----boundc-----------------------------------------------------------

  SUBROUTINE boundc(rhoin, deltax, deltay, deltaz, &
                    potbcxl, potbcxu, potbcyl, potbcyu, potbczl, potbczu, &
                    potini, ifinit)

    USE params, ONLY: DP, maxx, maxy, maxz, minx, miny, minz, one, zero, pi
    IMPLICIT NONE
! For given density 'rhoin'
! compute multipole moments 'mulmLM'
! and the Coulomb potential at the six boundary planes of the cube.
! The mesh spacing is entered through 'deltax,deltay,deltaz'
! and the mesh SIZE is prescribed IN 'all.inc'.

! Additionally, an initial guess for the potential is
! computed and returned on 'potini'.

! Note: This code is optimized for readability not yet for speed.

!INCLUDE 'work3D.inc'

    REAL(DP), INTENT(IN) :: rhoin(minx:maxx, miny:maxy, minz:maxz)
    REAL(DP), INTENT(IN) :: deltax
    REAL(DP), INTENT(IN) :: deltay
    REAL(DP), INTENT(IN) :: deltaz
    REAL(DP), INTENT(OUT) :: potbcxl(miny:maxy, minz:maxz)
    REAL(DP), INTENT(OUT) :: potbcxu(miny:maxy, minz:maxz)
    REAL(DP), INTENT(OUT) :: potbcyl(minx:maxx, minz:maxz)
    REAL(DP), INTENT(OUT) :: potbcyu(minx:maxx, minz:maxz)
    REAL(DP), INTENT(OUT) :: potbczl(minx:maxx, miny:maxy)
    REAL(DP), INTENT(OUT) :: potbczu(minx:maxx, miny:maxy)
    REAL(DP), INTENT(OUT) :: potini(minx:maxx, miny:maxy, minz:maxz)
    INTEGER, INTENT(IN) :: ifinit

! local variables

    LOGICAL :: testpr
    REAL(DP):: mulm00, mulm10, mulm11r, mulm11i, mulm20, mulm21r, mulm21i, mulm22r, mulm22i

    INTEGER :: ix, iy, iz
    REAL(DP):: fac, x, y, z, x2, y2, z2, dr2, r, r2
    REAL(DP):: xmin, xmax, ymin, ymax, zmin, zmax
    INTEGER:: mx, my, mz
!REAL(DP),DIMENSION(2):: xx2,yy2,zz2,xx,yy,zz,drr2,rr2,rr
!!$REAL(DP):: y00(minx:maxx,miny:maxy,minz:maxz),&
!!$ y10(minx:maxx,miny:maxy,minz:maxz),&
!!$ y11r(minx:maxx,miny:maxy,minz:maxz),&
!!$ y11i(minx:maxx,miny:maxy,minz:maxz),&
!!$ y20(minx:maxx,miny:maxy,minz:maxz),&
!!$ y21r(minx:maxx,miny:maxy,minz:maxz),&
!!$ y21i(minx:maxx,miny:maxy,minz:maxz),&
!!$ y22r(minx:maxx,miny:maxy,minz:maxz),&
!!$ y22i(minx:maxx,miny:maxy,minz:maxz)

! damping radius for initialization

    REAL(DP), PARAMETER :: rms2ini = 2.0D0

! forefactors of the spherical harmonics

!!$REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
!!$REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
!!$REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
!!$REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
!!$REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
!!$REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

!DATA testpr/.false./

!INTEGER :: maxmx,maxmy,maxmz

!maxmx = maxx-1
!maxmy = maxy-1
!maxmz = maxz-1

!maxmx = maxx
!maxmy = maxy
!maxmz = maxz

! the (REAL) spherical harmonics as internal statement functions

!y00(x,y,z) = pfy00
!y10(x,y,z) = pfy10*z
!y11r(x,y,z) = -pfy11*x
!y11i(x,y,z) = -pfy11*y
!y20(x,y,z) = pfy20*(2.0D0*(z*z)-(x*x)-(y*y))
!y21r(x,y,z) = -pfy21*(x*z)
!y21i(x,y,z) = -pfy21*(y*z)
!y22r(x,y,z) = pfy22*((x*x)-(y*y))
!y22i(x,y,z) = pfy22*2.0D0*(x*y)

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

!mz = (maxz-minz)/2
!my = (maxy-miny)/2
!mx = (maxx-minx)/2

    mz = (maxz)/2
    my = (maxy)/2
    mx = (maxx)/2

    DO iz = minz, maxz
      z = deltaz*(iz - mz)
      DO iy = miny, maxy
        y = deltay*(iy - my)
        DO ix = minx, maxx
          x = deltax*(ix - mx)
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

!WRITE(6,*) mx,my,mz
!WRITE(6,*) maxx,maxy,maxz
!WRITE(6,*) minx,miny,minz
!WRITE(6,*) x,y,z
!STOP'middle x,y,z'

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

!!$! z-boundaries
!!$
!!$zz(1) = deltaz*(maxz)
!!$zz(2) = deltaz*(minz)
!!$
!!$zz2 = zz*zz
!!$
!!$DO iy=miny,maxy
!!$ y = deltay*(iy)
!!$ drr2 = y*y+zz2
!!$ DO ix=minx,maxx
!!$ x = deltax*(ix)
!!$ rr2 = one/(x*x+drr2)
!!$ rr = SQRT(rr2)
!!$ potbczu(ix,iy) = rr(1)*(mulm00*y00(x,y,zz(1)) +rr2(1)*( mulm10*y10(x,y,zz(1)) &
!!$ +mulm11r*y11r(x,y,zz(1)) +mulm11i*y11i(x,y,zz(1)) &
!!$ +rr2(1)*( mulm20*y20(x,y,zz(1)) +mulm21r*y21r(x,y,zz(1)) &
!!$ +mulm21i*y21i(x,y,zz(1)) +mulm22r*y22r(x,y,zz(1)) &
!!$ +mulm22i*y22i(x,y,zz(1)) )))
!!$ potbczl(ix,iy) = rr(2)*(mulm00*y00(x,y,zz(2)) +rr2(2)*( mulm10*y10(x,y,zz(2)) &
!!$ +mulm11r*y11r(x,y,zz(2)) +mulm11i*y11i(x,y,zz(2)) &
!!$ +rr2(2)*( mulm20*y20(x,y,zz(2)) +mulm21r*y21r(x,y,zz(2)) &
!!$ +mulm21i*y21i(x,y,zz(2)) +mulm22r*y22r(x,y,zz(2)) &
!!$ +mulm22i*y22i(x,y,zz(2)) )))
!!$ END DO
!!$END DO
!!$
!!$! y-boundaries
!!$
!!$yy(1) = deltay*(maxy)
!!$yy(2) = deltay*(miny)
!!$
!!$yy2 = yy*yy
!!$
!!$DO iz=minz,maxz
!!$ z = deltaz*(iz)
!!$ drr2 = z*z+yy2
!!$ DO ix=minx,maxx
!!$ x = deltax*(ix)
!!$ rr2 = one/(x*x+drr2)
!!$ rr = SQRT(rr2)
!!$ potbcyu(ix,iz) = rr(1)*(mulm00*y00(x,yy(1),z) +rr2(1)*( mulm10*y10(x,yy(1),z) &
!!$ +mulm11r*y11r(x,yy(1),z) +mulm11i*y11i(x,yy(1),z) &
!!$ +rr2(1)*( mulm20*y20(x,yy(1),z) +mulm21r*y21r(x,yy(1),z) &
!!$ +mulm21i*y21i(x,yy(1),z) +mulm22r*y22r(x,yy(1),z) &
!!$ +mulm22i*y22i(x,yy(1),z) )))
!!$ potbcyl(ix,iz) = rr(2)*(mulm00*y00(x,yy(2),z) +rr2(2)*( mulm10*y10(x,yy(2),z) &
!!$ +mulm11r*y11r(x,yy(2),z) +mulm11i*y11i(x,yy(2),z) &
!!$ +rr2(2)*( mulm20*y20(x,yy(2),z) +mulm21r*y21r(x,yy(2),z) &
!!$ +mulm21i*y21i(x,yy(2),z) +mulm22r*y22r(x,yy(2),z) &
!!$ +mulm22i*y22i(x,yy(2),z) )))
!!$ END DO
!!$END DO
!!$! x-boundaries
!!$
!!$xx(1) = deltax*(maxx)
!!$xx(2) = deltax*(minx)
!!$xx2 = xx*xx
!!$DO iz=minz,maxz
!!$ z = deltaz*(iz-1)
!!$ drr2 = z*z+xx2
!!$ DO iy=miny,maxy
!!$ y = deltay*(iy)
!!$ rr2 = one/(y*y+drr2)
!!$ rr = SQRT(rr2)
!!$ potbcxu(iy,iz) = rr(1)*(mulm00*y00(xx(1),y,z) +rr2(1)*( mulm10*y10(xx(1),y,z) &
!!$ +mulm11r*y11r(xx(1),y,z) +mulm11i*y11i(xx(1),y,z) &
!!$ +rr2(1)*( mulm20*y20(xx(1),y,z) +mulm21r*y21r(xx(1),y,z) &
!!$ +mulm21i*y21i(xx(1),y,z) +mulm22r*y22r(xx(1),y,z) &
!!$ +mulm22i*y22i(xx(1),y,z) )))
!!$ potbcxl(iy,iz) = rr(2)*(mulm00*y00(xx(2),y,z) +rr2(2)*( mulm10*y10(xx(2),y,z) &
!!$ +mulm11r*y11r(xx(2),y,z) +mulm11i*y11i(xx(2),y,z) &
!!$ +rr2(2)*( mulm20*y20(xx(2),y,z) +mulm21r*y21r(xx(2),y,z) &
!!$ +mulm21i*y21i(xx(2),y,z) +mulm22r*y22r(xx(2),y,z) &
!!$ +mulm22i*y22i(xx(2),y,z) )))
!!$ END DO
!!$END DO

! z-boundaries

    z = deltaz*(maxz - mz)

    zmax = deltaz*(maxz - mz)
    zmin = deltaz*(minz - mz)

    z2 = z*z
    DO iy = miny, maxy
      y = deltay*(iy - my)
      dr2 = y*y + z2
      DO ix = minx, maxx
        x = deltax*(ix - mx)
        r2 = one/(x*x + dr2)
        r = SQRT(r2)
        potbczu(ix, iy) = r*(mulm00*y00(x, y, zmax) + r2*(mulm10*y10(x, y, zmax) &
                                                          + mulm11r*y11r(x, y, zmax) + mulm11i*y11i(x, y, zmax) &
                                                          + r2*(mulm20*y20(x, y, zmax) + mulm21r*y21r(x, y, zmax) &
                                                                + mulm21i*y21i(x, y, zmax) + mulm22r*y22r(x, y, zmax) &
                                                                + mulm22i*y22i(x, y, zmax))))
        potbczl(ix, iy) = r*(mulm00*y00(x, y, zmin) + r2*(mulm10*y10(x, y, zmin) &
                                                          + mulm11r*y11r(x, y, zmin) + mulm11i*y11i(x, y, zmin) &
                                                          + r2*(mulm20*y20(x, y, zmin) + mulm21r*y21r(x, y, zmin) &
                                                                + mulm21i*y21i(x, y, zmin) + mulm22r*y22r(x, y, zmin) &
                                                                + mulm22i*y22i(x, y, zmin))))
      END DO
    END DO

! y-boundaries

    y = deltay*(maxy - my)

    ymax = deltay*(maxy - my)
    ymin = deltay*(miny - my)

    y2 = y*y
    DO iz = minz, maxz
      z = deltaz*(iz - mz)
      dr2 = z*z + y2
      DO ix = minx, maxx
        x = deltax*(ix - mx)
        r2 = one/(x*x + dr2)
        r = SQRT(r2)
        potbcyu(ix, iz) = r*(mulm00*y00(x, ymax, z) + r2*(mulm10*y10(x, ymax, z) &
                                                          + mulm11r*y11r(x, ymax, z) + mulm11i*y11i(x, ymax, z) &
                                                          + r2*(mulm20*y20(x, ymax, z) + mulm21r*y21r(x, ymax, z) &
                                                                + mulm21i*y21i(x, ymax, z) + mulm22r*y22r(x, ymax, z) &
                                                                + mulm22i*y22i(x, ymax, z))))
        potbcyl(ix, iz) = r*(mulm00*y00(x, ymin, z) + r2*(mulm10*y10(x, ymin, z) &
                                                          + mulm11r*y11r(x, ymin, z) + mulm11i*y11i(x, ymin, z) &
                                                          + r2*(mulm20*y20(x, ymin, z) + mulm21r*y21r(x, ymin, z) &
                                                                + mulm21i*y21i(x, ymin, z) + mulm22r*y22r(x, ymin, z) &
                                                                + mulm22i*y22i(x, ymin, z))))
      END DO
    END DO

! x-boundaries

    x = deltax*(maxx - mx)

    xmax = deltax*(maxx - mx)
    xmin = deltax*(minx - mx)

    x2 = x*x
    DO iz = minz, maxz
      z = deltaz*(iz - mz)
      dr2 = z*z + x2
      DO iy = miny, maxy
        y = deltay*(iy - my)
        r2 = one/(y*y + dr2)
        r = SQRT(r2)
        potbcxu(iy, iz) = r*(mulm00*y00(xmax, y, z) + r2*(mulm10*y10(xmax, y, z) &
                                                          + mulm11r*y11r(xmax, y, z) + mulm11i*y11i(xmax, y, z) &
                                                          + r2*(mulm20*y20(xmax, y, z) + mulm21r*y21r(xmax, y, z) &
                                                                + mulm21i*y21i(xmax, y, z) + mulm22r*y22r(xmax, y, z) &
                                                                + mulm22i*y22i(xmax, y, z))))
        potbcxl(iy, iz) = r*(mulm00*y00(xmin, y, z) + r2*(mulm10*y10(xmin, y, z) &
                                                          + mulm11r*y11r(xmin, y, z) + mulm11i*y11i(xmin, y, z) &
                                                          + r2*(mulm20*y20(xmin, y, z) + mulm21r*y21r(xmin, y, z) &
                                                                + mulm21i*y21i(xmin, y, z) + mulm22r*y22r(xmin, y, z) &
                                                                + mulm22i*y22i(xmin, y, z))))
      END DO
    END DO

    potini(minx, miny:maxy, minz:maxz) = potbcxl(miny:maxy, minz:maxz)
    potini(maxx, miny:maxy, minz:maxz) = potbcxu(miny:maxy, minz:maxz)
    potini(minx:maxx, miny, minz:maxz) = potbcyl(minx:maxx, minz:maxz)
    potini(minx:maxx, maxy, minz:maxz) = potbcyu(minx:maxx, minz:maxz)
    potini(minx:maxx, miny:maxy, minz) = potbczl(minx:maxx, miny:maxy)
    potini(minx:maxx, miny:maxy, maxz) = potbczu(minx:maxx, miny:maxy)

!DO ix = minx,maxx
!WRITE(6,*) ix,2*potini(ix,miny,minz)
!END DO
!STOP'test boundc'

    IF (ifinit /= 1) RETURN

! initialization IN interior

    DO iz = minz, maxz
      z = deltaz*iz
      dr2 = z*z + rms2ini
      DO iy = miny, maxy
        y = deltay*iy
        DO ix = minx, maxx
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

! the (REAL) spherical harmonics as internal statement functions
  REAL(DP) FUNCTION y00(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y00 = pfy00
    RETURN
  ENDFUNCTION y00
  FUNCTION y10(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y10
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y10 = pfy10*z
  ENDFUNCTION y10
  FUNCTION y11r(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y11r
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y11r = -pfy11*x
  ENDFUNCTION y11r
  FUNCTION y11i(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y11i
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y11i = -pfy11*y
  ENDFUNCTION y11i
  FUNCTION y20(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y20
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y20 = pfy20*(2.0D0*(z*z) - (x*x) - (y*y))
  ENDFUNCTION y20
  FUNCTION y21r(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y21r
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y21r = -pfy21*(x*z)
  ENDFUNCTION y21r
  FUNCTION y21i(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y21i
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y21i = -pfy21*(y*z)
  ENDFUNCTION y21i
  FUNCTION y22r(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y22r
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y22r = pfy22*((x*x) - (y*y))
  ENDFUNCTION y22r
  FUNCTION y22i(x, y, z)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP) :: y22i
    REAL(DP), INTENT(IN) :: x, y, z

    REAL(DP), PARAMETER :: pfy00 = 0.282094791D0
    REAL(DP), PARAMETER :: pfy10 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy11 = 0.488602511D0
    REAL(DP), PARAMETER :: pfy20 = 0.315391565D0
    REAL(DP), PARAMETER :: pfy21 = 1.092548431D0
    REAL(DP), PARAMETER :: pfy22 = 0.546274215D0

    y22i = pfy22*2.0D0*(x*y)
  ENDFUNCTION y22i

END MODULE coulsolv
