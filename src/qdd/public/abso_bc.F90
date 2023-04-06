
! Package containing absorbing boundary conditions and
! subsequent analysis of observables from electron emission.


!-----absbc----------------------------------------------------------

SUBROUTINE absbc(psi, rho)

! Apply absorbing bounds, optionally accumulate absorbed density.
!
! Input/Output:
! psi = set of s.p. wavefunctions
! rho = local densities (spin-up, spin-down)
! Other I/O is handled through MODULE 'params'.

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)

  REAL(DP), DIMENSION(:), ALLOCATABLE :: w3

  LOGICAL :: firstcall = .true.
  INTEGER :: ind, nbe
!------------------------------------------------------------

  IF (nabsorb <= 0) RETURN

! OPTIONAL initialization of mask

  IF (firstcall) THEN
    IF (ispherabso == 1) THEN
      CALL init_spherabso()
    ELSE IF (ispherabso == 0) THEN
      CALL init_abso()
    ELSE IF (ispherabso == 2) THEN
      CALL init_ellipsabso()
    ELSE
      STOP ' this TYPE of absoring boundaries not yet implemented'
    END IF
    firstcall = .false.
  END IF

! SAVE old density

  IF (myn == 0) THEN
    ALLOCATE (w3(kdfull2))
    w3 = rho(1:kdfull2)
  END IF

! apply mask FUNCTION (and accumulate absorption per state)


  DO nbe = 1, nstate
    DO ind = 1, kdfull2
      IF (tgridabso(ind)) THEN
        psi(ind, nbe) = sphermask(ind)*psi(ind, nbe)
      END IF
    END DO
  END DO

! accumulate absorbed density

  CALL calcrho(rho, psi)
  IF (myn == 0) THEN
    DO ind = 1, kdfull2
      rhoabso(ind) = rhoabso(ind) + w3(ind) - rho(ind)
    END DO
    DEALLOCATE (w3)
  END IF

  RETURN
END SUBROUTINE absbc

!-----init_abs_accum----------------------------------------------------------

SUBROUTINE init_abs_accum()

! Initializes accumulator for absorbed densities with zero.

  USE params
  IMPLICIT NONE

  INTEGER :: ind
!--------------------------------------------------------------------

  IF (nabsorb <= 0) RETURN

  DO ind = 1, kdfull2
    rhoabso(ind) = 0D0
  END DO


  RETURN
END SUBROUTINE init_abs_accum

!-----init_absbc-----------------------------------------------------

SUBROUTINE init_absbc(rho)

! Initializes geometry parameters for absorbing bounds.
!
! Input:
! rho = local density

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)

  INTEGER :: i, ind, iphi, itheta, ix, iy, iz, jj
  REAL(DP) :: pp, rmin, tt, x1, y1, z1, xn, yn, zn

  INTEGER, EXTERNAL :: getnearestgridpoint

  IF (iangabso == 1) THEN ! origin is at center of density

    codx = 0D0
    cody = 0D0
    codz = 0D0

    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx

          ind = ind + 1
          codx = codx + rho(ind)*x1
          cody = cody + rho(ind)*y1
          codz = codz + rho(ind)*z1

        END DO
      END DO
    END DO

    codx = codx*dvol
    cody = cody*dvol
    codz = codz*dvol

    xango = codx
    yango = cody
    zango = codz

  ELSE IF (iangabso == 2) THEN ! origin is at center of cluster mass

    xango = 0D0
    yango = 0D0
    zango = 0D0

    DO i = 1, nion
      xango = xango + cx(i)
      yango = yango + cy(i)
      zango = zango + cz(i)
    END DO

    xango = xango/nion
    yango = yango/nion
    zango = zango/nion

  ELSE ! the origin is at the center of the box

    xango = 0D0
    yango = 0D0
    zango = 0D0

  END IF

  IF (iangabso /= 0 .AND. jmp > 0) THEN
! prepare arrays for photo-electron-spectra
! get distance from polar origin TO edge of box

    rmin = 1D10

    IF (ABS(xango - (minx - nxsh)*dx) - dx*(nabsorb + 1) < rmin) &
      rmin = ABS(xango - (minx - nxsh)*dx)
    IF (ABS(xango - (maxx - nxsh)*dx - dx*(nabsorb + 1)) < rmin) &
      rmin = ABS(xango - (maxx - nxsh)*dx)
    IF (ABS(yango - (miny - nysh)*dy - dy*(nabsorb + 1)) < rmin) &
      rmin = ABS(yango - (miny - nysh)*dy)
    IF (ABS(yango - (maxy - nysh)*dy - dy*(nabsorb + 1)) < rmin) &
      rmin = ABS(yango - (maxy - nysh)*dy)
    IF (ABS(zango - (minz - nzsh)*dz - dz*(nabsorb + 1)) < rmin) &
      rmin = ABS(zango - (minz - nzsh)*dz)
    IF (ABS(zango - (maxz - nzsh)*dz - dz*(nabsorb + 1)) < rmin) &
      rmin = ABS(zango - (maxz - nzsh)*dz)

! now make sphere of gridpoints with radius rmin
! these gridpoints shall serve as measure points for the
! outgoing waves

    jj = 0

    DO itheta = 1, nangtheta

      tt = (angthetah - angthetal)/(nangtheta - 1)*(itheta - 1) + angthetal

      DO iphi = 1, nangphi

        jj = jj + 1

        pp = (angphih - angphil)/(nangphi - 1)*(iphi - 1) + angphil

        xn = rmin*COS(pp)*SIN(tt)
        yn = rmin*SIN(pp)*SIN(tt)
        zn = rmin*COS(tt)

        indicesmp(jj) = getnearestgridpoint(xn, yn, zn)

      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE init_absbc

!------------------------------------------------------------

SUBROUTINE init_spherabso()

! Initializes mask FUNCTION for spherical boundary conditions.

  USE params
  USE util, ONLY: prifld
  IMPLICIT NONE

  INTEGER :: i, ind, ix, iy, iz
  REAL(DP) :: cosact, dist, dist2, dmin1, dmin2, dmin12, dmin22
  REAL(DP) :: rx, ry, rz, x1, y1, z1

!------------------------------------------------------------

  WRITE (6, *) 'x,y,zango', xango, yango, zango

  bcrad = nabsorb*dx
  dmin2 = 1D10

  IF (ABS((maxz - nzsh)*dz - zango) < dmin2) &
    dmin2 = ABS((maxz - nzsh)*dz - zango)
  IF (ABS((maxy - nysh)*dy - yango) < dmin2) &
    dmin2 = ABS((maxy - nysh)*dy - yango)
  IF (ABS((maxx - nxsh)*dx - xango) < dmin2) &
    dmin2 = ABS((maxx - nxsh)*dx - xango)

  WRITE (6, *) 'Setting spherical absorbing mask...'
  WRITE (6, *) 'Distance TO edge of box: ', dmin2

  DMIN1 = dmin2 - bcrad

  IF (DMIN1 < 0D0) STOP 'Error IN abso: dmin1<0'

  dmin22 = dmin2**2
  dmin12 = DMIN1**2

  ind = 0
  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = z1 - zango
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = y1 - yango
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = x1 - xango
        dist2 = rx*rx + ry*ry + rz*rz

        ind = ind + 1
        IF (dist2 <= dmin12) THEN
          sphermask(ind) = 1.0D0
          tgridabso(ind) = .false.
        ELSE IF (dist2 > dmin22) THEN
          sphermask(ind) = 0D0
          tgridabso(ind) = .true.
        ELSE
          dist = MAX(small, SQRT(dist2))
          cosact = COS((dist - DMIN1)*0.5D0*pi/bcrad)
          IF (cosact > 0D0) THEN
            sphermask(ind) = cosact**powabso
          ELSE
            sphermask(ind) = 0D0
          END IF
          tgridabso(ind) = .true.
        END IF
        spherloss(ind) = 1D0 - sphermask(ind)**2
      END DO
    END DO
  END DO

  CALL prifld(sphermask, ' mask ')

  DO i = 1, nmps
    WRITE (*, *) ' pt.,mask=', imps(i), sphermask(imps(i)), tgridabso(imps(i))
  END DO

  RETURN
END SUBROUTINE init_spherabso

!-----abso------------------------------------------------------------

SUBROUTINE init_abso()

! Initializes mask for rectangular absorbing boundaries conditions.

  USE params
  IMPLICIT NONE

  INTEGER :: ind, ix, iy, iz

  REAL(DP) :: xmask(nx2), ymask(ny2), zmask(nz2)

  LOGICAL, PARAMETER :: wflag = .true.

!-------------------------------------------------------------------------

! prepare mask functions IN each direction separately

  bcrad = nabsorb*dx

  DO iz = 1, nz2
    zmask(iz) = 1.0D0
  END DO
  DO iz = 1, nabsorb
    zmask(iz) = COS(pi*0.5D0*(nabsorb + 1.0D0 - iz)/nabsorb)**powabso
    zmask(nz2 + 1 - iz) = zmask(iz)
  END DO
  DO iy = 1, ny2
    ymask(iy) = 1.0D0
  END DO
  DO iy = 1, nabsorb
    ymask(iy) = COS(pi*0.5D0*(nabsorb + 1.0D0 - iy)/nabsorb)**powabso
    ymask(ny2 + 1 - iy) = ymask(iy)
  END DO
  DO ix = 1, nx2
    xmask(ix) = 1.0D0
  END DO
  DO ix = 1, nabsorb
    xmask(ix) = COS(pi*0.5D0*(nabsorb + 1.0D0 - ix)/nabsorb)**powabso
    xmask(nx2 + 1 - ix) = xmask(ix)
  END DO
  IF (wflag) THEN
    WRITE (6, '(a)') ' ZMASK:'
    WRITE (6, '(1x,5(1pg12.4))') zmask
    WRITE (6, '(a)') ' YMASK:'
    WRITE (6, '(1x,5(1pg12.4))') ymask
    WRITE (6, '(a)') ' XMASK:'
    WRITE (6, '(1x,5(1pg12.4))') xmask
  END IF

! compose TO one mask FUNCTION on all grid points

  ind = 0
  DO iz = minz, maxz
    DO iy = miny, maxy
      DO ix = minx, maxx
        ind = ind + 1
        sphermask(ind) = xmask(ix)*ymask(iy)*zmask(iz)
        tgridabso(ind) = sphermask(ind) < 0.999999999999D0
        spherloss(ind) = 1D0 - sphermask(ind)**2
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE init_abso

!------------------------------------------------------------

SUBROUTINE init_ellipsabso

! Initializes mask FUNCTION for ellipsoidal boundary conditions.

  USE params
  USE util, ONLY: printfield
  IMPLICIT NONE

  INTEGER :: i, ind, ix, iy, iz
  REAL(DP) :: cosact, ellips1, ellips2, fac
  REAL(DP) :: dmin1x, dmin1y, dmin1z, dmin2, dmin2x, dmin2y, dmin2z
  REAL(DP) :: rx, ry, rz, x1, y1, z1
!------------------------------------------------------------

  WRITE (6, *) 'x,y,zango', xango, yango, zango

  bcrad = nabsorb*dx
  dmin2 = 1D10

  dmin2x = ABS((maxx - nxsh)*dx - xango)
  dmin2y = ABS((maxy - nysh)*dy - yango)
  dmin2z = ABS((maxz - nzsh)*dz - zango)

  dmin2 = MIN(dmin2x, dmin2y, dmin2z)

  WRITE (6, *) 'Setting ellipsoidal absorbing mask...'
  WRITE (6, *) 'Minimum distance TO edge of box: ', dmin2

  dmin1x = dmin2x - dmin2x/dmin2*bcrad != dmin2x*(1.0D0-bcrad/dmin2)
  dmin1y = dmin2y - dmin2y/dmin2*bcrad
  dmin1z = dmin2z - dmin2z/dmin2*bcrad

  IF (dmin1x < 0D0 .OR. dmin1y < 0.0D0 .OR. dmin1z < 0.0D0) &
    STOP 'Error IN abso: dmin1<0'

  WRITE (6, *) dmin1x, dmin1y, dmin1z
  WRITE (6, *) dmin2x, dmin2y, dmin2z
  WRITE (6, *) dmin2

  ind = 0
  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = z1 - zango
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = y1 - yango
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = x1 - xango
        ellips1 = rx*rx/dmin1x/dmin1x + ry*ry/dmin1y/dmin1y &
                 + rz*rz/dmin1z/dmin1z
        ellips2 = rx*rx/dmin2x/dmin2x + ry*ry/dmin2y/dmin2y &
                 + rz*rz/dmin2z/dmin2z

        ind = ind + 1
        IF (ellips1 - 1.0D0 <= 0.0D0) THEN
          sphermask(ind) = 1.0D0
          tgridabso(ind) = .false.
        ELSE IF (ellips2 - 1.0D0 > 0.0D0) THEN
          sphermask(ind) = 0.0D0
          tgridabso(ind) = .false.
        ELSE
          fac = dsqrt(rx*rx/dmin2x/dmin2x + ry*ry/dmin2y/dmin2y &
               + rz*rz/dmin2z/dmin2z)
          fac = (fac - 1.0D0)/bcrad*dmin2 + 1.0D0
          cosact = COS(fac*0.5D0*pi)
          IF (cosact > 0D0) THEN
            sphermask(ind) = cosact**powabso
          ELSE
            sphermask(ind) = 0D0
          END IF
          tgridabso(ind) = .true.
        END IF
        spherloss(ind) = 1.0D0 - sphermask(ind)**2
      END DO
    END DO
  END DO

  OPEN (666, STATUS='unknown', FILE='pabsomask.'//outnam)
  CALL printfield(666, spherloss, 'mask')
  CLOSE (666)

  DO i = 1, nmps
    WRITE (*, *) ' pt.,mask=', imps(i), sphermask(imps(i)), tgridabso(imps(i))
  END DO

  RETURN
END SUBROUTINE init_ellipsabso

!------------------------------------------------------------

SUBROUTINE initmeasurepoints

! Initializes measuring points for collecting PES information.

  USE params
  IMPLICIT NONE
  INTEGER :: ik, ith, iph, impsact, impsx, impsy, impsz
  REAL(DP):: dmin2, p, r, t, x, y, z
  INTEGER, EXTERNAL :: getnearestgridpoint

  dmin2 = 1D10

  IF (ABS((maxz - nzsh)*dz - zango) < dmin2) &
    dmin2 = ABS((maxz - nzsh)*dz - zango)
  IF (ABS((maxy - nysh)*dy - yango) < dmin2) &
    dmin2 = ABS((maxy - nysh)*dy - yango)
  IF (ABS((maxx - nxsh)*dx - xango) < dmin2) &
    dmin2 = ABS((maxx - nxsh)*dx - xango)

  bcrad = nabsorb*dx
  r = dmin2 - bcrad
  WRITE (*, *) ' analyzing point at r,ir=', r, nint(r/dx)

  nmps = 0

  DO ith = 1, nmptheta + 1
    DO iph = 1, nmpphi

      p = (iph - 1)*2D0*pi/nmpphi
      t = (ith - 1)*pi/nmptheta

      nmps = nmps + 1

      IF (nmps > maxmps) STOP ' array for analyzing points exhausted'

      x = r*COS(p)*SIN(t)
      y = r*SIN(p)*SIN(t)
      z = r*COS(t)

      imps(nmps) = getnearestgridpoint(x, y, z)
      DO ik = 1, nmps - 1
        IF (imps(nmps) == imps(ik)) THEN
          nmps = nmps - 1
        END IF
      END DO

    END DO
  END DO

  DO ik = 1, nmps
    impsact = imps(ik)
    impsx = mod(impsact - 1, nx2) + 1
    impsy = mod(impsact/nx2, ny2) + 1
    impsz = mod(impsact/(nx2*ny2), nz2) + 1
    WRITE (*, *) ' analyzing pts: nmps,imps=', nmps, impsact
    WRITE (*, *) ' nx,ny,nz=', impsx, impsy, impsz
    WRITE (*, '(1a,4f12.5)') ' x,y,z,r=', x, y, z, dsqrt(x*x + y*y + z*z)
    WRITE (*, *) ' theta,phi=', t, p
  END DO

  RETURN
END SUBROUTINE initmeasurepoints

SUBROUTINE angular_distribution()

! PAD by trapezoidal interpolation from electron loss
! collected in 'rhoabso'.

  USE params
  IMPLICIT NONE

  REAL(8), ALLOCATABLE :: angdistp(:),angdist(:,:)
  REAL(8), ALLOCATABLE ::  wt(:), theta(:)
!  REAL(8), PARAMETER :: pi = 3.141592653589793D0
  INTEGER :: nth, nph, nrh, ind, nmaxrho, ix, iy, iz
  REAL(8) :: drho, phiang, rvecx, rvecy, rvecz, x1, y1, z1
  REAL(8) :: rxrel, ryrel, rzrel, rxdiff, rydiff, rzdiff
  REAL(8) :: total, wg, func, dmax2, deltatheta, rhoact

  ALLOCATE(angdist(nangtheta, nangphi))
  ALLOCATE(angdistp(nangtheta),wt(nangtheta),theta(nangtheta))

  angdist = 0D0
  angdistp = 0D0
  drho = dx/8D0 ! empirical choice
  dmax2 = MAX(kxbox, kybox, kzbox)/2D0*dx
  nmaxrho = dmax2/rhoact
  deltatheta = pi/(nangtheta - 1)

  WRITE (6, '(1a,1i4)') ' nmaxrho = ', nmaxrho
  WRITE (6, '(1a,1i4)') ' nangtheta = ', nangtheta
  WRITE (6, '(1a,1i4)') ' nangphi = ', nangphi
  WRITE (6, '(1a,1i4)') ' nabsorb = ', nabsorb
  WRITE (6, '(1a,1f8.4)') ' dx = ', dx
  WRITE (6, '(1a,1f8.4)') ' dmax2 = ', dmax2
  WRITE (6, '(1a,1f8.4)') ' drho = ', drho

! initialize theta and phiang
  DO nth = 1, nangtheta
    theta(nth) = deltatheta*(nth - 1)
    IF (nth == 1 .OR. nth == nangtheta) THEN
      wt(nth) = 2D0*pi*(1D0 - cos(deltatheta/2D0))
    ELSE
      wt(nth) = 2D0*pi*(dcos(theta(nth) - deltatheta/2.) &
                        - dcos(theta(nth) + deltatheta/2.))
    END IF
    wt(nth) = dabs(wt(nth))
  END DO

WRITE(*,*) 'WT:',wt

! calculate angdist
  DO nth = 1, nangtheta
    DO nph = 1, nangphi
      DO nrh = 1, nmaxrho
        phiang = 2D0*pi*(nph - 1)/nangphi
        rhoact = drho*nrh

        rvecx = rhoact*dcos(phiang)*dsin(theta(nth))
        rvecy = rhoact*dsin(phiang)*dsin(theta(nth))
        rvecz = rhoact*dcos(theta(nth))

        IF (rvecx > kxbox/2D0*dx .OR. rvecx < -(kxbox - 2)/2D0*dx) goto 23 ! OUT of box
        IF (rvecy > kybox/2D0*dx .OR. rvecy < -(kybox - 2)/2D0*dx) goto 23 ! OUT of box
        IF (rvecz > kzbox/2D0*dx .OR. rvecz < -(kzbox - 2)/2D0*dx) goto 23 ! OUT of box
        IF (icheckabsozone(rvecx, rvecy, rvecz) == 0) goto 23 ! OUT of abso zone

        ind = 0
        DO iz = 1, kzbox
          DO iy = 1, kybox
            DO ix = 1, kxbox
              z1 = (iz - kzbox/2)*dx
              y1 = (iy - kybox/2)*dx
              x1 = (ix - kxbox/2)*dx

              ind = ind + 1

              rxrel = x1!-xo
              ryrel = y1!-xo
              rzrel = z1!-xo

              rxdiff = rvecx - rxrel
              rydiff = rvecy - ryrel
              rzdiff = rvecz - rzrel

              func = gtent(rxdiff, rydiff, rzdiff) ! tent FUNCTION
              angdist(nth, nph) = angdist(nth, nph) &
                     + rhoabso(ind)*rhoact**2D0*func*drho
            END DO
          END DO
        END DO

23      CONTINUE
      END DO ! nmaxrho loop
    END DO ! nangphi loop
  END DO ! nangtheta loop

! WRITE OUT
  DO nth = 1, nangtheta
    DO nph = 1, nangphi
      IF (nth == 1) angdist(nth, nph) = angdist(1, 1)
      IF (nth == nangtheta) angdist(nth, nph) = angdist(nangtheta, 1)
      WRITE (47, '(2f17.4,1e17.7)') (nph - 1)*360D0/nangphi, &
         theta(nth)/pi*180D0, angdist(nth, nph)
    END DO
    WRITE (47, '(2f17.4,1e17.7)') 360D0, theta(nth)/pi*180D0, angdist(nth, 1)
    WRITE (47, *)
  END DO

! calculate phiang average and PRINT OUT pangdist3-2
  total = 0.0
  wg = 0.0
  DO nth = 1, nangtheta
    DO nph = 1, nangphi
      angdistp(nth) = angdistp(nth) + angdist(nth, nph)/nangphi
    END DO
    total = total + angdistp(nth)*wt(nth)
    wg = wg + wt(nth)
  END DO

  DO nth = 1, nangtheta
    WRITE (29, '(1f17.4,2e17.7)') theta(nth)/pi*180D0, &
      angdistp(nth)/(total/4D0/pi), angdistp(nth)/4D0/pi
  END DO

! finish
  WRITE (*, '(1a,1e17.7)') ' total = ', total
  WRITE (*, '(1a,1f17.14)') ' wg = ', wg

  DEALLOCATE(angdist)
  DEALLOCATE(angdistp)


  RETURN

CONTAINS

! cccccccccccccccccccccccccccccccccccccccccccccccccc

REAL(8) FUNCTION gbox(xx, yy, zz)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: xx, yy, zz

  gbox = 1D0
  IF (abs(xx) .gt. dx/2D0) gbox = 0D0
  IF (abs(yy) .gt. dx/2D0) gbox = 0D0
  IF (abs(zz) .gt. dx/2D0) gbox = 0D0
END FUNCTION gbox

REAL(8) FUNCTION gtent(xx, yy, zz)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: xx, yy, zz

  gtent = max(dx - dabs(xx), 0D0)*max(dx - dabs(yy), 0D0)&
         *max(dx - dabs(zz), 0D0)
  gtent = gtent/dx**3
END FUNCTION gtent

! cccccccccccccccccccccccccccccccccccccccccccccccccc

INTEGER FUNCTION icheckabsozone(xx, yy, zz)
  IMPLICIT REAL(KIND(1D0)) (A - H, O - Z)
  REAL(8), INTENT(IN) :: xx, yy, zz

  REAL(8) :: dmin2, dmin1, dist2, dmin12, dmin22, ellips1, ellips2
  REAL(8) :: dmin2x, dmin2y, dmin2z, dmin1x, dmin1y, dmin1z
  REAL(8) :: dmin22x, dmin22y, dmin22z, dmin12x, dmin12y, dmin12z

  icheckabsozone = 0

  IF (ispherabso == 0) THEN
    IF (xx < (-(kxbox - 2)/2D0 + nabsorb)*dx .OR. xx > (kxbox/2D0 - nabsorb)*dx) icheckabsozone = 1
    IF (yy < (-(kybox - 2)/2D0 + nabsorb)*dx .OR. yy > (kybox/2D0 - nabsorb)*dx) icheckabsozone = 1
    IF (zz < (-(kzbox - 2)/2D0 + nabsorb)*dx .OR. zz > (kzbox/2D0 - nabsorb)*dx) icheckabsozone = 1
  ELSE IF (ispherabso == 1) THEN
    dmin2 = MIN(kxbox, kybox, kzbox)/2D0*dx
    dmin1 = dmin2 - nabsorb*dx
    dist2 = xx*xx + yy*yy + zz*zz
    dmin12 = dmin1*dmin1
    dmin22 = dmin2*dmin2
    IF (dist2 > dmin12 .AND. dist2 < dmin22) icheckabsozone = 1
  ELSE IF (ispherabso == 2) THEN
    dmin2 = MIN(kxbox, kybox, kzbox)/2D0*dx
    dmin2x = kxbox/2D0*dx
    dmin2y = kybox/2D0*dx
    dmin2z = kzbox/2D0*dx
    dmin1x = dmin2x - dmin2x/dmin2*nabsorb*dx
    dmin1y = dmin2y - dmin2y/dmin2*nabsorb*dx
    dmin1z = dmin2z - dmin2z/dmin2*nabsorb*dx
    dmin12x = dmin1x*dmin1x
    dmin12y = dmin1y*dmin1y
    dmin12z = dmin1z*dmin1z
    dmin22x = dmin2x*dmin2x
    dmin22y = dmin2y*dmin2y
    dmin22z = dmin2z*dmin2z
    ellips1 = xx*xx/dmin12x + yy*yy/dmin12y + zz*zz/dmin12z
    ellips2 = xx*xx/dmin22x + yy*yy/dmin22y + zz*zz/dmin22z
    IF (ellips1 > 1D0 .AND. ellips2 <= 1D0) icheckabsozone = 1
  END IF

  RETURN
END FUNCTION icheckabsozone

END SUBROUTINE angular_distribution


!-----nesacpe------------------------------------------------

SUBROUTINE nescape(it, rho)

! Compute and PRINT total NUMBER of escaped electrons (on the fly).
!
! Input:
! rho = electron density
! it = nr. of time step IN calling routine

  USE params
  USE util, ONLY: safeopen
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: it
  REAL(DP), INTENT(IN) :: rho(kdfull2)

  REAL(DP) :: absosum, tinfs
!------------------------------------------------------------

  tinfs = it*dt1*0.0484D0/2D0/ame
  apnum = dvol*SUM(rho)
  absosum = dvol*SUM(rhoabso)

  CALL safeopen(23, it, jesc, 'pescel')
! OPEN(23,POSITION='append',FILE='pescel.'//outnam)
  WRITE (23, '(f13.5,3(1pg13.5))') tinfs, 1D0 - apnum/nelect, nelect - apnum, absosum
  FLUSH (23)

  RETURN
END SUBROUTINE nescape

!------------------------------------------------------------

SUBROUTINE evalmp(iunit, q0)

! Evaluate and store information for later analysis of PES.
!
! Input:
! iunit = wanted output UNIT
! q0 = set of s.p. wavefunctions

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iunit
  COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)

  LOGICAL :: topenf
  INTEGER :: i, ii, j, nbe
  REAL(DP) :: scal1, scal2, x1, y1, z1 !,fpulseinteg1,fpulseinteg2
  REAL(DP) :: rp(2*maxmps)
  COMPLEX(DP) :: q0phase


  IF (myn == 0) THEN
    INQUIRE (803, OPENED=topenf)
    IF (.NOT. topenf) &
      OPEN (803, POSITION='append', FILE='pMP.'//outnam)
  END IF


  DO nbe = 1, nstate
    j = 0
    DO i = 1, nmps
      ii = imps(i)
      CALL conv1to3(ii)
      x1 = (iindtmp(1) - nxsh)*dx
      y1 = (iindtmp(2) - nysh)*dy
      z1 = (iindtmp(3) - nzsh)*dz
      scal1 = x1*e1x + y1*e1y + z1*e1z
      scal2 = x1*e2x + y1*e2y + z1*e2z
      q0phase = q0(ii, nbe)* &
                EXP(CMPLX(0D0, -e0*(scal1*fpulseinteg1 + scal2*fpulseinteg2), DP))
      j = j + 1
      rp(j) = REAL(q0phase, DP)
      j = j + 1
      rp(j) = AIMAG(q0phase)
    END DO
    ! write phase corrections including field strengths
    WRITE (iunit, '(1f14.5,1000e18.8)') tfs, rp(1:2*nmps), &
       fpulse2integ1*e0,fpulse2integ1*e0_2,&
       f2pulse2integ1*e0**2,f2pulse2integ1*e0_2**2
  END DO


  IF (myn == 0) CLOSE (iunit)

  RETURN
END SUBROUTINE evalmp
!------------------------------------------------------------

