PROGRAM angdist5
  IMPLICIT REAL(KIND(1D0)) (A - H, O - Z)
  INTEGER, PARAMETER :: DP = KIND(1D0)
! plots the angular distribution of ionization as density plot;
! requires pescmask FILE
! collects electron loss by using grid point represented by
! box/tent distribution

  INTEGER, PARAMETER :: nstepstheta = 30
  INTEGER, PARAMETER :: nstepsphi = 60

  REAL(8) :: dx, bcrad, dum
  INTEGER :: kxbox, kybox, kzbox, nabsorb, ispherabso


  REAL(8), ALLOCATABLE :: rhoabso(:)
  REAL(8), PARAMETER :: pi = 3.141592653589793D0

! ---------------- Settings ----------------------------------------------------
  dx = 0.4120935 !0.588705D0 !dx = 0.824187016040328D0
  kxbox = 80
  kybox = 80
  kzbox = 80
  ispherabso = 1 ! (0) cartesian (1) spherical (2) ellipsoidal mask
  nabsorb = 6
!xo = 0D0
!yo = 0D0
!zo = 0D0

  kdfull2 = kxbox*kybox*kzbox
  bcrad = nabsorb*dx
  ALLOCATE (rhoabso(kdfull2))


! READ IN mask
  DO ind = 1, kdfull2
    READ (5, *) dum, dum, dum, rhoabso(ind)
  END DO


  OPEN (28, STATUS='unknown', FILE='2dPAD.res')
  OPEN (29, STATUS='unknown', FILE='1dPAD.res')
  CALL angular_distribution()

STOP "regular end"

CONTAINS

SUBROUTINE angular_distribution()
! cccccccccccccccccccccccccccccccccccccccccccccccccc

  REAL(8), ALLOCATABLE :: angdistp(:),angdist(:,:)
  REAL(8), ALLOCATABLE ::  wt(:), theta(:)
  REAL(8), PARAMETER :: pi = 3.141592653589793D0
  INTEGER :: nth, nph, nrh, ind, nmaxrho
  REAL(8) :: drho, phi, rvecx, rvecy, rvecz, x1, y1, z1
  REAL(8) :: rxrel, ryrel, rzrel, rxdiff, rydiff, rzdiff
  REAL(8) :: total, wg, func, dmax2, deltatheta

  ALLOCATE(angdist(nstepstheta, nstepsphi))
  ALLOCATE(angdistp(nstepstheta),wt(nstepstheta),theta(nstepstheta))

  angdist = 0D0
  angdistp = 0D0
  drho = dx/8D0 ! empirical choice
  dmax2 = MAX(kxbox, kybox, kzbox)/2D0*dx
  nmaxrho = dmax2/drho
  deltatheta = pi/(nstepstheta - 1)

  WRITE (6, '(1a,1i4)') ' nmaxrho = ', nmaxrho
  WRITE (6, '(1a,1i4)') ' nstepstheta = ', nstepstheta
  WRITE (6, '(1a,1i4)') ' nstepsphi = ', nstepsphi
  WRITE (6, '(1a,1i4)') ' nabsorb = ', nabsorb
  WRITE (6, '(1a,1f8.4)') ' dx = ', dx
  WRITE (6, '(1a,1f8.4)') ' dmax2 = ', dmax2
  WRITE (6, '(1a,1f8.4)') ' drho = ', drho

! initialize theta and phi
  DO nth = 1, nstepstheta
    theta(nth) = deltatheta*(nth - 1)
    IF (nth == 1 .OR. nth == nstepstheta) THEN
      wt(nth) = 2D0*pi*(1D0 - cos(deltatheta/2D0))
    ELSE
      wt(nth) = 2D0*pi*(dcos(theta(nth) - deltatheta/2.) - dcos(theta(nth) + deltatheta/2.))
    END IF
    wt(nth) = dabs(wt(nth))
  END DO

WRITE(*,*) 'WT:',wt

! calculate angdist
  DO nth = 1, nstepstheta
    DO nph = 1, nstepsphi
      DO nrh = 1, nmaxrho
        phi = 2D0*pi*(nph - 1)/nstepsphi
        rho = drho*nrh

        rvecx = rho*dcos(phi)*dsin(theta(nth))
        rvecy = rho*dsin(phi)*dsin(theta(nth))
        rvecz = rho*dcos(theta(nth))
!WRITE(*,*) 'nth,nph,nrho=',nth,nph,nrho,rvecx,rvecz,rvecz

        IF (rvecx > kxbox/2D0*dx .OR. rvecx < -(kxbox - 2)/2D0*dx) goto 23 ! OUT of box
        IF (rvecy > kybox/2D0*dx .OR. rvecy < -(kybox - 2)/2D0*dx) goto 23 ! OUT of box
        IF (rvecz > kzbox/2D0*dx .OR. rvecz < -(kzbox - 2)/2D0*dx) goto 23 ! OUT of box
! IF(icheckabsozone(rvecx,rvecy,rvecz) == 0) WRITE(777,'(3f,i)') rvecx,rvecy,rvecz,0 ! OUT of abso zone
        IF (icheckabsozone(rvecx, rvecy, rvecz) == 0) goto 23 ! OUT of abso zone

!WRITE(777,'(3f,i)') rvecx,rvecy,rvecz,1

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
              angdist(nth, nph) = angdist(nth, nph) + rhoabso(ind)*rho**2D0*func*drho
            END DO
          END DO
        END DO

23      CONTINUE
      END DO ! nmaxrho
!      IF (nth == 1 .OR. nth == nstepstheta) EXIT
    END DO ! nstepsphi
  END DO ! nstepstheta

! WRITE OUT
  DO nth = 1, nstepstheta
    DO nph = 1, nstepsphi
      IF (nth == 1) angdist(nth, nph) = angdist(1, 1)
      IF (nth == nstepstheta) angdist(nth, nph) = angdist(nstepstheta, 1)
      WRITE (28, '(2f17.4,1e17.7)') (nph - 1)*360D0/nstepsphi, theta(nth)/pi*180D0, angdist(nth, nph)
    END DO
    WRITE (28, '(2f17.4,1e17.7)') 360D0, theta(nth)/pi*180D0, angdist(nth, 1)
    WRITE (28, *)
  END DO

! calculate phi average and PRINT OUT pangdist3-2
  total = 0.0
  wg = 0.0
  DO nth = 1, nstepstheta
    DO nph = 1, nstepsphi
      angdistp(nth) = angdistp(nth) + angdist(nth, nph)/nstepsphi
    END DO
    total = total + angdistp(nth)*wt(nth)
    wg = wg + wt(nth)
  END DO

  DO nth = 1, nstepstheta
    WRITE (29, '(1f17.4,2e17.7)') theta(nth)/pi*180D0, angdistp(nth)/(total/4D0/pi), &
      angdistp(nth)/4D0/pi
  END DO

! finish
  WRITE (*, '(1a,1e17.7)') ' total = ', total
  WRITE (*, '(1a,1f17.14)') ' wg = ', wg

  DEALLOCATE(angdist)
  DEALLOCATE(angdistp)


  RETURN

END SUBROUTINE angular_distribution

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


END PROGRAM angdist5
