
#if(coudoub)
!-----main--of--falr Coulomb solver-------------------------------------

! SUBROUTINE pois(rhoc,wcoul)

SUBROUTINE falr(rhoc, wcoul)
  USE singleton, ONLY: fftn, resetfft
  USE params
  USE util, ONLY: prifld
  USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

  REAL, INTENT(IN) :: rhoc(kdfull2)
  REAL, INTENT(OUT) :: wcoul(kdfull2)
  INTEGER, PARAMETER :: kx2 = 2*kxbox
  INTEGER, PARAMETER :: ky2 = 2*kybox
  INTEGER, PARAMETER :: kz2 = 2*kzbox

  COMPLEX, PARAMETER :: czero = (0.0, 0.0)

! solve the Poisson equation TO obtain the Coulomb potential at
! collocation points.
! -d**2*wc=4*pi*e2*rhoc

! output:
! wcoul : Coulomb potential at collocation points.

  COMPLEX :: rho2(kx2, ky2, kz2) ! density on large grid
  COMPLEX :: qinv(kx2, ky2, kz2) ! inverse Green FUNCTION IN q space
  SAVE qinv
  LOGICAL :: tfirst
  DATA tfirst/.true./
  LOGICAL :: ttest/.false./

!-----------------------------------------------------------------------

! initialize q-field

  IF (tfirst) THEN
    DO i3 = 1, kz2
      IF (i3 <= kz2/2) THEN
        qz = (dz*(i3 - 1))**2
      ELSE
        qz = (dz*(i3 - kz2 - 1))**2
      END IF
      DO i2 = 1, ky2
        IF (i2 <= ky2/2) THEN
          qy = (dy*(i2 - 1))**2
        ELSE
          qy = (dy*(i2 - ky2 - 1))**2
        END IF
        DO i1 = 1, kx2
          IF (i1 <= kx2/2) THEN
            qx = (dx*(i1 - 1))**2
          ELSE
            qx = (dx*(i1 - kx2 - 1))**2
          END IF
          IF (i1 /= 1 .OR. i2 /= 1 .OR. i3 /= 1) THEN
            qinv(i1, i2, i3) = e2*(dx*dy*dz)/SQRT(qx + qy + qz)
          ELSE
            qinv(i1, i2, i3) = 0.0
          END IF
        END DO
      END DO
    END DO
!c kz0 = kzbox/2
!c ky0 = kybox/2
!c kx0 = kxbox/2
!c kzcut = 3*kz0
!c kycut = 3*ky0
!c kxcut = 3*kx0
!c DO i3=1,kz2
!c IF(i3.le.kzcut) THEN
!c qz=(dz*(i3-kz0))**2
!c ELSE
!c qz=(dz*(i3-kz0-kz2))**2
!c END IF
!c IF(ttest) WRITE(6,*) ' i3,qz=',i3,qz
!c DO i2=1,ky2
!c IF(i2.le.kycut) THEN
!c qy=(dy*(i2-ky0))**2
!c ELSE
!c qy=(dy*(i2-ky2-ky0))**2
!c END IF
!c DO i1=1,kx2
!c IF(i1.le.kxcut) THEN
!c qx=(dx*(i1-kx0))**2
!c ELSE
!c qx=(dx*(i1-kx2-kx0))**2
!c END IF
!c IF(i1.NE.kx0 .OR. i2.NE.ky0 .OR. i3.NE.kz0) THEN
!c qinv(i1,i2,i3) = e2*(dx*dy*dz)/sqrt(qx+qy+qz)
!c ELSE
!c qinv(i1,i2,i3) = 0.0
!c END IF
!c END DO
!c END DO
!c END DO
    IF (ttest) WRITE (6, '(a,2000(/i5,2g12.4))') &
      ' Qr=', (i1, qinv(i1, ky0, kz0), i1=1, kx2)
    CALL resetfft
    CALL fftn(qinv, shape(qinv), arrayout=qinv)
    IF (ttest) WRITE (6, '(a,20(/5g13.5))') ' Qq=', (qinv(i1, 1, 1), i1=1, nx2)
    CALL resetfft
    tfirst = .false.
  END IF

! reset auxiliary density field and fill relevant part

  DO i3 = 1, kz2
    DO i2 = 1, ky2
      DO i1 = 1, kx2
        rho2(i1, i2, i3) = czero
      END DO
    END DO
  END DO
  IF (ttest) acc = 0D0
  ind = 0
  DO i3 = 1, kzbox
    DO i2 = 1, kybox
      DO i1 = 1, kxbox
        ind = 1 + ind
        rho2(i1, i2, i3) = rhoc(ind)
        IF (ttest) acc = acc + rhoc(ind)
      END DO
    END DO
  END DO
  IF (ttest) WRITE (6, *) ' COUL: charge=', acc*dvol
  IF (ttest) WRITE (6, '(a,2000(/i5,2g12.4))') &
    ' rho2=', (i1, rho2(i1, ky0, kz0), i1=1, kx2)

! the central solver:
! forward FFT --> division by q**2 --> backward FFT

  CALL fftn(rho2, shape(rho2), arrayout=rho2)
  IF (ttest) WRITE (6, '(a,20(/5g13.5))') ' rho2=', (rho2(i1, 1, 1), i1=1, nx2)
  DO i3 = 1, kz2
    DO i2 = 1, ky2
      DO i1 = 1, kx2
        rho2(i1, i2, i3) = rho2(i1, i2, i3)*qinv(i1, i2, i3)
      END DO
    END DO
  END DO
  CALL fftn(rho2, shape(rho2), inv=.true., arrayout=rho2)
  IF (ttest) WRITE (6, '(a,2000(/i5,2g12.4))') &
    ' w2=', (i1, rho2(i1, ky0, kz0), i1=1, kx2)
  ind = 0
  DO i3 = 1, kzbox
    DO i2 = 1, kybox
      DO i1 = 1, kxbox
        ind = 1 + ind
        wcoul(ind) = rho2(i1, i2, i3)
      END DO
    END DO
  END DO
  IF (ttest) CALL prifld(wcoul, ' coulomb ')

  RETURN
END SUBROUTINE falr
#else

SUBROUTINE pois()
! dummy
  RETURN
END SUBROUTINE pois
#endif

