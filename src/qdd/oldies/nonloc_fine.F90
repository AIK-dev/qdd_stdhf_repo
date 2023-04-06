
! This version CONTAINS attempts TO represent projectors on a fine
! grid (commented OUT, at the END).
! It CONTAINS also the old calling SEQUENCE depending on the "row"
! of an element.

! *****************************

#ifdef REALSWITCH
SUBROUTINE nonlocalr(psi, aux)
#else
  SUBROUTINE nonlocalc(psi, aux, ionact)
#endif

! *****************************

    USE params
    IMPLICIT NONE
#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: psi(kdfull2)
    REAL(DP), INTENT(IN OUT) :: aux(kdfull2)
#else
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2)
    COMPLEX(DP), INTENT(IN OUT) :: aux(kdfull2)
    INTEGER, INTENT(IN) :: ionact
    INTEGER :: maxion, minion

#endif

    INTEGER :: i, ion

#ifdef REALSWITCH

    DO i = 1, kdfull2
      aux(i) = 0D0
    END DO
    DO ion = 1, nion
! npact = np(ion)
! IF(nrow(npact) <= 2) CALL effnl2r(psi,aux,ion)
! IF(nrow(npact) == 3) CALL effnl3r(psi,aux,ion)
! IF(nrow(npact) == 4) CALL effnl4r(psi,aux,ion)
      IF (tnonloc(ion)) CALL effnl4r(psi, aux, ion)
    END DO

#else

    DO i = 1, kdfull2
      aux(i) = CMPLX(0D0, 0D0, DP)
    END DO
    IF (ionact == 0) THEN
      minion = 1
      maxion = nion
    ELSE IF (ionact <= nion .AND. ionact > 0) THEN
      minion = ionact
      maxion = ionact
    ELSE
      STOP ' NONLOCALC: ion OUT of range'
    END IF
    DO ion = minion, maxion
! npact = np(ion)
! IF(nrow(npact) <= 2) CALL effnl2c(psi,aux,ion)
! IF(nrow(npact) == 3) CALL effnl3c(psi,aux,ion)
! IF(nrow(npact) == 4) CALL effnl4c(psi,aux,ion)
      IF (tnonloc(ion)) CALL effnl4c(psi, aux, ion)
    END DO
#endif
    RETURN
#ifdef REALSWITCH
  END SUBROUTINE nonlocalr
#else
END SUBROUTINE nonlocalc
#endif

! **************************

#ifdef REALSWITCH
SUBROUTINE effnl2r(psi, aux, ion)
#else
  SUBROUTINE effnl2c(psi, aux, ion)
#endif

! **************************

    USE params
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN OUT) :: aux(kdfull2)
    REAL(DP), INTENT(IN) :: psi(kdfull2)
#else
    COMPLEX(DP), INTENT(IN OUT) :: aux(kdfull2)
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2)
#endif
    INTEGER, INTENT(IN) :: ion

    INTEGER :: i, ii, nfin
    REAL(DP) :: h0_11, h1_11
#ifdef REALSWITCH
    REAL(DP) :: sumr01, sumr04, sumr05, sumr06
    REAL(DP) :: f0_11, f1_1x, f1_1y, f1_1z
#else
    COMPLEX(DP) :: sumr01, sumr04, sumr05, sumr06
    COMPLEX(DP) :: f0_11, f1_1x, f1_1y, f1_1z
#endif

    h0_11 = h0_11g(np(ion))
    h1_11 = h1_11g(np(ion))

! compute all the integrals over r'-dependent part of the
! wavefunction times the r'-dependent part of the pseudo:

! **** a) l=0 -> p0_1, p0_2 (no LS-coupling) ****

    sumr01 = 0D0
    sumr04 = 0D0
    sumr05 = 0D0
    sumr06 = 0D0

    nfin = ifin(ion)
    DO i = 1, nfin
      ii = icount(i, ion)
      sumr01 = sumr01 + psi(ii)*p0_1(i, ion)
    END DO
    sumr01 = sumr01*dvol

! **** b) l=1 -> p1_1x,y,z (including LS-coupling) ****

    DO i = 1, nfin
      ii = icount(i, ion)
      sumr04 = sumr04 + psi(ii)*p1_1x(i, ion)
      sumr05 = sumr05 + psi(ii)*p1_1y(i, ion)
      sumr06 = sumr06 + psi(ii)*p1_1z(i, ion)
    END DO
    sumr04 = sumr04*dvol
    sumr05 = sumr05*dvol
    sumr06 = sumr06*dvol

! Now, the RESULT of the r' integral is multiplied with the
! r-dependent part and summed up over all l-contributions as well
! as over all the ions:

    f0_11 = e2*(0.25D0/pi)*h0_11*sumr01
    IF (h1_11 == 0D0) THEN
      DO i = 1, nfin
        ii = icount(i, ion)
        aux(ii) = aux(ii) + p0_1(i, ion)*f0_11
      END DO
    ELSE
      f1_1x = e2*0.75D0/pi*h1_11*sumr04
      f1_1y = e2*0.75D0/pi*h1_11*sumr05
      f1_1z = e2*0.75D0/pi*h1_11*sumr06
      DO i = 1, nfin
        ii = icount(i, ion)
        aux(ii) = aux(ii) &
                  + p0_1(i, ion)*f0_11 &
                  + p1_1x(i, ion)*f1_1x + p1_1y(i, ion)*f1_1y + p1_1z(i, ion)*f1_1z
! + e2 * (1/(4D0*pi)*p0_1(i,ion)*h0_11*sumr01 &
! + 3D0/(4D0*pi)*h1_11*(p1_1x(i,ion)*sumr04+ p1_1y(i,ion)*sumr05+ &
! p1_1z(i,ion)*sumr06))
      END DO
    END IF

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE effnl2r
#else
END SUBROUTINE effnl2c
#endif
! **************************

#ifdef REALSWITCH
SUBROUTINE effnl3r(psi, aux, ion)
#else
  SUBROUTINE effnl3c(psi, aux, ion)
#endif

! **************************

    USE params
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN OUT) :: aux(kdfull2)
    REAL(DP), INTENT(IN) :: psi(kdfull2)
#else
    COMPLEX(DP), INTENT(IN OUT) :: aux(kdfull2)
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2)
#endif
    INTEGER, INTENT(IN) :: ion

    INTEGER :: i, ii, nfin
    REAL(DP) :: h0_11, h0_12, h0_21, h0_22, h1_11
#ifdef REALSWITCH
    REAL(DP) :: sumr01, sumr02, sumr04, sumr05, sumr06
    REAL(DP) :: f0_11, f0_22, f1_1x, f1_1y, f1_1z
#else
    COMPLEX(DP) :: sumr01, sumr02, sumr04, sumr05, sumr06
    COMPLEX(DP) :: f0_11, f0_22, f1_1x, f1_1y, f1_1z
#endif

    h0_11 = h0_11g(np(ion))
    h0_22 = h0_22g(np(ion))
    h1_11 = h1_11g(np(ion))

    IF (h0_12g(np(ion)) .gt. -1D10) THEN
      h0_12 = h0_12g(np(ion))
    ELSE
      h0_12 = -0.5D0*SQRT(3D0/5D0)*h0_22
    END IF
    h0_21 = h0_12

! compute all the integrals over r'-dependent part of the
! wavefunction times the r'-dependent part of the pseudo:

! **** a) l=0 -> p0_1, p0_2 (no LS-coupling) ****

    sumr01 = 0D0
    sumr02 = 0D0
    sumr04 = 0D0
    sumr05 = 0D0
    sumr06 = 0D0

    nfin = ifin(ion)
    DO i = 1, nfin
      ii = icount(i, ion)
      sumr01 = sumr01 + psi(ii)*p0_1(i, ion)
      sumr02 = sumr02 + psi(ii)*p0_2(i, ion)
    END DO
    sumr01 = sumr01*dvol
    sumr02 = sumr02*dvol

! **** b) l=1 -> p1_1x,y,z (including LS-coupling) ****

    DO i = 1, nfin
      ii = icount(i, ion)
      sumr04 = sumr04 + psi(ii)*p1_1x(i, ion)
      sumr05 = sumr05 + psi(ii)*p1_1y(i, ion)
      sumr06 = sumr06 + psi(ii)*p1_1z(i, ion)
    END DO
    sumr04 = sumr04*dvol
    sumr05 = sumr05*dvol
    sumr06 = sumr06*dvol
! Now, the RESULT of the r' integral is multiplied with the
! r-dependant part and summed up over all l-contributions as well
! as over all the ions:

    f0_11 = e2*(0.25D0/pi)*(h0_11*sumr01 + h0_12*sumr02)
    f0_22 = e2*(0.25D0/pi)*(h0_21*sumr01 + h0_22*sumr02)
    IF (h1_11 == 0D0) THEN
      DO i = 1, nfin
        ii = icount(i, ion)
        aux(ii) = aux(ii) &
                  + p0_1(i, ion)*f0_11 + p0_2(i, ion)*f0_22
      END DO
    ELSE
      f1_1x = e2*0.75D0/pi*h1_11*sumr04
      f1_1y = e2*0.75D0/pi*h1_11*sumr05
      f1_1z = e2*0.75D0/pi*h1_11*sumr06
      DO i = 1, nfin
        ii = icount(i, ion)
        aux(ii) = aux(ii) &
                  + p0_1(i, ion)*f0_11 + p0_2(i, ion)*f0_22 &
                  + p1_1x(i, ion)*f1_1x + p1_1y(i, ion)*f1_1y + p1_1z(i, ion)*f1_1z
! aux(ii)= aux(ii) + e2 * (1/(4*pi)*(p0_1(i,ion)*h0_11*sumr01+ &
! p0_1(i,ion)*h0_12*sumr02+ p0_2(i,ion)*h0_21*sumr01+ &
! p0_2(i,ion)*h0_22*sumr02) +3D0/(4D0*pi)*h1_11*(p1_1x(i,ion)*sumr04+ &
! p1_1y(i,ion)*sumr05+ p1_1z(i,ion)*sumr06))
      END DO
    END IF

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE effnl3r
#else
END SUBROUTINE effnl3c
#endif

! **************************

#ifdef REALSWITCH
SUBROUTINE effnl4r(psi, aux, ion)
#else
  SUBROUTINE effnl4c(psi, aux, ion)
#endif

! **************************

    USE params
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN OUT) :: aux(kdfull2)
    REAL(DP), INTENT(IN) :: psi(kdfull2)
#else
    COMPLEX(DP), INTENT(IN OUT) :: aux(kdfull2)
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2)
#endif
    INTEGER, INTENT(IN) :: ion

    INTEGER :: i, ii, nfin
    REAL(DP) :: h0_11, h0_12, h0_13, h0_21, h0_22, h0_23, h0_31, h0_32, h0_33, h1_11, h1_12, h1_21, h1_22, h2_11
#ifdef REALSWITCH
    REAL(DP) :: sumr01, sumr02, sumr03, sumr04, sumr05, sumr06
    REAL(DP) :: sumr07, sumr08, sumr09
    REAL(DP) :: sumr16, sumr17, sumr18, sumr19, sumr20
    REAL(DP) :: f0_1, f0_2, f0_3, f1_1x, f1_1y, f1_1z, f1_2x, f1_2y, f1_2z, &
                f2_xy, f2_xz, f2_yz, f2_xy2, f2_z2
#else
    COMPLEX(DP) :: sumr01, sumr02, sumr03, sumr04, sumr05, sumr06
    COMPLEX(DP) :: sumr07, sumr08, sumr09
    COMPLEX(DP) :: sumr16, sumr17, sumr18, sumr19, sumr20
    COMPLEX(DP) :: f0_1, f0_2, f0_3, f1_1x, f1_1y, f1_1z, f1_2x, f1_2y, f1_2z, &
                   f2_xy, f2_xz, f2_yz, f2_xy2, f2_z2
#endif

    REAL(DP), PARAMETER :: fac1_12 = -0.422577127364258D0 ! -0.5D0*SQRT(5D0/7D0)
    REAL(DP), PARAMETER :: fac0_23 = -0.629940788348712D0 ! -0.5D0*SQRT(100D0/63D0)
    REAL(DP), PARAMETER :: fac0_13 = 0.243975018237133D0 ! 0.5D0*SQRT(5D0/21D0)
    REAL(DP), PARAMETER :: fac0_12 = -0.387298334620742D0 ! -0.5D0*SQRT(3D0/5D0)

    h0_11 = h0_11g(np(ion))

    h0_22 = h0_22g(np(ion))
    h1_11 = h1_11g(np(ion))

    h0_33 = h0_33g(np(ion))
    h1_22 = h1_22g(np(ion))
    h2_11 = h2_11g(np(ion))

    IF (h0_12g(np(ion)) .gt. -1D10) THEN
      h0_12 = h0_12g(np(ion))
    ELSE
      h0_12 = fac0_12*h0_22
    END IF
    h0_21 = h0_12

    h0_13 = fac0_13*h0_33
    h0_31 = h0_13

    h0_23 = fac0_23*h0_33
    h0_32 = h0_23

    h1_12 = fac1_12*h1_22
    h1_21 = h1_12

! initialize the RESULT of the integration:

    sumr01 = 0D0
    sumr02 = 0D0
    sumr03 = 0D0
    sumr04 = 0D0
    sumr05 = 0D0
    sumr06 = 0D0
    sumr07 = 0D0
    sumr08 = 0D0
    sumr09 = 0D0
    sumr16 = 0D0
    sumr17 = 0D0
    sumr18 = 0D0
    sumr19 = 0D0
    sumr20 = 0D0

    nfin = ifin(ion)

! compute all the integrals over r'-dependant part of the
! wavefunction times the r'-dependent part of the pseudo:

    IF (ABS(h0_11) > small) THEN
      DO i = 1, nfin
        ii = icount(i, ion)
        sumr01 = sumr01 + psi(ii)*p0_1(i, ion)
      END DO
      sumr01 = sumr01*dvol
    ELSE IF (ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > small) THEN
      DO i = 1, nfin
        ii = icount(i, ion)
        sumr02 = sumr02 + psi(ii)*p0_2(i, ion)
        sumr04 = sumr04 + psi(ii)*p1_1x(i, ion)
        sumr05 = sumr05 + psi(ii)*p1_1y(i, ion)
        sumr06 = sumr06 + psi(ii)*p1_1z(i, ion)
      END DO
      sumr02 = sumr02*dvol
      sumr04 = sumr04*dvol
      sumr05 = sumr05*dvol
      sumr06 = sumr06*dvol
    ELSE IF (ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > small) THEN
      DO i = 1, nfin
        ii = icount(i, ion)
        sumr03 = sumr03 + psi(ii)*p0_3(i, ion)
        sumr07 = sumr07 + psi(ii)*p1_2x(i, ion)
        sumr08 = sumr08 + psi(ii)*p1_2y(i, ion)
        sumr09 = sumr09 + psi(ii)*p1_2z(i, ion)
        sumr16 = sumr16 + psi(ii)*p2_xy(i, ion)
        sumr17 = sumr17 + psi(ii)*p2_xz(i, ion)
        sumr18 = sumr18 + psi(ii)*p2_yz(i, ion)
        sumr19 = sumr19 + psi(ii)*p2_xy2(i, ion)
        sumr20 = sumr20 + psi(ii)*p2_z2(i, ion)
      END DO
      sumr03 = sumr03*dvol
      sumr07 = sumr07*dvol
      sumr08 = sumr08*dvol
      sumr09 = sumr09*dvol
      sumr16 = sumr16*dvol
      sumr17 = sumr17*dvol
      sumr18 = sumr18*dvol
      sumr19 = sumr19*dvol
      sumr20 = sumr20*dvol
    END IF

! Now, the RESULT of the r' integral is multiplied with the
! r-dependent part and summed up over all l-contributions

    IF (ABS(h0_11) > small) THEN
      f0_1 = e2*(h0_11*sumr01 + h0_12*sumr02 + h0_13*sumr03)/(4D0*pi)
! IF(ion==1) WRITE(6,'(a,4(1pg13.5))') ' norm of 1s projector,h0:', &
! dvol*SUM(p0_1(:,ion)**2)/(4D0*pi),(f0_1),abs(f0_1)
!CALL FLUSH(6)
! WRITE(*,*) ' f0_1=',f0_1,sumr01
      DO i = 1, nfin
        ii = icount(i, ion)
        aux(ii) = aux(ii) + p0_1(i, ion)*f0_1
      END DO
    ELSE IF (ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > small) THEN
      f0_2 = e2*(h0_21*sumr01 + h0_22*sumr02 + h0_23*sumr03)/(4D0*pi)
      f1_1x = e2*(h1_11*sumr04 + h1_12*sumr07)*3D0/(4D0*pi)
      f1_1y = e2*(h1_11*sumr05 + h1_12*sumr08)*3D0/(4D0*pi)
      f1_1z = e2*(h1_11*sumr06 + h1_12*sumr09)*3D0/(4D0*pi)
      DO i = 1, nfin
        ii = icount(i, ion)
        aux(ii) = aux(ii) + p0_2(i, ion)*f0_2 &
                  + (p1_1x(i, ion)*f1_1x + p1_1y(i, ion)*f1_1y + p1_1z(i, ion)*f1_1z)
      END DO
    ELSE IF (ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > small) THEN
      f0_3 = e2*(h0_31*sumr01 + h0_32*sumr02 + h0_33*sumr03)/(4D0*pi)
      f1_2x = e2*(h1_21*sumr04 + h1_22*sumr07)*3D0/(4D0*pi)
      f1_2y = e2*(h1_21*sumr05 + h1_22*sumr08)*3D0/(4D0*pi)
      f1_2z = e2*(h1_21*sumr06 + h1_22*sumr09)*3D0/(4D0*pi)
      f2_xy = e2*h2_11*sumr16*15D0/(4D0*pi)
      f2_xz = e2*h2_11*sumr17*15D0/(4D0*pi)
      f2_yz = e2*h2_11*sumr18*15D0/(4D0*pi)
      f2_xy2 = e2*h2_11*sumr19*5D0/(16D0*pi)
      f2_z2 = e2*h2_11*sumr20*5D0/(16D0*pi)
      DO i = 1, nfin
        ii = icount(i, ion)
        aux(ii) = aux(ii) + p0_3(i, ion)*f0_3 &
                  + (p1_2x(i, ion)*f1_2x + p1_2y(i, ion)*f1_2y + p1_2z(i, ion)*f1_2z) &
                  + (p2_xy(i, ion)*f2_xy + p2_xz(i, ion)*f2_xz + p2_yz(i, ion)*f2_yz) &
                  + (p2_xy2(i, ion)*f2_xy2 + p2_z2(i, ion)*f2_z2)
      END DO
    END IF

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE effnl4r
#else
END SUBROUTINE effnl4c
#endif

!~ ! *****************************
!~ #ifdef REALSWITCH
!~ SUBROUTINE nonlocalrfine(psi,aux)
!~ #else
!~ SUBROUTINE nonlocalcfine(psi,aux,ionact)
!~ #endif

!~ ! *****************************

!~ USE params
!~ IMPLICIT REAL(DP) (A-H,O-Z)

!~ #ifdef REALSWITCH
!~ REAL(DP), INTENT(IN ) :: psi(kdfull2fine)
!~ REAL(DP), INTENT(OUT ) :: aux(kdfull2)
!~ DO i=1,kdfull2fine
!~ aux(i)=0D0
!~ END DO
!~ !test WRITE(*,*) ' IN NONLOCALR'
!~ DO ion=1,nion
!~ ! npact = np(ion)
!~ ! IF(nrow(npact) <= 2) CALL effnl2r(psi,aux,ion)
!~ ! IF(nrow(npact) == 3) CALL effnl3r(psi,aux,ion)
!~ ! IF(nrow(npact) == 4) CALL effnl4r(psi,aux,ion)
!~ IF(tnonloc(ion)) CALL effnl4rfine(psi,aux,ion)
!~ END DO
!~ #else
!~ COMPLEX(DP), INTENT(IN) :: psi(kdfull2fine)
!~ COMPLEX(DP), INTENT(OUT) :: aux(kdfull2)
!~ DO i=1,kdfull2
!~ aux(i)=CMPLX(0D0,0D0,DP)
!~ END DO
!~ IF(ionact == 0) THEN
!~ minion = 1
!~ maxion = nion
!~ ELSE IF(ionact <= nion .AND. ionact > 0) THEN
!~ minion = ionact
!~ maxion = ionact
!~ ELSE
!~ STOP ' NONLOCALC: ion OUT of range'
!~ END IF
!~ DO ion=minion,maxion
!~ ! npact = np(ion)
!~ ! IF(nrow(npact) <= 2) CALL effnl2c(psi,aux,ion)
!~ ! IF(nrow(npact) == 3) CALL effnl3c(psi,aux,ion)
!~ ! IF(nrow(npact) == 4) CALL effnl4c(psi,aux,ion)
!~ IF(tnonloc(ion)) CALL effnl4cfine(psi,aux,ion)
!~ END DO
!~ #endif

!~ RETURN
!~ #ifdef REALSWITCH
!~ END SUBROUTINE nonlocalrfine
!~ #else
!~ END SUBROUTINE nonlocalcfine
!~ #endif

!~ ! **************************
!~ #ifdef REALSWITCH
!~ SUBROUTINE effnl4rfine(psi,aux,ion)
!~ #else
!~ SUBROUTINE effnl4cfine(psi,aux,ion)
!~ #endif

!~ ! **************************

!~ USE params
!~ IMPLICIT REAL(DP) (A-H,O-Z)
!~ #ifdef REALSWITCH
!~ REAL(DP), INTENT(IN OUT) :: aux(kdfull2)
!~ REAL(DP), INTENT(IN ) :: psi(kdfull2fine)
!~ REAL(DP) :: sumr01,sumr02,sumr03,sumr04,sumr05,sumr06
!~ REAL(DP) :: sumr07,sumr08,sumr09
!~ REAL(DP) :: sumr16,sumr17,sumr18,sumr19,sumr20
!~ REAL(DP) :: f0_1,f0_2,f0_3,f1_1x,f1_1y,f1_1z,f1_2x,f1_2y,f1_2z, &
!~ f2_xy,f2_xz,f2_yz,f2_xy2,f2_z2
!~ #else
!~ COMPLEX(DP), INTENT(IN OUT) :: aux(kdfull2)
!~ COMPLEX(DP), INTENT(IN ) :: psi(kdfull2fine)
!~ COMPLEX(DP) :: sumr01,sumr02,sumr03,sumr04,sumr05,sumr06
!~ COMPLEX(DP) :: sumr07,sumr08,sumr09
!~ COMPLEX(DP) :: sumr16,sumr17,sumr18,sumr19,sumr20
!~ COMPLEX(DP) :: f0_1,f0_2,f0_3,f1_1x,f1_1y,f1_1z,f1_2x,f1_2y,f1_2z, &
!~ f2_xy,f2_xz,f2_yz,f2_xy2,f2_z2
!~ #endif

!~ REAL(DP),PARAMETER :: fac1_12=-0.422577127364258 ! -0.5D0*SQRT(5D0/7D0)
!~ REAL(DP),PARAMETER :: fac0_23=-0.629940788348712 ! -0.5D0*SQRT(100D0/63D0)
!~ REAL(DP),PARAMETER :: fac0_13= 0.243975018237133 ! 0.5D0*SQRT(5D0/21D0)
!~ REAL(DP),PARAMETER :: fac0_12=-0.387298334620742 ! -0.5D0*SQRT(3D0/5D0)

!~ dvolfine=dx*dy*dz/8.0

!~ h0_11=h0_11g(np(ion))

!~ h0_22=h0_22g(np(ion))
!~ h1_11=h1_11g(np(ion))

!~ h0_33=h0_33g(np(ion))
!~ h1_22=h1_22g(np(ion))
!~ h2_11=h2_11g(np(ion))

!~ IF(h0_12g(np(ion)).gt.-1D10) THEN
!~ h0_12=h0_12g(np(ion))
!~ ELSE
!~ h0_12 = fac0_12*h0_22
!~ END IF
!~ h0_21 = h0_12

!~ h0_13 = fac0_13*h0_33
!~ h0_31 = h0_13

!~ h0_23 = fac0_23*h0_33
!~ h0_32 = h0_23

!~ h1_12 = fac1_12*h1_22
!~ h1_21 = h1_12

!~ ! initialize the RESULT of the integration:

!~ sumr01 = 0D0
!~ sumr02 = 0D0
!~ sumr03 = 0D0
!~ sumr04 = 0D0
!~ sumr05 = 0D0
!~ sumr06 = 0D0
!~ sumr07 = 0D0
!~ sumr08 = 0D0
!~ sumr09 = 0D0
!~ sumr16 = 0D0
!~ sumr17 = 0D0
!~ sumr18 = 0D0
!~ sumr19 = 0D0
!~ sumr20 = 0D0

!~ nfin=ifinfine(ion)

!~ ! compute all the integrals over r'-dependant part of the
!~ ! wavefunction times the r'-dependent part of the pseudo:

!~ IF(ABS(h0_11) > small) THEN
!~ DO i=1,nfin
!~ ii = icountfine(i,ion)
!~ sumr01 = sumr01 + psi(ii)*p0_1fine(i,ion)
!~ END DO
!~ sumr01=sumr01*dvolfine
!~ ELSE IF(ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > small) THEN
!~ DO i=1,nfin
!~ ii = icountfine(i,ion)
!~ sumr02 = sumr02 + psi(ii)*p0_2fine(i,ion)
!~ sumr04 = sumr04 + psi(ii)*p1_1xfine(i,ion)
!~ sumr05 = sumr05 + psi(ii)*p1_1yfine(i,ion)
!~ sumr06 = sumr06 + psi(ii)*p1_1zfine(i,ion)
!~ END DO
!~ sumr02=sumr02*dvolfine
!~ sumr04=sumr04*dvolfine
!~ sumr05=sumr05*dvolfine
!~ sumr06=sumr06*dvolfine
!~ ELSE IF(ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > small) THEN
!~ DO i=1,nfin
!~ ii = icountfine(i,ion)
!~ sumr03 = sumr03 + psi(ii)*p0_3fine(i,ion)
!~ sumr07 = sumr07 + psi(ii)*p1_2xfine(i,ion)
!~ sumr08 = sumr08 + psi(ii)*p1_2yfine(i,ion)
!~ sumr09 = sumr09 + psi(ii)*p1_2zfine(i,ion)
!~ sumr16 = sumr16 + psi(ii)*p2_xyfine(i,ion)
!~ sumr17 = sumr17 + psi(ii)*p2_xzfine(i,ion)
!~ sumr18 = sumr18 + psi(ii)*p2_yzfine(i,ion)
!~ sumr19 = sumr19 + psi(ii)*p2_xy2fine(i,ion)
!~ sumr20 = sumr20 + psi(ii)*p2_z2fine(i,ion)
!~ END DO
!~ sumr03=sumr03*dvolfine
!~ sumr07=sumr07*dvolfine
!~ sumr08=sumr08*dvolfine
!~ sumr09=sumr09*dvolfine
!~ sumr16=sumr16*dvolfine
!~ sumr17=sumr17*dvolfine
!~ sumr18=sumr18*dvolfine
!~ sumr19=sumr19*dvolfine
!~ sumr20=sumr20*dvolfine
!~ END IF

!~ ! Now, the RESULT of the r' integral is multiplied with the
!~ ! r-dependent part and summed up over all l-contributions

!~ IF(ABS(h0_11) > small) THEN
!~ f0_1 = e2*( h0_11*sumr01 + h0_12*sumr02 + h0_13*sumr03 )/(4D0*pi)
!~ ! WRITE(6,'(a,4(1pg13.5))') ' norm of 1s projector,h0:', &
!~ ! dvolfine*SUM(p0_1(:,ion)**2)/(4D0*pi),(f0_1),abs(f0_1)
!~ ! WRITE(6,*) h0_11,sumr01,h0_12,sumr02,h0_13,sumr03

!~ nfin=ifin(ion)
!~ !CALL FLUSH(6)
!~ !WRITE(*,*) ' f0_1=',f0_1,sumr01
!~ !amax=-1e10
!~ DO i=1,nfin
!~ ii = icount(i,ion)
!~ aux(ii)= aux(ii) + p0_1(i,ion)*f0_1
!~ ! IF(abs(aux(ii)).gt.amax) THEN
!~ ! amax=aux(ii)
!~ ! imax=ii
!~ ! END IF
!~ END DO
!~ !ii=imax
!~ ! i1=ii/nxyf
!~ ! ii2=mod(ii,nxyf)
!~ ! i2=ii2/nx2
!~ ! i3=mod(i2,nx2)
!~ ! WRITE(166,*) ion,i,ii,i1,i2,i3,amax
!~ !WRITE(166,*) ssum,ion,sumr08,sumr09,sumr10,sumr11,sumr12,sumr13,sumr14
!~ !WRITE(166,*) ssum,ion,sumr15,sumr16,sumr17,sumr18,sumr19,sumr20
!~ ELSE IF(ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > small) THEN
!~ f0_2 = e2*( h0_21*sumr01 + h0_22*sumr02 + h0_23*sumr03 )/(4D0*pi)
!~ f1_1x = e2*( h1_11*sumr04 + h1_12*sumr07 )*3D0/(4D0*pi)
!~ f1_1y = e2*( h1_11*sumr05 + h1_12*sumr08 )*3D0/(4D0*pi)
!~ f1_1z = e2*( h1_11*sumr06 + h1_12*sumr09 )*3D0/(4D0*pi)
!~ DO i=1,nfin
!~ ii = icount(i,ion)
!~ !WRITE(166,*) i,ii,aux(ii)
!~ aux(ii)= aux(ii) + p0_2(i,ion)*f0_2 &
!~ +(p1_1x(i,ion)*f1_1x + p1_1y(i,ion)*f1_1y + p1_1z(i,ion)*f1_1z)
!~ END DO
!~ ELSE IF(ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > small) THEN
!~ f0_3 = e2*( h0_31*sumr01 + h0_32*sumr02 + h0_33*sumr03 )/(4D0*pi)
!~ f1_2x = e2*( h1_21*sumr04 + h1_22*sumr07 )*3D0/(4D0*pi)
!~ f1_2y = e2*( h1_21*sumr05 + h1_22*sumr08 )*3D0/(4D0*pi)
!~ f1_2z = e2*( h1_21*sumr06 + h1_22*sumr09 )*3D0/(4D0*pi)
!~ f2_xy = e2* h2_11*sumr16*15D0/(4D0*pi)
!~ f2_xz = e2* h2_11*sumr17*15D0/(4D0*pi)
!~ f2_yz = e2* h2_11*sumr18*15D0/(4D0*pi)
!~ f2_xy2 = e2* h2_11*sumr19*5D0/(16D0*pi)
!~ f2_z2 = e2* h2_11*sumr20*5D0/(16D0*pi)
!~ DO i=1,nfin
!~ ii = icount(i,ion)
!~ aux(ii)= aux(ii) + p0_3(i,ion)*f0_3 &
!~ +(p1_2x(i,ion)*f1_2x + p1_2y(i,ion)*f1_2y + p1_2z(i,ion)*f1_2z) &
!~ +(p2_xy(i,ion)*f2_xy + p2_xz(i,ion)*f2_xz + p2_yz(i,ion)*f2_yz) &
!~ +(p2_xy2(i,ion)*f2_xy2 + p2_z2(i,ion)*f2_z2)
!~ END DO
!~ END IF

!~ RETURN
!~ #ifdef REALSWITCH
!~ END SUBROUTINE effnl4rfine
!~ #else
!~ END SUBROUTINE effnl4cfine
!~ #endif
