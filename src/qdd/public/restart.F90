
! **************************

#ifdef REALSWITCH
SUBROUTINE resume(psi, outna)
#else
  SUBROUTINE restart2(psi, outna, tstatin)
#endif

! **************************

! Reads data on wavefunctions, ions, and fields.
! 'resume' is the version for REAL wavefunctions (static).
! The variant 'restart2' produces COMPLEX wavefunctions,
! for 'trealin=.false.' from saved COMPLEX wavefunctions and
! for 'trealin=.true.' from REAL wavefunctions converting them
! to COMPLEX after reading. The parameter 'trealin' is read from
! the file 'RSAVE.*'.
! The list parameter 'tstatin' regulates reading of static or
! dynamic input.
! All the data are saved in one file called '(R)SAVE', also in the
! parallel case.

    USE params
    USE kinetic
!#ifdef REALSWITCH
!#endif
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
    LOGICAL, PARAMETER :: tstatin = .false.
#else
    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
    LOGICAL, INTENT(IN) :: tstatin
#endif

    CHARACTER(LEN=13), INTENT(IN) :: outna

    INTEGER :: iact, nstate_test, mynact, n, nb
    REAL(DP) :: dummy
    REAL(DP), ALLOCATABLE :: psiauxr(:)
    LOGICAL :: trealin
    LOGICAL, PARAMETER :: ttest = .TRUE.

#ifdef COMPLEXSWITCH
    INTEGER :: nbe
#endif
    LOGICAL :: topenf


    mynact = myn

#ifdef REALSWITCH

    IF (mynact == 0) &
      OPEN (UNIT=ifile, STATUS='old', FORM='UNFORMATTED', FILE='RSAVE.'//outna)

#else

    IF (mynact == 0) THEN
      IF (tstatin) THEN
        INQUIRE (ifile, OPENED=topenf)
        IF (.NOT. topenf) THEN
          OPEN (UNIT=ifile, STATUS='old', FORM='UNFORMATTED', FILE='RSAVE.'//outna)
          IF (TTEST) WRITE (*, *) ' RSAVE OPENED'
        ELSE
          REWIND (ifile)
          IF (TTEST) WRITE (*, *) ' unit ifile taken as is'
        END IF
      ELSE
        OPEN (UNIT=ifile, STATUS='old', FORM='UNFORMATTED', FILE='SAVE.'//outna)
      END IF
    END IF

#endif

    IF (mynact == 0) THEN
      ! read the iteration WHERE the data has been saved last:
      READ (ifile) iact, nstate_test, nelect, nion, nspdw, trealin
      IF (nstate_test /= nstate_all) &
        STOP ' RESTART: inconsistent nr. of states'
      IF (ttest) WRITE (*, *) 'READ iact etc:', iact, nstate_test, nelect, nion, nspdw, trealin
    END IF

#ifdef COMPLEXSWITCH
    IF (tstatin) THEN
      irest = 0
    ELSE
      irest = iact
    END IF
#endif

! read wavefunctions:

    IF (trealin) ALLOCATE (psiauxr(kdfull2))

    IF (nelect > 0) THEN
    DO nb = 1, nstate
    IF (trealin) THEN
      READ (ifile) occup(nb), psiauxr(1:nxyz)
      psi(1:nxyz, nb) = psiauxr(1:nxyz)
    ELSE
      READ (ifile) occup(nb), psi(1:nxyz, nb)
    END IF
    END DO
    IF (ttest) WRITE (*, *) ' READ occup:', occup(1:nstate)
    READ (ifile) amoy(1:nstate), epotsp(1:nstate), ekinsp(1:nstate)
    IF (ttest) THEN
      WRITE (*, *) ' READ amoy:', amoy(1:nstate)
      WRITE (*, *) ' READ epotsp:', epotsp(1:nstate)
      WRITE (*, *) ' READ ekinsp:', ekinsp(1:nstate)
    END IF
    END IF
    IF (trealin) DEALLOCATE (psiauxr)

    IF (mynact == 0) THEN

! read protonic coordinates and momenta
      IF (nion > 0) THEN
#ifdef REALSWITCH
        READ (ifile) dummy
#else
        IF (tstatin) THEN
          READ (ifile) dummy
        ELSE
          READ (ifile) cx(1:nion), cy(1:nion), cz(1:nion), &
            cpx(1:nion), cpy(1:nion), cpz(1:nion), np(1:nion)
        END IF
#endif
        IF (ttest) THEN
          WRITE (*, *) ' ionic positions/velocities:'
          DO n = 1, nion
            WRITE (*, *) cx(n), cy(n), cz(n), cpx(n), cpy(n), cpz(n)
          END DO
        END IF
      END IF

! read substrate coordinates and momenta

      IF (nelect > 0) THEN
        READ (ifile) qe(1:kmom), se(1:3)
        IF (ttest) WRITE (*, *) ' moments read in:', qe(1:4)
#ifdef COMPLEXSWITCH
        IF (nabsorb > 0) THEN
          IF (tstatin) THEN
            rhoabso = 0D0
          ELSE
            READ (ifile) rhoabso(1:kdfull2)
            IF (ttest) WRITE (*, *) ' rhoabso read in'
          END IF
        END IF
! reading accumulators for laser field
        IF (ttest) WRITE (*, *) ' before laser switch:', tstatin
        IF (.NOT. tstatin) THEN
          IF (ttest) WRITE (*, *) ' before reading laser'
          READ (ifile) ilas, dataold(1:7), datalaser(1:14)
          IF (ttest) WRITE (*, *) 'laser accumulators READ:', dataold(1:7)
        END IF
#endif
      END IF

    END IF




    IF (tstatin) THEN
      CLOSE (UNIT=ifile)
    ELSE
      CLOSE (UNIT=ifile, STATUS='keep')
    END IF


    RETURN

#ifdef REALSWITCH
  END SUBROUTINE resume
#else
END SUBROUTINE restart2
#endif

! **************************

#ifdef REALSWITCH
SUBROUTINE RSAVE(psi, isa, outna)
#else
  SUBROUTINE SAVE(psi, isa, outna)
#endif

! **************************

! writes out the data if MOD(iter,isave)=0
! all the data is saved in the same file called 'SAVE.outna' or 'RSAVE.outna', even in the parallel case

    USE params
    USE kinetic
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
#endif

    INTEGER, INTENT(IN) :: isa
    CHARACTER(LEN=13), INTENT(IN) :: outna
    INTEGER :: iact, mynact, nstate_test, nb
#ifdef COMPLEXSWITCH
    INTEGER :: nbe
#endif
    LOGICAL, PARAMETER :: ttest = .TRUE.
    LOGICAL :: trealin, topenf



    mynact = myn

    IF (ttest) WRITE (*, *) ' SAVE-BEFORE: myn=', myn

#ifdef REALSWITCH

    IF (mynact == 0) THEN
    IF (isaves > 0) THEN
      OPEN (UNIT=ifile, STATUS='unknown', FORM='UNFORMATTED', FILE='RSAVE.'//outna)
      WRITE (*, *) ' RSAVE OPENED'
    ELSE
      OPEN (UNIT=ifile, STATUS='scratch', FORM='UNFORMATTED')
      WRITE (*, *) ' scratch OPENED'
    END IF
    trealin = .true.
    END IF

#else

    IF (mynact == 0) THEN
      IF (isa < 0) THEN
        OPEN (UNIT=ifile, STATUS='unknown', FORM='UNFORMATTED', FILE='RSAVE.'//outna)
        WRITE (*, *) ' RSAVE OPENED for COMPLEX output'
        REWIND (ifile)
      ELSE
        OPEN (UNIT=ifile, STATUS='unknown', FORM='UNFORMATTED', FILE='SAVE.'//outna)
        WRITE (*, *) ' SAVE OPENED for COMPLEX output'
        REWIND (ifile)
      END IF
      trealin = .false.
    END IF

#endif

! WRITE iteration at which the data is saved
    IF (mynact == 0) THEN
      WRITE (ifile) isa, nstate_all, nelect, nion, nspdw, trealin
      IF (ttest) WRITE (*, *) 'WROTE iact etc:', iact, nstate_test, nelect, nion, nspdw, trealin
    END IF

! WRITE wavefunctions:
    IF (nelect > 0) THEN
      DO nb = 1, nstate
        WRITE (ifile) occup(nb), psi(1:nxyz, nb)
      END DO
      WRITE (ifile) amoy(1:nstate), epotsp(1:nstate), ekinsp(1:nstate)
      IF (TTEST) WRITE (6, *) ' SAVE: wavefunctions written'
    END IF


    IF (mynact == 0) THEN

! WRITE protonic coordinates and momenta:
      IF (nion > 0) &
        WRITE (ifile) cx(1:nion), cy(1:nion), cz(1:nion), &
        cpx(1:nion), cpy(1:nion), cpz(1:nion), np(1:nion)
      IF (TTEST) WRITE (*, *) 'ionic coordinates written'
! WRITE substrate coordinates and momenta

! WRITE dipole moment etc:
      IF (nelect > 0) THEN
        WRITE (ifile) qe(1:kmom), se(1:3)
        IF (TTEST) WRITE (*, *) 'electronic moments written'

#ifdef COMPLEXSWITCH
        IF (isa > 0 .AND. nabsorb > 0) THEN
          WRITE (ifile) rhoabso(1:kdfull2)
          WRITE (*, *) ' RHOABSO written'
        END IF
! writing cumulators for laser field
        IF (isa .GE. 0) THEN
          WRITE (ifile) ilas, dataold(1:7), datalaser(1:14)
          WRITE (*, *) 'laser cumulators written'
        END IF
#endif
      END IF
    END IF


    IF (mynact == 0 .AND. (isaves > 0 .OR. isaved > 0)) THEN
      INQUIRE (ifile, OPENED=topenf)
      IF(topenf) CLOSE (UNIT=ifile, STATUS='keep')
    END IF


    RETURN
#ifdef REALSWITCH
  END SUBROUTINE RSAVE
#else
END SUBROUTINE SAVE
#endif

#ifdef REALSWITCH

! **************************

SUBROUTINE restherm()

! **************************

  USE params
  USE kinetic
  IMPLICIT NONE
  INTEGER :: iact, ion

  OPEN (UNIT=20, STATUS='unknown', FORM='UNFORMATTED', FILE='therm')

  READ (20) iact
  IF (iact /= irest) THEN
    WRITE (7, *) 'iact=', iact
    STOP 'bad irest in for005dyn !'
  END IF
  DO ion = 1, nion
    READ (20) cx(ion), cy(ion), cz(ion)
    READ (20) cpx(ion), cpy(ion), cpz(ion)
  END DO
  CLOSE (UNIT=20, STATUS='keep')
  RETURN
END SUBROUTINE restherm

! **************************

SUBROUTINE addcluster(psi, outna)

! **************************

  USE params
  USE kinetic
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
  CHARACTER(LEN=13), INTENT(IN OUT) :: outna

  INTEGER :: i, iact, idum, ii, ion, k, nb, nelectt, niont, nspdwt, nstatet
  OPEN (UNIT=ifile, STATUS='unknown', FORM='UNFORMATTED', FILE='SAVE.'//outna)

! read the iteration where the data has been saved last:

  READ (ifile) iact, nstatet, nelectt, niont, nspdwt
  DO i = 1, nstatet
    READ (ifile) occup(i + nstate)
  END DO

  irest = iact

! read wavefunctions:
  IF (nelectt > 0) THEN
  DO nb = 1, nstatet
  DO i = 1, nxyz
    READ (ifile) psi(i, nb + nstate)
  END DO
  END DO
  END IF

  DO i = 1, ksttot
    IF (i <= nstatet) THEN
      READ (ifile) ispin(i + nstate), nrel2abs(i + nstate), nabs2rel(i + nstate)
    ELSE
      READ (ifile) idum, idum, idum
    END IF
  END DO

! READ protonic coordinates and momenta
  WRITE (6, *) nion, niont
  DO ion = 1, niont
    READ (ifile) cx(ion + nion), cy(ion + nion), cz(ion + nion), np(ion + nion)
    READ (ifile) cpx(ion + nion), cpy(ion + nion), cpz(ion + nion), np(ion + nion)
  END DO


  IF (nelectt > 0) THEN
  DO k = 1, kmom
    READ (ifile) qe(k)
  END DO
  DO i = 1, 3
    READ (ifile) se(i)
  END DO
  END IF

  IF (nabsorb > 0) THEN
    DO ii = 1, kdfull2
      READ (ifile) rhoabso(ii)
    END DO
  END IF

  CLOSE (UNIT=ifile, STATUS='keep')

  nstate = nstate + nstatet
  nelect = nelect + nelectt
  nion = nion + niont
  nspdw = nspdw + nspdwt
  nion2 = nion

  RETURN
END SUBROUTINE addcluster

#endif

