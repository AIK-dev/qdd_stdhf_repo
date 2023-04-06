
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
#if(twostsic)
    USE twostr, ONLY: vecsr, ndims
    USE twost, ONLY: vecs, expdabold, wfrotate
#endif
!#endif
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
    LOGICAL, PARAMETER :: tstatin = .false.
#else
    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
#if(mpi)
    REAL(DP), ALLOCATABLE :: rhoabsoorb_all(:, :)
#endif
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

#if(mpi)
    INCLUDE 'mpif.h'
    INTEGER :: is(mpi_status_size)
    INTEGER :: i, nba, nod, nstnod
    REAL(DP):: occupact
#ifdef REALSWITCH
    REAL(DP), ALLOCATABLE :: psiaux(:)
#else
    COMPLEX(DP), ALLOCATABLE :: psiaux(:)
#endif
    REAL(DP), DIMENSION(:), ALLOCATABLE :: amoya, epotspa, ekinspa
    REAL(DP), DIMENSION(:), ALLOCATABLE :: amoys, epotsps, ekinsps
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: nrel2absf

    CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)

    ALLOCATE (psiaux(kdfull2))
    ALLOCATE (amoya(ksttot), epotspa(ksttot), ekinspa(ksttot))
    ALLOCATE (amoys(kstate), epotsps(kstate), ekinsps(kstate))
    ALLOCATE (nrel2absf(kstate, 0:knode - 1))

    DO nod = 0, knode - 1
      DO i = 1, kstate
        nrel2absf(i, nod) = i + nstart_node(nod)
      END DO
    END DO

#endif

#if(mpi)
    mynact = 0
#else
    mynact = myn
#endif

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
#if(mpi)
    IF (knode > 1) THEN
      CALL mpi_bcast(iact, 1, mpi_integer, 0, mpi_comm_world, mpi_ierror)
      CALL mpi_bcast(nelect, 1, mpi_integer, 0, mpi_comm_world, mpi_ierror)
      CALL mpi_bcast(nion, 1, mpi_integer, 0, mpi_comm_world, mpi_ierror)
      CALL mpi_bcast(nspdw, 1, mpi_integer, 0, mpi_comm_world, mpi_ierror)
    END IF
#endif

#ifdef COMPLEXSWITCH
    IF (tstatin) THEN
      irest = 0
    ELSE
      irest = iact
    END IF
#endif

! read wavefunctions:

    IF (trealin) ALLOCATE (psiauxr(kdfull2))

#if(nompi)
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
#endif
#if(mpi)
    CALL mpi_barrier(mpi_comm_world, mpi_ierror)

    IF (nelect > 0) THEN
    DO nb = 1, nstate_all
      nod = nhome(nb)
      nba = nabs2rel(nb)
      IF (myn == 0) THEN
      IF (trealin) THEN
        READ (ifile) occupact, psiauxr(1:nxyz)
        psiaux = psiauxr
      ELSE
        READ (ifile) occupact, psiaux(1:nxyz)
      END IF
      IF (ttest) WRITE (*, *) ' READ: psiaux at myn,nb,occup=', myn, nb, occupact
      END IF
      CALL send_and_receive(occupact, occup(nba), 1, 0, nod)
#ifdef REALSWITCH
        CALL send_and_receive(psiaux, psi(1, nba), kdfull2, 0, nod)
#else
        CALL csend_and_receive(psiaux, psi(1, nba), kdfull2, 0, nod)
#endif
        IF (ttest .AND. myn == nod) WRITE (*, *) ' SENT to node,nb=', nod, nb
      END DO

      IF (myn == 0) READ (ifile) amoya(1:nstate_all), epotspa(1:nstate_all), &
        ekinspa(1:nstate_all)

      CALL mpi_barrier(mpi_comm_world, mpi_ierror)
      DO nod = 0, knode - 1
        CALL mpi_barrier(mpi_comm_world, mpi_ierror)
        IF (nod == myn) WRITE (*, *) 'myn ', myn, ' occup :', occup(:nstate)
        CALL flush (6)
      END DO

      DO nod = 0, knode - 1
        nstnod = nstate_node(nod)
        DO nb = 1, nstnod
          nba = nrel2absf(nb, nod)
          amoys(nb) = amoya(nba)
          epotsps(nb) = epotspa(nba)
          ekinsps(nb) = ekinspa(nba)
        END DO

        CALL send_and_receive(amoys, amoy, nstnod, 0, nod)
        CALL send_and_receive(epotsps, epotsp, nstnod, 0, nod)
        CALL send_and_receive(ekinsps, ekinsp, nstnod, 0, nod)
      END DO
    END IF

#endif
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
#if(raregas)
      IF (isurf /= 0) THEN
        READ (ifile) imobc(1:nc), xc(1:nc), yc(1:nc), zc(1:nc), &
          pxc(1:nc), pyc(1:nc), pzc(1:nc)
        READ (ifile) imobe(1:ne), xe(1:ne), ye(1:ne), ze(1:ne), &
          pxe(1:ne), pye(1:ne), pze(1:ne)
        READ (ifile) imobk(1:nk), xk(1:nk), yk(1:nk), zk(1:nk), &
          pxk(1:nk), pyk(1:nk), pzk(1:nk)
        READ (ifile) potfixedion(1:kdfull2)
        IF (ttest) WRITE (*, *) ' surface read in. nc,nk,ne=', nc, nk, ne
      END IF
#endif

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
#if(extended)
          IF (jescmaskorb /= 0) THEN
#if(mpi)
            ALLOCATE (rhoabsoorb_all(kdfull2, nstate_all))
            IF (tstatin) THEN
              rhoabsoorb_all = 0D0
            ELSE
              DO nbe = 1, nstate_all
                READ (ifile) rhoabsoorb_all(1:kdfull2, nbe)
              END DO
            END IF
#else
            IF (tstatin) THEN
              rhoabsoorb = 0D0
            ELSE
              DO nbe = 1, nstate
                READ (ifile) rhoabsoorb(1:kdfull2, nbe)
                IF (ttest) WRITE (*, *) ' rhoabsoorb read in'
              END DO
            END IF
#endif
          END IF
#endif
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

#if(mpi)
    IF (knode > 1) THEN
    IF (nelect > 0) THEN
      CALL mpi_bcast(qe(1:kmom), kmom, mpi_double_precision, 0, mpi_comm_world, mpi_ierror)
      CALL mpi_bcast(se(1:3), 3, mpi_double_precision, 0, mpi_comm_world, mpi_ierror)
    END IF
    IF (nabsorb > 0) &
      CALL mpi_bcast(rhoabso(1:kdfull2), kdfull2, &
                     mpi_double_precision, 0, mpi_comm_world, mpi_ierror)
    CALL mpi_bcast(dataold(1:7), 7, mpi_double_precision, 0, mpi_comm_world, mpi_ierror)
    CALL mpi_bcast(datalaser(1:14), 11, mpi_double_precision, 0, mpi_comm_world, mpi_ierror)
    CALL mpi_bcast(ilas, 1, mpi_integer, 0, mpi_comm_world, mpi_ierror)
    END IF

#ifdef COMPLEXSWITCH
    IF (nabsorb > 0 .AND. jescmaskorb /= 0) THEN
      DO nb = 1, nstate_all
        nod = nhome(nb)
        nba = nabs2rel(nb)
        CALL send_and_receive(rhoabsoorb_all(1:kdfull2, nb), rhoabsoorb(1:kdfull2, nba), kdfull2, 0, nod)
      END DO
      IF (myn == 0) DEALLOCATE (rhoabsoorb_all)
    END IF
#endif
#endif

#if(twostsic)
    IF (ifsicp >= 6) THEN
#ifdef REALSWITCH
      READ (ifile) vecsr(1:kstate, 1:kstate, 1:2), ndims(1:2)
      WRITE (*, *) ' READ vecsr:'
      DO n = 1, ndims(1)
        WRITE (*, '(8f10.6)') vecsr(1:ndims(1), n, 1)
      END DO
#endif
#ifdef COMPLEXSWITCH
      WRITE (*, *) ' before reading VECS'
      READ (ifile) vecs(1:kstate, 1:kstate, 1:2), ndims(1:2)
      WRITE (*, *) vecs(1, 1, 1)
      IF (.NOT. tstatin) THEN
        READ (ifile) ExpDABold(1:kstate, 1:kstate, 1:2)
        READ (ifile) wfrotate(1:kstate, 1:kstate, 1:2)
      END IF
      WRITE (*, *) ' READ vecs:', ndims
      DO n = 1, ndims(1)
        WRITE (*, '(8f10.6)') vecs(1:ndims(1), n, 1)
      END DO
      DO n = 1, ndims(2)
        WRITE (*, '(5(2f10.5,2x))') vecs(1:ndims(2), n, 2)
      END DO
#endif
    END IF
#endif

#if(extended)
#ifdef COMPLEXSWITCH
    IF (mynact == 0) THEN
      IF (.NOT. tstatin .AND. jattach /= 0) THEN
        WRITE (*, *) 'READ totintegprob: tstat,irest=', tstat, irest
        READ (ifile) totintegprob
        WRITE (*, *) totintegprob
        READ (ifile) reference_energy
        WRITE (*, *) reference_energy
      END IF
    END IF
#endif
#endif

    IF (tstatin) THEN
      CLOSE (UNIT=ifile)
    ELSE
      CLOSE (UNIT=ifile, STATUS='keep')
    END IF

#if(mpi)
    CALL comm_ionconfig()
    DEALLOCATE (psiaux)
    DEALLOCATE (amoya, epotspa, ekinspa)
    DEALLOCATE (amoys, epotsps, ekinsps)
    DEALLOCATE (nrel2absf)
#endif

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
#if(twostsic)
    USE twostr, ONLY: vecsr, ndims
    USE twost, ONLY: vecs, expdabold, wfrotate
#endif
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
#if(mpi)
    REAL(DP), ALLOCATABLE :: rhoabsoorb_all(:, :)
#endif
#endif

    INTEGER, INTENT(IN) :: isa
    CHARACTER(LEN=13), INTENT(IN) :: outna
    INTEGER :: iact, mynact, nstate_test, nb
#ifdef COMPLEXSWITCH
    INTEGER :: nbe
#endif
    LOGICAL, PARAMETER :: ttest = .TRUE.
    LOGICAL :: trealin, topenf

#if(twostsic)
    INTEGER :: n
#endif

#if(mpi)
    INCLUDE 'mpif.h'
    INTEGER :: is(mpi_status_size)
    INTEGER :: i, nba, nod, nstnod
#ifdef REALSWITCH
    REAL(DP), ALLOCATABLE :: psiaux(:)

#else
    COMPLEX(DP), ALLOCATABLE :: psiaux(:)
#endif

    REAL(DP) :: occupact(1)
    REAL(DP), DIMENSION(:), ALLOCATABLE :: amoya, epotspa, ekinspa
    REAL(DP), DIMENSION(:), ALLOCATABLE :: amoys, epotsps, ekinsps
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: nrel2absf

    CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)

    ALLOCATE (psiaux(kdfull2))
    ALLOCATE (amoya(ksttot), epotspa(ksttot), ekinspa(ksttot))
    ALLOCATE (amoys(kstate), epotsps(kstate), ekinsps(kstate))
    ALLOCATE (nrel2absf(kstate, 0:knode - 1))

    DO nod = 0, knode - 1
      DO i = 1, kstate
        nrel2absf(i, nod) = i + nstart_node(nod)
      END DO
    END DO

    CALL mpi_barrier(mpi_comm_world, mpi_ierror)

#endif

#if(mpi)
    mynact = 0
#else
    mynact = myn
#endif

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
#if(nompi)
    IF (nelect > 0) THEN
      DO nb = 1, nstate
        WRITE (ifile) occup(nb), psi(1:nxyz, nb)
      END DO
      WRITE (ifile) amoy(1:nstate), epotsp(1:nstate), ekinsp(1:nstate)
      IF (TTEST) WRITE (6, *) ' SAVE: wavefunctions written'
    END IF
#endif
#if(mpi)

    CALL mpi_barrier(mpi_comm_world, mpi_ierror)

    IF (nelect > 0) THEN
    DO nb = 1, nstate_all
      nod = nhome(nb)
      nba = nabs2rel(nb)
      IF (ttest .AND. myn == nod) WRITE (*, *) ' SENT from node,nb=', nod, nb
      CALL send_and_receive(occup(nba), occupact, 1, nod, 0)
#ifdef REALSWITCH
        CALL send_and_receive(psi(1, nba), psiaux, kdfull2, nod, 0)
#else
        CALL csend_and_receive(psi(1, nba), psiaux, kdfull2, nod, 0)
#endif
        IF (myn == 0) THEN
          WRITE (ifile) occupact, psiaux(1:nxyz)
          IF (ttest) WRITE (*, *) ' RECEIVED and written: node,nb=', nod, nb
        END IF
      END DO

      CALL mpi_barrier(mpi_comm_world, mpi_ierror)

      WRITE (*, *) ' before send ekins: myn=', myn
      DO nod = 0, knode - 1
        nstnod = nstate_node(nod)
        CALL send_and_receive(amoy, amoys, nstnod, nod, 0)
        CALL send_and_receive(epotsp, epotsps, nstnod, nod, 0)
        CALL send_and_receive(ekinsp, ekinsps, nstnod, nod, 0)
        DO nb = 1, nstnod
          nba = nrel2absf(nb, nod)
          amoya(nba) = amoys(nb)
          epotspa(nba) = epotsps(nb)
          ekinspa(nba) = ekinsps(nb)
        END DO
      END DO
      IF (myn == 0) WRITE (ifile) amoya(1:nstate_all), epotspa(1:nstate_all), &
        ekinspa(1:nstate_all)
      WRITE (*, *) ' after send ekins: myn=', myn

    END IF
#endif

#if(mpi)
#ifdef COMPLEXSWITCH
    IF (nabsorb > 0 .AND. jescmaskorb /= 0) THEN
      IF (myn == 0) ALLOCATE (rhoabsoorb_all(kdfull2, nstate_all))
      DO nb = 1, nstate_all
        nod = nhome(nb)
        nba = nabs2rel(nb)
        CALL send_and_receive(rhoabsoorb(1:kdfull2, nba), rhoabsoorb_all(1:kdfull2, nb), kdfull2, nod, 0)
      END DO
    END IF
#endif
#endif

    IF (mynact == 0) THEN

! WRITE protonic coordinates and momenta:
      IF (nion > 0) &
        WRITE (ifile) cx(1:nion), cy(1:nion), cz(1:nion), &
        cpx(1:nion), cpy(1:nion), cpz(1:nion), np(1:nion)
      IF (TTEST) WRITE (*, *) 'ionic coordinates written'
! WRITE substrate coordinates and momenta
#if(raregas)
      IF (isurf /= 0) THEN
        WRITE (ifile) imobc(1:nc), xc(1:nc), yc(1:nc), zc(1:nc), &
          pxc(1:nc), pyc(1:nc), pzc(1:nc)
        WRITE (ifile) imobe(1:ne), xe(1:ne), ye(1:ne), ze(1:ne), &
          pxe(1:ne), pye(1:ne), pze(1:ne)
        WRITE (ifile) imobk(1:nk), xk(1:nk), yk(1:nk), zk(1:nk), &
          pxk(1:nk), pyk(1:nk), pzk(1:nk)
        WRITE (ifile) potfixedion(1:kdfull2)
      END IF
#endif

! WRITE dipole moment etc:
      IF (nelect > 0) THEN
        WRITE (ifile) qe(1:kmom), se(1:3)
        IF (TTEST) WRITE (*, *) 'electronic moments written'

#ifdef COMPLEXSWITCH
        IF (isa > 0 .AND. nabsorb > 0) THEN
          WRITE (ifile) rhoabso(1:kdfull2)
#if(extended)
          IF (jescmaskorb /= 0) THEN
#if(mpi)
            DO nbe = 1, nstate_all
              WRITE (ifile) rhoabsoorb_all(1:kdfull2, nbe)
            END DO
#else
            DO nbe = 1, nstate
              WRITE (ifile) rhoabsoorb(1:kdfull2, nbe)
            END DO
#endif
          END IF
#endif
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

#if(twostsic)
#ifdef REALSWITCH
    IF (ifsicp >= 6) THEN
      WRITE (ifile) vecsr(1:kstate, 1:kstate, 1:2), ndims(1:2)
      WRITE (*, *) 'vecsr written'
      DO n = 1, ndims(1)
        WRITE (*, '(10f10.5)') vecsr(1:ndims(1), n, 1)
      END DO
      DO n = 1, ndims(2)
        WRITE (*, '(10f10.5)') vecsr(1:ndims(2), n, 2)
      END DO
    END IF
#endif
#ifdef COMPLEXSWITCH
    IF (ifsicp >= 6) THEN
      WRITE (ifile) vecs(1:kstate, 1:kstate, 1:2), ndims(1:2)
      WRITE (ifile) expdabold(1:kstate, 1:kstate, 1:2)
      WRITE (ifile) wfrotate(1:kstate, 1:kstate, 1:2)
      WRITE (*, *) 'vecs and expdabold written'
      WRITE (*, *) ndims
      DO n = 1, ndims(1)
        WRITE (*, '(5(2f10.5,2x))') vecs(1:ndims(1), n, 1)
      END DO
      DO n = 1, ndims(2)
        WRITE (*, '(5(2f10.5,2x))') vecs(1:ndims(2), n, 2)
      END DO
    END IF
#endif
#endif

#if(extended)
#ifdef COMPLEXSWITCH
    IF (jattach /= 0) THEN
      WRITE (ifile) totintegprob
      WRITE (ifile) reference_energy
    END IF
#endif
#endif
    IF (mynact == 0 .AND. (isaves > 0 .OR. isaved > 0)) THEN
      INQUIRE (ifile, OPENED=topenf)
      IF(topenf) CLOSE (UNIT=ifile, STATUS='keep')
    END IF

#if(mpi)
    DEALLOCATE (psiaux)
    DEALLOCATE (amoya, epotspa, ekinspa)
    DEALLOCATE (amoys, epotsps, ekinsps)
    DEALLOCATE (nrel2absf)
#endif

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

#if(raregas)
  IF (isurf /= 0) THEN

    DO i = 1, nc
      READ (ifile) imobc(i)
      READ (ifile) xc(i), yc(i), zc(i)
      READ (ifile) pxc(i), pyc(i), pzc(i)
    END DO
    DO i = 1, NE
      READ (ifile) imobe(i)
      READ (ifile) xe(i), ye(i), ze(i)
      READ (ifile) pxe(i), pye(i), pze(i)
    END DO
    DO i = 1, nk
      READ (ifile) imobk(i)
      READ (ifile) xk(i), yk(i), zk(i)
      READ (ifile) pxk(i), pyk(i), pzk(i)
    END DO

    DO i = 1, kdfull2
      READ (ifile) potfixedion(i)
    END DO

  END IF
#endif

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

#if(mpi)
#ifdef REALSWITCH
SUBROUTINE send_and_receive(instring, outstring, length, in_node, dest_node)

! Sends 'instring' from node 'in_node' TO
! 'outstring' on node 'out_node'.
! Both strings are DOUBLE PRECISION and have length 'length'

  USE params
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  REAL(DP) :: instring(length), outstring(length)
  INTEGER :: length, in_node, dest_node

  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)

  IF (dest_node == in_node .AND. myn == dest_node) THEN
    outstring = instring
    RETURN
  END IF

  IF (myn == dest_node) &
    CALL mpi_recv(outstring, length, mpi_double_precision, in_node, &
                  mpi_any_tag, mpi_comm_world, is, mpi_ierror)

  IF (myn == in_node) &
    CALL mpi_send(instring, length, mpi_double_precision, dest_node, 1, &
                  mpi_comm_world, mpi_ierror)

  RETURN

END SUBROUTINE send_and_receive

#else

SUBROUTINE csend_and_receive(instring, outstring, length, in_node, dest_node)

! Sends 'instring' from node 'in_node' to
! 'outstring' on node 'out_node'.
! Both strings are DOUBLE COMPLEX and have length 'length'

  USE params
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  COMPLEX(DP) :: instring(length), outstring(length)
  INTEGER :: length, in_node, dest_node

  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)

  IF (dest_node == in_node .AND. myn == dest_node) THEN
    outstring = instring
    RETURN
  END IF

  IF (myn == dest_node) &
    CALL mpi_recv(outstring, length, mpi_double_complex, in_node, &
                  mpi_any_tag, mpi_comm_world, is, mpi_ierror)

  IF (myn == in_node) &
    CALL mpi_send(instring, length, mpi_double_complex, dest_node, 1, &
                  mpi_comm_world, mpi_ierror)

  RETURN

END SUBROUTINE csend_and_receive
#endif
#endif
