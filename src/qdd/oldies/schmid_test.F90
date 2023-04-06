
#if(nompi)

! ******************************

SUBROUTINE schmidt(q0)

! ******************************

! Serial and OpenMP version of Schmidt orthogonalization.
!
! Input/Output:
! q0 = set of REAL s.p. wavefunctions TO be ortho-normalized

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)

  INTEGER :: i, imin, isav, n, nbe, nbes, ncc, ncs
  REAL(DP) :: cs, emin
  REAL(DP) :: eord(kstate)
  INTEGER :: isort(kstate)

  LOGICAL, PARAMETER :: tord = .true.

!*********************************************************

! sort the s.p. energies

  DO n = 1, nstate
    eord(n) = amoy(n)
    isort(n) = n
  END DO
  IF (tord) THEN
    DO n = 1, nstate
      emin = 1D32
      imin = n
      DO i = n, nstate
        IF (eord(i) < emin) THEN
          emin = eord(i)
          imin = i
        END IF
      END DO
      isav = isort(imin)
      eord(imin) = eord(n)
      isort(imin) = isort(n)
      eord(n) = emin
      isort(n) = isav
    END DO
  END IF

! Schmidt ortho-normalisation

#if(omp)

  DO nbes = 1, nstate
    nbe = isort(nbes)
    q0(1:nxyz, nbe) = q0(1:nxyz, nbe)/SQRT(SUM(q0(1:nxyz, nbe)**2)*dvol)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ncs,ncc,cs) SCHEDULE(STATIC)
    DO ncs = nbes + 1, nstate
      ncc = isort(ncs)
      IF ((ispin(nbe) == ispin(ncc))) THEN
        cs = SUM(q0(1:nxyz, nbe)*q0(1:nxyz, ncc))*dvol
        q0(1:nxyz, ncc) = q0(1:nxyz, ncc) - cs*q0(1:nxyz, nbe)
      END IF
    END DO
!$OMP END PARALLEL DO

  END DO

#else

  DO nbes = 1, nstate
    nbe = isort(nbes)

    DO ncs = 1, nbes - 1
      ncc = isort(ncs)
      IF ((ispin(nbe) == ispin(ncc))) THEN
        cs = SUM(q0(1:nxyz, nbe)*q0(1:nxyz, ncc))*dvol
        q0(1:nxyz, nbe) = q0(1:nxyz, nbe) - cs*q0(1:nxyz, ncc)
      END IF
    END DO

    q0(1:nxyz, nbe) = q0(1:nxyz, nbe)/SQRT(SUM(q0(1:nxyz, nbe)**2)*dvol)

  END DO

#endif

  RETURN
END SUBROUTINE schmidt

#endif

#if(mpi)

! ******************************

SUBROUTINE schmidt(q0)

! ******************************

! MPI parallele version of Schmidt orthogonalization.
!
! Input/Output:
! q0 = set of REAL s.p. wavefunctions TO be ortho-normalized

  USE params
  IMPLICIT NONE
  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)

  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  REAL(DP), ALLOCATABLE :: q3(:) ! workspace for foreign w.f.
  INTEGER :: iend, isp3, nbe1, nbe2, nod1, nod2

! fields for some RECURSIVE algorithm, this is a dummy
! operation which serves TO inhibit loop unrolling

  LOGICAL :: tsync
  LOGICAL, PARAMETER :: ttest = .false.

!*********************************************************

  ALLOCATE (q3(kdfull2))

! determine own node 'myn'

  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)

! loop over nodes and states: w.f. TO be normalized

! orthogonalize IN each next node with wfs sent from previous nodes

  IF (ttest) WRITE (*, *) 'SCHMID: myn=', myn
  DO nod2 = 0, myn - 1
    DO nbe2 = 1, nstate_node(nod2)
      IF (ttest) THEN
        WRITE (6, '(a,5i5)') 'SCHMID: wait myn,nbe2,nod2,tag=', &
          myn, nbe2, nod2, 2*nbe2 + 1
        WRITE (7, '(a,5i5)') 'SCHMID: wait myn,nbe2,nod2,tag=', &
          myn, nbe2, nod2, 2*nbe2 + 1
        CALL flush ()
      END IF
      CALL mpi_recv(q3, kdfull2, mpi_double_precision, nod2, nbe2*2 + 1, &
                    mpi_comm_world, is, mpi_ierror)
      isp3 = ispin_node(nbe2, nod2)
      IF (ttest) THEN
        WRITE (6, '(a,5i5)') 'SCHMID: received myn,nbe2,nod2,tag=', &
          myn, nbe2, nod2, nbe2*2 + 1
        WRITE (7, *) 'SCHMID: received myn,nbe2,nod2=', myn, nbe2, nod2
        CALL flush ()
      END IF
      DO nbe1 = 1, nstate_node(myn)
        IF (ispin(nrel2abs(nbe1)) == isp3) THEN
          CALL orthogonalize(q0(1, nbe1), q3)
          IF (ttest) THEN
            WRITE (7, '(a,4i5,3i3)') &
              ' orthog: myn,nbe1,nod2,nbe2,spin2,spin1=', &
              myn, nbe1, nod2, nbe2, isp3, ispin(nrel2abs(nbe1)), ispin_node(nbe1, myn)
            WRITE (6, '(a,4i5,3i3)') &
              ' orthog: myn,nbe1,nod2,nbe2,spin2,spin1=', &
              myn, nbe1, nod2, nbe2, isp3, ispin(nrel2abs(nbe1)), ispin_node(nbe1, myn)
          END IF
        END IF
      END DO
    END DO
  END DO

! ortho-normalize within node assuming orthogonality on
! previous nodes (established above)

  DO nbe1 = 1, nstate_node(myn)
    DO nbe2 = 1, nbe1 - 1
      IF (ispin(nrel2abs(nbe1)) == ispin(nrel2abs(nbe2))) THEN
        CALL orthogonalize(q0(1, nbe1), q0(1, nbe2))
        IF (ttest) WRITE (*, *) 'SCHMID orth. myn,nbe1,nbe2=', myn, nbe1, nbe2
        IF (ttest) WRITE (7, *) 'SCHMID orth. myn,nbe1,nbe2=', myn, nbe1, nbe2
      END IF
    END DO
    CALL normalize(q0(1, nbe1))
    IF (ttest) WRITE (*, *) 'SCHMID: norm. myn,nbe1=', myn, nbe1
    IF (ttest) WRITE (7, *) 'SCHMID: norm. myn,nbe1=', myn, nbe1

! distribute TO next nodes

    DO nod2 = myn + 1, knode - 1
      IF (ttest) THEN
        WRITE (6, '(a,5i5)') 'SCHMID: before send myn,nbe1,nod2,tag=', &
          myn, nbe1, nod2, nbe1*2 + 1
        WRITE (7, '(a,5i5)') 'SCHMID: before send myn,nbe1,nod2,tag=', &
          myn, nbe1, nod2, nbe1*2 + 1
        CALL flush ()
      END IF
      CALL mpi_ssend(q0(1, nbe1), kdfull2, mpi_double_precision, &
                     nod2, nbe1*2 + 1, mpi_comm_world, mpi_ierror)
      IF (ttest) THEN
        WRITE (6, '(a,5i5)') 'SCHMID: sent myn,nbe1,nod2,tag=', &
          myn, nbe1, nod2, nbe1*2 + 1
        WRITE (7, *) 'SCHMID: sent myn,nbe1,nod2=', myn, nbe1, nod2
        CALL flush ()
      END IF
    END DO

  END DO

! synchronize by sending and receiving termination signal
! from last node

  tsync = .false.
  IF (tsync) THEN
    IF (myn == knode - 1) THEN
      iend = 1
      DO nod2 = 0, knode - 2
        CALL mpi_bsend(iend, 1, mpi_integer, nod2, 1, mpi_comm_world, mpi_ierror)
        IF (ttest) WRITE (*, *) ' terminator sent myn=', myn
      END DO
    ELSE
      CALL mpi_recv(iend, 1, mpi_integer, knode - 1, &
                    mpi_any_tag, mpi_comm_world, is, mpi_ierror)
      IF (ttest) WRITE (*, *) ' terminator received myn,iend=', myn, iend
    END IF
  END IF

  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
  DEALLOCATE (q3)

  RETURN
END SUBROUTINE schmidt

#endif

! ******************************

SUBROUTINE normalize(qact)

! ******************************

! normalizes REAL wavefunction on 'qact'

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: qact(kdfull2)

  INTEGER :: i
  REAL(DP) :: cs

!*********************************************************

! compute norm

  cs = 0D0
  DO i = 1, nxyz
    cs = cs + qact(i)*qact(i)
  END DO
  cs = cs*dvol

! normalize

  cs = 1D0/SQRT(cs)
  DO i = 1, nxyz
    qact(i) = qact(i)*cs
  END DO

  RETURN

END SUBROUTINE normalize

! ******************************

SUBROUTINE orthogonalize(qact, qorth)

! ******************************

! orthogonalizes REAL wavefunction 'qact' on 'qorth'

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: qact(kdfull2)
  REAL(DP), INTENT(IN) :: qorth(kdfull2)

  INTEGER :: i
  REAL(DP) :: cs

!*********************************************************

! compute overlap

  cs = 0D0
  DO i = 1, nxyz
    cs = cs + qact(i)*qorth(i)
  END DO
  cs = cs*dvol

! orthogonalize

  DO i = 1, nxyz
    qact(i) = qact(i) - cs*qorth(i)
  END DO

  RETURN

END SUBROUTINE orthogonalize

! ******************************

SUBROUTINE cschmidt(q0)

! ******************************

! Serial version of Schmidt orthogonalization for COMPLEX wfs.
!
! Input/Output:
! q0 = set of COMPLEX s.p. wavefunctions TO be ortho-normalized

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)

  REAL(DP) :: eord(kstate)
  INTEGER :: isort(kstate)

  INTEGER :: i, imin, isav, n, nbe, nbes, ncc, ncs
  REAL(DP) :: emin
  COMPLEX(DP) :: cs

  LOGICAL, PARAMETER :: tord = .false.

!*********************************************************

! sort the s.p. energies

  DO n = 1, nstate
    eord(n) = amoy(n)
    isort(n) = n
  END DO
  IF (tord) THEN
    DO n = 1, nstate
      emin = 1D32
      imin = n
      DO i = n, nstate
        IF (eord(i) < emin) THEN
          emin = eord(i)
          imin = i
        END IF
      END DO
      isav = isort(imin)
      eord(imin) = eord(n)
      isort(imin) = isort(n)
      eord(n) = emin
      isort(n) = isav
    END DO
  END IF

! Schmidt ortho-normalisation

  DO nbes = 1, nstate
    nbe = isort(nbes)

    DO ncs = 1, nstate
      ncc = isort(ncs)
      IF ((ispin(nbe) == ispin(ncc)) .AND. ncc <= nbe) THEN
        cs = CMPLX(0D0, 0D0, DP)
        DO i = 1, nxyz
          cs = cs + CONJG(q0(i, ncc))*q0(i, nbe)
        END DO
        cs = cs*dvol
        IF (ncc == nbe) THEN
          DO i = 1, nxyz
            q0(i, nbe) = q0(i, nbe)/SQRT(REAL(cs, DP))
          END DO
        ELSE
          DO i = 1, nxyz
            q0(i, nbe) = q0(i, nbe) - cs*q0(i, ncc)
          END DO
        END IF
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE cschmidt

SUBROUTINE Loewdin(q0, aloc, tdiag)
! ******************************

! Serial version of orthogonalization with Loewdin scheme
!
! Input/Output:
! q0 = set of REAL s.p. wavefunctions TO be ortho-normalized
! tdiag = switch TO Hamiltonian diagonalization

  USE params
  USE util, ONLY: wfovlp
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  LOGICAL, INTENT(IN) :: tdiag

  INTEGER :: nbe, nbes
  REAL(DP) :: cs, emin
  REAL(DP) :: eord(kstate)
  INTEGER :: isort(kstate)

  REAL(DP), ALLOCATABLE :: omatr(:) ! tridiagonal storage
  REAL(DP), ALLOCATABLE :: oeigen(:)
  REAL(DP), ALLOCATABLE :: vect(:, :), vect2(:, :), hmatr(:, :)
  REAL(DP), ALLOCATABLE :: psistate(:), hpsi(:)
  INTEGER, ALLOCATABLE :: npoi(:, :)
  INTEGER :: ishift, ispact, ii, nbc, nbcs, ntridig, nlower, nupper, nstsp(2), ktridig

  LOGICAL, PARAMETER :: tprint = .TRUE.

!*********************************************************

  IF (tprint) WRITE (*, *) 'Loewdin entered: nstate=', nstate

  CALL sortwf_energ(q0)

! check that states comes IN contiguous blocks of spin
  ALLOCATE (npoi(nstate, 2))
  ispact = 1
  nstsp = 0
  npoi = 0
  DO nbe = 1, nstate
    IF (ispin(nbe) == ispact) THEN
      nstsp(ispact) = 1 + nstsp(ispact)
      npoi(nstsp(ispact), ispact) = nbe
    ELSE
      IF (ispact == 1) THEN
        ispact = 2
        nstsp(ispact) = 1 + nstsp(ispact)
        npoi(nstsp(ispact), ispact) = nbe
      ELSE
        STOP "LOEWDIN: states not sorted with spin"
      END IF
    END IF
  END DO
  IF (tprint) WRITE (*, *) 'Loewdin: nstsp,npoi=', nstsp, npoi

! ALLOCATE work spaces
  ktridig = (nstate + nstate*nstate)/2
  ALLOCATE (omatr(ktridig), oeigen(nstate))
  ALLOCATE (vect(nstate, nstate), vect2(nstate, nstate))
  ALLOCATE (psistate(nstate))
! IF(tdiag) ALLOCATE(hmatr(nstate,nstate),hpsi(kdfull2))
  ALLOCATE (hmatr(nstate, nstate), hpsi(kdfull2))

! orthogonlize spin-wise
  DO ispact = 1, 2

! set limits
    IF (nstsp(ispact) == 0) CYCLE
    IF (ispact == 1) THEN
      nlower = 1
      nupper = nstsp(1)
    ELSE
      nlower = nstsp(1) + 1
      nupper = nstsp(1) + nstsp(2)
    END IF

! compute the overlap matrix
    ntridig = 0
    DO nbe = nlower, nupper
      DO nbc = 1, nbe
        ntridig = 1 + ntridig
        omatr(ntridig) = wfovlp(q0(:, nbc), q0(:, nbe))
      END DO
    END DO
    IF (tprint) WRITE (*, '(a,200(1pg13.5))') 'nr. states, ntridig=', &
      nstsp(ispact), ntridig
    IF (tprint) WRITE (*, '(a,200(1pg13.5))') 'overlap matrix:', omatr(1:ntridig)

! diagonalize the overlap matrix
    CALL givens(omatr, oeigen, vect, nstsp(ispact), nstsp(ispact), nstate)
    IF (tprint) WRITE (*, '(a,200(1pg13.5))') 'eigenvalues:', &
      oeigen(1:nstsp(ispact))
    IF (tprint) THEN
      DO nbes = 1, nstsp(ispact)
        WRITE (*, '(a,200(1pg13.5))') 'eigenvectors:', &
          vect(1:nstsp(ispact), nbes)
      END DO
    END IF

! build the Loewdin matrix
    DO nbes = 1, nstsp(ispact)
      DO nbcs = 1, nstsp(ispact)
        vect(nbcs, nbes) = vect(nbcs, nbes)/(oeigen(nbes)**0.25D0)
      END DO
    END DO
    DO nbes = 1, nstsp(ispact)
      DO nbcs = 1, nstsp(ispact)
        vect2(nbcs, nbes) = SUM(vect(nbcs, 1:nstsp(ispact))*vect(nbes, 1:nstsp(ispact)))
      END DO
    END DO

    ii = 1
    IF (tdiag .AND. ii == 0) THEN

! now the branch for simultaneous Hamiltonian diagonalization

! compute Hamiltonian matrix
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
        hpsi = q0(:, nbe)
        CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 0)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          hmatr(nbcs, nbes) = wfovlp(q0(:, nbc), hpsi)
          IF (nbcs == nbes) hmatr(nbcs, nbes) = hmatr(nbcs, nbes) + amoy(nbes)
        END DO
      END DO

! transform Hamiltonian matrix TO orthonormalized basis
      DO nbes = 1, nstsp(ispact)
        DO nbcs = 1, nstsp(ispact)
          vect(nbcs, nbes) = &
            SUM(hmatr(nbcs, 1:nstsp(ispact))*vect2(1:nstsp(ispact), nbes))
        END DO
      END DO
      DO nbes = 1, nstsp(ispact)
        DO nbcs = 1, nstsp(ispact)
          hmatr(nbes, nbcs) = &
            SUM(vect2(1:nstsp(ispact), nbes)*vect(1:nstsp(ispact), nbcs))
        END DO
      END DO
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') ' new H matrix:'
      ntridig = 0
      DO nbes = 1, nstsp(ispact)
        IF (tprint) WRITE (*, '(200(1pg13.5))') hmatr(1:nstsp(ispact), nbes)
        DO nbcs = 1, nbes
          ntridig = 1 + ntridig
          omatr(ntridig) = hmatr(nbcs, nbes)
        END DO
      END DO

! diagonalize
      CALL givens(omatr, oeigen, vect, nstsp(ispact), nstsp(ispact), nstate)
      WRITE (*, '(a,200(1pg13.5))') ' new s.p.energies:', &
        oeigen(1:nstsp(ispact))
      amoy(1:nstsp(ispact)) = oeigen(1:nstsp(ispact))

! combine transformations, USE 'hmatr' as workspace
      DO nbes = 1, nstsp(ispact)
        DO nbcs = 1, nstsp(ispact)
          hmatr(nbcs, nbes) = &
            SUM(vect2(nbcs, 1:nstsp(ispact))*vect(1:nstsp(ispact), nbes))
        END DO
      END DO
      vect2 = hmatr

    END IF

! compose new wavefunctions
    DO ii = 1, nxyz
      psistate = 0D0
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          psistate(nbe) = psistate(nbe) + q0(ii, nbc)*vect2(nbcs, nbes)
        END DO
      END DO
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        q0(ii, nbe) = psistate(nbe)
      END DO
    END DO
    IF (tprint) THEN
      WRITE (*, '(a)') 'overl,after:'
      DO nbe = 1, nstate
        WRITE (*, '(100(1pg13.5))') (wfovlp(q0(:, nbc), q0(:, nbe)), nbc=1, nbe)
      END DO
    END IF

! test BLOCK
    IF (tdiag) THEN
      ntridig = 0
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
        hpsi = q0(:, nbe)
        CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 1)
        DO nbcs = 1, nbes
          nbc = npoi(nbcs, ispact)
          ntridig = 1 + ntridig
          omatr(ntridig) = wfovlp(q0(:, nbc), hpsi)
        END DO
      END DO
      CALL givens(omatr, oeigen, vect, nstsp(ispact), nstsp(ispact), nstate)
      WRITE (*, '(a,200(1pg13.5))') ' new s.p.energies:', &
        oeigen(1:nstsp(ispact))
      amoy(1:nstsp(ispact)) = oeigen(1:nstsp(ispact))
      DO ii = 1, nxyz
        psistate = 0D0
        DO nbes = 1, nstsp(ispact)
          nbe = npoi(nbes, ispact)
          DO nbcs = 1, nstsp(ispact)
            nbc = npoi(nbcs, ispact)
            psistate(nbe) = psistate(nbe) + q0(ii, nbc)*vect(nbcs, nbes)
          END DO
        END DO
        DO nbes = 1, nstsp(ispact)
          nbe = npoi(nbes, ispact)
          q0(ii, nbe) = psistate(nbe)
        END DO
      END DO
    END IF

    ! compute Hamiltonian matrix
! IF(tprint .AND. tdiag) THEN
    IF (.TRUE.) THEN
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
        hpsi = q0(:, nbe)
        CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 0)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          hmatr(nbcs, nbes) = wfovlp(q0(:, nbc), hpsi)
! IF(nbcs==nbes) hmatr(nbcs,nbes) = hmatr(nbcs,nbes) + amoy(nbes)
        END DO
        WRITE (*, *) ' variance:', nbe, SQRT(wfovlp(hpsi, hpsi) - hmatr(nbes, nbes)**2)
      END DO
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') ' H matrix after:'
      DO nbes = 1, nstsp(ispact)
        IF (tprint) WRITE (*, '(200(1pg13.5))') hmatr(1:nstsp(ispact), nbes)
      END DO
    END IF

  END DO

  DEALLOCATE (omatr, oeigen, vect)
  DEALLOCATE (npoi)
  DEALLOCATE (psistate)
! IF(tdiag) DEALLOCATE(hmatr,hpsi)
  DEALLOCATE (hmatr, hpsi)

END SUBROUTINE Loewdin

! ******************************

SUBROUTINE sortwf_energ(q0)

! ******************************

! Sorts wavefunctions accorting TO energies
!
! Input/Output:
! q0 = set of REAL s.p. wavefunctions TO be ortho-normalized

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)

  INTEGER :: i, imin, isav, n, nbe, nbes, ncc, ncs, ispact, ii
  REAL(DP) :: cs, emin
  REAL(DP) :: eord(kstate)
  REAL(DP) :: psistate(kstate)
  INTEGER :: isort(kstate)

!*********************************************************

! sort the s.p. energies

  DO n = 1, nstate
    eord(n) = amoy(n)
    isort(n) = n
  END DO

  DO ispact = 1, numspin
    DO n = 1, nstate
      IF (ispin(n) == ispact) THEN
        emin = 1D32
        imin = n
        DO i = n, nstate
          IF (ispin(i) == ispact .AND. eord(i) < emin) THEN
            emin = eord(i)
            imin = i
          END IF
        END DO
        isav = isort(imin)
        eord(imin) = eord(n)
        isort(imin) = isort(n)
        eord(n) = emin
        isort(n) = isav
      END IF
    END DO
  END DO

  WRITE (*, *) 'isort:', isort(1:nstate)
  WRITE (*, *) 'amoy:', amoy(1:nstate)
  WRITE (*, *) 'eord:', eord(1:nstate)

  amoy(1:nstate) = eord(1:nstate)
  DO ii = 1, nxyz
    psistate = 0D0
    DO nbe = 1, nstate
      psistate(nbe) = q0(ii, isort(nbe))
    END DO
    DO nbe = 1, nstate
      q0(ii, nbe) = psistate(nbe)
    END DO
  END DO

END SUBROUTINE sortwf_energ

SUBROUTINE schmid2(q0, aloc, tdiag)
! ******************************

! Serial version of orthogonalization with Loewdin scheme
!
! Input/Output:
! q0 = set of REAL s.p. wavefunctions TO be ortho-normalized
! tdiag = switch TO Hamiltonian diagonalization

  USE params
  USE util, ONLY: wfovlp
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  LOGICAL, INTENT(IN) :: tdiag

  INTEGER :: nbe, nbes
  REAL(DP) :: cs, emin
  REAL(DP) :: eord(kstate)
  INTEGER :: isort(kstate)

  REAL(DP), ALLOCATABLE :: omatr(:) ! tridiagonal storage
  REAL(DP), ALLOCATABLE :: oeigen(:)
  REAL(DP), ALLOCATABLE :: vect(:, :), vect2(:, :), hmatr(:, :)
  REAL(DP), ALLOCATABLE :: psistate(:), hpsi(:)
  INTEGER, ALLOCATABLE :: npoi(:, :)
  INTEGER :: ishift, ispact, ii, nbc, nbcs, ntridig, nlower, nupper, nstsp(2), ktridig
  INTEGER :: ncc, ncs, n

  LOGICAL, PARAMETER :: tprint = .TRUE.

!*********************************************************

  IF (tprint) WRITE (*, *) 'Loewdin entered: nstate=', nstate

  CALL sortwf_energ(q0)
  DO n = 1, nstate
    isort(n) = n
  END DO

! check that states comes IN contiguous blocks of spin
  ALLOCATE (npoi(nstate, 2))
  ispact = 1
  nstsp = 0
  npoi = 0
  DO nbe = 1, nstate
    IF (ispin(nbe) == ispact) THEN
      nstsp(ispact) = 1 + nstsp(ispact)
      npoi(nstsp(ispact), ispact) = nbe
    ELSE
      IF (ispact == 1) THEN
        ispact = 2
        nstsp(ispact) = 1 + nstsp(ispact)
        npoi(nstsp(ispact), ispact) = nbe
      ELSE
        STOP "LOEWDIN: states not sorted with spin"
      END IF
    END IF
  END DO
  IF (tprint) WRITE (*, *) 'Loewdin: nstsp,npoi=', nstsp, npoi

  DO nbes = 1, nstate
    nbe = isort(nbes)

    DO ncs = 1, nbes - 1
      ncc = isort(ncs)
      IF ((ispin(nbe) == ispin(ncc))) THEN
        cs = SUM(q0(1:nxyz, nbe)*q0(1:nxyz, ncc))*dvol
        q0(1:nxyz, nbe) = q0(1:nxyz, nbe) - cs*q0(1:nxyz, ncc)
      END IF
    END DO

    q0(1:nxyz, nbe) = q0(1:nxyz, nbe)/SQRT(SUM(q0(1:nxyz, nbe)**2)*dvol)

  END DO

! ALLOCATE work spaces
  ktridig = (nstate + nstate*nstate)/2
  ALLOCATE (omatr(ktridig), oeigen(nstate))
  ALLOCATE (vect(nstate, nstate), vect2(nstate, nstate))
  ALLOCATE (psistate(nstate))
! IF(tdiag) ALLOCATE(hmatr(nstate,nstate),hpsi(kdfull2))
  ALLOCATE (hmatr(nstate, nstate), hpsi(kdfull2))

! orthogonlize spin-wise
  DO ispact = 1, 2

! set limits
    IF (nstsp(ispact) == 0) CYCLE
    IF (ispact == 1) THEN
      nlower = 1
      nupper = nstsp(1)
    ELSE
      nlower = nstsp(1) + 1
      nupper = nstsp(1) + nstsp(2)
    END IF

    IF (tdiag) THEN

! now the branch for simultaneous Hamiltonian diagonalization

! compute Hamiltonian matrix
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
        hpsi = q0(:, nbe)
        CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 0)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          hmatr(nbcs, nbes) = wfovlp(q0(:, nbc), hpsi)
          IF (nbcs == nbes) hmatr(nbcs, nbes) = hmatr(nbcs, nbes) !+ amoy(nbes)
        END DO
        WRITE (*, *) ' variance b:', nbe, SQRT(wfovlp(hpsi, hpsi) - hmatr(nbes, nbes)**2)
      END DO

! IF(tprint) WRITE(*,*) 'H matrix before:'
      ntridig = 0
      cs = 0D0
      DO nbes = 1, nstsp(ispact)
! IF(tprint) WRITE(*,'(200(1pg13.5))') hmatr(1:nstsp(ispact),nbes)
        DO nbcs = 1, nbes
          IF (nbcs .NE. nbes) cs = hmatr(nbcs, nbes)**2 + cs
          ntridig = 1 + ntridig
          omatr(ntridig) = hmatr(nbcs, nbes)
        END DO
      END DO
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') '# H before:', &
        cs/(nstsp(ispact)*nstsp(ispact) - nstsp(ispact)), &
        (hmatr(nbes, nbes), nbes=1, nstsp(ispact))

! diagonalize
      CALL givens(omatr, oeigen, vect2, nstsp(ispact), nstsp(ispact), nstate)
      WRITE (*, '(a,200(1pg13.5))') ' new s.p.energies:', &
        oeigen(1:nstsp(ispact))
      amoy(1:nstsp(ispact)) = oeigen(1:nstsp(ispact))

! compose new wavefunctions
      DO ii = 1, nxyz
        psistate = 0D0
        DO nbes = 1, nstsp(ispact)
          nbe = npoi(nbes, ispact)
          DO nbcs = 1, nstsp(ispact)
            nbc = npoi(nbcs, ispact)
            psistate(nbe) = psistate(nbe) + q0(ii, nbc)*vect2(nbcs, nbes)
          END DO
        END DO
        DO nbes = 1, nstsp(ispact)
          nbe = npoi(nbes, ispact)
          q0(ii, nbe) = psistate(nbe)
        END DO
      END DO

    END IF

    ! compute Hamiltonian matrix
    IF (tprint .AND. tdiag) THEN
      cs = 0D0
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
        hpsi = q0(:, nbe)
        CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 0)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          hmatr(nbcs, nbes) = wfovlp(q0(:, nbc), hpsi)
          IF (nbcs .NE. nbes) cs = hmatr(nbcs, nbes)**2 + cs
! IF(nbcs==nbes) hmatr(nbcs,nbes) = hmatr(nbcs,nbes) + amoy(nbes)
        END DO
        WRITE (*, *) ' variance a:', nbe, SQRT(wfovlp(hpsi, hpsi) - hmatr(nbes, nbes)**2)
      END DO
! IF(tprint) WRITE(*,'(a,200(1pg13.5))') ' H matrix after:'
! DO nbes=1,nstsp(ispact)
! IF(tprint) WRITE(*,'(200(1pg13.5))') hmatr(1:nstsp(ispact),nbes)
! END DO
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') '# H after: ', &
        cs/(nstsp(ispact)*nstsp(ispact) - nstsp(ispact)), &
        (hmatr(nbes, nbes), nbes=1, nstsp(ispact))
    END IF

  END DO

  DEALLOCATE (omatr, oeigen, vect)
  DEALLOCATE (npoi)
  DEALLOCATE (psistate)
! IF(tdiag) DEALLOCATE(hmatr,hpsi)
  DEALLOCATE (hmatr, hpsi)

END SUBROUTINE schmid2

SUBROUTINE schmid3(q0, aloc, tdiag)
! ******************************

! Serial version of orthogonalization with Loewdin scheme
!
! Input/Output:
! q0 = set of REAL s.p. wavefunctions TO be ortho-normalized
! tdiag = switch TO Hamiltonian diagonalization

  USE params
  USE util, ONLY: wfovlp
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  LOGICAL, INTENT(IN) :: tdiag

  INTEGER :: nbe, nbes
  REAL(DP) :: cs, emin
  REAL(16) :: qs
  REAL(DP) :: eord(kstate)
  INTEGER :: isort(kstate)

  REAL(DP), ALLOCATABLE :: omatr(:) ! tridiagonal storage
  REAL(DP), ALLOCATABLE :: oeigen(:)
  REAL(DP), ALLOCATABLE :: vect(:, :), vect2(:, :), hmatr(:, :)
  REAL(16), ALLOCATABLE :: ovect(:, :, :), ovmatr(:, :, :)
  REAL(DP), ALLOCATABLE :: psistate(:), hpsi(:)
  INTEGER, ALLOCATABLE :: npoi(:, :)
  INTEGER :: ishift, ispact, ii, nbc, nbcs, ntridig, nlower, nupper, nstsp(2), ktridig
  INTEGER :: ncc, ncs, n

  LOGICAL, PARAMETER :: tprint = .TRUE.

!*********************************************************

  IF (tprint) WRITE (*, *) 'Loewdin entered: nstate=', nstate

  CALL sortwf_energ(q0)
  DO n = 1, nstate
    isort(n) = n
  END DO

! check that states comes IN contiguous blocks of spin
  ALLOCATE (npoi(nstate, 2))
  ispact = 1
  nstsp = 0
  npoi = 0
  DO nbe = 1, nstate
    IF (ispin(nbe) == ispact) THEN
      nstsp(ispact) = 1 + nstsp(ispact)
      npoi(nstsp(ispact), ispact) = nbe
    ELSE
      IF (ispact == 1) THEN
        ispact = 2
        nstsp(ispact) = 1 + nstsp(ispact)
        npoi(nstsp(ispact), ispact) = nbe
      ELSE
        STOP "LOEWDIN: states not sorted with spin"
      END IF
    END IF
  END DO
  IF (tprint) WRITE (*, *) 'Loewdin: nstsp,npoi=', nstsp, npoi

  ALLOCATE (ovmatr(nstate, nstate, 2), ovect(nstate, nstate, 2))
  ALLOCATE (hmatr(nstate, nstate), hpsi(kdfull2))
  ALLOCATE (vect(nstate, nstate), vect2(nstate, nstate))
  ALLOCATE (psistate(nstate))

  hmatr = 0D0
  DO nbes = 1, nstate
    nbe = isort(nbes)
    DO ncs = 1, nbes
      ncc = isort(ncs)
      IF ((ispin(nbe) == ispin(ncc))) THEN
        ovmatr(ncs, nbes, 1) = SUM(q0(1:nxyz, nbe)*q0(1:nxyz, ncc))*dvol
      END IF
      ovmatr(nbes, ncs, 1) = ovmatr(ncs, nbes, 1)
    END DO
  END DO

  DO ii = 1, 1
    IF (tprint) THEN
      WRITE (*, *) 'Norm overlap before: iteration=', ii
      DO nbes = 1, nstate
        WRITE (*, '(200(1pg13.5))') (ovmatr(ncs, nbes, ii), ncs=1, nstate)
      END DO
    END IF

    ovect(:, :, ii) = 0D0
    DO nbes = 1, nstate
      DO ncs = nbes + 1, nstate
        ovect(ncs, nbes, ii) = -ovmatr(ncs, nbes, ii)/ovmatr(ncs, ncs, ii)
      END DO
      ovect(nbes, nbes, ii) = 1D0
    END DO
    IF (tprint) THEN
      WRITE (*, *) 'Gram-Schmid matrix before:'
      DO nbes = 1, nstate
        WRITE (*, '(200(1pg13.5))') (ovect(ncs, nbes, ii), ncs=1, nstate)
      END DO
    END IF
    DO nbes = 1, nstate
      qs = 0D0
      DO ncs = nbes, nstate
        qs = qs + ovect(ncs, nbes, ii)* &
             SUM(ovmatr(ncs, nbes:nstate, ii)*ovect(nbes:nstate, nbes, ii))
      END DO
      psistate(nbes) = 1D0/SQRT(qs)
    END DO
    WRITE (*, *) 'qs=', qs
    DO nbes = 1, nstate
      DO ncs = nbes, nstate
        ovect(ncs, nbes, ii) = ovect(ncs, nbes, ii)*psistate(nbes)
      END DO
    END DO
    IF (tprint) THEN
      WRITE (*, *) 'Gram-Schmid matrix after:'
      DO nbes = 1, nstate
        WRITE (*, '(200(1pg13.5))') (ovect(ncs, nbes, ii), ncs=1, nstate)
      END DO
      WRITE (*, *) 'Norm overlap intermediate:'
      DO nbes = 1, nstate
        DO ncs = 1, nstate
          qs = 0D0
          DO n = 1, nstate
            qs = qs + ovect(n, ncs, ii)* &
                 SUM(ovmatr(n, 1:nstate, ii)*ovect(1:nstate, nbes, ii))
          END DO
          vect2(ncs, nbes) = qs
        END DO
        WRITE (*, '(200(1pg13.5))') vect2(1:nstate, nbes)
      END DO
      IF (ii == 1) ovmatr(:, :, 2) = vect2
    END IF
    IF (ii == 2) THEN
      DO nbes = 1, nstate
        DO ncs = 1, nstate
          vect2(ncs, nbes) = SUM(ovect(ncs, 1:nstate, 1)*ovect(1:nstate, nbes, 2))
        END DO
      END DO
    ELSE
      vect2 = ovect(:, :, 1)
    END IF
  END DO
! ALLOCATE work spaces
  ktridig = (nstate + nstate*nstate)/2
  ALLOCATE (omatr(ktridig), oeigen(nstate))
! IF(tdiag) ALLOCATE(hmatr(nstate,nstate),hpsi(kdfull2))

! orthogonlize spin-wise
  DO ispact = 1, 2

! set limits
    IF (nstsp(ispact) == 0) CYCLE
    IF (ispact == 1) THEN
      nlower = 1
      nupper = nstsp(1)
    ELSE
      nlower = nstsp(1) + 1
      nupper = nstsp(1) + nstsp(2)
    END IF

    IF (tdiag) THEN

! now the branch for simultaneous Hamiltonian diagonalization

! compute Hamiltonian matrix
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
        hpsi = q0(:, nbe)
        CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 0)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          hmatr(nbcs, nbes) = wfovlp(q0(:, nbc), hpsi)
          IF (nbcs == nbes) hmatr(nbcs, nbes) = hmatr(nbcs, nbes) !+ amoy(nbes)
        END DO
        WRITE (*, *) ' variance b:', nbe, SQRT(wfovlp(hpsi, hpsi) - hmatr(nbes, nbes)**2)
      END DO

! IF(tprint) WRITE(*,*) 'H matrix before:'
      ntridig = 0
      cs = 0D0
      DO nbes = 1, nstsp(ispact)
! IF(tprint) WRITE(*,'(200(1pg13.5))') hmatr(1:nstsp(ispact),nbes)
        DO nbcs = 1, nbes
          IF (nbcs .NE. nbes) cs = hmatr(nbcs, nbes)**2 + cs
          ntridig = 1 + ntridig
          omatr(ntridig) = hmatr(nbcs, nbes)
        END DO
      END DO
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') '# H before:', &
        cs/(nstsp(ispact)*nstsp(ispact) - nstsp(ispact)), &
        (hmatr(nbes, nbes), nbes=1, nstsp(ispact))

! diagonalize
      CALL givens(omatr, oeigen, vect2, nstsp(ispact), nstsp(ispact), nstate)
      WRITE (*, '(a,200(1pg13.5))') ' new s.p.energies:', &
        oeigen(1:nstsp(ispact))
      amoy(1:nstsp(ispact)) = oeigen(1:nstsp(ispact))

    END IF

! compose new wavefunctions
    DO ii = 1, nxyz
      psistate = 0D0
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          psistate(nbe) = psistate(nbe) + q0(ii, nbc)*vect2(nbcs, nbes)
        END DO
      END DO
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        q0(ii, nbe) = psistate(nbe)
      END DO
    END DO

    ! test overlap
    IF (tprint) THEN
      hmatr = 0D0
      DO nbes = 1, nstate
        nbe = isort(nbes)
        DO ncs = 1, nbes
          ncc = isort(ncs)
          IF ((ispin(nbe) == ispin(ncc))) THEN
            hmatr(ncs, nbes) = SUM(q0(1:nxyz, nbe)*q0(1:nxyz, ncc))*dvol
            hmatr(nbes, ncs) = hmatr(ncs, nbes)
          END IF
        END DO
      END DO
      WRITE (*, *) 'Norm overlap after:'
      DO nbes = 1, nstate
        WRITE (*, '(200(1pg13.5))') (hmatr(ncs, nbes), ncs=1, nstate)
      END DO
    END IF

    ! compute Hamiltonian matrix
! IF(tprint .AND. tdiag) THEN
    IF (tprint) THEN
      cs = 0D0
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
        hpsi = q0(:, nbe)
        CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 0)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          hmatr(nbcs, nbes) = wfovlp(q0(:, nbc), hpsi)
          IF (nbcs .NE. nbes) cs = hmatr(nbcs, nbes)**2 + cs
! IF(nbcs==nbes) hmatr(nbcs,nbes) = hmatr(nbcs,nbes) + amoy(nbes)
        END DO
! cs = 00D0
! DO nbcs=1,nstsp(ispact)
! cs=cs+hmatr(nbes,nbcs)*hmatr(nbcs,nbes)
! END DO
        WRITE (*, *) ' variance a:', nbe, &
          SQRT(wfovlp(hpsi, hpsi) - hmatr(nbes, nbes)**2)
! SQRT(wfovlp(hpsi,hpsi)-cs)
      END DO
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') ' H matrix after:'
      DO nbes = 1, nstsp(ispact)
        IF (tprint) WRITE (*, '(a,200(1pg13.5))') '# H after: ', &
          cs/(nstsp(ispact)*nstsp(ispact) - nstsp(ispact)) !,&
! (hmatr(nbes,nbes),nbes=1,nstsp(ispact))
        IF (tprint) WRITE (*, '(200(1pg13.5))') hmatr(1:nstsp(ispact), nbes)
      END DO
    END IF

  END DO

  DEALLOCATE (ovmatr, ovect)
  DEALLOCATE (omatr, oeigen, vect)
  DEALLOCATE (npoi)
  DEALLOCATE (psistate)
! IF(tdiag) DEALLOCATE(hmatr,hpsi)
  DEALLOCATE (hmatr, hpsi)

END SUBROUTINE schmid3

