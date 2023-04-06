
SUBROUTINE hamdiag(q0, aloc)
! ******************************

! Diagonalization of actual mean field
!
! Input/Output:
! q0 = set of REAL s.p. wavefunctions, already ortho-normalized
! aloc = given local potential
! (non-local part indirectly through 'rhpsi')
! tdiag = switch TO Hamiltonian diagonalization

  USE params
  USE util, ONLY: wfovlp
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)

  INTEGER :: nbe, nbes
  REAL(DP) :: cs, emin
  REAL(DP) :: eord(kstate)
  INTEGER :: isort(kstate)

  REAL(DP), ALLOCATABLE :: tmatr(:) ! tridiagonal storage
  REAL(DP), ALLOCATABLE :: eigen(:)
  REAL(DP), ALLOCATABLE :: vect(:, :), hmatr(:, :)
  REAL(DP), ALLOCATABLE :: vect0(:, :), hmatr0(:, :)
  REAL(DP), ALLOCATABLE :: psistate(:), hpsi(:)
  INTEGER, ALLOCATABLE :: npoi(:, :)
  INTEGER :: ishift, ispact, ii, nbc, nbcs, ntridig, nlower, nupper, nstsp(2), ktridig
  INTEGER :: ncc, ncs, n

  LOGICAL, PARAMETER :: tprint = .TRUE.
  LOGICAL, PARAMETER :: tsimp_diag = .FALSE. ! choice of diagonalization scheme

!*********************************************************

  IF (tprint) THEN
    WRITE (*, *) 'HAMDIAG entered: nstate=', nstate
    OPEN (199, POSITION='append', FILE='sspvariances.'//outnam)
    WRITE (199, *) '# HAMDIAG'
    CLOSE (199)
  END IF

  CALL sortwf_energ(q0)
  DO n = 1, nstate
    isort(n) = n
  END DO

! check that states come IN contiguous blocks of spin
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
        STOP "HAMDIAG: states not sorted with spin"
      END IF
    END IF
  END DO
  IF (tprint) WRITE (*, *) 'HAMDIAG: nstsp,npoi=', nstsp, npoi

! ALLOCATE work spaces
  ktridig = (nstate + nstate*nstate)/2
  ALLOCATE (tmatr(ktridig), eigen(nstate))
  ALLOCATE (vect(nstate, nstate), hmatr(nstate, nstate), psistate(nstate))
  ALLOCATE (hpsi(kdfull2))

! work up spin-wise
  DO ispact = 1, numspin

! set limits
    IF (nstsp(ispact) == 0) CYCLE
    IF (ispact == 1) THEN
      nlower = 1
      nupper = nstsp(1)
    ELSE
      nlower = nstsp(1) + 1
      nupper = nstsp(1) + nstsp(2)
    END IF

! compute Hamiltonian matrix
    DO nbes = 1, nstsp(ispact)
      nbe = npoi(nbes, ispact)
      ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK
      hpsi = q0(:, nbe)
      CALL rhpsi(hpsi, aloc(1 + ishift), nbe, 0)
      DO nbcs = 1, nstsp(ispact)
        nbc = npoi(nbcs, ispact)
        hmatr(nbcs, nbes) = wfovlp(q0(:, nbc), hpsi)
! IF(nbcs==nbes) hmatr(nbcs,nbes) = hmatr(nbcs,nbes) !+ amoy(nbes)
      END DO
      IF (tprint) WRITE (*, *) ' variance b:', nbe, &
        SQRT(ABS(wfovlp(hpsi, hpsi) - hmatr(nbes, nbes)**2))
    END DO

    IF (tsimp_diag) THEN
! diaginalization by simple gradient iteration IN mattrix space
      nbe = nstsp(ispact)
      ALLOCATE (vect0(nbe, nbe), hmatr0(nbe, nbe))
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          hmatr0(nbcs, nbes) = hmatr(nbcs, nbes)
        END DO
      END DO
      CALL diagsimp(hmatr0, vect0, eigen, nstsp(ispact))
      DO nbes = 1, nstsp(ispact)
        nbe = npoi(nbes, ispact)
        DO nbcs = 1, nstsp(ispact)
          nbc = npoi(nbcs, ispact)
          vect(nbcs, nbes) = vect0(nbcs, nbes)
        END DO
      END DO
      DEALLOCATE (vect0, hmatr0)
    ELSE
! diagonalization by Givens-Householder algorithm
! map TO upper triangular storage
      ntridig = 0
      cs = 0D0
      DO nbes = 1, nstsp(ispact)
        DO nbcs = 1, nbes
          IF (nbcs .NE. nbes) cs = hmatr(nbcs, nbes)**2 + cs
          ntridig = 1 + ntridig
          tmatr(ntridig) = hmatr(nbcs, nbes)
        END DO
      END DO
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') '# H variance, s.p.energies:', &
        cs/(nstsp(ispact)*nstsp(ispact) - nstsp(ispact)), &
        (hmatr(nbes, nbes), nbes=1, nstsp(ispact))
! diagonalize
      CALL givens(tmatr, eigen, vect, nstsp(ispact), nstsp(ispact), nstate)
      IF (tprint) WRITE (*, '(a,200(1pg13.5))') ' new s.p.energies:', &
        eigen(1:nstsp(ispact))
!
    END IF
    amoy(nlower:nupper) = eigen(1:nstsp(ispact))

! compose new wavefunctions
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

! test Hamiltonian matrix
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

  DEALLOCATE (tmatr, eigen, vect, hpsi)
  DEALLOCATE (npoi)
  DEALLOCATE (psistate)

END SUBROUTINE hamdiag

!---------------------------------------------------------------------------
! Employs a simple gradient step method for approximate diagonalization of
! the Hamiltonian matrix.
!
! List variables:
! hmatr = Hamiltonian matrix.
! eigvec = eigenvectors of the Hamiltonian Matrix.
! eigval = eigenvalues of the Hamiltonian Matrix.
! ndimact = DIMENSION of matries.
!
!---------------------------------------------------------------------------

SUBROUTINE diagsimp(hmatr, eigvec, eigval, ndimact)

  USE params, ONLY: small
  IMPLICIT NONE

  REAL(8), INTENT(IN OUT) :: hmatr(ndimact, ndimact)
  REAL(8), INTENT(OUT) :: eigvec(ndimact, ndimact)
  REAL(8), INTENT(OUT) :: eigval(ndimact)
  INTEGER, INTENT(IN) :: ndimact

  INTEGER, PARAMETER :: itdiagmax = 500
  REAL(8), PARAMETER :: diagstep = 0.05D0
  REAL(8), PARAMETER :: diagprecis = 1D-6
  LOGICAL, PARAMETER :: tprint = .FALSE.

  REAL(8), ALLOCATABLE :: workv1(:, :), workv2(:, :), propag(:, :)
  REAL(8) :: variance, sumeig, hmin, sumeig2
  INTEGER :: na, nb, itdiag

!----------------------------------------------------------------------

  IF (ndimact == 1) THEN
    eigvec(1, 1) = 1D0
    RETURN
  END IF

  ALLOCATE (workv1(ndimact, ndimact), workv2(ndimact, ndimact))
  ALLOCATE (propag(ndimact, ndimact))
  workv1 = 0D0
  workv2 = 0D0

  IF (tprint) THEN
    WRITE (6, '(a)') 'H matrix:'
    DO na = 1, ndimact
      WRITE (6, '(10(1pg13.5))') (hmatr(na, nb), nb=1, ndimact)
    END DO
  END IF
!***********************************************************
! prepare stepping matrix and eigenvectors
!***********************************************************
!search for minimum diagonal ENTRY
  hmin = 1D99
  DO na = 1, ndimact
    hmin = MIN(hmatr(na, na), hmin)
  END DO

  hmatr = -diagstep*hmatr
  eigvec = 0.0d0
  propag = hmatr
  DO na = 1, ndimact
    propag(na, na) = 1.0D0 + hmatr(na, na) + diagstep*hmin
    eigvec(na, na) = 1.0D0
  END DO
! iterate
  iterations: DO itdiag = 1, itdiagmax
! analysis of variance
    workv1 = MATMUL(hmatr(:, :), eigvec(:, :))
    workv2 = MATMUL(hmatr(:, :), workv1(:, :))
    variance = 0D0
    DO na = 1, ndimact
      sumeig = SUM(eigvec(:, na)*workv1(:, na))
      sumeig2 = SUM(eigvec(:, na)*workv2(:, na))
      variance = variance + (sumeig2 - sumeig**2)
      eigval(na) = -sumeig/diagstep
    END DO
    variance = SQRT(MAX(ABS(variance), small)/ndimact)
    IF (tprint) WRITE (6, '(a,2i5,200(1pg12.4))') &
      'diag: BLOCK,iter,variance:', ndimact, itdiag, variance, eigval
    IF (variance < diagprecis) EXIT iterations
! stepping
! eigvec = workv1
    eigvec = MATMUL(propag(:, :), eigvec(:, :))
    CALL orthmatr(eigvec, ndimact)
  END DO iterations

  DEALLOCATE (workv1, workv2)
  RETURN
END SUBROUTINE diagsimp

SUBROUTINE orthmatr(eigvec, ndimact)

!---------------------------------------------------------------------------
!
! Employs the Gram-Schmidt algorithm TO orthonormalize matrix
! of eigenvectors.
!
! List variables:
! eigvec = input non-orthonormal eigenvectors and returns the
! orthonormal ones.
! ndimact = DIMENSION of matrix.
!---------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(8), INTENT(IN OUT) :: eigvec(ndimact, ndimact)
  INTEGER, INTENT(IN) :: ndimact

  LOGICAL, PARAMETER :: ttest = .false.

  INTEGER :: na, nb
  REAL(8) :: ov
!----------------------------------------------------------------------

  DO na = 1, ndimact
    DO nb = 1, na - 1
      ov = sum(eigvec(:, nb)*eigvec(:, na))
      eigvec(:, na) = eigvec(:, na) - ov*eigvec(:, nb)
    END DO
    ov = 1.0D0/SQRT(sum(eigvec(:, na)**2))
    eigvec(:, na) = ov*eigvec(:, na)
  END DO

! test
  IF (.NOT. ttest) RETURN

  WRITE (6, '(a)') ' test orthonormalization:'
  DO na = 1, ndimact
    DO nb = 1, ndimact
      ov = sum(eigvec(:, nb)*eigvec(:, na))
      WRITE (6, '(2i5,1pg13.5)') na, nb, ov
    END DO
  END DO

  RETURN
END SUBROUTINE orthmatr
