
!-----tstep_exp-------------------------------------------------------------

SUBROUTINE tstep_exp(q0, aloc, rho, it, qwork, timagtime)

! One electronic time step by exponential evolution, consisting
! of a half step TO produce an intermediate mean field followed
! by a full with using this mean field.
!
! q0 = s.p. wavefunctions TO be propagated
! aloc = array for local mean field
! rho = array for local density
! it = nr. of time step
! qwork = work space for s.p. wavefunctions
! timagtime = LOGICAL variable, switch TO imaginary time step

  USE params
  USE twost, ONLY: tnearest

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT) :: akv(kdfull2)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it
  COMPLEX(DP), INTENT(OUT) :: qwork(kdfull2, kstate)
  LOGICAL, INTENT(IN) :: timagtime

  INTEGER :: nb, nterms
  COMPLEX(DP) :: q1(kdfull2)
  COMPLEX(DP) :: cdtact

! The PARAMETER 'tnorotate' de-activates the subtraction of the
! Lagrangian matrix IN the SIC step. The version of exponential
! evolution with subtraction of the Lagrangian matrix is found
! IN 'exp_evolp'. The strategy needs yet testing. It was not
! really beneficial so far. (01/2013)
!
! what is the current STATUS of this issue ? F.L. (03/2017)
  LOGICAL, PARAMETER :: tnorotate = .true.

!----------------------------------------------------------------------

#if(mpi)
  STOP 'exponential evolution not yet MPI parallelized'
#else
  myn = 0
#endif

!STOP

#if(raregas)
  IF (nc + NE + nk > 0) STOP 'TSTEP_EXP not appropriate for rare gas'
#endif

! one half time step TO define new mean field
! USE exponential evolution TO second order

  IF (ifsicp == 5) psisavex = q0

  IF (.NOT. timagtime) THEN
    cdtact = CMPLX(dt1/2D0, 0D0, DP)
    IF (tnorotate .OR. ifsicp < 8) THEN
! WRITE(*,*) 'EXPEVOL: half step, dt=',cdtact
      DO nb = 1, nstate
        qwork(:, nb) = q0(:, nb)
        CALL exp_evol(qwork(:, nb), aloc, nb, 4, cdtact, q1)
      END DO
    ELSE
      qwork = q0
      CALL exp_evolp(qwork, aloc, 4, cdtact, q1, q0)
    END IF

    IF (tnearest .AND. ifsicp .GE. 8) CALL eval_unitrot(qwork, q0)

! compute mean field at half time step
    CALL dyn_mfield(rho, aloc, qwork, dt1*0.5D0, it)

  END IF

! full time step TO next wavefunctions
! USE exponential evolution TO fourth order

  IF (timagtime) THEN
    cdtact = CMPLX(0D0, -dt1, DP)
  ELSE
    cdtact = CMPLX(dt1, 0D0, DP)
  END IF

  nterms = 4

  IF (tnorotate .OR. ifsicp < 8) THEN
    DO nb = 1, nstate
      CALL exp_evol(q0(:, nb), aloc, nb, nterms, cdtact, q1)
    END DO
  ELSE
    qwork = q0
    CALL exp_evolp(q0, aloc, nterms, cdtact, q1, qwork)
  END IF

  IF (tnearest .AND. ifsicp .GE. 8 .AND. .NOT. timagtime) CALL eval_unitrot(q0, qwork)

! switch of absorption after times step 'ntref'
  IF (ntref > 0 .AND. it > ntref) nabsorb = 0

! compute mean field at new time

  CALL dyn_mfield(rho, aloc, q0, dt1, it)

  RETURN
END SUBROUTINE tstep_exp

!-----exp_evol-------------------------------------------------------------

SUBROUTINE exp_evol(qact, aloc, nbe, norder, dtact, qwork)

! Propagation of the wavefunction of state 'nbe' by Taylor
! expanded evolution:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies IN momentum space
! nbe = NUMBER of state
! norder = order of expansion (4 recommended for full step))
! dtact = time step

! Note: The propagation uses the ACTION of the Hamiltonian
! WHERE the diagonal element (s.p.energy) is subtracted.
! That diagonal element is evaluated IN the first order
! CALL 'nterm=1'.

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
  INTEGER, INTENT(IN) :: nbe
  INTEGER, INTENT(IN) :: norder
  COMPLEX(DP), INTENT(IN) :: dtact
  COMPLEX(DP), INTENT(OUT) :: qwork(kdfull2)

  COMPLEX(DP) :: dti, cfac
  INTEGER :: i, ilocbas, isig, nterm

!----------------------------------------------------------------------

  IF (ispin(nrel2abs(nbe)) == 1) THEN
    ilocbas = 1
  ELSE IF (ispin(nrel2abs(nbe)) == 2) THEN
    ilocbas = nxyz + 1
  ELSE
    STOP " EXPEVOL: spin index must be 1 or 2"
  END IF
  IF (ABS(IMAG(dtact)) > 1D-10) THEN
    isig = -1
  ELSE
    isig = 1
  END IF
!isig=-1
  dti = dtact*CMPLX(0D0, 1D0, DP)
  cfac = CMPLX(1D0, 0D0, DP)
  DO i = 1, nxyz
    qwork(i) = qact(i)
  END DO
  DO nterm = 1, norder
    CALL hpsi(qwork, aloc(ilocbas), nbe, isig*nterm)
    cfac = -dti/nterm*cfac
    DO i = 1, nxyz
      qact(i) = qact(i) + cfac*qwork(i)
    END DO
  END DO
  RETURN
END SUBROUTINE exp_evol

!-----exp_evolp-------------------------------------------------------------

SUBROUTINE exp_evolp(qact, aloc, norder, dtact, qwork, psi)

! Propagation of all s.p. wavefuntions by Taylor expanded
! exponential evolution (version needed for SIC):
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies IN momentum space
! norder = order of epxansion (4 recommended for full step))
! dtact = time step
! qwork = work space for COMPLEX wavefunction
! psi = set of wavefunctions before the step

! Note: The propagation uses the ACTION of the Hamiltonian
! WHERE the diagonal element (s.p.energy) is subtracted.
! That diagonal element is evalutaed IN the first order
! CALL 'nterm=1'.

  USE params
  USE util, ONLY: wfovlp
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
  INTEGER, INTENT(IN) :: norder
  COMPLEX(DP), INTENT(IN) :: dtact
  COMPLEX(DP), INTENT(OUT) :: qwork(kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

  COMPLEX(DP), ALLOCATABLE :: chmatrix(:, :)

  COMPLEX(DP) :: dti, cfac, cacc(kstate)
  INTEGER :: i, ilocbas
  INTEGER :: na, nbe, ncs, nterm

!----------------------------------------------------------------------

!WRITE(*,*) 'entering EXP_EVOLP'

  ALLOCATE (chmatrix(kstate, kstate))

  dti = dtact*CMPLX(0D0, 1D0, DP)

! compute H-matrix, store h*psi wavefunctions

  DO nbe = 1, nstate
    ilocbas = 1 + (ispin(nrel2abs(nbe)) - 1)*nxyz
    CALL hpsi(qact(1, nbe), aloc(ilocbas), nbe, 1)
    DO ncs = 1, nstate
      IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(ncs))) THEN
        chmatrix(ncs, nbe) = wfovlp(psi(:, ncs), qact(:, nbe))
      ELSE
        chmatrix(ncs, nbe) = CMPLX(0D0, 0D0, DP)
      END IF
    END DO
  END DO
! symmetrize H-matrix
  DO nbe = 1, nstate; DO ncs = 1, nbe - 1
      chmatrix(ncs, nbe) = (chmatrix(ncs, nbe) + CONJG(chmatrix(nbe, ncs)))/2D0
      chmatrix(nbe, ncs) = CONJG(chmatrix(ncs, nbe))
    END DO; END DO

! now the Taylor expansion (recycle stored h*psi IN first step)

  DO nbe = 1, nstate
    ilocbas = 1 + (ispin(nrel2abs(nbe)) - 1)*nxyz
    cfac = CMPLX(1D0, 0D0, DP)
    DO nterm = 1, norder
      ! H*psi and prepare subtraction factors
      IF (nterm == 1) THEN
        qwork(:) = qact(:, nbe)
        qact(:, nbe) = psi(:, nbe) ! restore original wavefunctions
        cacc(:) = chmatrix(:, nbe)
      ELSE
        DO ncs = 1, nstate
          IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(ncs))) THEN
            cacc(ncs) = CMPLX(0D0, 0D0, DP)
            DO na = 1, nstate
              IF (ispin(nrel2abs(na)) == ispin(nrel2abs(ncs))) THEN
                cacc(ncs) = cacc(ncs) + chmatrix(ncs, na)*wfovlp(psi(:, na), qwork)
              END IF
            END DO
          END IF
        END DO
        CALL hpsi(qwork, aloc(ilocbas), nbe, 2)
      END IF
      ! project H-matrix
      DO ncs = 1, nstate
        IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(ncs))) THEN
          qwork(:) = qwork(:) - psi(:, ncs)*cacc(ncs)
        END IF
      END DO
      ! accumulate TO Taylor series
      cfac = -dti/nterm*cfac
      qact(:, nbe) = qact(:, nbe) + cfac*qwork(:)
    END DO
  END DO

  DEALLOCATE (chmatrix)

  RETURN
END SUBROUTINE exp_evolp

!-----hpsi -------------------------------------------------------------

SUBROUTINE hpsi(qact, aloc, nbe, itpri)

! ACTION of mean-field Hamiltonian on one s.p. wavefunction:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies IN momentum space
! nbe = NUMBER of state
! itpri = switch for computing s.p. energies (for ABVS(itpri)=1)
! <0 switches TO subtract mean-value of s.p. energy

  USE params
  USE util, ONLY: wfovlp
  USE kinetic
  USE twost

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN) :: akv(kdfull2)
  INTEGER, INTENT(IN) :: nbe
  INTEGER, INTENT(IN) :: itpri
  COMPLEX(DP), ALLOCATABLE :: qex(:)

! workspaces
  COMPLEX(DP), ALLOCATABLE :: q1(:), q2(:), q1fine(:), q2fine(:), qactfine(:)
  COMPLEX(DP), ALLOCATABLE :: qarray(:, :, :), qarrayfine(:, :, :)
  REAL(DP) :: wnorm
  LOGICAL :: tpri
  LOGICAL, PARAMETER :: tsubmean = .TRUE.
  LOGICAL, PARAMETER :: ttest = .FALSE.
  INTEGER :: i, is, na
  COMPLEX(DP) :: cf

!----------------------------------------------------------------------

  tpri = ABS(itpri) == 1

  ALLOCATE (q1(kdfull2), q2(kdfull2))
  q1 = (0D0, 0D0)
  q2 = (0D0, 0D0)
! ACTION of kinetic energy

#if(netlib_fft|fftw_cpu)
  CALL fftf(qact, q1)
  DO i = 1, nxyz
    q1(i) = akv(i)*q1(i)
  END DO
  CALL fftback(q1, q2)
#else
  CALL ckin3d(qact, q1)
#endif

! ACTION of potential and non-local PsP (optionally)

  IF (ipsptyp == 1) THEN
    CALL nonlocalc(qact, q1, 0)

    IF (tpri) enonlo(nbe) = wfovlp(qact, q1)
    DO i = 1, nxyz
      q1(i) = q1(i) + qact(i)*aloc(i)
    END DO
  ELSE
    DO i = 1, nxyz
      q1(i) = qact(i)*aloc(i)
    END DO
  END IF

  IF (ifsicp == 5) THEN
    ALLOCATE (qex(kdfull2))
    CALL exchg(qact, qex, nbe)
    q1 = q1 + qex
    DEALLOCATE (qex)
  END IF

! subtract SIC potential for state NBE
  IF (ifsicp .GE. 8) THEN
    is = ispin(nrel2abs(nbe))
    DO na = 1, nstate
      IF (ispin(nrel2abs(na)) == is) THEN
        cf = wfovlp(psiut(:, na), qact)
        DO i = 1, nxyz
          q1(i) = q1(i) - qnewut(i, na)*cf
        END DO
      END IF
    END DO

  END IF

  IF (tpri) THEN
    epotsp(nbe) = wfovlp(qact, q1)
    amoy(nbe) = ekinsp(nbe) + epotsp(nbe)
    q2 = q1 + q2
    spvariance(nbe) = SQRT(MAX(REAL(wfovlp(q2, q2), DP) - ABS(wfovlp(qact, q2))**2, 1D-99))
    is = ispin(nrel2abs(nbe))
    IF (ttest) WRITE (*, '(a,2i4,5(1pg13.5))') &
      ' HPSI: nbe,is,esp,var=', nbe, is, amoy(nbe), spvariance(nbe), &
      ekinsp(nbe), epotsp(nbe), REAL(wfovlp(q2, q2), DP)
    CALL flush (6)
  ELSE
    q2 = q1 + q2
  END IF

  IF (itpri < 0) THEN
    qact = q2 - amoy(nbe)*qact
  ELSE
    qact = q2
  END IF

  DEALLOCATE (q1, q2)

  RETURN
END SUBROUTINE hpsi

SUBROUTINE hpsi_boostinv(qact, aloc, current, rho, nbe)

! ACTION of boost-invariant Hamiltonian on one s.p. wavefunction:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! current = average local momentum IN x,y, and z-direction
! nbe = NUMBER of state
! The routine requires that 'current' has been accumulated before.

  USE params
  USE kinetic
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  REAL(DP), INTENT(IN) :: current(kdfull2, 3)
  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: nbe

! workspaces
  COMPLEX(DP), ALLOCATABLE :: q1(:), q2(:)

  COMPLEX(DP) :: cf

!----------------------------------------------------------------------

! check availability of FFT
  IF (.NOT. ALLOCATED(akv)) STOP "HPSI_BOOSTINVARIANT requires FFT"
  ALLOCATE (q1(kdfull2), q2(kdfull2))

  q1 = qact
  CALL hpsi(qact, aloc, nbe, 0)

  CALL xgradient_rspace(q1, q2)
  qact = qact - h2m*current(:, 1)*q2
  q2 = current(:, 1)*q1
  CALL xgradient_rspace(q2, q2)
  qact = qact - h2m*q2

  CALL ygradient_rspace(q1, q2)
  qact = qact - h2m*current(:, 2)*q2
  q2 = current(:, 2)*q1
  CALL ygradient_rspace(q2, q2)
  qact = qact - h2m*q2

  CALL zgradient_rspace(q1, q2)
  qact = qact - eye*h2m*current(:, 3)*q2
  q2 = current(:, 3)*q1
  CALL zgradient_rspace(q2, q2)
  qact = qact - eye*h2m*q2

  qact = qact + h2m* &
         (current(:, 1)**2 + current(:, 2)**2 + current(:, 3)**2)

  WRITE (*, *) ' HPSI_BOOSTINV: E_coll=', dvol*h2m*SUM(rho(:)* &
                                                       (current(:, 1)**2 + current(:, 2)**2 + current(:, 3)**2))

  RETURN
END SUBROUTINE hpsi_boostinv
