
#ifdef REALSWITCH
MODULE twostr
  USE params
  USE kinetic
  USE twost_util
  USE orthmat
  IMPLICIT NONE

  INTERFACE spmomsmatrix
! automated choice between REAL and COMPLEX versions of spmomsmatrix,
! with spmomsmatrix_r, spmomsmatrix_c both defined at first compilation of twostr.
    MODULE PROCEDURE spmomsmatrix_r, spmomsmatrix_c
  END INTERFACE spmomsmatrix

  INTERFACE dalphabeta
! same as before, for the REAL, COMPLEX and REAL+COMPLEX versions of dalphabeta.
    MODULE PROCEDURE dalphabeta_r, dalphabeta_rc, dalphabeta_c
  END INTERFACE dalphabeta

  LOGICAL, PARAMETER :: tconv = .true. ! TO PRINT convergence

  COMPLEX(DP), PRIVATE, ALLOCATABLE :: qnewc(:, :)
  COMPLEX(DP), PRIVATE, ALLOCATABLE :: psicut(:, :)
  COMPLEX(DP), ALLOCATABLE :: rExpDABoldc(:, :, :)
  COMPLEX(DP), ALLOCATABLE, SAVE :: vecsc(:, :, :) ! searched eigenvectors

  REAL(DP), PRIVATE, ALLOCATABLE :: qnewr(:, :)
  REAL(DP), PRIVATE, ALLOCATABLE :: psirut(:, :)
  REAL(DP), ALLOCATABLE :: rExpDABold(:, :, :)!MV added
  REAL(DP), ALLOCATABLE, SAVE :: vecsr(:, :, :) ! searched eigenvectors

!
!REAL(DP),ALLOCATABLE :: usicall(:,:)
!COMMON /twost/ qnewr,qnew,psirut,psiut

  REAL(DP) :: step = 0.1D0, precis = 1D-6, precisfact = 1D-3, phiini = 0D0
  REAL(DP) :: ener_2st(2)
  INTEGER :: symutbegin = 100, itut
!COMMON /twostut/ step,precis,symutbegin,radmaxsym,itut
  REAL(DP) :: dampopt = 0.5D0, steplow = 0.01D0, steplim = 1.2D0 ! optimized stepsize
  LOGICAL :: toptsicstep = .true. ! switch TO optimized step

  INTEGER, PRIVATE :: kdim
! matrices of radial moments

#endif

#ifdef COMPLEXSWITCH
  MODULE twost
    USE params
    USE kinetic
    USE twostr, ONLY: step, precis, symutbegin, precisfact, itut, vecsr, tconv, ener_2st, toptsicstep, steplim, dampopt
    USE orthmat
    IMPLICIT NONE

    COMPLEX(DP), ALLOCATABLE :: qnewut(:, :)
    COMPLEX(DP), ALLOCATABLE :: psiut(:, :)

!INTEGER,PRIVATE :: kdim=kstate
    INTEGER, PRIVATE :: kdim

! matrices of radial moments

    COMPLEX(DP), ALLOCATABLE :: ExpDABold(:, :, :), wfrotate(:, :, :)
    COMPLEX(DP), ALLOCATABLE, SAVE :: vecs(:, :, :) ! searched eigenvectors

! This PARAMETER 'tnearest' activates the computation of the
! rotation amongst occupuied states. The RESULT ist stored on
! 'wfrotate' and so transferred TO the SIC routines.
    LOGICAL, PARAMETER :: tnearest = .true.

! The switch 'texpo' activates the "extrapolation" of the unitary
! transformation 'vecs' usign the array 'expdabold'. This array is
! initialzed as UNIT matrix and propagated IN parallel TO 'vecs'.
! This strategy seems TO help. But it is not so obvious why.
    LOGICAL, PARAMETER :: texpo = .true.

#endif

  CONTAINS

#ifdef REALSWITCH
!-----init_fsicr------------------------------------------------

    SUBROUTINE init_fsicr()

! Iinitializes fields for FSIC etc

      USE params
      USE kinetic
      IMPLICIT NONE

      INTEGER :: i, j, is, na
      INTEGER :: ifcmplxin, ni, nim, nip ! oly used IN COMPLEX SIC
      REAL(DP) :: vecexist
      COMPLEX(DP) :: ccr, csi

!
! NAMELIST FSIC2 CONTAINS numerical parameters for the solution
! of the symmetry condition:
! step = step SIZE (actually 'step/radmaxysm' is the step SIZE)
! precis = final limit on PRECISION,
! early stages USE actual variance*1D-4 as limit
! symutbegin = NUMBER of iterations IN symmetry condition
! toptsicstep = switch TO tuning of stepsize by quadratic extrapolation
! of energy TO a maximum using the last three steps;
! this option is presently ONLY relevant for the static CASE;
! the following parameters are relevant ONLY
! IF this option is set TO .true.
! precisfact = factor applied TO variance TO compute actual PRECISION limit
! dampopt = damping factor on the optimal step (DEFAULT 0.7)
! steplow = minimum step SIZE (tentative DEFAULT 0.01)
! steplim = maximum step SIZE (tentative DEFAULT 1.2)
!
      NAMELIST /fsic2/ step, precis, symutbegin, precisfact, &
        dampopt, steplow, steplim, toptsicstep, phiini !!! UT parameters

!----------------------------------------------------------------

      IF (ifsicp < 6) RETURN

! READ SIC specific parameters

!      OPEN (5, STATUS='old', FORM='FORMATTED', FILE='for005.'//outnam)
      READ (5, fsic2)
      WRITE (6, '(a,4(1pg13.5))') &
        ' SIC running with step,precis,SymUtBegin=', &
        step, precis, symutbegin
      CLOSE (5)

! dimensions of spin sub-spaces

      ndims(1) = 0
      ndims(2) = 0
      DO na = 1, nstate
        IF (ispin(nrel2abs(na)) == 1) THEN
          ndims(1) = na
        ELSE IF (ispin(nrel2abs(na)) == 2) THEN
          ndims(2) = na - ndims(1)
        END IF
      END DO
      WRITE (6, '(a,2i5)') ' DIMENSION of sub-matrices:', ndims

! initialize work space
      IF (ifsicp == 8) THEN
        ALLOCATE (qnewr(kdfull2, kstate))
        ALLOCATE (psirut(kdfull2, kstate))
        ALLOCATE (rExpDABold(kstate, kstate, 2))
        ALLOCATE (vecsr(kstate, kstate, 2)) ! searched eigenvectors
      ELSE IF (ifsicp == 9) THEN
        ALLOCATE (qnewc(kdfull2, kstate))
        ALLOCATE (psicut(kdfull2, kstate))
        ALLOCATE (rExpDABoldc(kstate, kstate, 2))
        ALLOCATE (vecsc(kstate, kstate, 2)) ! searched eigenvectors
      END IF

      kdim = kstate

! initialize unitary matrix

      IF (ifsicp == 8) THEN
        vecexist = SUM(ABS(vecsr))
      ELSE IF(ifsicp == 9) THEN
        vecexist = SUM(ABS(vecsc))
      END IF
      IF (.NOT.tstat .OR. vecexist<1D-30) THEN
        DO is = 1, 2
          IF (ifsicp == 8) THEN
            DO i = 1, ndims(is)
            DO j = 1, ndims(is)
              IF (i == j) THEN
                vecsr(i, j, is) = 1D0
              ELSE
                vecsr(i, j, is) = 0D0
              END IF
            END DO
            END DO
          ELSE IF (ifsicp == 9) THEN
            DO i = 1, ndims(is)
            DO j = 1, ndims(is)
              IF (i == j) THEN
                vecsc(i, j, is) = CMPLX(1D0, 0D0)
              ELSE
                vecsc(i, j, is) = 0D0
              END IF
            END DO
            END DO
            ifcmplxin = 2
            IF (ifcmplxin == 1) THEN
              ni = ndims(is)
              nim = ni - 1
              vecsc(ni, ni, is) = CMPLX(COS(phiini), 0D0, DP)
              vecsc(nim, nim, is) = vecsr(ni, ni, is)
              vecsc(ni, nim, is) = CMPLX(0D0, SIN(phiini), DP)
              vecsc(nim, ni, is) = CMPLX(0D0, -SIN(phiini), DP)
            ELSE IF (ifcmplxin == 2) THEN
              csi = CMPLX(0D0, SIN(phiini)/SQRT(2D0), DP)
              ccr = CMPLX(SQRT(1D0 - 2D0*ABS(csi)**2), 0D0, DP)
              DO ni = 1, ndims(is)
                vecsc(ni, ni, is) = ccr
                nip = MOD(ni, ndims(is)) + 1
                vecsc(ni, nip, is) = csi
                nim = ni - 1
                IF (nim == 0) nim = ndims(is)
                vecsc(ni, nim, is) = -csi
              END DO
              CALL orthnorm(vecsc(:, :, is), ndims(is), kdim)
            END IF
          END IF
        END DO
      END IF

! initialize protocol FILE
      IF (tconv) THEN
        OPEN (353, FILE='2st-STAT-conv-1.res')
        REWIND (353)
        WRITE (353, '(a)') '# protocol of convergence of symmetry condition for ispin=1'
        CLOSE (353)
        OPEN (353, FILE='2st-STAT-conv-2.res')
        REWIND (353)
        WRITE (353, '(a)') '# protocol of convergence of symmetry condition for ispin=2'
        CLOSE (353)
      END IF

      RETURN
    END SUBROUTINE init_fsicr

!-----end_fsicr------------------------------------------------
    SUBROUTINE end_fsicr()

! Frees workspace used for static SIC

      USE params
      USE kinetic
      IMPLICIT NONE

      IF (ifsicp == 8) THEN
        DEALLOCATE (rExpDABold)
        DEALLOCATE (vecsr)
        DEALLOCATE (qnewr)
        DEALLOCATE (psirut)
      ELSE IF (ifsicp == 9) THEN
        DEALLOCATE (rExpDABoldc)
        DEALLOCATE (vecsc)
        DEALLOCATE (qnewc)
        DEALLOCATE (psicut)
      END IF

      RETURN
    END SUBROUTINE end_fsicr

#endif

#ifdef COMPLEXSWITCH
!-----init_fsic------------------------------------------------

    SUBROUTINE init_fsic()

! Initializes fields for dynamical FSIC etc

      USE params
      USE kinetic
      IMPLICIT NONE

!NAMELIST /fsic/step,precis,symutbegin !!! UT parameters

!----------------------------------------------------------------

      kdim = kstate
      ALLOCATE (ExpDABold(kstate, kstate, 2), wfrotate(kstate, kstate, 2))
      ALLOCATE (vecs(kstate, kstate, 2)) ! searched eigenvectors

      ALLOCATE (qnewut(kdfull2, kstate))
      ALLOCATE (psiut(kdfull2, kstate))

      IF (tconv) THEN
        OPEN (353, FILE='2st-STAT-conv-1.res')
        REWIND (353)
        WRITE (353, '(a)') '# protocol of convergence of symmetry condition for ispin=1'
        CLOSE (353)
        OPEN (353, FILE='2st-STAT-conv-2.res')
        REWIND (353)
        WRITE (353, '(a)') '# protocol of convergence of symmetry condition for ispin=2'
        CLOSE (353)
      END IF

      RETURN
    END SUBROUTINE init_fsic

!-----ExpDabvol_rotatewf_init----------------------------------

    SUBROUTINE expdabvol_rotate_init

! Initializes matrices used later IN rotation.

      IMPLICIT NONE
      INTEGER::is
      DO is = 1, 2
        CALL MatUnite(ExpDabOld(:, :, is), kstate, ndims(is))
        CALL MatUnite(wfrotate(:, :, is), kstate, ndims(is))
      END DO
    END SUBROUTINE expdabvol_rotate_init

!-----end_fsic------------------------------------------------

    SUBROUTINE end_fsic()

! Terminates fields used for FSIC etc

      USE params
      USE kinetic
      IMPLICIT NONE

      kdim = kstate
      DEALLOCATE (ExpDABold, wfrotate, vecs)
      DEALLOCATE (qnewut, psiut)

      IF (tconv) CLOSE (353)

      RETURN
    END SUBROUTINE end_fsic

    SUBROUTINE eval_unitrot(qact, qold)

! Computes the piece of rotation amongst occupied states contained
! IN the tansformation from 'qold' TO 'qact'. The resulting unitary
! transformation is transferred via 'wfrotate'.

      USE params
      USE util, ONLY: wfovlp
!USE twost
      USE orthmat, ONLY: matdorth
      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN) :: qact(kdfull2, kstate)
      COMPLEX(DP), INTENT(IN) :: qold(kdfull2, kstate)

      COMPLEX(DP) :: ovl
      REAL(DP) :: rmo
      INTEGER :: is, ni, na, nb, naeff, nbeff

!----------------------------------------------------------------------

      wfrotate = 0D0
      DO is = 1, 2
        ni = ndims(is)

! evaluate matrix from overlaps
        DO na = 1, nstate
          IF (ispin(nrel2abs(na)) == is) THEN
            naeff = na - (is - 1)*ndims(1)
            DO nb = 1, nstate
              IF (ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
                nbeff = nb - (is - 1)*ndims(1)
                wfrotate(naeff, nbeff, is) = wfovlp(qold(:, na), qact(:, nb))
              END IF
            END DO
          END IF
        END DO

! ortho-normalize
        DO nb = 1, ni
          DO na = 1, nb - 1
            ovl = SUM(CONJG(wfrotate(1:ni, na, is))*wfrotate(1:ni, nb, is))
            wfrotate(:, nb, is) = wfrotate(:, nb, is) - wfrotate(:, na, is)*ovl
          END DO
          ovl = SUM(CONJG(wfrotate(1:ni, nb, is))*wfrotate(1:ni, nb, is))
          wfrotate(1:ni, nb, is) = wfrotate(1:ni, nb, is)*CMPLX(1D0/SQRT(REAL(ovl, DP)), 0D0, DP)
        END DO

! transpose for further USE
        wfrotate(1:ni, 1:ni, is) = CONJG(TRANSPOSE(wfrotate(1:ni, 1:ni, is)))

        rmo = matdorth(wfrotate(:, :, is), kstate, ni)
        WRITE (*, *) 'is,unitarity wfrotate:', is, rmo
        IF (ABS(rmo) > 1D-10) WRITE (*, *) ' WFROTATE:', wfrotate(1:ni, 1:ni, is)

      END DO

    END SUBROUTINE eval_unitrot

#endif

! ******************************

#ifdef COMPLEXSWITCH
!-----init_vecs-------------------------------------------------

    SUBROUTINE init_vecs()

! Initializes COMPLEX superposition coefficients for SIC.

      USE params
      USE kinetic
      IMPLICIT NONE

      INTEGER :: i, ii, j, jj
!----------------------------------------------------------------

      WRITE (6, *) 'vecsr initvecs'
      WRITE (6, '(6f12.6)') ((vecsr(ii, jj, 1), ii=1, 3), jj=1, 3)

      DO i = 1, kstate
        DO j = 1, kstate
          vecs(i, j, 1) = vecsr(i, j, 1)
          vecs(i, j, 2) = vecsr(i, j, 2)
        END DO
      END DO

      WRITE (6, *) 'vecs fin initvecs'
      WRITE (6, '(6f12.6)') ((vecs(ii, jj, 1), ii=1, 3), jj=1, 3)

      RETURN
    END SUBROUTINE init_vecs
#endif

#ifdef REALSWITCH

!-----static_sicfield------------------------------------------------

    SUBROUTINE static_sicfield(rho, aloc, psir, iter1)

! The Coulomb part of the SIC mean field.
!
! Input:
! rho = electron density
! psir = REAL wavefunctions
! iter1 = nr. of iteration or time steps
! Output:
! aloc = local mean-field potential

      USE params
      USE kinetic
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: rho(2*kdfull2)
      REAL(DP), INTENT(OUT) :: aloc(2*kdfull2)
      REAL(DP), INTENT(IN) :: psir(kdfull2, kstate)
      INTEGER, INTENT(IN) :: iter1

!----------------------------------------------------------------

      IF (ifsicp < 7) RETURN

! computation of the mean fields

      IF (ifsicp == 8) THEN
        CALL calc_utwfr(psir, psirut, iter1)
      ELSE IF (ifsicp == 9) THEN
        CALL calc_utwfrc(psir, psicut, iter1)
      END IF

      IF (ifsicp == 8) THEN
        CALL calc_fullsicr(psirut, qnewr)
      ELSE IF (ifsicp == 9) THEN
        CALL calc_fullsic(psicut, qnewc)
      END IF

      RETURN
    END SUBROUTINE static_sicfield

!-----infor_sic------------------------------------------------

    SUBROUTINE infor_sic(psir)

! Compute and PRINT variances, static CASE.
!
! Input:
! psir = set of s.p. wavefunctions

      USE params
      USE kinetic
      USE util, ONLY: wfovlp
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: psir(kdfull2, kstate)

      INTEGER :: is, na, nb
      COMPLEX(DP) :: cacc
      REAL(DP) :: acc

!----------------------------------------------------------------

      IF (ifsicp < 7) RETURN

! TESTS

      WRITE (6, '(a)') 'DIAGONAL STATES :'
      CALL spmomsmatrix(psir, 1) !!! TO PRINT the total variance

      IF (ifsicp == 8) THEN !!! TO calculate the total
        WRITE (6, '(a)') 'LOCALIZED STATES :'
        CALL spmomsmatrix(psirut, 1)
        DO is = 1, 2 !!! violation of symmetry condition
          acc = 0D0
          DO na = 1, nstate
          IF (ispin(na) == is .AND. occup(na) > 0.5D0) THEN
          DO nb = 1, nstate
          IF (ispin(nb) == is .AND. occup(nb) > 0.5D0) THEN
            acc = (wfovlp(psirut(:, na), qnewr(:, nb)) - &
                   wfovlp(psirut(:, nb), qnewr(:, na)))**2 + acc
          END IF
          END DO
          END IF
          END DO
          WRITE (6, '(a,i3,a,2(1pg12.4))') &
            'For spin', is, ' Total violation of SymCond', SQRT(acc)
        END DO
      ELSE IF (ifsicp == 9) THEN !!! TO calculate the total
        WRITE (6, '(a)') 'LOCALIZED STATES :'
        CALL spmomsmatrix(psicut, 1)
        DO is = 1, 2 !!! violation of symmetry condition
          cacc = 0D0
          DO na = 1, nstate
            IF (ispin(na) == is) THEN
              DO nb = 1, nstate
                IF (ispin(nb) == is) THEN
                  cacc = (wfovlp(psicut(:, na), qnewc(:, nb)) - &
                          wfovlp(psicut(:, nb), qnewc(:, na)))**2 + cacc
                END IF
              END DO
            END IF
          END DO
          WRITE (6, '(a,i3,a,2(1pg12.4))') &
            'For spin', is, ' Total violation of SymCond', SQRT(cacc)
        END DO
      END IF

      RETURN
    END SUBROUTINE infor_sic

!-----diag_lagr------------------------------------------------

    SUBROUTINE diag_lagr(psir)

! Diagonalize the matrix of Lagrange parameters and PRINT
! information from that step.
!
! Input/Output:
! psir = set of s.p. wavefunctions

      USE params
      USE kinetic
!USE twost_util, ONLY:superpose_state
      IMPLICIT NONE

      REAL(DP), INTENT(IN OUT) :: psir(kdfull2, kstate)

! vect1,2 must be declared IN each spin subspace for correct
! diagonalization of SIC with multipliers

      INTEGER :: i, j, is, na, nb, nbeff
      REAL(DP) :: a(nstate*(nstate + 1)/2), root(kdim, 2), &
                  vect1(nstate - nspdw, nstate - nspdw), vect2(nspdw, nspdw)

!----------------------------------------------------------------

      IF (ifsicp < 7) RETURN

      DO is = 1, 2 !!! Diagonalization on the Lagrange matrix for FSIC

        i = 1
        DO na = 1, nstate
          IF (ispin(nrel2abs(na)) == is) THEN
            DO nb = 1, na
              IF (ispin(nrel2abs(nb)) == is) THEN
                a(i) = hmatrix(nb, na)
                i = i + 1
              END IF
            END DO
          END IF
        END DO

! Lagrange mult matrix diag (needs vect1,2 IN each sp subs, instead bugs)
        IF (is == 1) THEN
          CALL givens(a, root(:, is), vect1(:, :), ndims(is), ndims(is), ndims(is))
        ELSE IF (is == 2) THEN
          CALL givens(a, root(:, is), vect2(:, :), ndims(is), ndims(is), ndims(is))
        END IF

        DO na = 1, ndims(is)
          DO nb = 1, ndims(is)
            IF (is == 1) THEN
              vecsr(na, nb, is) = vect1(nb, na)
            ELSE IF (is == 2) THEN
              vecsr(na, nb, is) = vect2(nb, na)
            END IF
          END DO
        END DO

        DO nb = 1, nstate
          IF (is == ispin(nrel2abs(nb))) THEN
            nbeff = nb - (is - 1)*ndims(1)
            IF (is == 1) THEN
              CALL superpose_state(psir(:, nb), vect1(1:kstate, nbeff), psirut, is)
            ELSE IF (is == 2) THEN
              CALL superpose_state(psir(:, nb), vect2(1:kstate, nbeff), psirut, is)
            END IF
          END IF
        END DO

        WRITE (6, *) 'Diag of the Lagrange matrix for spin subs', is
        DO j = 1, ndims(is)
          WRITE (6, *) 'state', j, 'single energies', root(j, is)
        END DO

      END DO

      WRITE (6, *) 'DIAGONAL STATES :'
      CALL spmomsmatrix(psir, 1) !!! TO PRINT the total variance
      WRITE (6, *) 'LOCALIZED STATES :'
      CALL spmomsmatrix(psirut, 1)

      RETURN
    END SUBROUTINE diag_lagr

!-----subtr_sicpot------------------------------------------------

    SUBROUTINE subtr_sicpot(q1, nbe)

! Subtract the SIC potential acting on the actual 'q1'.
! COMPLEX version.
!
! Input:
! nbe = nr. of state
! Input/Output:
! q1 = wavefunction TO be modified

      USE params
      USE kinetic
      IMPLICIT NONE
      REAL(DP), INTENT(IN OUT) :: q1(kdfull2)
      INTEGER, INTENT(IN) :: nbe

      INTEGER :: is, na, nb, nae
!----------------------------------------------------------------

      IF (ifsicp < 8) RETURN

      is = ispin(nbe)
      nb = nbe - (is - 1)*ndims(1)
      IF (ifsicp == 8) THEN
        DO na = 1, ndims(is)
          nae = na + (is - 1)*ndims(1)
          q1(:) = q1(:) - qnewr(:, nae)*vecsr(nb, na, is)
        END DO
      ELSE IF (ifsicp == 9) THEN
        DO na = 1, ndims(is)
          nae = na + (is - 1)*ndims(1)
          q1(:) = q1(:) - qnewc(:, nae)*CONJG(vecsc(nb, na, is))
        END DO
      END IF
      RETURN
    END SUBROUTINE subtr_sicpot

#endif

#ifdef COMPLEXSWITCH

    SUBROUTINE calc_utwfc(q0, q0ut, iter1)

! Computes localized wavefunctions. COMPLEX version.
!
! Input:
! q0 = set of s.p. wavefunctions
! iter1 = nr. of time step
! Output:
! q0ut = set of localized wavefunctions

      USE params
      USE kinetic
      USE util, ONLY: wfovlp
      USE twost_util, ONLY: superpose_state
      IMPLICIT NONE

      COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
      COMPLEX(DP), INTENT(OUT) :: q0ut(kdfull2, kstate)
      INTEGER, INTENT(IN) :: iter1

      INTEGER :: is, nb, nbeff

!------------------------------------------------------------------

! Compute matrix of radial moments and determine
! optimally localizing transformation.
! Results is unitary matrix 'vecs' communicated via
! COMMON /radmatrix/.

      DO is = 1, 2
        IF (ndims(is) > 1) THEN
          CALL utgradstepc(is, 0, q0, iter1)
        END IF
      END DO

      DO nb = 1, nstate
        is = ispin(nrel2abs(nb))
        nbeff = nb - (is - 1)*ndims(1)
        CALL superpose_state(q0ut(:, nb), vecs(:, nbeff, is), q0, is)
      END DO

      RETURN
    END SUBROUTINE calc_utwfc
#endif

!-----calc_utwf--------------------------------------------------

#ifdef REALSWITCH
    SUBROUTINE calc_utwfr(q0, q0ut, iter1)

! Computes localized wavefunctions. REAL version.
!
! Input:
! q0 = set of s.p. wavefunctions
! iter1 = nr. of time step
! Output:
! q0ut = set of localized wavefunctions

      USE params
      USE kinetic
      USE twost_util, ONLY: superpose_state

      IMPLICIT NONE

      INTEGER :: is, nb, nbeff
! basic arrays and workspace

!#ifdef REALSWITCH
      REAL(DP) :: q0(kdfull2, kstate)
      REAL(DP) :: q0ut(kdfull2, kstate)
      INTEGER, INTENT(IN) :: iter1

!------------------------------------------------------------------

! Compute matrix of radial moments and determine
! optimally localizing transformation.
! Results is unitary matrix 'vecs' communicated via
! COMMON /radmatrix/.

      DO is = 1, 2
        IF (ndims(is) > 1) THEN
          CALL utgradstepr(is, 0, q0, iter1)
        END IF
      END DO

      DO nb = 1, nstate
        is = ispin(nrel2abs(nb))
        nbeff = nb - (is - 1)*ndims(1)
        CALL superpose_state(q0ut(:, nb), vecsr(:, nbeff, is), q0, is)
      END DO

      RETURN
    END SUBROUTINE calc_utwfr

    SUBROUTINE calc_utwfrc(q0, q0ut, iter1)

! Computes localized wavefunctions. REAL/COMPLEX version.
!
! Input:
! q0 = set of s.p. wavefunctions
! iter1 = nr. of time step
! Output:
! q0ut = set of localized wavefunctions

      USE params
      USE kinetic
      USE twost_util, ONLY: superpose_state

      IMPLICIT NONE

      INTEGER :: is, nb, nbeff
! basic arrays and workspace

      REAL(DP) :: q0(kdfull2, kstate)
      COMPLEX(DP) :: q0ut(kdfull2, kstate)
      INTEGER, INTENT(IN) :: iter1

!------------------------------------------------------------------

! Compute matrix of radial moments and determine
! optimally localizing transformation.
! Results is unitary matrix 'vecs' communicated via
! COMMON /radmatrix/.

      DO is = 1, 2
        IF (ndims(is) > 1) THEN
          CALL utgradsteprc(is, 0, q0, iter1) !!!!! new vecs (gradient method)
        END IF
      END DO

      DO nb = 1, nstate
        is = ispin(nrel2abs(nb))
        nbeff = nb - (is - 1)*ndims(1)
        CALL superpose_state(q0ut(:, nb), vecsc(:, nbeff, is), q0, is)
      END DO

      RETURN
    END SUBROUTINE calc_utwfrc
#endif

!-----utgradstepc------------------------------------------------------------------------

#ifdef COMPLEXSWITCH
    SUBROUTINE utgradstepc(is, iprint, q0, iter1)

! Nonlinear gradient iteration TO optimally localized states.
! COMPLEX version.
!
! Input:
! q0 = set of s.p. wavefunctions
! is = isospin
! iprint = PRINT level: <0 --> no PRINT at all
! 0 --> ONLY final RESULT
! >0 --> PRINT modulus
! iter1 = nr. of time step

      USE params
      USE kinetic
      USE twostr, ONLY: dalphabeta
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: is
      INTEGER, INTENT(IN) :: iprint
      COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
      INTEGER, INTENT(IN) :: iter1

      INTEGER :: iter
      REAL(DP) :: e1der, e2der, e1old, enold_2st, ERR_c, stepnew, steplow
      COMPLEX(DP) :: xlambda(kdim, kdim)
      COMPLEX(DP) :: acc
      COMPLEX(DP) :: dab(kdim, kdim), expdab(kdim, kdim)
      COMPLEX(DP) :: dabsto(kdim, kdim)

      INTEGER :: itmax2, ni
      LOGICAL, PARAMETER :: ttest = .false.

      REAL(DP) :: norm ! variance of step, erreur du produit
      REAL(DP) :: actstep
      REAL(DP) :: enstore(3), stepstore(3) ! storage for optimized step

!-------------------------------------------------------

      itmax2 = symutbegin
      ni = ndims(is)

! update unitary transformation 'vecs' with previous exponentional
      IF (tnearest) vecs(1:ni, 1:ni, is) = &
        MATMUL(wfrotate(1:ni, 1:ni, is), vecs(1:ni, 1:ni, is))

      IF (texpo) THEN
        vecs(1:ni, 1:ni, is) = MATMUL(vecs(1:ni, 1:ni, is), expdabold(1:ni, 1:ni, is))
      END IF

      actstep = step

      IF (tconv) THEN
        IF (is == 1) OPEN (353, FILE='2st-STAT-conv-1.res', POSITION='append')
        IF (is == 2) OPEN (353, FILE='2st-STAT-conv-2.res', POSITION='append')
        WRITE (353, *) '# convergence symmetry condition. timestep,is=', iter1, is
        WRITE (353, '(a)') '# iter Ortho , variance, erreur , actstep'
      ELSE
        WRITE (6, '(a)') ' iter Ortho , variance, non-linearity, actstep'
      END IF

      enold_2st = 0D0
      DO iter = 1, itmax2
        CALL dalphabeta(is, dab, vecs, q0) !actual symmetry condition on 'dab'
        IF (toptsicstep) THEN
          IF (iter .LE. 3) THEN
            enstore(iter) = ener_2st(is)
            IF (iter == 1) THEN
              stepstore(1) = 0D0
            ELSE
              stepstore(iter) = stepstore(iter - 1) + actstep
            END IF
          ELSE
            enstore(1) = enstore(2)
            enstore(2) = enstore(3)
            enstore(3) = ener_2st(is)
            stepstore(1) = stepstore(2)
            stepstore(2) = stepstore(3)
            stepstore(3) = stepstore(3) + actstep
            e1der = (enstore(3) - enstore(2))/(stepstore(3) - stepstore(2))
            e1old = (enstore(2) - enstore(1))/(stepstore(2) - stepstore(1))
            e2der = (e1der - e1old)*2D0/(stepstore(3) - stepstore(1))
            stepnew = -dampopt*e1der/e2der
            IF (stepnew < 0D0) stepnew = step
            actstep = max(steplow, min(stepnew, steplim))
            IF (ttest) THEN
              WRITE (*, *) ' stepnew,steplim=', stepnew, steplim, actstep
              WRITE (*, '(a,30(1pg15.7))') ' enstore=', enstore
              WRITE (*, '(a,30(1pg15.7))') ' stepstore=', stepstore
              WRITE (6, '(a,3(1pg15.7))') 'e1d,e1o,e2d=', e1der, e1old, e2der
              WRITE (6, '(a,2(1pg15.7))') 'act/newstep=', actstep, stepnew
            END IF
          END IF
        END IF
        norm = SQRT(SUM(dab(1:ni, 1:ni)*CONJG(dab(1:ni, 1:ni))))
        dab(1:ni, 1:ni) = CMPLX(-actstep, 0D0)*dab(1:ni, 1:ni)
        CALL matexp(dab, expdab, kdim, ni) ! exponentiate 'dab' TO 'ExpDab'
        vecs(1:ni, 1:ni, is) = MATMUL(vecs(1:ni, 1:ni, is), expdab(1:ni, 1:ni))
        IF (texpo) THEN
          expdabold(1:ni, 1:ni, is) = MATMUL(expdabold(1:ni, 1:ni, is), expdab(1:ni, 1:ni))
          dab(1:ni, 1:ni) = CONJG(expdabold(1:ni, 1:ni, is))
          dabsto(1:ni, 1:ni) = MATMUL(expdabold(1:ni, 1:ni, is), dab(1:ni, 1:ni))
          ERR_c = SQRT(ABS(SUM(dabsto(1:ni, 1:ni)*CONJG(dabsto(1:ni, 1:ni)) - ndims(is))))
        END IF
        IF (tconv) THEN
          WRITE (353, '(i4,6(1pg13.5))') &
            iter, matdorth(vecs(:, :, is), kdim, ndims(is)), ABS(norm), &
            actstep*norm**2/(ener_2st(is) - enold_2st), actstep, &
            ener_2st(is) - enold_2st
          CALL FLUSH (353)

        ELSE
          WRITE (6, '(a,4f16.13)') ' Ortho , variance, erreur, actstep', &
            matdorth(vecs(:, :, is), kdim, ndims(is)), ABS(norm), ABS(ERR_c), actstep
        END IF
        IF (iter .GE. 1) enold_2st = ener_2st(is)
        IF (ABS(norm) < precis) EXIT
      END DO

      IF (tconv) CLOSE (353)
      RETURN
    END SUBROUTINE utgradstepc
#endif

#ifdef REALSWITCH
    SUBROUTINE utgradstepr(is, iprint, q0, iter1)

! Nonlinear gradient iteration TO optimally localized states.
! REAL version
!
! Input:
! q0 = set of s.p. wavefunctions
! is = isospin
! iprint = PRINT level: <0 --> no PRINT at all
! 0 --> ONLY final RESULT
! >0 --> PRINT modulus
! iter1 = nr. of time step

      USE params
      USE kinetic
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
      REAL(DP) :: dab(kdim, kdim), expdab(kdim, kdim)
      REAL(DP) :: dabsto(kdim, kdim)
      INTEGER, INTENT(IN) :: is, iprint, iter1

      INTEGER :: i, iter, itmax2, jj, ni
      LOGICAL, PARAMETER :: ttest = .FALSE.

      REAL(DP) :: e1der, e1old, e2der, ERR_r, enorm, norm

      REAL(DP) :: actstep, stepnew, enold_2st, actprecis
      REAL(DP) :: enstore(3), stepstore(3) ! storage for optimized step

! The optimized step computes a quadratic from for the SIC energy
! from the previous three iterations and extrapolates this FORM
! TO the maximum.
! The step SIZE is limited from above by the PARAMETER 'steplim'
! and from below by 'steplow'.
! IF a negative step SIZE emerges, the initial step SIZE is used.
!

!-------------------------------------------------------

! tentative change of step SIZE during iteration
      IF (iter1 < 100000) THEN
        itmax2 = symutbegin
      ELSE
        itmax2 = 1
      END IF
      WRITE(*,*) 'utgradstepr entered:',iter1,symutbegin,itmax2

      ni = ndims(is)

      IF (ttest) THEN
        WRITE (6, *) 'entree utgradstepr: is=', is
        DO jj = 1, ni
          WRITE (6, '(6(1pg13.5))') vecsr(:, jj, is)
        END DO
      END IF

! update unitary transformation 'vecs' with previous exponentional
      actstep = step ! radmaxsym obsolete, set TO 1D0
      actprecis = max(precis, precisfact*sumvar2)
      IF (tconv) THEN
        IF (is == 1) OPEN (353, FILE='2st-STAT-conv-1.res', POSITION='append')
        IF (is == 2) OPEN (353, FILE='2st-STAT-conv-2.res', POSITION='append')
        WRITE (353, *) '# convergence symmetry condition. iter1,is=', iter1, is
        WRITE (353, '(a)') '# Ortho , variance, erreur , actstep, Spin'
        IF (ttest) THEN
          WRITE (353, *) ' vecsr before: is=', is
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') vecsr(1:ni, i, is)
          END DO
          WRITE (353, *) ' rexpdabold before: is=', is
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') rexpdabold(1:ni, i, is)
          END DO
          WRITE (353, '(a,i4,1pg13.5)') &
            ' Iter,Ortho,variance,erreur,energy. Spin,precis=', is, actprecis
        END IF
        CALL FLUSH (353)
      END IF
      WRITE (6, '(a,i4,1pg13.5)') &
        ' Iter,Ortho,variance,erreur,energy. Spin,precis=', is, actprecis

      enold_2st = 0D0
      DO iter = 1, itmax2
      IF(ttest) WRITE(*,*) 'inside iteration: iter=',iter
        CALL dalphabeta(is, dab, vecsr, q0) !new DAB
        IF (iter == 1) THEN
          WRITE (353, *) ' dab before: is=', is
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') dab(1:ni, i)
          END DO
        END IF
        IF (toptsicstep) THEN
          IF (iter .LE. 3) THEN
            enstore(iter) = ener_2st(is)
            IF (iter == 1) THEN
              stepstore(1) = 0D0
            ELSE
              stepstore(iter) = stepstore(iter - 1) + actstep
            END IF
          ELSE
            enstore(1) = enstore(2)
            enstore(2) = enstore(3)
            enstore(3) = ener_2st(is)
            stepstore(1) = stepstore(2)
            stepstore(2) = stepstore(3)
            stepstore(3) = stepstore(3) + actstep
            e1der = (enstore(3) - enstore(2))/(stepstore(3) - stepstore(2))
            e1old = (enstore(2) - enstore(1))/(stepstore(2) - stepstore(1))
            e2der = (e1der - e1old)*2D0/(stepstore(3) - stepstore(1))
            stepnew = -dampopt*e1der/e2der
            IF (stepnew < 0D0) stepnew = step
            actstep = max(steplow, min(stepnew, steplim))
            IF (ttest) THEN
              WRITE (*, *) ' stepnew,steplim=', stepnew, steplim, actstep
              WRITE (*, '(a,30(1pg15.7))') ' enstore=', enstore
              WRITE (*, '(a,30(1pg15.7))') ' stepstore=', stepstore
              WRITE (6, '(a,3(1pg15.7))') 'e1d,e1o,e2d=', e1der, e1old, e2der
              WRITE (6, '(a,2(1pg15.7))') 'act/newstep=', actstep, stepnew
            END IF
          END IF
        END IF
        norm = SQRT(SUM(dab(1:ni, 1:ni)**2)) ! rmatnorme(dab,kdim,ni)
        dab(1:ni, 1:ni) = -actstep*dab(1:ni, 1:ni)
        CALL matexp(dab, expdab, kdim, ni) ! MV exp IN ExpDab
        rexpdabold(1:ni, 1:ni, is) = MATMUL(rexpdabold(1:ni, 1:ni, is), expdab(1:ni, 1:ni))
        vecsr(1:ni, 1:ni, is) = MATMUL(vecsr(1:ni, 1:ni, is), expdab(1:ni, 1:ni))
        dab(1:ni, 1:ni) = rexpdabold(1:ni, 1:ni, is)
        dabsto(1:ni, 1:ni) = MATMUL(rexpdabold(1:ni, 1:ni, is), dab(1:ni, 1:ni))
        ERR_r = SQRT(ABS(SUM(dab(1:ni, 1:ni)**2) - ndims(is)))
        IF (tconv) THEN
          WRITE (353, '(i4,8(1pg13.5))') &
            iter, matdorth(vecsr(:, :, is), kdim, ndims(is)), &
            norm, ERR_r, ener_2st(is) - enold_2st, actstep, &
            actstep*norm**2/(ener_2st(is) - enold_2st), ener_2st(is)
          CALL FLUSH (353)
        ELSE
          WRITE (6, '(i4,5(1pg13.5))') &
            iter, matdorth(vecsr(:, :, is), kdim, ndims(is)), &
            norm, ERR_r, ener_2st(is) - enold_2st, actstep
          CALL FLUSH (6)
        END IF
        IF (iter .GE. 1) enold_2st = ener_2st(is)

        IF (iter > 10 .AND. ABS(norm) < actprecis) EXIT
!        IF (iter > 0 .AND. ABS(norm) < actprecis) EXIT

      END DO

      IF (tconv) THEN
        IF (ttest) THEN
          WRITE (353, *) ' after: vecsr='
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') vecsr(1:ni, i, is)
          END DO
          WRITE (353, *) ' after: rexpdabold='
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') rexpdabold(1:ni, i, is)
          END DO
        END IF
        WRITE (353, '(1x/1x)')
        CLOSE (353)
      END IF

      IF (iprint >= 0) THEN
        WRITE (6, '(2(a,i4),a,1pg12.4)') 'UT cv is reached at it nr.', iter, &
          ' for spin =', is, ' with variance=', norm
      END IF
      RETURN
    END SUBROUTINE utgradstepr
    SUBROUTINE utgradsteprc(is, iprint, q0, iter1)

! Nonlinear gradient iteration TO optimally localized states.
! REAL/COMPLEX version
!
! Input:
! q0 = set of s.p. wavefunctions
! is = isospin
! iprint = PRINT level: <0 --> no PRINT at all
! 0 --> ONLY final RESULT
! >0 --> PRINT modulus
! iter1 = nr. of time step

      USE params
      USE kinetic
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
      COMPLEX(DP) :: dab(kdim, kdim), expdab(kdim, kdim)
      COMPLEX(DP) :: dabsto(kdim, kdim)
      INTEGER, INTENT(IN) :: is, iprint, iter1

      INTEGER :: i, iter, itmax2, jj, ni
      LOGICAL, PARAMETER :: ttest = .false.

      REAL(DP) :: e1der, e1old, e2der, ERR_r, enorm, norm

      REAL(DP) :: actstep, stepnew, enold_2st, actprecis
      REAL(DP) :: enstore(3), stepstore(3) ! storage for optimized step

! The optimized step computes a quadratic from for the SIC energy
! from the previous three iterations and extrapolates this FORM
! TO the maximum.
! The step SIZE is limited from above by the PARAMETER 'steplim'
! and from below by 'steplow'.
! IF a negative step SIZE emerges, the initial step SIZE is used.
!

!-------------------------------------------------------

! tentative change of step SIZE during iteration
      IF (iter1 < 100000) THEN
        itmax2 = symutbegin
      ELSE
        itmax2 = 1
      END IF

      ni = ndims(is)

      IF (ttest) THEN
        WRITE (6, *) 'entree utgradstepr: is=', is
        DO jj = 1, ni
          WRITE (6, '(6(1pg13.5))') vecsc(:, jj, is)
        END DO
      END IF

! update unitary transformation 'vecs' with previous exponentional
      actstep = step ! radmaxsym obsolete, set TO 1D0
      actprecis = max(precis, precisfact*sumvar2)
      IF (tconv) THEN
        IF (is == 1) OPEN (353, FILE='2st-STAT-conv-1.res', POSITION='append')
        IF (is == 2) OPEN (353, FILE='2st-STAT-conv-2.res', POSITION='append')
        WRITE (353, *) '# convergence symmetry condition. iter1,is=', iter1, is
        WRITE (353, '(a)') '# Ortho , variance, erreur , actstep, Spin'
        IF (ttest) THEN
          WRITE (353, *) ' vecsr before: is=', is
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') vecsc(1:ni, i, is)
          END DO
          WRITE (353, *) ' rexpdabold before: is=', is
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') rexpdaboldc(1:ni, i, is)
          END DO
          WRITE (353, '(a,i4,1pg13.5)') &
            ' Iter,Ortho,variance,erreur,energy. Spin,precis=', is, actprecis
        END IF
        CALL FLUSH (353)
      END IF
      WRITE (6, '(a,i4,1pg13.5)') &
        ' Iter,Ortho,variance,erreur,energy. Spin,precis=', is, actprecis

      enold_2st = 0D0
      DO iter = 1, itmax2

        CALL dalphabeta(is, dab, vecsc, q0) !new DAB
        IF (iter == 1) THEN
          WRITE (353, *) ' dab before: is=', is
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') dab(1:ni, i)
          END DO
        END IF
        IF (toptsicstep) THEN
          IF (iter .LE. 3) THEN
            enstore(iter) = ener_2st(is)
            IF (iter == 1) THEN
              stepstore(1) = 0D0
            ELSE
              stepstore(iter) = stepstore(iter - 1) + actstep
            END IF
          ELSE
            enstore(1) = enstore(2)
            enstore(2) = enstore(3)
            enstore(3) = ener_2st(is)
            stepstore(1) = stepstore(2)
            stepstore(2) = stepstore(3)
            stepstore(3) = stepstore(3) + actstep
            e1der = (enstore(3) - enstore(2))/(stepstore(3) - stepstore(2))
            e1old = (enstore(2) - enstore(1))/(stepstore(2) - stepstore(1))
            e2der = (e1der - e1old)*2D0/(stepstore(3) - stepstore(1))
            stepnew = -dampopt*e1der/e2der
            IF (stepnew < 0D0) stepnew = step
            actstep = max(steplow, min(stepnew, steplim))
            IF (ttest) THEN
              WRITE (*, *) ' stepnew,steplim=', stepnew, steplim, actstep
              WRITE (*, '(a,30(1pg15.7))') ' enstore=', enstore
              WRITE (*, '(a,30(1pg15.7))') ' stepstore=', stepstore
              WRITE (6, '(a,3(1pg15.7))') 'e1d,e1o,e2d=', e1der, e1old, e2der
              WRITE (6, '(a,2(1pg15.7))') 'act/newstep=', actstep, stepnew
            END IF
          END IF
        END IF
        norm = SQRT(SUM(dab(1:ni, 1:ni)*CONJG(dab(1:ni, 1:ni)))) ! rmatnorme(dab,kdim,ni)
        dab(1:ni, 1:ni) = -actstep*dab(1:ni, 1:ni)
        CALL matexp(dab, expdab, kdim, ni) ! MV exp IN ExpDab
        rexpdaboldc(1:ni, 1:ni, is) = MATMUL(rexpdaboldc(1:ni, 1:ni, is), expdab(1:ni, 1:ni))
        vecsc(1:ni, 1:ni, is) = MATMUL(vecsc(1:ni, 1:ni, is), expdab(1:ni, 1:ni))
        dab(1:ni, 1:ni) = CONJG(rexpdaboldc(1:ni, 1:ni, is))
        dabsto(1:ni, 1:ni) = MATMUL(rexpdaboldc(1:ni, 1:ni, is), dab(1:ni, 1:ni))
        ERR_r = SQRT(ABS(SUM(dab(1:ni, 1:ni)*CONJG(dab(1:ni, 1:ni)) - ndims(is))))
        IF (tconv) THEN
          WRITE (353, '(i4,8(1pg13.5))') &
            iter, matdorth(vecsc(:, :, is), kdim, ndims(is)), &
            norm, ERR_r, ener_2st(is) - enold_2st, actstep, &
            actstep*norm**2/(ener_2st(is) - enold_2st), ener_2st(is)
          CALL FLUSH (353)
        ELSE
          WRITE (6, '(i4,5(1pg13.5))') &
            iter, matdorth(vecsc(:, :, is), kdim, ndims(is)), &
            norm, ERR_r, ener_2st(is) - enold_2st, actstep
          CALL FLUSH (6)
        END IF
        IF (iter .GE. 1) enold_2st = ener_2st(is)

        IF (iter > 0 .AND. ABS(norm) < actprecis) EXIT

      END DO

      IF (tconv) THEN
        IF (ttest) THEN
          WRITE (353, *) ' after: vecsr='
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') vecsc(1:ni, i, is)
          END DO
          WRITE (353, *) ' after: rexpdabold='
          DO i = 1, ni
            WRITE (353, '(20(1pg13.5))') rexpdaboldc(1:ni, i, is)
          END DO
        END IF
        WRITE (353, '(1x/1x)')
        CLOSE (353)
      END IF

      IF (iprint >= 0) THEN
        WRITE (6, '(2(a,i4),a,1pg12.4)') 'UT cv is reached at it nr.', iter, &
          ' for spin =', is, ' with variance=', norm
      END IF
      RETURN
    END SUBROUTINE utgradsteprc

    SUBROUTINE test_symmcond(is, vecact, q0)

! Tests the symmetry condition.
!
! Input:
! is = spin
! vecact = actual transformation vector
! q0 = set of s.p. wavefunctions

      USE params
      USE kinetic
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
      REAL(DP), INTENT(IN) :: vecact(kdim, kdim)
      INTEGER, INTENT(IN) :: is

      LOGICAL, PARAMETER :: ttest = .false.

      REAL(DP) :: actstep, enold_2st, norm, enorm
      INTEGER :: ni, i

      REAL(DP), PARAMETER :: dactstep = 0.1D0
      INTEGER, PARAMETER :: numteststep = 40

      REAL(DP), ALLOCATABLE :: vecsav(:, :), dabstep(:, :), dab(:, :), expdab(:, :)

!-------------------------------------------------------

      ALLOCATE (vecsav(kdim, kdim), dabstep(kdim, kdim), dab(kdim, kdim), expdab(kdim, kdim))

      ni = ndims(is)
      vecsav(1:ni, 1:ni) = vecsr(1:ni, 1:ni, is)
      CALL dalphabeta(is, dab, vecsr, q0)
      norm = SQRT(SUM(dab(1:ni, 1:ni)**2))
      enold_2st = ener_2st(is)

      WRITE (353, '(a,i2)') ' test linearity of unitary transformation for is=', is
      WRITE (353, '(a)') ' step-SIZE energy old-energy energy-ratio '

      DO i = 0, numteststep
        actstep = MAX(i*dactstep, 2D-3)
        dabstep(1:ni, 1:ni) = -actstep*dab(1:ni, 1:ni)
        CALL matexp(dabstep, expdab, kdim, ni)
        vecsr(1:ni, 1:ni, is) = MATMUL(vecsav(1:ni, 1:ni), expdab(1:ni, 1:ni))
        CALL dalphabeta(is, dabstep, vecsr, q0)
        WRITE (353, '(f8.3,6(1pg13.5))') &
          actstep, ener_2st(is), ener_2st(is) - enold_2st, &
          actstep*norm**2/(ener_2st(is) - enold_2st), actstep*norm**2
      END DO

      vecsr(1:ni, 1:ni, is) = vecsav(1:ni, 1:ni)

      DEALLOCATE (vecsav, dabstep, dab, expdab)

    END SUBROUTINE test_symmcond
#endif

#ifdef REALSWITCH
    SUBROUTINE dalphabeta_rc(is, dab, vec, q0)

! Compose matrix 'dab=<psi_a|h_a-h_b|psi_b>. REAL/COMPLEX version.
!
! Input:
! vec = actual states IN configuration space
! q0 = set of s.p. wavefunctions IN coordinate space
! is = spin
! Output:
! dab = matrix of symmetry condition

      USE params
      USE kinetic
      USE util, ONLY: wfovlp
      USE twost_util, ONLY: superpose_state
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: is
      COMPLEX(DP), INTENT(OUT) :: dab(kdim, kdim) ! DAB
      COMPLEX(DP), INTENT(IN) :: vec(kstate, kstate, 2)
      REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)

      LOGICAL, PARAMETER :: ttest = .false.

      INTEGER :: ind, ishift, na, naa, nb, nbeff
      REAL(DP) :: save1, save2, acc2

      REAL(DP), ALLOCATABLE :: usicsp(:), rhosp(:)
      COMPLEX(DP), ALLOCATABLE :: qsym(:, :), utcond(:, :)
      COMPLEX(DP), ALLOCATABLE :: uqsym(:, :), symcond(:, :)

!-------------------------------------------------------

      ALLOCATE (usicsp(2*kdfull2), rhosp(2*kdfull2))
      ALLOCATE (qsym(kdfull2, kstate), utcond(kdim, kdim))
      ALLOCATE (uqsym(kdfull2, kstate), symcond(kstate, kstate))

      ishift = (is - 1)*nxyz

      ener_2st(is) = 0D0
      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)
          CALL superpose_state(qsym(:, nb), vec(:, nbeff, is), q0, is)
          save1 = enrear !!! TO deactivate cumulation of enrear, enerpw
          save2 = enerpw !!! (ELSE wrong total energy)
          CALL calc_sicsp(rhosp, usicsp, qsym(1, nb), nb)
          ener_2st(is) = ener_2st(is) + (encoulsp + enerpw)*occup(nb)
          enrear = save1
          enerpw = save2

          DO ind = 1, nxyz
            uqsym(ind, nb) = usicsp(ind + ishift)*qsym(ind, nb)
          END DO
        END IF
      END DO

! variational RESULT of the UT constraint

      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)
          DO na = 1, nstate
            IF (ispin(nrel2abs(na)) == is) THEN
              naa = na - (is - 1)*ndims(1)
              utcond(naa, nbeff) = -wfovlp(qsym(:, na), uqsym(:, nb))
            END IF
          END DO
        END IF
      END DO

      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)
          DO na = 1, nstate
            IF (ispin(nrel2abs(na)) == is) THEN
              naa = na - (is - 1)*ndims(1)
              dab(naa, nbeff) = utcond(naa, nbeff) - CONJG(utcond(nbeff, naa))
            END IF
          END DO
        END IF
      END DO

      DEALLOCATE (usicsp, rhosp)
      DEALLOCATE (qsym, utcond)
      DEALLOCATE (uqsym, symcond)

      RETURN
    END SUBROUTINE dalphabeta_rc

!-----------DAlphaBetar-------------------

    SUBROUTINE dalphabeta_r(is, dab, vec, q0)

! Compose matrix 'dab=<psi_a|h_a-h_b|psi_b>. REAL version.
!
! Input:
! vec = actual states IN configuration space
! q0 = set of s.p. wavefunctions IN coordinate space
! is = spin
! Output:
! dab = matrix of symmetry condition

      USE params
      USE kinetic
      USE util, ONLY: wfovlp
      USE twost_util, ONLY: superpose_state
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: is
      REAL(DP), INTENT(OUT) :: dab(kdim, kdim)
      REAL(DP), INTENT(IN) :: vec(kstate, kstate, 2)
      REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)

      LOGICAL, PARAMETER :: ttest = .false.
      INTEGER :: ind, ishift, na, naa, nb, nbeff
      REAL(DP) :: save1, save2, acc2

      REAL(DP), ALLOCATABLE :: usicsp(:), rhosp(:)
      REAL(DP), ALLOCATABLE :: uqsym(:, :), symcond(:, :)
      REAL(DP), ALLOCATABLE :: utcond(:, :), qsym(:, :)

!-------------------------------------------------------

      ALLOCATE (usicsp(2*kdfull2), rhosp(2*kdfull2))
      ALLOCATE (uqsym(kdfull2, kstate), symcond(kstate, kstate))
      ALLOCATE (utcond(kdim, kdim), qsym(kdfull2, kstate))

      ishift = (is - 1)*nxyz
      dab = 0D0
      ener_2st(is) = 0D0
      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)

          CALL superpose_state(qsym(:, nb), vec(:, nbeff, is), q0, is)
          save1 = enrear !!! TO deactivate cumulation of enrear, enerpw
          save2 = enerpw !!! (ELSE wrong total energy)
          CALL calc_sicspr(rhosp, usicsp, qsym(1, nb), nb)
! CALL printfield(6,rhosp,'rhosp')
! CALL printfield(6,usicsp,'usicsp')
IF(ttest) THEN
WRITE(*,*) 'vec:',vec(:, nbeff, is)
WRITE(*,*) 'is,nb,encoulsp,enerpw=',is,nb,encoulsp,enerpw,occup(nb)
WRITE(*,*) 'is,nb,norm-wf=',is,nb,SUM(qsym(:,nb)**2)*dvol,&
 SUM(rhosp),SUM(usicsp)
END IF
        ener_2st(is) = ener_2st(is) + (encoulsp + enerpw)*occup(nb)
        enrear = save1
        enerpw = save2

          DO ind = 1, nxyz
            uqsym(ind, nb) = usicsp(ind + ishift)*qsym(ind, nb)
          END DO
        END IF
      END DO !nb

! variational RESULT of the UT constraint

      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)
          DO na = 1, nstate
            IF (ispin(nrel2abs(na)) == is) THEN
              naa = na - (is - 1)*ndims(1)
              utcond(naa, nbeff) = -wfovlp(qsym(:, na), uqsym(:, nb))
            END IF !ispin
          END DO !na
        END IF
      END DO !nb
! WRITE(6,*) ' utcond: is=',is
! DO i=1,ndims(is)
! WRITE(6,'(20(1pg13.5))') utcond(1:ndims(is),i)
! END DO
! CALL rMatPrint ('utCond', utcond, kstate, ndims(is))
      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)
          DO na = 1, nstate
            IF (ispin(nrel2abs(na)) == is) THEN
              naa = na - (is - 1)*ndims(1)
              dab(naa, nbeff) = utcond(naa, nbeff) - utcond(nbeff, naa)
            END IF
          END DO
        END IF
      END DO

      IF (ttest) THEN
        WRITE (*, *) 'dab:'
        WRITE (*, *) dab(1:ndims(is), 1:ndims(is))
      END IF

      DEALLOCATE (usicsp, rhosp)
      DEALLOCATE (uqsym, symcond)
      DEALLOCATE (utcond, qsym)

      RETURN
    END SUBROUTINE dalphabeta_r

!-----------DAlphaBeta-------------------

    SUBROUTINE dalphabeta_c(is, dab, vec, q0)

! Compose matrix 'dab=<psi_a|h_a-h_b|psi_b>. COMPLEX version.
!
! Input:
! vec = actual states IN configuration space
! q0 = set of s.p. wavefunctions IN coordinate space
! is = spin
! Output:
! dab = matrix of symmetry condition

      USE params
      USE kinetic
      USE util, ONLY: wfovlp
      USE twost_util, ONLY: superpose_state
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: is
      COMPLEX(DP), INTENT(OUT) :: dab(kdim, kdim) ! DAB
      COMPLEX(DP), INTENT(IN) :: vec(kstate, kstate, 2)
      COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)

      LOGICAL, PARAMETER :: ttest = .false.
      INTEGER :: ind, ishift, na, naa, nb, nbeff
      REAL(DP) :: save1, save2, acc2

      REAL(DP), ALLOCATABLE :: usicsp(:), rhosp(:)
      COMPLEX(DP), ALLOCATABLE :: qsym(:, :), utcond(:, :)
      COMPLEX(DP), ALLOCATABLE :: uqsym(:, :), symcond(:, :)

!-------------------------------------------------------

      ALLOCATE (usicsp(2*kdfull2), rhosp(2*kdfull2))
      ALLOCATE (qsym(kdfull2, kstate), utcond(kdim, kdim))
      ALLOCATE (uqsym(kdfull2, kstate), symcond(kstate, kstate))

      ishift = (is - 1)*nxyz

      ener_2st(is) = 0D0
      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)

          CALL superpose_state(qsym(:, nb), vec(:, nbeff, is), q0, is)
          save1 = enrear !!! TO deactivate cumulation of enrear, enerpw
          save2 = enerpw !!! (ELSE wrong total energy)
! WRITE(*,*) 'DALPHAETA: nb,wfnorm=',nb,wfnorm(qsym(1,nb))
          CALL calc_sicsp(rhosp, usicsp, qsym(1, nb), nb)
! WRITE(*,*) ' nb,energs=',nb,encoulsp,enerpw,occup(nb)
          ener_2st(is) = ener_2st(is) + (encoulsp + enerpw)*occup(nb)
          enrear = save1
          enerpw = save2

          DO ind = 1, nxyz
            uqsym(ind, nb) = usicsp(ind + ishift)*qsym(ind, nb)
          END DO! nxyz
        END IF
      END DO !nb

! variational RESULT of the UT constraint

      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)
          DO na = 1, nstate
            IF (ispin(nrel2abs(na)) == is) THEN
              naa = na - (is - 1)*ndims(1)
              utcond(naa, nbeff) = -wfovlp(qsym(:, na), uqsym(:, nb))
            END IF !ispin
          END DO !na
        END IF
      END DO !nb
!WRITE(*,*) ' utcond:',utcond(1:nstate,a:nstate)
! CALL MatPrint ('utCond', utcond, kstate, ndims(is))
      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == is) THEN
          nbeff = nb - (is - 1)*ndims(1)
          DO na = 1, nstate
            IF (ispin(nrel2abs(na)) == is) THEN
              naa = na - (is - 1)*ndims(1)
              dab(naa, nbeff) = utcond(naa, nbeff) - CONJG(utcond(nbeff, naa))
            END IF
          END DO
        END IF
      END DO
! CALL MatPrint ('DAB', DAB, kstate, ndims(is))

      DEALLOCATE (usicsp, rhosp)
      DEALLOCATE (qsym, utcond)
      DEALLOCATE (uqsym, symcond)

      RETURN
    END SUBROUTINE dalphabeta_c

!----------------------------------------------------------------------
! REAL version
!----------------------------------------------------------------------
    SUBROUTINE spmomsmatrix_r(wfr, PRINT)

! Matrix of spatial moments between single-particle states
! from REAL wf's:
! wfr = set of single particle wavefunctions
! PRINT = switch TO PRINT variance
! The RESULT is stored within the MODULE variables.

      USE params
      USE kinetic
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: PRINT
      REAL(DP), INTENT(IN) :: wfr(kdfull2, kstate)

      INTEGER, PARAMETER :: iunit = 0 ! set zero TO disable testprint

      LOGICAL :: tfirst = .TRUE.
      INTEGER :: ind, is, ix, iy, iz, ma, mb, na, nac, nb, noff
      REAL(DP) :: var, x1, y1, z1, x2, y2, z2
      REAL(DP) :: s, wfmom, xmom, ymom, zmom, xxmom, yymom, zzmom
      REAL(DP), ALLOCATABLE :: rrmatr(:, :, :) ! matrix of r**2
      REAL(DP), ALLOCATABLE :: xxmatr(:, :, :) ! matrix of x**2
      REAL(DP), ALLOCATABLE :: yymatr(:, :, :) ! matrix of y**2
      REAL(DP), ALLOCATABLE :: zzmatr(:, :, :) ! matrix of z**2
      REAL(DP), ALLOCATABLE :: xmatr(:, :, :) ! matrix of x
      REAL(DP), ALLOCATABLE :: ymatr(:, :, :) ! matrix of y
      REAL(DP), ALLOCATABLE :: zmatr(:, :, :) ! matrix of z
!----------------------------------------------------------------------

      kdim = kstate
      ALLOCATE (rrmatr(kdim, kdim, 2)) ! matrix of r**2
      ALLOCATE (xxmatr(kdim, kdim, 2)) ! matrix of x**2
      ALLOCATE (yymatr(kdim, kdim, 2)) ! matrix of y**2
      ALLOCATE (zzmatr(kdim, kdim, 2)) ! matrix of z**2
      ALLOCATE (xmatr(kdim, kdim, 2)) ! matrix of x
      ALLOCATE (ymatr(kdim, kdim, 2)) ! matrix of y
      ALLOCATE (zmatr(kdim, kdim, 2)) ! matrix of z
!
      IF (iunit > 0) OPEN (iunit, POSITION='append', FILE='pstatmomsmatrix.'//outnam)

      noff = ndims(1)
      WRITE (*, '(a,2i5)') 'NDIMS=', ndims

! compute matrix elements and store

      IF (iunit > 0) THEN
        WRITE (iunit, '(a)') 'matrix of s.p. moments:', &
          ' na nb x y z xx yy zz'
        tfirst = .false.
      END IF
      DO na = 1, nstate

        DO nb = 1, na
          IF (ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
            wfmom = 0D0
            xmom = 0D0
            ymom = 0D0
            zmom = 0D0
            xxmom = 0D0
            yymom = 0D0
            zzmom = 0D0
            ind = 0
            DO iz = minz, maxz
              z1 = (iz - nzsh)*dz
              z2 = z1*z1
              DO iy = miny, maxy
                y1 = (iy - nysh)*dy
                y2 = y1*y1
                DO ix = minx, maxx
                  ind = ind + 1
                  IF ((ix /= nx2) .AND. (iy /= ny2) .AND. (iz /= nz2)) THEN
                    x1 = (ix - nxsh)*dx
                    x2 = x1*x1
                    s = wfr(ind, na)*wfr(ind, nb)
                    wfmom = wfmom + s ! monopole
                    xmom = xmom + s*x1 ! dipole
                    ymom = ymom + s*y1
                    zmom = zmom + s*z1
                    xxmom = xxmom + s*x2 ! quadrupole
                    yymom = yymom + s*y2
                    zzmom = zzmom + s*z2
                  END IF
                END DO
              END DO
            END DO

            wfmom = dvol*wfmom
            xmom = dvol*xmom
            ymom = dvol*ymom
            zmom = dvol*zmom
            xxmom = dvol*xxmom
            yymom = dvol*yymom
            zzmom = dvol*zzmom

            IF (iunit > 0) WRITE (iunit, '(3i4,6(1pg13.5))') &
              na, nb, ispin(nrel2abs(nb)), xmom, ymom, zmom, xxmom, yymom, zzmom

            is = ispin(nrel2abs(nb))

            IF (is == 1) THEN
              ma = na
              mb = nb
            ELSE
              ma = na - noff
              mb = nb - noff
            END IF
            IF (iunit > 0) WRITE (iunit, *) 'ma,mb=', ma, mb, na, nb
            xmatr(ma, mb, is) = xmom
            ymatr(ma, mb, is) = ymom
            zmatr(ma, mb, is) = zmom
            xxmatr(ma, mb, is) = xxmom
            yymatr(ma, mb, is) = yymom
            zzmatr(ma, mb, is) = zzmom
            xmatr(mb, ma, is) = xmom
            ymatr(mb, ma, is) = ymom
            zmatr(mb, ma, is) = zmom
            xxmatr(mb, ma, is) = xxmom

            yymatr(mb, ma, is) = yymom
            zzmatr(mb, ma, is) = zzmom
          END IF
        END DO
      END DO

! initialize eigen-vectors

      DO is = 1, 2
        DO na = 1, ndims(is)
          DO nb = 1, ndims(is)
            rrmatr(na, nb, is) = xxmatr(na, nb, is) + yymatr(na, nb, is) &
                                 + zzmatr(na, nb, is)
          END DO
        END DO
      END DO

      IF (PRINT == 1) THEN
        DO is = 1, 2
          var = 0D0
          WRITE (6, '(a)') ' na x y z xx yy zz'
          DO na = 1, ndims(is)
            var = var + rrmatr(na, na, is) - xmatr(na, na, is)**2 - &
                  ymatr(na, na, is)**2 - zmatr(na, na, is)**2
            nac = na + (is - 1)*noff
            WRITE (6, '(i5,6(1pg13.5))') &
              na, REAL(xmatr(na, na, is), DP), REAL(ymatr(na, na, is), DP), &
              REAL(zmatr(na, na, is), DP), REAL(xxmatr(na, na, is), DP), &
              REAL(yymatr(na, na, is), DP), REAL(zzmatr(na, na, is), DP)
          END DO
          WRITE (6, *) 'For spin =', is, ' Total spacial variance ', var
        END DO
      END IF

      IF (iunit > 0) CLOSE (iunit)
      DEALLOCATE (rrmatr) ! matrix of r**2
      DEALLOCATE (xxmatr) ! matrix of x**2
      DEALLOCATE (yymatr) ! matrix of y**2
      DEALLOCATE (zzmatr) ! matrix of z**2
      DEALLOCATE (xmatr) ! matrix of x
      DEALLOCATE (ymatr) ! matrix of y
      DEALLOCATE (zmatr) ! matrix of z

      RETURN
    END SUBROUTINE spmomsmatrix_r

!-----------------------------------------------------------------------
! COMPLEX version
!-----------------------------------------------------------------------

    SUBROUTINE spmomsmatrix_c(wfr, PRINT)

! Matrix of spatial moments between single-particle states
! from COMPLEX wf's:
! wfr = set of single particle wavefunctions
! PRINT = switch TO PRINT variance
! The RESULT is stored within the MODULE variables.

      USE params
      USE kinetic
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: PRINT
      COMPLEX(DP), INTENT(IN) :: wfr(kdfull2, kstate)

      INTEGER, PARAMETER :: iunit = 0 ! set zero TO disable testprint

      LOGICAL :: tfirst = .TRUE.
      INTEGER :: ind, is, ix, iy, iz, ma, mb, na, nac, nb, noff
      REAL(DP) :: var, x1, y1, z1, x2, y2, z2
      COMPLEX(DP) :: s, wfmom, xmom, ymom, zmom, xxmom, yymom, zzmom
      COMPLEX(DP), ALLOCATABLE :: rrmatr(:, :, :) ! matrix of r**2
      COMPLEX(DP), ALLOCATABLE :: xxmatr(:, :, :) ! matrix of x**2
      COMPLEX(DP), ALLOCATABLE :: yymatr(:, :, :) ! matrix of y**2
      COMPLEX(DP), ALLOCATABLE :: zzmatr(:, :, :) ! matrix of z**2
      COMPLEX(DP), ALLOCATABLE :: xmatr(:, :, :) ! matrix of x
      COMPLEX(DP), ALLOCATABLE :: ymatr(:, :, :) ! matrix of y
      COMPLEX(DP), ALLOCATABLE :: zmatr(:, :, :) ! matrix of z

!----------------------------------------------------------------------

      kdim = kstate
      ALLOCATE (rrmatr(kdim, kdim, 2)) ! matrix of r**2
      ALLOCATE (xxmatr(kdim, kdim, 2)) ! matrix of x**2
      ALLOCATE (yymatr(kdim, kdim, 2)) ! matrix of y**2
      ALLOCATE (zzmatr(kdim, kdim, 2)) ! matrix of z**2
      ALLOCATE (xmatr(kdim, kdim, 2)) ! matrix of x
      ALLOCATE (ymatr(kdim, kdim, 2)) ! matrix of y
      ALLOCATE (zmatr(kdim, kdim, 2)) ! matrix of z
!
      IF (iunit > 0) OPEN (iunit, POSITION='append', FILE='pstatmomsmatrix.'//outnam)

      noff = ndims(1)

      WRITE (*, '(a,2i5)') 'NDIMS=', ndims

! compute matrix elements and store

      IF (iunit > 0) THEN
        WRITE (iunit, '(a)') 'matrix of s.p. moments:', &
          ' na nb x y z xx yy zz'
        tfirst = .false.
      END IF
      DO na = 1, nstate

        DO nb = 1, na
          IF (ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
            wfmom = CMPLX(0D0, 0D0, DP)
            xmom = CMPLX(0D0, 0D0, DP)
            ymom = CMPLX(0D0, 0D0, DP)
            zmom = CMPLX(0D0, 0D0, DP)
            xxmom = CMPLX(0D0, 0D0, DP)
            yymom = CMPLX(0D0, 0D0, DP)
            zzmom = CMPLX(0D0, 0D0, DP)

            ind = 0
            DO iz = minz, maxz
              z1 = (iz - nzsh)*dz
              z2 = z1*z1
              DO iy = miny, maxy
                y1 = (iy - nysh)*dy
                y2 = y1*y1
                DO ix = minx, maxx
                  ind = ind + 1
                  IF ((ix /= nx2) .AND. (iy /= ny2) .AND. (iz /= nz2)) THEN
                    x1 = (ix - nxsh)*dx
                    x2 = x1*x1
                    s = CONJG(wfr(ind, na))*wfr(ind, nb)
                    wfmom = wfmom + s ! monopole
                    xmom = xmom + s*x1 ! dipole
                    ymom = ymom + s*y1
                    zmom = zmom + s*z1
                    xxmom = xxmom + s*x2 ! quadrupole
                    yymom = yymom + s*y2
                    zzmom = zzmom + s*z2
                  END IF
                END DO
              END DO
            END DO

            wfmom = dvol*wfmom
            xmom = dvol*xmom
            ymom = dvol*ymom
            zmom = dvol*zmom
            xxmom = dvol*xxmom
            yymom = dvol*yymom
            zzmom = dvol*zzmom

            IF (iunit > 0) WRITE (iunit, '(3i4,6(1pg13.5))') &
              na, nb, ispin(nrel2abs(nb)), xmom, ymom, zmom, xxmom, yymom, zzmom

            is = ispin(nrel2abs(nb))

            IF (is == 1) THEN
              ma = na
              mb = nb
            ELSE
              ma = na - noff
              mb = nb - noff
            END IF
            IF (iunit > 0) WRITE (iunit, *) 'ma,mb=', ma, mb, na, nb
            xmatr(ma, mb, is) = xmom
            ymatr(ma, mb, is) = ymom
            zmatr(ma, mb, is) = zmom
            xxmatr(ma, mb, is) = xxmom
            yymatr(ma, mb, is) = yymom
            zzmatr(ma, mb, is) = zzmom
            xmatr(mb, ma, is) = CONJG(xmom)
            ymatr(mb, ma, is) = CONJG(ymom)
            zmatr(mb, ma, is) = CONJG(zmom)
            xxmatr(mb, ma, is) = CONJG(xxmom)
            yymatr(mb, ma, is) = CONJG(yymom)
            zzmatr(mb, ma, is) = CONJG(zzmom)
          END IF
        END DO
      END DO

! initialize eigen-vectors

      DO is = 1, 2
        DO na = 1, ndims(is)
          DO nb = 1, ndims(is)
            rrmatr(na, nb, is) = xxmatr(na, nb, is) + yymatr(na, nb, is) &
                                 + zzmatr(na, nb, is)
          END DO
        END DO
      END DO

      IF (PRINT == 1) THEN
        DO is = 1, 2
          var = CMPLX(0D0, 0D0, DP)
          WRITE (6, '(a)') ' na x y z xx yy zz'
          DO na = 1, ndims(is)
            var = var + rrmatr(na, na, is) - xmatr(na, na, is)**2 - &
                  ymatr(na, na, is)**2 - zmatr(na, na, is)**2
            nac = na + (is - 1)*noff
            WRITE (6, '(i5,6(1pg13.5))') &
              na, REAL(xmatr(na, na, is), DP), REAL(ymatr(na, na, is), DP), &
              REAL(zmatr(na, na, is), DP), REAL(xxmatr(na, na, is), DP), &
              REAL(yymatr(na, na, is), DP), REAL(zzmatr(na, na, is), DP)
          END DO
          WRITE (6, *) 'For spin =', is, ' Total spacial variance ', var
        END DO
      END IF

      IF (iunit > 0) CLOSE (iunit)
      DEALLOCATE (rrmatr) ! matrix of r**2
      DEALLOCATE (xxmatr) ! matrix of x**2
      DEALLOCATE (yymatr) ! matrix of y**2
      DEALLOCATE (zzmatr) ! matrix of z**2
      DEALLOCATE (xmatr) ! matrix of x
      DEALLOCATE (ymatr) ! matrix of y
      DEALLOCATE (zmatr) ! matrix of z

      RETURN
    END SUBROUTINE spmomsmatrix_c

#endif

!#endif

#ifdef REALSWITCH
  END MODULE twostr
#endif

#ifdef COMPLEXSWITCH
END MODULE twost
#endif

#ifdef REALSWITCH
SUBROUTINE calc_fullsicr(q0, qsic)
#else
  SUBROUTINE calc_fullsic(q0, qsic)
#endif

! ******************************

! Full SIC:
! input is set of wavefunctions on 'q0'
! output are SIC s.p. wavefunctions on 'qsic'

    USE params
    USE kinetic
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(OUT) :: qsic(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
    COMPLEX(DP), INTENT(OUT) :: qsic(kdfull2, kstate)
#endif

    LOGICAL, PARAMETER :: ttest = .TRUE.

    INTEGER :: ind, idx, ishift, nb
    REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd

! workspaces
    REAL(DP), ALLOCATABLE :: usicsp(:), rhosp(:)

    IF (numspin .NE. 2) STOP "CALC_FULLSIC requires full spin"

    enrearsave = enrear
    enerpwsave = enerpw
    enrear1 = 0D0
    enrear2 = 0D0
    enpw1 = 0D0
    enpw2 = 0D0
    encadd = 0D0

! compute ACTION of SIC potential
    ALLOCATE (usicsp(2*kdfull2), rhosp(2*kdfull2))

    DO nb = 1, nstate
      IF (occup(nb) > small) THEN
        ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 IN upper BLOCK
#ifdef REALSWITCH
        IF (ifsicp == 9) THEN
          CALL calc_sicsp(rhosp, usicsp, q0(1, nb), nb)
        ELSE IF (ifsicp == 8) THEN
          CALL calc_sicspr(rhosp, usicsp, q0(1, nb), nb)
        ELSE
          STOP "invalid IFSICP IN CALC_FULLSIC"
        END IF
#else
        CALL calc_sicsp(rhosp, usicsp, q0(1, nb), nb)
#endif

        IF (ispin(nrel2abs(nb)) == 1) THEN
          enrear1 = enrear1 + enrear*occup(nb)
          IF (directenergy) THEN
            enpw1 = enpw1 + enerpw*occup(nb)
          END IF
        ELSE
          enrear2 = enrear2 + enrear*occup(nb)
          IF (directenergy) THEN
            enpw2 = enpw2 + enerpw*occup(nb)
          END IF
        END IF
        IF (directenergy) THEN
          encadd = encadd + encoulsp*occup(nb)
        END IF

        IF (ispin(nrel2abs(nb)) == 1) THEN
          DO ind = 1, nxyz
            qsic(ind, nb) = usicsp(ind)*q0(ind, nb)
          END DO
        ELSE
          DO ind = 1, nxyz
            idx = ind + nxyz
            qsic(ind, nb) = usicsp(idx)*q0(ind, nb)
          END DO
        END IF
      END IF
    END DO

    enrear = enrearsave - enrear1 - enrear2
    IF (directenergy) THEN
      enerpw = enerpwsave - enpw1 - enpw2 - encadd
    END IF
    DEALLOCATE (usicsp, rhosp)

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE calc_fullsicr
#else
END SUBROUTINE calc_fullsic
#endif

