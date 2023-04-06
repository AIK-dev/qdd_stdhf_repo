
#if(symmcond&&nompi)
MODULE symcond
  USE params
  USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

  SAVE

  REAL(DP) :: usicall(kdfull2, kstate)
!COMMON /sicsav/usicall(kdfull2,kstate)

CONTAINS
! ******************************

  SUBROUTINE symmcond_step(q0, qaux, symmstep, symcond_add)

! ******************************

! iteration of symmetry condition for the CASE
! of REAL wavefunctions.

! list parameters:
! q0 = set of wavefunctions TO be propagated
! qaux = workspace of same complexity as 'q0'
! symmstep = step SIZE
! symcond_add = accumulated symmetry condition (at ENTRY)

! the routine requires that the SIC potentials have been
! computed before IN 'calc_fullsic'.

    USE params

    REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(OUT) :: qaux(kdfull2, kstate)
    REAL(DP), INTENT(IN) :: symmstep
    REAL(DP), INTENT(IN OUT) :: symcond_add

    REAL(DP) :: qaux2(kdfull2)
    REAL(DP) :: xlambda(kstate, kstate)

! the SIC potentials come IN from 'calc_fullsic'
! through that COMMON BLOCK

!*********************************************************

! outer loop over state TO be propagated,
! the stepped step is SAVE temporarily on 'qaux',
! the symmetry condition is accumulated simultaneously.

    symcond_add = 0D0
    DO nc = 1, nstate

! initialize accumulator

      DO i = 1, nxyz
        qaux(i, nc) = q0(i, nc)
      END DO
      DO nb = 1, nstate
        xlambda(nc, nb) = 0D0
      END DO

! inner loop over partner states IN symmetry condition,
! compute symmetry condition when new,
! and step.

      DO na = 1, nstate
        IF (ispin(na) == ispin(nc)) THEN
          DO i = 1, nxyz
            qaux2(i) = (usicall(i, nc) - usicall(i, na))*q0(i, na)
          END DO
          symac = rwfovlp(qaux2, q0(1, nc))
          symcond_add = symcond_add + symac**2

! project onto occupied states

          DO nb = 1, nstate
            IF (ispin(nb) == ispin(nc)) THEN
              ovabc = rwfovlp(qaux2, q0(1, nb))
              factor = symac*ovabc
              xlambda(nc, nb) = factor + xlambda(nc, nb)
              factor = -factor*symmstep
              DO i = 1, nxyz
                qaux(i, nc) = q0(i, nb)*factor + qaux(i, nc)
              END DO
            END IF
          END DO
        END IF
      END DO
    END DO

! symmetrize Lagrangian multipliers and subtract

    DO nc = 1, nstate
      DO nb = 1, nc - 1
        xlambda(nc, nb) = 0.5D0*(xlambda(nc, nb) + xlambda(nb, nc))
        xlambda(nb, nc) = xlambda(nc, nb)
      END DO
    END DO
    DO nc = 1, nstate
      DO nb = 1, nstate
        IF (ispin(nb) == ispin(nc)) THEN
          factor = -symmstep*xlambda(nc, nb)
          DO i = 1, nxyz
            qaux(i, nc) = qaux(i, nc) - factor*q0(i, nb)
          END DO
        END IF
      END DO
    END DO

! store RESULT and so update 'q0'

    DO nc = 1, nstate
      DO i = 1, nxyz
        q0(i, nc) = qaux(i, nc)
      END DO
    END DO

! ortho-normalize

    CALL schmidt(q0)

    RETURN
  END SUBROUTINE symmcond_step

! ******************************

  SUBROUTINE symmcond_emin(q0, qaux, symmstep, symcond_add)

! ******************************

! Iteration of symmetry condition for the CASE
! of REAL wavefunctions.
! The energy gradient is used.

! list parameters:
! q0 = set of wavefunctions TO be propagated
! qaux = workspace of same complexity as 'q0'
! symmstep = step SIZE
! symcond_add = accumulated symmetry condition (at ENTRY)

! the routine requires that the SIC potentials have been
! computed before IN 'calc_fullsic'.

    USE params
    USE kinetic
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(OUT) :: qaux(kdfull2, kstate)
    REAL(DP), INTENT(IN) :: symmstep
    REAL(DP), INTENT(IN OUT) :: symcond_add

    REAL(DP) :: qaux2(kdfull2)
    REAL(DP) :: xlambda(kstate, kstate)

! the SIC potentials come IN from 'calc_fullsic'
! through that COMMON BLOCK
!COMMON /sicsav/usicall(kdfull2,kstate)

! the matrices for vector and curvature of stepping

    REAL(DP) :: symcond(kstate, kstate)

    LOGICAL, PARAMETER :: ttest = .true.

!*********************************************************

! the symmetry condition = constrained energy gradient

    sum = 0D0
    DO na = 1, nstate
      DO nb = 1, nstate
        IF (ispin(na) == ispin(nb)) THEN
          DO i = 1, nxyz
            qaux2(i) = (usicall(i, na) - usicall(i, nb))*q0(i, nb)
          END DO
          symcond(na, nb) = rwfovlp(q0(1, na), qaux2)
          sum = symcond(na, nb)**2 + sum
        ELSE
          symcond(na, nb) = 0.0D0
        END IF
      END DO
    END DO
    symcond_add = sum
    IF (ttest) THEN
      WRITE (6, '(a,1pg13.5)') ' SYMMCON-DIRECT:', SQRT(sum/nstate)
    END IF

! the step

    DO na = 1, nstate
      DO i = 1, nxyz
        qaux(i, na) = 0D0
      END DO
      DO nb = 1, nstate
        IF (ispin(ni) == ispin(nc)) THEN
          factor = symmstep*symcond(na, nb)
          DO i = 1, nxyz
            qaux(i, na) = qaux(i, na) + factor*q0(i, nb)
          END DO
        END IF
      END DO
      IF (ttest) THEN
        WRITE (7, '(a,i5,1pg13.5)') ' SYMMCON-EMIN: na, norm of changes=', ny, &
          rwfovlp(qaux(1, na), qaux(1, na))
      END IF
    END DO

! store RESULT and so update 'q0'

    DO na = 1, nstate
      DO i = 1, nxyz
        q0(i, na) = qaux(i, na) + q0(i, na)
      END DO
      IF (ttest) THEN
        WRITE (7, '(a,i5,1pg13.5)') &
          ' SYMMCON-EMIN: na,norm=', na, rwfovlp(q0(1, na), q0(1, na))
      END IF
    END DO

! ortho-normalize

    CALL schmidt(q0)

    RETURN
  END SUBROUTINE symmcond_emin

! ******************************

  SUBROUTINE symmcond_direct(q0, qaux, symmstep, symcond_add)

! ******************************

! Iteration of symmetry condition for the CASE
! of REAL wavefunctions.
! The DIRECT Newton-Raphson step is used.

! list parameters:
! q0 = set of wavefunctions TO be propagated
! qaux = workspace of same complexity as 'q0'
! symmstep = step SIZE
! symcond_add = accumulated symmetry condition (at ENTRY)

! the routine requires that the SIC potentials have been
! computed before IN 'calc_fullsic'.

    USE params
    USE kinetic
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(OUT) :: qaux(kdfull2, kstate)
    REAL(DP), INTENT(IN) :: symmstep
    REAL(DP), INTENT(IN OUT) :: symcond_add

    REAL(DP) :: qaux2(kdfull2)
    REAL(DP) :: xlambda(kstate, kstate)

! the SIC potentials come IN from 'calc_fullsic'
! through that COMMON BLOCK
!COMMON /sicsav/usicall(kdfull2,kstate)

! the matrices for vector and curvature of stepping

    INTEGER, PARAMETER :: ksymcond = kstate*kstate
    REAL(DP) :: symcondvec(ksymcond)
    REAL(DP) :: stepvec(ksymcond)
    REAL(DP) :: symcondmat(ksymcond, ksymcond)
    INTEGER :: iptsymcond(2, ksymcond)

    LOGICAL, PARAMETER :: ttest = .false.

!*********************************************************

! the condition vector

    sum = 0D0
    ipt = 0
    DO na = 1, nstate
      DO nb = 1, nstate
        ipt = 1 + ipt
        iptsymcond(1, ipt) = na
        iptsymcond(2, ipt) = nb
        IF (na < nb .AND. ispin(na) == ispin(nb)) THEN
          DO i = 1, nxyz
            qaux2(i) = (usicall(i, na) - usicall(i, nb))*q0(i, nb)
          END DO
          symcondvec(ipt) = rwfovlp(q0(1, na), qaux2)
          sum = symcondvec(ipt)**2 + sum
        ELSE
          symcondvec(ipt) = 0.0D0
        END IF
      END DO
    END DO
    iptmax = ipt
    sum = 2.0*sum
    symcond_add = sum
    IF (ttest) THEN
      WRITE (6, '(a,1pg13.5)') ' SYMMCON-DIRECT:', SQRT(sum/nstate)
      WRITE (6, '(a)') ' symcondvec:'
      DO ni = 1, nstate
        WRITE (6, '(20(1pg13.5))') (symcondvec(ipt), ipt=(ni - 1)*nstate + 1, ni*nstate)
      END DO
    END IF

! the condition matrix

    DO ivec1 = 1, iptmax
      DO ivec2 = 1, iptmax
        ni = iptsymcond(1, ivec1)
        nc = iptsymcond(2, ivec1)
        na = iptsymcond(1, ivec2)
        nb = iptsymcond(2, ivec2)
        IF (ispin(na) == ispin(nb) .AND. ispin(na) == ispin(nc) &
            .AND. ispin(na) == ispin(ni)) THEN
          IF (na /= nc) THEN
            symcondmat(ivec1, ivec2) = 0.0D0
          ELSE
            IF (na < nb) THEN
              DO i = 1, nxyz
                qaux2(i) = (usicall(i, na) - usicall(i, nb))*q0(i, nb)
              END DO
              symcondmat(ivec1, ivec2) = rwfovlp(q0(1, ni), qaux2)
            ELSE
              IF (ni == nb) THEN
                symcondmat(ivec1, ivec2) = 1.0D0
              ELSE
                symcondmat(ivec1, ivec2) = 0.0D0
              END IF
            END IF
          END IF
        ELSE
          symcondmat(ivec1, ivec2) = 0.0D0
        END IF
      END DO
    END DO
    IF (ttest) THEN
      WRITE (6, '(a)') ' symcondmat:'
      DO ivec2 = 1, iptmax
        WRITE (6, '(20(1pg13.5))') (symcondmat(ivec1, ivec2), ivec1=1, iptmax)
      END DO
    END IF

! invert matrix

    CALL matin2(symcondmat, iptmax, ksymcond, nerror, determ)

! compose step-vector

    DO ivec1 = 1, iptmax
      sum = 0D0
      DO ivec2 = 1, iptmax
        sum = sum + symcondvec(ivec2)*symcondmat(ivec2, ivec1)
      END DO
      stepvec(ivec1) = sum
    END DO
    IF (ttest) THEN
      WRITE (6, '(a,i10,1pg12.5)') ' nerror,determ=', nerror, determ
      WRITE (6, '(a)') ' inverse symcondmat:'
      DO ivec2 = 1, iptmax
        WRITE (6, '(20(1pg13.5))') (symcondmat(ivec1, ivec2), ivec1=1, iptmax)
      END DO
      WRITE (6, '(a)') ' stepvec:'
      DO ni = 1, nstate
        WRITE (6, '(20(1pg13.5))') (stepvec(ipt), ipt=(ni - 1)*nstate + 1, ni*nstate)
      END DO
    END IF

! compose step

    DO nc = 1, nstate
      DO i = 1, nxyz
        qaux(i, nc) = q0(i, nc)
      END DO
    END DO
    ipt = 0
    DO ni = 1, nstate
      DO nc = 1, nstate
        ipt = 1 + ipt
        IF (ispin(ni) == ispin(nc)) THEN
          factor = -symmstep*stepvec(ipt)
          DO i = 1, nxyz
            qaux(i, nc) = qaux(i, nc) + factor*q0(i, ni)
          END DO
        END IF

      END DO
    END DO

! store RESULT and so update 'q0'

    DO nc = 1, nstate
      DO i = 1, nxyz
        q0(i, nc) = qaux(i, nc)
      END DO
      IF (ttest) THEN
        WRITE (7, '(a,1pg13.5)') &
          ' SYMMCON-DIRECT: overlap=', rwfovlp(q0(1, nc), q0(1, nc))
      END IF
    END DO

! ortho-normalize

    CALL schmidt(q0)

    RETURN
  END SUBROUTINE symmcond_direct

!-----matin2----------------------------------------------------------

  SUBROUTINE matin2(a, nmax, ndim, nerror, determ)

! include 'params.inc'
    USE params, ONLY: DP

    REAL(DP), INTENT(OUT) :: a(*)
    INTEGER, INTENT(IN) :: nmax
    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(OUT) :: nerror
    REAL(DP), INTENT(OUT) :: determ
    INTEGER, PARAMETER :: kindex = 900

! inverts the matrix 'c'.

! the parameters on the i/o list are:

! a = at input: matrix TO be inverted
! at output: inverse matrix
! nmax = NUMBER of actual rows and columns
! ndim = DIMENSION of matrix IN the calling routine
! nerror = error flag
! determ = determinant of input matrix 'a'

! the routine is taken from 'cernlib' and modified for this
! package.

    INTEGER :: INDEX(kindex)

!-----------------------------------------------------------------------------

    IF (nmax > kindex) STOP ' too large DIMENSION IN matin2'

    deter = 1D0
    n = nmax
! iemat=n+n2
    iemat = n + 0
    idim = nmax
    nmin1 = n - 1

! umsetzen des nmax*nmax-anteils der matrix a IN gepackte
! ein-dimensionale FORM.

    DO i = 2, idim
      DO j = 1, idim
        ipack = (i - 1)*nmax + j
        imatr = (i - 1)*ndim + j
        a(ipack) = a(imatr)
      END DO
    END DO
! the routine does its own evaluation for DOUBLE subscripting of
! array a.
    ipivc = 1 - idim
! main loop TO invert the matrix
    DO main = 1, n
      pivot = 0D0
      ipivc = ipivc + idim
! search for next pivot IN column main.
      ipivc1 = ipivc + main - 1
      ipivc2 = ipivc + nmin1
      DO i1 = ipivc1, ipivc2
        IF (ABS(a(i1)) - ABS(pivot) <= 0D0) CYCLE
        pivot = a(i1)
        lpiv = i1
      END DO
! is pivot different from zero
      IF (pivot == 0D0) THEN
        GO TO 15
      END IF
! get the pivot-line indicator and swap lines IF necessary
3     icol = lpiv - ipivc + 1
      INDEX(main) = icol
      IF (icol - main > 0) THEN
        GO TO 4
      ELSE
        GO TO 6
      END IF
! complement the determinant
4     deter = -deter
! POINTER TO line pivot found
      icol = icol - idim
! POINTER TO exact pivot line
      i3 = main - idim
      DO i = 1, iemat
        icol = icol + idim
        i3 = i3 + idim
        swap = a(i3)
        a(i3) = a(icol)
        a(icol) = swap
      END DO
! compute determinant
6     deter = deter*pivot
      pivot = 1./pivot
! transform pivot column
      i3 = ipivc + nmin1
      DO i = ipivc, i3
        a(i) = -a(i)*pivot
      END DO
      a(ipivc1) = pivot
! pivot element transformed

! now convert rest of the matrix
      i1 = main - idim
! POINTER TO pivot line elements
      icol = 1 - idim
! general column POINTER
      DO i = 1, iemat
        icol = icol + idim
        i1 = i1 + idim
! pointers moved
        IF (i - main == 0) CYCLE
! pivot column excluded
        jcol = icol + nmin1
        swap = a(i1)
        i3 = ipivc - 1
        DO i2 = icol, jcol
          i3 = i3 + 1
          a(i2) = a(i2) + swap*a(i3)
        END DO
        a(i1) = swap*pivot
      END DO
    END DO
! now rearrange the matrix TO get right invers
    DO i1 = 1, n
      main = n + 1 - i1
      lpiv = INDEX(main)
      IF (lpiv - main == 0) CYCLE
12    icol = (lpiv - 1)*idim + 1
      jcol = icol + nmin1
      ipivc = (main - 1)*idim + 1 - icol
      DO i2 = icol, jcol
        i3 = i2 + ipivc
        swap = a(i2)
        a(i2) = a(i3)
        a(i3) = swap
      END DO
    END DO
    determ = deter
    nerror = 0
! ruecksetzen der gepackten eindimensionalen FORM IN den
! nmax*nmax-anteil der matrix c

    DO i = idim, 2, -1
      DO j = idim, 1, -1
        ipack = (i - 1)*nmax + j
        imatr = (i - 1)*ndim + j
        a(imatr) = a(ipack)
      END DO
    END DO

    RETURN
15  nerror = main
    determ = deter
    RETURN
  END SUBROUTINE matin2

END MODULE symcond

#else

FUNCTION symmdummy()
  RETURN
END FUNCTION symmdummy

#endif
