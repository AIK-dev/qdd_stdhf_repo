!-----givens------------------------------------------------------------------

SUBROUTINE givens(a, root, vect, nx, nrootx, njx)

! 62.3 givens -eigenvalues and eigenvectors by the givens method.
! by Franklin Prosser, Indiana University.
! September, 1967
! calculates eigenvalues and eigenvectors of REAL symmetric matrix
! stored in packed upper triangular FORM.

! thanks are due to F. E. Harris (Stanford University) and H. H.
! Michels (United Aircraft Research Laboratories) for excellent
! work on numerical difficulties with earlier versions of this
! program.

! the parameters for the routine are...
! nx order of matrix
! nrootx number of roots wanted. the nrootx smallest (most
! negative) roots will be calculated. IF no vectors
! are wanted, make this number negative.
! njx row DIMENSION of vect array. see 'vect' below.
! njx must be not less than nx.
! a matrix stored by columns in packed upper triangular
! FORM, i.e. occupying nx*(nx+1)/2 consecutive
! locations.
! b scratch array used by givens. must be at least
! nx*5 cells.
! root array to hold the eigenvalues. must be at least
! nrootx cells long. the nrootx smallest roots are
! ordered largest first in this array.
! vect eigenvector array. each column will hold an
! eigenvector for the corresponding root. must be
! dimensioned with 'njx' rows and at least 'nrootx'
! columns, unless no vectors
! are requested (negative nrootx). in this latter
! CASE, the argument vect is just a dummy, and the
! storage is not used.

! the arrays a and b are destroyed by the computation. the results
! appear in root and vect.

! for proper functioning of this routine, the RESULT of a floating
! point underflow should be a zero.

! to convert this routine to DOUBLE PRECISION (e.g. on ibm 360
! machines), be sure that all REAL variables and function
! references are properly made DOUBLE PRECISION.

! the original reference to the givens technique is in oak ridge
! report number ornl 1574 (physics), by wallace givens.
! the method as presented in this program consists of four steps,
! all modifications of the original method...
! first, the input matrix is reduced to tridiagonal FORM by the
! householder technique (J. H. Wilkinson, Comp. J. 3, 23 (1960)).
! the roots are then located by the Sturm sequence method (j. m.
! Ortega (see reference below). the vectors of the tridiagonal
! FORM are then evaluated (J. H. Wilkinson, Comp. J. 1, 90 (1958)),
! and last the tridiagonal vectors are rotated to vectors of the
! original array (first reference).
! vectors for degenerate (or near-degenerate) roots are forced
! to be orthogonal, using a method suggested by B. Garbow, Argonne
! National Labs (private communication, 1964). the Gram-Schmidt
! process is used for the orthogonalization.

! an excellent presentation of the givens technique is found IN
! J. M. Ortega's article in 'mathematics for digital computers,'
! volume 2, ed. by Ralston and Wilf, Wiley (1967), page 94.

! ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! ** users please note...
! ** the following two parameters, eta and theta, should be adjusted
! ** by the user for his particular machine.
! ** eta is an indication of the PRECISION of the floating point
! ** representation on the computer being used (roughly 10**(-m),
! ** WHERE m is the number of decimals of PRECISION ).
! ** theta is an indication of the range of numbers that can be
! ** expressed in the floating point representation (roughly the
! ** largest number or the reciprocal of the smallest NUMBER,
! ** whichever is smaller).
! ** some recommended values follow.
! ** for ibm 7094, univac 1108, etc. (27-bit binary fraction, 8-bit
! ** binary exponent), eta=1.e-8, theta=1.e37.
! ** for control DATA 3600 (36-bit binary fraction, 11-bit binary
! ** exponent), eta=1.e-11, theta=1.e307.
! ** for control DATA 6600 (48-bit binary fraction, 11-bit binary
! ** exponent), eta=1.e-14, theta=1.e295.
! ** for ibm 360/50 and 360/65 DOUBLE PRECISION (56-bit hexadecimal
! ** fraction, 7-bit hexadecimal exponent), eta=1.e-16, theta=1.e75.
! **
! IMPLICIT DOUBLE PRECISION (a-h,o-z)

  USE params, ONLY: DP
  IMPLICIT NONE
!INCLUDE 'com.inc'

  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: nrootx
  INTEGER, INTENT(IN) :: njx
  REAL(DP), INTENT(IN OUT) :: a(*)
  REAL(DP), INTENT(OUT) :: root(nrootx)
  REAL(DP), INTENT(OUT) :: vect(njx, nrootx)

  INTEGER, PARAMETER :: n1 = 300

  INTEGER :: i, ia, ic, id, ii, im, iter, j, j1, jj, jp2, jump, k, l, n, nroot, nsize, nm1, nm2, nom, nomtch
  REAL(DP) :: acc, ak, anorm, alimit, aroot, bb, del1, delbig, delta
  REAL(DP) :: elim1, elim2, f0, factor, small, temp, toler, theta1, trial
  REAL(DP) :: rand1, rootl, rootx, rpower, rpow1
  REAL(DP) :: b(n1, 5)

  REAL(DP) :: eta = 5D-15, theta = 1D88

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! check size

  IF (nx > n1) STOP ' givens was called with nx too large '

  del1 = eta/100D0
  delta = eta**2*100D0
  small = eta**2/100D0
  delbig = theta*delta/1000D0
  theta1 = 1000D0/theta
! toler is a factor used to determine if two roots are close
! enough to be considered degenerate for purposes of orthogonali-
! zing their vectors. for the matrix normed to unity, if the
! difference between two roots is less than toler, then
! orthogonalization will occur.
  toler = eta*100D0

! initial value for pseudorandom number generator... (2**23)-3
  rpower = 8388608D0
  rpow1 = rpower/2D0
  rand1 = rpower - 3D0

  n = nx
  nroot = ABS(nrootx)
  IF (nroot == 0) RETURN
! for n<1 or n=1 : no solution or trivial solution
  IF (n - 1 < 0) THEN
    RETURN
  ELSE IF (n - 1 == 0) THEN
    root(1) = a(1)
    IF (nrootx > 0) vect(1, 1) = 1D0
    RETURN
  END IF

! nsize number of elements in the packed array
  nsize = (n*(n + 1))/2
  nm1 = n - 1
  nm2 = n - 2

! scale matrix to euclidean norm of 1. scale factor is anorm.
  factor = 0D0
  DO i = 1, nsize
    factor = MAX(factor, ABS(a(i)))
  END DO
  IF (ABS(factor) <= 0D0) THEN
! null matrix. fix up roots and vectors, THEN EXIT.
    DO i = 1, nroot
      IF (nrootx >= 0) THEN
        DO j = 1, n
          vect(j, i) = 0D0
        END DO
        vect(i, i) = 1D0
      END IF
      root(i) = 0D0
    END DO
    RETURN
  END IF

  anorm = 0D0
  j = 1
  k = 1
  DO i = 1, nsize
    IF (i == j) THEN
      anorm = anorm + (a(i)/factor)**2/2D0
      k = k + 1
      j = j + k
      CYCLE
    END IF
81  anorm = anorm + (a(i)/factor)**2
  END DO
83 anorm = SQRT(anorm*2D0)*factor
  DO i = 1, nsize
    a(i) = a(i)/anorm
  END DO
  alimit = 1D0

! tridia section.
! tridiagonalization of symmetric matrix
  id = 0
  ia = 1
  IF (nm2 /= 0) THEN
    DO j = 1, nm2
! j counts row of a-matrix to be diagonalized
! ia start of non-codiagonal elements in the row
! id index of codiagonal element on row being codiagonalized.
      ia = ia + j + 2
      id = id + j + 1
      jp2 = j + 2
! sum squares of non-codiagonal elements in row j
      ii = ia
      acc = 0D0
      DO i = jp2, n
        acc = acc + a(ii)**2
        ii = ii + i
      END DO
      temp = a(id)
      IF (acc <= small) THEN
! no transformation necessary IF all the non-codiagonal
! elements are tiny.
        b(j, 1) = temp
        a(id) = 0D0
        CYCLE
      END IF
! now complete the sum of off-diagonal squares
      acc = SQRT(acc + temp**2)
! new codiagonal element
      b(j, 1) = -SIGN(acc, temp)
! first non-zero element of this w-vector
      b(j + 1, 2) = SQRT((1D0 + ABS(temp)/acc)/2D0)
! form rest of the w-vector elements
      temp = SIGN(0.5D0/(b(j + 1, 2)*acc), temp)
      ii = ia
      DO i = jp2, n
        b(i, 2) = a(ii)*temp
        ii = ii + i
      END DO
! form p-vector and scalar. p-vector = a-matrix*w-vector.
! scalar = w-vector*p-vector.
      ak = 0D00
! ic location of next diagonal element
      ic = id + 1
      j1 = j + 1
      DO i = j1, n
        jj = ic
        temp = 0D0
        DO ii = j1, n
! i runs over the non-zero p-elements
! ii runs over elements of w-vector
          temp = temp + b(ii, 2)*a(jj)
! change incrementing mode at the diagonal elements.
          IF (ii >= i) THEN
            jj = jj + ii
            CYCLE
          END IF
          jj = jj + 1
        END DO
! build up the k-scalar (ak)
        ak = ak + temp*b(i, 2)
        b(i, 1) = temp
! move ic to top of next a-matrix 'row'
        ic = ic + i
      END DO
! form the q-ve+tor
      DO i = j1, n
        b(i, 1) = b(i, 1) - ak*b(i, 2)
      END DO
! transform the rest of the a-matrix
! jj start-1 of the rest of the a-matrix
      jj = id
! move w-vector into the old a-matrix locations to SAVE space
! i runs over the significant elements of the w-vector
      DO i = j1, n
        a(jj) = b(i, 2)
        DO ii = j1, i
          jj = jj + 1
          a(jj) = a(jj) - 2D0*(b(i, 1)*b(ii, 2) + b(i, 2)*b(ii, 1))
        END DO
        jj = jj + j
      END DO
    END DO
  END IF
! move last codiagonal element OUT into its proper place
  b(nm1, 1) = a(nsize - 1)
  a(nsize - 1) = 0D0
  b(n, 1) = 0D0

! Sturm section.
! Sturm sequence iteration to obtain roots of tridiagonal form.
! move diagonal elements into second n elements of b-vector.
! this is a more convenient indexing position.
! also, put square of codiagonal elements in third n elements.
  jump = 1
  DO j = 1, n
    b(j, 2) = a(jump)
    b(j, 3) = b(j, 1)**2
    jump = jump + j + 1
  END DO
  DO i = 1, nroot
    root(i) = +alimit
  END DO
  rootl = -alimit
! isolate the roots. the nroot lowest roots are found, lowest first
  DO i = 1, nroot
! find current 'best' upper bound
    rootx = +alimit
    DO j = i, nroot
      rootx = MIN(rootx, root(j))
    END DO
    root(i) = rootx
! get improved trial root
    DO
      trial = (rootl + root(i))*0.5D0
! IF (trial.eq.rootl.or.trial.eq.root(i)) go to 330
      IF (.NOT. ABS(trial - rootl) > 0D0 .OR. .NOT. ABS(trial - root(i)) > 0D0) EXIT
! form Sturm sequence ratios, using Ortega's algorithm (modified).
! nomtch is the number of roots less than the trial value.
      nomtch = n
      j = 1
      jloop: DO WHILE (j <= n)
        f0 = b(j, 2) - trial
        DO WHILE (ABS(f0) >= theta1)
          IF (f0 >= 0D0) nomtch = nomtch - 1
          j = j + 1
          IF (j > n) EXIT jloop
! since matrix is normed to unity, magnitude of b(j,3) is less than
! one, so overflow is not possible at the division step, since
! f0 is greater than theta1.
! f0 = b(j,2) - trial - b(j-1,3)/f0
          bb = b(j - 1, 3)/f0
! IF (f0.eq.0) bb=0D0
          IF (.NOT. ABS(f0) > 0D0) bb = 0D0
          f0 = b(j, 2) - trial - bb
        END DO
        j = j + 2
        nomtch = nomtch - 1
      END DO jloop
! fix new bounds on roots
      IF (nomtch < i) THEN
        rootl = trial
      ELSE
        root(i) = trial
        nom = MIN0(nroot, nomtch)
        root(nom) = trial
      END IF
    END DO
  END DO
! reverse the order of the eigenvalues, since custom dictates
! 'largest first'. this section may be removed IF desired without
! affecting the remainder of the routine.
! nrt = nroot/2
! DO 10 i=1,nrt
! SAVE = root(i)
! nmip1 = nroot - i + 1
! root(i) = root(nmip1)
! 10 root(nmip1) = SAVE

! trivec section.
! eigenvectors of codiagonal form

! quit now IF no vectors were requested.
  IF (nrootx < 0) THEN
    root(1:nroot) = anorm*root(1:nroot)
    RETURN
  END IF
! initialize vector array.
  DO i = 1, n
    DO j = 1, nroot
      vect(i, j) = 1D0
    END DO
  END DO
  DO i = 1, nroot
    aroot = root(i)
! orthogonalize if roots are close.
! the absolute value in the next test is to assure that the trivec
! section is independent of the order of the eigenvalues.
    IF (i /= 1) THEN
      IF (ABS(root(i - 1) - aroot) < toler) THEN
        ia = ia + 1
      ELSE
        ia = 0
      END IF
    ELSE
      ia = 0
    END IF

    elim1 = a(1) - aroot
    elim2 = b(1, 1)
    jump = 1
    DO j = 1, nm1
      jump = jump + j + 1
! get the correct pivot equation for this step.
      IF (ABS(elim1) > ABS(b(j, 1))) THEN
! first (elim1) equation is the pivot this time. CASE 1.
        b(j, 2) = elim1
        b(j, 3) = elim2
        b(j, 4) = 0D0
        temp = b(j, 1)/elim1
        elim1 = a(jump) - aroot - temp*elim2
        elim2 = b(j + 1, 1)
      ELSE
! second equation is the pivot this time. CASE 2.
        b(j, 2) = b(j, 1)
        b(j, 3) = a(jump) - aroot
        b(j, 4) = b(j + 1, 1)
        temp = 1D0
        IF (ABS(b(j, 1)) > theta1) temp = elim1/b(j, 1)
        elim1 = elim2 - temp*b(j, 3)
        elim2 = -temp*b(j + 1, 1)
      END IF
! save factor for the second iteration.
      b(j, 5) = temp
    END DO
    b(n, 2) = elim1
    b(n, 3) = 0D0
    b(n, 4) = 0D0
    b(nm1, 4) = 0D0
    iter = 1
    IF (ia /= 0) THEN
! produce a random vector
      DO j = 1, n
! generate pseudorandom numbers with uniform distribution in (-1,1).
! this random number scheme is of the form...
! rand1 = amod((2**12+3)*rand1,2**23)
! it has a period of 2**21 numbers.
        rand1 = MOD(4099D0*rand1, rpower)
        vect(j, i) = rand1/rpow1 - 1D0
      END DO
    END IF
! back substitute to get this vector.
    eternal: DO
      l = n + 1
      DO j = 1, n
        l = l - 1
        infinite: DO
          IF (l < nm1) THEN
            elim1 = vect(l, i) - vect(l + 1, i)*b(l, 3) - vect(l + 2, i)*b(l, 4)
          ELSE IF (l == nm1) THEN
            elim1 = vect(l, i) - vect(l + 1, i)*b(l, 3)
          ELSE
            elim1 = vect(l, i)
          END IF
          ! IF overflow is conceivable, scale the vector down.
          ! this approach is used to avoid machine-dependent and system-
          ! dependent calls to overflow routines.
          IF (ABS(elim1) > delbig) THEN
            ! vector is too big. scale it down.
            DO k = 1, n
              vect(k, i) = vect(k, i)/delbig
            END DO
          ELSE
            temp = b(l, 2)
            IF (ABS(b(l, 2)) < delta) temp = delta
            vect(l, i) = elim1/temp
            EXIT infinite
          END IF
        END DO infinite
      END DO
      IF (iter == 2) EXIT eternal
      ! second iteration. (both iterations for repeated-root vectors).
      iter = iter + 1
      elim1 = vect(1, i)
      DO j = 1, nm1
        IF (ABS(b(j, 2) - b(j, 1)) > 0D0) THEN
          ! CASE one.
          vect(j, i) = elim1
          elim1 = vect(j + 1, i) - elim1*b(j, 5)
        ELSE
          ! CASE two.
          vect(j, i) = vect(j + 1, i)
          elim1 = elim1 - vect(j + 1, i)*temp
        END IF
      END DO
      vect(n, i) = elim1
    END DO eternal

! orthogonalize this repeated-root vector to others with this root.
    IF (ia /= 0) THEN
      DO j1 = 1, ia
        k = i - j1
        temp = 0D0
        DO j = 1, n
          temp = temp + vect(j, i)*vect(j, k)
        END DO
        DO j = 1, n
          vect(j, i) = vect(j, i) - temp*vect(j, k)
        END DO
      END DO
    END IF
! normalize the vector
    elim1 = 0D0
    DO j = 1, n
      elim1 = MAX(ABS(vect(j, i)), elim1)
    END DO
    temp = 0D0
    DO j = 1, n
      elim2 = vect(j, i)/elim1
      temp = temp + elim2**2
    END DO
    temp = 1D0/(SQRT(temp)*elim1)
    DO j = 1, n
      vect(j, i) = vect(j, i)*temp
      IF (ABS(vect(j, i)) < del1) vect(j, i) = 0D0
    END DO
  END DO

! simvec section.
! rotate codiagonal vectors into vectors of original array
! loop over all the transformation vectors
  IF (nm2 == 0) THEN
    root = anorm*root
    RETURN
  END IF

  jump = nsize - (n + 1)
  im = nm1
  DO i = 1, nm2
    j1 = jump
! move a transformation vector OUT into better indexing POSITION.
    DO j = im, n
      b(j, 2) = a(j1)
      j1 = j1 + j
    END DO
! modify all requested vectors.
    DO k = 1, nroot
      temp = 0D0
! form scalar product of transformation vector with eigenvector
      DO j = im, n
        temp = temp + b(j, 2)*vect(j, k)
      END DO
      temp = temp + temp
      DO j = im, n
        vect(j, k) = vect(j, k) - temp*b(j, 2)
      END DO
    END DO
    jump = jump - im
    im = im - 1
  END DO
! restore roots to their proper SIZE.
  root = anorm*root
  RETURN
END SUBROUTINE givens
