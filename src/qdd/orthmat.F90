
MODULE orthmat

! Subroutines and functions for matrices and vectors overlaps, as well
! as their orthonormalization.
! Most procedures exist in real and complex version, depending on the
! TYPE of their arguments (regulated by INTERFACE).
  USE params
  IMPLICIT NONE

  INTERFACE matdorth
    MODULE PROCEDURE matdorth_r, matdorth_c
  END INTERFACE matdorth

  INTERFACE matunite
    MODULE PROCEDURE matunite_r, matunite_c
  END INTERFACE matunite

  INTERFACE matexp
    MODULE PROCEDURE matexp_r, matexp_c
  END INTERFACE matexp

  INTERFACE orthnorm
    MODULE PROCEDURE orthnorm_r, orthnorm_c
  END INTERFACE orthnorm

  INTERFACE vecovlp
    MODULE PROCEDURE vecovlp_r, vecovlp_c
  END INTERFACE vecovlp

  INTERFACE vecnorm
    MODULE PROCEDURE vecnorm_r, vecnorm_c
  END INTERFACE vecnorm

  INTERFACE ovlpmatrix
    MODULE PROCEDURE ovlpmatrix_r, ovlpmatrix_c
  END INTERFACE ovlpmatrix

CONTAINS

!_________________________MatdOrth___________________________
!
! Matrix norm: Sum_ij [A(i,:)*A(j,:)-delta_ij]

!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  REAL(DP) FUNCTION matdorth_r(aa, n, ndim)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: aa(n, n)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: i, j
    matdorth_r = 0D0
    DO i = 1, ndim
      DO j = i, ndim
        IF (i == j) THEN
          matdorth_r = matdorth_r + SUM(aa(i, 1:ndim)*aa(j, 1:ndim)) - 1D0
        ELSE
          matdorth_r = matdorth_r + SUM(aa(i, 1:ndim)*aa(j, 1:ndim))
        END IF
      END DO
    END DO
  END FUNCTION matdorth_r

!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
  REAL(DP) FUNCTION matdorth_c(aa, n, ndim)

    INTEGER, INTENT(IN) :: n
    COMPLEX(DP), INTENT(IN) :: aa(n, n)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: i, j
    matdorth_c = 0D0
    DO i = 1, ndim
      DO j = i, ndim
        IF (i == j) THEN
          matdorth_c = matdorth_c + SUM(aa(i, 1:ndim)*CONJG(aa(j, 1:ndim))) - 1D0
        ELSE
          matdorth_c = matdorth_c + SUM(aa(i, 1:ndim)*CONJG(aa(j, 1:ndim)))
        END IF
      END DO
    END DO
  END FUNCTION matdorth_c

!____________________________rMat Unite______________________________
!
! Set matrix TO UNIT matrix: AA_ij=delta_ij
!
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  SUBROUTINE matunite_r(aa, n, ndim)

    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(OUT) :: aa(n, n)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: i
    aa(1:ndim, 1:ndim) = 0D0
    DO i = 1, ndim
      aa(i, i) = 1D0
    END DO

  END SUBROUTINE matunite_r

!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
  SUBROUTINE matunite_c(aa, n, ndim)

    INTEGER, INTENT(IN) :: n
    COMPLEX(DP), INTENT(OUT) :: aa(n, n)
    INTEGER, INTENT(IN) :: ndim
    INTEGER :: i
    aa(1:ndim, 1:ndim) = CMPLX(0D0, 0D0, DP)
    DO i = 1, ndim
      aa(i, i) = CMPLX(1D0, 0D0, DP)
    END DO

  END SUBROUTINE matunite_c

!__________________________rMat Exp_____________________________
!
! Matrix exponential: BB=Exp(AA)
!
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  SUBROUTINE matexp_r(aa, bb, n, ndim)

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: ndim
    REAL(DP), INTENT(IN OUT) :: aa(n, n)
    REAL(DP), INTENT(IN OUT) :: bb(n, n)
    REAL(DP) cc(n, n), dd(n, n)

    REAL(DP) eps, delta, rn

    CALL matunite(bb, n, ndim)
    CALL matunite(cc, n, ndim)
    eps = 1D-20
    rn = 1D0
    delta = 1D0
    DO WHILE (delta > eps)
      dd(1:ndim, 1:ndim) = MATMUL(aa(1:ndim, 1:ndim), cc(1:ndim, 1:ndim))
      cc(1:ndim, 1:ndim) = dd(1:ndim, 1:ndim)/rn
      delta = SQRT(SUM(cc(1:ndim, 1:ndim)**2))
      bb(1:ndim, 1:ndim) = cc(1:ndim, 1:ndim) + bb(1:ndim, 1:ndim)
      rn = rn + 1D0
    END DO

  END SUBROUTINE matexp_r

!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
  SUBROUTINE matexp_c(aa, bb, n, ndim)

    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN OUT) :: ndim
    COMPLEX(DP), INTENT(IN OUT) :: aa(n, n)
    COMPLEX(DP), INTENT(IN OUT) :: bb(n, n)

    COMPLEX(DP) :: cc(n, n), dd(n, n)
    REAL(DP) eps, delta, rn

    eps = 1D-20
    CALL matunite(bb, n, ndim)
    CALL matunite(cc, n, ndim)
    rn = CMPLX(1D0, 0D0, DP)
    delta = 1D0

    DO WHILE (delta > eps)
      dd(1:ndim, 1:ndim) = MATMUL(aa(1:ndim, 1:ndim), cc(1:ndim, 1:ndim))
      cc(1:ndim, 1:ndim) = dd(1:ndim, 1:ndim)/rn
      delta = SQRT(SUM(cc(1:ndim, 1:ndim)*CONJG(cc(1:ndim, 1:ndim))))
      bb(1:ndim, 1:ndim) = cc(1:ndim, 1:ndim) + bb(1:ndim, 1:ndim)
      rn = rn + 1
    END DO
  END SUBROUTINE matexp_c

!________________________OrthNorm__________________________
!
! Ortho-normalizes system of vectors 'vecs' with actual
! dimension 'ndim' ('kdim' is dimension in calling routine).
!
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  SUBROUTINE orthnorm_r(vecs, ndim, kdim)
    USE params, ONLY: DP

    INTEGER, INTENT(IN) :: ndim, kdim
    REAL(DP), INTENT(IN OUT) :: vecs(kdim, kdim)

    INTEGER :: i, j
    REAL(DP) :: acc2
    REAL(DP) :: acc

!-------------------------------------------------------
    DO i = 1, ndim
      DO j = 1, i - 1
        acc = vecovlp(vecs(:, j), vecs(:, i), ndim)
        vecs(1:ndim, i) = vecs(1:ndim, i) - acc*vecs(1:ndim, j)
      END DO
      acc2 = 1D0/SQRT(vecnorm(vecs(:, i), ndim))
      vecs(1:ndim, i) = vecs(1:ndim, i)*acc2
    END DO

    RETURN
  END SUBROUTINE orthnorm_r
!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
  SUBROUTINE orthnorm_c(vecs, ndim, kdim)
    USE params, ONLY: DP

    INTEGER, INTENT(IN) :: ndim, kdim
    COMPLEX(DP), INTENT(IN OUT) :: vecs(kdim, kdim)

    INTEGER :: i, j
    REAL(DP) :: acc2
    COMPLEX(DP) :: acc

!-------------------------------------------------------
    DO i = 1, ndim
      DO j = 1, i - 1
        acc = vecovlp(vecs(:, j), vecs(:, i), ndim)
        vecs(1:ndim, i) = vecs(1:ndim, i) - acc*vecs(1:ndim, j)
      END DO
      acc2 = 1D0/SQRT(vecnorm(vecs(:, i), ndim))
      vecs(1:ndim, i) = vecs(1:ndim, i)*acc2
    END DO

    RETURN
  END SUBROUTINE orthnorm_c

!___________________vecovlp__________________________________
!
! Overlap 'vec1' with 'vec2' having active dimension 'ndim'.
!
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  REAL(DP) FUNCTION vecovlp_r(vec1, vec2, ndim)

    USE params, ONLY: DP
    INTEGER, INTENT(IN) :: ndim
    REAL(DP), INTENT(IN) :: vec1(ndim), vec2(ndim)

    vecovlp_r = SUM(vec1(1:ndim)*vec2(1:ndim))

    RETURN
  END FUNCTION vecovlp_r
!-------------------------------------------------------
!COMPLEX version
!-------------------------------------------------------
  COMPLEX(DP) FUNCTION vecovlp_c(vec1, vec2, ndim)

    USE params, ONLY: DP
    INTEGER, INTENT(IN) :: ndim
    COMPLEX(DP), INTENT(IN) :: vec1(ndim), vec2(ndim)

    vecovlp_c = SUM(CONJG(vec1(1:ndim))*vec2(1:ndim))

    RETURN
  END FUNCTION vecovlp_c

!_______________________vecnorm__________________________
!
! Norm of vector 'vec' with active dimension 'ndim'.
!
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  REAL(DP) FUNCTION vecnorm_r(vec, ndim)
    USE params, ONLY: DP
    INTEGER, INTENT(IN) :: ndim
    REAL(DP), INTENT(IN) :: vec(ndim)

    vecnorm_r = SUM(vec(1:ndim)**2)

    RETURN
  END FUNCTION vecnorm_r
!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
  REAL(DP) FUNCTION vecnorm_c(vec, ndim)

    USE params, ONLY: DP
    INTEGER, INTENT(IN) :: ndim
    COMPLEX(DP), INTENT(IN) :: vec(ndim)

    vecnorm_c = SUM(ABS(vec(1:ndim)**2))

    RETURN
  END FUNCTION vecnorm_c

!_____________________ovlpmatrix_________________________
!
! Transition matrix elmenent <vec1|A|vec2> for matrix A
! and vectors vec1 and vec2.
!
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  REAL(DP) FUNCTION ovlpmatrix_r(a, vec1, vec2, ndim, kdim)
    USE params, ONLY: DP

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(IN) :: kdim
    REAL(DP), INTENT(IN) :: a(kdim, kdim)
    REAL(DP), INTENT(IN) :: vec1(kdim)
    REAL(DP), INTENT(IN) :: vec2(kdim)
    REAL(DP) :: ovlp
    INTEGER :: i, j
!-------------------------------------------------------

    ovlp = 0D0
    DO i = 1, ndim
      DO j = 1, ndim
        ovlp = ovlp + vec1(j)*a(j, i)*vec2(i)
      END DO
    END DO
    ovlpmatrix_r = ovlp

    RETURN
  END FUNCTION ovlpmatrix_r
!-------------------------------------------------------
!COMPLEX version
!-------------------------------------------------------
  COMPLEX(DP) FUNCTION ovlpmatrix_c(a, vec1, vec2, ndim, kdim)
    USE params, ONLY: DP

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(IN) :: kdim
    COMPLEX(DP), INTENT(IN) :: a(kdim, kdim)
    COMPLEX(DP), INTENT(IN) :: vec1(kdim)
    COMPLEX(DP), INTENT(IN) :: vec2(kdim)
    COMPLEX(DP) :: ovlp
    INTEGER :: i, j
!-------------------------------------------------------

    ovlp = CMPLX(0D0, 0D0, DP)
    DO i = 1, ndim
      DO j = 1, ndim
        ovlp = ovlp + CONJG(vec1(j))*a(j, i)*vec2(i)
      END DO
    END DO
    ovlpmatrix_c = ovlp

    RETURN
  END FUNCTION ovlpmatrix_c

END MODULE orthmat
