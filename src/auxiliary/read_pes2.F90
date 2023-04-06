PROGRAM test
  IMPLICIT REAL*8(A - H, O - Z)
  INTEGER, PARAMETER :: nlines = 100000000
  INTEGER, PARAMETER :: nheader = 2 ! commentary lines
  INTEGER, PARAMETER :: ncomment = 9 ! commentary lines
  INTEGER, PARAMETER :: np = 18 ! the spacing between 0 TO 180 degrees
  INTEGER, DIMENSION(0:np) :: rc
  REAL*8, DIMENSION(0:np) :: pes, pes2
  REAL*8, DIMENSION(0:np) :: radius
  CHARACTER :: symbol

  OPEN (9, FILE="c60_e136_nz1001.ekinet") ! raw DATA
  OPEN (10, FILE="pes_e136_final.c60") ! with two head lines
  OPEN (11, FILE="pes_e136_rescaling.ek") ! new ouput for ARPES

  DO i0 = 1, ncomment
    READ (9, *)
  END DO

  DO i1 = 0, np
    READ (9, *) symbol, &
      xx, xx, xx, rc(i1), radius(i1), xx, xx
  END DO

  DO i = 1, nheader
    READ (10, *)
  END DO
  DO j = 1, nlines

    READ (10, *, ERR=100, END=100) omega, pes
    DO k = 0, np
      pes2(k) = pes(k)/(rc(k)*8.0) ! 1/(rc(k)*8.0) is the scaling factor
      IF ((k == 0) .OR. (k == np)) pes2(k) = pes(k) ! no scaling for polar points
    END DO
    WRITE (11, '(1pg14.6,18(1pg14.4))') omega, pes2
  END DO

100 CONTINUE
END
