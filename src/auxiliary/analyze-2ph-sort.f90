IMPLICIT REAL(8) (a - h, o - z)

INTEGER, PARAMETER :: kmax = 5510
REAL(8) :: e(kmax), v(3, kmax), vs(3), vs0(2), vall(3), vall0(2)
INTEGER :: idum(8)

READ (5, '(1x/1x)')
DO i = 1, kmax
  READ (5, *, ERR=19, END=19) idum, ein, v1, v2, v2b
  e(i) = ein
  v(1, i) = v1
  v(2, i) = v2
  v(3, i) = v2b
  nmax = i
END DO

19 CONTINUE

DO i = 1, nmax
  DO j = i, nmax
    IF (e(j) < e(i)) THEN
      es = e(j)
      e(j) = e(i)
      e(i) = es
      vs = v(:, j)
      v(:, j) = v(:, i)
      v(:, i) = vs
    END IF
  END DO
END DO

nall = 0
nall0 = 0
iacc = 0
iacc0 = 0
es = 0D0
vs(:) = 0D0
vs0 = 0D0
vall = 0D0
vall0 = 0D0
DO i = 1, nmax
! WRITE(6,'(4(1pg13.5))') e(i),v(:,i)
  IF (es == e(i)) THEN
    vs(:) = vs(:) + v(:, i)
    IF (v(2, i) > 1D-8) THEN
      vs0 = vs0 + v(1:2, i)
      iacc0 = 1 + iacc0
    END IF
    iacc = 1 + iacc
  ELSE IF (i > 1) THEN
    WRITE (6, '(6(1pg13.5))') es, vs/iacc, vs0/iacc0
    iacc = 1
    es = e(i)
    vs(:) = v(:, i)
    IF (v(2, i) > 1D-8) THEN
      vs0 = v(1:2, i)
      iacc0 = 1
    END IF
  ELSE
    vs(:) = vs(:) + v(:, i)
    IF (v(2, i) > 1D-8) THEN
      vs0 = vs0 + v(1:2, i)
      iacc0 = 1 + iacc0
    END IF
    es = e(i)
    iacc = 1 + iacc
  END IF
  vall = vall + v(:, i)
  IF (v(2, i) > 1D-8) THEN
    vall0 = vall0 + v(1:2, i)
    nall0 = 1 + nall0
  END IF
  nall = 1 + nall
END DO
vall = vall/nall
vall0 = vall0/nall0
WRITE (6, '(a,5(1pg13.5)/a,2(1pg13.5))') '#', vall, vall0, &
  '#', vall(1)/vall(2), vall0(1)/vall0(2)

STOP
END
