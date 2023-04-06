! integrate 1/r IN rectangular box

IMPLICIT DOUBLE PRECISION(a - h, o - z)

DO iter = 2, 11
  ngrid = 2**iter

  acc = 0D0
  dx = 1D0/ngrid
  DO ix = 1, ngrid
    x2 = (ix - 0.5D0)*dx
    x2 = x2*x2
    DO iy = 1, ngrid
      y2 = (iy - 0.5D0)*dx
      y2 = y2*y2
      DO iz = 1, ngrid
        z = (iz - 0.5D0)*dx
        acc = acc + 1D0/SQRT(x2 + y2 + z*z)
      END DO
    END DO
  END DO
  WRITE (*, *) ' ngrid,integral=', ngrid, acc*dx**3

END DO
STOP
END
