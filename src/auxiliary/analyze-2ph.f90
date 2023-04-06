IMPLICIT REAL(8) (a - h, o - z)

REAL(8), PARAMETER :: emin = 0.21, emax = 0.67, dele = 0.1, wide = 3*dele
REAL(8), PARAMETER :: rmin = 1D-2, rmax = 1D2, delrl = 0.1, widrl = 3*delrl
REAL(8), PARAMETER :: vmin = 1D-7, vmax = 1D-2, delv = 0.1, widv = 3*delv

REAL(8), ALLOCATABLE :: dens(:, :), dv(:, :), dens0(:, :), dv0(:, :)
INTEGER :: idum(8)

rlmin = LOG10(rmin)
rlmax = LOG10(rmax)
vlmin = LOG10(vmin)
vlmax = LOG10(vmax)

nenerg = NINT((emax - emin)/dele)
ndens = NINT((rlmax - rlmin)/delrl)
nv = NINT((vlmax - vlmin)/delv)

ALLOCATE (dens(ndens, nenerg), dens0(ndens, nenerg))
ALLOCATE (dv(nv, nv), dv0(nv, nv))
dens = 0D0
dens0 = 0D0
dv = 0D0
dv0 = 0D0

READ (5, '(1x/1x)')
DO i = 1, 999999 !4092
  READ (5, *, ERR=19, END=19) idum, e, v1, v2, v2b

! ratio versus energy
  IF (v2 > 1D-10 .AND. v1 > 1D-10) THEN
    ql = LOG10(v1/v2)
    DO ie = 1, nenerg
      en = emin + (ie - 1)*dele
      DO id = 1, ndens
        rl = rlmin + delrl*(id - 1)
        dens(id, ie) = dens(id, ie) + EXP(-((en - e)/wide)**2 - ((rl - ql)/widrl)**2)
      END DO
    END DO
  END IF
! v2 versus v1
  IF (v2 > 1D-10 .AND. v1 > 1D-10) THEN
    DO i1 = 1, nv
      vl1 = vlmin + (i1 - 1)*delv
      DO i2 = 1, nv
        vl2 = vlmin + (i2 - 1)*delv
        dv(i2, i1) = dv(i2, i1) + &
                     EXP(-((vl1 - LOG10(v1))/widv)**2 - ((vl2 - LOG10(v2))/widv)**2)
      END DO
    END DO
  END IF
! ratio TO v2b versus energy
  IF (v2b > 1D-10 .AND. v1 > 1D-10) THEN
    ql = LOG10(v1/v2b)
    DO ie = 1, nenerg
      en = emin + (ie - 1)*dele
      DO id = 1, ndens
        rl = rlmin + delrl*(id - 1)
        IF (v2 == 0D0) dens0(id, ie) = dens0(id, ie) + EXP(-((en - e)/wide)**2 - ((rl - ql)/widrl)**2)
      END DO
    END DO
  END IF
! v2b versus v1
  IF (v2b > 1D-10 .AND. v1 > 1D-10) THEN
    DO i1 = 1, nv
      vl1 = vlmin + (i1 - 1)*delv
      DO i2 = 1, nv
        vl2 = vlmin + (i2 - 1)*delv
        IF (v2 == 0D0) dv0(i2, i1) = dv0(i2, i1) + &
                                     EXP(-((vl1 - LOG10(v1))/widv)**2 - ((vl2 - LOG10(v2b))/widv)**2)
      END DO
    END DO
  END IF
END DO

19 CONTINUE

WRITE (6, *) '# ratio VYuk/Vdelta versus energy'
DO ie = 1, nenerg
  en = emin + (ie - 1)*dele
  DO id = 1, ndens
    rl = rlmin + delrl*(id - 1)
    WRITE (6, '(3(1pg13.5))') rl, en, dens(id, ie)
  END DO
  WRITE (6, '(1x)')
END DO

WRITE (6, '(1x)')

WRITE (6, *) '# Vdelta versus Vyuk'
DO i1 = 1, nv
  vl1 = vlmin + (i1 - 1)*delv
  DO i2 = 1, nv
    vl2 = vlmin + (i2 - 1)*delv
    WRITE (6, '(3(1pg13.5))') vl1, vl2, dv(i2, i1)
  END DO
  WRITE (6, '(1x)')
END DO

WRITE (6, '(1x)')

WRITE (6, *) '# ratio VYuk/Vdelta0 versus energy'
DO ie = 1, nenerg
  en = emin + (ie - 1)*dele
  DO id = 1, ndens
    rl = rlmin + delrl*(id - 1)
    WRITE (6, '(3(1pg13.5))') rl, en, dens0(id, ie)
  END DO
  WRITE (6, '(1x)')
END DO

WRITE (6, '(1x)')

WRITE (6, *) '# Vdelta0 versus Vyuk'
DO i1 = 1, nv
  vl1 = vlmin + (i1 - 1)*delv
  DO i2 = 1, nv
    vl2 = vlmin + (i2 - 1)*delv
    WRITE (6, '(3(1pg13.5))') vl1, vl2, dv0(i2, i1)
  END DO
  WRITE (6, '(1x)')
END DO

STOP
END

