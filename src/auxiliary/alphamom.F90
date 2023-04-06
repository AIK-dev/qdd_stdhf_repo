!--alphamom--------------------------------------------------

! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16 Time: 12:34:36

!
PROGRAM alphamom

  IMPLICIT NONE

! computes dimensionless quadrupole moments from <x**2>
! <y**2> and <z**2>

  REAL :: pi
  PARAMETER(pi=3.1415926E0)

  REAL :: xx, yy, zz, rms2, q20, q22, fac20, fac22, alp20, alp22, beta, gamma

! READ spatial variances from standard input

  WRITE (6, *) ' enter the three moments <x**2>, <y**2> and <z**2>:'
  READ (5, *) xx, yy, zz

! reshuffle TO multipole moments

  rms2 = xx + yy + zz
  q20 = 2.0*zz - xx - yy
  q22 = xx - yy

! rephrase into rms and dimensionless quadrupole moments

  fac20 = (4.0*pi/5.0)*SQRT(5.0/(16.0*pi)) !/PI
  fac22 = (4.0*pi/5.0)*SQRT(15.0/(32.0*pi)) !/PI
  alp20 = fac20*q20/rms2
  alp22 = fac22*q22/rms2
  beta = SQRT(alp20*alp20 + 2.0*alp22*alp22)
  gamma = ABS(ATAN(SQRT(2.0)*alp22/alp20)*180.0/pi)
  IF (q20 < 0.0E0) THEN
    gamma = 180.0E0 - gamma
    gamma = ABS(gamma - 120.0E0)
  END IF
  IF (gamma > 60.0E0) THEN
    gamma = 60.0E0 - (gamma - 60.0E0)
  END IF

  WRITE (6, *) ' alpha20,alpha22,beta2,gamma=', alp20, alp22, beta, gamma

  STOP
END PROGRAM alphamom

