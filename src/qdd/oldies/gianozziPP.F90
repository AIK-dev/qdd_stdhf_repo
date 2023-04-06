
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16 Time: 12:36:26

DOUBLE PRECISION FUNCTION vgian(r)

  DOUBLE PRECISION, INTENT(IN) :: r
  IMPLICIT NONE

  DOUBLE PRECISION :: erf
  EXTERNAL erf

! returns the Gianozzi hydrogen PP IN Rydberg units.
! Reference: F. Gygi, PRB 48, 11692 (1993).
! This PP was used for the calculations IN
! K"ummel, Kronik, Perdew, PRL 93, 213002 (2004)

  DOUBLE PRECISION ::
  DOUBLE PRECISION, PARAMETER :: e2 = 2.d0
  DOUBLE PRECISION, PARAMETER :: rc1 = 0.25D0
  DOUBLE PRECISION, PARAMETER :: rc2 = 0.284D0
  DOUBLE PRECISION, PARAMETER :: a = -1.9287D0*e2
  DOUBLE PRECISION, PARAMETER :: b = 0.3374D0*e2

  vgian = -e2*erf(r/rc1)/r + (a + b*r**2)*EXP(-(r/rc2)**2)

END FUNCTION vgian
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

DOUBLE PRECISION FUNCTION erf(x)

  DOUBLE PRECISION, INTENT(IN) :: x
  IMPLICIT NONE

! **************************************************************

! TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
! ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING CDD ,
! WHERE DD IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
! DIGITS REPRESENTED BY THE MACHINE
! DELETE THE ILLEGAL DUMMY STATEMENT OF THE FORM
! * EXPANSION (DATA) *

! **************************************************************
! .. Scalar Arguments ..

  INTEGER :: ifail
! .. Local Scalars ..
  DOUBLE PRECISION :: bj, bjp1, bjp2, half, one, sqrtpi, three, twenty, two, &
    x2, xup, xv, zero
  INTEGER :: j, ncfc, ncfd
! .. Local Arrays ..
!07 DOUBLE PRECISION C(8), D(8)
!12 DOUBLE PRECISION C(11), D(12)
!14 DOUBLE PRECISION C(15), D(15)
  DOUBLE PRECISION :: c(18), d(17)
! .. INTRINSIC Functions ..
  INTRINSIC ABS, EXP, SIGN
! .. DATA statements ..
! * EXPANSION (DATA) *
!07 DATA NCFC,NCFD/8,8/,XUP/4.0D0/,SQRTPI/1.772454D0/
!07 A,C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8)
!07 A/1.944907D0,4.2019D-2,-1.8687D-2,5.129D-3,-1.068D-3
!07 A,1.74D-4,-2.1D-5,2.0D-6/
!07 A,D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8)
!07 A/1.483110D0,-3.01071D-1,6.8995D-2,-1.3916D-2,2.421D-3
!07 A,-3.66D-4,4.9D-5,-6.0D-6/

!12 DATA NCFC,NCFD/11,12/,XUP/5.0D0/,SQRTPI/1.7724538509D0/,C(1),
!12 * C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11)
!12 * /1.9449071068D0,4.20186582D-2,-1.86866104D-2,5.1281062D-3,
!12 * -1.0683107D-3,1.744738D-4,-2.15642D-5,1.7283D-6,-2.D-8,-1.65D-8,
!12 * 2.D-9/,D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),
!12 * D(10),D(11),D(12)/1.4831105641D0,-3.010710734D-1,6.89948307D-2,
!12 * -1.39162713D-2,2.4207995D-3,-3.658640D-4,4.86210D-5,-5.7493D-6,
!12 * 6.113D-7,-5.90D-8,5.2D-9,-4.D-10/

!14 DATA NCFC,NCFD/15,15/,XUP/5.75D0/,SQRTPI/1.7724538509055D0/
!14 A,C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10)
!14 A,C(11),C(12),C(13),C(14),C(15)
!14 A/1.9449071068179D0,4.20186582324D-2,-1.86866103977D-2
!14 A,5.1281061839D-3,-1.0683107462D-3,1.744737872D-4
!14 A,-2.15642056D-5,1.7282658D-6,-2.00479D-8,-1.64782D-8
!14 A,2.0008D-9,2.58D-11,-3.06D-11,1.9D-12,4.0D-13/
!14 A,D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),D(10)
!14 A,D(11),D(12),D(13),D(14),D(15)
!14 A/1.4831105640848D0,-3.010710733866D-1,6.89948306898D-2
!14 A,-1.39162712647D-2,2.4207995224D-3,-3.658639686D-4
!14 A,4.86209844D-5,-5.7492565D-6,6.113243D-7,-5.89910D-8
!14 A,5.2070D-9,-4.233D-10,3.19D-11,-2.2D-12,1.0D-13/

  DATA ncfc, ncfd/18, 17/, xup/6.25D0/, sqrtpi/1.7724538509055160D0/ &
    , c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10), c(11) &
    , c(12), c(13), c(14), c(15), c(16), c(17), c(18) &
    /1.9449071068178803D0, 4.20186582324414D-2, -1.86866103976769D-2 &
    , 5.1281061839107D-3, -1.0683107461726D-3, 1.744737872522D-4 &
    , -2.15642065714D-5, 1.7282657974D-6, -2.00479241D-8 &
    , -1.64782105D-8, 2.0008475D-9, 2.57716D-11, -3.06343D-11 &
    , 1.9158D-12, 3.703D-13, -5.43D-14, -4.0D-15, 1.2D-15/ &
    , d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), d(9), d(10), d(11) &
    , d(12), d(13), d(14), d(15), d(16), d(17) &
    /1.4831105640848036D0, -3.010710733865950D-1, 6.89948306898316D-2 &
    , -1.39162712647222D-2, 2.4207995224335D-3, -3.658639685849D-4 &
    , 4.86209844323D-5, -5.7492565580D-6, 6.113243578D-7 &
    , -5.89910153D-8, 5.2070091D-9, -4.232976D-10, 3.18811D-11 &
    , -2.2361D-12, 1.467D-13, -9.0D-15, 5.0D-16/

  DATA zero, one, two, three, twenty, &
    half/0.0D0, 1.0D0, 2.0D0, 3.0D0, 20.0D0, 0.5D0/
! .. Executable Statements ..

! NO FAILURE EXITS
  ifail = 0
  xv = ABS(x)
  IF (xv >= xup) GO TO 120
  IF (xv <= two) GO TO 60
  x2 = two - twenty/(xv + three)

! SUMMATION
  bjp2 = zero
  bjp1 = c(ncfc)
  j = ncfc - 1
20 bj = x2*bjp1 - bjp2 + c(j)
  IF (j == 1) GO TO 40
  bjp2 = bjp1
  bjp1 = bj
  j = j - 1
  GO TO 20
40 x2 = half*(bj - bjp2)/xv*EXP(-x*x)/sqrtpi
  erf = (one - x2)*SIGN(one, x)
  GO TO 140

60 x2 = x*x - two
! SUMMATION
  bjp2 = zero
  bjp1 = d(ncfd)
  j = ncfd - 1
80 bj = x2*bjp1 - bjp2 + d(j)
  IF (j == 1) GO TO 100
  bjp2 = bjp1
  bjp1 = bj
  j = j - 1
  GO TO 80
100 erf = half*(bj - bjp2)*x
  GO TO 140

120 erf = SIGN(one, x)
140 RETURN
END FUNCTION erf
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

