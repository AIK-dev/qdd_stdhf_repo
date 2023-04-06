
!------------------------------------------------------------

SUBROUTINE addfunctofield(field, func, x0, y0, z0, fact)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! adds function func(x-x0,y-y0,z-z0) to field(x,y,z)
! version with zero parameters in func

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: x0
  REAL(DP), INTENT(IN) :: y0
  REAL(DP), INTENT(IN) :: z0
  REAL(DP), INTENT(IN) :: fact
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r)
      USE params
      REAL(DP), INTENT(IN) :: r
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz, max
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

  ind = 0

  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = z1 - z0
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = y1 - y0
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = x1 - x0
        rr = SQRT(rx*rx + ry*ry + rz*rz)

        rr = MAX(rr, small)

        ind = ind + 1

        field(ind) = field(ind) + func(rr)*fact

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE addfunctofield
!------------------------------------------------------------

SUBROUTINE addfunctofield1(field, func, x0, y0, z0, fact, para)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! adds function func(x-x0,y-y0,z-z0) to field(x,y,z)
! version with 1 parameter in func

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: x0
  REAL(DP), INTENT(IN) :: y0
  REAL(DP), INTENT(IN) :: z0
  REAL(DP), INTENT(IN) :: fact
  REAL(DP), INTENT(IN) :: para
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, s)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: s
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) ::rr, rx, ry, rz, x1, y1, z1

  ind = 0

  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = z1 - z0
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = y1 - y0
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = x1 - x0
        rr = SQRT(rx*rx + ry*ry + rz*rz)

        rr = MAX(rr, small)

        ind = ind + 1

        field(ind) = field(ind) + func(rr, para)*fact

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE addfunctofield1

SUBROUTINE addfunctofield3(field, func, x0, y0, z0, fact, para1, para2, para3)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! adds function func(x-x0,y-y0,z-z0) to field(x,y,z)
! version with 3 parameter in func

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: x0
  REAL(DP), INTENT(IN) :: y0
  REAL(DP), INTENT(IN) :: z0
  REAL(DP), INTENT(IN) :: fact
  REAL(DP), INTENT(IN) :: para1
  REAL(DP), INTENT(IN) :: para2
  REAL(DP), INTENT(IN) :: para3
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, a, b, c)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: a
      REAL(DP), INTENT(IN) :: b
      REAL(DP), INTENT(IN) :: c
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1
  ind = 0

  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = z1 - z0
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = y1 - y0
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = x1 - x0
        rr = SQRT(rx*rx + ry*ry + rz*rz)

        rr = MAX(rr, small)

        ind = ind + 1

        field(ind) = field(ind) + func(rr, para1, para2, para3)*fact

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE addfunctofield3

!------------------------------------------------------------

SUBROUTINE addfunctofieldonsubgrid(field, func, x0, y0, z0, fact, nsgsize)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! same as above routine but with subgrids

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: x0
  REAL(DP), INTENT(IN) :: y0
  REAL(DP), INTENT(IN) :: z0
  REAL(DP), INTENT(IN) :: fact
  INTEGER, INTENT(IN) :: nsgsize
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r)
      USE params
      REAL(DP), INTENT(IN) :: r
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

  INTEGER, EXTERNAL :: getnearestgridpoint
  INTEGER, EXTERNAL :: conv3to1

  ind = getnearestgridpoint(x0, y0, z0)

  CALL conv1to3(ind)

  DO iz = iindtmp(3) - nsgsize, iindtmp(3) + nsgsize
    IF (iz >= minz .AND. iz <= maxz) THEN
      z1 = (iz - nzsh)*dz
      rz = z1 - z0

      DO iy = iindtmp(2) - nsgsize, iindtmp(2) + nsgsize
        IF (iy >= miny .AND. iy <= maxy) THEN
          y1 = (iy - nysh)*dy
          ry = y1 - y0

          DO ix = iindtmp(1) - nsgsize, iindtmp(1) + nsgsize

            IF (ix >= minx .AND. ix <= maxx) THEN

              x1 = (ix - nxsh)*dx
              rx = x1 - x0

              rr = SQRT(rx*rx + ry*ry + rz*rz)

              rr = MAX(rr, small) ! avoid zero

              ind = conv3to1(ix, iy, iz)

              field(ind) = field(ind) + func(rr)*fact

            END IF
          END DO
        END IF
      END DO
    END IF

  END DO

  RETURN
END SUBROUTINE addfunctofieldonsubgrid

!------------------------------------------------------------

SUBROUTINE addfunctofieldonsubgrid1(field, func, x0, y0, z0, fact, para, nsgsize)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! same as above routine but with subgrids

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: x0
  REAL(DP), INTENT(IN) :: y0
  REAL(DP), INTENT(IN) :: z0
  REAL(DP), INTENT(IN) :: fact
  REAL(DP), INTENT(IN) :: para
  INTEGER, INTENT(IN) :: nsgsize
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, s)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: s
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

  INTEGER, EXTERNAL :: getnearestgridpoint
  INTEGER, EXTERNAL :: conv3to1

  ind = getnearestgridpoint(x0, y0, z0)

  CALL conv1to3(ind)

  DO iz = iindtmp(3) - nsgsize, iindtmp(3) + nsgsize
    IF (iz >= minz .AND. iz <= maxz) THEN
      z1 = (iz - nzsh)*dz
      rz = z1 - z0

      DO iy = iindtmp(2) - nsgsize, iindtmp(2) + nsgsize
        IF (iy >= miny .AND. iy <= maxy) THEN
          y1 = (iy - nysh)*dy
          ry = y1 - y0

          DO ix = iindtmp(1) - nsgsize, iindtmp(1) + nsgsize

            IF (ix >= minx .AND. ix <= maxx) THEN

              x1 = (ix - nxsh)*dx
              rx = x1 - x0

              rr = SQRT(rx*rx + ry*ry + rz*rz)

              rr = MAX(rr, small) ! avoid zero

              ind = conv3to1(ix, iy, iz)

              field(ind) = field(ind) + func(rr, para)*fact

            END IF
          END DO
        END IF
      END DO
    END IF

  END DO

  RETURN
END SUBROUTINE addfunctofieldonsubgrid1

!------------------------------------------------------------

SUBROUTINE addfunctofieldonsubgrid3(field, func, x0, y0, z0, fact, &
                                    para1, para2, para3, nsgsize)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! same as above routine but with subgrids

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: x0
  REAL(DP), INTENT(IN) :: y0
  REAL(DP), INTENT(IN) :: z0
  REAL(DP), INTENT(IN) :: fact
  REAL(DP), INTENT(IN) :: para1
  REAL(DP), INTENT(IN) :: para2
  REAL(DP), INTENT(IN) :: para3
  INTEGER, INTENT(IN) :: nsgsize
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, a, b, c)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: a
      REAL(DP), INTENT(IN) :: b
      REAL(DP), INTENT(IN) :: c
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

  INTEGER, EXTERNAL :: getnearestgridpoint
  INTEGER, EXTERNAL :: conv3to1

  ind = getnearestgridpoint(x0, y0, z0)

  CALL conv1to3(ind)

  DO iz = iindtmp(3) - nsgsize, iindtmp(3) + nsgsize
    IF (iz >= minz .AND. iz <= maxz) THEN
      z1 = (iz - nzsh)*dz
      rz = z1 - z0

      DO iy = iindtmp(2) - nsgsize, iindtmp(2) + nsgsize
        IF (iy >= miny .AND. iy <= maxy) THEN
          y1 = (iy - nysh)*dy
          ry = y1 - y0

          DO ix = iindtmp(1) - nsgsize, iindtmp(1) + nsgsize

            IF (ix >= minx .AND. ix <= maxx) THEN

              x1 = (ix - nxsh)*dx
              rx = x1 - x0

              rr = SQRT(rx*rx + ry*ry + rz*rz)

              rr = MAX(rr, small) ! avoid zero

              ind = conv3to1(ix, iy, iz)

              field(ind) = field(ind) + func(rr, para1, para2, para3)*fact

            END IF
          END DO
        END IF
      END DO
    END IF

  END DO

  RETURN
END SUBROUTINE addfunctofieldonsubgrid3
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE dintfieldfunc(field, func, xx, yy, zz, param)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! folds function around xx,yy,zz

! calculates the integral

! \int d^3r field(\vec r) \times \nabla_R func(|\vec r-\vec R|)

! where \vec R = (xx,yy,zz)

  REAL(DP), INTENT(IN) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz
  REAL(DP), INTENT(IN) :: param
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, s)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: s
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: sumx, sumy, sumz, radial, rr, rx, ry, rz, x1, y1, z1

  ind = 0
  sumx = 0D0
  sumy = 0D0
  sumz = 0D0

  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = zz - z1
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = yy - y1
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = xx - x1
        rr = SQRT(rx*rx + ry*ry + rz*rz)
        rr = MAX(rr, small) ! avoid zero
        ind = ind + 1

        radial = func(rr, param)/rr

        sumx = sumx + field(ind)*radial*rx
        sumy = sumy + field(ind)*radial*ry
        sumz = sumz + field(ind)*radial*rz

      END DO
    END DO
  END DO

  rvectmp(1) = sumx*dvol
  rvectmp(2) = sumy*dvol
  rvectmp(3) = sumz*dvol

  RETURN
END SUBROUTINE dintfieldfunc

!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE dintfieldfunconsubgrid(field, func, xx, yy, zz, param)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! folds function around xx,yy,zz using a subgrid ONLY

! calculates the integral

! \int d^3r field(\vec r) \times \nable_R func(|\vec r-\vec R|)

! where \vec R = (xx,yy,zz)

  REAL(DP), INTENT(IN) :: field(kdfull2)
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, s)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: s
    END FUNCTION func
  END INTERFACE

  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz
  REAL(DP), INTENT(IN) :: param

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: sumx, sumy, sumz, radial, rr, rx, ry, rz, x1, y1, z1

  INTEGER, EXTERNAL :: getnearestgridpoint
  INTEGER, EXTERNAL :: conv3to1

  ind = getnearestgridpoint(xx, yy, zz)

  CALL conv1to3(ind)

  sumx = 0D0
  sumy = 0D0
  sumz = 0D0

  DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
    z1 = (iz - nzsh)*dz
    rz = zz - z1

    DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
      y1 = (iy - nysh)*dy
      ry = yy - y1

      DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg
        x1 = (ix - nxsh)*dx
        rx = xx - x1

        rr = SQRT(rx*rx + ry*ry + rz*rz)
        rr = MAX(rr, small)

        ind = conv3to1(ix, iy, iz)

        radial = func(rr, param)/rr

        sumx = sumx + field(ind)*radial*rx
        sumy = sumy + field(ind)*radial*ry
        sumz = sumz + field(ind)*radial*rz

      END DO
    END DO
  END DO
  rvectmp(1) = sumx*dvol
  rvectmp(2) = sumy*dvol
  rvectmp(3) = sumz*dvol

  RETURN
END SUBROUTINE dintfieldfunconsubgrid

!GB

!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE updatesubgrids
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! if a particle has moved so far that we have to move the
! subgrid as well, then do so!


END SUBROUTINE updatesubgrids
!------------------------------------------------------------

!------------------------------------------------------------

INTEGER FUNCTION conv3to1(k, j, i)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN):: i, j, k

  conv3to1 = (i - 1)*2*nx*2*ny + (j - 1)*2*nx + k

  RETURN
END FUNCTION conv3to1
!------------------------------------------------------------

!------------------------------------------------------------

INTEGER FUNCTION getnearestgridpoint(xpos, ypos, zpos)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: xpos
  REAL(DP), INTENT(IN) :: ypos
  REAL(DP), INTENT(IN) :: zpos
!------------------------------------------------------------
! getNearestGridPoint

  INTEGER :: ind, ix0, iy0, iz0
  REAL(DP) :: x0, y0, z0


  x0 = xpos/dx
  ix0 = NINT(x0)

  y0 = ypos/dy
  iy0 = NINT(y0)

  z0 = zpos/dz
  iz0 = NINT(z0)

  ix0 = ix0 + nx
  iy0 = iy0 + ny
  iz0 = iz0 + nz

  IF (ix0 > maxx) ix0 = maxx ! avoid grid point OUT of box
  IF (ix0 < minx) ix0 = minx
  IF (iy0 > maxy) iy0 = maxy
  IF (iy0 < miny) iy0 = miny
  IF (iz0 > maxz) iz0 = maxz
  IF (iz0 < minz) iz0 = minz

! convert triple indices to single index:

  ind = (iz0 - 1)*2*nx*2*ny + (iy0 - 1)*2*nx + ix0

  getnearestgridpoint = ind

  IF (ind > kdfull2 .OR. ind <= 0) STOP 'invalid grid point in getnearestgridpoint'

  RETURN
END FUNCTION getnearestgridpoint

!------------------------------------------------------------
INTEGER FUNCTION getnearestgridpoint2(xpos, ypos, zpos)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

!------------------------------------------------------------
! getNearestGridPoint

  REAL(DP), INTENT(IN) :: xpos
  REAL(DP), INTENT(IN) :: ypos
  REAL(DP), INTENT(IN) :: zpos

  INTEGER :: ind, ix0, iy0, iz0
  REAL(DP) :: x0, y0, z0

  x0 = xpos/dx
  ix0 = INT(x0 + nxsh)

  y0 = ypos/dy
  iy0 = INT(y0 + nysh)

  z0 = zpos/dz
  iz0 = INT(z0 + nzsh)

  IF (ix0 > maxx) ix0 = maxx ! avoid grid point out of box
  IF (ix0 < minx) ix0 = minx
  IF (iy0 > maxy) iy0 = maxy
  IF (iy0 < miny) iy0 = miny
  IF (iz0 > maxz) iz0 = maxz
  IF (iz0 < minz) iz0 = minz

  WRITE (*, *) ' ix0,iy0,iz0=', ix0, iy0, iz0

! convert triple indices to single index:

  ind = (iz0 - 1)*nx2*ny2 + (iy0 - 1)*nx2 + ix0

  getnearestgridpoint2 = ind

  IF (ind > kdfull2 .OR. ind <= 0) STOP 'invalid grid point in getnearestgridpoint2'

  RETURN
END FUNCTION getnearestgridpoint2
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE conv1to3(ind)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! single index to triple index:

  INTEGER, INTENT(IN) :: ind

  INTEGER :: ind1, indx, indy, indz, iupper
  REAL(DP) :: rupper

  ind1 = ind

  rupper = REAL(ind1, DP)/(2*nx*2*ny)

  IF (rupper == INT(rupper)) THEN
    iupper = INT(rupper)
  ELSE
    iupper = INT(rupper) + 1
  END IF

  indz = iupper

  ind1 = ind1 - (indz - 1)*2*nx*2*ny

  rupper = REAL(ind1, DP)/(2*nx)
  IF (rupper == INT(rupper)) THEN
    iupper = INT(rupper)
  ELSE
    iupper = INT(rupper) + 1
  END IF

  indy = iupper
  indx = ind1 - (indy - 1)*2*nx

  iindtmp(1) = indx
  iindtmp(2) = indy
  iindtmp(3) = indz

  RETURN
END SUBROUTINE conv1to3

!----------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getxval(ind)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind

  CALL conv1to3(ind)

  getxval = (iindtmp(1) - nx)*dx
  RETURN
END FUNCTION getxval
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getyval(ind)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind

  CALL conv1to3(ind)

  getyval = (iindtmp(2) - ny)*dy
  RETURN
END FUNCTION getyval
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getzval(ind)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind

  CALL conv1to3(ind)

  getzval = (iindtmp(3) - nz)*dz
  RETURN
END FUNCTION getzval
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE renormsg(field, isgc, fromv, tov)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  INTEGER, INTENT(IN) :: isgc
  REAL(DP), INTENT(IN) :: fromv
  REAL(DP), INTENT(IN) :: tov

  INTEGER :: ind, ix, iy, iz

  INTEGER, EXTERNAL :: conv3to1

  ind = isgc

  CALL conv1to3(ind)

  DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
    DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
      DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg

        ind = conv3to1(ix, iy, iz)

        field(ind) = field(ind)*tov/fromv

      END DO
    END DO
  END DO

END SUBROUTINE renormsg
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION foldfunc(field, func, xx, yy, zz, param)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! folds function around xx,yy,zz

! calculates the integral

! \int d^3r field(\vec r) \times func(|\vec r-\vec R|)

! where \vec R = (xx,yy,zz)

  REAL(DP), INTENT(IN) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz
  REAL(DP), INTENT(IN) :: param
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, s)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: s
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: acc, rr, rx, ry, rz, x1, y1, z1

  ind = 0
  acc = 0D0

  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = z1 - zz
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = y1 - yy
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = x1 - xx
        rr = SQRT(rx*rx + ry*ry + rz*rz)
        rr = MAX(rr, small) ! avoid zero
        ind = ind + 1

        acc = acc + field(ind)*func(rr, param)

      END DO
    END DO
  END DO

  foldfunc = acc*dvol

  RETURN
END FUNCTION foldfunc
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE foldgradfunc(field, func, xx, yy, zz, param)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! folds function around xx,yy,zz

! calculates the integral

! \int d^3r field(\vec r) \times \nabla_R func(|\vec r-\vec R|)

! where \vec R = (xx,yy,zz)

  REAL(DP), INTENT(IN) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz
  REAL(DP), INTENT(IN) :: param
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, s)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: s
    END FUNCTION func
  END INTERFACE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: rder, sumx, sumy, sumz
  REAL(DP):: radial, rr, rx, ry, rz, x1, y1, z1

  ind = 0
  sumx = 0D0
  sumy = 0D0
  sumz = 0D0
  rder = 1.0D-5

  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    rz = zz - z1
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      ry = yy - y1
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        rx = xx - x1
        rr = SQRT(rx*rx + ry*ry + rz*rz)
        rr = MAX(rr, small) ! avoid zero
        ind = ind + 1

        radial = func(rr + rder, param) - func(rr - rder, param)

        sumx = sumx + field(ind)*radial*rx/rr
        sumy = sumy + field(ind)*radial*ry/rr
        sumz = sumz + field(ind)*radial*rz/rr

      END DO
    END DO
  END DO

  rvectmp(1) = sumx*dvol/(rder + rder)
  rvectmp(2) = sumy*dvol/(rder + rder)
  rvectmp(3) = sumz*dvol/(rder + rder)

  RETURN
END SUBROUTINE foldgradfunc
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE foldgradfunconsubgrid(field, func, xx, yy, zz, param)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! folds function around xx,yy,zz using a subgrid ONLY

! calculates the integral

! \int d^3r field(\vec r) \times \nable_R func(|\vec r-\vec R|)

! where \vec R = (xx,yy,zz)

  REAL(DP), INTENT(IN) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz
  REAL(DP), INTENT(IN) :: param
! function name passed as an argument :
  INTERFACE
    REAL(DP) FUNCTION func(r, s)
      USE params
      REAL(DP), INTENT(IN) :: r
      REAL(DP), INTENT(IN) :: s
    END FUNCTION func
  END INTERFACE

  INTEGER, EXTERNAL :: getnearestgridpoint
  INTEGER, EXTERNAL :: conv3to1

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: radial, rr, rder, rx, ry, rz, sumx, sumy, sumz, x1, y1, z1

  ind = getnearestgridpoint(xx, yy, zz)

  CALL conv1to3(ind)

  sumx = 0.0D0
  sumy = 0.0D0
  sumz = 0.0D0
  rder = 1.0D-5

  DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
    IF (iz >= minz .AND. iz <= maxz) THEN
      z1 = (iz - nzsh)*dz
      rz = zz - z1

      DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
        IF (iy >= miny .AND. iy <= maxy) THEN
          y1 = (iy - nysh)*dy
          ry = yy - y1

          DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg
            IF (ix >= minx .AND. ix <= maxx) THEN
              x1 = (ix - nxsh)*dx
              rx = xx - x1

              rr = SQRT(rx*rx + ry*ry + rz*rz)
              rr = MAX(rr, small)

              ind = conv3to1(ix, iy, iz)

              radial = func(rr + rder, param) - func(rr - rder, param)

              sumx = sumx + field(ind)*radial*rx/rr
              sumy = sumy + field(ind)*radial*ry/rr
              sumz = sumz + field(ind)*radial*rz/rr

            END IF
          END DO
        END IF
      END DO
    END IF
  END DO

  rvectmp(1) = sumx*dvol/(rder + rder)
  rvectmp(2) = sumy*dvol/(rder + rder)
  rvectmp(3) = sumz*dvol/(rder + rder)

  RETURN
END SUBROUTINE foldgradfunconsubgrid
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION dintprod(field1, field2)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: field1(kdfull2)
  REAL(DP), INTENT(IN) :: field2(kdfull2)

  INTEGER :: ind, ix, iy, iz
  REAL(DP):: acc

  ind = 0
  acc = 0.0D0
  DO iz = minz, maxz
    DO iy = miny, maxy
      DO ix = minx, maxx
        ind = ind + 1
        acc = acc + field1(ind)*field2(ind)
      END DO
    END DO
  END DO

  dintprod = acc*dvol

  RETURN
END FUNCTION dintprod
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION dintprodsubgrid(field1, field2, xx, yy, zz)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! calculates the integral of the product of two scalar
! fields not over the whole box, but ONLY over a subgrid
! which is centered at xx,yy,zz

  REAL(DP), INTENT(IN) :: field1(kdfull2)
  REAL(DP), INTENT(IN) :: field2(kdfull2)
  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: acc

  INTEGER, EXTERNAL :: getnearestgridpoint
  INTEGER, EXTERNAL :: conv3to1

  ind = getnearestgridpoint(xx, yy, zz)

  CALL conv1to3(ind)

  acc = 0.0D0
  DO iz = iindtmp(3) - nzsg, iindtmp(3) + nzsg
    DO iy = iindtmp(2) - nysg, iindtmp(2) + nysg
      DO ix = iindtmp(1) - nxsg, iindtmp(1) + nxsg

        ind = conv3to1(ix, iy, iz)

        acc = acc + field1(ind)*field2(ind)

      END DO
    END DO
  END DO

  dintprodsubgrid = acc*dvol

  RETURN
END FUNCTION dintprodsubgrid
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION dintfieldtimesgauss(field, xx, yy, zz, sigm, ccharg)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! calculates the integral:

! \int d^3r f(\vec{r})\cdot A exp(-|\vec{r}-\vec{R}|/(2\sigma^2))
! with A being the proper factor normalising the Gaussian to ccharg

! ATTENTION: be aware that the factor e2 is NOT included here,
! but must be attached in the calling routine

! xx,yy,zz are the center of the Gaussian, sigm its width
! according to the convention in this code:
! exp(-r**2/(2.*sigm**2)) !!!

  REAL(DP), INTENT(IN) :: field(kdfull2)
  REAL(DP), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: yy
  REAL(DP), INTENT(IN) :: zz
  REAL(DP), INTENT(IN) :: sigm
  REAL(DP), INTENT(IN) :: ccharg

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: a, acc, sigm2
  REAL(DP) :: r2, rx, ry, rz, x1, y1, z1

  INTEGER, EXTERNAL :: getnearestgridpoint
  INTEGER, EXTERNAL :: conv3to1

  a = ccharg/(pi**1.5D0*2D0**1.5D0*sigm**3)
  sigm2 = sigm*sigm
  acc = 0D0


    ind = 0

    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      rz = z1 - zz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        ry = y1 - yy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          rx = x1 - xx

          r2 = rx*rx + ry*ry + rz*rz

          acc = acc + field(ind)*EXP(-r2/2D0/sigm2)

        END DO
      END DO
    END DO


  acc = acc*a*dvol

  dintfieldtimesgauss = acc

  RETURN
END FUNCTION dintfieldtimesgauss
!------------------------------------------------------------

!------------------------------------------------------------


