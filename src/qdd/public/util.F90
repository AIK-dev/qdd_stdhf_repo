MODULE util

  IMPLICIT NONE

  INTERFACE wfovlp
    MODULE PROCEDURE r_wfovlp, c_wfovlp
  END INTERFACE wfovlp

  INTERFACE wfnorm
    MODULE PROCEDURE r_wfnorm, c_wfnorm
  END INTERFACE wfnorm

  INTERFACE project
    MODULE PROCEDURE r_project, c_project
  END INTERFACE project

CONTAINS

!------------------------------------------------------------
  SUBROUTINE shiftfield(q0, shix, shiy, shiz)
!------------------------------------------------------------
    USE params
    IMPLICIT NONE

! Shifts array 'q0' by shix,shiy,shiz along grid points. Shift values
! shix,shiy,shiz must be positive
!
! Input/Output:
! q0 = REAL array for one spatial wavefunction to be shifted
! Input:
! shix = shift distance along x direction, nearest multiple of
! grid spacing is actually taken; must be > 0.
! shiy = as shix, but for y-direction
! shiz = as shix, but for y-direction

    REAL(DP), INTENT(IN OUT) :: q0(kdfull2)
    REAL(DP), INTENT(IN) :: shix
    REAL(DP), INTENT(IN) :: shiy
    REAL(DP), INTENT(IN) :: shiz

    REAL(DP), ALLOCATABLE :: q1(:)

    INTEGER :: ind, ind1, ix, iy, iz, ix1, iy1, iz1, nshift
    INTEGER, EXTERNAL :: conv3to1

    ALLOCATE (q1(kdfull2))

    q1 = q0

    nshift = nint(shix/dx)
    IF (nshift < 0) STOP ' SHIX must be >= 0'

    IF (nshift > 0) THEN
      DO ix = minx + nshift, maxx
        DO iy = miny, maxy
          DO iz = minz, maxz

            ix1 = ix - nshift

            ind = conv3to1(ix, iy, iz)
            ind1 = conv3to1(ix1, iy, iz)

            q0(ind) = q1(ind1)

          END DO
        END DO
      END DO

      DO ix = minx, minx + nshift - 1
        DO iy = miny, maxy
          DO iz = minz, maxz

            ind = conv3to1(ix, iy, iz)

            q0(ind) = 0D0

          END DO
        END DO
      END DO
    END IF

!::::::::::::::::::::::::::::::

    q1 = q0

    nshift = nint(shiy/dy)
    IF (nshift < 0) STOP ' SHIY must be >= 0'

    IF (nshift > 0) THEN
      DO iy = miny + nshift, maxy
        DO ix = minx, maxx
          DO iz = minz, maxz

            iy1 = iy - nshift

            ind = conv3to1(ix, iy, iz)
            ind1 = conv3to1(ix, iy1, iz)

            q0(ind) = q1(ind1)

          END DO
        END DO
      END DO

      DO iy = miny, miny + nshift - 1
        DO ix = minx, maxx
          DO iz = minz, maxz

            ind = conv3to1(ix, iy, iz)

            q0(ind) = 0D0

          END DO
        END DO
      END DO
    END IF

!::::::::::::::::::::::::::::::

    q1 = q0

    nshift = nint(shiz/dz)
    IF (nshift < 0) STOP ' SHIZ must be >= 0'

    IF (nshift > 0) THEN
      DO iz = minz + nshift, maxz
        DO ix = minx, maxx
          DO iy = miny, maxy

            iz1 = iz - nshift

            ind = conv3to1(ix, iy, iz)
            ind1 = conv3to1(ix, iy, iz1)

            q0(ind) = q1(ind1)

          END DO
        END DO
      END DO

      DO iz = minz, minz + nshift - 1
        DO ix = minx, maxx
          DO iy = miny, maxy

            ind = conv3to1(ix, iy, iz)

            q0(ind) = 0D0

          END DO
        END DO
      END DO
    END IF

    DEALLOCATE (q1)

    RETURN
  END SUBROUTINE shiftfield
!------------------------------------------------------------

!------------------------------------------------------------
  CHARACTER*9 FUNCTION inttostring(inumber)
!------------------------------------------------------------

! translates INTEGER number to corresponding CHARACTER variable

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: inumber

    SELECT CASE (inumber)
    CASE (0:999999999)
      WRITE (inttostring, '(i9)') inumber
    CASE DEFAULT
      STOP 'ERROR in intToString : input argument outside range 0 to 999999999'
    END SELECT

    RETURN
  END FUNCTION inttostring
!------------------------------------------------------------

!------------------------------------------------------------
  SUBROUTINE pm3dcut(iunit, iax1, iax2, value, field)
!------------------------------------------------------------

! Prepares 3D plotting with 'pm3d' in 'gnuplot.
!
! Input:
! iunit = write unit number-
! iax1,iax2 = choices of the two axes (from 1..3) along which is plotted.
! value = coordinate on third axis (complement of iax1,iax2) along
! which field is plotted.
! field = 3D array from which 2D cut is retrieved.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    INTEGER, INTENT(IN) :: iax1
    INTEGER, INTENT(IN) :: iax2
    REAL(DP), INTENT(IN) :: value
    REAL(DP), INTENT(IN) :: field(kdfull2)

    INTEGER :: ind, ix, iy, iz, x1, y1, z1

    IF (iax1 == 1 .AND. iax2 == 2) THEN ! xy-plane

      ind = 0
      DO iz = 1, maxz
        z1 = (iz - nzsh)*dz

        DO iy = 1, maxy
          y1 = (iy - nysh)*dy
          DO ix = 1, maxx
            x1 = (ix - nxsh)*dx

            ind = ind + 1

            IF (z1 == value) THEN
              WRITE (iunit, '(3f10.3,1e15.6)') x1, y1, z1, field(ind)
            END IF
          END DO
          IF (z1 == value) THEN
            WRITE (iunit, *)
          END IF
        END DO

      END DO

    ELSE IF (iax1 == 1 .AND. iax2 == 3) THEN ! xz-plane

      ind = 0
      DO iz = 1, maxz
        z1 = (iz - nzsh)*dz

        DO iy = 1, maxy
          y1 = (iy - nysh)*dy

          DO ix = 1, maxx
            x1 = (ix - nxsh)*dx

            ind = ind + 1
            IF (y1 == value) THEN
              WRITE (iunit, '(3f10.3,1e15.6)') x1, y1, z1, field(ind)
            END IF

          END DO
          IF (y1 == value) THEN
            WRITE (iunit, *)
          END IF

        END DO

      END DO

    ELSE IF (iax1 == 2 .AND. iax2 == 3) THEN ! yz-plane

      ind = 0
      DO iz = 1, maxz
        z1 = (iz - nzsh)*dz

        DO iy = 1, maxy
          y1 = (iy - nysh)*dy

          DO ix = 1, maxx
            x1 = (ix - nxsh)*dx

            ind = ind + 1
            IF (x1 == value) THEN
              WRITE (iunit, '(3f10.3,1e15.6)') x1, y1, z1, field(ind)
            END IF

          END DO

        END DO
        WRITE (iunit, *)
      END DO

    END IF

    RETURN
  END SUBROUTINE pm3dcut
!------------------------------------------------------------

!------------------------------------------------------------
  COMPLEX(DP) FUNCTION orbitaloverlap(q1, q2)
!------------------------------------------------------------
    USE params

    IMPLICIT NONE

! Overlap <q1|q2> = \int\d^3r q_1^* q_2 between the two COMPLEX
! wavefunctions q1 and q2.

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(IN) :: q2(kdfull2)

    INTEGER :: ind
    REAL(DP) :: sumr, sumi
    COMPLEX(DP) :: csum

    sumr = 0D0
    sumi = 0D0

    DO ind = 1, kdfull2
      sumr = sumr + REAL(q1(ind), DP)*REAL(q2(ind), DP) + AIMAG(q1(ind))*AIMAG(q2(ind))
      sumi = sumi + REAL(q1(ind), DP)*AIMAG(q2(ind)) - REAL(q2(ind), DP)*AIMAG(q1(ind))
    END DO

    csum = CMPLX(sumr*dx*dy*dz, sumi*dx*dy*dz, DP)

    orbitaloverlap = csum

    RETURN
  END FUNCTION orbitaloverlap
!------------------------------------------------------------

!------------------------------------------------------------
  REAL(DP) FUNCTION realoverlap(q1, q2)
!------------------------------------------------------------
    USE params
    IMPLICIT NONE

! returns \int\d^3r q_1^* q_2

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(IN) :: q2(kdfull2)

    INTEGER :: ind
    REAL(DP) :: sumr

    sumr = 0D0

    DO ind = 1, kdfull2
      sumr = sumr + REAL(q1(ind), DP)*REAL(q2(ind), DP) + AIMAG(q1(ind))*AIMAG(q2(ind))
    END DO

    realoverlap = sumr*dvol

    RETURN
  END FUNCTION realoverlap
!------------------------------------------------------------

!-------------------------------------------------------------
  REAL(DP) FUNCTION realovsubgrid(q1, q2, ion)

! REAL part of overlap between wavefunctions 'q1' and 'q2'
! accumulated on the subgrid of ion 'ion' defined by
! the non-local PsP.

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(IN) :: q2(kdfull2)
    INTEGER, INTENT(IN) :: ion

    INTEGER :: i, ii
    REAL(DP) :: sumr
!-------------------------------------------------------------

    sumr = 0D0
    DO i = 1, ifin(ion)
      ii = icount(i, ion)
      sumr = REAL(q2(ii), DP)*REAL(q1(ii), DP) + AIMAG(q2(ii))*AIMAG(q1(ii)) + sumr
    END DO

    realovsubgrid = sumr*dvol

    RETURN
  END FUNCTION realovsubgrid
!------------------------------------------------------------

!------------------------------------------------------------
  SUBROUTINE gettemperature(iflag,rkt)
!------------------------------------------------------------

! Compute ionic kinetic energy and PRINT corresponding ionic temperature
! on unit 175.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iflag
    REAL(8), INTENT(OUT) :: rkt

    INTEGER :: i
    REAL(DP) :: ek, ekcm, pcmx, pcmy, pcmz, rm, rmcm, tt

    IF (iflag == 4) THEN

      pcmx = 0D0
      pcmy = 0D0
      pcmz = 0D0

      DO i = 1, nion
        pcmx = pcmx + cpx(i)
        pcmy = pcmy + cpy(i)
        pcmz = pcmz + cpz(i)
      END DO

      rmcm = nion*amu(np(1))*1836D0*0.5D0
      rm = amu(np(1))*1836D0*0.5D0
      ekcm = pcmx**2 + pcmy**2 + pcmz**2
      ekcm = ekcm/2D0/rmcm

      ek = 0D0
      DO i = 1, nion
        ek = ek + (cpx(i)**2 + cpy(i)**2 + cpz(i)**2)/2D0/rm
      END DO

      rkt = 2D0/(3D0*nion - 3D0)*(ek - ekcm) ! thermal energy kT in Ry

! Temperature in Kelvin



    END IF

    RETURN
  END SUBROUTINE gettemperature
!------------------------------------------------------------

  SUBROUTINE densbelow(field)
!------------------------------------------------------------

! Integrate 'field' ober x-y plane, accumulate along z axis and print
! on unit 650.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: field(kdfull2)

    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: acc, x1, y1, z1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vs(:)

    ALLOCATE (vs(maxz + 1))
    ind = 0

    acc = 0D0

    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          acc = acc + field(ind)
        END DO
      END DO
      vs(iz) = acc*dx*dy*dz
      WRITE (650, '(2f10.5,1e20.9)') tfs, z1, vs(iz)
    END DO

    DEALLOCATE (vs)

    RETURN
  END SUBROUTINE densbelow
!------------------------------------------------------------

!-----priCM-------------------------------------------------------
  SUBROUTINE pricm(rho)
!------------------------------------------------------------
    USE params
    IMPLICIT NONE

! calculates center of density rho
! returns result on rVecTmp(1:3) via module params

    REAL(DP), INTENT(IN) :: rho(kdfull2)

    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: acc, xcm, ycm, zcm, x1, y1, z1

!------------------------------------------------------------

    xcm = 0D0
    ycm = 0D0
    zcm = 0D0
    acc = 0D0

    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1

          xcm = xcm + rho(ind)*x1
          ycm = ycm + rho(ind)*y1
          zcm = zcm + rho(ind)*z1
          acc = acc + rho(ind)
        END DO
      END DO
    END DO

    rvectmp(1) = xcm/acc
    rvectmp(2) = ycm/acc
    rvectmp(3) = zcm/acc

  END SUBROUTINE pricm
!------------------------------------------------------------

!-----priCM_state-------------------------------------------------------

  SUBROUTINE pricm_state(psir)
!------------------------------------------------------------
    USE params
    IMPLICIT NONE

! calculates center of density for states in wavefunctions 'psir'
! and prints results on standard output

    REAL(DP), INTENT(IN) :: psir(kdfull2, kstate)

    INTEGER :: i1, i2, i3, ind, nbr
    REAL(DP) :: xx, yy, zz

!------------------------------------------------------------

    DO nbr = 1, nstate
      xx = 0D0
      yy = 0D0
      zz = 0D0
      ind = 0
      DO i3 = minz, maxz
        DO i2 = miny, maxy
          DO i1 = minx, maxx
            ind = ind + 1
            xx = xx + psir(ind, nbr)*psir(ind, nbr)*(i1 - nxsh)*dx
            yy = yy + psir(ind, nbr)*psir(ind, nbr)*(i2 - nysh)*dy
            zz = zz + psir(ind, nbr)*psir(ind, nbr)*(i3 - nzsh)*dz
          END DO
        END DO
      END DO
      xx = xx*dvol
      yy = yy*dvol
      zz = zz*dvol
      WRITE (6, '(a,i5,a,i5,a,3(1pg13.5))') 'node: myn=', &
        myn, ', state: n=', nbr, '<x>,<y>,<z>=', xx, yy, zz
    END DO

    RETURN
  END SUBROUTINE pricm_state

!------------------------------------------------------------

  SUBROUTINE rotatevec(x, y, z, iax, alpha)
!------------------------------------------------------------
! rotation of vector '(x,y,z)' about axis 'iax' with angle 'alpha'
! returning the result on vector 'rvectmp' via module 'params'..
    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN)::x
    REAL(DP), INTENT(IN)::y
    REAL(DP), INTENT(IN)::z
    INTEGER, INTENT(IN)::iax
    REAL(DP), INTENT(IN)::alpha

    IF (iax == 1) THEN ! rotate along x-axis
      rvectmp(1) = x
      rvectmp(2) = COS(alpha)*y - SIN(alpha)*z
      rvectmp(3) = SIN(alpha)*y + COS(alpha)*z
    ELSE IF (iax == 2) THEN
      rvectmp(1) = COS(alpha)*x + SIN(alpha)*z
      rvectmp(2) = y
      rvectmp(3) = -SIN(alpha)*x + COS(alpha)*z
    ELSE IF (iax == 3) THEN
      rvectmp(1) = COS(alpha)*x - SIN(alpha)*y
      rvectmp(2) = SIN(alpha)*x + COS(alpha)*y
      rvectmp(3) = z
    ELSE
      STOP 'Error in SUBROUTINE rotateVec'
    END IF

    RETURN
  END SUBROUTINE rotatevec
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE rotatevec3D(vecin, vecout, vecalpha)
!------------------------------------------------------------
! rotation of 3D vector 'vecin' by rotation vector 'vecalpha'
! returning result on 'vecout'.
! Angle to be given in radian.
!
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: vecin(3)
    REAL(DP), INTENT(IN) :: vecalpha(3)
    REAL(DP), INTENT(OUT) :: vecout(3)

    REAL(DP) :: absalpha, sa, c, c2a2

! prepare auxiliary variables

    absalpha = SQRT(SUM(vecalpha**2))
    c = COS(absalpha)
    sa = SIN(absalpha)/absalpha
    c2a2 = (c - 1D0)/(absalpha*absalpha)

! compute transformation

    vecout(1) = (c - c2a2*vecalpha(1)**2)*vecin(1) &
                - (sa*vecalpha(3) + c2a2*vecalpha(1)*vecalpha(2))*vecin(2) &
                + (sa*vecalpha(2) - c2a2*vecalpha(1)*vecalpha(3))*vecin(3)

    vecout(2) = (sa*vecalpha(3) - c2a2*vecalpha(1)*vecalpha(2))*vecin(1) &
                + (c - c2a2*vecalpha(2)**2)*vecin(2) &
                - (sa*vecalpha(1) + c2a2*vecalpha(2)*vecalpha(3))*vecin(3)

    vecout(3) = -(sa*vecalpha(2) + c2a2*vecalpha(1)*vecalpha(3))*vecin(1) &
                + (sa*vecalpha(1) - c2a2*vecalpha(2)*vecalpha(3))*vecin(2) &
                + (c - c2a2*vecalpha(3)**2)*vecin(3)

    RETURN
  END SUBROUTINE rotatevec3D
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE swaparr(arr1, arr2, n)
!------------------------------------------------------------
    USE params
    IMPLICIT NONE

! swaps arrays arr1 and arr2

    REAL(DP), INTENT(IN OUT) :: arr1(kdfull2)
    REAL(DP), INTENT(IN OUT) :: arr2(kdfull2)
    INTEGER, INTENT(IN) :: n

    INTEGER :: i
    REAL(DP), ALLOCATABLE :: arrt(:)

    ALLOCATE (arrt(kdfull2))
    DO i = 1, n
      arrt(i) = arr1(i)
      arr1(i) = arr2(i)
      arr2(i) = arrt(i)
    END DO
    DEALLOCATE (arrt)

    RETURN
  END SUBROUTINE swaparr
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE cparr(arr1, arr2, n)
!------------------------------------------------------------
    USE params
    IMPLICIT NONE

! copies array arr1 to arr2

    REAL(DP), INTENT(IN) :: arr1(kdfull2)
    REAL(DP), INTENT(OUT) :: arr2(kdfull2)
    INTEGER, INTENT(IN) :: n

    INTEGER :: i

    DO i = 1, n
      arr2(i) = arr1(i)
    END DO

    RETURN
  END SUBROUTINE cparr
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE makexbsfile(ionlymob, iwithsh)
!------------------------------------------------------------

! Prepares input for plotting ionic configuration with the chemical
! grafical freeware 'xbs'.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ionlymob
    INTEGER, INTENT(IN) :: iwithsh


    OPEN (789, STATUS='unknown', FILE='xbs.Na_MgO')


    WRITE (789, *) 'spec O 0.6 0.7 0.3 0.3'
    WRITE (789, *) 'spec Mg 0.6 1.0 0.0 0.0'
    WRITE (789, *) 'spec Osh 1.4 0.0 0.0 1.0'
    WRITE (789, *) 'tmat 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0'
    WRITE (789, *) 'dist 51.629'
    WRITE (789, *) 'inc 5.000'
    WRITE (789, *) 'scale 15.671'
    WRITE (789, *) 'rfac 1.00'
    WRITE (789, *) 'bfac 1.00'
    WRITE (789, *) 'pos 0.000 0.000'
    WRITE (789, *) 'switches 1 0 1 0 0 1 1 0 0'

    CLOSE (789)

    RETURN
  END SUBROUTINE makexbsfile
!------------------------------------------------------------

  SUBROUTINE printfieldx(iunit, field, y, z)
!------------------------------------------------------------

! PRINT 'field' along x-axis for grid points y and z.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    REAL(DP), INTENT(IN) :: field(kdfull2)
    REAL(DP), INTENT(IN) :: y
    REAL(DP), INTENT(IN) :: z

    INTEGER :: ind, ix
    REAL(DP) :: x
    INTEGER, EXTERNAL :: getnearestgridpoint

    DO ix = minx, maxx
      x = (ix - nxsh)*dx

      ind = getnearestgridpoint(x, y, z)

      WRITE (iunit, '(3f12.4,2e17.7)') x, y, z, field(ind)

    END DO

    RETURN
  END SUBROUTINE printfieldx
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE printfieldy(iunit, field, x, z)
!------------------------------------------------------------

! print 'field' along y-axis for grid points x and z.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    REAL(DP), INTENT(IN) :: field(kdfull2)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: z

    INTEGER :: ind, iy
    REAL(DP) :: y
    INTEGER, EXTERNAL :: getnearestgridpoint

    DO iy = miny, maxy
      y = (iy - nysh)*dy

      ind = getnearestgridpoint(x, y, z)

      WRITE (iunit, '(3f12.4,2e17.7)') x, y, z, field(ind)

    END DO

    RETURN
  END SUBROUTINE printfieldy
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE printfieldz(iunit, field, x, y)
!------------------------------------------------------------

! print 'field' along z-axis for grid points x and y.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    REAL(DP), INTENT(IN) :: field(kdfull2)
    REAL(DP), INTENT(IN) :: x
    REAL(DP), INTENT(IN) :: y

    INTEGER :: ind, iz
    REAL(DP) :: z
    INTEGER, EXTERNAL :: getnearestgridpoint

    DO iz = minz, maxz
      z = (iz - nzsh)*dz

      ind = getnearestgridpoint(x, y, z)

      WRITE (iunit, '(3f12.4,2e17.7)') x, y, z, field(ind)

    END DO

    RETURN
  END SUBROUTINE printfieldz
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE printfield(iunit, field, comment)
!------------------------------------------------------------

! print full 3D field with associated coordinates

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    REAL(DP), INTENT(IN) :: field(kdfull2)
    CHARACTER(LEN=*), INTENT(IN) :: comment

    LOGICAL :: topenunit
    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: x1, y1, z1

    INQUIRE (iunit, OPENED=topenunit)
    IF (.NOT. topenunit) THEN
      OPEN (iunit, FILE=comment//'.'//outnam)
      WRITE (iunit, '(a)') '# x y z '//comment
    END IF

    ind = 0
    DO iz = 1, maxz
      z1 = (iz - nzsh)*dz
      DO iy = 1, maxy
        y1 = (iy - nysh)*dy
        DO ix = 1, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          WRITE (iunit, '(3f12.4,1e17.7)') x1, y1, z1, field(ind)
        END DO
      END DO
    END DO

    IF (.NOT. topenunit) CLOSE (iunit)

    RETURN
  END SUBROUTINE printfield
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE printCubeFile(iunit, field)
!------------------------------------------------------------

! print a Gaussian-style Cube FILE with ion
! structure and electron density. Details on this FORMAT
! can be obtained: http://paulbourke.net/dataformats/cube/

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    REAL(DP), INTENT(IN) :: field(kdfull2)

    LOGICAL :: topenunit
    INTEGER :: index, ix, iy, iz, ion
    REAL(DP) :: Ox, Oy, Oz
    REAL(DP), ALLOCATABLE :: density(:, :, :)

    ALLOCATE (density(nx2, ny2, nz2))

    INQUIRE (iunit, OPENED=topenunit)
    IF (.NOT. topenunit) THEN
      OPEN (iunit, FILE=TRIM(outnam)//'.cube')
    END IF

    ! Building 3D array OUT of 1D array by using FUNCTION for index conversion
    DO index = 1, nxyz
      CALL conv1to3(index)
      density(iindtmp(1), iindtmp(2), iindtmp(3)) = field(index)
    END DO

    WRITE (iunit, '(A13)') "QDD CUBE FILE"
    WRITE (iunit, '(A44)') "Electron number-density"

    Ox = (1 - nx)*dx
    Oy = (1 - ny)*dy
    Oz = (1 - nz)*dz

    WRITE (iunit, '(I3,3(2X, F9.5))') nion, Ox, Oy, Oz
    WRITE (iunit, '(I3,3(2X, F9.7))') nx2, dx, 0D0, 0D0
    WRITE (iunit, '(I3,3(2X, F9.7))') ny2, 0D0, dy, 0D0
    WRITE (iunit, '(I3,3(2X, F9.7))') nz2, 0D0, 0D0, dz

    DO ion = 1, nion
      WRITE (iunit, '(I3, 2X, I1, 3(2X,F11.7))') np(ion), 0, cx(ion), cy(ion), cz(ion)
    END DO

    index = 1
    DO ix = 1, maxx
      DO iy = 1, maxy
        DO iz = 1, maxz
          WRITE (iunit, '(E13.7,2X)', ADVANCE='no') density(ix, iy, iz)
          IF (MOD(index, 6) == 0) THEN
            WRITE (iunit, *)
          END IF
          index = index + 1
        END DO
      END DO
    END DO

    IF (.NOT. topenunit) CLOSE (iunit)

    DEALLOCATE (density)
    RETURN
  END SUBROUTINE printCubeFile
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE printfield2(iunit, field1, field2)
!------------------------------------------------------------

! print two full 3D fields with associated coordinates

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    REAL(DP), INTENT(IN) :: field1(kdfull2)
    REAL(DP), INTENT(IN) :: field2(kdfull2)

    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: ft1, ft2, x1, y1, z1

    ind = 0
    DO iz = 1, maxz
      z1 = (iz - nzsh)*dz
      DO iy = 1, maxy
        y1 = (iy - nysh)*dy
        DO ix = 1, maxx
          x1 = (ix - nxsh)*dx

          ind = ind + 1

          ft1 = field1(ind)
          ft2 = field2(ind)

          IF (ABS(ft1) < 1D-90) ft1 = 0D0
          IF (ABS(ft2) < 1D-90) ft2 = 0D0

          WRITE (iunit, '(3f12.4,2e17.7)') x1, y1, z1, ft1, ft2

        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE printfield2
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE printforces(iflag, iterat)
!------------------------------------------------------------

! print forces on ions

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iflag
    INTEGER, INTENT(IN) :: iterat

    INTEGER :: i

    IF (iterat > 0) THEN

      WRITE (6, *) 'CHECKING FORCES !!', iflag

      IF (iflag == 4) THEN
        DO i = 1, nion
          WRITE (6, '(3e17.8)') fx(i), fy(i), fz(i)
        END DO
      END IF

    END IF

    RETURN
  END SUBROUTINE printforces

!------------------------------------------------------------

  SUBROUTINE getcm(iflag, iflagc, iflagk)
!------------------------------------------------------------

! Calculates the center of mass and stores it in the vector
! rVecTmp(1:3) which is communicated via module 'params'.

    USE params
    IMPLICIT NONE
    INTEGER, INTENT(IN)::iflag
    INTEGER, INTENT(IN)::iflagc
    INTEGER, INTENT(IN)::iflagk

    INTEGER :: i
    REAL(DP) :: summ, sumx, sumy, sumz

    IF (nion2 == 0) THEN
      rvectmp = 0D0
      RETURN
    END IF

    summ = 0D0
    sumx = 0D0
    sumy = 0D0
    sumz = 0D0

! the COMMON scaling factor 1836.0*ame can be skipped here!

    IF (iflag /= 0) THEN
      DO i = 1, nion
        summ = summ + amu(np(i))
        sumx = sumx + amu(np(i))*cx(i)
        sumy = sumy + amu(np(i))*cy(i)
        sumz = sumz + amu(np(i))*cz(i)
      END DO
    END IF


    rvectmp(1) = sumx/summ
    rvectmp(2) = sumy/summ
    rvectmp(3) = sumz/summ

    RETURN
  END SUBROUTINE getcm
!------------------------------------------------------------

!-----emoms-------------------------------------------------------emoms

! Various multipole moments for the given electronic distribution
! 'rho'. result is array of moments 'qe' communicated via module 'params'.

  SUBROUTINE emoms(rho)
    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rho(2*kdfull2)

    INTEGER :: ind, ix, iy, iz, k
    REAL(DP) :: updens, dodens
    REAL(DP) :: sux, sdx, suy, sdy, suz, sdz, s, s1, s2
    REAL(DP) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
!----------------------------------------------------------------------

    sux = 0D0
    sdx = 0D0
    suy = 0D0
    sdy = 0D0
    suz = 0D0
    sdz = 0D0
    s1 = 0D0
    s2 = 0D0

    nrmom = 35
    IF (nrmom > kmom) STOP ' too many moments in EMOMS'

    DO k = 1, nrmom
      qe(k) = 0D0
    END DO

! switch for calculating moments relative to center of mass (1)
! or center of box (0)

    rvectmp = 0D0
    IF (tmoms_rel_cm .AND. nion2 > 0) CALL getcm(1, 0, 0)

    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      z1 = z1 - rvectmp(3)

      z2 = z1*z1
      z3 = z2*z1
      z4 = z2*z2
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        y1 = y1 - rvectmp(2)

        y2 = y1*y1
        y3 = y2*y1
        y4 = y2*y2
        DO ix = minx, maxx
          ind = ind + 1
          IF ((ix <= nx2) .AND. (iy <= ny2) .AND. (iz <= nz2)) THEN
            x1 = (ix - nxsh)*dx
            x1 = x1 - rvectmp(1)

            x2 = x1*x1
            x3 = x2*x1
            x4 = x2*x2
            s = rho(ind)
! monopole
            qe(1) = qe(1) + s
! dipole
            qe(2) = qe(2) + s*x1
            qe(3) = qe(3) + s*y1
            qe(4) = qe(4) + s*z1
! quadrupole
            qe(5) = qe(5) + s*x2
            qe(6) = qe(6) + s*y2
            qe(7) = qe(7) + s*z2
            qe(8) = qe(8) + s*x1*y1
            qe(9) = qe(9) + s*z1*x1
            qe(10) = qe(10) + s*z1*y1
! octupole
            qe(11) = qe(11) + s*x3
            qe(12) = qe(12) + s*x2*y1
            qe(13) = qe(13) + s*z1*x2
            qe(14) = qe(14) + s*x1*y2
            qe(15) = qe(15) + s*x1*y1*z1
            qe(16) = qe(16) + s*z2*x1
            qe(17) = qe(17) + s*y3
            qe(18) = qe(18) + s*z1*y2
            qe(19) = qe(19) + s*z2*y1
            qe(20) = qe(20) + s*z3
! hexadecupole
            qe(21) = qe(21) + s*x4
            qe(22) = qe(22) + s*x3*y1
            qe(23) = qe(23) + s*x3*z1
            qe(24) = qe(24) + s*x2*y2
            qe(25) = qe(25) + s*z1*x2*y1
            qe(26) = qe(26) + s*x2*z2
            qe(27) = qe(27) + s*x1*z3
            qe(28) = qe(28) + s*z1*x1*y2
            qe(29) = qe(29) + s*z2*x1*y1
            qe(30) = qe(30) + s*z3*x1
            qe(31) = qe(31) + s*y4
            qe(32) = qe(32) + s*z1*y3
            qe(33) = qe(33) + s*z2*y2
            qe(34) = qe(34) + s*z3*y1
            qe(35) = qe(35) + s*z4
! spin dipole
            updens = 0.5D0*(rho(ind)*rho(ind + nxyz) + rho(ind))
            sux = sux + updens*x1
            suy = suy + updens*y1
            suz = suz + updens*z1
            s1 = s1 + updens
            dodens = -0.5D0*(rho(ind)*rho(ind + nxyz) - rho(ind))
            sdx = sdx + dodens*x1
            sdy = sdy + dodens*y1
            sdz = sdz + dodens*z1
            s2 = s2 + dodens !integral over downdensity
          END IF
        END DO
      END DO
    END DO

    s1 = s1*dvol
    s2 = s2*dvol


    se(1) = sux/(s1) - sdx/(s2) !spindipole-x
    se(2) = suy/(s1) - sdy/(s2) !spindipole-y
    se(3) = suz/(s1) - sdz/(s2) !spindipole-z
    se(4) = s1 ! total number spin-up
    se(5) = s2 ! total number spin-down

    DO k = 1, nrmom
      qe(k) = qe(k)*dvol
    END DO

    DO k = 2, nrmom
      qe(k) = qe(k)/qe(1) !normalization
    END DO

    DO k = 1, 3
      se(k) = se(k)*dvol
    END DO

    RETURN
  END SUBROUTINE emoms

! *********************

  SUBROUTINE laserp(vlaser, rho)

! **********************

! computes a laser potential corresponding to an electric
! field polarized in unit direction given by ux,uy,uz

! laser potential is :

! vlaser (x,y,z) = e_0 (x*ux + y*uy + z*uz) f(t) sin(omega*t)

! caution vlaser is stored as a 1-d vector (as any other grid field)

! the energy absorbed from the laser field is accumulated
! on the variable 'elaser' (stored in module 'params')

! photon pulsation is omega
! electric field intensity is e_0

! the laser pulse is characterized by
! its length deltat
! and peak time tpeak
! and begins at tnode

! pulse profile is given by f(t) (t is time)
! 4 options are presently implemented :

! itft = 1 --> ramp laser pulse, sine switching on/off up TO/from f(t) = 1

! tstart = tnode + tpeak
! tend = tstart + deltat
! tvend = tend + tpeak

! 0. < t < tnode f(t) = 0.
! tnode < t < tstart f(t) = sin(.5 *PI * (t-tnode) / tpeak)
! tstart < t < tend f(t) = 1.
! tend < t < tvend f(t) = sin(.5 *PI* (t - tend) / tpeak)
! tvend < t < ... f(t) = 0.

! itft = 2 --> Gaussian laser pulse
! tmax = tnode + tpeak
! 0. < t < tnode f(t) = 0.
! tnode < t < ... f(t) = exp(-((t-tmax)/deltat)**2)

! itft = 3 --> cos**2 envelope
! tnode = start of pulse
! tmax = tnode + deltat = END of pulse
! deltat = pulse duration
!
! this pulse allows a second cos**2 pulse (probe) switched by
! start time tstart2 and length 2*tpeak2.

! itft = 4 --> cos**4 envelope
! tnode = start of pulse
! tmax = tnode + deltat = END of pulse
! deltat = pulse duration

! The routine accumulates the integrated pulse profile in 'fpulseinteg1/2'.
! Trapezoidal integration is used and the last pulse value is
! saved in 'foft1/2old'.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: vlaser(kdfull2)
    REAL(DP), INTENT(IN) :: rho(2*kdfull2)

    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: acc1, acc2, accXUV
    REAL(DP) :: foft1, foft2, foftXUV
!REAL(DP) :: ppower,elaser,fpulseinteg1,fpulseinteg2,fXUVinteg
!REAL(DP) :: timeold,foft1old,foft2old,foftXUVold,acc1old,acc2old,accXUVold
    REAL(DP) :: tstart, tend, tmax, tmax2=0D0, tpulse2, tvend, deltime
    REAL(DP) :: scal1, scal2, snorm
!REAL(DP) :: ex1, ey1, ez1, ex2, ey2, ez2, exXUV, eyXUV, ezXUV
    REAL(DP) :: x1, y1, z1
    REAL(DP) :: tAPTstart, gaussenv, scalXUV
    INTEGER :: itrain, i

    IF (ABS(e0) <= 1D-20) THEN
      ind = 0
      DO iz = minz, maxz
        DO iy = miny, maxy
          DO ix = minx, maxx
            ind = ind + 1
            vlaser(ind) = 0.0D0
          END DO
        END DO
      END DO
      RETURN
    END IF

    IF (ilas == 0) THEN ! switch 'ilas' is maintained in 'params.F90'
      ilas = 1

!     normalize photon pulse orientations to uni vectors
      snorm = (e1x*e1x + e1y*e1y + e1z*e1z)
      snorm = SQRT(snorm)
      e1x = e1x/snorm
      e1y = e1y/snorm
      e1z = e1z/snorm
      IF(ABS(e0_2)>1D-10) THEN
        snorm = (e2x*e2x + e2y*e2y + e2z*e2z)
        IF(snorm == 0D0) STOP 'zero vector for orientation of 2. photon pulse'
        snorm = SQRT(snorm)
        e2x = e2x/snorm
        e2y = e2y/snorm
        e2z = e2z/snorm
      END IF
      datalaser(1:5) = 0D0 !ppower,elaser,fpulseinteg1,fpulseinteg2,fXUVinteg
      dataold = 0D0 !timeold,foft1old,foft2old,foftXUVold,acc1old,acc2old,accXUVold

    END IF

! prepare time profile

    foft2 = 0D0 ! standard CASE is: no second pulse
    foftXUV = 0D0 ! standard CASE: no attopulse train


    IF (itft == 1) THEN
      tstart = tnode + tpeak
      tend = tstart + deltat
      tvend = tend + tpeak
      tmax = (tnode + tvend)/2D0

      IF (tfs <= tnode .OR. tfs >= tvend) foft1 = 0D0
      IF (tfs > tnode .AND. tfs <= tstart) foft1 = SIN(0.5D0*pi*(tfs - tnode)/tpeak)
      IF (tfs > tstart .AND. tfs <= tend) foft1 = 1D0
      IF (tfs > tend .AND. tfs < tvend) foft1 = SIN(0.5D0*pi*(tvend - tfs)/tpeak)
    END IF

    IF (itft == 2) THEN
      tmax = tpeak
      foft1 = EXP(-((tfs - tpeak)/deltat)**2)
    END IF

    IF (itft == 3) THEN
      tpeak = deltat/2D0 ! is that needed ???
      tmax = deltat/2D0 + tnode
      tvend = tnode + deltat
      IF (tfs <= tnode .OR. tfs >= tvend) THEN
        foft1 = 0D0
      ELSE
        foft1 = COS((-0.5D0 + (tfs - tnode)/deltat)*pi)**2
      END IF

      IF (ABS(e0_2) > small) THEN
        tmax2 = deltat2/2D0 + tstart2
        IF (tfs >= tstart2 .AND. tfs <= tstart2 + deltat2) THEN
          foft2 = COS((-0.5D0 + (tfs - tstart2)/deltat2)*pi)**2
        END IF
      END IF


    END IF

    IF (itft == 4) THEN
      tmax = deltat/2D0 + tnode
      tvend = tnode + deltat
      IF (tfs <= tnode .OR. tfs >= tvend) THEN
        foft1 = 0D0
      ELSE
        foft1 = COS((-0.5D0 + (tfs - tnode)/deltat)*pi)**4
      END IF
    END IF

    IF (itft <= 0 .OR. itft >= 5) THEN
      STOP ' this pulse profile not yet implemented'
    END IF

! implement time profile in space

    foft1 = COS(omega*(tfs-tmax)/0.0484D0 + phi)*foft1
    foft2 = COS(omega2*(tfs-tmax2)/0.0484D0 + phase2)*foft2
    deltime = (tfs/0.0484D0 - timeold)
    acc1 = fpulseinteg1         ! save old value
    fpulseinteg1 = fpulseinteg1 + (foft1 + foft1old)*0.5D0*deltime
    acc2 = fpulse2integ1
    fpulse2integ1 = fpulse2integ1 + (fpulseinteg1+acc1)*0.5D0*deltime
    f2pulse2integ1 = f2pulse2integ1 + (fpulseinteg1**2+acc1**2)*0.5D0*deltime
    acc1 = fpulseinteg2         ! save old value
    fpulseinteg2 = fpulseinteg2 + (foft2 + foft2old)*0.5D0*deltime
    acc2 = fpulse2integ2
    fXUVinteg = fXUVinteg + (foftXUV + foftXUVold)*0.5D0*deltime
    fpulse2integ2 = fpulse2integ2 + (fpulseinteg2+acc1)*0.5D0*deltime
    f2pulse2integ2 = f2pulse2integ2 + (fpulseinteg2**2+acc1**2)*0.5D0*deltime

    WRITE (7, '(a,2f9.3,3(1pg12.4))') &
      ' tfs,tRy,foft1,foft2,foftXUV=', tfs, tfs/0.048D0, foft1, foft2, foftXUV

    ind = 0
    acc1 = 0D0
    acc2 = 0D0
    accXUV = 0D0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          scal1 = x1*e1x + y1*e1y + z1*e1z
          scal2 = x1*e2x + y1*e2y + z1*e2z
          vlaser(ind) = -e0*scal1*foft1 - e0_2*scal2*foft2
          acc1 = acc1 + e0*scal1*rho(ind)
          acc2 = acc2 + e0_2*scal2*rho(ind)
        END DO
      END DO
    END DO

    acc1 = acc1*dvol
    acc2 = acc2*dvol
    accXUV = accXUV*dvol
    elaser = elaser + (acc1 - acc1old)*(foft1 + foft1old)*0.5D0 &
             + (acc2 - acc2old)*(foft2 + foft2old)*0.5D0 &
             + (accXUV - accXUVold)*(foftXUV + foftXUVold)*0.5D0


    dataold(1) = tfs/0.0484D0 !timeold = tfs/0.0484D0
    dataold(2) = foft1 !foft1old = foft1
    dataold(3) = foft2 !foft2old = foft2
    dataold(4) = foftXUV !foftXUVold=foftXUV
    dataold(5) = acc1 !acc1old = acc1
    dataold(6) = acc2 !acc2old = acc2
    dataold(7) = accXUV !accXUVold = accXUV

    ex1 = e1x*foft1
    ey1 = e1y*foft1
    ez1 = e1z*foft1
    ppower = e0*e0*(ex1*ex1 + ey1*ey1 + ez1*ez1)
    IF (e0_2 /= 0D0) THEN
      ex2 = e2x*foft2
      ey2 = e2y*foft2
      ez2 = e2z*foft2
      ppower = ppower + e0_2*e0_2*(ex2*ex2 + ey2*ey2 + ez2*ez2)
    END IF


    RETURN
  END SUBROUTINE laserp

!-----projectp-----------------------------------------------------

  SUBROUTINE projectp(Vproj)

! Calculates the potential 'Vproj' from the point charge projectile
! on the valence electrons.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: Vproj(kdfull2)

    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: x1, y1, z1
!-----------------------------------------------------------------

    IF (ABS(projcharge) < 1D-3) THEN
      ind = 0
      DO iz = minz, maxz
        DO iy = miny, maxy
          DO ix = minx, maxx
            ind = ind + 1
            Vproj(ind) = 0D0
          END DO
        END DO
      END DO
    ELSE
      ind = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz
        DO iy = miny, maxy
          y1 = (iy - nysh)*dy
          DO ix = minx, maxx
            x1 = (ix - nxsh)*dx
            ind = ind + 1
            Vproj(ind) = -projcharge*e2/sqrt &
                         ((x1 - (projvelx*tfs + projinix))**2 + &
                          (y1 - (projvely*tfs + projiniy))**2 + &
                          (z1 - (projvelz*tfs + projiniz))**2)
          END DO
        END DO
      END DO
    END IF

    RETURN
  END SUBROUTINE projectp

! **************************

  REAL(DP) FUNCTION c_wfnorm(psi)

! Norm of wavefunction 'psi'. COMPLEX psi version

    USE params
    IMPLICIT NONE
    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2)

    INTEGER :: ii
    REAL(DP) :: acc
    acc = 0D0
    DO ii = 1, nxyz
      acc = REAL(psi(ii), DP)*REAL(psi(ii), DP) + AIMAG(psi(ii))*AIMAG(psi(ii)) + acc
    END DO
    c_wfnorm = acc*dvol
    RETURN
  END FUNCTION c_wfnorm

! **************************

  REAL(DP) FUNCTION r_wfnorm(psi1)

! Norm of wavefunction 'psi'. REAL psi version

    USE params
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: psi1(kdfull2)

    INTEGER :: ii
    REAL(DP):: acc

    acc = 0D0
    DO ii = 1, nxyz
      acc = (psi1(ii))*psi1(ii) + acc
    END DO
    r_wfnorm = acc*dvol
    RETURN
  END FUNCTION r_wfnorm

! **************************

  COMPLEX(DP) FUNCTION c_wfovlp(psi1, psi2)

! *************************

! Overlap <psi1|psi2>
! COMPLEX version

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: psi1(kdfull2)
    COMPLEX(DP), INTENT(IN) :: psi2(kdfull2)

    INTEGER :: ii
    COMPLEX(DP) :: csum

    csum = 0D0
    DO ii = 1, nxyz
      csum = CONJG(psi1(ii))*psi2(ii) + csum
    END DO
    c_wfovlp = csum*dvol
    RETURN
  END FUNCTION c_wfovlp
! **************************

  REAL(DP) FUNCTION r_wfovlp(psi1, psi2)

! *************************

! Overlap <psi1|psi2>
! REAL version

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: psi1(kdfull2)
    REAL(DP), INTENT(IN) :: psi2(kdfull2)

    INTEGER :: ii
    REAL(DP) :: acc

    acc = 0D0
    DO ii = 1, nxyz
      acc = (psi1(ii))*psi2(ii) + acc
    END DO
    r_wfovlp = acc*dvol

    RETURN
  END FUNCTION r_wfovlp

  SUBROUTINE rhointxy(rho, itime)

! **********************

! provides data file to analyse time-evolution of density.
! density 'rho' is integrated over x and y .
! the integrated density versus z is written on FILE 29
! for later analysis.
! at present a print modulus of 10 is hardwired.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rho(2*kdfull2)
    INTEGER, INTENT(IN) :: itime

    INTEGER :: ind, ix, iy, iz
    REAL(DP) ::acc
    REAL(DP) :: rhoint(minz:maxz)
    LOGICAL :: tfirst = .true.

    IF (jdensity1d == 0) THEN
      STOP 'irhoint must be larger than zero'
    END IF
    IF (MOD(itime, jdensity1d) /= 0) RETURN
    time = itime*dt1
    tfs = time*0.0484D0
    IF (tfirst) THEN
      OPEN (UNIT=28, FORM='UNFORMATTED', FILE='rho1Dz.'//outnam)
      WRITE (28) minz, maxz, nzsh, dz, centfz, dt1, itime, jdensity1d
      tfirst = .false.
      WRITE (6, *) ' output for integrated rho initialised'
    END IF

    ind = 0
    DO iz = minz, maxz
      acc = 0D0
      DO iy = miny, maxy
        DO ix = minx, maxx
          ind = ind + 1
          acc = rho(ind) + acc
        END DO
      END DO
      rhoint(iz) = acc*dx*dy
    END DO
    WRITE (28) tfs, rhoint

    RETURN
  END SUBROUTINE rhointxy

! *********************

  SUBROUTINE rhointyz(rho, itime)

! **********************

! provides data file to analyse time-evolution of density.
! density 'rho' is integrated over x and y .
! the integrated density versus z is written on FILE 29
! for later analysis.
! at present a print modulus of 10 is hardwired.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rho(2*kdfull2)
    INTEGER, INTENT(IN) :: itime

    LOGICAL :: tfirst
    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: acc
    REAL(DP) :: rhoint(minx:maxx)

    DATA tfirst/.true./

    IF (MOD(itime, jdensity1d) /= 0) RETURN
    time = itime*dt1
    tfs = time*0.0484D0
    IF (tfirst) THEN
      OPEN (UNIT=29, FORM='UNFORMATTED', FILE='rho1Dx.'//outnam)
      WRITE (29) minx, maxx, nxsh, dx, centfx, dt1, itime, jdensity1d
      tfirst = .false.
      WRITE (6, *) ' output for integrated rho initialised'
    END IF

    ind = 0
    DO ix = minx, maxx
      acc = 0D0
      DO iz = minz, maxz
        DO iy = miny, maxy
          ind = (iz - 1)*nxyf + (iy - 1)*nyf + ix
          acc = rho(ind) + acc
        END DO
      END DO
      rhoint(ix) = acc*dz*dy
    END DO
    WRITE (29) tfs, rhoint

    RETURN
  END SUBROUTINE rhointyz

! *********************

  SUBROUTINE rhointxz(rho, itime)

! **********************

! provides data file to analyse time-evolution of density.
! density 'rho' is integrated over x and y .
! the integrated density versus z is written on FILE 29
! for later analysis.
! at present a print modulus of 10 is hardwired.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rho(2*kdfull2)
    INTEGER, INTENT(IN) :: itime

    LOGICAL :: tfirst
    INTEGER :: ind, ix, iy, iz
    REAL(DP) :: acc
    REAL(DP) :: rhoint(miny:maxy)

    DATA tfirst/.true./

    IF (MOD(itime, jdensity1d) /= 0) RETURN
    time = itime*dt1
    tfs = time*0.0484D0
    IF (tfirst) THEN
      OPEN (UNIT=27, FORM='UNFORMATTED', FILE='rho1Dy.'//outnam)
      WRITE (27) miny, maxy, nysh, dy, centfy, dt1, itime, jdensity1d
      tfirst = .false.
      WRITE (6, *) ' output for integrated rho initialised'
    END IF

    ind = 0
    DO iy = miny, maxy
      acc = 0D0
      DO iz = minz, maxz
        DO ix = minx, maxx
          ind = (iz - 1)*nxyf + (iy - 1)*nyf + ix
          acc = rho(ind) + acc
        END DO
      END DO
      rhoint(iy) = acc*dx*dz
    END DO
    WRITE (27) tfs, rhoint

    RETURN
  END SUBROUTINE rhointxz

! **************************

  SUBROUTINE prifld(field, comment)

! print 'field' along x-axis for y=0 and z=0.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: field(kdfull2)
    CHARACTER(LEN=*), INTENT(IN) :: comment

    LOGICAL :: tnopri
    INTEGER :: ind, jx, jy, jz

    DATA tnopri/.false./

!c print field along x-axis

    IF (tnopri) RETURN
    WRITE (6, '(/2a)') 'field along x:', comment
    ind = 0
    DO jz = minz, maxz
      DO jy = miny, maxy
        DO jx = minx, maxx
          ind = 1 + ind
          IF (jz == nzsh .AND. jy == nysh) &
            WRITE (6, '(1x,f6.2,g13.5)') (jx - nxsh)*dx, field(ind)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE prifld

  SUBROUTINE prifldz(field, comment)

! print 'field' along z-axis for x=0 and y=0.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: field(kdfull2)
    CHARACTER(LEN=*), INTENT(IN) :: comment

    LOGICAL :: tnopri
    INTEGER :: ind, jx, jy, jz
    DATA tnopri/.false./

!c print field along z-axis

    IF (tnopri) RETURN
    WRITE (6, '(/2a)') 'field along z:', comment
    ind = 0
    DO jz = minz, maxz
      DO jy = miny, maxy
        DO jx = minx, maxx
          ind = 1 + ind
          IF (jx == nxsh .AND. jy == nysh) &
            WRITE (6, '(1x,f6.2,g13.5)') (jz - nzsh)*dx, field(ind)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE prifldz

! **************************

  SUBROUTINE prifld2(iunit, field, comment)

! print 'field' along y-axis for z=0 and y=0.
! 'iunit' is the file unit to which output is written.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    REAL(DP), INTENT(IN) :: field(kdfull2)
    CHARACTER(LEN=*), INTENT(IN) :: comment

    LOGICAL :: tnopri
    INTEGER :: ind, jx, jy, jz
    DATA tnopri/.false./

!c print field along x-axis

    IF (tnopri) RETURN
    WRITE (iunit, '(/2a)') 'field along x:', comment
    ind = 0
    DO jz = minz, maxz
      DO jy = miny, maxy
        DO jx = minx, maxx
          ind = 1 + ind
          IF (jz == nzsh .AND. jy == nysh) &
            WRITE (iunit, '(1x,f6.2,g13.5)') (jx - nxsh)*dx, field(ind)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE prifld2

! ******************************

  SUBROUTINE orthogwfr(wfr, isp, q0)

! ******************************

! Orthogonalizes REAL wavefunction 'wfr' with spin label 'isp'
! on the set of occupied states in 'q0'.
! Occupation numbers must be 0 or 1, i.e. 'temp=0.0' required.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(IN OUT) :: wfr(kdfull2)
    INTEGER, INTENT(IN) :: isp

    INTEGER :: i, nbe
    REAL(DP) :: overlap
!*********************************************************

! orthogonalization

    IF (temp /= 0D0) STOP ' ORTHOGWFR requires temperature 0'
    DO nbe = 1, nstate
      IF (ispin(nbe) == isp .AND. occup(nbe) > 0.999D0) THEN
        overlap = wfovlp(wfr, q0(:, nbe))
        DO i = 1, nxyz
          wfr(i) = wfr(i) - overlap*q0(i, nbe)
        END DO
      END IF
    END DO

    RETURN
  END SUBROUTINE orthogwfr

!------------------------------------------------------------
  REAL(DP) FUNCTION ran1(idum)
!------------------------------------------------------------
    USE params, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT)::idum
    INTEGER :: ia, im, iq, ir, ntab, ndiv
    REAL(DP) :: am, eps, rnmx
    PARAMETER(ia=16807, im=2147483647, am=1D0/im, iq=127773, ir=2836, &
              ntab=32, ndiv=1 + (im - 1)/ntab, eps=1.2D-7, rnmx=1D0 - eps)

! Minimal random number generator of Park and Miller with Bays-Durham
! shuffle and added safeguards. Returns a uniform random deviate
! between 0.0 and 1.0 (exclusive of the endpoints). CALL with idum
! a negative INTEGER to initialize; thereafter, do not alter
! idum between successive deviates in a sequence. RNMX should
! approximate the largest floating value that is less than 1.

    INTEGER :: j, k, iv(ntab), iy
    SAVE iv, iy
    DATA iv/ntab*0/, iy/0/

    IF (idum <= 0 .OR. iy == 0) THEN
      idum = MAX(-idum, 1)
      DO j = ntab + 8, 1, -1
        k = idum/iq
        idum = ia*(idum - k*iq) - ir*k
        IF (idum < 0) idum = idum + im
        IF (j <= ntab) iv(j) = idum
      END DO
      iy = iv(1)
    END IF
    k = idum/iq
    idum = ia*(idum - k*iq) - ir*k
    IF (idum < 0) idum = idum + im
    j = 1 + iy/ndiv
    iy = iv(j)
    iv(j) = idum
    ran1 = MIN(am*iy, rnmx)

    RETURN
  END FUNCTION ran1
!------------------------------------------------------------

!------------------------------------------------------------
  REAL(DP) FUNCTION gasdev(idum)

! Random Gaussian distribution

    USE params, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN OUT) :: idum

    INTEGER :: iset
    REAL(DP) :: fac, gset, rsq, v1, v2
    SAVE iset, gset

    DATA iset/0/
    IF (iset == 0) THEN
      rsq = 2D0 ! some value higher than 1.0 to start to enter the loop
      DO WHILE (rsq >= 1D0 .OR. rsq == 0D0)
        v1 = 2D0*ran1(idum) - 1D0
        v2 = 2D0*ran1(idum) - 1D0
        rsq = v1**2 + v2**2
      END DO
      fac = SQRT(-2D0*LOG(rsq)/rsq)
      gset = v1*fac
      gasdev = v2*fac
      iset = 1
    ELSE
      gasdev = gset
      iset = 0
    END IF

    RETURN
  END FUNCTION gasdev
!------------------------------------------------------------

!------------------------------------------------------------
  SUBROUTINE givetemperature(pxt, pyt, pzt, nteil, temperat, masse, ipflag)
!------------------------------------------------------------

! Produces thermal distribution of ionic velocities.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: pxt(:)
    REAL(DP), INTENT(OUT) :: pyt(:)
    REAL(DP), INTENT(OUT) :: pzt(:)
    INTEGER, INTENT(IN) :: nteil
    REAL(DP), INTENT(IN) :: temperat
    REAL(DP), INTENT(IN) :: masse
    INTEGER, INTENT(IN) :: ipflag

! temper given in Kelvin, convert to kT in Ryd now

    INTEGER :: i, irand, ntmob
    REAL(DP) :: ekhaben, eksoll, fac, sumx, sumy, sumz, sc, temper

    WRITE (6, *) 'ENTERING GIVETEMPERATURE'

    temper = temperat*6.507D-6

    fac = SQRT(temper/masse)

    sumx = 0D0
    sumy = 0D0
    sumz = 0D0
    irand = -4

    DO i = 1, nteil
      pxt(i) = gasdev(irand)*fac
      pyt(i) = gasdev(irand)*fac
      pzt(i) = gasdev(irand)*fac
      sumx = sumx + pxt(i)
      sumy = sumy + pyt(i)
      sumz = sumz + pzt(i)
    END DO

    DO i = 1, nteil
      pxt(i) = pxt(i) - sumx/nteil
      pyt(i) = pyt(i) - sumy/nteil
      pzt(i) = pzt(i) - sumz/nteil
    END DO

    ekhaben = 0D0

    ntmob = 0

    IF (ipflag == 4) THEN
      IF (ionmdtyp /= 0) ntmob = nion
    ELSE
      STOP 'Error in routine giveTemperature'
    END IF

    eksoll = (3D0*ntmob - 3D0)/2D0*temper

! WRITE(6,*) 'eksoll: ',eksoll

    DO i = 1, nteil
      IF (ipflag == 4 .AND. ionmdtyp /= 0) THEN
        ekhaben = ekhaben + pxt(i)**2 + pyt(i)**2 + pzt(i)**2
      END IF
    END DO

    ekhaben = ekhaben/2D0/masse

    sc = SQRT(eksoll/ekhaben)

    DO i = 1, nteil
      pxt(i) = pxt(i)*sc
      pyt(i) = pyt(i)*sc
      pzt(i) = pzt(i)*sc
    END DO

    ekhaben = 0D0

    DO i = 1, nteil
      IF (ipflag == 4 .AND. ionmdtyp /= 0) THEN
        ekhaben = ekhaben + pxt(i)**2 + pyt(i)**2 + pzt(i)**2
      END IF
    END DO

    ekhaben = ekhaben/2D0/masse

    RETURN
  END SUBROUTINE givetemperature
!------------------------------------------------------------

!-----safeopen-----------------------------------------------

  SUBROUTINE safeopen(iunit, iact, icond, filename)

! Opens file unit 'iunit' and asks for status before.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit, iact, icond
    CHARACTER(LEN=*), INTENT(IN) :: filename

    LOGICAL :: topenf

    IF (icond > 0) THEN
      IF (MOD(iact, icond) == 0) THEN
        INQUIRE (iunit, OPENED=topenf)
        IF (.NOT. topenf) OPEN (iunit, POSITION='append', &
                                FILE=filename//'.'//outnam)
      END IF
    END IF

    RETURN
  END SUBROUTINE safeopen
! ******************************

! ******************************
  SUBROUTINE cleanfile(iunit)

! Closes FILE unit 'iunit' IF it was OPEN.

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iunit
    LOGICAL :: topen

    INQUIRE (iunit, OPENED=topen)
    IF (topen) CLOSE (iunit)

    RETURN
  END SUBROUTINE cleanfile
! ******************************

! ******************************

  SUBROUTINE r_project(qin, qout, ispact, q0)

! ******************************
! projects all occupied states 'q0' OUT of 'qin'.
! q0 = set of s.p. wavefunctions (REAL)
! qin = wavefunction from which 'q0' are to be removed (REAL)
! qout = resulting wavefunction
! ispact = spin of 'qin'

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(IN) :: qin(kdfull2)
    REAL(DP), INTENT(OUT) :: qout(kdfull2)
    INTEGER, INTENT(IN) :: ispact

    INTEGER :: i, nbe
    REAL(DP) :: ovl
!*********************************************************

    DO i = 1, nxyz
      qout(i) = qin(i)
    END DO

    DO nbe = 1, nstate
      IF (ispin(nbe) == ispact .AND. occup(nbe) > 0.99999999D0) THEN
        ovl = 0D0
        DO i = 1, nxyz
          ovl = ovl + q0(i, nbe)*qout(i)
        END DO
        ovl = ovl*dvol
        DO i = 1, nxyz
          qout(i) = qout(i) - ovl*q0(i, nbe)
        END DO
      END IF
    END DO

    RETURN
  END SUBROUTINE r_project
! ******************************

  SUBROUTINE c_project(qin, qout, ispact, q0)

! ******************************
! projects all occupied states 'q0' out of 'qin'.
! q0 = set of s.p. wavefunctions (COMPLEX)
! qin = wavefunction from which 'q0' are to be removed (COMPLEX)
! qout = resulting wavefunction
! ispact = spin of 'qin'

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
    COMPLEX(DP), INTENT(IN) :: qin(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: qout(kdfull2)
    INTEGER, INTENT(IN) :: ispact

    INTEGER :: nbe
    COMPLEX(DP) :: ovl

!*********************************************************

    qout = qin

    DO nbe = 1, nstate
      IF (ispin(nbe) == ispact .AND. occup(nbe) > 0.99999999D0) THEN
        ovl = dvol*SUM(CONJG(q0(:, nbe))*qout)
        qout(:) = qout(:) - q0(:, nbe)*ovl
      END IF
    END DO

    RETURN
  END SUBROUTINE c_project


!-----pair ------------ part of the Hartree-Fock package --------------

  SUBROUTINE pair(e, gw, ph, nz, nmax, gp, eferm, delta, partnm, &
                  iter, ipair, eps, iab, kstate)
    USE params, ONLY: DP, one, zero
    IMPLICIT NONE

! * * * * * * * * * * pairing iteration * * * * * * * * * *

! determines the Fermi-energy 'eferm' and the gap 'delta'
! such that the "gap-equation" and the "particle number condition"
! are fulfilled simultaneously. this is achieved by iteration
! with a gradient scheme (newton's tangential formula)
! where the gradient is evaluated by finite differences.

! there are two schemes possible:
! 1. the "constant force" approach, where the pairing force is given
! in 'gp'.
! 2. the "constant gap" approach, where the gap is kept fixed and
! is given in 'delta'.
! there is only 'eferm' to be determined in the "constant gap"
! approach. the gap-equation is used at the END TO determine an
! according force 'gp'.

! input via list:
! e = array(1:nmax) of single particle energies.
! gw = array(1:nmax) of occupation probabilities,
! input is used as initial guess for the iteration.
! ph = array(1:nmax) of the multiplicities of the single
! particle levels.
! nz = desired particle number
! nmax = maximum number of single particle states
! (needs to be smaller than the PARAMETER 'kstate'
! and smaller than 46, for dimensioning reasons).
! gp = size of pairing-force, is input in CASE of
! "constant force" approach.
! eferm = Fermi-energy, input is used as initial guess;
! IF eferm=0.0 is input THEN a suitable initial guess
! is computed from input 'gw'.
! delta = pairing gap, is input in CASE of "constant gap" approach
! iter = number of pairing iterations,
! abs(iter) is the number of iterations
! ipair = switch for pairing TYPE
! ipair=1 --> "constant force" approach, using
! pairing force 'gp'
! ipair=0 --> "constant gap" approach, using
! given gap 'delta'
! ipair=4 --> thermal occupation with temperature
! given by 'delta'
! eps = end condition for pairing iteration;
! iteration terminates if error in particle number
! is smaller than 'eps' and if 'eferm' and/or 'delta'
! are smaller than 'eps'.
! iab = print level; significant values are -2, -1, 0 and 1;
! larger value prints more.

! output via list:
! gw = array(1..nmax) of the occupation probabilities
! finally achieved.
! gp = size of pairing-force, is output in case of
! "constant gap" approach.
! delta = pairing gap, is output in case of
! "constant force" approach.
! eferm = Fermi-energy.
! sum = particle number, from adding 'gw' with degeneracy 'ph'.

! some auxiliary fields are:
! ph = array(1:46) of the multiplicities of the single
! particle levels. array is ordered according to
! the harmonic oscillator scheme, with rising shell
! number and lowering angular momentum within each shell.
! dparn = array(1:3) of variation in the particle number as
! consequence of 'deferm' and 'ddelta'.
! dgap = array(1:3) of variation in the gap equation as
! consequence of 'deferm' and 'ddelta'.
! gr = double array of variation of particle number and of
! gap equation, used for construction of the tangential
! step

    INTEGER, INTENT(IN) :: kstate
    REAL(DP), INTENT(IN) :: e(kstate)
    REAL(DP), INTENT(IN OUT) :: gw(kstate)
    REAL(DP), INTENT(IN OUT) :: ph(kstate)
    INTEGER, INTENT(IN) :: nz
    INTEGER, INTENT(IN) :: nmax
    REAL(DP), INTENT(IN OUT) :: gp
    REAL(DP), INTENT(IN OUT) :: eferm
    REAL(DP), INTENT(IN OUT) :: delta
    REAL(DP), INTENT(IN OUT) :: partnm
    INTEGER, INTENT(IN) :: iter
    INTEGER, INTENT(IN) :: ipair
    REAL(DP), INTENT(IN) ::eps
    INTEGER, INTENT(IN) :: iab

    LOGICAL :: converged = .false. ! will be true when convergence is reached
    INTEGER :: i, isum, it, itmaxp, itsum, k, ka
    REAL(DP):: dampin, deferm, ddelta, delti, delact, delold, del2
    REAL(DP):: elam, equasi, efold, efstep, parnum, gapeq, gphalf, grsav
    REAL(DP):: emax, emin
    REAL(Dp):: eferm1, eferm2, eferm3, parn1, parn2, parn3
    REAL(DP):: efmax, dlmax, det, detmin, defmin, ddlmin, dlstep
    REAL(DP):: xmaxlw, x10lw
    REAL(DP):: dparn(3), dgap(3), gr(2, 2)

    DATA gr, dparn, dgap/10*0D0/

! further numerical variables are:
! deferm = step in the Fermi-energy, to evaluate
! numerically the derivatives.
! ddelta = step in the gap, to evaluate numerically the
! derivatives; used only in the "constant force" approach
! dampin = damping factor for the tangential step
! efmax,dlmax,detmin,defmin,ddlmin = limitations for the step

! note: these parameters are adjusted for an application with
! the nuclear shell model. it is assumed that the single
! particle energies are given in units of [mev] and display
! a reasonable nuclear spectrum.

    DATA dampin, deferm, ddelta/0.9D0, 0.1D0, 0.1D0/
    DATA efmax, dlmax, detmin, defmin, ddlmin/3D0, .6D0, .1D-3, .1D-3, .1D-3/
    DATA xmaxlw, x10lw/1.0D-30, 1.0D-10/

!-----------------------------------------------------------------------

    IF (iab > 1) THEN

      WRITE (7, '(/a,2i4,3(a,f9.4)/a,10(/1x,10f8.3))') &
        ' *** pair entered with: iter,ipair=', iter, ipair, ' delta=', delta, &
        ' eferm=', eferm, ' gp=', gp, &
        ' single particle energies:', (e(i), i=1, nmax)
      WRITE (6, '(/a,2i4,3(a,f9.4)/a,10(/1x,10f8.3))') &
        ' *** pair entered with: iter,ipair=', iter, ipair, ' delta=', delta, &
        ' eferm=', eferm, ' gp=', gp, &
        ' single particle energies:', (e(i), i=1, nmax)
      WRITE (6, *) 'kstate', kstate
      WRITE (7, *) 'kstate', kstate
      WRITE (6, *) 'nmax', nmax
      WRITE (7, *) 'nmax', nmax
      WRITE (6, '(a,500f6.2)') 'gw', (gw(i), i=1, nmax)
      WRITE (6, '(a,500f6.2)') 'ph', (ph(i), i=1, nmax)
      WRITE (6, *) 'nz', nz
      WRITE (7, *) 'gw', (gw(i), i=1, nmax)
      WRITE (7, *) 'ph', (ph(i), i=1, nmax)
      WRITE (7, *) 'nz', nz
      WRITE (6, *) partnm, ipair, eps, iab
      WRITE (7, *) partnm, ipair, eps, iab
    END IF

! prepare Fermi-energy 'eferm', if not yet initialised.

    IF (eferm == zero) THEN
      isum = 0
      DO it = 1, nmax
        itsum = it
        isum = INT(ph(it) + .000001D0) + isum
        IF (isum + isum - nz >= 0) EXIT
      END DO
      IF (it > nmax) STOP ' not enough states to reach particle number'
      eferm = e(itsum) + deferm
      IF (iab >= 1) WRITE (7, '(1x,a,i5,(a,g11.4),a,i5)') &
        'initialistion: isum=', isum, ' eferm=', eferm, ' nz=', nz

    END IF

    IF (ipair == 1) THEN

!-----------------------------------------------------------------------

! case of constant pairing-force, given by 'gp'.
! gap 'delta', Fermi-energy 'eferm', and occupations 'gw'
! are computed.

!-----------------------------------------------------------------------

      delti = delta
      IF (ABS(delti) < x10lw) delti = one
      DO it = 1, ABS(iter)

! gap-equation & particle-number at 'eferm' & 'delti' (for k=1)
! and at 'eferm' varied (k=2) and at 'delta' varied (k=3).

        DO ka = 1, 3
! actual Fermi-energy and gap
          k = 4 - ka
          IF (k == 1) THEN
            elam = eferm
            delact = delti
          ELSE IF (k == 2) THEN
            elam = eferm + deferm
            delact = delti
          ELSE
            elam = eferm
            delact = delti + ddelta
          END IF
! compute actual quasiparticle e
! and occupation-weights,
! accumulate gap-equation and
! particle number
          parnum = zero
          gapeq = zero
          DO i = 1, nmax
            equasi = SQRT((e(i) - elam)*(e(i) - elam) + delact*delact)
            gw(i) = 0.5D0 - 0.5D0*(e(i) - elam)/equasi
            gapeq = ph(i)/equasi + gapeq
            parnum = gw(i)*ph(i) + parnum
          END DO
! store actual part.number and
! gap-equation
          dgap(k) = 0.25D0*gp*gapeq - one !?? extra 1/2 because of w
          dparn(k) = parnum - nz
        END DO

! construction of the gradient-matrix and its inverse
! ( determinant 'det' is limited in value to avoid singularities

        gr(1, 1) = (dgap(2) - dgap(1))/deferm
        gr(1, 2) = (dparn(2) - dparn(1))/deferm
        gr(2, 1) = (dgap(3) - dgap(1))/ddelta
        gr(2, 2) = (dparn(3) - dparn(1))/ddelta
        det = gr(1, 1)*gr(2, 2) - gr(1, 2)*gr(2, 1)
        det = SIGN(MAX(ABS(det), detmin), det)
        IF (ABS(det) < xmaxlw) det = SIGN(xmaxlw, det)

        IF (iab >= 1) THEN
          WRITE (7, '(a,i4,3(a,g11.3)/a,3g11.3,a,3g11.3/a,4g11.3)') &
            ' iteration at it=', it, ' : delti=', delti, ' eferm=', eferm, &
            ' det=', det, ' dgap=', dgap, ' , dparn=', dparn, &
            ' gr=', gr
          WRITE (7, '(a,10f7.4,10(/6x,10f7.4))') ' gw=', gw
        END IF

        grsav = gr(1, 1)
        gr(1, 1) = gr(2, 2)/det
        gr(1, 2) = -gr(1, 2)/det
        gr(2, 1) = -gr(2, 1)/det
        gr(2, 2) = grsav/det

! tangential step (including some step limiters)

        efold = eferm
        delold = delti
        efstep = dampin*(gr(1, 1)*dgap(1) + gr(2, 1)*dparn(1))
        dlstep = dampin*(gr(1, 2)*dgap(1) + gr(2, 2)*dparn(1))
        eferm = eferm - SIGN(MIN(ABS(efstep), efmax), efstep)
        delti = ABS(delti - SIGN(MIN(ABS(dlstep), dlmax), dlstep))
        deferm = efold - eferm
        ddelta = delold - delti
        IF (ABS(deferm) < defmin) deferm = defmin
        IF (ABS(ddelta) < ddlmin) ddelta = ddlmin

        IF (iab == 0) WRITE (7, '(a,i4,2(a,g11.3))') &
          ' pairing iteration it=', it, ' : delti=', delti, ' eferm=', eferm

        IF (ABS(dgap(1)) < eps .AND. ABS(dparn(1)) < eps) THEN
          converged = .true.
          EXIT
        END IF
      END DO

      delta = delti

    ELSE

!-----------------------------------------------------------------------

! case of constant paring-gap, given by 'delta'.
! using now the safer secant step.

!-----------------------------------------------------------------------

! determine upper and lower bound for secant step

      IF (iab > 1) WRITE (*, *) ' bracketing branch'
      converged = .FALSE.
      emin = +1.0D30
      emax = -1.0D30
      DO i = 1, nmax
        emin = MIN(e(i), emin)
        emax = MAX(e(i), emax)
      END DO
      itmaxp = ABS(iter)
      IF (itmaxp < 0) THEN
        eferm1 = eferm
        eferm2 = eferm
        eferm3 = eferm
        parn1 = parnm(e, gw, ph, nmax, delta, eferm, ipair)
        parn2 = parn1
        parn3 = parn1
        IF (parn1 - nz > zero) THEN
          deferm = (eferm - emin)/itmaxp
          IF (iab >= 1) WRITE (7, '(a,3(1pg12.4))') &
            ' search brackets: parn2,eferm,emin=', parn2, eferm, emin
          DO it = 1, itmaxp
            eferm1 = eferm1 - deferm
            parn1 = parnm(e, gw, ph, nmax, delta, eferm1, ipair)
            IF (iab >= 1) WRITE (7, '(a,i5,3(1pg12.4))') &
              ' search brackets: it,parn1,eferm1=', it, parn1, eferm1
            IF (parn1 - nz < zero) EXIT
            eferm2 = eferm1
          END DO
          IF (it > itmaxp) STOP ' particle number not embraced in pair'
        ELSE IF (parn1 - nz < zero) THEN
          deferm = (emax - eferm)/itmaxp
          IF (iab >= 1) WRITE (7, '(a,3(1pg12.4))') &
            ' search brackets: parn1,eferm,emax=', parn1, eferm, emax
          DO it = 1, itmaxp
            eferm2 = eferm2 + deferm
            parn2 = parnm(e, gw, ph, nmax, delta, eferm2, ipair)
            IF (iab >= 1) WRITE (7, '(a,i5,3(1pg12.4))') &
              ' search brackets: it,parn2,eferm2=', it, parn2, eferm2
            IF (parn2 - nz > zero) EXIT
            eferm1 = eferm2
          END DO
          IF (it > itmaxp) STOP ' particle number not embraced in pair'
        ELSE
          IF (iab >= 1) WRITE (7, '(a)') ' solution immediately found '
          converged = .true. ! jump to END
        END IF

        IF (.NOT. converged) THEN
! Perform the secant search
          IF (iab >= 1) WRITE (7, '(a)') ' brackets found '

          DO it = 1, itmaxp
            eferm3 = eferm1 + (nz - parn1)*(eferm2 - eferm1)/(parn2 - parn1)
            parn3 = parnm(e, gw, ph, nmax, delta, eferm3, ipair)
            IF (iab >= 0) WRITE (7, '(a,i4,2(a,g11.3))') &
              ' secant iteration it=', it, ' : eferm=', eferm3, ' dparn(1)=', parn3 - nz
            IF (ABS(parn3 - nz) < eps) THEN
              converged = .true.
              EXIT
            ELSE IF (parn3 - nz > zero) THEN
              eferm2 = eferm3
              parn2 = parn3
            ELSE
              eferm1 = eferm3
              parn1 = parn3
            END IF
          END DO
        END IF
      END IF

      IF (.NOT. converged) THEN
! the secant step did not succeed: use bisection
        eferm1 = emin - 4*delta
        eferm2 = emax + 4*delta
        DO it = 1, itmaxp
          eferm3 = 0.5D0*(eferm2 + eferm1)
          parn3 = parnm(e, gw, ph, nmax, delta, eferm3, ipair)
          IF (iab >= 0) WRITE (7, '(a,i4,2(a,g11.3))') &
            ' bisection it=', it, ' : eferm=', eferm3, ' dparn(1)=', parn3 - nz
          IF (ABS(parn3 - nz) < eps) THEN
            converged = .true.
            EXIT
          ELSE IF (parn3 - nz < zero) THEN
            eferm1 = eferm3
            parn1 = parn3
          ELSE IF (parn3 - nz > zero) THEN
            eferm2 = eferm3
            parn2 = parn3
          END IF
        END DO
      END IF

      eferm = eferm3
      dparn(1) = parn3 - nz
      dgap(1) = zero
      partnm = parn3

! determine force 'gp' according to given 'delta' and the
! equilibrium reached.

      IF (ipair == 0) THEN
        gapeq = zero
        DO i = 1, nmax
          gapeq = ph(i)/SQRT((e(i) - eferm)*(e(i) - eferm) + delta*delta) + gapeq
        END DO
        gp = 4D0/gapeq !?? account for wate=2.
      END IF

! END of big switch between cases

    END IF

!-----------------------------------------------------------------------

! finally the gap-equation and particle-number are checked and print

!-----------------------------------------------------------------------

! in case of incomplete convergence print a warning

    IF (iab >= -1 .AND. .NOT. converged) WRITE (7, '(2(a,g11.3),2(a,i3))') &
      ' ---> no convergence in pair! ERR(gap-eq)=', dgap(1), &
      ' ERR(part.nr)=', dparn(1), ' nz=', nz, ' nmzx=', nmax

! check gap equation

    IF (iab >= 0) THEN
      IF (ipair /= 4) THEN
        gphalf = gp*0.5D0
        del2 = MAX(xmaxlw, delta*delta)
        gapeq = zero
        partnm = zero
        DO i = 1, nmax
          partnm = gw(i)*ph(i) + partnm
          gapeq = gphalf*ph(i)/SQRT(del2 + (e(i) - eferm)**2) + gapeq
        END DO
      ELSE
        gapeq = zero
      END IF

      WRITE (7, '(a,i8,2(a,g13.5))') &
        ' pair finished with ', it, ' iterations: gap-eq.=', gapeq, ' part.nr=', partnm
    END IF

    RETURN
  END SUBROUTINE pair
!-----parnm------------ part of the Hartree-Fock package ---------------

  REAL(DP) FUNCTION parnm(e, gw, ph, nmax, delta, elam, ipair)

! the include file 'params.inc' communicates the PRECISION
! transfering a line 'IMPLICIT DOUBLE PRECISION (a-h,o-z)' or not.
! furthermore it carries the COMMON BLOCK /constn/ from which
! the pseudo-constant one=1.0 is taken here.

! the function computes the particle number in a pairing or
! temperature distribution. the parameters are:
! e = array of s.p. energies
! gw = array of occupation weights (is workspace and output!)
! ph = array of degeneracy factors
! nmax = number of states in the above three arrays
! delta = pairing gap or temperature
! elam = Fermi energy
! ipair = switch to pairing (0) or temperature (4).

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: e(ksttot)
    REAL(DP), INTENT(OUT) :: gw(ksttot)
    REAL(DP), INTENT(IN) :: ph(ksttot)
    INTEGER, INTENT(IN) :: nmax
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), INTENT(IN) :: elam
    INTEGER, INTENT(IN) :: ipair

    INTEGER :: i
    REAL(DP) :: acc, equasi

!----------------------------------------------------------------------

    acc = zero
    IF (ipair == 0) THEN
      DO i = 1, nmax
        equasi = SQRT((e(i) - elam)*(e(i) - elam) + delta*delta)
        gw(i) = 0.5D0 - 0.5D0*(e(i) - elam)/equasi
        acc = gw(i)*ph(i) + acc
      END DO
    ELSE IF (ipair == 4) THEN ! case of temperature
      DO i = 1, nmax
        IF (ph(i) == 0D0) CYCLE
        equasi = (e(i) - elam)/delta
        IF (equasi > 50D0) THEN
          gw(i) = 0D0
        ELSE IF (equasi < -50D0) THEN
          gw(i) = 1D0
        ELSE
          gw(i) = 1D0/(1D0 + EXP(equasi))
        END IF
        acc = gw(i)*ph(i) + acc
      END DO
    END IF

    parnm = acc

    RETURN
  END FUNCTION parnm


!-----timer------------------------------------------------------

  SUBROUTINE timer(iaction)

! Initializes timer and retrieves CPU times
! iaction = 1 for initial call
! 2 for subsequent call
! The relative time is taken between step 2 and 1.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iaction

    INTEGER :: i
    REAL(DP), ALLOCATABLE, SAVE :: cputimenew(:), cputimeold(:)

    LOGICAL, SAVE :: tfirst = .true.

!-----------------------------------------------------------------

    IF (tfirst) THEN
      ALLOCATE (cputimenew(0:knode), cputimeold(0:knode))
      tfirst = .FALSE.
    END IF

    IF (iaction == 1) THEN
      CALL cpu_time(cputimeold(0))
    ELSE IF (iaction == 2) THEN
      CALL cpu_time(cputimenew(0))
      IF (myn == 0) WRITE (6, '(a,1pg20.10)') 'time for step in node 0=', &
        cputimenew(0) - cputimeold(0)
      DO i = 0, knode
        cputimeold(i) = cputimenew(i)
      END DO
    ELSE
      STOP 'TIMER: this IACTION is not valid'
    END IF

    RETURN
  END SUBROUTINE timer

!-----stimer------------------------------------------------------

  SUBROUTINE stimer(iaction, text)

! Initializes timer and retrieves CPU times
! iaction = 1 for initial call
! 2 for subsequent acll
! The relative time is taken between step 2 and 1.
! text is ussed in print

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iaction
    CHARACTER(*), INTENT(IN) :: text

    INTEGER, SAVE :: itimeold = 0, itimenew = 0, irate = 0

!-----------------------------------------------------------------

    IF (irate == 0) CALL system_clock(itimeold, irate)

    IF (iaction == 1) THEN
      CALL system_clock(itimeold)
    ELSE IF (iaction == 2) THEN
      CALL system_clock(itimenew)
      IF (myn == 0) WRITE (6, '(a,1pg13.5)') 'systime '//text//'=', &
        (itimenew - itimeold)*1D0/irate
      IF (myn == 0) WRITE (7, '(a,1pg13.5)') 'systime '//text//'=', &
        (itimenew - itimeold)*1D0/irate
      itimeold = itimenew
    ELSE
      STOP 'SYS.TIMER: this action is not valid'
    END IF

    RETURN
  END SUBROUTINE stimer

!---probab----------------------------------------------------------

  SUBROUTINE probab(psitmp)

! Calculates the probabilities to find a certain charge stae from
! the distribution of loss of norms of the s.p. wavefunctions
! given set of wavefunctions in input array 'psitmp'.

    USE params
    IMPLICIT NONE
    COMPLEX(DP), INTENT(IN) :: psitmp(kdfull2, kstate)

    INTEGER :: k
    REAL(DP) :: prob(nstate + 1), pold(nstate + 1)
    REAL(DP) :: rtmpuse(nstate)

    INTEGER :: i, it
    COMPLEX(DP) :: cscal
!-----------------------------------------------------------------

    DO i = 1, nstate
      cscal = orbitaloverlap(psitmp(:, i), psitmp(:, i))
      rtmpuse(i) = SQRT(REAL(cscal, DP)**2 + AIMAG(cscal)**2)*occup(i)
      prob(i) = 0D0
    END DO

    prob(nstate + 1) = 0D0

    prob(1) = rtmpuse(1)
    prob(2) = 1D0 - rtmpuse(1)

    DO i = 2, nstate
      DO k = 1, i
        pold(k) = prob(k)
      END DO

      DO k = 1, i
        prob(k) = pold(k)*rtmpuse(i)
      END DO

      DO k = 2, i + 1
        prob(k) = prob(k) + pold(k - 1)*(1D0 - rtmpuse(i))
      END DO
    END DO

    CALL safeopen(808, it, jnorms, 'pproba')
    DO k = nstate + 1, 1, -1
      WRITE (808, '(f12.6,i4,1pg13.5)') tfs, nstate + 1 - k, prob(k)
    END DO
    WRITE (808, '(1x)')
    FLUSH (808)

    RETURN
  END SUBROUTINE probab

!-----------------------------------------------------------------

!-----stateoverl---------------------------------------------------------

  SUBROUTINE stateoverl(psi1, psi2)

! Computes overlap of two Slater states with possibly
! different sets of s.p. states.
! Both states must have the same sequence of occupation numbers.
! The set of s.p. states is given by 'psi1' and 'psi2'.
! The occupation numbers 'occup' are communicated via 'params'.
! All states should have occupation 1.

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: psi1(kdfull2, kstate)
    COMPLEX(DP), INTENT(IN) :: psi2(kdfull2, kstate)

    INTEGER :: nbe, nb2
    COMPLEX(DP), ALLOCATABLE :: spoverlaps(:, :)

    LOGICAL, SAVE :: tfirst = .true.
    LOGICAL, PARAMETER :: ttest = .false.

    IF (TFIRST) THEN
      OPEN (382, FILE='povlp.'//outnam)
      WRITE (382, '(a)') &
        '# overlap of static state with dynamically propagated one:', &
        '# time (fs) overlap'
      DO nbe = 1, nstate
        IF (occup(nbe) < 0.9999999999999D0) &
          STOP ' STATEOVERL requires full occupation'
      END DO
      tfirst = .false.
    END IF

    ALLOCATE (spoverlaps(nstate, nstate))

    DO nbe = 1, nstate
      DO nb2 = 1, nstate
        IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(nb2))) THEN
          spoverlaps(nb2, nbe) = SUM(CONJG(psi1(:, nbe))*psi2(:, nb2))*dvol
        ELSE
          spoverlaps(nb2, nbe) = 0D0
        END IF
      END DO
    END DO

    WRITE (382, '(f10.4,2(1pg13.5))') time, determinant(spoverlaps, nstate, nstate)
    IF (ttest) THEN
      WRITE (382, '(a)') ' matrix up:'
      DO nbe = 1, nstate/2
        WRITE (382, '(16(1pg13.5))') (spoverlaps(nb2, nbe), nb2=1, nstate/2)
      END DO
      WRITE (382, '(a)') ' matrix down:'
      DO nbe = nstate/2 + 1, nstate
        WRITE (382, '(16(1pg13.5))') (spoverlaps(nb2, nbe), nb2=nstate/2 + 1, nstate)
      END DO
    END IF
    FLUSH (382)

    DEALLOCATE (spoverlaps)

    RETURN

  END SUBROUTINE stateoverl

!------------------------------------------------------------

  SUBROUTINE shiftsmall(q0, shix, shiy, shiz)
!------------------------------------------------------------
    USE params
    USE kinetic
    IMPLICIT NONE

! Shifts COMPLEX wavefunction array 'q0' by 'shix', 'shiy', 'shiz' IN
! x-, y-, and z-direction.
! The shift is produced by an exponential in Fourier space.

    COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2)
    REAL(DP), INTENT(IN) :: shix
    REAL(DP), INTENT(IN) :: shiy
    REAL(DP), INTENT(IN) :: shiz

    INTEGER :: ind
    REAL(DP) :: dkx, dky, dkz
    COMPLEX(DP), ALLOCATABLE :: q1(:)

    IF (.NOT. ALLOCATED(akv)) STOP "SHIFTSMALL requires FFT"
    ALLOCATE (q1(kdfull2))

    CALL fftf(q0, q1)

    dkx = pi/(dx*REAL(nx, DP))
    dky = pi/(dy*REAL(ny, DP))
    dkz = pi/(dz*REAL(nz, DP))

    DO ind = 1, kdfull2
      q1(ind) = q1(ind)*EXP((shix*akx(ind) + shiy*aky(ind) + shiz*akz(ind)))
    END DO

    CALL fftback(q1, q0)

    DEALLOCATE (q1)

    RETURN
  END SUBROUTINE shiftsmall
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE dipole_qp(q0, shix, shiy, shiz, bostx, bosty, bostz)
!------------------------------------------------------------
    USE params
    USE kinetic
    IMPLICIT NONE

! Initialization of a dipole deformation in space and momentum.
! Shifts COMPLEX wavefunction array 'q0' by 'shix', 'shiy', 'shiz' IN
! x-, y-, and z-direction and boosts it by 'bostx', 'bosty' and 'bostz'.

    INTEGER :: ind, ix, iy, iz, nb
    REAL(DP) :: arg, x1, y1, z1
    COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(IN) :: shix, shiy, shiz, bostx, bosty, bostz

    LOGICAL, PARAMETER :: ttest = .true.
    REAL(DP), ALLOCATABLE :: rhotest(:)

    IF (ttest) THEN
      ALLOCATE (rhotest(kdfull2))
      CALL calcrho(rhotest, q0)
      WRITE (*, *) ' test DIPOLE_QP'
      CALL prifld(rhotest, 'rho before ')
    END IF

    DO nb = 1, nstate

! shift

      CALL shiftsmall(q0(1, nb), shix, shiy, shiz)

! boost

      ind = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz
        DO iy = miny, maxy
          y1 = (iy - nysh)*dy
          DO ix = minx, maxx
            x1 = (ix - nxsh)*dx
            ind = ind + 1
            arg = x1*bostx + y1*bosty + z1*bostz
            q0(ind, nb) = CMPLX(COS(arg), SIN(arg), DP)*q0(ind, nb)
          END DO
        END DO
      END DO

    END DO

    IF (ttest) THEN
      CALL calcrho(rhotest, q0)
      CALL prifld(rhotest, 'rho after ')
      DEALLOCATE (rhotest)
    END IF

    RETURN

  END SUBROUTINE dipole_qp

  SUBROUTINE testcurrent(wfin, rho, iteract)

! ******************************

! Compute currents and prints them on standard output.
!

    USE params
    USE kinetic
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: wfin(kdfull2, kstate)
    INTEGER, INTENT(IN) :: iteract
    REAL(DP), INTENT(IN) :: rho(2*kdfull2)

    REAL(DP), ALLOCATABLE :: curr0(:, :)
    REAL(DP) :: ekincoll
    INTEGER :: icomp,ind

    ALLOCATE (curr0(kdfull2, 3))
    CALL calc_current(curr0, wfin)
    IF (iteract == 0) THEN
      CALL prifld(curr0(:, 1), 'j_1')
      CALL prifld(curr0(:, 2), 'j_2')
      CALL prifld(curr0(:, 3), 'j_3')
    END IF
    WRITE (*, '(a,i5,a,3(1pg13.5))') 'step nr.:', iteract, &
      ', integrated current: ', (SUM(curr0(:, icomp)), icomp=1, 3)
    WRITE (*, '(a,i5,a,3(1pg13.5))') 'step nr.:', iteract, &
      ', integrated |current|:', (SUM(ABS(curr0(:, icomp))), icomp=1, 3)

    ekincoll = 0D0
    DO ind=1,kdfull2
      ekincoll = ekincoll + &
                 (curr0(ind,1)**2+curr0(ind,2)**2+curr0(ind,3)**2)/&
                 MAX(rho(ind),1D-20)
    END DO
    ekincoll = ekincoll*dvol*h2m
    WRITE(*,'(a,1pg13.5)')  'step nr.:', iteract, &
      ', collective kinetic energy=',ekincoll
    WRITE(7,'(a,1pg13.5)')  'step nr.:', iteract, &
      ', collective kinetic energy=',ekincoll
    DEALLOCATE (curr0)

  END SUBROUTINE testcurrent

! ******************************

  SUBROUTINE testgradient(wfin)

! ******************************

! Tests the routines for gradients in x, y, and z direction.
! The test-field 'wfin' is a COMPLEX field.
!

    USE params
    USE kinetic
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN OUT) :: wfin(kdfull2)

    COMPLEX(DP), ALLOCATABLE :: wftest(:)

    REAL(DP) :: ekintestx, ekintesty, ekintestz, ekintot

    IF (.NOT. ALLOCATED(akv)) STOP "TESTGRADIENT requires FFT"
    ALLOCATE (wftest(kdfull2))

    CALL xgradient_rspace(wfin, wftest)
    CALL xgradient_rspace(wftest, wftest)
    ekintestx = dvol*SUM(wfin*wftest)

    CALL ygradient_rspace(wfin, wftest)
    CALL ygradient_rspace(wftest, wftest)
    ekintesty = dvol*SUM(wfin*wftest)

    CALL zgradient_rspace(wfin, wftest)
    CALL zgradient_rspace(wftest, wftest)
    ekintestz = dvol*SUM(wfin*wftest)

    CALL fftf(wfin, wftest)
    wftest = akv*wftest
    CALL fftback(wftest, wftest)
    ekintot = dvol*SUM(wfin*wftest)

    WRITE (6, '(a,5(1pg13.5))') ' test gradients:', &
      ekintestx, ekintesty, ekintestz, ekintestx + ekintesty + ekintestz, &
      ekintot
    WRITE (7, '(a,5(1pg13.5))') ' test gradients:', &
      ekintestx, ekintesty, ekintestz, ekintestx + ekintesty + ekintestz, &
      ekintot

    DEALLOCATE (wftest)

  END SUBROUTINE testgradient

  COMPLEX(DP) FUNCTION determinant(a, n, np)
! Determinant of a COMPLEX matrix 'a'.
! This function is merely an INTERFACE for cludcmp.
    USE params, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: np
    COMPLEX(DP), INTENT(IN OUT) :: a(np, np)
    INTEGER, INTENT(IN) :: n

    REAL(DP) :: d
    COMPLEX(DP) :: det
    INTEGER :: indx(n)
    INTEGER :: ierror

    CALL cludcmp(a, n, np, indx, d, det, ierror)
    determinant = det

    RETURN
  END FUNCTION determinant

  SUBROUTINE cludcmp(a, n, np, indx, d, det, ierror)

! Lower upper decomposition for a COMPLEX matrix:
! a matrix to be decomposed, at the end decomposed matrix
! n actual dimension of matrix
! np physical dimension of matrix
! indx array keeping the permutations from pivoting
! d sign of number of row interchanges
! det determinant of 'a'
    USE params, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: np
    COMPLEX(DP), INTENT(IN OUT) :: a(np, np)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: indx(n)
    REAL(DP), INTENT(OUT) :: d
    COMPLEX(DP), INTENT(OUT) :: det
    INTEGER, INTENT(OUT) :: ierror

    INTEGER, PARAMETER :: nmax = 500
    DOUBLE PRECISION, PARAMETER :: tiny = 1.0D-20
    INTEGER :: i, imax, j, k
    COMPLEX(DP) :: csum, cdum
    REAL(DP) :: dum, aamax, vv(nmax)

    IF (np > nmax) STOP 'LUDCMP: too large matrix'
    ierror = 0

    d = 1.d0
    DO i = 1, n
      aamax = 0.d0
      DO j = 1, n
        IF (ABS(a(i, j)) > aamax) aamax = ABS(a(i, j))
      END DO
      IF (aamax == CMPLX(0D0, 0D0, DP)) THEN
        STOP 'singular matrix in cludcmp. Press ENTER to CONTINUE.'
        ierror = 99
      END IF
      vv(i) = 1.d0/aamax
    END DO
    DO j = 1, n
      DO i = 1, j - 1
        csum = a(i, j)
        DO k = 1, i - 1
          csum = csum - a(i, k)*a(k, j)
        END DO
        a(i, j) = csum
      END DO
      aamax = 0.d0
      DO i = j, n
        csum = a(i, j)
        DO k = 1, j - 1
          csum = csum - a(i, k)*a(k, j)
        END DO
        a(i, j) = csum
        dum = vv(i)*ABS(csum)
        IF (dum >= aamax) THEN
          imax = i
          aamax = dum
        END IF
      END DO
      IF (j /= imax) THEN
        DO k = 1, n
          cdum = a(imax, k)
          a(imax, k) = a(j, k)
          a(j, k) = cdum
        END DO
        d = -d
        vv(imax) = vv(j)
      END IF
      indx(j) = imax
      IF (a(j, j) == CMPLX(0D0, 0D0, DP)) a(j, j) = tiny
      IF (j /= n) THEN
        cdum = 1.d0/a(j, j)
        DO i = j + 1, n
          a(i, j) = a(i, j)*cdum
        END DO
      END IF
    END DO

!calculate determinant
    det = CMPLX(d, 0D0, DP)
    DO i = 1, n
      det = det*a(i, i)
    END DO
    RETURN
  END SUBROUTINE cludcmp

!------------------------------------------------------------

  SUBROUTINE getclustergeometry

! Compute and print geometry parameters of ionic configuration.

    USE params
    IMPLICIT NONE
! calculates characteristic observables for describing
! the geometry of the cluster ions

    INTEGER :: i, ico, ion, j
    REAL(DP) :: dist, rr, sumx, sumy, sumz
    REAL(DP) :: cw(ng, 3)

! get center of mass first

    sumx = 0D0
    sumy = 0D0
    sumz = 0D0

    DO i = 1, nion
      sumx = sumx + cx(i)
      sumy = sumy + cy(i)
      sumz = sumz + cz(i)
    END DO

    comx = sumx/nion
    comy = sumy/nion
    comz = sumz/nion

! now use center of mass as origin

    DO i = 1, nion
      cw(i, 1) = cx(i) - comx
      cw(i, 2) = cy(i) - comy
      cw(i, 3) = cz(i) - comz
    END DO

! ! calculate quadrupole tensor of mass distribution

    DO i = 1, 3
      DO j = i, 3
        qtion(i, j) = 0D0
        DO ion = 1, nion
          qtion(i, j) = qtion(i, j) + 3*cw(ion, i)*cw(ion, j)
          IF (i == j) THEN
            rr = cw(ion, 1)**2 + cw(ion, 2)**2 + cw(ion, 3)**2
            qtion(i, j) = qtion(i, j) - rr
          END IF
        END DO
        qtion(i, j) = qtion(i, j)/nion
        qtion(j, i) = qtion(i, j)
      END DO
    END DO

! ! calculate root mean square radius

    rmsion = 0D0

    DO ion = 1, nion
      rmsion = rmsion + cw(ion, 1)**2 + cw(ion, 2)**2 + cw(ion, 3)**2
    END DO

    rmsion = SQRT(rmsion/nion)

! determine maximum distance; allows to see IF
! cluster is dissociating

    dmdistion = 0D0
    DO i = 1, nion - 1
      DO j = i + 1, nion
        dist = 0D0
        DO ico = 1, 3
          dist = dist + (cw(i, ico) - cw(j, ico))**2
        END DO
        IF (dist > dmdistion) dmdistion = dist
      END DO
    END DO

    dmdistion = SQRT(dmdistion)

    RETURN
  END SUBROUTINE getclustergeometry
!------------------------------------------------------------

!------------------------------------------------------------

  SUBROUTINE getelectrongeometry(q0, nbe)

! Calculates geometry observables of electron orbital density
! number nbe.
!
! Input:
! q0 = set of s.p. wavefunctions
! nbe = state which is wanted
! for total electronic density, set nbe=0
! for total spin-up density, set nbe=-1
! for total spin-down density, set nbe=-2

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
    INTEGER, INTENT(IN) :: nbe

    INTEGER :: i, ind, ix, iy, iz, j, k, ncount
    REAL(DP) :: rr, sum, x1, y1, z1
    REAL(DP) :: xic(3)
    REAL(DP), ALLOCATABLE :: q1(:)

    ALLOCATE (q1(kdfull2))

    q1 = 0D0

    IF (nbe == 0) THEN ! total density
      ncount = nelect
      DO i = 1, nstate
        DO ind = 1, kdfull2
          q1(ind) = q1(ind) + q0(ind, i)*CONJG(q0(ind, i))
        END DO
      END DO
    ELSE IF (nbe == -1) THEN ! total spin-up density
      ncount = nelect  - nspdw
      STOP 'analyse.F: to be implemented'
    ELSE IF (nbe == -2) THEN ! total spin-down density
      ncount = nspdw
      STOP 'analyse.F: to be implemented'
    ELSE ! single orbital density
      ncount = 1
      DO ind = 1, kdfull2
        q1(ind) = q1(ind) + q0(ind, nbe)*CONJG(q0(ind, nbe))
      END DO
    END IF

! get center of density

    codx = 0D0
    cody = 0D0
    codz = 0D0

    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx

          ind = ind + 1
          codx = codx + q1(ind)*x1
          cody = cody + q1(ind)*y1
          codz = codz + q1(ind)*z1

        END DO
      END DO
    END DO

    codx = codx*dvol/ncount
    cody = cody*dvol/ncount
    codz = codz*dvol/ncount

! get quadrupole tensor

    DO i = 1, 3
      DO j = i, 3
        qtel(i, j) = 0D0

        ind = 0
        DO iz = minz, maxz
          xic(3) = (iz - nzsh)*dz - codz
          DO iy = miny, maxy
            xic(2) = (iy - nysh)*dy - cody
            DO ix = minx, maxx
              xic(1) = (ix - nxsh)*dx - codx

              ind = ind + 1
              sum = 3*xic(i)*xic(j)
              IF (i == j) THEN
                rr = 0.
                DO k = 1, 3
                  rr = xic(k)**2 + rr
                END DO
                sum = sum - rr
              END IF

              qtel(i, j) = qtel(i, j) + sum*q1(ind)

            END DO
          END DO
        END DO

        qtel(i, j) = qtel(i, j)*dvol/ncount
        qtel(j, i) = qtel(i, j)

      END DO
    END DO

! get root mean square radius

    ind = 0
    rmsel = 0D0
    DO iz = minz, maxz
      xic(3) = (iz - nzsh)*dz - codz
      DO iy = miny, maxy
        xic(2) = (iy - nysh)*dy - cody
        DO ix = minx, maxx
          xic(1) = (ix - nxsh)*dx - codx

          ind = ind + 1
          rr = 0
          DO k = 1, 3
            rr = rr + xic(k)*xic(k)
          END DO
          rmsel = rmsel + q1(ind)*rr

        END DO
      END DO
    END DO

    rmsel = SQRT(rmsel*dvol/ncount)

    DEALLOCATE (q1)

    RETURN
  END SUBROUTINE getelectrongeometry
!------------------------------------------------------------
!------------------------------------------------------------

  SUBROUTINE calcchargdist(field)

! Computes and prints charge radial charge distribution by
! integrating charge up to a given radius.
! The radii of the analyzing spheres are scanned on an equidstant
! mesh.
!
! Input:
! field = 3D array of (charge-)density
!
! Output:
! Array 'chfld' contains integrated charge in equidistant
! sequence of radii and is printed to unit 323.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: field(kdfull2)

    INTEGER :: ii, ind, ix, iy, iz, n
    REAL(DP) :: chtotal, r, rr, x1, y1, z1
    REAL(DP) :: chfld(100)

    chtotal = 0.0D0
    DO ii = 1, INT(nzsh*dz/drcharges)
      chfld(ii) = 0.0D0
    END DO

    ind = 0

    CALL getcm(1, 0, 0)

    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx

          rr = (x1 - rvectmp(1))**2 + (y1 - rvectmp(2))**2 + (z1 - rvectmp(3))**2

          ind = ind + 1

          chtotal = chtotal + field(ind)

          DO ii = 1, INT(nzsh*dz/drcharges)
            r = ii*drcharges
            IF (rr <= r*r) THEN
              chfld(ii) = chfld(ii) + field(ind)
            END IF
          END DO

        END DO
      END DO
    END DO

    DO ii = 1, INT(nzsh*dz/drcharges)
      chfld(ii) = chfld(ii)*dvol
    END DO

    chtotal = chtotal*dvol

    WRITE (323, '(a,f15.5,1pg13.5)') '# time,total charge=', tfs, chtotal
    WRITE (323, '(2(0pf10.3),1pg13.5)') &
      (tfs, n*drcharges, (chfld(n) - chfld(n - 1))/drcharges, n=2, INT(nzsh*dz/drcharges))
    WRITE (323, '(1x)')
    FLUSH (323)

    RETURN
  END SUBROUTINE calcchargdist

  SUBROUTINE sort_energy(eord, isort, nmax)

    USE params, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nmax
    REAL(DP), INTENT(IN OUT) :: eord(nmax) ! energy array to be sorted
    INTEGER, INTENT(OUT) :: isort(nmax) ! pointer to position before sorting

    REAL(DP) :: emin
    INTEGER :: n, i, isav, imin

! standard ordering before sorting
    DO n = 1, nmax
      isort(n) = n
    END DO

! sort
    DO n = 1, nmax
      emin = 1D32
      imin = n
      DO i = n, nmax
        IF (eord(i) < emin) THEN
          emin = eord(i)
          imin = i
        END IF
      END DO
      isav = isort(imin)
      eord(imin) = eord(n)
      isort(imin) = isort(n)
      eord(n) = emin
      isort(n) = isav
    END DO

  END SUBROUTINE sort_energy

  PURE LOGICAL FUNCTION tmodbuf(it,irep)
!
! This functions evaluates the modulo function only if the
! the denomintor non zero.
!
  INTEGER, INTENT(IN) :: it, irep

    tmodbuf = irep>0 .AND. MOD(it,MAX(irep,1))==0

  END FUNCTION tmodbuf

END MODULE util
