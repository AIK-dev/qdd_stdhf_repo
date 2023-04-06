!
!------falr Coulomb solver-------------------------------------------------

MODULE coulsolv_f

! Package for the FALR Coulomb solver (G. Lauritsch et al,
! Int. J. Mod. Phys. C 5 (1994) 65).
! A model density which reproduces the first multipole moments
! of the given charge density is subtracted. The Coulomb problem
! for the remaining density yields a short-range Coulomb potential
! which is obtained by FFT forward and backward.
! The FALR method is faster THEN the exact Coulomb solver on a
! doubled grid (see 'coulex.F90'). But looses precision if too
! much density comes close to  the bounds of the box.

  USE FFTW
  USE, INTRINSIC :: iso_c_binding
  USE params, ONLY: DP, PI
  IMPLICIT NONE

  PUBLIC
  INTEGER, PRIVATE :: kxmax, kymax, kzmax, ksmax
  INTEGER, PRIVATE :: kdfull, kdred, kddoub
  INTEGER, PRIVATE :: kmom = 25
  INTEGER, PRIVATE :: kfft, kfftx, kffty, kfftz
  INTEGER, PRIVATE :: kdcorf, kfft2
  INTEGER, PARAMETER, PRIVATE :: kqlabel = 7 ! nr. of quadrants minus 1
  REAL(DP), PARAMETER, PRIVATE :: zero = 0D0

  REAL(DP), PRIVATE, ALLOCATABLE :: akv2(:)
  INTEGER, PRIVATE, ALLOCATABLE :: ikm(:, :), indfc(:)
  REAL(DP), PRIVATE :: dx, dy, dz, dxsp, grnorm, fnorm
  REAL(DP), PRIVATE :: dkx, dky, dkz, akmax, dksp
  INTEGER, PRIVATE :: nx, ny, nz, nx1, ny1, nz1, nxr, nxi, nyr, nyi, nzr, nzi, nxy1, nxyz
  INTEGER, PRIVATE :: nxhigh, nxlow, nyhigh, nylow, nzhigh, nzlow
  INTEGER, PRIVATE :: nxk, nxklo, nxkhi, nksp, nkxyz, ecut

  REAL(DP), PRIVATE, ALLOCATABLE :: xval(:), yval(:), zval(:)
  REAL(DP), PRIVATE, ALLOCATABLE :: xt2(:), yt2(:), zt2(:)

! include BLOCK: fields
  REAL(DP), PRIVATE, ALLOCATABLE :: potc(:), rho(:), rhozw(:)
  REAL(DP), PRIVATE, ALLOCATABLE :: rq(:, :)

  REAL(DP), PRIVATE, ALLOCATABLE :: qlcows(:)
  REAL(DP), PRIVATE, ALLOCATABLE, TARGET :: rokall(:, :), pcoall(:, :)
  REAL(DP), PRIVATE :: ax2(0:4), bdiv(0:4), radpow
  INTEGER, PRIVATE :: lpow


  TYPE(C_PTR), PRIVATE :: pforwx, pforwy, pforwz, pbackx, pbacky, pbackz
  INTEGER(C_INT), PRIVATE :: wisdomtest
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, PRIVATE :: fftax(:), fftay(:), fftb(:, :)

CONTAINS

!-------------------------------------------------------------------
  SUBROUTINE init_coul_f(dx0, dy0, dz0, nx0, ny0, nz0)

! Initializes basic arrays.
!
! Input:
! dx0,dy0,dz0 = grid spacings
! nx0,ny0,nz0 = box sizes

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx0, ny0, nz0
    REAL(DP), INTENT(IN):: dx0, dy0, dz0

! compute array parameters
    kxmax = nx0; kymax = ny0; kzmax = nz0
    ksmax = kxmax
! x DIMENSION must be largest
    IF (kymax > kxmax) STOP ' x DIMENSION must be largest'
    IF (kzmax > kxmax) STOP ' x DIMENSION must be largest'

    kdfull = kxmax*kymax*kzmax
    kdred = kdfull
    kddoub = kdfull
    kmom = 25
    kfft = ksmax; kfftx = kxmax; kffty = kymax; kfftz = kzmax
    kdcorf = (kxmax/2 + 1)*(kymax/2 + 1)*(kzmax/2 + 1)
    kfft2 = kfft*2

    ALLOCATE (xval(kxmax), yval(kymax), zval(kzmax))
    ALLOCATE (xt2(kxmax), yt2(kymax), zt2(kzmax))
!ALLOCATE(xrow(kxmax),yrow(kymax),zrow(kzmax))
    ALLOCATE (akv2(kdred))
    ALLOCATE (ikm(kxmax, kymax), indfc(kdfull))
    ALLOCATE (potc(kdfull), rho(kddoub), rhozw(kddoub))
    ALLOCATE (qlcows(kmom), rq(kmom, 0:kqlabel))
    ALLOCATE (rokall(kmom, kdcorf), pcoall(kmom, kdcorf))

! ALLOCATE work arrays for NETLIB

! ALLOCATE work arrays for FFTW3
    ALLOCATE (fftax(kxmax), fftay(kymax), fftb(kzmax, kxmax))

! READ grid parameters from FILE or simply initialize them
    nx = nx0/2
    ny = ny0/2
    nz = nz0/2
    dx = dx0
    dy = dy0
    dz = dz0
    ecut = 100.1D0

! call input routine fftinp, which initializes the grid and fft tabl
    CALL fftinp

! call coucor, which initializes the correction densities and potent
! this has to  be done only once, even if you call coufou (see below)
! more than one time. of course, your grid must stay the same.
    CALL coucor

    RETURN
  END SUBROUTINE init_coul_f

  SUBROUTINE solv_poisson_f(rhoinp, chpfalr, kdf)

! Driver for the solution of the Poisson equation.
!
! Input:
! rhoinp = charge density
! kdf = linear DIMENSION of 3D array
! Output:
! chpfalr = resulting Coulomb potential

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kdf
    REAL(DP), INTENT(IN) :: rhoinp(kdf)
    REAL(DP), INTENT(OUT) :: chpfalr(kdf)

! CALL a routine written by you which writes your density field
! on the array rho.
! remember not to  send your original density array to  the fcs.
! IN this CASE we have a homogeneously charged sphere .
    CALL rhofld(rhoinp, kdf)

! CALL coufou, which contains the fcs procedure.
    CALL coufou2

! CALL a routine written by you which outputs the results of the fcs
! and maybe some other things to  an output file or the screen.
    CALL RESULT(chpfalr, kdf)

  END SUBROUTINE solv_poisson_f


!-----rhofld------------------------------------------------------------

  SUBROUTINE rhofld(rhoinp, kdf)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kdf
    REAL(DP), INTENT(IN) :: rhoinp(kdf)

    rho(1:nxyz) = rhoinp(1:nxyz)

! save the original rho data in rhozw, because rho will be overwritten
    rhozw(1:nxyz) = rho(1:nxyz)

    RETURN
  END SUBROUTINE rhofld

!-----RESULT------------------------------------------------------------

  SUBROUTINE RESULT(chpfalr, kdf)
    USE params, ONLY: e2
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kdf
    REAL(DP), INTENT(OUT) :: chpfalr(kdf)
    INTEGER:: nmax

! inclusion of e2

    nmax = nxi*nyi*nzi
    chpfalr(1:nmax) = e2*potc(1:nmax)

    RETURN
  END SUBROUTINE RESULT

!-----cofows------------------------------------------------------------

  SUBROUTINE coufou2

! The FALR solver. I/O had been prepared in the driver
! routine 'solv_poisson_f' and is communicated via module
! variables.

    IMPLICIT NONE
    INTEGER :: i, i1, i2, i3, icase, ik, ikzero, j
    REAL(DP), ALLOCATABLE :: rhokr(:), rhoki(:)
    INTEGER :: iprho(3)
    REAL(DP):: q00, q10, q11r, q11i, q20, q21r, q21i, &
               q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
               q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2
    REAL(DP)::x00, x10, x11r, x11i, x20, x21r, x21i, &
               x22r, x22i, x30, x31r, x31i, x32r, x32i, x33r, x33i, &
               x40, x41r, x41i, x42r, x42i, x43r, x43i, x44r, x44i, xr2
    REAL(DP) :: ftpi, potcor
    LOGICAL, PARAMETER :: tprint = .false.
    LOGICAL, PARAMETER :: rqplot = .false.

    DATA iprho/1, 1, 1/

!------------------------------------------------------------------------------

    ALLOCATE (rhokr(kdred), rhoki(kdred))

! calculation of the multipole-moments

    CALL mulmws(rho, q00, q10, q11r, q11i, q20, q21r, q21i, &
                q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, tprint)

! effective moments 'rq...' in ordering as fields 'rkoall', 'pcoall'
    IF (ABS(qlcows(1)) > 0D0) THEN
      rq(1, 0) = q00/qlcows(1)
    ELSE
      rq(1, 0) = 0D0
    END IF
    rq(2, 0) = q10/qlcows(2)
    rq(3, 0) = q11r/qlcows(3)
    rq(4, 0) = q11i/qlcows(4)
    rq(5, 0) = q20/qlcows(5)
    rq(6, 0) = q21r/qlcows(6)
    rq(7, 0) = q21i/qlcows(7)
    rq(8, 0) = q22r/qlcows(8)
    rq(9, 0) = q22i/qlcows(9)
    rq(10, 0) = q30/qlcows(10)
    rq(11, 0) = q31r/qlcows(11)
    rq(12, 0) = q31i/qlcows(12)
    rq(13, 0) = q32r/qlcows(13)
    rq(14, 0) = q32i/qlcows(14)
    rq(15, 0) = q33r/qlcows(15)
    rq(16, 0) = q33i/qlcows(16)
    rq(17, 0) = q40/qlcows(17)
    rq(18, 0) = q41r/qlcows(18)
    rq(19, 0) = q41i/qlcows(19)
    rq(20, 0) = q42r/qlcows(20)
    rq(21, 0) = q42i/qlcows(21)
    rq(22, 0) = q43r/qlcows(22)
    rq(23, 0) = q43i/qlcows(23)
    rq(24, 0) = q44r/qlcows(24)
    rq(25, 0) = q44i/qlcows(25)

    DO icase = 1, 7
      rq(:, icase) = rq(:, 0)
    END DO

! set sign for quadrants

    rq(2, 1) = -rq(2, 0)
    rq(6, 1) = -rq(6, 0)
    rq(7, 1) = -rq(7, 0)
    rq(10, 1) = -rq(10, 0)
    rq(13, 1) = -rq(13, 0)
    rq(14, 1) = -rq(14, 0)
    rq(18, 1) = -rq(18, 0)
    rq(19, 1) = -rq(19, 0)
    rq(22, 1) = -rq(22, 0)
    rq(23, 1) = -rq(23, 0)

    rq(4, 2) = -rq(4, 0)
    rq(7, 2) = -rq(7, 0)
    rq(9, 2) = -rq(9, 0)
    rq(12, 2) = -rq(12, 0)
    rq(14, 2) = -rq(14, 0)
    rq(16, 2) = -rq(16, 0)
    rq(19, 2) = -rq(19, 0)
    rq(21, 2) = -rq(21, 0)
    rq(23, 2) = -rq(23, 0)
    rq(25, 2) = -rq(25, 0)

    rq(3, 3) = -rq(3, 0)
    rq(6, 3) = -rq(6, 0)
    rq(9, 3) = -rq(9, 0)
    rq(11, 3) = -rq(11, 0)
    rq(14, 3) = -rq(14, 0)
    rq(15, 3) = -rq(15, 0)
    rq(18, 3) = -rq(18, 0)
    rq(21, 3) = -rq(21, 0)
    rq(22, 3) = -rq(22, 0)
    rq(25, 3) = -rq(25, 0)

    rq(3, 4) = -rq(3, 0)
    rq(4, 4) = -rq(4, 0)
    rq(6, 4) = -rq(6, 0)
    rq(7, 4) = -rq(7, 0)
    rq(11, 4) = -rq(11, 0)
    rq(12, 4) = -rq(12, 0)
    rq(15, 4) = -rq(15, 0)
    rq(16, 4) = -rq(16, 0)
    rq(18, 4) = -rq(18, 0)
    rq(19, 4) = -rq(19, 0)
    rq(22, 4) = -rq(22, 0)
    rq(23, 4) = -rq(23, 0)

    rq(2, 5) = -rq(2, 0)
    rq(3, 5) = -rq(3, 0)
    rq(7, 5) = -rq(7, 0)
    rq(9, 5) = -rq(9, 0)
    rq(10, 5) = -rq(10, 0)
    rq(11, 5) = -rq(11, 0)
    rq(13, 5) = -rq(13, 0)
    rq(15, 5) = -rq(15, 0)
    rq(19, 5) = -rq(19, 0)
    rq(21, 5) = -rq(21, 0)
    rq(23, 5) = -rq(23, 0)
    rq(25, 5) = -rq(25, 0)

    rq(2, 6) = -rq(2, 0)
    rq(4, 6) = -rq(4, 0)
    rq(6, 6) = -rq(6, 0)
    rq(9, 6) = -rq(9, 0)
    rq(10, 6) = -rq(10, 0)
    rq(12, 6) = -rq(12, 0)
    rq(13, 6) = -rq(13, 0)
    rq(16, 6) = -rq(16, 0)
    rq(18, 6) = -rq(18, 0)
    rq(21, 6) = -rq(21, 0)
    rq(22, 6) = -rq(22, 0)
    rq(25, 6) = -rq(25, 0)

    rq(2, 7) = -rq(2, 0)
    rq(3, 7) = -rq(3, 0)
    rq(4, 7) = -rq(4, 0)
    rq(10, 7) = -rq(10, 0)
    rq(11, 7) = -rq(11, 0)
    rq(12, 7) = -rq(12, 0)
    rq(13, 7) = -rq(13, 0)
    rq(14, 7) = -rq(14, 0)
    rq(15, 7) = -rq(15, 0)
    rq(16, 7) = -rq(16, 0)

    IF (rqplot) THEN
      WRITE (7, '(a,e12.4)') ' rq(1,0) =', rq(1, 0)
      WRITE (7, '(3(a,e12.4))') ' rq(2,0) =', rq(2, 0),&
        ' rq(3,0)=', rq(3, 0), ' rq(4,0)=', rq(4, 0)
      WRITE (7, '(a,e12.4/4(a,e12.4))') &
        ' rq(5,0) =', rq(5, 0), ' rq(6,0)=', rq(6, 0), ' rq(7,0)=', rq(7, 0), &
        ' rq(8,0)=', rq(8, 0), ' rq(9,0)=', rq(9, 0)
      WRITE (7, '(3(a,e12.4)/4(a,e12.4))') &
        ' rq(10,0) =', rq(10, 0), ' rq(11,0)=', rq(11, 0), &
        ' rq(12,0)=', rq(12, 0), &
        ' rq(13,0)=', rq(13, 0), ' rq(14,0)=', rq(14, 0), &
        ' rq(15,0)=', rq(15, 0), ' rq(16,0)=', rq(16, 0)
      WRITE (7, '(a,e12.4,4(a,e12.4)/4(a,e12.4))') ' rq(17,0) =', rq(17, 0), &
        ' rq(18,0)=', rq(18, 0), ' rq(19,0)=', rq(19, 0), &
        ' rq(20,0)=', rq(20, 0), ' rq(21,0)=', rq(21, 0), &
        ' rq(22,0)=', rq(22, 0), ' rq(23,0)=', rq(21, 0), &
        ' rq(24,0)=', rq(24, 0), ' rq(25,0)=', rq(25, 0)
    END IF

    i = 0
    DO i3 = 0, nz
      DO i2 = 0, ny
        DO i1 = 0, nx
          i = i + 1
!+++
          j = (i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (i1 + nx)
          rho(j) = rho(j) - SUM(rq(:, 0)*rokall(:, i))
!++-
          IF (i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (i1 + nx)
            rho(j) = rho(j) - SUM(rq(:, 1)*rokall(:, i))
          END IF
!+-+
          IF (i2 /= 0 .AND. i2 /= ny) THEN
            j = (i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (i1 + nx)
            rho(j) = rho(j) - SUM(rq(:, 2)*rokall(:, i))
          END IF
!-++
          IF (i1 /= 0 .AND. i1 /= nx) THEN
            j = (i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (-i1 + nx)
            rho(j) = rho(j) - SUM(rq(:, 3)*rokall(:, i))
          END IF
!--+
          IF (i1 /= 0 .AND. i1 /= nx .AND. i2 /= 0 .AND. i2 /= ny) THEN
            j = (i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (-i1 + nx)
            rho(j) = rho(j) - SUM(rq(:, 4)*rokall(:, i))
          END IF
!-+-
          IF (i1 /= 0 .AND. i1 /= nx .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (-i1 + nx)
            rho(j) = rho(j) - SUM(rq(:, 5)*rokall(:, i))
          END IF
!+--
          IF (i2 /= 0 .AND. i2 /= ny .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (i1 + nx)
            rho(j) = rho(j) - SUM(rq(:, 6)*rokall(:, i))
          END IF
!---
          IF (i1 /= 0 .AND. i1 /= nx .AND. i2 /= 0 .AND. i2 /= ny &
              .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (-i1 + nx)
            rho(j) = rho(j) - SUM(rq(:, 7)*rokall(:, i))
          END IF
        END DO
      END DO
    END DO

    CALL mulmws(rho, x00, x10, x11r, x11i, x20, x21r, x21i, &
                x22r, x22i, x30, x31r, x31i, x32r, x32i, x33r, x33i, &
                x40, x41r, x41i, x42r, x42i, x43r, x43i, x44r, x44i, xr2, tprint)

! Fourier transformation of the density
    CALL fourf(rho, rhokr, rhoki, iprho)

! calculation of the Coulomb field (writing on the density field)

    ftpi = 4D0*pi
    ikzero = nxy1*(nz - 1) + nxi*(ny - 1) + nx
    ikzero = indfc(ikzero)
    DO ik = 1, nksp
      IF (ik /= ikzero) THEN
        rhokr(ik) = ftpi/akv2(ik)*rhokr(ik)
        rhoki(ik) = ftpi/akv2(ik)*rhoki(ik)
      ELSE
! special treatment of the singularity at k=0:
! because of the special norm of the Fourier fields one has to
! divide the analytic term by sqrt(2*pi)**3
        rhokr(ik) = -qr2/SQRT(8D0)/(3D0*SQRT(2D0*pi))
        rhoki(ik) = 0D0
      END IF
    END DO

! Fourier back transformation

    CALL fourb(potc, rhokr, rhoki, iprho)
    potcor = -potc(nxyz)
    DO i = 1, nxyz
      potc(i) = potcor + potc(i)
    END DO

! addition of the correction field

    i = 0
    DO i3 = 0, nz
      DO i2 = 0, ny
        DO i1 = 0, nx
          i = i + 1
!+++
          j = (i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (i1 + nx)
          potc(j) = potc(j) + SUM(rq(:, 0)*pcoall(:, i))
!++-CASE 1
          IF (i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (i1 + nx)
            potc(j) = potc(j) + SUM(rq(:, 1)*pcoall(:, i))
          END IF
!+-+CASE 2
          IF (i2 /= 0 .AND. i2 /= ny) THEN
            j = (i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (i1 + nx)
            potc(j) = potc(j) + SUM(rq(:, 2)*pcoall(:, i))
          END IF
!-++CASE 3
          IF (i1 /= 0 .AND. i1 /= nx) THEN
            j = (i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (-i1 + nx)
            potc(j) = potc(j) + SUM(rq(:, 3)*pcoall(:, i))
          END IF
!--+CASE 4
          IF (i1 /= 0 .AND. i1 /= nx .AND. i2 /= 0 .AND. i2 /= ny) THEN
            j = (i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (-i1 + nx)
            potc(j) = potc(j) + SUM(rq(:, 4)*pcoall(:, i))
          END IF
!-+-CASE 5
          IF (i1 /= 0 .AND. i1 /= nx .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (-i1 + nx)
            potc(j) = potc(j) + SUM(rq(:, 5)*pcoall(:, i))
          END IF
!+--CASE 6
          IF (i2 /= 0 .AND. i2 /= ny .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (i1 + nx)
            potc(j) = potc(j) + SUM(rq(:, 6)*pcoall(:, i))
          END IF
!---CASE 7
          IF (i1 /= 0 .AND. i1 /= nx .AND. i2 /= 0 .AND. i2 /= ny &
              .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (-i1 + nx)
            potc(j) = potc(j) + SUM(rq(:, 7)*pcoall(:, i))
          END IF
        END DO
      END DO
    END DO

    DEALLOCATE (rhokr, rhoki)

    RETURN

  END SUBROUTINE coufou2

!-----fx1-------------------------------------------------------------------fx1

  REAL(DP) FUNCTION fx1(x)

! Model function for FALR solver.

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: x

    REAL(DP), PARAMETER :: eighty = 80.0D0
    REAL(DP) :: expa, x2
!------------------------------------------------------------------------------

    x2 = x*x
    expa = x2*ax2(lpow)
    IF (expa <= eighty) THEN
      expa = EXP(-expa)
    ELSE
      expa = zero
    END IF

    IF (lpow > 0) THEN
      IF (radpow == zero) THEN
        fx1 = expa/(bdiv(lpow) + x**lpow)
      ELSE IF (x /= zero) THEN
        fx1 = x**radpow*expa/(bdiv(lpow) + x**lpow)
      ELSE
        fx1 = zero
      END IF
    ELSE
      IF (radpow == zero) THEN
        fx1 = expa
      ELSE IF (x /= zero) THEN
        fx1 = x**radpow*expa
      ELSE
        fx1 = zero
      END IF
    END IF

    RETURN
  END FUNCTION fx1

!-----cocows-------------------------------------------------------------cocows

  SUBROUTINE coucor

! initializes the multipole correction fields

! these correction fields have the following parities:
! lab field x y z
! 1 00 + + +
! 2 10 + + -
! 3 11r - + +
! 4 11i + - +
! 5 20 + + +
! 6 21r - + -
! 7 21i + - -
! 8 22r + + +
! 9 22i - - +
! 10 30 + + -
! 11 31r - + +
! 12 31i + - +
! 13 32r + + -
! 14 32i - - -
! 15 33r - + +
! 16 33i + - +
! 17 40 + + +
! 18 41r - + -
! 19 41i + - -
! 20 42r + + +
! 21 42i - - +
! 22 43r - + -
! 23 43i + - -
! 24 44r + + +
! 25 44i - - +

    IMPLICIT NONE

    LOGICAL, PARAMETER :: mumpri = .false.
    LOGICAL, PARAMETER :: refmom = .true.

    INTEGER :: i, i1, i2, i3, l
    REAL(DP) :: one, two, xmax, ymax, zmax, xyzmin, xyzmax
    REAL(DP) :: agen, binup, rad0, fzrm, qg, qgint, qr2, rad, gfc
    REAL(DP) :: pfy00, pfy10, pfy11, pfy20, pfy21, pfy22, pfy30, pfy31, &
                pfy32, pfy33, pfy40, pfy41, pfy42, pfy43, pfy44
    REAL(DP) :: y00, y10, y11r, y11i, y20, y21r, y21i, y22r, y22i, y30, y31r, y31i, &
                y32r, y32i, y33r, y33i, y40, y41r, y41i, y42r, y42i, y43r, y43i, y44r, y44i
    REAL(DP) :: q00, q10, q11r, q11i, q20, q21r, q21i, q22r, q22i, q30, q31r, q31i, &
                q32r, q32i, q33r, q33i, q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i
    REAL(DP) :: t1, t2
    REAL(DP) :: x, y, z, rr, xx, yy, zz

!------------------------------------------------------------------------------

    xmax = nx1*dx
    ymax = ny1*dy
    zmax = nz1*dz
    xyzmin = MIN(xmax, ymax, zmax)
    xyzmax = MAX(xmax, ymax, zmax)
    binup = SQRT(xmax*xmax + ymax*ymax + zmax*zmax)
    rad0 = xyzmax
    fzrm = xyzmin/5.5D0

    agen = 6D0*SQRT(LOG(2D0))/xyzmin
    ax2(0) = agen*agen
    bdiv(0) = 1D0
    agen = SQRT(LOG(2D0))/fzrm
    agen = agen*agen
    DO l = 1, 4
      ax2(l) = agen
      bdiv(l) = rad0**l
    END DO

    pfy00 = 1D0/SQRT(4D0*pi)
    pfy10 = SQRT(3D0/(4D0*pi))
    pfy11 = SQRT(3D0/(4D0*pi))
    pfy20 = SQRT(5D0/(16D0*pi))
    pfy21 = SQRT(15D0/(4D0*pi))
    pfy22 = SQRT(15D0/(16D0*pi))
    pfy30 = SQRT(7D0/(16D0*pi))
    pfy31 = SQRT(21D0/(32D0*pi))
    pfy32 = SQRT(105D0/(16D0*pi))
    pfy33 = SQRT(35D0/(32D0*pi))
    pfy40 = 3D0/(16D0*SQRT(pi))
    pfy41 = SQRT(45D0/(32D0*pi))
    pfy42 = (3D0/8D0)*SQRT(5D0/(1D0*pi))
    pfy43 = (3D0/4D0)*SQRT(35D0/(2D0*pi))
    pfy44 = (3D0/8D0)*SQRT(35D0/(4D0*pi))

    DO l = 0, 4

      WRITE (6, '(a,i2,a)') 'l =', l, ':'
      WRITE (7, '(a,i2,a)') 'l =', l, ':'

      lpow = l

      radpow = 2D0*l + 2D0
! CALL qgaus(fx1,zero,binup,qgint)
      CALL qgaus_fx1(zero, binup, qgint)
      IF (qgint /= 0D0) THEN
        qg = 1/qgint
      ELSE
        qg = 0D0
      END IF

      i = 0
      DO i3 = nz, nzi
        DO i2 = ny, nyi
          DO i1 = nx, nxi
            i = i + 1

            x = xval(i1)
            y = yval(i2)
            z = zval(i3)
            xx = xt2(i1)
            yy = yt2(i2)
            zz = zt2(i3)

            rr = xx + yy + zz
            rad = SQRT(rr)

            radpow = zero
            gfc = qg*fx1(rad)

            IF (l == 0) THEN
              y00 = pfy00
              rokall(1, i) = gfc*y00
            ELSE IF (l == 1) THEN
              y10 = pfy10*z
              y11r = pfy11*(-x)
              y11i = pfy11*(-y)
              rokall(2, i) = gfc*y10
              rokall(3, i) = gfc*y11r
              rokall(4, i) = gfc*y11i
            ELSE IF (l == 2) THEN
              y20 = pfy20*(3D0*zz - rr)
              y21r = pfy21*(-x*z)
              y21i = pfy21*(-y*z)
              y22r = pfy22*(xx - yy)
              y22i = pfy22*2D0*x*y
              rokall(5, i) = gfc*y20
              rokall(6, i) = gfc*y21r
              rokall(7, i) = gfc*y21i
              rokall(8, i) = gfc*y22r
              rokall(9, i) = gfc*y22i
            ELSE IF (l == 3) THEN
              y30 = pfy30*(5D0*zz - 3D0*rr)*z
              y31r = pfy31*(-5D0*zz + rr)*x
              y31i = pfy31*(-5D0*zz + rr)*y
              y32r = pfy32*(xx - yy)*z
              y32i = pfy32*2D0*x*y*z
              y33r = pfy33*(-xx + 3D0*yy)*x
              y33i = pfy33*(-3D0*xx + yy)*y
              rokall(10, i) = gfc*y30
              rokall(11, i) = gfc*y31r
              rokall(12, i) = gfc*y31i
              rokall(13, i) = gfc*y32r
              rokall(14, i) = gfc*y32i
              rokall(15, i) = gfc*y33r
              rokall(16, i) = gfc*y33i
            ELSE IF (l == 4) THEN
              y40 = pfy40*(35D0*zz*zz - 30D0*rr*zz + 3D0*rr*rr)
              y41r = pfy41*(-7D0*zz + 3D0*rr)*x*z
              y41i = pfy41*(-7D0*zz + 3D0*rr)*y*z
              y42r = pfy42*(xx - yy)*(7D0*zz - rr)
              y42i = pfy42*2D0*x*y*(7D0*zz - rr)
              y43r = pfy43*(-xx + 3D0*yy)*x*z
              y43i = pfy43*(-3D0*xx + yy)*y*z
              y44r = pfy44*(xx*xx - 6D0*xx*yy + yy*yy)
              y44i = pfy44*(xx - yy)*4D0*x*y
              rokall(17, i) = gfc*y40
              rokall(18, i) = gfc*y41r
              rokall(19, i) = gfc*y41i
              rokall(20, i) = gfc*y42r
              rokall(21, i) = gfc*y42i
              rokall(22, i) = gfc*y43r
              rokall(23, i) = gfc*y43i
              rokall(24, i) = gfc*y44r
              rokall(25, i) = gfc*y44i
            END IF

            IF (rad < 1.d-2) THEN
              t1 = zero
            ELSE
              radpow = 2D0*l + 2D0
! CALL qgaus(fx1,zero,rad,t1)
              CALL qgaus_fx1(zero, rad, t1)
              t1 = t1/(rad**(2*l + 1))
            END IF
            radpow = 1D0
            CALL qgaus_fx1(rad, binup, t2)
            gfc = 4D0*pi/(2D0*l + 1D0)*qg*(t1 + t2)

            IF (l == 0) THEN
              pcoall(1, i) = gfc*y00
            ELSE IF (l == 1) THEN
              pcoall(2, i) = gfc*y10
              pcoall(3, i) = gfc*y11r
              pcoall(4, i) = gfc*y11i
            ELSE IF (l == 2) THEN
              pcoall(5, i) = gfc*y20
              pcoall(6, i) = gfc*y21r
              pcoall(7, i) = gfc*y21i
              pcoall(8, i) = gfc*y22r
              pcoall(9, i) = gfc*y22i
            ELSE IF (l == 3) THEN
              pcoall(10, i) = gfc*y30
              pcoall(11, i) = gfc*y31r
              pcoall(12, i) = gfc*y31i
              pcoall(13, i) = gfc*y32r
              pcoall(14, i) = gfc*y32i
              pcoall(15, i) = gfc*y33r
              pcoall(16, i) = gfc*y33i
            ELSE IF (l == 4) THEN
              pcoall(17, i) = gfc*y40
              pcoall(18, i) = gfc*y41r
              pcoall(19, i) = gfc*y41i
              pcoall(20, i) = gfc*y42r
              pcoall(21, i) = gfc*y42i
              pcoall(22, i) = gfc*y43r
              pcoall(23, i) = gfc*y43i
              pcoall(24, i) = gfc*y44r
              pcoall(25, i) = gfc*y44i
            END IF

          END DO
        END DO
      END DO

    END DO

    IF (refmom) THEN

      WRITE (7, '(a)') 'y00 density:'
      CALL expand(rokall(1, :), potc, 1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(1) = q00

      WRITE (7, '(a)') 'y10 density:'
      CALL expand(rokall(2, :), potc, 1, 1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(2) = q10

      WRITE (7, '(a)') 'y11r density:'
      CALL expand(rokall(3, :), potc, -1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(3) = q11r

      WRITE (7, '(a)') 'y11i density:'
      CALL expand(rokall(4, :), potc, 1, -1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(4) = q11i

      WRITE (7, '(a)') 'y20 density:'
      CALL expand(rokall(5, :), potc, 1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(5) = q20

      WRITE (7, '(a)') 'y21r density:'
      CALL expand(rokall(6, :), potc, -1, 1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(6) = q21r

      WRITE (7, '(a)') 'y21i density:'
      CALL expand(rokall(7, :), potc, 1, -1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(7) = q21i

      WRITE (7, '(a)') 'y22r density:'
      CALL expand(rokall(8, :), potc, 1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(8) = q22r

      WRITE (7, '(a)') 'y22i density:'
      CALL expand(rokall(9, :), potc, -1, -1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(9) = q22i

      WRITE (7, '(a)') 'y30 density:'
      CALL expand(rokall(10, :), potc, 1, 1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(10) = q30

      WRITE (7, '(a)') 'y31r density:'
      CALL expand(rokall(11, :), potc, -1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(11) = q31r

      WRITE (7, '(a)') 'y31i density:'
      CALL expand(rokall(12, :), potc, 1, -1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(12) = q31i

      WRITE (7, '(a)') 'y32r density:'
      CALL expand(rokall(13, :), potc, 1, 1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(13) = q32r

      WRITE (7, '(a)') 'y32i density:'
      CALL expand(rokall(14, :), potc, -1, -1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(14) = q32i

      WRITE (7, '(a)') 'y33r density:'
      CALL expand(rokall(15, :), potc, -1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(15) = q33r

      WRITE (7, '(a)') 'y33i density:'
      CALL expand(rokall(16, :), potc, 1, -1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(16) = q33i

      WRITE (7, '(a)') 'y40 density:'
      CALL expand(rokall(17, :), potc, 1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(17) = q40

      WRITE (7, '(a)') 'y41r density:'
      CALL expand(rokall(18, :), potc, -1, 1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(18) = q41r

      WRITE (7, '(a)') 'y41i density:'
      CALL expand(rokall(19, :), potc, 1, -1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(19) = q41i

      WRITE (7, '(a)') 'y42r density:'
      CALL expand(rokall(20, :), potc, 1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(20) = q42r

      WRITE (7, '(a)') 'y42i density:'
      CALL expand(rokall(21, :), potc, -1, -1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(21) = q42i

      WRITE (7, '(a)') 'y43r density:'
      CALL expand(rokall(22, :), potc, -1, 1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(22) = q43r

      WRITE (7, '(a)') 'y43i density:'
      CALL expand(rokall(23, :), potc, 1, -1, -1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(23) = q43i

      WRITE (7, '(a)') 'y44r density:'
      CALL expand(rokall(24, :), potc, 1, 1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(24) = q44r

      WRITE (7, '(a)') 'y44i density:'
      CALL expand(rokall(25, :), potc, -1, -1, 1)
      CALL mulmws(potc, q00, q10, q11r, q11i, q20, q21r, q21i, &
                  q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                  q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)
      qlcows(25) = q44i

    ELSE

      qlcows(1) = 1D0/pfy00
      qlcows(2) = 1D0/pfy10
      qlcows(3) = 1D0/pfy11
      qlcows(4) = 1D0/pfy11
      qlcows(5) = 1D0/pfy20
      qlcows(6) = 1D0/pfy21
      qlcows(7) = 1D0/pfy21
      qlcows(8) = 1D0/pfy22
      qlcows(9) = 1D0/pfy22
      qlcows(10) = 1D0/pfy30
      qlcows(11) = 1D0/pfy31
      qlcows(12) = 1D0/pfy31
      qlcows(13) = 1D0/pfy32
      qlcows(14) = 1D0/pfy32
      qlcows(15) = 1D0/pfy33
      qlcows(16) = 1D0/pfy33
      qlcows(17) = 1D0/pfy40
      qlcows(18) = 1D0/pfy41
      qlcows(19) = 1D0/pfy41
      qlcows(20) = 1D0/pfy42
      qlcows(21) = 1D0/pfy42
      qlcows(22) = 1D0/pfy43
      qlcows(23) = 1D0/pfy43
      qlcows(24) = 1D0/pfy44
      qlcows(25) = 1D0/pfy44

    END IF

    WRITE (7, '(a,10(/1x,5g12.4))') 'qlcor:', qlcows

    RETURN
  END SUBROUTINE coucor
!-----tstbnd-------------------------------------------------------------tstbnd

  SUBROUTINE tstbnd(rhotst)

! tests boundary values of 'rhotst'

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rhotst(kdfull)
    INTEGER :: i1, i2, i3, ii, numb
    REAL(DP):: bndmax, bndmin, rhp

!------------------------------------------------------------------------------

! preoccupation

    bndmax = -1000D0
    bndmin = +1000D0

    ii = 0
    numb = 0
    DO i3 = 1, nzi
      DO i2 = 1, nyi
        DO i1 = 1, nxi
          ii = ii + 1
          IF (i3 == 1 .OR. i3 == nzi .OR. i2 == 1 .OR. i2 == nyi &
              .OR. i1 == 1 .OR. i1 == nxi) THEN
            rhp = rhotst(ii)
            bndmax = MAX(bndmax, rhp)
            bndmin = MIN(bndmin, rhp)
            numb = 1 + numb
          END IF
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE tstbnd

!-----mulmws-------------------------------------------------------------mulmws

  SUBROUTINE mulmws(rhomul, q00, q10, q11r, q11i, q20, q21r, q21i, &
                    q22r, q22i, q30, q31r, q31i, q32r, q32i, q33r, q33i, &
                    q40, q41r, q41i, q42r, q42i, q43r, q43i, q44r, q44i, qr2, mumpri)

! calculates the multipolemoments of rhomul

! output: q00 (monopole)
! q10,q11r,q11i (dipole)
! q20,q2#r,q2#i (quadrupole)
! q30,q3#r,q3#i (octupole)
! q40,q4#r,q4#i (hexadecupole)
! qr2 (second moment)

! input : rhomul (spatial array of density)
! mumpri (PRINT option)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rhomul(kdfull)
    REAL(DP), INTENT(OUT) :: q00
    REAL(DP), INTENT(OUT) :: q10
    REAL(DP), INTENT(OUT) :: q11r
    REAL(DP), INTENT(OUT) :: q11i
    REAL(DP), INTENT(OUT) :: q20
    REAL(DP), INTENT(OUT) :: q21r
    REAL(DP), INTENT(OUT) :: q21i
    REAL(DP), INTENT(OUT) :: q22r
    REAL(DP), INTENT(OUT) :: q22i
    REAL(DP), INTENT(OUT) :: q30
    REAL(DP), INTENT(OUT) :: q31r
    REAL(DP), INTENT(OUT) :: q31i
    REAL(DP), INTENT(OUT) :: q32r
    REAL(DP), INTENT(OUT) :: q32i
    REAL(DP), INTENT(OUT) :: q33r
    REAL(DP), INTENT(OUT) :: q33i
    REAL(DP), INTENT(OUT) :: q40
    REAL(DP), INTENT(OUT) :: q41r
    REAL(DP), INTENT(OUT) :: q41i
    REAL(DP), INTENT(OUT) :: q42r
    REAL(DP), INTENT(OUT) :: q42i
    REAL(DP), INTENT(OUT) :: q43r
    REAL(DP), INTENT(OUT) :: q43i
    REAL(DP), INTENT(OUT) :: q44r
    REAL(DP), INTENT(OUT) :: q44i
    REAL(DP), INTENT(OUT) :: qr2
    LOGICAL, INTENT(IN) :: mumpri

    INTEGER :: i1, i2, i3, ii
    REAL(DP) ::rhp, rr, x, y, z, xx, yy, zz

!------------------------------------------------------------------------------

    IF (mumpri) WRITE (7, '(a)') ' begin mulmws'

! preoccupation

    q00 = 0D0
    q10 = 0D0
    q11r = 0D0
    q11i = 0D0
    q20 = 0D0
    q21r = 0D0
    q21i = 0D0
    q22r = 0D0
    q22i = 0D0
    q30 = 0D0
    q31r = 0D0
    q31i = 0D0
    q32r = 0D0
    q32i = 0D0
    q33r = 0D0
    q33i = 0D0
    q40 = 0D0
    q41r = 0D0
    q41i = 0D0
    q42r = 0D0
    q42i = 0D0
    q43r = 0D0
    q43i = 0D0
    q44r = 0D0
    q44i = 0D0
    qr2 = 0D0

    ii = 0
    DO i3 = 1, nzi
      DO i2 = 1, nyi
        DO i1 = 1, nxi
          ii = ii + 1

          x = xval(i1)
          y = yval(i2)
          z = zval(i3)
          xx = xt2(i1)
          yy = yt2(i2)
          zz = zt2(i3)

          rr = xx + yy + zz
          rhp = rhomul(ii)

          qr2 = qr2 + rhp*rr

          q00 = q00 + rhp

          q10 = q10 + rhp*z
          q11r = q11r + rhp*(-x)
          q11i = q11i + rhp*(-y)

          q20 = q20 + rhp*(3D0*zz - rr)
          q21r = q21r + rhp*(-x*z)
          q21i = q21i + rhp*(-y*z)
          q22r = q22r + rhp*(xx - yy)
          q22i = q22i + rhp*2D0*x*y

          q30 = q30 + rhp*(5D0*zz - 3D0*rr)*z
          q31r = q31r + rhp*(-5D0*zz + rr)*x
          q31i = q31i + rhp*(-5D0*zz + rr)*y
          q32r = q32r + rhp*(xx - yy)*z
          q32i = q32i + rhp*2D0*x*y*z
          q33r = q33r + rhp*(-xx + 3D0*yy)*x
          q33i = q33i + rhp*(-3D0*xx + yy)*y

          q40 = q40 + rhp*(35D0*zz*zz - 30D0*rr*zz + 3D0*rr*rr)
          q41r = q41r + rhp*(-7D0*zz + 3D0*rr)*x*z
          q41i = q41i + rhp*(-7D0*zz + 3D0*rr)*y*z
          q42r = q42r + rhp*(xx - yy)*(7D0*zz - rr)
          q42i = q42i + rhp*2D0*x*y*(7D0*zz - rr)
          q43r = q43r + rhp*(-xx + 3D0*yy)*x*z
          q43i = q43i + rhp*(-3D0*xx + yy)*y*z
          q44r = q44r + rhp*(xx*xx - 6D0*xx*yy + yy*yy)
          q44i = q44i + rhp*(xx - yy)*4D0*x*y

        END DO
      END DO
    END DO

    q00 = q00*dxsp
    q10 = q10*dxsp
    q11r = q11r*dxsp
    q11i = q11i*dxsp
    q20 = q20*dxsp
    q21r = q21r*dxsp
    q21i = q21i*dxsp
    q22r = q22r*dxsp
    q22i = q22i*dxsp
    q30 = q30*dxsp
    q31r = q31r*dxsp
    q31i = q31i*dxsp
    q32r = q32r*dxsp
    q32i = q32i*dxsp
    q33r = q33r*dxsp
    q33i = q33i*dxsp
    q40 = q40*dxsp
    q41r = q41r*dxsp
    q41i = q41i*dxsp
    q42r = q42r*dxsp
    q42i = q42i*dxsp
    q43r = q43r*dxsp
    q43i = q43i*dxsp
    q44r = q44r*dxsp
    q44i = q44i*dxsp
    qr2 = qr2*dxsp

    IF (mumpri) THEN
      WRITE (7, '(2(a,e12.4))') ' qr2 =', qr2, ' q00 =', q00
      WRITE (7, '(3(a,e12.4))') ' q10 =', q10, ' q11r=', q11r, ' q11i=', q11i
      WRITE (7, '(a,e12.4/4(a,e12.4))') ' q20 =', q20, ' q21r=', q21r, ' q21i=', q21i, &
        ' q22r=', q22r, ' q22i=', q22i
      WRITE (7, '(3(a,e12.4)/4(a,e12.4))') &
        ' q30 =', q30, ' q31r=', q31r, ' q31i=', q31i, ' q32r=', q32r, ' q32i=', q32i, &
        ' q33r=', q33r, ' q33i=', q33i
      WRITE (7, '(a,e12.4,4(a,e12.4)/4(a,e12.4))') ' q40 =', q40, &
        ' q41r=', q41r, ' q41i=', q41i, ' q42r=', q42r, ' q42i=', q42i, &
        ' q43r=', q43r, ' q43i=', q43i, ' q44r=', q44r, ' q44i=', q44i
    END IF

    IF (mumpri) THEN
      WRITE (7, '(a)') ' END mulmws'
    END IF

    RETURN

  END SUBROUTINE mulmws

!-----expand-------------------------------------------------------------cofows

  SUBROUTINE expand(rhoin, rhoout, ipx, ipy, ipz)

! Expands a compressed field into full 3d
! rhoin = input field in upper octant
! rhoout = output field in full 3d box
! ipx,y,z = parities IN x,y,z

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rhoin(kdcorf)
    REAL(DP), INTENT(OUT) :: rhoout(kdfull)
    INTEGER, INTENT(IN) :: ipx
    INTEGER, INTENT(IN) :: ipy
    INTEGER, INTENT(IN) :: ipz

!INTEGER :: i,i1,i2,i3,ii,j
    INTEGER :: i, i1, i2, i3, j
!------------------------------------------------------------------------------

    i = 0
    DO i3 = 0, nz
      DO i2 = 0, ny
        DO i1 = 0, nx
          i = i + 1
!+++
          j = (i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (i1 + nx)
          rhoout(j) = rhoin(i)
!++-
          IF (i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (i1 + nx)
            rhoout(j) = ipz*rhoin(i)
          END IF
!+-+
          IF (i2 /= 0 .AND. i2 /= ny) THEN
            j = (i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (i1 + nx)
            rhoout(j) = ipy*rhoin(i)
          END IF
!-++
          IF (i1 /= 0 .AND. i1 /= nx) THEN
            j = (i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (-i1 + nx)
            rhoout(j) = ipx*rhoin(i)
          END IF
!--+
          IF (i1 /= 0 .AND. i1 /= nx .AND. i2 /= 0 .AND. i2 /= ny) THEN
            j = (i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (-i1 + nx)
            rhoout(j) = ipx*ipy*rhoin(i)
          END IF
!-+-
          IF (i1 /= 0 .AND. i1 /= nx .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (i2 + ny - 1)*nxi + (-i1 + nx)
            rhoout(j) = ipx*ipz*rhoin(i)
          END IF
!+--
          IF (i2 /= 0 .AND. i2 /= ny .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (i1 + nx)
            rhoout(j) = ipy*ipz*rhoin(i)
          END IF
!---
          IF (i1 /= 0 .AND. i1 /= nx .AND. i2 /= 0 .AND. i2 /= ny &
              .AND. i3 /= 0 .AND. i3 /= nz) THEN
            j = (-i3 + nz - 1)*nxy1 + (-i2 + ny - 1)*nxi + (-i1 + nx)
            rhoout(j) = ipx*ipy*ipz*rhoin(i)
          END IF
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE expand

!----qgaus-------------------------------------------------------------

! Gauss quadrature adapted from numerical recipes routine qgaus (20 point)
!
  SUBROUTINE qgaus_fx1(a, b, ss)
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a
    REAL(DP), INTENT(IN) :: b
    REAL(DP), INTENT(OUT) :: ss
    REAL(DP) :: x(10), w(10)

    INTEGER :: j
    REAL(DP) :: xm, xr

    DATA x/.076526521133497D0, .227785851141645D0, .373706088715419D0, &
      .510867001950827D0, .636053680726515D0, .746331906460150D0, &
      .839116971822218D0, .912234428251325D0, .963971927277913D0, &
      .993128599185094D0/
    DATA w/.152753387130725D0, .149172986472603D0, .142096109318382D0, &
      .131688638449176D0, .118194531961518D0, .101930119817240D0, &
      .083276741576704D0, .062672036334109D0, .040601429800386D0, &
      .017614007139152D0/

    xm = 0.5D0*(b + a)
    xr = 0.5D0*(b - a)
    ss = 0.0D0
    DO j = 1, 10
      dx = xr*x(j)
      ss = ss + w(j)*(fx1(xm + dx) + fx1(xm - dx))
    END DO
    ss = xr*ss

    RETURN
  END SUBROUTINE qgaus_fx1

!-----fftinp------------------------------------------------------------

  SUBROUTINE fftinp
    IMPLICIT NONE

! Initializes parameters for Coulomb solver (mostly for FFT).

!-----------------------------------------------------------------------

    INTEGER :: i0, i1, i2, i3, ii, ix, iy, iz
    REAL(DP) :: ak2, sqak2
    REAL(DP) :: xx1, xy1, xz1, xx2, xy2, xz2

! initialize grid IN coordinate space
    nx1 = nx + 1
    ny1 = ny + 1
    nz1 = nz + 1
    nxr = nx + nx
    nyr = ny + ny
    nzr = nz + nz
    nxi = nxr
    nyi = nyr
    nzi = nzr
    nxy1 = nxi*nyi
    nxyz = nxi*nyi*nzi

! grid lengths must match with parameters IN incs
    IF (kxmax < nxi) THEN
      WRITE (6, '(a)') 'ERROR: parameter kxmax too small '
      STOP ' error in parameter'
    ELSE IF (kymax < nyi) THEN
      WRITE (6, '(a)') 'ERROR: parameter kymax too small '
      STOP ' error in parameter'
    ELSE IF (kzmax < nzi) THEN
      WRITE (6, '(a)') 'ERROR: parameter kzmax too small '
      STOP ' error in parameter'
    END IF

! calculate "real" coordinates of the gridpoints

    xx1 = -REAL(nx)*dx
    DO ix = 1, nxi
      xx1 = xx1 + dx
      xval(ix) = xx1
      xt2(ix) = xx1*xx1
    END DO

    xy1 = -REAL(ny)*dy
    DO iy = 1, nyi
      xy1 = xy1 + dy
      yval(iy) = xy1
      yt2(iy) = xy1*xy1
    END DO

    xz1 = -REAL(nz)*dz
    DO iz = 1, nzi
      xz1 = xz1 + dz
      zval(iz) = xz1
      zt2(iz) = xz1*xz1
    END DO

    nkxyz = nxi*nyi*nzi

! initialize grid in Fourier space

    dkx = pi/(dx*REAL(nx, DP))
    dky = pi/(dy*REAL(ny, DP))
    dkz = pi/(dz*REAL(nz, DP))

    dxsp = dx*dy*dz
    dksp = dkx*dky*dkz

    grnorm = SQRT(dxsp/dksp)
    fnorm = 1D0/SQRT(REAL(nx*ny*nz, DP))
    akmax = SQRT(3D0*(nx*nx)*dx*dx) + 2D0
    nxk = INT(akmax/dkx) + 1
    IF (nxk > nx1) nxk = nx1

    xy1 = -dky
    DO i2 = 1, ny1
      xy1 = xy1 + dky
      xy2 = xy1*xy1
      xx1 = -dkx
      DO i1 = 1, nx1
        xx1 = xx1 + dkx
        xx2 = xx1*xx1
        ak2 = xx2 + xy2
        sqak2 = SQRT(ak2)
        ikm(nx + i1 - 1, ny + i2 - 1) = 0
        IF ((i1 < nx1) .AND. (i2 < ny1)) THEN
          ikm(nx + i1 - 1, ny - i2 + 1) = 0
          ikm(nx - i1 + 1, ny - i2 + 1) = 0
          ikm(nx - i1 + 1, ny + i2 - 1) = 0
        END IF
        IF (sqak2 <= akmax) THEN
          ikm(nx + i1 - 1, ny + i2 - 1) = INT(SQRT(akmax**2 - ak2)/dkz) + 1
          IF (ikm(nx + i1 - 1, ny + i2 - 1) > nz1) ikm(nx + i1 - 1, ny + i2 - 1) = nz1
          IF ((i1 < nx1) .AND. (i2 < ny1)) THEN
            ikm(nx + i1 - 1, ny - i2 + 1) = ikm(nx + i1 - 1, ny + i2 - 1)
            ikm(nx - i1 + 1, ny - i2 + 1) = ikm(nx + i1 - 1, ny + i2 - 1)
            ikm(nx - i1 + 1, ny + i2 - 1) = ikm(nx + i1 - 1, ny + i2 - 1)
          END IF
        END IF
      END DO
    END DO

    ii = 0
    i0 = 0
    xz1 = -nz*dkz
    DO i3 = 1, nzi
      xz1 = xz1 + dkz
      xz2 = xz1*xz1
      xy1 = -ny*dky
      DO i2 = 1, nyi
        xy1 = xy1 + dky
        xy2 = xy1*xy1
        xx1 = -nx*dkx
        DO i1 = 1, nxi
          xx1 = xx1 + dkx
          xx2 = xx1*xx1
          i0 = i0 + 1
          ak2 = xx2 + xy2 + xz2
          sqak2 = SQRT(ak2)
          indfc(i0) = 0
          IF (sqak2 <= akmax) THEN
            ii = ii + 1
            IF (ii > kdred) THEN
              WRITE (6, '(a,a)') ' e r r o r IN i n c l u d e ', &
                ' parameter kdred too small'
              STOP 'kdred to  small, reduce ecut'
            END IF
          END IF
          indfc(i0) = ii
          akv2(ii) = ak2
        END DO
      END DO
    END DO
    nksp = ii

    RETURN
  END SUBROUTINE fftinp

!-----fourf-------------------------------------------------------fourf

  SUBROUTINE fourf(psx, pskr, pski, ipar)
    IMPLICIT NONE

! Fourier forward transformation
!
! Input:
! psx input wave-function
! ipar parity in x- and y-direction
! Output:
! pskr real part of the wave-function
! pski imaginary part of the wave-function

    REAL(DP), INTENT(IN OUT) :: psx(kdfull)
    REAL(DP), INTENT(OUT) :: pskr(kdred)
    REAL(DP), INTENT(OUT) :: pski(kdred)
    INTEGER, INTENT(IN) :: ipar(3)

    INTEGER :: i1, i
    REAL(DP) :: a(kdfull)
    REAL(DP) :: mxini, myini, mzini
    REAL(DP) :: nzzh, nyyh, sqh, tnorm, xp, yp, zp

    DATA mxini, myini, mzini/0, 0, 0/ ! flag for initialization

!----------------------------------------------------------------------

! initialization NETLIB
! END initialization NETLIB

! initialization FFTW3
    IF (mxini == 0) THEN
#if(fftwnomkl)
      wisdomtest = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'wisdom_fftw.dat not found, creating it'
        WRITE (7, *) 'wisdom_fftw.dat not found, creating it'
      END IF
      CALL dfftw_plan_dft_1d(pforwx, kfftx, fftax, fftax, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      CALL dfftw_plan_dft_1d(pbackx, kfftx, fftax, fftax, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      mxini = kfftx
      WRITE (7, '(a)') ' x-fft initialized '
    ELSE IF (mxini /= kfftx) THEN
      STOP ' nx2 in four3d not as initialized!'
    END IF
    IF (myini == 0) THEN
      CALL dfftw_plan_dft_1d(pforwy, kfftx, fftax, fftax, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      CALL dfftw_plan_dft_1d(pbacky, kfftx, fftax, fftax, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      myini = kffty
      WRITE (7, '(a)') ' y-fft initialized '
    ELSE IF (myini /= kffty) THEN
      STOP ' ny2 in four3d not as initialized!'
    END IF
    IF (mzini == 0) THEN
      CALL dfftw_plan_dft_1d(pforwz, kfftx, fftax, fftax, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      CALL dfftw_plan_dft_1d(pbackz, kfftx, fftax, fftax, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      mzini = kfftz
#if(fftwnomkl)
      wisdomtest = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'Error exporting wisdom to  FILE wisdom_fftw.dat'
        WRITE (7, *) 'Error exporting wisdom to  FILE wisdom_fftw.dat'
      END IF
      WRITE (7, '(a)') ' z-fft initialized '
    ELSE IF (mzini /= kfftz) THEN
      STOP ' nz2 in four3d not as initialized!'
    END IF
! END initialization FFTW3

    nzzh = (nzi - 1)*nxy1
    nyyh = (nyi - 1)*nxi
    sqh = SQRT(0.5D0)
    tnorm = grnorm*fnorm

    xp = ipar(1)
    yp = ipar(2)
    zp = ipar(3)

    a = 0D0

    CALL fftx(psx, a)
    CALL ffty(psx, a)
    CALL fftz(psx, a)

    DO i1 = 1, nkxyz
      i = indfc(i1)
      IF (i > 0) THEN
        pskr(i) = tnorm*psx(i1)
        pski(i) = tnorm*a(i1)
      END IF
    END DO

    RETURN
  END SUBROUTINE fourf

!-----fourb------------------------------------------------------------

  SUBROUTINE fourb(psx, pskr, pski, ipar)

! Fourier backward transformation
!
! Input:
! pskr real part of the wave-function
! pski imaginary part of the wave-function
! ipar parity in x- and y-direction
! Output:
! psx input wave-function

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: psx(kdfull)
    REAL(DP), INTENT(IN) :: pskr(kdred)
    REAL(DP), INTENT(IN) :: pski(kdred)
    INTEGER, INTENT(IN) :: ipar(3)

    INTEGER :: i, i1
    REAL(DP) :: nzzh, nyyh, sq2, tnorm, xp, yp, zp
    REAL(DP) :: a(kdfull)

!----------------------------------------------------------------------

    nzzh = (nzi - 1)*nxy1
    nyyh = (nyi - 1)*nx1
    sq2 = SQRT(2D0)
    tnorm = fnorm/(8D0*grnorm)

    xp = ipar(1)
    yp = ipar(2)
    zp = ipar(3)

    DO i1 = 1, nkxyz
      i = indfc(i1)
      psx(i1) = 0D0
      a(i1) = 0D0
      IF (i > 0) THEN
        psx(i1) = tnorm*pskr(i)
        a(i1) = tnorm*pski(i)
      END IF
    END DO

    CALL ffbz(psx, a)
    CALL ffby(psx, a)
    CALL ffbx(psx, a)

    RETURN
  END SUBROUTINE fourb

!-----fftx-------------------------------------------------------------

  SUBROUTINE fftx(psxr, psxi)

! Performs the Fourier-transformation in x-direction.
! The input-wave-function (psxr,psxi) (i.e. real and imaginary part)
! is overwritten by the Fourier-transformed wave-function.

    IMPLICIT NONE

    REAL(DP), INTENT(IN OUT) :: psxr(kdfull)
    REAL(DP), INTENT(OUT) :: psxi(kdfull)
    INTEGER :: i0, i1, i2, i3, i30, ii, nx11

!----------------------------------------------------------------------


    nx11 = nx1 + 1
    i30 = -nxy1
    DO i3 = 1, nzi
      i30 = i30 + nxy1
      i0 = i30 - nxi
      DO i2 = 1, nyi
        i0 = i0 + nxi
! composition of the wave-function
! positive space
        ii = i0 + nx - 1
        DO i1 = 1, nx1
          ii = ii + 1
          fftax(i1) = CMPLX(psxr(ii), 0D0, DP)
        END DO
! negative space
        ii = i0
        DO i1 = nx11, nxi
          ii = ii + 1
          fftax(i1) = CMPLX(psxr(ii), 0D0, DP)
        END DO
! execution of the Fourier-transformation
        CALL fftw_execute_dft(pforwx, fftax, fftax)
! decomposition of the wave-function
! positive space
        ii = i0 + nx - 1
        DO i1 = 1, nx1
          ii = ii + 1
          psxr(ii) = REAL(fftax(i1), DP)
          psxi(ii) = AIMAG(fftax(i1))
        END DO
! negative space
        ii = i0
        DO i1 = nx11, nxi
          ii = ii + 1
          psxr(ii) = REAL(fftax(i1), DP)
          psxi(ii) = AIMAG(fftax(i1))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE fftx

!-----ffty--------------------------------------------------------------

  SUBROUTINE ffty(psxr, psxi)

! Performs the Fourier-transformation in y-direction.
! The input-wave-function (psxr,psxi) (i.e. real and imaginary part)
! is overwritten by the Fourier-transformed wave-function

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: psxr(kdfull)
    REAL(DP), INTENT(IN OUT) :: psxi(kdfull)
    INTEGER :: i0, i1, i2, i3, i30, ii, ny11

!----------------------------------------------------------------------


    ny11 = ny1 + 1
    IF (nxk < nx1) THEN
      nxklo = nx - nxk + 1
    ELSE
      nxklo = 1
    END IF
    nxkhi = nx + nxk - 1

    i30 = -nxy1
    DO i3 = 1, nzi
      i30 = i30 + nxy1
      DO i1 = nxklo, nxkhi
        i0 = i1 + i30
! composition of the wave-function
! positive space
        ii = i0 + nxi*(ny - 2)
        DO i2 = 1, ny1
          ii = ii + nxi
          fftay(i2) = CMPLX(psxr(ii), psxi(ii), DP)
        END DO
! negative space
        ii = i0 - nxi
        DO i2 = ny11, nyi
          ii = ii + nxi
          fftay(i2) = CMPLX(psxr(ii), psxi(ii), DP)
        END DO
! execution of the Fourier-transformation
        CALL fftw_execute_dft(pforwy, fftay, fftay)

! decomposition of the wave-function
! positive space
        ii = i0 + nxi*(ny - 2)
        DO i2 = 1, ny1
          ii = ii + nxi
          psxr(ii) = REAL(fftay(i2), DP)
          psxi(ii) = AIMAG(fftay(i2))
        END DO
! negative space
        ii = i0 - nxi
        DO i2 = ny11, nyi
          ii = ii + nxi
          psxr(ii) = REAL(fftay(i2), DP)
          psxi(ii) = AIMAG(fftay(i2))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE ffty

!-----fftz-------------------------------------------------------------
  SUBROUTINE fftz(psxr, psxi)

! Performs the Fourier-transformation in z-direction.
! The input-wave-function (psxr,psxi) (i.e. real and imaginary part)
! is overwritten by the Fourier-transformed wave-function

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: psxr(kdfull)
    REAL(DP), INTENT(IN OUT) :: psxi(kdfull)
    INTEGER :: nxyf
    INTEGER :: nyf
    INTEGER :: nzh
    INTEGER :: i1, i2, i3, i3m, ind

!----------------------------------------------------------------------
    nxyf = kfftx*kffty
    nyf = kfftx
    nzh = kfftz/2


    DO i2 = 1, kffty
      DO i3 = 1, kfftz
        i3m = MOD(i3 + nzh, kfftz) + 1
        DO i1 = 1, kfftx
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftb(i3m, i1) = CMPLX(psxr(ind), psxi(ind), DP)
        END DO
      END DO
      DO i1 = 1, kfftx
        CALL fftw_execute_dft(pforwz, fftb(1, i1), fftb(1, i1))
      END DO
      DO i3 = 1, kfftz
        i3m = MOD(i3 + nzh, kfftz) + 1
        DO i1 = 1, kfftx
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          psxr(ind) = REAL(fftb(i3m, i1), DP)
          psxi(ind) = AIMAG(fftb(i3m, i1))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE fftz

!-----ffbz-------------------------------------------------------------
  SUBROUTINE ffbz(psxr, psxi)

! Performs the Fourier-backtransformation in z-direction.
! The input/output-wave-function (psxr,psxi) (i.e. real and imaginary part)
! is overwritten by the Fourier-transformed wave-function

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: psxr(kdfull)
    REAL(DP), INTENT(IN OUT) :: psxi(kdfull)
    INTEGER :: nxyf
    INTEGER :: nyf
    INTEGER :: nzh
    INTEGER :: i1, i2, i3, i3m, ind

!----------------------------------------------------------------------

    nxyf = kfftx*kffty
    nyf = kfftx
    nzh = kfftz/2


    DO i2 = 1, kffty
      DO i3 = 1, kfftz
        i3m = MOD(i3 + nzh, kfftz) + 1
        DO i1 = 1, kfftx
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftb(i3m, i1) = CMPLX(psxr(ind), psxi(ind), DP)
        END DO
      END DO
      DO i1 = 1, kfftx
        CALL fftw_execute_dft(pbackz, fftb(1, i1), fftb(1, i1))
      END DO
      DO i3 = 1, kfftz ! copy back
        i3m = MOD(i3 + nzh, kfftz) + 1
        DO i1 = 1, kfftx
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          psxr(ind) = REAL(fftb(i3m, i1), DP)
          psxi(ind) = AIMAG(fftb(i3m, i1))
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE ffbz

!-----ffby-------------------------------------------------------------
  SUBROUTINE ffby(psxr, psxi)

! Performs the Fourier-backtransformation in y-direction.
! The input/output-wave-function (psxr,psxi) (i.e. real and imaginary part)
! is overwritten by the Fourier-transformed wave-function

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: psxr(kdfull)
    REAL(DP), INTENT(IN OUT) :: psxi(kdfull)
    INTEGER :: i0, i1, i2, i3, i30, ii, ny11

!----------------------------------------------------------------------

    ny11 = ny1 + 1

    IF (nxk < nx1) THEN
      nxklo = nx - nxk + 1
    ELSE
      nxklo = 1
    END IF
    nxkhi = nx + nxk - 1


    i30 = -nxy1
    DO i3 = 1, nzi
      i30 = i30 + nxy1
      DO i1 = nxklo, nxkhi
        i0 = i1 + i30
! composition of the complex wave-function
! positive space
        ii = i0 + nxi*(ny - 2)
        DO i2 = 1, ny1
          ii = ii + nxi
          fftay(i2) = CMPLX(psxr(ii), psxi(ii), DP)
        END DO
! negative space
        ii = i0 - nxi
        DO i2 = ny11, nyi
          ii = ii + nxi
          fftay(i2) = CMPLX(psxr(ii), psxi(ii), DP)
        END DO
! execution
        CALL fftw_execute_dft(pbacky, fftay, fftay)
! decomposition of the inverse transformed wave-function
! positive space
        ii = i0 + nxi*(ny - 2)
        DO i2 = 1, ny1
          ii = ii + nxi
          psxr(ii) = REAL(fftay(i2), DP)
          psxi(ii) = AIMAG(fftay(i2))
        END DO
! negative space
        ii = i0 - nxi
        DO i2 = ny11, nyi
          ii = ii + nxi
          psxr(ii) = REAL(fftay(i2), DP)
          psxi(ii) = AIMAG(fftay(i2))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE ffby

!-----ffbx-------------------------------------------------------------
  SUBROUTINE ffbx(psxr, psxi)

! Performs the Fourier-backtransformation in x-direction.
! The input/output-wave-function (psxr,psxi) (i.e. real and imaginary part)
! is overwritten by the Fourier-transformed wave-function

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: psxr(kdfull)
    REAL(DP), INTENT(IN OUT) :: psxi(kdfull)
    INTEGER :: i0, i1, i2, i3, i30, ii, nx11

!----------------------------------------------------------------------

    nx11 = nx1 + 1

    i30 = -nxy1


    DO i3 = 1, nzi
      i30 = i30 + nxy1
      i0 = i30 - nxi
      DO i2 = 1, nyi
        i0 = i0 + nxi
! composition of the complex wave-function
! positive space
        ii = i0 + nx - 1
        DO i1 = 1, nx1
          ii = ii + 1
          fftax(i1) = CMPLX(psxr(ii), psxi(ii), DP)
        END DO
! negative space
        ii = i0
        DO i1 = nx11, nxi
          ii = ii + 1
          fftax(i1) = CMPLX(psxr(ii), psxi(ii), DP)
        END DO
! execution
        CALL fftw_execute_dft(pbackx, fftax, fftax)
! decomposition of the inverse transformed wave-function
! positive space
        ii = i0 + nx - 1
        DO i1 = 1, nx1
          ii = ii + 1
          psxr(ii) = REAL(fftax(i1), DP)
          psxi(ii) = AIMAG(fftax(i1))
        END DO
! negative space
        ii = i0
        DO i1 = nx11, nxi
          ii = ii + 1
          psxr(ii) = REAL(fftax(i1), DP)
          psxi(ii) = AIMAG(fftax(i1))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE ffbx

  SUBROUTINE coulsolv_end()

! Epilogie for FFTW3

    CALL fftw_destroy_plan(pforwx)
    CALL fftw_destroy_plan(pforwy)
    CALL fftw_destroy_plan(pforwz)
    CALL fftw_destroy_plan(pbackx)
    CALL fftw_destroy_plan(pbacky)
    CALL fftw_destroy_plan(pbackz)

    DEALLOCATE (xval, yval, zval)
    DEALLOCATE (xt2, yt2, zt2)
    DEALLOCATE (fftax, fftay, fftb)
    DEALLOCATE (ikm, indfc)
    DEALLOCATE (akv2)
    DEALLOCATE (potc, rho, rhozw)
    DEALLOCATE (qlcows, rq)
    DEALLOCATE (rokall, pcoall)

    RETURN
  END SUBROUTINE coulsolv_end

END MODULE coulsolv_f
