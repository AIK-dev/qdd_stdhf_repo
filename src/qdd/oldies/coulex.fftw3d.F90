
MODULE coulsolv
  USE, INTRINSIC :: iso_c_binding
  USE params, ONLY: DP
  USE FFTW
  IMPLICIT REAL(DP) (A - H, O - Z)

  SAVE
  INTEGER, PRIVATE :: kxmax, kymax, kzmax, ksmax
! kxamx must be the largest
  INTEGER, PRIVATE :: kdfull
  INTEGER, PRIVATE :: kdred
  INTEGER, PRIVATE :: kfft2
!INTEGER,PARAMETER,PRIVATE :: kddoub=kdfull
  INTEGER, PRIVATE :: kfft, kfftx, kffty, kfftz
!INTEGER,PARAMETER,PRIVATE :: kdcorf=(kxmax/2+1)*(kymax/2+1)*(kzmax/2+1)
! include BLOCK: xkgrid
  REAL(DP), ALLOCATABLE, PRIVATE :: xval(:), yval(:), zval(:)
  REAL(DP), ALLOCATABLE, PRIVATE :: xt2(:), yt2(:), zt2(:)
  REAL(DP), PRIVATE :: dx, dy, dz, dxsp, grnorm, fnorm
  INTEGER, PRIVATE :: nx, ny, nz, nx1, ny1, nz1, nxr, nxi, nyr, nyi, nzr, nzi, nxy1, nxyz
  INTEGER, PRIVATE :: nxhigh, nxlow, nyhigh, nylow, nzhigh, nzlow
!COMMON /xgrid/ xval,yval,zval,xt2,yt2,zt2,dx,dy,dz,dxsp,grnorm,fnorm, &
! nx,ny,nz,nx1,ny1,nz1, nxr,nxi,nyr,nyi,nzr,nzi,nxy1,nxyz, &
! nxhigh,nxlow,nyhigh,nylow,nzhigh,nzlow

  REAL(DP), ALLOCATABLE, PRIVATE :: akv2r(:), akv2i(:)
  INTEGER, ALLOCATABLE, PRIVATE :: ikm(:, :)
  REAL(DP), PRIVATE :: dkx, dky, dkz, akmax, dksp, ecut
  INTEGER, PRIVATE :: nxk, nxklo, nxkhi, nksp, nkxyz
!COMMON /kgrid/ akv2r,akv2i, dkx,dky,dkz,akmax,dksp, &
! ikm,nxk,nxklo,nxkhi,nksp,nkxyz,ecut

  TYPE(C_PTR), PRIVATE :: pforw, pback
  INTEGER(C_INT), PRIVATE :: wisdomtest
!COMMON /fftini/wrkx,wrky,wrkz,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz

! include BLOCK: option
  REAL(DP), PARAMETER, PRIVATE :: zero = 0D0
  REAL(DP), PARAMETER, PRIVATE :: pi = 3.141592653589793D0

!COMPLEX(C_DOUBLE_COMPLEX),ALLOCATABLE,PRIVATE :: fftax(:),fftay(:),fftb(:,:),ffta(:,:,:)
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, PRIVATE :: ffta(:, :, :)
!COMMON /fftcom/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax)

CONTAINS

  SUBROUTINE init_coul(dx0, dy0, dz0, nx0, ny0, nz0)

!-----------------------------------------------------------------------

! READ grid parameters from FILE or simply initialize them
! note that the Coulomb solver doubles the grid internally
    kxmax = 2*nx0; kymax = 2*ny0; kzmax = 2*nz0; ksmax = kxmax
    kdfull = nx0*ny0*nz0
    kdred = kxmax*kymax*kzmax
    kfft = ksmax; kfftx = kxmax; kffty = kymax; kfftz = kzmax
    kfft2 = kfft*2

    nx = nx0 !/2
    ny = ny0 !/2
    nz = nz0 !/2
    dx = dx0
    dy = dy0
    dz = dz0

    ALLOCATE (xval(kxmax), yval(kymax), zval(kzmax))
    ALLOCATE (xt2(kxmax), yt2(kymax), zt2(kzmax))
!ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
    ALLOCATE (ffta(kxmax, kymax, kzmax))
    ALLOCATE (akv2r(kdred), akv2i(kdred))
    ALLOCATE (ikm(kxmax, kymax))

! CALL input routine fftinp, which initializes the grid and fft tabl

    CALL fftinp

    RETURN
  END SUBROUTINE init_coul

!-----fftinp------------------------------------------------------------

  SUBROUTINE fftinp
    IMPLICIT REAL(DP) (A - H, O - Z)

! initializes work tables for FFT

! grid parameters nx,ny,nz,dx,dy,dz,ecut must have been READ or
! initialized before !

!-----------------------------------------------------------------------
    INTEGER, SAVE :: fini = 0
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
    nkxyz = nxi*nyi*nzi

! grid lengths must match with parameters IN incs

    IF (kxmax < nxi) THEN
      WRITE (6, '(a)') ' ERROR: PARAMETER kxmax too small'
      STOP ' error IN PARAMETER: KXMAX IN COULEX too small'
    ELSE IF (kymax < nyi) THEN
      WRITE (6, '(a)') ' ERROR: PARAMETER kymax too small'
      STOP ' error IN PARAMETER: KYMAX IN COULEX too small'
    ELSE IF (kzmax < nzi) THEN
      WRITE (6, '(a)') ' ERROR: PARAMETER kzmax too small'
      STOP ' error IN PARAMETER: KZMAX IN COULEX too small'
    END IF

! initialize grid IN fourier space

    dkx = pi/(dx*REAL(nx))
    dky = pi/(dy*REAL(ny))
    dkz = pi/(dz*REAL(nz))

    dxsp = dx*dy*dz
    dksp = dkx*dky*dkz
    WRITE (*, *) ' dkx,dky,dkz,dksp=', dkx, dky, dkz, dksp

    grnorm = SQRT(dxsp/dksp)
    fnorm = 1.0/SQRT(REAL(nx*ny*nz))
!test akmax=sqrt(3*(nx*nx)*dx*dx)+2.0
!test nxk=int(akmax/dkx)+1
!test IF(nxk.gt.nx1) nxk=nx1
    nxk = nx1

! built Greens FUNCTION IN Fourier space
! by Fourier transformation from REAL space

    ikzero = nxy1*(nz - 1) + nxi*(ny - 1) + nx
    WRITE (*, *) ' nzi,nyi,nxi,nx,ny,nz,ikzero=', nzi, nyi, nxi, nx, ny, nz, ikzero
    ii = 0
    xz1 = -nz*dz
    DO i3 = 1, nzi
      xz1 = xz1 + dz
      xz2 = xz1*xz1
      xy1 = -ny*dy
      DO i2 = 1, nyi
        xy1 = xy1 + dy
        xy2 = xy1*xy1
        xx1 = -nx*dx
        DO i1 = 1, nxi
          xx1 = xx1 + dx
          xx2 = xx1*xx1
          ak2 = xx2 + xy2 + xz2
          ii = ii + 1
! WRITE(*,*) ' i1,i2,i3,ii=',i1,i2,i3,ii
          IF (ii /= ikzero) THEN
            akv2r(ii) = 1D0/SQRT(ak2)
          ELSE
! akv2r(ii) = (6D0*pi/(dx*dy*dz))**(1D0/3D0) ! spherical approx
! akv2r(ii) = 1.19003868*(dx*dy*dz)**(-1D0/3D0)
            akv2r(ii) = 2.34D0*1.19003868D0*(dx*dy*dz)**(-1D0/3D0) ! empirical
          END IF
          akv2i(ii) = 0D0
        END DO
      END DO
    END DO
    nksp = ii

    CALL fourf(akv2r(1), akv2i(1))

!STOP

    RETURN
  END SUBROUTINE fftinp

!-------------------------------------------------------------------

  SUBROUTINE falr(rhoinp, chpfalr, nxdum, nydum, nzdum, kdum)

    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN) :: rhoinp(kdfull)
    REAL(DP), INTENT(OUT) :: chpfalr(kdfull)
    INTEGER, INTENT(IN) :: nxdum
    INTEGER, INTENT(IN) :: nydum
    INTEGER, INTENT(IN) :: nzdum
    INTEGER, INTENT(IN) :: kdum

    REAL(DP), ALLOCATABLE :: rhokr(:), rhoki(:)

    ALLOCATE (rhokr(kdred), rhoki(kdred))

! CALL a routine written by you which writes your density field
! on the array rho.
! remember not TO send your original density array TO the fcs.
! IN this CASE we have a homogeneously charged sphere .

    CALL rhofld(rhoinp, rhokr, rhoki)

! CALL coufou, which CONTAINS the fcs PROCEDURE.

    CALL coufou2(rhokr, rhoki)

! CALL a routine written by you which outputs the results of the fcs
! and maybe some other things TO an output FILE or the screen.

    CALL RESULT(chpfalr, rhokr, rhoki)

    DEALLOCATE (rhokr, rhoki)

  END SUBROUTINE falr

!-----rhofld------------------------------------------------------------

  SUBROUTINE rhofld(rhoinp, rhokr, rhoki)
    IMPLICIT REAL(DP) (A - H, O - Z)

! copy density on COMPLEX array of DOUBLE extnesion IN x,y,z

    REAL(DP), INTENT(IN) :: rhoinp(kdfull)
    REAL(DP), INTENT(OUT) :: rhokr(kdred)
    REAL(DP), INTENT(OUT) :: rhoki(kdred)

    ii = 0
    i0 = 0
    DO i3 = 1, nzi
      DO i2 = 1, nyi
        DO i1 = 1, nxi
          ii = ii + 1
          IF (i3 <= nz .AND. i2 <= ny .AND. i1 <= nx) THEN
            i0 = i0 + 1
            rhokr(ii) = rhoinp(i0)
          ELSE
            rhokr(ii) = 0D0
          END IF
          rhoki(ii) = 0D0
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE rhofld

!-----RESULT------------------------------------------------------------

  SUBROUTINE RESULT(chpfalr, rhokr, rhoki)

    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(OUT) :: chpfalr(kdfull)
    REAL(DP), INTENT(IN) :: rhokr(kdred)
    REAL(DP), INTENT(IN OUT) :: rhoki(kdred)

! copy Coulomb field back TO standard grid

    ii = 0
    i0 = 0
    DO i3 = 1, nzi
      DO i2 = 1, nyi
        DO i1 = 1, nxi
          ii = ii + 1
          IF (i3 <= nz .AND. i2 <= ny .AND. i1 <= nx) THEN
            i0 = i0 + 1
            chpfalr(i0) = 2D0*rhokr(ii)
          END IF
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE RESULT

!-----cofows------------------------------------------------------------

  SUBROUTINE coufou2(rhokr, rhoki)

    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN OUT) :: rhokr(kdred)
    REAL(DP), INTENT(IN OUT) :: rhoki(kdred)

!INTEGER,SAVE :: fini=0
    LOGICAL, PARAMETER :: tprint = .false.
    LOGICAL, PARAMETER :: rqplot = .false.

!------------------------------------------------------------------------------

! fourier transformation of the density

    CALL fourf(rhokr, rhoki)

! calculation of the coulomb field (writing on the density field)

    DO ik = 1, kdred
      SAVE2 = akv2r(ik)*rhokr(ik) + akv2i(ik)*rhoki(ik)
      rhoki(ik) = akv2r(ik)*rhoki(ik) + akv2i(ik)*rhokr(ik)
      rhokr(ik) = SAVE2
    END DO

! fourier back transformation

    CALL fourb(rhokr, rhoki)

    RETURN

  END SUBROUTINE coufou2

!-----fourf-------------------------------------------------------fourf

  SUBROUTINE fourf(pskr, pski)
    USE, INTRINSIC :: iso_c_binding
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(OUT) :: pskr(kdred)
    REAL(DP), INTENT(OUT) :: pski(kdred)

! fourier forward transformation
! I/O: pskr REAL part of the wave-FUNCTION
! pski imaginary part of the wave-FUNCTION

    DATA mini/0/ ! flag for initialization

!----------------------------------------------------------------------

! check initialization

    IF (mini == 0) THEN
      wisdomtest = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'wisdom_fftw.dat not found, creating it'
        WRITE (7, *) 'wisdom_fftw.dat not found, creating it'
      END IF
      pforw = fftw_plan_dft_3d(kzmax, kymax, kxmax, ffta, ffta, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      pback = fftw_plan_dft_3d(kzmax, kymax, kxmax, ffta, ffta, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      mini = kxmax*kymax*kzmax
      wisdomtest = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'Error exporting wisdom TO FILE wisdom_fftw.dat'
        WRITE (7, *) 'Error exporting wisdom TO FILE wisdom_fftw.dat'
      END IF
      WRITE (7, '(a)') ' fft initialized '
    ELSE IF (mini /= kxmax*kymax*kzmax) THEN
      STOP ' n2 IN four3d not as initialized!'
    END IF

    nzzh = (nzi - 1)*nxy1
    nyyh = (nyi - 1)*nxi
!test sqh=sqrt(0.5)
    tnorm = grnorm*fnorm

    CALL fft(pskr, pski)

    DO i1 = 1, nkxyz
      pskr(i1) = tnorm*pskr(i1)
      pski(i1) = tnorm*pski(i1)
    END DO

    RETURN
  END SUBROUTINE fourf
!-----fourb------------------------------------------------------------

  SUBROUTINE fourb(pskr, pski)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(OUT) :: pskr(kdred)
    REAL(DP), INTENT(OUT) :: pski(kdred)

! fourier backward transformation
! I/O: pskr REAL part of the wave-FUNCTION
! pski imaginary part of the wave-FUNCTION
!----------------------------------------------------------------------

    nzzh = (nzi - 1)*nxy1
    nyyh = (nyi - 1)*nx1
! sq2=sqrt(2.0)
    tnorm = fnorm/(8D0*grnorm)*pi**1.5D0
! tnorm=fnorm/(16D0*grnorm)

    DO i = 1, nkxyz
      pskr(i) = tnorm*pskr(i)
      pski(i) = tnorm*pski(i)
    END DO

    CALL ffb(pskr, pski)

    RETURN
  END SUBROUTINE fourb

!-----fft--------------------------------------------------------------

  SUBROUTINE fft(psxr, psxi)

    USE, INTRINSIC :: iso_c_binding
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN OUT) :: psxr(kdred)
    REAL(DP), INTENT(IN OUT) :: psxi(kdred)

    CALL copyr1dto3d(psxr, ffta, kfftx, kffty, kfftz)

    CALL fftw_execute_dft(pforw, ffta, ffta)

    CALL copy3dto1d(ffta, psxr, psxi, kfftx, kffty, kfftz)

    RETURN
  END SUBROUTINE fft

!-----ffb--------------------------------------------------------------

  SUBROUTINE ffb(psxr, psxi)

    USE, INTRINSIC :: iso_c_binding
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN OUT) :: psxr(kdred)
    REAL(DP), INTENT(IN OUT) :: psxi(kdred)

!----------------------------------------------------------------------

    CALL secopy1dto3d(psxr, psxi, ffta, kfftx, kffty, kfftz)

    CALL fftw_execute_dft(pback, ffta, ffta)

    CALL copyr3dto1d(ffta, psxr, psxi, kfftx, kffty, kfftz)

    RETURN
  END SUBROUTINE ffb

! ******************************

  SUBROUTINE copyr1dto3d(psxr, ffta, nbx2, nby2, nbz2)

! ******************************

    USE params
    USE, INTRINSIC :: iso_c_binding

    REAL(DP), INTENT(IN) :: psxr(kdred)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ffta(nbx2, nby2, nbz2)
    INTEGER :: nxyfn, nyfn, nnx, nny, nnz

    nxyfn = nbx2*nbz2
    nyfn = kfftx
    nnx2 = nx + nx
    nny2 = ny + ny
    nnz2 = nz + nz

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyfn + (i2 - 1)*nyfn + i1
          ffta(MOD(i1 + nnx2, nbx2) + 1, MOD(i2 + nny2, nby2) + 1, MOD(i3 + nnz2, nbz2) + 1) = CMPLX(psxr(ind), 0D0, DP)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copyr1dto3d

! ******************************

  SUBROUTINE copy3dto1d(ffta, psxr, psxi, nbx2, nby2, nbz2)

! ******************************

    USE params
    USE, INTRINSIC :: iso_c_binding

    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: ffta(nbx2, nby2, nbz2)
    REAL(DP), INTENT(OUT) :: psxr(kdred), psxi(kdred)
    INTEGER :: nxyfn, nyfn, nnx, nny, nnz

    nxyfn = nbx2*nbz2
    nyfn = kfftx
    nnx2 = nx + nx
    nny2 = ny + ny
    nnz2 = nz + nz

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyfn + (i2 - 1)*nyfn + i1
          psxr(ind) = REAL(ffta(MOD(i1 + nnx2, nbx2) + 1, MOD(i2 + nny2, nby2) + 1, MOD(i3 + nnz2, nbz2) + 1))
          psxi(ind) = AIMAG(ffta(MOD(i1 + nnx2, nbx2) + 1, MOD(i2 + nny2, nby2) + 1, MOD(i3 + nnz2, nbz2) + 1))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copy3dto1d

! ******************************

  SUBROUTINE secopy1dto3d(psxr, psxi, ffta, nbx2, nby2, nbz2)

! ******************************

    USE params
    USE, INTRINSIC :: iso_c_binding

    REAL(DP), INTENT(IN) :: psxr(kdred), psxi(kdred)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ffta(nbx2, nby2, nbz2)
    INTEGER :: nxyfn, nyfn, nnx, nny, nnz

    nxyfn = nbx2*nbz2
    nyfn = kfftx
    nnx2 = nx + nx
    nny2 = ny + ny
    nnz2 = nz + nz

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyfn + (i2 - 1)*nyfn + i1
          ffta(MOD(i1 + nnx2, nbx2) + 1, MOD(i2 + nny2, nby2) + 1, MOD(i3 + nnz2, nbz2) + 1) = CMPLX(psxr(ind), psxi(ind), DP)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE secopy1dto3d

! ******************************

  SUBROUTINE copyr3dto1d(ffta, psxr, psxi, nbx2, nby2, nbz2)

! ******************************

    USE params
    USE, INTRINSIC :: iso_c_binding

    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: ffta(nbx2, nby2, nbz2)
    REAL(DP), INTENT(OUT) :: psxr(kdred), psxi(kdred)
    INTEGER :: nxyfn, nyfn, nnx, nny, nnz

    nxyfn = nbx2*nbz2
    nyfn = kfftx
    nnx2 = nx + nx
    nny2 = ny + ny
    nnz2 = nz + nz

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyfn + (i2 - 1)*nyfn + i1
          psxr(ind) = REAL(ffta(MOD(i1 + nnx2, nbx2) + 1, MOD(i2 + nny2, nby2) + 1, MOD(i3 + nnz2, nbz2) + 1))
          psxi(ind) = AIMAG(ffta(MOD(i1 + nnx2, nbx2) + 1, MOD(i2 + nyy2, nby2) + 1, MOD(i3 + nnz2, nbz2) + 1))
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copyr3dto1d

! ******************************

  SUBROUTINE coulex_end()

    CALL fftw_destroy_plan(pforw)
    CALL fftw_destroy_plan(pback)

    DEALLOCATE (xval, yval, zval)
    DEALLOCATE (xt2, yt2, zt2)
    DEALLOCATE (ffta)
    DEALLOCATE (akv2r, akv2i)
    DEALLOCATE (ikm)

  END SUBROUTINE coulex_end

! ******************************

END MODULE coulsolv
