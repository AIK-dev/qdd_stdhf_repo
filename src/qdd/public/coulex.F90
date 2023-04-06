
MODULE coulsolv_e

! Exact Coulomb solver using doubled grid to treat long-range
! part exactly.

  USE, INTRINSIC :: iso_c_binding
  USE FFTW
  USE kinetic, ONLY: FFTW_planflag
  USE params, ONLY: DP, PI, numthr, e2, zero, nx2, ny2, nz2, systime_factor
  IMPLICIT NONE

  SAVE
  INTEGER, PRIVATE :: kxmax, kymax, kzmax, ksmax
! kxmax must be the largest
  INTEGER, PRIVATE :: kdfull
  INTEGER, PRIVATE :: kdred
  INTEGER, PRIVATE :: kfft2
!INTEGER,PARAMETER,PRIVATE :: kddoub=kdfull
  INTEGER, PRIVATE :: kfft, kfftx, kffty, kfftz
!INTEGER,PARAMETER,PRIVATE :: kdcorf=(kxmax/2+1)*(kymax/2+1)*(kzmax/2+1)
! include BLOCK: xkgrid
  REAL(DP), ALLOCATABLE, PRIVATE :: xt2(:), yt2(:), zt2(:)
  REAL(DP), PRIVATE :: dx, dy, dz, dxsp, grnorm, fnorm
  INTEGER, PRIVATE :: nx, ny, nz, nx1, ny1, nz1, nxi, nyi, nzi, nxy1, nxyz
  INTEGER, PRIVATE :: nxhigh, nxlow, nyhigh, nylow, nzhigh, nzlow
  INTEGER, ALLOCATABLE, PRIVATE :: ikm(:, :)
  REAL(DP), PRIVATE :: dkx, dky, dkz, akmax, dksp, ecut
  INTEGER, PRIVATE :: nxk, nxklo, nxkhi, nksp, nkxyz, iret


  COMPLEX(C_DOUBLE_COMPLEX), POINTER, PRIVATE :: akv2(:, :, :)
  INTEGER, PRIVATE :: nini = 0
  INTEGER(C_INT), PRIVATE :: wisdomtest
  TYPE(C_PTR), PRIVATE :: pforwx, pforwy, pforwz, pbackx, pbacky, pbackz
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, PRIVATE :: fftax(:), fftay(:), fftb(:, :)
  TYPE(C_PTR), PRIVATE :: pforw, pback, pforwc
  COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, POINTER :: ffta(:, :, :)
  REAL(C_DOUBLE), PRIVATE, POINTER :: rffta(:, :, :)
  REAL(C_DOUBLE), PRIVATE, POINTER :: rfftb(:, :, :)
  TYPE(C_PTR) :: pakv2, pffta, prffta, prfftb

CONTAINS

  SUBROUTINE init_coul_e(dx0, dy0, dz0, nx0, ny0, nz0)

! Initialize grid parameters, basic arrays, and FFTW3 plans.
!
! Input:
! dx0,dy0,dz0 = grid spacings
! nx0,ny0,nz0 = box sizes
!

    REAL(DP), INTENT(IN)::dx0, dy0, dz0
    INTEGER, INTENT(IN)::nx0, ny0, nz0
!-----------------------------------------------------------------------
    INTEGER :: ii, i1, i2, i3, kdum
    LOGICAL, PARAMETER :: tcoultest = .false.
    REAL(DP) :: charge
    REAL(DP), ALLOCATABLE :: rhotest(:), ctest(:)

    kxmax = 2*nx0; kymax = 2*ny0; kzmax = 2*nz0; ksmax = kxmax
    kdfull = nx0*ny0*nz0
    kdred = kxmax*kymax*kzmax
    kfft = ksmax; kfftx = kxmax; kffty = kymax; kfftz = kzmax
    kfft2 = kfft*2

    nx = nx0
    ny = ny0
    nz = nz0
    dx = dx0
    dy = dy0
    dz = dz0

    ALLOCATE (xt2(kxmax), yt2(kymax), zt2(kzmax))
    ALLOCATE (ikm(kxmax, kymax))


! start NFFTW3 BLOCK
    ALLOCATE (fftax(kxmax), fftay(kymax), fftb(kzmax, kxmax))
    pffta = fftw_alloc_complex(int(kxmax*kymax*kzmax, C_SIZE_T))
    CALL c_f_pointer(pffta, ffta, [kxmax, kymax, kzmax])
    pakv2 = fftw_alloc_complex(int(kxmax*kymax*kzmax, C_SIZE_T))
    CALL c_f_pointer(pakv2, akv2, [kxmax, kymax, kzmax])
    prffta = fftw_alloc_real(int(kxmax*kymax*kzmax, C_SIZE_T))
    CALL c_f_pointer(prffta, rffta, [kxmax, kymax, kzmax])
    prfftb = fftw_alloc_real(int(kxmax*kymax*kzmax, C_SIZE_T))
    CALL c_f_pointer(prfftb, rfftb, [kxmax, kymax, kzmax])
!ALLOCATE(ffta(kxmax,kymax,kzmax),akv2(kxmax,kymax,kzmax))
!ALLOCATE(rffta(kxmax,kymax,kzmax))
#if(omp)
    CALL dfftw_init_threads(iret)
    CALL dfftw_plan_with_nthreads(numthr)
    WRITE (*, *) ' init Coul FFTW threads: iret=', iret, ', nr. of threads=', numthr
#endif
    IF (nini == 0) THEN
#if(fftwnomkl)
      wisdomtest = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw_coul.dat'//C_NULL_CHAR)
#endif
      IF (wisdomtest == 0) THEN
        wisdomtest = fftw_import_system_wisdom()
        IF (wisdomtest == 0) THEN
          WRITE (6, *) 'wisdom_fftw.dat not found, creating it'
          WRITE (7, *) 'wisdom_fftw.dat not found, creating it'
        ELSE
          WRITE (*, *) 'Coulex: wisdom from system'
        END IF
      END IF
      pforw = fftw_plan_dft_r2c_3d(kzmax, kymax, kxmax, rffta, ffta, FFTW_planflag)
      pback = fftw_plan_dft_c2r_3d(kzmax, kymax, kxmax, ffta, rfftb, FFTW_planflag)
! pforwc=fftw_plan_dft_3d(kzmax,kymax,kxmax,ffta,ffta,FFTW_FORWARD,FFTW_planflag)
      nini = kxmax*kymax*kzmax
      WRITE (*, *) ' Coul-Solv initialized nini=', nini, kxmax, kymax, kzmax
    ELSE IF (nini /= kxmax*kymax*kzmax) THEN
      WRITE (*, *) ' nini,nx2,ny2,nz2=', nini, nx2, ny2, nz2
      STOP ' nx2, ny2 or/and nz2 in four3d not as initialized!'
    END IF
#if(fftwnomkl)
    wisdomtest = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw_coul.dat'//C_NULL_CHAR)
#endif
    IF (wisdomtest == 0) THEN
      WRITE (*, *) ' export wisdom_fftw_coul.dat failed'
    ELSE
      WRITE (*, *) ' export wisdom_fftw_coul.dat successfull'
    END IF
    CALL fftw_forget_wisdom
! END FFTW3 BLOCK

! CALL input routine fftinp, which initializes the grid and fft table
    CALL fftinp

! test section
    IF (tcoultest) THEN
      kdum = nx*ny*nz
      ALLOCATE (rhotest(nx*ny*nz), ctest(nx*ny*nz))
      rhotest = 0D0
      ii = 0
      DO i3 = 1, nz; DO i2 = 1, ny; DO i1 = 1, nx
            ii = ii + 1
            IF (i3 == nz/2 .AND. i2 == ny/2 .AND. i1 == nx/2) &
              rhotest(ii) = 1D0/(dx*dy*dz)
          END DO; END DO; END DO
      charge = SUM(rhotest)*(dx*dy*dz)
      WRITE (*, *) '# test Coulomb for point charge:', charge
! CALL falr(rhotest,ctest,kdum)
      CALL solv_poisson_e(rhotest, ctest, kdum)
      ii = 0
      DO i3 = 1, nz; DO i2 = 1, ny; DO i1 = 1, nx
            ii = ii + 1
            IF (i3 == nz/2 .AND. i2 == ny/2) &
              WRITE (*, *) (i1 - nx/2)*dx, rhotest(ii),ctest(ii) ! *(i1-nx/2)*dx/2D0
          END DO; END DO; END DO
      STOP "Coulomb test finished"
      DEALLOCATE (rhotest, ctest)
    END IF

    RETURN
  END SUBROUTINE init_coul_e

!-----fftinp------------------------------------------------------------

  SUBROUTINE fftinp
    IMPLICIT NONE

! Initializes work tables for FFT

! Grid parameters nx,ny,nz,dx,dy,dz,ecut must have been
! initialized before !

!-----------------------------------------------------------------------

    INTEGER :: ii, i1, i2, i3, ind, ikzero
    REAL(DP) :: ak2, xx1, xx2, xy1, xy2, xz1, xz2
    REAL(DP) :: factor
! initialize grid in coordinate space

    nx1 = nx + 1
    ny1 = ny + 1
    nz1 = nz + 1
    nxi = nx + nx
    nyi = ny + ny
    nzi = nz + nz
    nxy1 = nxi*nyi
    nxyz = nxi*nyi*nzi
    nkxyz = nxi*nyi*nzi

! grid lengths must match with parameters in incs

    IF (kxmax < nxi) THEN
      WRITE (6, '(a)') ' ERROR: parameter kxmax too small'
      STOP ' error in parameter: KXMAX in COULEX too small'
    ELSE IF (kymax < nyi) THEN
      WRITE (6, '(a)') ' ERROR: parameter kymax too small'
      STOP ' error in parameter: KYMAX in COULEX too small'
    ELSE IF (kzmax < nzi) THEN
      WRITE (6, '(a)') ' ERROR: parameter kzmax too small'
      STOP ' error in parameter: KZMAX in COULEX too small'
    END IF

! initialize grid in Fourier space

    dkx = pi/(dx*REAL(nx, DP))
    dky = pi/(dy*REAL(ny, DP))
    dkz = pi/(dz*REAL(nz, DP))

    dxsp = dx*dy*dz
    dksp = dkx*dky*dkz
    WRITE (*, *) ' dkx,dky,dkz,dksp=', dkx, dky, dkz, dksp

    grnorm = SQRT(dxsp/dksp)
    fnorm = 1.0D0/SQRT(REAL(nx*ny*nz, DP))
    nxk = nx1

! built Greens function in Fourier space
! by Fourier transformation from REAL space

    ikzero = nxy1*(nz - 1) + nxi*(ny - 1) + nx
    WRITE (*, *) ' nzi,nyi,nxi,nx,ny,nz,ikzero=', nzi, nyi, nxi, nx, ny, nz, ikzero

! start FFTW3 BLOCK
    factor = e2*(dx*dy*dz)/(nx*ny*nz*8D0)
    ii = 0
    DO i3 = 1, nzi
      IF (i3 <= nz) THEN
        xz1 = (i3 - 1)*dz
      ELSE
        xz1 = (i3 - nzi - 1)*dz
      END IF
      xz2 = xz1*xz1
      DO i2 = 1, nyi
        IF (i2 <= ny) THEN
          xy1 = (i2 - 1)*dy
        ELSE
          xy1 = (i2 - nyi - 1)*dy
        END IF
        xy2 = xy1*xy1
        DO i1 = 1, nxi
          IF (i1 <= nx) THEN
            xx1 = (i1 - 1)*dx
          ELSE
            xx1 = (i1 - nxi - 1)*dx
          END IF
          xx2 = xx1*xx1
          ak2 = xx2 + xy2 + xz2
          ii = ii + 1
          IF (ak2 > dx**2/10D0) THEN
            akv2(i1, i2, i3) = factor/SQRT(ak2)
          ELSE
            akv2(i1, i2, i3) = factor*2.34D0*1.19003868D0*(dx*dy*dz)**(-1D0/3D0)
          END IF
          rffta(i1, i2, i3) = akv2(i1, i2, i3)
        END DO
      END DO
    END DO
    nksp = ii
    CALL fftw_execute_dft_r2c(pforw, rffta, ffta)
    DO i3 = 1, nzi
      DO i2 = 1, nyi
        DO i1 = 1, nxi
          akv2(i1, i2, i3) = ffta(i1, i2, i3)
        END DO
      END DO
    END DO
    WRITE (6, *) '1/k**2 along x'
    DO i1 = 1, nxi
      WRITE (6, '(f8.2,2(1pg13.5))') (i1 - nx)*dx, akv2(i1, 1, 1)
    END DO
    WRITE (6, *) '1/k**2 along z'
    DO i3 = 1, nzi
      WRITE (6, '(f8.2,2(1pg13.5))') (i3 - nz)*dz, akv2(1, 1, i3)
    END DO
    rffta = 0D0
! END FFTW3 BLOCK



    RETURN
  END SUBROUTINE fftinp

!start FFTW3 BLOCK
!-------------------------------------------------------------------

  SUBROUTINE solv_poisson_e(rhoinp, chpfalr, kdf)
    IMPLICIT NONE

! Coulomb solver using FFTW
!
! Input:
! rhoinp = charge density
! kdf = size of array
! Output:
! chpfalr = resulting Coulomb potential

    INTEGER, INTENT(IN) :: kdf
    REAL(DP), INTENT(IN) :: rhoinp(kdf)
    REAL(DP), INTENT(OUT) :: chpfalr(kdf)

    LOGICAL, PARAMETER :: tmultdirect = .FALSE.
    LOGICAL, PARAMETER :: ttimestop = .FALSE.
    INTEGER :: it1, it2
    INTEGER :: i0, i1, i2, i3
    REAL(DP) ::factor

! map density to 2**3 larger grid, immediately on FFTW3 work space

!rffta = 0D0
    i0 = 0
    DO i3 = 1, nz
      DO i2 = 1, ny
        DO i1 = 1, nx
          i0 = i0 + 1
          rffta(i1, i2, i3) = rhoinp(i0)
        END DO
      END DO
    END DO

! FFT forward, multiply with Green's function, FFT backward

    CALL fftw_execute_dft_r2c(pforw, rffta, ffta)

    IF (ttimestop) CALL system_clock(it1)
    IF (tmultdirect) THEN
      CALL mult_direct3(ffta, akv2)
    ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i1,i2,i3) SCHEDULE(STATIC)
      DO i3 = 1, kzmax
        DO i2 = 1, kymax
          DO i1 = 1, kxmax
            ffta(i1, i2, i3) = akv2(i1, i2, i3)*ffta(i1, i2, i3)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF
    IF (ttimestop) THEN
      CALL system_clock(it2)
      WRITE (*, '(a,f8.3)') 'systime5=', (it2 - it1)*systime_factor
    END IF
!ffta = akv2*ffta

    CALL fftw_execute_dft_c2r(pback, ffta, rfftb)

! map back to standard grid, augment with normalization factor

!facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(kxmax*kymax*kzmax,DP)) * 2D0
    i0 = 0
    DO i3 = 1, nz
      DO i2 = 1, ny
        DO i1 = 1, nx
          i0 = i0 + 1
          chpfalr(i0) = rfftb(i1, i2, i3) ! *facnr ! all scaling in 'akv2'
        END DO
      END DO
    END DO

  END SUBROUTINE solv_poisson_e

  SUBROUTINE mult_direct3(asub, bsub)

    COMPLEX(8), INTENT(IN) :: bsub(kxmax, kymax, kzmax)
    COMPLEX(8), INTENT(IN OUT) :: asub(kxmax, kymax, kzmax)

    INTEGER :: i1, i2, i3, it9

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i1,i2,i3) SCHEDULE(STATIC)
    DO i3 = 1, kzmax
      DO i2 = 1, kymax
        DO i1 = 1, kxmax
          asub(i1, i2, i3) = bsub(i1, i2, i3)*asub(i1, i2, i3)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    RETURN

  END SUBROUTINE mult_direct3

  SUBROUTINE coulsolv_end()

! Epilogue for Coulomb solver

    USE FFTW

    CALL fftw_destroy_plan(pforwx)
    CALL fftw_destroy_plan(pforwy)
    CALL fftw_destroy_plan(pforwz)
    CALL fftw_destroy_plan(pbackx)
    CALL fftw_destroy_plan(pbacky)
    CALL fftw_destroy_plan(pbackz)

    DEALLOCATE (xt2, yt2, zt2)
    DEALLOCATE (fftax, fftay, fftb)
    DEALLOCATE (ikm)

    RETURN
  END SUBROUTINE coulsolv_end

! END of NETLIB BLOCK

END MODULE coulsolv_e
