
MODULE coulsolv
  USE params, ONLY: DP
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
  INTEGER, PRIVATE :: nxhigh, nxlow, nyhigh, nylow, nzhigh, nzlow, nx2, ny2, nz2
!COMMON /xgrid/ xval,yval,zval,xt2,yt2,zt2,dx,dy,dz,dxsp,grnorm,fnorm, &
! nx,ny,nz,nx1,ny1,nz1, nxr,nxi,nyr,nyi,nzr,nzi,nxy1,nxyz, &
! nxhigh,nxlow,nyhigh,nylow,nzhigh,nzlow

  REAL(DP), ALLOCATABLE, PRIVATE :: akv2r(:)
  COMPLEX(DP), ALLOCATABLE, PRIVATE :: akv2c(:)
  INTEGER, ALLOCATABLE, PRIVATE :: ikm(:, :)
  REAL(DP), PRIVATE :: dkx, dky, dkz, akmax, dksp, ecut
  INTEGER, PRIVATE :: nxk, nxklo, nxkhi, nksp, nkxyz
!COMMON /kgrid/ akv2r,akv2i, dkx,dky,dkz,akmax,dksp, &
! ikm,nxk,nxklo,nxkhi,nksp,nkxyz,ecut

  REAL(DP), ALLOCATABLE, PRIVATE :: wrkx(:), wrky(:), wrkz(:)
  REAL(DP), ALLOCATABLE, PRIVATE :: wsavex(:), wsavey(:), wsavez(:)
  INTEGER, ALLOCATABLE, PRIVATE :: ifacx(:), ifacy(:), ifacz(:)
!COMMON /fftini/wrkx,wrky,wrkz,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz

! include BLOCK: option
  REAL(DP), PARAMETER, PRIVATE :: zero = 0D0
  REAL(DP), PARAMETER, PRIVATE :: pi = 3.141592653589793D0

  COMPLEX(DP), ALLOCATABLE, PRIVATE :: fftax(:), fftay(:), fftb(:, :)
!COMMON /fftcom/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax)

CONTAINS

  SUBROUTINE init_coul(dx0, dy0, dz0, nx0, ny0, nz0)

!-----------------------------------------------------------------------

! READ grid parameters from FILE or simply initialize them
! note that the Coulomb solver doubles the grid internally
    kxmax = 2*nx0; kymax = 2*ny0; kzmax = 2*nz0; ksmax = kxmax
    kdfull = nx0*ny0*nz0
    kdred = kxmax*kymax*kzmax
    kfft = 2*ksmax; kfftx = kxmax; kffty = kymax; kfftz = kzmax
    kfft2 = kfft*2 + 1

    nx = nx0 !/2
    ny = ny0 !/2
    nz = nz0 !/2
    dx = dx0
    dy = dy0
    dz = dz0

    ALLOCATE (xval(kxmax), yval(kymax), zval(kzmax))
    ALLOCATE (xt2(kxmax), yt2(kymax), zt2(kzmax))
    ALLOCATE (wrkx(kfft2), wrky(kfft2), wrkz(kfft2))
    ALLOCATE (wsavex(kfft2), wsavey(kfft2), wsavez(kfft2))
    ALLOCATE (ifacx(kfft2), ifacy(kfft2), ifacz(kfft2))
    ALLOCATE (fftax(kxmax), fftay(kymax), fftb(kzmax, kxmax))
    ALLOCATE (akv2r(kdred), akv2c(kdred))
    ALLOCATE (ikm(kxmax, kymax))

! CALL input routine fftinp, which initializes the grid and fft tabl

    CALL fftinp

    RETURN
  END SUBROUTINE init_coul

!-----fftinp------------------------------------------------------------

  SUBROUTINE fftinp
    IMPLICIT REAL(DP) (A - H, O - Z)
    REAL(DP), ALLOCATABLE :: akvpr(:)

! initializes work tables for FFT

! grid parameters nx,ny,nz,dx,dy,dz,ecut must have been READ or
! initialized before !

!-----------------------------------------------------------------------

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
    nx2 = nxr
    ny2 = nyr
    nz2 = nzr
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

! initialize grid IN Fourier space

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
        END DO
      END DO
    END DO
    nksp = ii

    WRITE (*, *) ' akv2r ready'

    CALL rftf2(akv2r(1), akv2c(1))
    ALLOCATE (akvpr(kdred))
    CALL rfftback2(akv2c(1), akvpr(1))

! akvpr = REAL(akv2c)
! CALL prifld(akvpr,'k**2 for Co')
! akvpr = AIMAG(akv2c)
! CALL prifld(akvpr,'k**2 for Co')
    nxyf = nx2*ny2
    nyf = nx2
    ind = nz*nxyf + ny*nyf
!WRITE(6,*) 'k**2 along x'
    WRITE (6, *) ' norms:', SUM(akv2r**2), SUM(akvpr**2)
    WRITE (6, *) 'test repro'
    DO i1 = 1, nx2
      WRITE (6, '(f8.2,2(1pg13.5))') (i1 - nx)*dx, akv2r(i1 + ind), akvpr(i1 + ind)
    END DO
    DEALLOCATE (akvpr)

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
!WRITE(*,*) 'FALR entered'

    CALL rhofld(rhoinp, rhokr, rhoki)

!WRITE(*,*) 'after RHOFLD'

! CALL coufou, which CONTAINS the fcs PROCEDURE.

    CALL coufou2(rhokr)
!WRITE(*,*) 'after COUFOU'

! CALL a routine written by you which outputs the results of the fcs
! and maybe some other things TO an output FILE or the screen.

    CALL RESULT(chpfalr, rhokr)
!WRITE(*,*) 'after RESULT'

    DEALLOCATE (rhokr, rhoki)

  END SUBROUTINE falr

!-----rhofld------------------------------------------------------------

  SUBROUTINE rhofld(rhoinp, rhokr, rhoki)
    IMPLICIT REAL(DP) (A - H, O - Z)

! copy density on COMPLEX array of DOUBLE extension IN x,y,z

    REAL(DP), INTENT(IN) :: rhoinp(kdfull)
    REAL(DP), INTENT(OUT) :: rhokr(kdred)
    REAL(DP), INTENT(OUT) :: rhoki(kdred)

    rhokr = 0D0
    rhoki = 0D0
    i0 = 0
    DO i3 = 1, nz
      DO i2 = 1, ny
        ii = (i3 - 1)*nxi*nyi + (i2 - 1)*nxi
        DO i1 = 1, nx
          ii = ii + 1
          i0 = i0 + 1
          rhokr(ii) = rhoinp(i0)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE rhofld

!-----RESULT------------------------------------------------------------

  SUBROUTINE RESULT(chpfalr, rhokr)

    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(OUT) :: chpfalr(kdfull)
    REAL(DP), INTENT(IN) :: rhokr(kdred)

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

  SUBROUTINE coufou2(rhokr)

    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN OUT) :: rhokr(kdred)
    COMPLEX(DP), ALLOCATABLE :: rhok(:)

    LOGICAL, PARAMETER :: tprint = .false.
    LOGICAL, PARAMETER :: rqplot = .false.

!------------------------------------------------------------------------------

    ALLOCATE (rhok(kdred))
! Fourier transformation of the density

    CALL rftf2(rhokr, rhok)

! calculation of the Coulomb field (writing on the density field)

    rhok(:) = akv2c(:)*rhok(:)

! Fourier back transformation

    CALL rfftback2(rhok, rhokr)

    DEALLOCATE (rhok)

    RETURN

  END SUBROUTINE coufou2

! ******************************

  SUBROUTINE rftf2(q1, q2)

! ******************************

    USE params
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN OUT) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax), &
! wrkx(kfft2),wrky(kfft2),wrkz(kfft2), &
! wsavex(kfft2),wsavey(kfft2),wsavez(kfft2), &
! ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

    INTEGER, SAVE :: nxini = 0, nyini = 0, nzini = 0 ! flag for initialization

! check initialization

! nx=nx2/2
! ny=ny2/2
! nz=nz2/2
    tnorm = 1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2, DP))
!WRITE(*,*) 'TNORM=',tnorm
    IF (nxini == 0) THEN
      CALL dcfti1(nx2, wsavex, ifacx)
      nxini = nx2
! WRITE(6,'(a)') ' x-fft initialized '
    ELSE IF (nxini /= nx2) THEN
      STOP ' nx2 IN four3d not as initialized!'
    END IF
    IF (nyini == 0) THEN
      CALL dcfti1(ny2, wsavey, ifacy)
      nyini = ny2
! WRITE(6,'(a)') ' y-fft initialized '
    ELSE IF (nyini /= ny2) THEN
      STOP ' ny2 IN four3d not as initialized!'
    END IF
    IF (nzini == 0) THEN
      CALL dcfti1(nz2, wsavez, ifacz)
      nzini = nz2
! WRITE(6,'(a)') ' z-fft initialized '
    ELSE IF (nzini /= nz2) THEN
      STOP ' nz2 IN four3d not as initialized!'
    END IF

! transformation IN x-direction

    nxyf = nx2*ny2
    nyf = nx2
    DO i3 = 1, nz2
      DO i2 = 1, ny2
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftax(MOD(i1 + nx, nx2) + 1) = CMPLX(q1(ind), 0D0, DP) ! copy TO workspace
        END DO
        CALL dcftf1(nx2, fftax, wrkx, wsavex, ifacx) ! basic fft
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = fftax(i1) ! copy back IN strange order
        END DO
      END DO
    END DO

! transformation IN y-direction

    DO i3 = 1, nz2
      DO i1 = 1, nx2
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftay(MOD(i2 + ny, ny2) + 1) = q2(ind)
        END DO
        CALL dcftf1(ny2, fftay, wrky, wsavey, ifacy)
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = fftay(i2)
        END DO
      END DO
    END DO

!oldc
!oldc transformation IN z-direction
!oldc
!old DO i2=1,ny2
!old DO i3=1,nz2
!old DO i1=1,nx2
!old ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old fftb(mod(i3+nz,nz2)+1,i1) = q2(ind)
!old END DO
!old END DO
!old DO i1=1,nx2
!old CALL dcftf1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
!old END DO
!old DO i3=1,nz2
!old DO i1=1,nx2
!old ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old q2(ind)= fftb(i3,i1)*tnorm
!old END DO
!old END DO
!old END DO

! transformation IN z-direction

    DO i2 = 1, ny2
      DO i3 = 1, nz2
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftb(MOD(i3 + nz, nz2) + 1, i1) = q2(ind)
        END DO
      END DO
      DO i1 = 1, nx2
        CALL dcftf1(nz2, fftb(1, i1), wrkz, wsavez, ifacz)
      END DO
      DO i3 = 1, nz2
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = fftb(i3, i1)*tnorm
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE rftf2

! ******************************

  SUBROUTINE rfftback2(q1, q3)
!SUBROUTINE rfftback(q1,q2)

! ******************************

    USE params
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
!REAL(DP), INTENT(OUT) :: q2(kdfull2)
    REAL(DP), INTENT(OUT) :: q3(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax), &
! wrkx(kfft2),wrky(kfft2),wrkz(kfft2), &
! wsavex(kfft2),wsavey(kfft2),wsavez(kfft2), &
! ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

    COMPLEX(DP), ALLOCATABLE :: q2(:)

! DATA nxini,nyini,nzini/0,0,0/ ! flag for initialization
    nxyf = nx2*ny2
    nyf = nx2
! nx=nx2/2
! ny=ny2/2
! nz=nz2/2
    facnr = SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2, DP))
!WRITE(*,*) 'FACNR=',facnr
!oldc
!oldc transformation IN z-direction
!oldc
!old DO i2=1,ny2
!old DO i1=1,nx2
!old DO i3=1,nz2
!old ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old ffta(i3) = q1(ind)*facnr
!old END DO
!old CALL dcftb1 (nz2,ffta,wrkz,wsavez,ifacz) ! basic fft
!old DO i3=1,nz2 ! copy back
!old ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old q2(ind)= ffta(mod(i3+nz,nz2)+1)
!old END DO
!old END DO
!old END DO

    ALLOCATE (q2(kdfull2))

! transformation IN z-direction

    DO i2 = 1, ny2
      DO i3 = 1, nz2
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftb(i3, i1) = q1(ind)*facnr
        END DO
      END DO
      DO i1 = 1, nx2
        CALL dcftb1(nz2, fftb(1, i1), wrkz, wsavez, ifacz) ! basic fft
      END DO
      DO i3 = 1, nz2 ! copy back
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = fftb(MOD(i3 + nz, nz2) + 1, i1)
        END DO
      END DO
    END DO

! transformation IN y-direction

    DO i3 = 1, nz2
      DO i1 = 1, nx2
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftay(i2) = q2(ind)
        END DO
        CALL dcftb1(ny2, fftay, wrky, wsavey, ifacy)
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = fftay(MOD(i2 + ny, ny2) + 1)
        END DO
      END DO
    END DO

! transformation IN x-direction

    DO i3 = 1, nz2
      DO i2 = 1, ny2
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftax(i1) = q2(ind)
        END DO
        CALL dcftb1(nx2, fftax, wrkx, wsavex, ifacx)
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
! q2(ind)= REAL(fftax(MOD(i1+nx,nx2)+1),DP)
          q3(ind) = REAL(fftax(MOD(i1 + nx, nx2) + 1), DP)
        END DO
      END DO
    END DO

    DEALLOCATE (q2)

    RETURN
  END SUBROUTINE rfftback2

END MODULE coulsolv
