
MODULE kinetic
  USE, INTRINSIC :: iso_c_binding
  USE params, ONLY: DP
  IMPLICIT REAL(DP) (A - H, O - Z)

  SAVE
! arrays for kinetic energy and electronic propagation
! kinetic energy coefficients IN strange ordered fourier space
! akv = fouier-field for 0.5*k^2
! ak = fourier-field for exp(i*dt*(h^2/2m)*k^2)
  COMPLEX(DP), ALLOCATABLE :: ak(:)
  REAL(DP), ALLOCATABLE :: akv(:)
  COMPLEX(DP), ALLOCATABLE, PRIVATE :: akpropx(:), akpropy(:), akpropz(:)
  COMPLEX(DP), PARAMETER, PRIVATE :: eye = (0D0, 1D0)
  REAL(DP), PARAMETER, PRIVATE :: PI = 3.141592653589793D0

  INTEGER, PRIVATE :: kfft, kfft2, kdfull2

  COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: fftax(:), fftay(:), fftaz(:), fftb(:, :), ffta(:, :, :)
  TYPE(C_PTR), PRIVATE :: pforwx, pforwy, pforwz, pforwz1, pbackx, pbacky, pbackz, pbackz1
  TYPE(C_PTR), PRIVATE :: pforw, pback
  INTEGER(C_INT), PRIVATE :: wisdomtest

CONTAINS
!-----init_grid_fft-----------------------------------------------------

  SUBROUTINE init_grid_fft(dx0, dy0, dz0, nx0, ny0, nz0, dt1, h2m)

! initialize details for FFT

    USE FFTW

    REAL(DP) :: dt1, h2m
    INTEGER, SAVE :: nxini = 0, nyini = 0, nzini = 0, nini = 0 ! flag for initialization

    nx2 = nx0; ny2 = ny0; nz2 = nz0
    kxmax = nx0; kymax = ny0; kzmax = nz0
    nx = nx2/2; ny = ny2/2; nz = nz2/2

    kfft = 2*kxmax
    kfft2 = kfft*2 + 1
    kdfull2 = kxmax*kymax*kzmax

    dkx = pi/(dx0*nx)
    dky = pi/(dy0*ny)
    dkz = pi/(dz0*nz)

    ALLOCATE (ak(kdfull2), akv(kdfull2))
    ALLOCATE (akpropx(kxmax), akpropy(kymax), akpropz(kzmax))
    ALLOCATE (fftax(kxmax), fftay(kymax), fftb(kzmax, kxmax), ffta(kxmax, kymax, kzmax))

    WRITE (7, *) 'h bar squared over two m electron', h2m
    WRITE (7, *) ' testprint EYE=', eye
    WRITE (7, *) ' dkx,dky,dkz=', dkx, dky, dkz
    WRITE (7, *) ' testprint: nx2,ny2,nz2=', nx2, ny2, nz2
    WRITE (7, *) ' testprint: kdfull2,kfft2=', kdfull2, kfft2

! prepare k**2 and kinetic propagation factor IN 3D momentum space

    ind = 0
    DO i3 = 1, nz2
      IF (i3 >= (nz + 1)) THEN
        zkz = (i3 - nz2 - 1)*dkz
      ELSE
        zkz = (i3 - 1)*dkz
      END IF
      DO i2 = 1, ny2
        IF (i2 >= (ny + 1)) THEN
          zky = (i2 - ny2 - 1)*dky
        ELSE
          zky = (i2 - 1)*dky
        END IF
        DO i1 = 1, nx2
          IF (i1 >= (nx + 1)) THEN
            zkx = (i1 - nx2 - 1)*dkx
          ELSE
            zkx = (i1 - 1)*dkx
          END IF
          ind = ind + 1
          ak(ind) = EXP(-eye*dt1*(zkx**2 + zky**2 + zkz**2)*h2m)
          akv(ind) = (zkx**2 + zky**2 + zkz**2)*h2m
        END DO
      END DO
    END DO

! prepare kinetic propagation factors IN 1D momentum spaces

    DO i3 = 1, nz2
      IF (i3 >= (nz + 1)) THEN
        zkz = (i3 - nz2 - 1)*dkz
      ELSE
        zkz = (i3 - 1)*dkz
      END IF
      akpropz(i3) = EXP(-eye*dt1*zkz**2*h2m)
    END DO

    DO i2 = 1, ny2
      IF (i2 >= (ny + 1)) THEN
        zky = (i2 - ny2 - 1)*dky
      ELSE
        zky = (i2 - 1)*dky
      END IF
      akpropy(i2) = EXP(-eye*dt1*zky**2*h2m)
    END DO

    DO i1 = 1, nx2
      IF (i1 >= (nx + 1)) THEN
        zkx = (i1 - nx2 - 1)*dkx
      ELSE
        zkx = (i1 - 1)*dkx
      END IF
      akpropx(i1) = EXP(-eye*dt1*zkx**2*h2m)
    END DO

    IF (nini == 0) THEN
      wisdomtest = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'wisdom_fftw.dat not found, creating it'
        WRITE (7, *) 'wisdom_fftw.dat not found, creating it'
      END IF
      pforw = fftw_plan_dft_3d(nz2, ny2, nx2, ffta, ffta, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      nini = nx2*ny2*nz2
      pback = fftw_plan_dft_3d(nz2, ny2, nx2, ffta, ffta, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      nini = nx2*ny2*nz2
      wisdomtest = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'Error exporting wisdom TO FILE wisdom_fftw.dat'
        WRITE (7, *) 'Error exporting wisdom TO FILE wisdom_fftw.dat'
      END IF
      nini = nx2*ny2*nz2
    ELSE IF (nini /= nx2*ny2*nz2) THEN
      STOP ' nx2, ny2 or/and nz2 IN four3d not as initialized!'
    END IF
    IF (nxini == 0) THEN
      wisdomtest = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'wisdom_fftw.dat not found, creating it'
        WRITE (7, *) 'wisdom_fftw.dat not found, creating it'
      END IF
      pforwx = fftw_plan_dft_1d(nx2, fftax, fftax, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      pbackx = fftw_plan_dft_1d(nx2, fftax, fftax, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      nxini = nx2
! WRITE(6,'(a)') ' x-fft initialized '
    ELSE IF (nxini /= nx2) THEN
      STOP ' nx2 IN four3d not as initialized!'
    END IF
    IF (nyini == 0) THEN
      pforwy = fftw_plan_dft_1d(ny2, fftay, fftay, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      pbacky = fftw_plan_dft_1d(ny2, fftay, fftay, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      nyini = ny2
! WRITE(6,'(a)') ' y-fft initialized '
    ELSE IF (nyini /= ny2) THEN
      STOP ' ny2 IN four3d not as initialized!'
    END IF
    IF (nzini == 0) THEN
      pforwz = fftw_plan_dft_1d(nz2, fftb(1, i1), fftb(1, i1), FFTW_FORWARD, FFTW_EXHAUSTIVE)
      pbackz = fftw_plan_dft_1d(nz2, fftb(1, i1), fftb(1, i1), FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      pforwz1 = fftw_plan_dft_1d(nz2, fftaz, fftaz, FFTW_FORWARD, FFTW_EXHAUSTIVE)
      pbackz1 = fftw_plan_dft_1d(nz2, fftaz, fftaz, FFTW_BACKWARD, FFTW_EXHAUSTIVE)
      nzini = nz2
      wisdomtest = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'Error exporting wisdom TO FILE wisdom_fftw.dat'
        WRITE (7, *) 'Error exporting wisdom TO FILE wisdom_fftw.dat'
      END IF
! WRITE(6,'(a)') ' z-fft initialized '
    ELSE IF (nzini /= nz2) THEN
      STOP ' nz2 IN four3d not as initialized!'
    END IF
  END SUBROUTINE init_grid_fft

! ******************************

  SUBROUTINE kinprop(q1, q2)

! ******************************

! propagation with exp(-i*dt*e_kin)

    USE params
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN OUT) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)

! check initialization

!tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

!fnorm=1D0/(nx2*ny2*nz2)

! propagation IN x-direction

    xfnorm = 1D0/nx2
    DO i3 = 1, nz2
      DO i2 = 1, ny2
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftax(MOD(i1 + nx, nx2) + 1) = q1(ind) ! copy TO workspace
        END DO
        CALL fftw_execute_dft(pforwx, fftax, fftax)
        DO i1 = 1, nx2
          fftax(i1) = akpropx(i1)*fftax(i1)
        END DO
        CALL fftw_execute_dft(pbackx, fftax, fftax)
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = fftax(MOD(i1 + nx, nx2) + 1)*xfnorm
        END DO
      END DO
    END DO

! transformation IN y-direction

    yfnorm = 1D0/ny2
    DO i3 = 1, nz2
      DO i1 = 1, nx2
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftay(MOD(i2 + ny, ny2) + 1) = q2(ind)
        END DO
        CALL fftw_execute_dft(pforwy, fftay, fftay)
        DO i2 = 1, ny2
          fftay(i2) = akpropy(i2)*fftay(i2)
        END DO
        CALL fftw_execute_dft(pbacky, fftay, fftay)
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = fftay(MOD(i2 + ny, ny2) + 1)*yfnorm
        END DO
      END DO
    END DO

! propagation IN z-direction

    zfnorm = 1D0/nz2
    DO i2 = 1, ny2
      DO i3 = 1, nz2
        i3m = MOD(i3 + nz, nz2) + 1
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftb(i3m, i1) = q2(ind)
        END DO
      END DO
      DO i1 = 1, nx2
        CALL fftw_execute_dft(pforwz, fftb(1, i1), fftb(1, i1))
        DO i3 = 1, nz2
          fftb(i3, i1) = akpropz(i3)*fftb(i3, i1)
        END DO
        CALL fftw_execute_dft(pbackz, fftb(1, i1), fftb(1, i1))
      END DO
      DO i3 = 1, nz2
        i3m = MOD(i3 + nz, nz2) + 1
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q1(ind) = fftb(i3m, i1)*zfnorm
        END DO
      END DO
    END DO

!IF (nini==0) THEN
! CALL dfftw_plan_dft_3d(pkinforw,nx2,ny2,nz2,ffta,ffta, &
! FFTW_FORWARD,FFTW_PATIENT)
! CALL dfftw_plan_dft_3d(pkinback,nx2,ny2,nz2,ffta,ffta, &
! FFTW_BACKWARD,FFTW_PATIENT)
! nini = nx2*ny2*nz2
! WRITE(6,'(a)') ' fft initialized '
!ELSE IF(nini /= nx2*ny2*nz2) THEN
! STOP ' nx2 IN four3d not as initialized!'
!END IF

!CALL copy1dto3d(q1,ffta,nx2,ny2,nz2)

!CALL dfftw_execute_dft(pkinforw,ffta,ffta)

!DO i3=1,nz2
! DO i2=1,ny2
! DO i1=1,nx2
! ind=(i3-1)*nxyf+(i2-1)*nyf+i1
! ffta(i1,i2,i3)=ffta(i1,i2,i3)*ak(ind)
! END DO
! END DO
!END DO

!CALL dfftw_execute_dft(pkinback,ffta,ffta)

!CALL secopy3dto1d(ffta,q2,fnorm,nx2,ny2,nz2)

    RETURN
  END SUBROUTINE kinprop

! ******************************

  SUBROUTINE gradient(fin, gradfout, idirec)

! ******************************

! The gradient of COMPLEX field 'fin' IN direction
! 'idirec' (x =1, y=2, z=3).
! The fields are given IN Fourier space and the gradient is
! applied as product with 'kx', 'ky', 'kz'.
!

    USE params
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN OUT) :: fin(kdfull2)
    COMPLEX(DP), INTENT(IN OUT) :: gradfout(kdfull2)
    INTEGER, INTENT(IN) :: idirec

! ************************************************************

    IF (idirec == 1) THEN
! x-derivative

      dkx = pi/(dx*REAL(nx, DP))
      ind = 0
      DO i3 = 1, nz2
        DO i2 = 1, ny2
          DO i1 = 1, nx2
            IF (i1 >= (nx + 1)) THEN
              zkx = (i1 - nx2 - 1)*dkx
            ELSE
              zkx = (i1 - 1)*dkx
            END IF
            ind = ind + 1
            gradfout(ind) = eye*zkx*fin(ind)
          END DO
        END DO
      END DO

    ELSEIF (idirec == 2) THEN
! y-derivative

      dky = pi/(dy*REAL(ny, DP))
      ind = 0
      DO i3 = 1, nz2
        DO i2 = 1, ny2
          IF (i2 >= (ny + 1)) THEN
            zky = (i2 - ny2 - 1)*dky
          ELSE
            zky = (i2 - 1)*dky
          END IF
          DO i1 = 1, nx2
            ind = ind + 1
            gradfout(ind) = eye*zky*fin(ind)
          END DO
        END DO
      END DO

    ELSEIF (idirec == 3) THEN
! z-derivative

      ind = 0
      dkz = pi/(dz*REAL(nz, DP))
      DO i3 = 1, nz2
        IF (i3 >= (nz + 1)) THEN
          zkz = (i3 - nz2 - 1)*dkz
        ELSE
          zkz = (i3 - 1)*dkz
        END IF
        DO i2 = 1, ny2
          DO i1 = 1, nx2
            ind = ind + 1
            gradfout(ind) = eye*zkz*fin(ind)
          END DO
        END DO
      END DO

    ELSE
      STOP ' RGRADIENT called with invalid IDIREC'
    END IF

    RETURN
  END SUBROUTINE gradient

! ******************************

  SUBROUTINE xgradient_rspace(fin, gradfout)

! ******************************

! The gradient of the COMPLEX field 'fin' IN x-direction.
! The fields are given IN coordinate space. The gradient is
! evaluated IN k_x space.
!

    USE params
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN OUT) :: fin(kdfull2)
    COMPLEX(DP), INTENT(IN OUT) :: gradfout(kdfull2)

! ************************************************************

    dkx = pi/(dx*nx)
    DO i3 = 1, nz2
      DO i2 = 1, ny2
! forward transform along x
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftax(MOD(i1 + nx, nx2) + 1) = fin(ind) ! copy TO workspace
        END DO
        CALL fftw_execute_dft(pforwx, fftax, fftax) ! basic fft
! multiply by k_x factor IN k_x space
        DO i1 = 1, nx2
          IF (i1 >= (nx + 1)) THEN
            zkx = (i1 - nx2 - 1)*dkx
          ELSE
            zkx = (i1 - 1)*dkx
          END IF
          fftax(i1) = fftax(i1)*eye*zkx
        END DO
! backward transform along x
        CALL fftw_execute_dft(pbackx, fftax, fftax)
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          gradfout(ind) = fftax(MOD(i1 + nx, nx2) + 1)/nx2
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE xgradient_rspace

! ******************************

  SUBROUTINE ygradient_rspace(fin, gradfout)

! ******************************

! The gradient of the COMPLEX field 'fin' IN y-direction.
! The fields are given IN coordinate space. The gradient is
! evaluated IN k_y space.
!

    USE params
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN OUT) :: fin(kdfull2)
    COMPLEX(DP), INTENT(IN OUT) :: gradfout(kdfull2)

! ************************************************************

    dky = pi/(dy*ny)
    DO i3 = 1, nz2
      DO i1 = 1, nx2
! forward transform along y
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftay(MOD(i2 + ny, ny2) + 1) = fin(ind)
        END DO
        CALL fftw_execute_dft(pforwy, fftay, fftay)
! multiply by k_y factor IN k_y space
        DO i2 = 1, ny2
          IF (i2 >= (ny + 1)) THEN
            zky = (i2 - ny2 - 1)*dky
          ELSE
            zky = (i2 - 1)*dky
          END IF
          fftay(i2) = fftay(i2)*eye*zky
        END DO
! backward transform along y
        CALL fftw_execute_dft(pbacky, fftay, fftay)
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          gradfout(ind) = fftay(MOD(i2 + ny, ny2) + 1)/ny2
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE ygradient_rspace

! ******************************

  SUBROUTINE zgradient_rspace(fin, gradfout)

! ******************************

! The gradient of the COMPLEX field 'fin' IN z-direction.
! The fields are given IN coordinate space. The gradient is
! evaluated IN k_z space.
!

    USE params
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN OUT) :: fin(kdfull2)
    COMPLEX(DP), INTENT(IN OUT) :: gradfout(kdfull2)

! ************************************************************

!ALLOCATE(fftaz(nz2))

    dkz = pi/(dz*nz)
    DO i2 = 1, ny2
      DO i1 = 1, nx2
! forward transform along z
        DO i3 = 1, nz2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          i3m = MOD(i3 + nz, nz2) + 1
          fftaz(i3m) = fin(ind)
        END DO
        CALL fftw_execute_dft(pforwz, fftaz, fftaz)
! multiply by k_z factor IN k_z space
        DO i3 = 1, nz2
          IF (i3 >= (nz + 1)) THEN
            zkz = (i3 - nz2 - 1)*dkz
          ELSE
            zkz = (i3 - 1)*dkz
          END IF
          fftaz(i3) = fftaz(i3)*eye*zkz
        END DO
! backward transform along z
        CALL fftw_execute_dft(pforwz, fftaz, fftaz)
        DO i3 = 1, nz2 ! copy back
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          gradfout(ind) = fftaz(MOD(i3 + nz, nz2) + 1)/nz2
        END DO
!
      END DO
    END DO

!DEALLOCATE(fftaz)

    RETURN

  END SUBROUTINE zgradient_rspace

! ******************************

  SUBROUTINE fftf(q1, q2)

! ******************************

    USE params
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax), &
! wrkx(kfft2),wrky(kfft2),wrkz(kfft2), &
! wsavex(kfft2),wsavey(kfft2),wsavez(kfft2), &
! ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)
    INTEGER, SAVE :: nini = 0
! check initialization

! nx=nx2/2
! ny=ny2/2
! nz=nz2/2

    tnorm = 1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2, DP))

    CALL copy1dto3d(q1, ffta, nx2, ny2, nz2)

    CALL fftw_execute_dft(pforw, ffta, ffta)

    CALL copy3dto1d(ffta, q2, tnorm, nx2, ny2, nz2)

    RETURN
  END SUBROUTINE fftf

! ******************************

  SUBROUTINE fftback(q1, q2)

! ******************************

    USE params
    USE FFTW
    IMPLICIT REAL(DP) (A - H, O - Z)

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax), &
! wrkx(kfft2),wrky(kfft2),wrkz(kfft2), &
! wsavex(kfft2),wsavey(kfft2),wsavez(kfft2), &
! ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

! nxyf=nx2*ny2
! nyf=nx2
! nx=nx2/2
! ny=ny2/2
! nz=nz2/2
!
!
    facnr = SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2, DP))

    CALL secopy1dto3d(q1, ffta, nx2, ny2, nz2)

    CALL fftw_execute_dft(pback, ffta, ffta)

    CALL secopy3dto1d(ffta, q2, facnr, nx2, ny2, nz2)

    RETURN
  END SUBROUTINE fftback

! ******************************

  SUBROUTINE rftf(q1, q2)

! ******************************

    USE params
    USE FFTW
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

! check initialization

! nx=nx2/2
! ny=ny2/2
! nz=nz2/2

    tnorm = 1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2, DP))

    CALL copyr1dto3d(q1, ffta, nx2, ny2, nz2)

    CALL fftw_execute_dft(pforw, ffta, ffta)

    CALL copy3dto1d(ffta, q2, tnorm, nx2, ny2, nz2)

    RETURN
  END SUBROUTINE rftf

! ******************************

  SUBROUTINE rfftback(q1, q3)
!SUBROUTINE rfftback(q1,q2)

! ******************************

    USE params
    USE FFTW
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

! nxyf=nx2*ny2
! nyf=nx2
! nx=nx2/2
! ny=ny2/2
! nz=nz2/2

    facnr = SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2, DP))

    CALL secopy1dto3d(q1, ffta, nx2, ny2, nz2)

    CALL fftw_execute_dft(pback, ffta, ffta)

    CALL copyr3dto1d(ffta, q3, facnr, nx2, ny2, nz2)

    RETURN
  END SUBROUTINE rfftback

! ******************************

  SUBROUTINE copy1dto3d(q1, ffta, nbx2, nby2, nbz2)

! ******************************

    USE params

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ffta(nbx2, nby2, nbz2)

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          ffta(MOD(i1 + nx, nbx2) + 1, MOD(i2 + ny, nby2) + 1, MOD(i3 + nz, nbz2) + 1) = q1(ind)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copy1dto3d

! ******************************

  SUBROUTINE copy3dto1d(ffta, q2, coef, nbx2, nby2, nbz2)

! ******************************

    USE params

    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: ffta(nbx2, nby2, nbz2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = coef*ffta(i1, i2, i3)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copy3dto1d

! ******************************

  SUBROUTINE copyr1dto3d(q1, ffta, nbx2, nby2, nbz2)

! ******************************

    USE params

    REAL(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ffta(nbx2, nby2, nbz2)

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          ffta(MOD(i1 + nx, nbx2) + 1, MOD(i2 + ny, nby2) + 1, MOD(i3 + nz, nbz2) + 1) = CMPLX(q1(ind), 0D0, DP)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copyr1dto3d

! ******************************

  SUBROUTINE copyr3dto1d(ffta, q2, coef, nbx2, nby2, nbz2)

! ******************************

    USE params

    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: ffta(nbx2, nby2, nbz2)
    REAL(DP), INTENT(OUT) :: q2(kdfull2)

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = REAL(ffta(MOD(i1 + nx, nbx2) + 1, MOD(i2 + ny, nby2) + 1, MOD(i3 + nz, nbz2) + 1), DP)*coef
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copyr3dto1d

! ******************************

  SUBROUTINE secopy1dto3d(q1, ffta, nbx2, nby2, nbz2)

! ******************************

    USE params

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: ffta(nbx2, nby2, nbz2)

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          ffta(i1, i2, i3) = q1(ind)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE secopy1dto3d

! ******************************

  SUBROUTINE secopy3dto1d(ffta, q2, coef, nbx2, nby2, nbz2)

! ******************************

    USE params

    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: ffta(nbx2, nby2, nbz2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          q2(ind) = ffta(MOD(i1 + nx, nbx2) + 1, MOD(i2 + ny, nby2) + 1, MOD(i3 + nz, nbz2) + 1)*coef
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE secopy3dto1d

! ******************************

  SUBROUTINE fft_end()

! ******************************

    USE FFTW

    CALL fftw_destroy_plan(pforw)
    CALL fftw_destroy_plan(pback)
    CALL fftw_destroy_plan(pforwx)
    CALL fftw_destroy_plan(pforwy)
    CALL fftw_destroy_plan(pforwz)
    CALL fftw_destroy_plan(pforwz1)
    CALL fftw_destroy_plan(pbackx)
    CALL fftw_destroy_plan(pbacky)
    CALL fftw_destroy_plan(pbackz)
    CALL fftw_destroy_plan(pbackz1)

    DEALLOCATE (ak, akv)
    DEALLOCATE (akpropx, akpropy, akpropz)
    DEALLOCATE (fftax, fftay, fftb, ffta)

    RETURN
  END SUBROUTINE fft_end

! ******************************

END MODULE kinetic
