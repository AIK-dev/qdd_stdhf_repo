
MODULE kinetic

! This module provides basic variables and arrays for FFT definition
! of kinetic energand Coulomb solver. Valid for FFTw3 as well as NETLIB.

  USE, INTRINSIC :: iso_c_binding
  USE params, ONLY: DP, numthr, PI
  IMPLICIT NONE

  SAVE
! Arrays for kinetic energy and electronic propagation.
! Kinetic energy coefficients in strange ordered Fourier space
! akv = Fourier-field for 0.5*k^2
! ak = Fourier-field for exp(i*dt*(h^2/2m)*k^2)

  INTEGER, PRIVATE, ALLOCATABLE :: modx(:), mody(:), modz(:)
  COMPLEX(DP), PARAMETER, PRIVATE :: eye = (0D0, 1D0)

  INTEGER, PRIVATE :: kfft, kfft2, kdfull2
  INTEGER, PRIVATE :: kxmax, kymax, kzmax
  INTEGER, PRIVATE :: iret
  INTEGER, PRIVATE :: ithr

  REAL(DP), ALLOCATABLE :: akv(:), akv3D(:, :, :)
  COMPLEX(DP), ALLOCATABLE :: akx(:), aky(:), akz(:) !DB
  COMPLEX(DP), ALLOCATABLE :: akpropx(:), akpropy(:), akpropz(:), akprop(:, :, :)
  COMPLEX(DP), ALLOCATABLE :: ak(:)


  INTEGER, PUBLIC :: FFTW_planflag
  COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: fftax(:,:), fftay(:,:), fftaz(:,:)
!  COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: fftb(:, :)
  COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: ffta(:, :, :, :)
  TYPE(C_PTR), PRIVATE, ALLOCATABLE :: pforwx(:), pforwy(:), pforwz(:)
  TYPE(C_PTR), PRIVATE, ALLOCATABLE :: pbackx(:), pbacky(:), pbackz(:)
!  TYPE(C_PTR), PRIVATE :: pforwz1, pbackz1
  TYPE(C_PTR), PRIVATE, ALLOCATABLE :: pforw(:), pback(:)
  INTEGER(C_INT), PRIVATE :: wisdomtest
  INTEGER, PRIVATE :: nacthr

CONTAINS
!-----init_grid_fft-----------------------------------------------------

  SUBROUTINE init_grid_fft(dx0, dy0, dz0, nx0, ny0, nz0, dt1, h2m)

! Initialize details for FFT, plans for FFTW3.
!
! Input:
! dx0,dy0,dz0 = grid spacing
! nx0,ny0,nz0 = nr. of grid points
! dt1 = electronic time step
! h2m = electronic hbar*+2/2m

    USE FFTW

    REAL(DP), INTENT(IN):: dx0, dy0, dz0
    INTEGER, INTENT(IN):: nx0, ny0, nz0
    REAL(DP), INTENT(IN):: dt1, h2m

    INTEGER:: nx, nx2, ny, ny2, nz, nz2
    INTEGER:: i1, i2, i3
    REAL(DP):: dkx, dky, dkz, facnr
    REAL(DP):: zkx, zky, zkz

    INTEGER, SAVE :: nxini = 0, nyini = 0, nzini = 0, nini = 0 ! flag for initialization
    INTEGER :: i

    INTEGER:: ind

#if(omp)
#if(!dynomp)
    INTEGER :: omp_get_num_threads
#endif
#endif
    nx2 = nx0; ny2 = ny0; nz2 = nz0
    kxmax = nx0; kymax = ny0; kzmax = nz0
    nx = nx2/2; ny = ny2/2; nz = nz2/2

    kfft = 2*kxmax
    kfft2 = kfft*2 + 1
    kdfull2 = kxmax*kymax*kzmax

    dkx = pi/(dx0*nx)
    dky = pi/(dy0*ny)
    dkz = pi/(dz0*nz)

    ALLOCATE (modx(kxmax), mody(kymax), modz(kzmax))
    ALLOCATE (ak(kdfull2), akv(kdfull2))
    ALLOCATE (akx(kdfull2), aky(kdfull2), akz(kdfull2)) !DB
    ALLOCATE (akpropx(kxmax), akpropy(kymax), akpropz(kzmax))
    ALLOCATE (akprop(kxmax, kymax, kzmax))
    ALLOCATE (akv3D(kxmax, kymax, kzmax))
!    ALLOCATE (fftb(kzmax, kxmax))
#if(omp)
    nacthr = numthr - 1
#else
    nacthr = 0
#endif
    WRITE (*, *) ' ALLOCATE with: kxmax,kymax,kzmax,nacthr=',&
        kxmax, kymax, kzmax, nacthr
    ALLOCATE (fftax(kxmax,0:nacthr), fftay(kymax,0:nacthr), fftaz(kzmax,0:nacthr))
    ALLOCATE (ffta(kxmax, kymax, kzmax, 0:nacthr), pforw(0:nacthr), pback(0:nacthr))
    ALLOCATE (pforwx(0:nacthr), pbackx(0:nacthr))
    ALLOCATE (pforwy(0:nacthr), pbacky(0:nacthr))
    ALLOCATE (pforwz(0:nacthr), pbackz(0:nacthr))
    WRITE (*, *) ' FFTA allocated with NACTHR=', nacthr
!
! central setting of FFTW planning expense --> edit here
!
! FFTW_planflag = FFTW_MEASURE
    FFTW_planflag = FFTW_PATIENT+FFTW_UNALIGNED
! FFTW_planflag = FFTW_EXHAUSTIVE

    WRITE (7, *) 'h bar squared over two m electron', h2m
    WRITE (7, *) ' testprint EYE=', eye
    WRITE (7, *) ' dkx,dky,dkz=', dkx, dky, dkz
    WRITE (7, *) ' testprint: nx2,ny2,nz2=', nx2, ny2, nz2
    WRITE (7, *) ' testprint: kdfull2,kfft2=', kdfull2, kfft2
    WRITE (6, *) ' testprint: nx2,ny2,nz2=', nx2, ny2, nz2
    WRITE (6,'(a,3i4)') 'testprint: kxmax,kymax,kzmax=',kxmax,kymax,kzmax


! prepare k**2 and kinetic propagation factor in 3D momentum space
    facnr = 1D0/(nx2*ny2*nz2)
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
          akx(ind) = zkx!*facnr
          aky(ind) = zky!*facnr
          akz(ind) = zkz!*facnr
          akprop(i1, i2, i3) = ak(ind)*facnr
          akv3D(i1, i2, i3) = akv(ind)*facnr
        END DO
      END DO
    END DO

! prepare kinetic propagation factors in 1D momentum spaces

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

! compose to kinetic propagation factors in 3D momentum spaces

    DO i3 = 1, nz2; DO i2 = 1, ny2; DO i1 = 1, nx2
          akprop(i1, i2, i3) = akpropx(i1)*akpropy(i2)*akpropz(i3)
        END DO; END DO; END DO
    akprop(:, :, :) = akprop(:, :, :)*facnr

! book-keeping arrays
    DO i1 = 1, nx2
      modx(i1) = MOD(i1 + nx, nx2) + 1
    END DO
    DO i2 = 1, ny2
      mody(i2) = MOD(i2 + ny, ny2) + 1
    END DO
    DO i3 = 1, nz2
      modz(i3) = MOD(i3 + nz, nz2) + 1
    END DO
    WRITE(*,'(a,200i4)') 'MODX:',modx


! FFTW3 switch, preparation of 'plans'

#if(omp)
    nacthr = numthr - 1
#else
    nacthr = 0
#endif
    !
    ! central setting of FFTW planning expense --> edit here
    !
    FFTW_planflag = FFTW_MEASURE
    ! FFTW_planflag = FFTW_PATIENT
    ! FFTW_planflag = FFTW_EXHAUSTIVE

    IF (wisdomtest == 0) THEN
      wisdomtest = fftw_import_system_wisdom()
      IF (wisdomtest == 0) THEN
        WRITE (6, *) 'wisdom_fftw.dat not found, creating it'
        WRITE (7, *) 'wisdom_fftw.dat not found, creating it'
      ELSE
        WRITE (*, *) 'wisdom from system'
      END IF
    END IF



    IF (nini == 0) THEN
#if(fftwnomkl)
      wisdomtest = fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    END IF

!   initialize plans, always in sequential mode
    WRITE(*,*) 'FFTINIT: nacthr=',nacthr
    DO i = 0, nacthr
      WRITE(*,*) 'FFTINIT: i=',i
      IF (nxini == 0) THEN
        pforwx(i) = fftw_plan_dft_1d(nx2, fftax(:,i), fftax(:,i), FFTW_FORWARD, FFTW_planflag)
        pbackx(i) = fftw_plan_dft_1d(nx2, fftax(:,i), fftax(:,i), FFTW_BACKWARD, FFTW_planflag)
        IF(i==nacthr) nxini = nx2
        fftax(:,i) = 1D0
        CALL fftw_execute_dft(pforwx(i), fftax(:,i), fftax(:,i))
      ELSE IF (nxini /= nx2) THEN
        STOP ' nx2 in four3d not as initialized!'
      END IF
      IF (nyini == 0) THEN
        pforwy(i) = fftw_plan_dft_1d(ny2, fftay(:,i), fftay(:,i), FFTW_FORWARD, FFTW_planflag)
        pbacky(i) = fftw_plan_dft_1d(ny2, fftay(:,i), fftay(:,i), FFTW_BACKWARD, FFTW_planflag)
        IF(i==nacthr) nyini = ny2
      ELSE IF (nyini /= ny2) THEN
        STOP ' ny2 in four3d not as initialized!'
      END IF
      IF (nzini == 0) THEN
!      pforwz = fftw_plan_dft_1d(nz2, fftb(1, nx2), fftb(1, nx2), FFTW_FORWARD, FFTW_planflag)
!      pbackz = fftw_plan_dft_1d(nz2, fftb(1, nx2), fftb(1, nx2), FFTW_BACKWARD, FFTW_planflag)
        pforwz(i) = fftw_plan_dft_1d(nz2, fftaz(:,i), fftaz(:,i), FFTW_FORWARD, FFTW_planflag)
        pbackz(i) = fftw_plan_dft_1d(nz2, fftaz(:,i), fftaz(:,i), FFTW_BACKWARD, FFTW_planflag)
        IF(i==nacthr) nzini = nz2
      ELSE IF (nzini /= nz2) THEN
        STOP ' nz2 in four3d not as initialized!'
      END IF
    END DO

#if(omp)
#if(!dynomp)
    CALL dfftw_init_threads(iret)
    WRITE (*, *) ' dfftw_init_threads: iret=', iret
    CALL dfftw_plan_with_nthreads(numthr)
    WRITE (*, *) ' init fft FFTW threads: nr. of threads=', numthr
#endif
#endif

! initialitze 3D plans, optionally for OpenMP
    IF (nini == 0) THEN
      DO i = 0, nacthr
        pforw(i) = fftw_plan_dft_3d(nz2, ny2, nx2, ffta(:, :, :, i), ffta(:, :, :, i), &
                                    FFTW_FORWARD, FFTW_planflag)
        pback(i) = fftw_plan_dft_3d(nz2, ny2, nx2, ffta(:, :, :, i), ffta(:, :, :, i), &
                                    FFTW_BACKWARD, FFTW_planflag)
      END DO
      nini = nx2*ny2*nz2
      WRITE (*, *) ' initialized nini=', nini, nx2, ny2, nz2
    ELSE IF (nini /= nx2*ny2*nz2) THEN
      WRITE (*, *) ' nini,nx2,ny2,nz2=', nini, nx2, ny2, nz2
      STOP ' nx2, ny2 or/and nz2 in four3d not as initialized!'
    END IF

#if(fftwnomkl)
    wisdomtest = fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    IF (wisdomtest == 0) THEN
      WRITE (6, *) 'Error exporting wisdom to FILE wisdom_fftw.dat'
      WRITE (7, *) 'Error exporting wisdom to FILE wisdom_fftw.dat'
    ELSE
      WRITE (*, *) ' successfull export of wisdom to wisdom_fftw.dat'
    END IF
    CALL fftw_forget_wisdom



! END of FFTW3 switch

    WRITE(*,*) 'END INIT: nx2,kxmax=',nx2,kxmax

  END SUBROUTINE init_grid_fft

! ******************************

  SUBROUTINE kinprop(q1)

! ******************************

! Propagation with kinetic energy exp(-i*dt*e_kin).
!
! Input/Output:
! q1 = coordinate-space array for propagated wavefunction

    USE params
    USE FFTW

    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN OUT) :: q1(kxbox, kybox, kzbox)
!COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)

    COMPLEX(DP), ALLOCATABLE :: ffttax(:), ffttay(:), ffttaz(:), ffttb(:, :)
!REAL(DP)::facnr
!REAL(DP):: tnorm, xfnorm, yfnorm, zfnorm

    LOGICAL, PARAMETER :: ttimestop = .FALSE.
    REAL(DP) :: acctime(0:2) = (/0D0, 0D0, 0D0/)
    INTEGER :: it1, it2, it3, it4, it5, it6

!tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

    ithr = 0
#if(omp)
#if(dynomp)
    ithr = OMP_GET_THREAD_NUM()
    IF (ithr > nacthr) THEN
      WRITE (*, *) ' in kinprop: ithr,nacthr=', ithr, nacthr
      STOP "too large ITHR"
    END IF
#endif
#endif
!facnr = 1D0/(nx2*ny2*nz2)

    IF (ttimestop) CALL system_clock(it1)
    IF (ttimestop) CALL system_clock(it2)
    CALL fftw_execute_dft(pforw(ithr), q1, q1)
    IF (ttimestop) CALL system_clock(it3)
    q1 = akprop*q1
    IF (ttimestop) CALL system_clock(it4)
    CALL fftw_execute_dft(pback(ithr), q1, q1)
    IF (ttimestop) CALL system_clock(it5)
    IF (ttimestop) THEN
      CALL system_clock(it6)
      acctime(0) = acctime(0) + 1D0
      acctime(1) = acctime(1) + (it4 - it3)*systime_factor
      acctime(2) = acctime(2) + (it3 - it2 + it5 - it4)*systime_factor
      IF (MOD(acctime(0), 8D0) == 0D0) WRITE (6, '(a,10(1pg11.3))') &
        'systime4: prop,fft=', acctime(1)/acctime(0), acctime(2)/acctime(0)
    END IF

    RETURN

  END SUBROUTINE kinprop
! ******************************

  SUBROUTINE kinetic_energy(q1, q2)

! ******************************

! ACTION of the kinetic energy
!
! Input/Output:
! q1 = coordinate-space array for propagated wavefunction

    USE params
    USE FFTW

    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)


    INTEGER :: ithr
    COMPLEX(DP), ALLOCATABLE :: qwork(:, :, :)


    ithr = 0
#if(omp)
#if(dynomp)
    ithr = OMP_GET_THREAD_NUM()
    IF (ithr > nacthr) THEN
      WRITE (*, *) ' in kinetic_energy: ithr,nacthr=', ithr, nacthr
      STOP "too large ITHR"
    END IF
#endif
#endif

    ALLOCATE (qwork(kxbox, kybox, kzbox))
    CALL copy1dto3d(q1, qwork, nx2, ny2, nz2)
    CALL fftw_execute_dft(pforw(ithr), qwork, qwork)
    qwork = akv3D*qwork
    CALL fftw_execute_dft(pback(ithr), qwork, qwork)
    CALL secopy3dto1d(qwork, q2, 1D0, nx2, ny2, nz2)
    DEALLOCATE (qwork)

    RETURN

  END SUBROUTINE kinetic_energy

  SUBROUTINE calc_ekin(psin, ekinout)

! Calculates kinetic energy for single particle state with
! complex wavefunction 'psin'.

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: psin(kdfull2)
    REAL(DP), INTENT(OUT) :: ekinout

    REAL(DP) :: sum0
    REAL(DP) :: sumk, sum0ex
    INTEGER ::ii
    REAL(DP) :: vol
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: psi2

!------------------------------------------------------------------

    ALLOCATE (psi2(kdfull2))

    CALL fftf(psin, psi2)
    sum0 = 0D0
    sumk = 0D0
    DO ii = 1, kdfull2
      vol = REAL(psi2(ii), DP)*REAL(psi2(ii), DP) + AIMAG(psi2(ii))*AIMAG(psi2(ii))
      sum0 = vol + sum0
      sumk = vol*akv(ii) + sumk
    END DO
    sum0ex = 1D0/((2D0*PI)**3*dx*dy*dz)
    ekinout = sumk/sum0ex

    DEALLOCATE (psi2)

    RETURN
  END SUBROUTINE calc_ekin

! ******************************

  SUBROUTINE gradient(fin, gradfout, idirec)

! The gradient of complex field 'fin' in direction 'idirec'
! (x =1, y=2, z=3).
! The fields are given in Fourier space and the gradient is
! applied as product with 'kx', 'ky', 'kz'.
!

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: fin(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: gradfout(kdfull2)
    INTEGER, INTENT(IN) :: idirec

! ************************************************************

    IF (idirec == 1) THEN
! x-derivative
      gradfout = -akx*fin
    ELSEIF (idirec == 2) THEN
! y-derivative
      gradfout = -aky*fin
    ELSEIF (idirec == 3) THEN
! z-derivative
      gradfout = -akz*fin
    ELSE
      STOP ' RGRADIENT called with invalid IDIREC'
    END IF

    RETURN
  END SUBROUTINE gradient

! BLOCK of subroutines for FFTW3
! ******************************

  SUBROUTINE xgradient_rspace(fin, gradfout)

! The gradient of the complex field 'fin' in x-direction.
! The fields are given in coordinate space. The gradient is
! evaluated in k_x space.
!

    USE params !, ONLY: kdfull2,pi,dx,nxyf,nyf,nx2,ny2,nz2,nx,omp_get_thread_num
    USE FFTW
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: fin(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: gradfout(kdfull2)

    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: work(:)

    INTEGER:: i1, i2, i3, ind
    INTEGER :: ithra
    REAL(DP):: dkx, zkx
! ************************************************************

!    ALLOCATE(work(kxbox))
    ithra = 0
#if(omp)
#if(dynomp)
    ithra = OMP_GET_THREAD_NUM()
    IF (ithra > nacthr) THEN
      WRITE (*, *) ' in xgradient: ithra,nacthr=', ithra, nacthr
      STOP "too large ITHRA"
    END IF
!    WRITE(*,*) 'XGRADIENT:',ithra,nacthr,nthr
!    WRITE(*,'(a,200i4)') 'MODX:',modx
!    CALL flush(6)
#endif
#endif

    dkx = pi/(dx*nx)
#if(omp)
#if(!dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ithra,i3,i2,i1,ind,zkx) SCHEDULE(STATIC)
#endif
#endif
    DO i3 = 1, nz2
#if(omp)
#if(!dynomp)
      ithra = OMP_GET_THREAD_NUM()
#endif
#endif
      DO i2 = 1, ny2
! forward transform along x
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftax(modx(i1),ithra) = fin(ind) ! copy to workspace
!          work(modx(i1)) = fin(ind) ! copy to workspace
        END DO
        CALL fftw_execute_dft(pforwx(ithra), fftax(:,ithra), fftax(:,ithra))
!        CALL fftw_execute_dft(pforwx(ithra), work, work)
! multiply by k_x factor in k_x space
        DO i1 = 1, nx2
          IF (i1 == (nx + 1)) THEN
            zkx = 0D0
          ELSE IF (i1 > (nx + 1)) THEN
            zkx = (i1 - nx2 - 1)*dkx
          ELSE
            zkx = (i1 - 1)*dkx
          END IF
          fftax(i1,ithra) = fftax(i1,ithra)*eye*zkx
!          work(i1) = work(i1)*eye*zkx
        END DO
! backward transform along x
        CALL fftw_execute_dft(pbackx(ithra), fftax(:,ithra), fftax(:,ithra))
!        CALL fftw_execute_dft(pbackx(ithra), work, work)
        DO i1 = 1, nx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          gradfout(ind) = fftax(modx(i1),ithra)/nx2
!          gradfout(ind) = work(modx(i1))/nx2
        END DO
      END DO
    END DO
#if(omp)
#if(!dynomp)
!$OMP END PARALLEL DO
#endif
#endif
!    DEALLOCATE(work)

    RETURN

  END SUBROUTINE xgradient_rspace

  SUBROUTINE ygradient_rspace(fin, gradfout)

! The gradient of the complex field 'fin' in y-direction.
! The fields are given in coordinate space. The gradient is
! evaluated in k_y space.
!

    USE params
    USE FFTW
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: fin(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: gradfout(kdfull2)

    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: work(:)

    INTEGER:: i1, i2, i3, ind
    INTEGER :: ithra
    REAL(DP):: dky, zky
! ************************************************************

!    ALLOCATE(work(kybox))
    ithra = 0
#if(omp)
#if(dynomp)
    ithra = OMP_GET_THREAD_NUM()
    IF (ithra > nacthr) THEN
      WRITE (*, *) ' in ygradient: ithra,nacthr=', ithra, nacthr
      STOP "too large ITHRA"
    END IF
#endif
#endif

    dky = pi/(dy*ny)
#if(omp)
#if(!dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ithra,i3,i2,i1,ind,zky) SCHEDULE(STATIC)
#endif
#endif
    DO i3 = 1, nz2
#if(omp)
#if(!dynomp)
      ithra = OMP_GET_THREAD_NUM()
#endif
#endif
      DO i1 = 1, nx2
! forward transform along y
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          fftay(mody(i2),ithra) = fin(ind)
!          work(mody(i2)) = fin(ind)
        END DO
        CALL fftw_execute_dft(pforwy(ithra), fftay(:,ithra), fftay(:,ithra))
!        CALL fftw_execute_dft(pforwy(ithra), work(:), work(:))
! multiply by k_y factor in k_y space
        DO i2 = 1, ny2
          IF (i2 == (ny + 1)) THEN
            zky = 0D0
          ELSE IF (i2 > (ny + 1)) THEN
            zky = (i2 - ny2 - 1)*dky
          ELSE
            zky = (i2 - 1)*dky
          END IF
          fftay(i2,ithra) = fftay(i2,ithra)*eye*zky
!          work(i2) = work(i2)*eye*zky
        END DO
! backward transform along y
        CALL fftw_execute_dft(pbacky(ithra), fftay(:,ithra), fftay(:,ithra))
!        CALL fftw_execute_dft(pbacky(0), work(:), work(:))
        DO i2 = 1, ny2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          gradfout(ind) = fftay(mody(i2),ithra)/ny2
!          gradfout(ind) = work(mody(i2))/ny2
        END DO
      END DO
    END DO
#if(omp)
#if(!dynomp)
!$OMP END PARALLEL DO
#endif
#endif
!    DEALLOCATE(work)

    RETURN

  END SUBROUTINE ygradient_rspace

  SUBROUTINE zgradient_rspace(fin, gradfout)

! The gradient of the complex field 'fin' in z-direction.
! The fields are given in coordinate space. The gradient is
! evaluated in k_z space.
!

    USE params
    USE FFTW
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: fin(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: gradfout(kdfull2)

    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: work(:)

    INTEGER:: i1, i2, i3, i3m, ind
    INTEGER :: ithra
    REAL(DP):: dkz, zkz

! ************************************************************

!    ALLOCATE(work(kzbox))
    ithra = 0
#if(omp)
#if(dynomp)
    ithra = OMP_GET_THREAD_NUM()
    IF (ithra > nacthr) THEN
      WRITE (*, *) ' in zgradient: ithra,nacthr=', ithra, nacthr
      STOP "too large ITHRA"
    END IF
#endif
#endif

    dkz = pi/(dz*nz)
#if(omp)
#if(!dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ithra,i3,i2,i1,ind,zkz) SCHEDULE(STATIC)
#endif
#endif
    DO i2 = 1, ny2
#if(omp)
#if(!dynomp)
      ithra = OMP_GET_THREAD_NUM()
#endif
#endif
      DO i1 = 1, nx2
! forward transform along z
        DO i3 = 1, nz2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          i3m = MOD(i3 + nz, nz2) + 1
! i3m = modz(i3)
          fftaz(i3m,ithra) = fin(ind)
!          work(i3m) = fin(ind)
        END DO
        CALL fftw_execute_dft(pforwz(ithra), fftaz(:,ithra), fftaz(:,ithra))
!        CALL fftw_execute_dft(pforwz(ithra), work(:), work(:))
! multiply by k_z factor in k_z space
        DO i3 = 1, nz2
          IF (i3 == (nz + 1)) THEN
            zkz = 0D0
          ELSE IF (i3 > (nz + 1)) THEN
            zkz = (i3 - nz2 - 1)*dkz
          ELSE
            zkz = (i3 - 1)*dkz
          END IF
          fftaz(i3,ithra) = fftaz(i3,ithra)*eye*zkz
!          work(i3) = work(i3)*eye*zkz
        END DO
! backward transform along z
        CALL fftw_execute_dft(pbackz(ithra), fftaz(:,ithra), fftaz(:,ithra))
!        CALL fftw_execute_dft(pbackz(ithra), work(:), work(:))
        DO i3 = 1, nz2 ! copy back
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          i3m = MOD(i3 + nz, nz2) + 1
          gradfout(ind) = fftaz(modz(i3),ithra)/nz2
!          gradfout(ind) = work(modz(i3))/nz2
        END DO
!
      END DO
    END DO
#if(omp)
#if(!dynomp)
!$OMP END PARALLEL DO
#endif
#endif
!    DEALLOCATE(work)

    RETURN

  END SUBROUTINE zgradient_rspace

  SUBROUTINE fftf(q1, q2)

! Forward 3D FFT from 'q1' in r-space to 'q2' in k-space.

    USE params
    USE FFTW
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)

    REAL(DP):: tnorm
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: work(:, :, :)
    INTEGER :: ithra

!    ALLOCATE(work(nx2,ny2,nz2))

    tnorm = 1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2, DP))

    ithra = 0
#if(omp)
#if(dynomp)
    ithra = OMP_GET_THREAD_NUM()
    IF (ithra > nacthr) THEN
      WRITE (*, *) ' in FFTF: ithra,nacthr=', ithra, nacthr
      STOP "too large ITHRA"
    END IF
#endif
#endif
    CALL copy1dto3d(q1, ffta(:, :, :, ithra), nx2, ny2, nz2)
    CALL fftw_execute_dft(pforw(ithra), ffta(:,:,:,ithra), ffta(:,:,:,ithra))
    CALL copy3dto1d(ffta(:, :, :, ithra), q2, tnorm, nx2, ny2, nz2)
!    CALL copy1dto3d(q1, work(:, :, :), nx2, ny2, nz2)
!    CALL fftw_execute_dft(pforw(ithra), work(:,:,:), work(:,:,:))
!    CALL copy3dto1d(work(:, :, :), q2, tnorm, nx2, ny2, nz2)
!    DEALLOCATE(work)

    RETURN
  END SUBROUTINE fftf

  SUBROUTINE fftback(q1, q2)

! Backward 3D FFT from 'q1' in k-space to 'q2' in r-space.

    USE params
    USE FFTW
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)
    REAL(DP):: facnr
    COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: work(:, :, :)
    INTEGER :: ithra

!    ALLOCATE(work(nx2,ny2,nz2))
    facnr = SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2, DP))

    ithra = 0
#if(omp)
#if(dynomp)
    ithra = OMP_GET_THREAD_NUM()
    IF (ithra > nacthr) THEN
      WRITE (*, *) ' in FFTF: ithra,nacthr=', ithra, nacthr
      STOP "too large ITHRA"
    END IF
#endif
#endif
    CALL secopy1dto3d(q1, ffta(:, :, :, ithra), nx2, ny2, nz2)
    CALL fftw_execute_dft(pback(ithra), ffta(:, :, :, ithra), ffta(:, :, :, ithra))
    CALL secopy3dto1d(ffta(:, :, :, ithra), q2, facnr, nx2, ny2, nz2)
!    CALL secopy1dto3d(q1, work(:, :, :), nx2, ny2, nz2)
!    CALL fftw_execute_dft(pback(ithra), work(:, :, :), work(:, :, :))
!    CALL secopy3dto1d(work(:, :, :), q2, facnr, nx2, ny2, nz2)

!    DEALLOCATE(work)

    RETURN
  END SUBROUTINE fftback

  SUBROUTINE rftf(q1, q2)

! Forward 3D FFT from real r-space array 'q1' to complex k-space 'q2'.

    USE params
    USE FFTW
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2)

    REAL(DP):: tnorm

    tnorm = 1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2, DP))

    CALL copyr1dto3d(q1, ffta(:, :, :, 0), nx2, ny2, nz2)
    CALL fftw_execute_dft(pforw(0), ffta(:, :, :, 0), ffta(:, :, :, 0))
    CALL copy3dto1d(ffta(:, :, :, 0), q2, tnorm, nx2, ny2, nz2)

    RETURN
  END SUBROUTINE rftf

  SUBROUTINE rfftback(q1, q3)

! Backward 3D FFT from complex k-space 'q1' to real r-space array 'q3'.

    USE params
    USE FFTW
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
    REAL(DP), INTENT(OUT) :: q3(kdfull2)

    REAL(DP):: facnr


    facnr = SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2, DP))

    CALL secopy1dto3d(q1, ffta(:, :, :, 0), nx2, ny2, nz2)
    CALL fftw_execute_dft(pback(0), ffta(:, :, :, 0), ffta(:, :, :, 0))
    CALL copyr3dto1d(ffta(:, :, :, 0), q3, facnr, nx2, ny2, nz2)

    RETURN
  END SUBROUTINE rfftback

  SUBROUTINE copy1dto3d(vec1d, vec3d, nbx2, nby2, nbz2)

! Copies 3D array from linear storage to 3D storage (both COMPLEX).
!
! Input:
! vec1d = array in linear storage
! nbx2,nby2,nbz2 = dimensions of 3D array
! Output:
! vec3d = array in 3D storage

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbx2
    INTEGER, INTENT(IN) :: nby2
    INTEGER, INTENT(IN) :: nbz2
    COMPLEX(DP), INTENT(IN) :: vec1d(kdfull2)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: vec3d(nbx2, nby2, nbz2)

    INTEGER:: i1, i2, i3, ind

    ind = 0
    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = ind + 1
          vec3d(i1, i2, i3) = vec1d(ind)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copy1dto3d

! ******************************

  SUBROUTINE copyr1dto3d(vec1d, vec3d, nbx2, nby2, nbz2)

! Copies 3D array from real & linear storage to complex & 3D storage.
!
! Input:
! vec1d = array in linear storage
! nbx2,nby2,nbz2 = dimensions of 3D array
! Output:
! vec3d = array in 3D storage

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbx2
    INTEGER, INTENT(IN) :: nby2
    INTEGER, INTENT(IN) :: nbz2
    REAL(DP), INTENT(IN) :: vec1d(kdfull2)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: vec3d(nbx2, nby2, nbz2)
    INTEGER:: i1, i2, i3, ind

    ind = 0
    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = ind + 1
          vec3d(i1, i2, i3) = CMPLX(vec1d(ind), 0D0, DP)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copyr1dto3d

! ******************************

  SUBROUTINE secopy1dto3d(vec1d, vec3d, nbx2, nby2, nbz2)

! Copies 3D array from linear storage to 3D storage (both COMPLEX).
!
! Input:
! vec1d = array in linear storage
! nbx2,nby2,nbz2 = dimensions of 3D array
! Output:
! vec3d = array in 3D storage
!
! Does the same as 'copy1dto3d' but with faster (?) index computation.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbx2
    INTEGER, INTENT(IN) :: nby2
    INTEGER, INTENT(IN) :: nbz2
    COMPLEX(DP), INTENT(IN) :: vec1d(kdfull2)
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: vec3d(nbx2, nby2, nbz2)
    INTEGER:: i1, i2, i3, ind

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          vec3d(i1, i2, i3) = vec1d(ind)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE secopy1dto3d

! ******************************

  SUBROUTINE copy3dto1d(vec3d, vec1d, coef, nbx2, nby2, nbz2)

! Copies 3D array from 3D storage to linear storage (both COMPLEX).
!
! Input:
! vec3d = array in 3D storage
! nbx2,nby2,nbz2 = dimensions of 3D array
! Output:
! vec1d = array in linear storage
!

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) ::coef
    INTEGER, INTENT(IN) :: nbx2
    INTEGER, INTENT(IN) :: nby2
    INTEGER, INTENT(IN) :: nbz2
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: vec3d(nbx2, nby2, nbz2)
    COMPLEX(DP), INTENT(OUT) :: vec1d(kdfull2)
    INTEGER:: i1, i2, i3, ind

    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = (i3 - 1)*nxyf + (i2 - 1)*nyf + i1
          vec1d(ind) = coef*vec3d(i1, i2, i3)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copy3dto1d

! ******************************

  SUBROUTINE copyr3dto1d(vec3d, vec1d, coef, nbx2, nby2, nbz2)

! Copies 3D array from 3D storage to linear storage.
!
! Input:
! vec3d = complex array in 3D storage
! nbx2,nby2,nbz2 = dimensions of 3D array
! Output:
! vec1d = real array in linear storage
!

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) ::coef
    INTEGER, INTENT(IN) :: nbx2
    INTEGER, INTENT(IN) :: nby2
    INTEGER, INTENT(IN) :: nbz2
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: vec3d(nbx2, nby2, nbz2)
    REAL(DP), INTENT(OUT) :: vec1d(kdfull2)

    INTEGER:: i1, i2, i3, ind

    ind = 0
    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = ind + 1
          vec1d(ind) = REAL(vec3d(i1, i2, i3), DP)*coef
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE copyr3dto1d

! ******************************

  SUBROUTINE secopy3dto1d(vec3d, vec1d, coef, nbx2, nby2, nbz2)

! Copies 3D array from 3D storage to linear storage (both COMPLEX).
!
! Input:
! vec3d = array in 3D storage
! nbx2,nby2,nbz2 = dimensions of 3D array
! Output:
! vec1d = array in linear storage
!
! Does the same as 'copy1dto3d' but with faster (?) index computation.

    USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) ::coef
    INTEGER, INTENT(IN) :: nbx2
    INTEGER, INTENT(IN) :: nby2
    INTEGER, INTENT(IN) :: nbz2
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: vec3d(nbx2, nby2, nbz2)
    COMPLEX(DP), INTENT(OUT) :: vec1d(kdfull2)

    INTEGER:: i1, i2, i3, ind

    ind = 0
    DO i3 = 1, nbz2
      DO i2 = 1, nby2
        DO i1 = 1, nbx2
          ind = ind + 1
          vec1d(ind) = vec3d(i1, i2, i3)*coef
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE secopy3dto1d

! ******************************

  SUBROUTINE fft_end()

! FFTW epilogue.

    USE FFTW
    IMPLICIT NONE

    CALL fftw_destroy_plan(pforw(0))
    CALL fftw_destroy_plan(pback(0))
    CALL fftw_destroy_plan(pforwx(0))
    CALL fftw_destroy_plan(pforwy(0))
    CALL fftw_destroy_plan(pforwz(0))
!    CALL fftw_destroy_plan(pforwz1)
    CALL fftw_destroy_plan(pbackx(0))
    CALL fftw_destroy_plan(pbacky(0))
    CALL fftw_destroy_plan(pbackz(0))
!    CALL fftw_destroy_plan(pbackz1)

    DEALLOCATE (ak, akv)
    DEALLOCATE (akpropx, akpropy, akpropz)
    DEALLOCATE (fftax, fftay, fftaz, ffta)
!    DEALLOCATE (fftb)

    RETURN
  END SUBROUTINE fft_end

! END BLOCK of subroutines for FFTW3

! ******************************

END MODULE kinetic
