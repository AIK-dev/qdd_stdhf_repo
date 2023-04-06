MODULE cukinetic

USE params, ONLY: DP,numthr,PI,kdfull2
USE kinetic, ONLY: akv
USE cudafor

IMPLICIT NONE

SAVE

INTEGER, PRIVATE :: kxmax, kymax, kzmax

integer,  public :: CUFFT_FORWARD = -1
integer,  public :: CUFFT_INVERSE =  1
!integer, public :: CUFFT_C2C = Z'29'
integer,  public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex
COMPLEX(DP),device,ALLOCATABLE :: akprop_d(:,:,:)
REAL(DP) , device , ALLOCATABLE :: akv_d(:)
type(dim3) :: gridfft3d,tBlockfft3d !!! (real(8),real(8),real(8))
integer, parameter, public :: fp_kind = kind(0.0d0)
complex(DP), device , allocatable :: ffta_d(:,:,:)
complex(DP), device , allocatable :: fullfftx_d(:),fullfftay_d(:),fullfftaz_d(:)
integer, allocatable,save :: fftsize(:),inembed(:),onembed(:)
type(c_ptr) ,save :: planmany1dx,planmany1dy,planmany1dz
integer, allocatable,save :: sizefftx(:),xinembed(:),xonembed(:)
integer, allocatable,save :: sizeffty(:),yinembed(:),yonembed(:)
integer, allocatable,save :: sizefftz(:),zinembed(:),zonembed(:)
type(c_ptr) ,save :: planmany,plan1dx,plan3d
  interface cufftDestroy
     subroutine cufftDestroy(plan) bind(C,name='cufftDestroy') 
       use iso_c_binding
       type(c_ptr),value:: plan
     end subroutine cufftDestroy
  end interface cufftDestroy

  interface cufftExec
     subroutine cufftExecZ2Z(plan, idata, odata, direction) &
          bind(C,name='cufftExecZ2Z') 
       use iso_c_binding
       type(c_ptr),value:: plan
       integer(c_int),value:: direction
       !pgi$ ignore_tr idata,odata
       complex(8),device:: idata(*),odata(*)
     end subroutine cufftExecZ2Z
  end interface cufftExec
  interface cufftPlanMany
     subroutine cufftPlanMany(plan, rank, n, inembed, istride, idist, & 
          onembed, ostride, odist,  &
          type, batch) bind(C,name='cufftPlanMany')
       use iso_c_binding
       implicit none
       !pgi$ ignore_tkr n, inembed, onembed       
       type(c_ptr) :: plan
       integer(c_int) :: n, inembed, onembed
       integer(c_int), value:: rank, istride, ostride, idist, odist, type, batch
     end subroutine cufftPlanMany
 end interface cufftPlanMany
  interface cufftPlan1d
     subroutine cufftPlan1d(plan, nx, type, batch) &
          bind(C,name='cufftPlan1d') 
       use iso_c_binding
       type(c_ptr):: plan
       integer(c_int),value:: nx, batch,type
     end subroutine cufftPlan1d
  end interface cufftPlan1d
  interface cufftPlan3d
     subroutine cufftPlan3d(plan, nx, ny, nz, type, batch) &
          bind(C,name='cufftPlan3d') 
       use iso_c_binding
       type(c_ptr):: plan
       integer(c_int),value:: nx, ny, nz, batch,type
     end subroutine cufftPlan3d
  end interface cufftPlan3d

  CONTAINS

  SUBROUTINE init_cuda_fft(dx0, dy0, dz0, nx0, ny0, nz0, dt1, h2m)
use params ,only:nstate
USE kinetic, ONLY: akprop,akv
implicit none

REAL(DP), INTENT(IN):: dx0, dy0, dz0
INTEGER, INTENT(IN):: nx0, ny0, nz0
REAL(DP), INTENT(IN):: dt1, h2m


kxmax = nx0; kymax = ny0; kzmax = nz0
tBlockfft3d = dim3(16,16,4)
gridfft3d   = dim3(ceiling(real(kxmax)/tBlockfft3d%x),&
		   ceiling(real(kymax)/tBlockfft3d%y),&
		   ceiling(real(kzmax)/tBlockfft3d%z))
allocate(akv_d(kdfull2))
allocate(ffta_d(kxmax,kymax,kzmax))
allocate(fullfftx_d(kxmax*kymax*kzmax),fullfftay_d(kxmax*kymax*kzmax),fullfftaz_d(kxmax*kymax*kzmax))
allocate(fftsize(3),inembed(3),onembed(3))
fftsize(1) = kxmax
fftsize(2) = kymax
fftsize(3) = kzmax
inembed(:) = kxmax*kymax*kzmax*nstate
onembed(:) = kxmax*kymax*kzmax*nstate
write(*,*) "init_cuda_fft"
write(*,*) "kxmax,kymax,kzmax,nstate",kxmax,kymax,kzmax,nstate
 call cufftPlanMany(planmany, 3, fftsize, inembed, &
          1,kxmax*kymax*kzmax, &
          onembed, &
          1,kxmax*kymax*kzmax,  &
          CUFFT_Z2Z, nstate)
 call cufftPlan1d(plan1dx,kxmax,CUFFT_Z2Z,1)
 call cufftPlan3d(plan3d,kzmax,kymax,kxmax,CUFFT_Z2Z,1)
write(*,*) "allocate akprop_d"
 ALLOCATE(akprop_d(kxmax,kymax,kzmax))
akprop_d = akprop
akv_d    = akv
!---- planmany1dx
allocate(sizefftx(1),xinembed(1),xonembed(1))
sizefftx(:) = kxmax
xinembed(:) = kxmax*kymax*kzmax
xonembed(:) = kxmax*kymax*kzmax
call cufftPlanMany(planmany1dx,1,sizefftx,xinembed, &
        1,kxmax,    &
        xonembed,   &
        1,kxmax,    &
        CUFFT_Z2Z,kymax*kzmax)
!---- planmany1dy
allocate(sizeffty(1),yinembed(1),yonembed(1))
sizeffty(:) = kymax
yinembed(:) = kxmax*kymax*kzmax
yonembed(:) = kxmax*kymax*kzmax
call cufftPlanMany(planmany1dy,1,sizeffty,yinembed, &
        1,kymax,    &
        yonembed,   &
        1,kymax,    &
        CUFFT_Z2Z,kxmax*kzmax)
!---- planmany1dz
allocate(sizefftz(1),zinembed(1),zonembed(1))
sizefftz(:) = kzmax
zinembed(:) = kxmax*kymax*kzmax
zonembed(:) = kxmax*kymax*kzmax
call cufftPlanMany(planmany1dz,1,sizefftz,zinembed, &
        1,kzmax,    &
        zonembed,   &
        1,kzmax,    &
        CUFFT_Z2Z,kxmax*kymax)
END SUBROUTINE init_cuda_fft

END MODULE cukinetic