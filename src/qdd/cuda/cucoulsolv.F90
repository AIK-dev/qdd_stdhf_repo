!This file is a part of PW-TELEMAN project.
!PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
!Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
!Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.
!
!PW-Teleman is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!PW-Teleman is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.

MODULE cucoulsolv
    USE cudafor
    USE params, ONLY: DP,PI,numthr,e2,zero
    IMPLICIT NONE
    
    SAVE
    INTEGER,PRIVATE :: kxmax,kymax,kzmax, ksmax
    ! kxmax must be the largest
    INTEGER,PRIVATE :: kdfull
    INTEGER,PRIVATE :: kdred
    INTEGER,PRIVATE :: kfft,kfftx,kffty,kfftz
    !INTEGER,PARAMETER,PRIVATE :: kdcorf=(kxmax/2+1)*(kymax/2+1)*(kzmax/2+1)
    ! include block: xkgrid
    REAL(DP),PRIVATE :: dx,dy,dz,dxsp,grnorm,fnorm
    INTEGER,PRIVATE :: nx,ny,nz,nx1,ny1,nz1,nxi,nyi,nzi,nxy1,nxyz

    REAL(DP),ALLOCATABLE,PRIVATE :: akv2r(:),akv2i(:)
    REAL(DP),PRIVATE :: dkx,dky,dkz,akmax,dksp,ecut
    INTEGER,PRIVATE :: nxk,nxklo,nxkhi,nksp,nkxyz
  
    
    REAL(DP),DEVICE ,ALLOCATABLE :: rhokr_d(:),rhoki_d(:)
    REAL(DP), ALLOCATABLE :: rhointest(:)
    REAL(DP), ALLOCATABLE :: testrhokr(:), testrhoki(:)
    type(dim3) :: gridrhofld,tBlockrhofld,gridfftz,tBlockfftz
    type(dim3) :: gridresult,tBlockresult
    COMPLEX(DP),DEVICE,ALLOCATABLE :: fftax_d(:),fftay_d(:),fftb_d(:),fullfftax_d(:),fftayback_d(:)
    REAL(DP),DEVICE,ALLOCATABLE    :: akv2r_d(:),akv2i_d(:)
    type(c_ptr) ,save :: plan1dcoulsolvx,plan1dcoulsolvy,plan1dcoulsolvz,plan1dcoulsolvyback
    integer, allocatable,save :: fftsizex(:),inembedx(:),onembedx(:)
    integer, allocatable,save :: fftsizey(:),inembedy(:),onembedy(:)
    integer, allocatable,save :: fftsizeby(:),inembedby(:),onembedby(:)
    integer, allocatable,save :: fftsizez(:),inembedz(:),onembedz(:)
    
    CONTAINS
    
    SUBROUTINE cuinit_coul(dx0,dy0,dz0,nx0,ny0,nz0)
    USE CUKINETIC
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::dx0,dy0,dz0
    INTEGER,INTENT(IN)::nx0,ny0,nz0
    !-----------------------------------------------------------------------
    INTEGER :: ii,i1,i2,i3, kdum
    LOGICAL,PARAMETER :: tcoultest=.false.
    REAL(DP) :: charge
    REAL(DP),ALLOCATABLE :: rhotest(:),ctest(:)
    
    
    !     read grid parameters from file or simply initialize them
    !     note that the Coulomb solver doubles the grid internally
    
    kxmax=2*nx0;kymax=2*ny0;kzmax=2*nz0;ksmax=kxmax
    kdfull=nx0*ny0*nz0
    kdred=kxmax*kymax*kzmax
    kfft=ksmax;kfftx=kxmax;kffty=kymax;kfftz=kzmax
    
    nx=nx0  !/2
    ny=ny0  !/2
    nz=nz0  !/2
    dx=dx0
    dy=dy0
    dz=dz0
    
    
    !--------------------------------
    

    ALLOCATE(akv2r(kdred),akv2i(kdred))

    !     call input routine fftinp, which initializes the grid and fft table
    
    CALL cufftinp
    
    allocate(testrhokr(kdred),testrhoki(kdred))
    tBlockrhofld = dim3(16,16,4)
    gridrhofld   = dim3(ceiling(real(nx)/tBlockrhofld%x),ceiling(real(ny)/tBlockrhofld%y), ceiling(real(nz)/tBlockrhofld%z))
    tBlockresult = dim3(16,16,4)
    gridresult   = dim3(ceiling(real(nxi)/tBlockresult%x),ceiling(real(nyi)/tBlockresult%y), ceiling(real(nzi)/tBlockresult%z))
    !call cufftPlan1d(plan1dcoulsolvx, kfftx, CUFFT_Z2Z, 1)
    allocate(fftsizex(1),inembedx(1),onembedx(1))
    fftsizex(:) = kfftx
    inembedx(:) = kfftx*nzi*nyi
    onembedx(:) = kfftx*nzi*nyi
    call cufftPlanMany(plan1dcoulsolvx,1,fftsizex,inembedx, &
            1,kfftx,    &
            onembedx,   &
            1,kfftx,    &
            CUFFT_Z2Z,nzi*nyi)
    !call cufftPlan1d(plan1dcoulsolvy, kffty, CUFFT_Z2Z, 1)
    write(*,*) kffty,nzi,nx1
    allocate(fftsizey(1),inembedy(1),onembedy(1))
    fftsizey(:) = kffty
    inembedy(:) = kffty*nzi*nx1
    onembedy(:) = kffty*nzi*nx1
    call cufftPlanMany(plan1dcoulsolvy,1,fftsizey,inembedy, &
            1,kffty,    &
            onembedy,   &
            1,kffty,    &
            CUFFT_Z2Z,nzi*nx1)
    !backward
    allocate(fftsizeby(1),inembedby(1),onembedby(1))
    fftsizeby(:) = kffty
    inembedby(:) = kffty*nz1*nx1
    onembedby(:) = kffty*nz1*nx1
    call cufftPlanMany(plan1dcoulsolvyback,1,fftsizeby,inembedby, &
            1,kffty,    &
            onembedby,   &
            1,kffty,    &
            CUFFT_Z2Z,nz1*nx1)
    !call cufftPlan1d(plan1dcoulsolvz, kfftz, CUFFT_Z2Z, 1)
    allocate(fftsizez(1),inembedz(1),onembedz(1))
    fftsizez(:) = kfftz
    inembedz(:) = kfftz*nyi*nx1
    onembedz(:) = kfftz*nyi*nx1
    call cufftPlanMany(plan1dcoulsolvz,1,fftsizez,inembedz, &
            1,kfftz,    &
            onembedz,   &
            1,kfftz,    &
            CUFFT_Z2Z,nyi*nx1)
    
    ALLOCATE(fftax_d(kxmax),fftay_d(nzi*nx1*nyi),fftb_d(kfftz*nyi*nx1))
    ALLOCATE(fftayback_d(nz1*nx1*nyi))
    ALLOCATE(fullfftax_d(kdred))

    RETURN
    END SUBROUTINE cuinit_coul
    

    !-----fftinp------------------------------------------------------------
    
    SUBROUTINE cufftinp
    IMPLICIT NONE
    
    !     initializes work tables for FFT
    
    !     grid parameters nx,ny,nz,dx,dy,dz,ecut must have been read or
    !     initialized before !
    
    !-----------------------------------------------------------------------
    
    INTEGER ::  ii, i1, i2, i3, ind, ikzero
    REAL(DP) :: ak2, xx1, xx2, xy1, xy2, xz1, xz2

    ALLOCATE(rhokr_d(kdred),rhoki_d(kdred))
    ALLOCATE(akv2r_d(kdred),akv2i_d(kdred))

    !     initialize grid in coordinate space
    
    nx1=nx+1
    ny1=ny+1
    nz1=nz+1
    nxi=nx+nx
    nyi=ny+ny
    nzi=nz+nz
    nxy1=nxi*nyi
    nxyz=nxi*nyi*nzi
    nkxyz=nxi*nyi*nzi
    
    !     grid lengths must match with parameters in incs
    
    IF(kxmax < nxi) THEN
      WRITE(6,'(a)') ' ERROR: parameter   kxmax   too small'
      STOP ' error in parameter: KXMAX in COULEX too small'
    ELSE IF(kymax < nyi) THEN
      WRITE(6,'(a)') ' ERROR: parameter   kymax   too small'
      STOP ' error in parameter: KYMAX in COULEX too small'
    ELSE IF(kzmax < nzi) THEN
      WRITE(6,'(a)') ' ERROR: parameter   kzmax   too small'
      STOP ' error in parameter: KZMAX in COULEX too small'
    END IF
    
    !     initialize grid in Fourier space
    
    dkx=pi/(dx*REAL(nx,DP))
    dky=pi/(dy*REAL(ny,DP))
    dkz=pi/(dz*REAL(nz,DP))
    
    dxsp=dx*dy*dz
    dksp=dkx*dky*dkz
    WRITE(*,*) ' dkx,dky,dkz,dksp=',dkx,dky,dkz,dksp
    
    grnorm=SQRT(dxsp/dksp)
    fnorm=1.0D0/SQRT(REAL(nx*ny*nz,DP))
    !test      akmax=sqrt(3*(nx*nx)*dx*dx)+2.0
    !test      nxk=int(akmax/dkx)+1
    !test      if(nxk.gt.nx1) nxk=nx1
    nxk=nx1
    
    !     built Greens function in Fourier space
    !     by Fourier transformation from real space
    
    ikzero = nxy1*(nz-1)+nxi*(ny-1)+nx
    write(*,*) ' nzi,nyi,nxi,nx,ny,nz,ikzero=',nzi,nyi,nxi,nx,ny,nz,ikzero
    ii=0
    xz1=-nz*dz
    DO i3=1,nzi
      xz1=xz1+dz
      xz2=xz1*xz1
      xy1=-ny*dy
      DO i2=1,nyi
        xy1=xy1+dy
        xy2=xy1*xy1
        xx1=-nx*dx
        DO i1=1,nxi
          xx1=xx1+dx
          xx2=xx1*xx1
          ak2=xx2+xy2+xz2
          ii=ii+1
    !        write(*,*) ' i1,i2,i3,ii=',i1,i2,i3,ii
          IF(ii /= ikzero) THEN
            akv2r(ii) =  1D0/SQRT(ak2)
          ELSE
    !              akv2r(ii) = (6D0*pi/(dx*dy*dz))**(1D0/3D0)  ! spherical approx
    !              akv2r(ii) = 1.19003868*(dx*dy*dz)**(-1D0/3D0)
            akv2r(ii) = 2.34D0*1.19003868D0*(dx*dy*dz)**(-1D0/3D0)  ! empirical
          END IF
          akv2i(ii) = 0D0
        END DO
      END DO
    END DO
    nksp=ii

    rhokr_d = akv2r
    rhoki_d = akv2i
    
    CALL cufourf

    akv2r_d = rhokr_d
    akv2i_d = rhoki_d
    
    RETURN
    END SUBROUTINE cufftinp
    
        
    attributes(host) subroutine cucoufou2()
    use cutools
    use params , only: kdfull2
    implicit none
    
    !real(DP) :: realtest(kdred),aimagtest(kdred)
    !integer  :: ot
    INTEGER :: cudastatus
    real(DP) :: save2,tendcufourf,tcufoufor
    
    
    !call rhofld(rhointest,testrhokr,testrhoki)
    
    !write(*,*) 'cufourf'; call flush(6)
    !CALL cpu_time(tcufoufor)
    call cufourf
    !CALL cpu_time(tendcufourf)
    !write(*,*) "time cufourf",tendcufourf-tcufoufor
    
    
    
    !CALL fftx(testrhokr,testrhoki)
    !CALL ffty(testrhokr,testrhoki)
    !CALL fftz(testrhokr,testrhoki)
    
    !DO ot=1,kdred
    !  testrhokr(ot)=grnorm*fnorm*testrhokr(ot)
    !  testrhoki(ot)=grnorm*fnorm*testrhoki(ot)
    !END DO
    
    !rhokr_d = testrhokr
    !rhoki_d = testrhoki
    !write(*,*) 'cukineticfactor'; call flush(6)
    call cukineticfactor<<<ceiling(real(kdred)/256),256>>>(rhokr_d,rhoki_d,akv2r_d,akv2i_d)
    
    
    !DO ot=1,kdred
    ! save2     = akv2r(ot)*testrhokr(ot)+akv2i(ot)*testrhoki(ot)
    !  testrhoki(ot) = akv2r(ot)*testrhoki(ot)+akv2i(ot)*testrhokr(ot)
    !  testrhokr(ot) = SAVE2
    !END DO
    
    !DO ot=1,kdred
    !  testrhokr(ot)=fnorm/(8D0*grnorm)*pi**1.5D0*testrhokr(ot)
    !  testrhoki(ot)=fnorm/(8D0*grnorm)*pi**1.5D0*testrhoki(ot)
    !END DO
    
    !write(*,*) 'cufourb'; call flush(6)
    !CALL cpu_time(tcufoufor)
    call cufourb
    !CALL cpu_time(tendcufourf)
    !write(*,*) "time cufourb",tendcufourf-tcufoufor
    
    
    !CALL ffbz(testrhokr,testrhoki)
    !CALL ffby(testrhokr,testrhoki)
    !CALL ffbx(testrhokr,testrhoki)
    !write(*,*) 'endcufourb'; call flush(6)
    ! cudastatus = cudaDeviceSynchronize()
    !realtest = rhokr_d
    !aimagtest = rhoki_d
    !do ot = 1 , kdred
    !if(abs(testrhokr(ot))>1.d-7) then
    !if (ABS((realtest(ot)-testrhokr(ot))/testrhokr(ot))>1.d-5) &
    !	write(*,*) "error rhokr", realtest(ot),testrhokr(ot)
    !end if 
    !end do
    
    !do ot = 1 , kdred
    !if(abs(testrhoki(ot))>1.d-7) then
    !if (ABS((aimagtest(ot)-testrhoki(ot))/testrhoki(ot))>1.d-5) &
    !	write(*,*) "error rhoki", aimagtest(ot),testrhoki(ot)
    !end if
    !end do
    !STOP "cucoufou2"
    
    
    end subroutine cucoufou2
    
    
    attributes(host) subroutine reinitrhokr()
    use cutools
    implicit none
    call cumultarrays<<<ceiling(real(kdred)/256),256>>>(rhokr_d,rhoki_d,0.D0)
    end subroutine reinitrhokr
    
    
    
    attributes(host) subroutine cufourf()
    use cutools
    implicit none
    REAL(DP) :: tnorm
    REAL(DP) :: t1,t2
    
    tnorm=grnorm*fnorm
    !CALL cpu_time(t1)
    call cufftx
    !CALL cpu_time(t2)
    !write(*,*) "cufftx_time: ",t2-t1
    !CALL cpu_time(t1)
    call cuffty
    !CALL cpu_time(t2)
    !write(*,*) "cuffty_time: ",t2-t1
    !CALL cpu_time(t1)
    call cufftz
    !CALL cpu_time(t2)
    !write(*,*) "cufftz_time: ",t2-t1
    call cumultarrays<<<ceiling(real(kdred)/256),256>>>(rhokr_d,rhoki_d,tnorm)
    
    end subroutine cufourf
    
    
    attributes(host) subroutine cufourb()
    use cutools
    implicit none
    REAL(DP) :: tnorm
    REAL(DP) :: t1,t2
    tnorm=fnorm/(8D0*grnorm)*pi**1.5D0
    
    call cumultarrays<<<ceiling(real(kdred)/256),256>>>(rhokr_d,rhoki_d,tnorm)
    !CALL cpu_time(t1)
    call cubftz
    !CALL cpu_time(t2)
    !write(*,*) "cubftz_time: ",t2-t1
    !CALL cpu_time(t1)
    call cubfty
    !CALL cpu_time(t2)
    !write(*,*) "cubfty_time: ",t2-t1
    !CALL cpu_time(t1)
    call cubftx
    !CALL cpu_time(t2)
    !write(*,*) "cubftx_time: ",t2-t1
    
    end subroutine cufourb
    
    attributes(host) subroutine cufftx()
        use cukinetic
        use cutools
        use cudafor
        implicit none
        !real(DP) :: actrho(kdred),curho(kdred)
        INTEGER :: i0, i1, i2, i3, i30, ii, nx11
        !INTEGER :: cudastatus
        !real(DP) :: time_begin,time_end
            !actrho = rhokr_d
        !cudastatus = cudaDeviceSynchronize()
    
            !CALL cpu_time(time_begin)
        !nx11=nx1+1
        !i30=-nxy1
        !DO i3=1,nzi
        !  i30=i30+nxy1
        !  i0=i30-nxi
        !  DO i2=1,nyi
        !    i0=i0+nxi
    
        !         composition of the wave-function
        !         positive space
        !    ii=i0+nx-1
        !    DO i1=1,nx1
        !      ii=ii+1
        !      fftax(i1) = CMPLX(psxr(ii),0D0,DP)
        !	write(*,*) i1,ii,i0 + nx -1+i1,actrho(i0 + nx -1+i1)!,curho(i0 + nx -1+i1)
        !    END DO
        !         negative space
        !    ii=i0
        !    DO i1=nx11,nxi
        !      ii=ii+1
        !      fftax(i1) = CMPLX(psxr(ii),0D0,DP)
        !	write(*,*) i1,ii,i0 - nx -1+i1,actrho(i0 - nx -1+i1)!,curho(i0 - nx -1+i1)
        !    END DO
        !cudastatus = cudaDeviceSynchronize()
        !write(*,*) 'compXFFT'; call flush(6)
        CALL compXFFT<<<dim3(nzi,nyi,1),dim3(nxi,1,1)>>>(rhokr_d,fullfftax_d,nx,nxy1)
        !cudastatus = cudaDeviceSynchronize()
        !write(*,*) 'cufftExec'; call flush(6)
        CALL cufftExec(plan1dcoulsolvx,fullfftax_d,fullfftax_d,CUFFT_FORWARD)
        !cudastatus = cudaDeviceSynchronize()
        !write(*,*) 'decompXFFT'; call flush(6)
        CALL decompXFFT<<<dim3(nzi,nyi,1),dim3(nxi,1,1)>>>(rhokr_d,rhoki_d,fullfftax_d,nx,nxy1)
        !    CALL fftw_execute_dft(pforwx,fftax,fftax)
        !!        decomposition of the wave-function
        !!        positive space
        !    ii=i0+nx-1
        !    DO i1=1,nx1
        !      ii=ii+1
        !      psxr(ii) = REAL(fftax(i1),DP)
        !      psxi(ii) = AIMAG(fftax(i1))
        !    END DO
        !!        negative space
        !    ii=i0
        !    DO  i1=nx11,nxi
        !      ii=ii+1
        !      psxr(ii) = REAL(fftax(i1),DP)
        !      psxi(ii) = AIMAG(fftax(i1))
        !    END DO
        !  END DO
        !END DO
            !cudastatus = cudaDeviceSynchronize()
            !CALL cpu_time(time_end)
            !write(*,*) "end cufftx",time_end-time_begin; call flush(6)
    end subroutine cufftx
    
    attributes(host) subroutine cuffty
        use cukinetic
        use cutools
        implicit none
        INTEGER :: i0, i1, i2, i3, i30, ii, ny11
    
        !ny11=ny1+1
        !IF(nxk < nx1) THEN
        !  nxklo=nx-nxk+1
        !ELSE
        !  nxklo=1
        !END IF
        !nxkhi=nx+nxk-1
        !i30=-nxy1
    
        !DO i3=1,nzi
        !  i30=i30+nxy1
        !  DO i1=nx,nxi ! 1,nx1
        !    i0=i1+i30
        !!         composition of the wave-function
        !!         positive space
        !    ii=i0+nxi*(ny-2)
        !    DO i2=1,ny1
        !      ii=ii+nxi
        !      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
            !    END DO
        !!         negative space
        !    ii=i0-nxi
        !    DO i2=ny11,nyi
        !      ii=ii+nxi
        !      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
        !    END DO
        !write(*,*) "compffty"; call flush(6)
            CALL compYFFT<<<dim3(nzi,nx1,1),dim3(nyi,1,1)>>>(rhokr_d,rhoki_d,fftay_d,ny,nxi,nxy1)
        !write(*,*) "cufftExec"; call flush(6)
            CALL cufftExec(plan1dcoulsolvy,fftay_d,fftay_d,CUFFT_FORWARD)
            !write(*,*) "decompffty"; call flush(6)
        CALL decompYFFT<<<dim3(nzi,nx1,1),dim3(nyi,1,1)>>>(rhokr_d,rhoki_d,fftay_d,ny,nxi,nxy1)
        !!         decomposition of the wave-function
        !!         positive space
        !    ii=i0+nxi*(ny-2)
        !    DO i2=1,ny1
        !      ii=ii+nxi
        !      psxr(ii)=REAL(fftay(i2),DP)
        !      psxi(ii)=AIMAG(fftay(i2))
        !    END DO
        !!         negative space
        !    ii=i0-nxi
        !    DO i2=ny11,nyi
        !      ii=ii+nxi
        !      psxr(ii)=REAL(fftay(i2),DP)
        !      psxi(ii)=AIMAG(fftay(i2))
        !    END DO
        !  END DO
        !END DO
            
    end subroutine cuffty
    
    
    attributes(host) subroutine cufftz
        use cukinetic
        use cutools
        implicit none
    
        INTEGER :: lnxyf 
        INTEGER :: lnyf  
        INTEGER :: lnzh ! ,custat
        INTEGER :: i1, i2, i3, i3m, ind,ot
        !complex(DP) , allocatable :: actfftb(:)
        !complex(DP) , allocatable :: reffftb(:)
        !real(DP) , allocatable :: psxr(:),psxi(:)
    
    
        !allocate(actfftb(nx1*kffty*kfftz),reffftb(nx1*kffty*kfftz))
        !allocate(psxr(kdred),psxi(kdred))
        !psxr = rhokr_d
        !psxi = rhoki_d
        lnxyf = kfftx*kffty
        lnyf  = kfftx
        lnzh  = kfftz/2
    !	write(*,*) "allocated" ; call flush(6)
    !	DO i2=1,kffty
    !	  DO i1=nx,nxi ! 1,nx1
    !	    DO i3=1,kfftz
    !	    i3m   = MOD(i3+nzh,kfftz)+1
    !	      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
        !      fftb(i3m,i1) = CMPLX(psxr(ind),psxi(ind),DP)
    !	    reffftb(((i2-1)*nx1+i1-nx)*kfftz+i3m) = &
    !		CMPLX(psxr(ind),psxi(ind),DP)
    !	    END DO
    !	  END DO
    !	ENDDO !test 
           call compZFFT<<<dim3(nx1,kffty,1),dim3(kfftz,1,1)>>>(rhokr_d,rhoki_d,fftb_d,lnxyf,lnyf,lnzh,kfftz)
    
    !		actfftb = fftb_d
    !do ot = 1 , kfftz!nx1*kffty*kfftz
    !if(abs(reffftb(ot))>1.d-5) then
    !if (ABS((reffftb(ot)-actfftb(ot))/reffftb(ot))>1.d-5) &
    !	write(*,*) "error cufftz", abs(reffftb(ot)),abs(actfftb(ot))
    !end if 
    !end do
    ! call flush(6)
    
    !	DO i2 = 1 ,kffty
    !	DO i1 = nx,nxi
    !CALL fftw_execute_dft(pforwz,&
    !	reffftb( ((i2-1)*nx1+i1-nx)*kfftz:((i2-1)*nx1+i1-nx)*kfftz+kfftz),&
    !	reffftb(((i2-1)*nx1+i1-nx)*kfftz:((i2-1)*nx1+i1-nx)*kfftz+kfftz))
    !	END DO
    !	END DO
    
    
    
           call cufftExec(plan1dcoulsolvz,fftb_d,fftb_d,CUFFT_FORWARD)
    !	custat  = cudaDeviceSynchronize()
    !		actfftb = fftb_d
    !do ot = 1 , kfftz!nx1*kffty*kfftz
    !if(abs(reffftb(ot))>1.d-5) then
    !if (ABS((reffftb(ot)-actfftb(ot))/reffftb(ot))>1.d-5) &
    !	write(*,*) "error cufftz", abs(reffftb(ot)),abs(actfftb(ot))
    !end if 
    !end do
    ! call flush(6)
    !fftb_d = reffftb
    !custat  = cudaDeviceSynchronize()
    !write(*,*) "here" ; call flush(6)
          call decompZFFT<<<dim3(nx1,kffty,1),dim3(kfftz,1,1)>>>(rhokr_d,rhoki_d,fftb_d,lnxyf,lnyf,lnzh,kfftz)
        !  DO i3=1,kfftz
        !    i3m   = MOD(i3+nzh,kfftz)+1
        !    DO i1=nx,nxi ! 1,nx1
        !      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
        !      psxr(ind)= REAL(fftb(i3m,i1),DP)
        !      psxi(ind)= AIMAG(fftb(i3m,i1))
        !    END DO
        !  END DO
        !END DO
    !write(*,*) "here" ; call flush(6)
    !deallocate(actfftb,reffftb,psxr,psxi)
            
    end subroutine  cufftz
    
    attributes(host) subroutine cubftz
        use cukinetic
        use cutools
        implicit none
    
        INTEGER :: lnxyf 
        INTEGER :: lnyf  
        INTEGER :: lnzh  
        INTEGER :: i1, i2, i3, i3m, ind
        complex(DP) :: fftest(nzi)
    
        lnxyf = kfftx*kffty
        lnyf  = kfftx
        lnzh  = kfftz/2
    
        !DO i2=1,kffty
        !  DO i1=nx,nxi ! 1,nx1
        !    DO i3=1,kfftz
        !    i3m   = MOD(i3+nzh,kfftz)+1
        !      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
        !      fftb(i3m,i1) = CMPLX(psxr(ind),psxi(ind),DP)
        !    END DO
        !  END DO
           call compZFFT<<<dim3(nx1,kffty,1),dim3(kfftz,1,1)>>>(rhokr_d,rhoki_d,fftb_d,lnxyf,lnyf,lnzh,kfftz)
           call cufftExec(plan1dcoulsolvz,fftb_d,fftb_d,CUFFT_INVERSE)
           call decompZFFT<<<dim3(nx1,kffty,1),dim3(kfftz,1,1)>>>(rhokr_d,rhoki_d,fftb_d,lnxyf,lnyf,lnzh,kfftz)
        !  DO i3=1,kfftz
        !    i3m   = MOD(i3+nzh,kfftz)+1
        !    DO i1=nx,nxi ! 1,nx1
        !      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
        !      psxr(ind)= REAL(fftb(i3m,i1),DP)
        !      psxi(ind)= AIMAG(fftb(i3m,i1))
        !    END DO
        !  END DO
        !END DO
        
    end subroutine  cubftz
    
    attributes(host) subroutine cubfty
        use cukinetic
        use cutools
        implicit none
        INTEGER :: i0, i1, i2, i3, i30, ii, ny11
    
        ny11=ny1+1
    
        !IF(nxk < nx1) THEN
        !  nxklo=nx-nxk+1
        !ELSE
        !  nxklo=1
        !END IF
        !nxkhi=nx+nxk-1
    
        !i30=-nxy1
    
        !DO i3=1,nz1
        !  i30=i30+nxy1
        !  DO i1=nx,nxi ! 1,nx1
        !    i0=i1+i30
        !!         composition of the complex wave-function
        !!         positive space
        !    ii=i0+nxi*(ny-2)
        !    DO i2=1,ny1
        !      ii=ii+nxi
        !      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
        !    END DO
        !!         negative space
        !    ii=i0-nxi
        !    DO i2=ny11,nyi
        !      ii=ii+nxi
        !      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
        !    END DO
        CALL compYBFT<<<dim3(nz1,nx1,1),dim3(nyi,1,1)>>>(rhokr_d,rhoki_d,fftayback_d,ny,nxi,nxy1)
        CALL cufftExec(plan1dcoulsolvyback,fftayback_d,fftayback_d,CUFFT_INVERSE)
        CALL decompYBFT<<<dim3(nz1,nx1,1),dim3(nyi,1,1)>>>(rhokr_d,rhoki_d,fftayback_d,ny,nxi,nxy1)
        !!         decomposition of the wave-function
        !!         positive space
        !    ii=i0+nxi*(ny-2)
        !    DO i2=1,ny1
        !      ii=ii+nxi
        !      psxr(ii)=REAL(fftay(i2),DP)
        !      psxi(ii)=AIMAG(fftay(i2))
        !    END DO
        !!         negative space
        !    ii=i0-nxi
        !    DO i2=ny11,nyi
        !      ii=ii+nxi
        !      psxr(ii)=REAL(fftay(i2),DP)
        !      psxi(ii)=AIMAG(fftay(i2))
        !    END DO
        !  END DO
        !END DO
    
    
    
    end subroutine cubfty
    
    attributes(host) subroutine cubftx()
        use cukinetic
        use cutools
        use cudafor
        implicit none
        INTEGER :: i0, i1, i2, i3, i30, ii, nx11
    
    !nx11=nx1+1
    !i30=-nxy1
    
    
    !DO i3=1,nz1
    !  i30=i30+nxy1
    !  i0=i30-nxi
    !  DO i2=1,ny1
    !    i0=i0+nxi
    !!         composition of the complex wave-function
    !!        positive space
    !    ii=i0+nx-1
    !    DO i1=1,nx1
    !      ii=ii+1
    !      fftax(i1)=CMPLX(psxr(ii),psxi(ii),DP)
    !    END DO
    !!        negative space
    !    DO i1=nx11,nxi
    !      fftax(i1)=CONJG(fftax(nxi-i1+2))
    !    END DO
        CALL compXBFT<<<dim3(nz1,ny1,1),dim3(nxi,1,1)>>>(rhokr_d,rhoki_d,fullfftax_d,nx,nxy1)
        CALL cufftExec(plan1dcoulsolvx,fullfftax_d,fullfftax_d,CUFFT_INVERSE)
        CALL decompXFFT<<<dim3(nz1,ny1,1),dim3(nxi,1,1)>>>(rhokr_d,rhoki_d,fullfftax_d,nx,nxy1)
    !!         decomposition of the inverse transformed wave-function
    !!         positive space
    !   ii=i0+nx-1
    !    DO i1=1,nx1
    !      ii=ii+1
    !      psxr(ii)=REAL(fftax(i1),DP)
    !      psxi(ii)=AIMAG(fftax(i1))
    !    END DO
    !!         negative space
    !    ii=i0
    !    DO i1=nx11,nxi
    !      ii=ii+1
    !      psxr(ii)=REAL(fftax(i1),DP)
    !      psxi(ii)=AIMAG(fftax(i1))
    !    END DO
    !  END DO  
    !END DO
    
    end subroutine cubftx
    
    
    
    SUBROUTINE cucoulsolv_end()
    
    deallocate(akv2r,akv2i)
    deallocate(testrhokr,testrhoki)
    deallocate(fftsizex,inembedx,onembedx)
    deallocate(fftsizey,inembedy,onembedy)
    deallocate(fftsizez,inembedz,onembedz)
    deallocate(fftax_d,fftay_d,fftb_d)
    deallocate(fftayback_d)
    deallocate(fullfftax_d)
    deallocate(rhokr_d,rhoki_d)
    deallocate(akv2r_d,akv2i_d)
    RETURN
    END SUBROUTINE cucoulsolv_end
    
    END MODULE cucoulsolv
    