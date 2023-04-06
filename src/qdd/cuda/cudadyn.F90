
MODULE cudadyn
  use params, only: DP,kdfull2,nstate
  use cukinetic
  use cudafor
  implicit none
  !--------------------------------------------
  complex(DP) ,device ,allocatable :: psi_d(:,:)
  complex(DP) ,device ,allocatable :: kpsi_d(:,:)
  complex(DP) ,device ,allocatable :: psiwork_d(:)
  real(DP) ,device ,allocatable :: rho_d(:)
  real(DP) ,device ,allocatable :: aloc_d(:)
  real(DP) ,device ,allocatable :: chpcoul_d(:)
  real(DP) ,device ,allocatable :: chpdft_d(:)
  real(DP) ,device ,allocatable :: rhon_d(:)
  real(DP) ,device ,allocatable :: rhojel_d(:)
  real(DP) ,device ,allocatable :: occup_d(:)
  real(DP) ,device ,allocatable :: ekin_d(:)
  real(DP) ,device ,allocatable :: epot_d(:)
  real(DP) ,device ,allocatable :: enonloc_d(:)
  real(DP) ,device ,allocatable :: potion_d(:)
  !--------------------------------------------
  !REAL(DP),device,allocatable   :: rhosp_d(:),chpdftsp_d(:)
  REAL(DP),device,allocatable   :: coulsum_d(:),couldif_d(:)
  REAL(DP),device,allocatable   :: rho1_d(:),rho2_d(:)
  REAL(DP),device,allocatable   :: rhospu_d(:),rhospd_d(:)
  REAL(DP),device,allocatable   :: chpdftspu_d(:),chpdftspd_d(:)
  !--------------------------------------------
  REAL(Dp),device		:: enrear_d,ecback_d,ecrho_d
  type(dim3) :: grid, tBlock
  contains
  
  subroutine init_device_psi()
          use params, only: rhojel,potion,nion2,dvol,numspin,ispin,dx, dy, dz, nx2, ny2, nz2, dt1, h2m
          use cuparams, only: dvol_d,ispin_d
          use cukinetic
          use cucoulsolv
          implicit none

    allocate(psi_d(kdfull2,nstate))
    allocate(kpsi_d(kdfull2,nstate))
    allocate(psiwork_d(kdfull2))
    allocate(ekin_d(nstate),epot_d(nstate),enonloc_d(nstate))
    allocate(rho_d(2*kdfull2))
    allocate(aloc_d(2*kdfull2))
    allocate(chpcoul_d(kdfull2))
    allocate(chpdft_d(2*kdfull2))
    allocate(rhon_d(kdfull2))
    allocate(rhojel_d(kdfull2))
    allocate(occup_d(nstate))
    allocate(coulsum_d(kdfull2),couldif_d(kdfull2))
    allocate(rho1_d(kdfull2),rho2_d(kdfull2))
    allocate(rhospu_d(2*kdfull2),rhospd_d(2*kdfull2))
    allocate(chpdftspu_d(2*kdfull2),chpdftspd_d(2*kdfull2))
    if (nion2 /=0) then
      allocate(potion_d(kdfull2)) 
      potion_d = potion	
    end if
    !------
    tBlock = dim3(1024,1,1) !warp size max = 32
    grid   = dim3(ceiling(real(kdfull2)/tBlock%x), 1, 1)
    rhojel_d = rhojel
    allocate(ispin_d(nstate))
    ispin_d = ispin(:nstate)
    write(*,*) "init_cuda_fft ",nstate; call flush(6)
    call init_cuda_fft(dx, dy, dz, nx2, ny2, nz2, dt1, h2m)
    write(*,*) "cuinit_coul ",nstate; call flush(6)
    call cuinit_coul(dx, dy, dz, nx2, ny2, nz2)
    !------
  end subroutine init_device_psi
  
  attributes(host) subroutine cudyn_mfield(actpsi_d)
    use params!, only: DP,kdfull2,kstate,occup
    use cutools
    use coulsolv
    implicit none
    real(DP)          ,allocatable  :: rho(:),rhocu(:) !test
    real(DP)	  ,allocatable  :: aloc(:),aloccu(:) !test
    complex(DP)  , optional , device   :: actpsi_d(kdfull2,nstate)
          complex(Dp) , allocatable :: psi(:,:) !test
    integer				:: i,custat
          REAL(DP)        ::   tout,tin
    !----------------------------------------------------
    allocate(rho(2*kdfull2),rhocu(2*kdfull2))
    allocate(aloc(2*kdfull2),aloccu(2*kdfull2))
    allocate(psi(kdfull2,kstate))
    !IF(numspin/=2) THEN
    !	STOP "numspin must be equal 2 for this cuda implementation"
    !END IF
    !ALLOCATE(rhokr(kdred),rhoki(kdred))
          custat  = cudaDeviceSynchronize()
          !CALL cpu_time(tin)
    IF (present(actpsi_d)) THEN
      call cucalcrho(actpsi_d,occup_d,rho_d)
    ELSE
      call cucalcrho(psi_d,occup_d,rho_d)
    END IF
    IF (present(actpsi_d)) THEN
    	psi(:,:nstate) = actpsi_d(:,:nstate)
    ELSE
    	psi(:,:nstate) = psi_d(:,:nstate)
    END IF
    !custat  = cudaDeviceSynchronize()
          !CALL cpu_time(tout)
          !write(*,*) "time cucalcrho: ",tout-tin
          !CALL cpu_time(tin)
          CALL calcrho(rho,psi)
    !CALL cpu_time(tout)
          !write(*,*) "time calcrho: ",tout-tin
          !custat  = cudaDeviceSynchronize()
          !CALL cpu_time(tin)
          !custat  = cudaDeviceSynchronize()
          call cucoul_mfield()
          !custat  = cudaDeviceSynchronize()
          !CALL cpu_time(tout)
          !write(*,*) "time cucoul_mfield: ",tout-tin
          !CALL cpu_time(tin)
          call coul_mfield(rho)
          !CALL cpu_time(tout)
          !write(*,*) "time coul_mfield: ",tout-tin
          !CALL cpu_time(tin)
          call cucalclocal()
          !CALL cpu_time(tout)
          !write(*,*) "time cucalclocal: ",tout-tin
          CALL calclocal(rho,aloc)




          !CALL cpu_time(tin)
  IF(ifsicp > 0 .AND.ifsicp <= 6) THEN
    CALL calc_sic(rho,aloc,psi)
    CALL cucalc_sic()
  #if(twostsic)
  ELSE IF(ifsicp >= 7)THEN
    STOP "twostsic ifsicp >= 7 not for cuda"
  #endif
  ELSE IF(ifsicp == 6) THEN
    STOP ' that kind of SIC not valid for dynamics'
  END IF
          !CALL cpu_time(tout)
          !write(*,*) "time cucalc_sic: ",tout-tin
    !custat  = cudaDeviceSynchronize()
    !rhocu = rho_d
    aloccu = aloc_d
    !custat  = cudaDeviceSynchronize()
    !i=1
    !do while (i<=size(rho,1))
    !if (ABS(rhocu(i)-rho(i))>1.d-8) then
    !	write(*,'(a,I12,2x,3F12.8)') 'error rho ',i,rhocu(i),rho(i),ABS(rhocu(i)-rho(i))
    !end if
    !i = i+1
    !end do
    i=1
    do while (i<=size(aloc,1))
    if (ABS(aloccu(i)-aloc(i))>1.d-8) then
    	write(*,'(a,I12,2x,3F12.8)') 'error ',i,aloccu(i),aloc(i),ABS(aloccu(i)-aloc(i))
    end if
    i = i+1
    end do
    !STOP "Bazinga !"
  end subroutine cudyn_mfield
  
  subroutine cuinfo(eout,it)
  USE params
  USE cutools
  !----------------------
  COMPLEX(DP), INTENT(IN)		  :: eout(kstate,kstate)
  INTEGER, INTENT(IN)                      :: it
  !----------------------
  REAL(DP) , DEVICE , allocatable :: spnorm_d(:)
  INTEGER :: nb
  !---------------------
  
  allocate(spnorm_d(nstate))
  CALL cudyn_mfield()
  
  CALL cuwfovlptot<<<nstate,256,256*128>>>(psi_d,psi_d,spnorm_d)
  CALL cucalc_ekin()
  CALL cucalc_epot()
  epotsp(:nstate) = epot_d(:nstate)
  ekinsp(:nstate) = ekin_d(:nstate)
  enonlo(:nstate) = enonloc_d(:nstate)
  !DO nb=1,nstate
  !WRITE(6,'(2(a,i2),5(a,f9.5))') 'level:',nrel2abs(nb), &
  !'  ispin=',ispin(nrel2abs(nb)),' occup=',occup(nb), &
  !'  ekin='  &
  !,ekin,'  epot=',ehilf,'  esp=',amoy(nb) ,'  enonlo=', enonlo(nb)
  !END DO
  
  !IF(jesc > 0 .AND. jnorms>0 .AND. MOD(it,jnorms) == 0) THEN
      spnorm(:nstate) = spnorm_d(:nstate) 
     ! CALL safeopen(806,it,jnorms,'pescOrb')
      WRITE(*,*) "spnorm", it,1D0-spnorm(1:nstate);
      WRITE(*,*) "ekinsp", it,ekinsp(1:nstate);
      WRITE(*,*) "epotsp", it,epotsp(1:nstate);
      WRITE(*,*) "enonlo", it,enonlo(1:nstate);
     ! CALL flush(806)
        CALL flush(6)
  !END IF
  !  CALL safeopen(808,it,jnorms,'pproba')
  !  sCALL probab(psi)
  !  CALL flush(808)
  
  tstinf = jstinf > 0 .AND. MOD(it,jstinf)==0
  !tstinf=.FALSE.
  IF(tstinf) then  !! pour calculer spvariance
    STOP "not implemented for cuda" 
  END IF
  
  #if(extended)
  IF(jstboostinv>0 .AND. MOD(it,jstboostinv)==0) THEN
    STOP "not implemented for cuda in extended"
  END IF
  #endif
  
  !call calc_estar(psi,iss,estar(iss),estarETF(iss))
  
  !---------rearrangement and background Coulomb energy
  !ecback=0D0
  !ecrho=0D0
  IF(nion2 /= 0) THEN
    !DO ind=1,nxyz
    !  ecback=ecback-rho(ind)*potion(ind)
     call cugetecbackion<<<1,256,256*64>>>(rho_d,potion_d,ecback_d)
     ! ecrho=ecrho+rho(ind)*(chpcoul(ind)-potion(ind))
     call cugetecrhoion<<<1,256,256*64>>>(rho_d,chpcoul_d,potion_d,ecrho_d)
   ! END DO
  !  WRITE(*,*) ' ECBACK loop:',ecback*dvol/2.0
  ELSE ! jellium case
   ! DO ind=1,nxyz
    !  ecback=ecback-rhojel(ind)*chpcoul(ind)
    call cugetecback<<<1,256,256*64>>>(rhojel_d,chpcoul_d,ecback_d)
      !ecrho=ecrho+rho(ind)*chpcoul(ind)
    call cugetecrho<<<1,256,256*64>>>(rho_d,chpcoul_d,ecrho_d)
   ! END DO
  END IF
  !ecback=ecback_d
  !ecback=ecback*dvol/2.0D0
  !ecrho=ecrho_d
  !ecrho=ecrho*dvol/2.0D0
  !eshell=eshell/2.0D0
  
  end subroutine cuinfo
  
  subroutine cucalc_epot()
  ! compute epot_d and enonloc_d inside the GPU
  ! the called kernels depend on the psp type and if numspin==2 or 1
  USE params, only: dvol,nstate,numspin,ipsptyp
  USE cutools
    IF(ipsptyp == 1) THEN
      call cutotnonlocalc(psi_d,kpsi_d,0)
      call cucalc_epotnonloc<<<nstate,256,256*16>>>(psi_d,kpsi_d,aloc_d,epot_d)
      call cugetnonlocenergies<<<nstate,256,256*8>>>(psi_d,kpsi_d,enonloc_d)
    ELSE
        call cucalc_epotloc<<<nstate,256,256*16>>>(psi_d,aloc_d,epot_d)
    END IF
  
  end subroutine cucalc_epot
  
  subroutine cucalcrho(inpsi_d,inoccup_d,outrho_d)
  use params 
  use cutools
  implicit none 
  COMPLEX(DP), DEVICE ,INTENT(IN) :: inpsi_d(kdfull2,nstate)
  REAL(DP), DEVICE ,INTENT(IN) :: inoccup_d(nstate)
  REAL(DP), DEVICE ,INTENT(IN OUT) :: outrho_d(2*kdfull2)
  !----------------------------
  
  IF (numspin==2) THEN
    call cucalcrhospindep<<<kdfull2/256,256>>>(inpsi_d,inoccup_d,outrho_d)
  ELSE 
    call cucalcrhonospin<<<kdfull2/256,256>>>(inpsi_d,inoccup_d,outrho_d) 
  END IF
  
  end subroutine cucalcrho
  
  
  subroutine cucalc_current(current_d)
  use params 
  use cutools
  implicit none 
  REAL(DP), DEVICE, INTENT(INOUT) :: current_d(kdfull2,3)
  !----------------------------
  COMPLEX(DP)	     , DEVICE, ALLOCATABLE :: dq0(:,:)
  !----------------------------
  INTEGER :: nb
  !----------------------------
  allocate(dq0(kdfull2,3))
  !----------------------------
  call reset2dr<<<dim3(kdfull2/256,3,1),dim3(256,1,1)>>>(current_d)
  DO nb=1,nstate
    call cuxyzgradient_rspace(psi_d(:,nb),dq0)
    call cuaddgradtocurr<<<dim3(kdfull2/256,3,1),dim3(256,1,1)>>>(current_d,dq0,psi_d(:,nb),occup_d(nb))
  END DO
  end subroutine cucalc_current
  
  subroutine cuxyzgradient_rspace(inpsi,outgrad)
  use params 
  use cutools
  use cukinetic
  implicit none
  COMPLEX(DP), INTENT(IN) , DEVICE :: inpsi(kdfull2)
  COMPLEX(DP), INTENT(OUT), DEVICE :: outgrad(kdfull2,3)
  !----------------------------
  REAL(DP)   :: dkx,dky,dkz
  !----------------------------
  dkx=pi/(dx*nx)
  dky=pi/(dy*ny)
  dkz=pi/(dz*nz)
  !----------------------------
  call compXYZgrad<<<dim3(nz2,ny2,1),dim3(nx2,1,1)>>>(inpsi,fullfftx_d,fullfftay_d,fullfftaz_d)
  CALL cufftExec(planmany1dx,fullfftx_d,fullfftx_d,CUFFT_FORWARD)
  CALL cufftExec(planmany1dy,fullfftay_d,fullfftay_d,CUFFT_FORWARD)
  CALL cufftExec(planmany1dz,fullfftaz_d,fullfftaz_d,CUFFT_FORWARD)
  call cuapplygrad<<<dim3(nz2,ny2,1),dim3(nx2,1,1)>>>(fullfftx_d,fullfftay_d,fullfftaz_d,dkx,dky,dkz)
  CALL cufftExec(planmany1dx,fullfftx_d,fullfftx_d,CUFFT_INVERSE)
  CALL cufftExec(planmany1dy,fullfftay_d,fullfftay_d,CUFFT_INVERSE)
  CALL cufftExec(planmany1dz,fullfftaz_d,fullfftaz_d,CUFFT_INVERSE)
  call decompXYZgrad<<<dim3(nz2,ny2,1),dim3(nx2,1,1)>>>(outgrad,fullfftx_d,fullfftay_d,fullfftaz_d)
  end subroutine cuxyzgradient_rspace
  
  subroutine cu3gradient_rspace(inpsi,outgrad)
  use params 
  use cutools
  use cukinetic
  implicit none
  COMPLEX(DP), INTENT(IN) , DEVICE :: inpsi(kdfull2,3)
  COMPLEX(DP), INTENT(OUT), DEVICE :: outgrad(kdfull2,3)
  !----------------------------
  REAL(DP)   :: dkx,dky,dkz
  !----------------------------
  dkx=pi/(dx*nx)
  dky=pi/(dy*ny)
  dkz=pi/(dz*nz)
  !----------------------------
  call comp3grad<<<dim3(nz2,ny2,1),dim3(nx2,1,1)>>>(inpsi,fullfftx_d,fullfftay_d,fullfftaz_d)
  CALL cufftExec(planmany1dx,fullfftx_d,fullfftx_d,CUFFT_FORWARD)
  CALL cufftExec(planmany1dy,fullfftay_d,fullfftay_d,CUFFT_FORWARD)
  CALL cufftExec(planmany1dz,fullfftaz_d,fullfftaz_d,CUFFT_FORWARD)
  call cuapplygrad<<<dim3(nz2,ny2,1),dim3(nx2,1,1)>>>(fullfftx_d,fullfftay_d,fullfftaz_d,dkx,dky,dkz)
  CALL cufftExec(planmany1dx,fullfftx_d,fullfftx_d,CUFFT_INVERSE)
  CALL cufftExec(planmany1dy,fullfftay_d,fullfftay_d,CUFFT_INVERSE)
  CALL cufftExec(planmany1dz,fullfftaz_d,fullfftaz_d,CUFFT_INVERSE)
  call decompXYZgrad<<<dim3(nz2,ny2,1),dim3(nx2,1,1)>>>(outgrad,fullfftx_d,fullfftay_d,fullfftaz_d)
  end subroutine cu3gradient_rspace
  
  attributes(host) subroutine cufftf()
  use params, only: kxbox,kybox,kybox
  use cukinetic
  use cutools
  implicit none
  real(DP)   :: tnorm
  integer    :: ii
  !------------------------
  tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(kxbox*kybox*kybox,DP))
  DO ii = 1 ,nstate
  !write(*,*) "go cucopy1dto3d"; call flush(6)
  call cucopy1dto3d<<<gridfft3d,tBlockfft3d>>>(psi_d,ffta_d,ii)
  !write(*,*) "go cufftExecZ2Z"; call flush(6)
  call cufftExecZ2Z(plan3d, ffta_d(:,:,:), ffta_d(:,:,:), CUFFT_FORWARD)
  !write(*,*) "go cucopy3dto1d"; call flush(6)
  call cucopy3dto1d<<<gridfft3d,tBlockfft3d>>>(kpsi_d,ffta_d,tnorm,ii) 
  END DO
  end subroutine cufftf
  
  attributes(host) subroutine cufftb()
  use params, only: kxbox,kybox,kybox
  use cukinetic
  use cutools
  implicit none
  real(DP)   :: tnorm
  integer    :: ii
  !------------------------
  tnorm =SQRT(8D0*pi*pi*pi)/SQRT(REAL(kxbox*kybox*kybox,DP))
  DO ii = 1 ,nstate
  call cucopy1dto3d<<<gridfft3d,tBlockfft3d>>>(kpsi_d,ffta_d,ii)
  call cufftExecZ2Z(plan3d, ffta_d(:,:,:), ffta_d(:,:,:), CUFFT_INVERSE)
  call cucopy3dto1d<<<gridfft3d,tBlockfft3d>>>(kpsi_d,ffta_d,tnorm,ii) 
  END DO
  end subroutine cufftb
  
  attributes(host) subroutine cucalc_ekin()
  use params, only: ekinsp,dx,dy,dz,PI
  use cutools
  use cukinetic, only: akv_d
  implicit none
  integer :: nb
  REAL(DP) :: sum0ex
  !COMPLEX(DP) :: psi(kdfull2),kpsi(kdfull2)
  
  

  CALL cufftf()
  !do nb = 1 , nstate
  !psi(:) = psi_d(:,nb)
  !CALL fftf(psi,kpsi)
  !kpsi_d(:,nb) = kpsi(:)
  !end do
  !!sum0 = 0D0
  !sumk = 0D0
  !DO ii=1,kdfull2
  !  vol   = REAL(psi2(ii),DP)*REAL(psi2(ii),DP) +AIMAG(psi2(ii))*AIMAG(psi2(ii))
    !!sum0  = vol + sum0
    !sumk  = vol*akv(ii) + sumk
  !END DO
  sum0ex = 1D0/((2D0*PI)**3*dx*dy*dz)
  call cukineticenerg<<<nstate,256,256*64>>>(kpsi_d,akv_d,ekin_d,sum0ex)
  !ekinout = sumk/sum0ex
  !ekinsp(:nstate) = ekin_d(:nstate)
  
  !WRITE(6,*) ' sum0,sum0ex=',sum0,sum0ex
  #if(findiff|numerov)
  STOP "not for CUDA"
  #endif
  
  end subroutine cucalc_ekin
  
  
  attributes(host) subroutine cucoul_mfield()
    USE cucoulsolv
    USE cutools
    USE params, only: e2,nx2,ny2,nz2
    implicit none
  
    call reinitrhokr
    !if(.not.allocated(rhointest)) allocate(rhointest(kdfull2))
    !rhointest(:kdfull2) = rho_d(:kdfull2)
    call curhofld<<<gridrhofld,tBlockrhofld>>>(rho_d,rhokr_d,rhoki_d,nx2,ny2,nz2)
    call cufalr
    call curesult<<<gridresult,tBlockresult>>>(chpcoul_d,rhokr_d,e2,nx2,ny2,nz2)

  #if(findiff|numerov)
    !call solv_poisson(rho,chpcoul,kdfull2)
    STOP "solv_poisson not implemented for cuda..."
  #endif
  end subroutine cucoul_mfield
  
  attributes(host) subroutine cucalc_sic()
  USE params
  IMPLICIT NONE

  #if(extended)
  IF (ifreezekspot == 1 .AND. tfs > 0) RETURN
  #endif
  SELECT CASE(ifsicp)
    CASE(1)
     ! CALL calc_sicgam(rho,aloc,q0)
    STOP "calc_sicgam not yet implemented for cuda"
    CASE(2)
      !CALL calc_adsic(rho,aloc)
          call cucalc_adsic()
    CASE(3)
     ! CALL calc_slater(rho,aloc,q0)
    STOP "calc_slater not yet implemented for cuda"
    CASE(4)
     ! CALL calc_sickli(rho,aloc,q0)
    STOP "calc_sickli not yet implemented for cuda"
    CASE DEFAULT
      RETURN
  END SELECT
  
  end subroutine cucalc_sic
  
  
  attributes(host) subroutine cucalc_adsic
  USE params
  USE cutools
  USE cucoulsolv
  IMPLICIT NONE
  
  
  INTEGER  :: npartdw, npartup, npartto,ind,i,idx
  REAL(DP) :: facuph, facdwh ,facdw, facup, fac
  REAL(DP) :: enrear1,enrear2
  REAL(DP) ,device :: enrear1_d,enrear2_d
  !test
  !real(DP) ,allocatable :: rho(:),aloc(:)
  !REAL(DP),DIMENSION(:),ALLOCATABLE :: rhosp,chpdftsp,coulsum,couldif
  !REAL(DP),DIMENSION(:),ALLOCATABLE :: rho1,rho2
  !REAL(DP),DIMENSION(:),ALLOCATABLE :: rhospu,rhospd,chpdftspu,chpdftspd
  !test
  
  
  CALL act_part_num(npartup,npartdw,npartto)
   ! allocate(rho(2*kdfull2))
   ! allocate(aloc(2*kdfull2))
   ! ALLOCATE(rhospu(2*kdfull2),rhospd(2*kdfull2),  &
   !                 chpdftspu(2*kdfull2),chpdftspd(2*kdfull2))
  
   !rho = rho_d
   !aloc = aloc_d
  IF(numspin==2) THEN
  
  
    IF(npartup > 0) THEN
      facuph = 0.5D0/npartup
    call constructrhospu<<<ceiling(real(nxyz)/256),256>>>(rhospu_d,rho_d,facuph)
  
    !  DO ind=1,nxyz
    !    rhospu(ind)=rho(ind)*(1D0+rho(ind+nxyz))*facuph
    !    rhospu(ind+nxyz)=1D0
    !  END DO
    END IF
  
    IF(npartdw > 0) THEN
      facdwh = 0.5D0/npartdw
    call constructrhospd<<<ceiling(real(nxyz)/256),256>>>(rhospd_d,rho_d,facdwh)
     ! DO ind=1,nxyz
     !   rhospd(ind)=rho(ind)*(1D0-rho(ind+nxyz))*facdwh
     !   rhospd(ind+nxyz)=-1D0
     ! END DO
    END IF
  
    IF(npartup > 0) THEN
      CALL cucalc_lda<<<1,256,256*8>>>(rhospu_d,chpdftspu_d,dvol,enrear1_d)
      !CALL calc_lda(rhospu,chpdftspu)
      ! enrear1=enrear
    END IF
    IF(npartdw > 0) THEN
      CALL cucalc_lda<<<1,256,256*8>>>(rhospd_d,chpdftspd_d,dvol,enrear2_d)
      !CALL calc_lda(rhospd,chpdftspd)
      !enrear2=enrear
    END IF
   !write(*,*) enrear1,enrear2
    enrear1  = enrear1_d
    enrear2  = enrear2_d
   !write(*,*) enrear1,enrear2
    enrear   = enrear-enrear1*npartup-enrear2*npartdw
    call sicalocaverdensity<<<ceiling(real(nxyz)/256),256>>>(aloc_d,chpdftspu_d,chpdftspd_d)
    !IF(npartup > 0) THEN
    !  DO ind=1,nxyz
    !    aloc(ind) = aloc(ind) - chpdftspu(ind)
    !  END DO
    !END IF
    !IF(npartdw > 0) THEN
    !  DO idx=nxyz+1,2*nxyz
    !    aloc(idx) = aloc(idx) - chpdftspd(idx)
    !  END DO
    !END IF
  ELSE
      facuph = 1D0/npartto
    call constructrhospdeg<<<ceiling(real(nxyz)/256),256>>>(rhospu_d,rho_d,facuph)
      !DO ind=1,nxyz
      !rhosp(ind)=rho(ind)*factotal
      !rhosp(ind+nxyz)=1D0
      !END DO
  
      CALL cucalc_lda<<<1,256,256*8>>>(rhospu_d,chpdftspu_d,dvol,enrear1_d)
      !CALL calc_lda(rhosp,chpdftsp)
      
      call sicalocaverdensitydeg<<<ceiling(real(nxyz)/256),256>>>(aloc_d,chpdftspu_d)
      !DO ind=1,size
      !aloc(ind) = aloc(ind) - chpdftsp(ind)
      !END DO
  
      enrear1  = enrear1_d
      enrear   = enrear-enrear1*npartto
  END IF
  
  !     correct Coulomb potential by 1/N
  
  IF(numspin==2) THEN
    !ALLOCATE(coulsum(kdfull2))
    !ALLOCATE(couldif(kdfull2))
    !ALLOCATE(rho1(kdfull2))
    !ALLOCATE(rho2(kdfull2))
  
  
    call compcoulombadsic<<<ceiling(real(nxyz)/256),256>>>(rho1_d,rho2_d,rho_d)
    !DO ind=1,nxyz
    !  rho2(ind)=  rho(ind)
    !  rho1(ind)=  rho(ind)*rho(ind+nxyz)
    !END DO
  
    call reinitrhokr
    call curhofld<<<gridrhofld,tBlockrhofld>>>(rho1_d,rhokr_d,rhoki_d,nx2,ny2,nz2)
    !rhointest(:kdfull2) = rho1(:kdfull2)
    call cucoufou2()
    call curesult<<<gridresult,tBlockresult>>>(couldif_d,rhokr_d,e2,nx2,ny2,nz2)
    !CALL falr(rho1,couldif,kdfull2)
#if(raregas)
    IF (idielec == 0 .AND. nion2 > 0) THEN
#else
      IF (nion2 > 0) THEN
#endif
      !STOP "Nooooo inadsic !!!"
    coulsum_d = chpcoul_d
    ELSE
    !  CALL falr(rho2,coulsum,kdfull2)
    call reinitrhokr
    call curhofld<<<gridrhofld,tBlockrhofld>>>(rho2_d,rhokr_d,rhoki_d,nx2,ny2,nz2)
    !rhointest(:kdfull2) = rho2(:kdfull2)
    call cucoufou2()
    call curesult<<<gridresult,tBlockresult>>>(coulsum_d,rhokr_d,e2,nx2,ny2,nz2)
    !CALL falr(rho2,coulsum,kdfull2)
    END IF
  
  #if(findiff|numerov)
   STOP "Nooooo !!!"  
  #endif
    facup = 1D0/npartup
    facdw = 1D0/npartdw
     call alocfromcoulsic<<<ceiling(real(nxyz)/256),256>>>(aloc_d,coulsum_d,couldif_d,facup,facdw)
    !DO ind=1,nxyz
    !  aloc(ind) = aloc(ind) - 0.5D0*(coulsum(ind)+couldif(ind))*facup
    !END DO
    !DO ind=1,nxyz
    !    aloc(ind+nxyz) = aloc(ind+nxyz) - facdw*0.5D0*(coulsum(ind)-couldif(ind))
    !END DO
  
  ELSE
    IF (nion == 0) THEN
        call reinitrhokr
        call curhofld<<<gridrhofld,tBlockrhofld>>>(rho_d(:kdfull2),rhokr_d,rhoki_d,nx2,ny2,nz2)
        call cucoufou2()
        call curesult<<<gridresult,tBlockresult>>>(chpcoul_d,rhokr_d,e2,nx2,ny2,nz2)
        !CALL falr(rho(1:kdfull2),chpcoul,kdfull2)
    END IF
     
    fac = 1D0/npartto
  
    call alocfromcoulsicsp1<<<ceiling(real(nxyz)/256),256>>>(aloc_d,chpcoul_d,fac) 
    !DO ind=1,size
    !  aloc(ind) = aloc(ind) - fac*chpcoul(ind)
    !END DO
  END IF
  
   !rho = aloc_d
   !i=1
   !do while (i<=size(aloc,1))
   !!if (ABS(rho(i)-aloc(i))>1.d-3) then
   !	write(*,'(a,I12,2x,3F12.8)') 'error ',i,rho(i),aloc(i),ABS(rho(i)-aloc(i))
   !end if
   !i = i+1
   !end do
  
  
  end subroutine cucalc_adsic
  
  
  attributes(host) subroutine cufalr()
    USE cucoulsolv
    USE cutools
    USE params, only: e2
  implicit none
  
  
  
  call cucoufou2()
  
  
  end subroutine cufalr
  
  
  
  attributes(host) subroutine cucalclocal()
  use params
  use cutools
  use coulsolv
  use cucoulsolv
  implicit none
  
  real(DP) :: rhotest(2*kdfull2)
  real(DP) :: chpdft(2*kdfull2),localresult(2*kdfull2)
  real(DP) :: rhon(kdfull2)
  !real(DP) :: rhokr(kdred),rhoki(kdred),localrhokr(kdred)
  !integer  :: ind
  integer :: i,custatus
  !real(dp) :: enreartestcu
  
  #if(extended)
  IF (ifreezekspot == 1 .AND. tfs > 0D0) RETURN
  #endif

  rhotest = rho_d
  
  !     the lda part
  
  IF(ifsicp /= 5) THEN
    call calc_lda(rhotest,chpdft)
   ! write(*,*) "enrear", enrear
    CALL cucalc_lda<<<1,256,256*8>>>(rho_d,chpdft_d,dvol,enrear_d)
   ! enreartestcu = enrear_d
     enrear = enrear_d 
   ! write(*,*) "enreartestcu", enreartestcu
  ELSE
    STOP "ifsicp == 5 not implemented for cuda"
    !chpdft = 0D0
  END IF
  !custatus =  cudaDeviceSynchronize()
  localresult = chpdft_d
  !custatus =  cudaDeviceSynchronize()
  
  do i = 1, 2*kdfull2
  if (ABS(localresult(i)-chpdft(i))>1.d-8) then
  		write(*,'(a,I12,2x,3F12.8)') 'error chpdft !',i,localresult(i),chpdft(i),ABS(localresult(i)-chpdft(i))
  end if
  end do
  
  
  !      the netto charge density
  
  
  
    IF(nion2 == 0) THEN
      DO i=1,kdfull2
      rhon(i)=rhotest(i)-rhojel(i)   
      END DO           
      call cunettocharg<<<ceiling(real(kdfull2)/256),256>>>(rhon_d,rho_d,rhojel_d)
   ELSE IF(nion2 /= 0) THEN
      DO i=1,kdfull2
      rhon(i)=rhotest(i) 
      END DO              
      call cuaffectcharg<<<ceiling(real(kdfull2)/256),256>>>(rhon_d,rho_d)
   END IF
  
  !     the Coulombic part
  !     warning : counet inserts the esquar factor
  !custatus =  cudaDeviceSynchronize()
  localresult(:kdfull2) = rhon_d
  !custatus =  cudaDeviceSynchronize()
  do i = 1, kdfull2
  if (ABS(localresult(i)-rhon(i))/ABS(rhon(i))>1.d-8) then
  		write(*,'(a,I12,2x,3F12.8)') 'error rhon ',i,localresult(i),rhon(i),ABS(localresult(i)-rhon(i))
  end if
  end do
  
  CALL solv_poisson(rhon,chpcoul,kdfull2)
  

  IF (nion2 == 0) THEN
  call reinitrhokr
  call curhofld<<<gridrhofld,tBlockrhofld>>>(rhon_d,rhokr_d,rhoki_d,nx2,ny2,nz2)
  custatus =  cudaDeviceSynchronize()
  call cufalr()
  call curesult<<<gridresult,tBlockresult>>>(chpcoul_d,rhokr_d,e2,nx2,ny2,nz2)
  END IF

  !custatus =  cudaDeviceSynchronize()
  !localresult(:kdfull2) = chpcoul_d
  !custatus =  cudaDeviceSynchronize()
  !i=1
  !do while (i<=kdfull2)
  !if (ABS(localresult(i)-chpcoul(i))/ABS(chpcoul(i))>1.d-8) then
  !		write(*,'(a,I12,2x,3F12.8)') 'error chpcoul ',i,localresult(i),chpcoul(i),ABS(localresult(i)-chpcoul(i))
  !end if
  !i=i+1
  !end do
  chpcoul_d = chpcoul
  !STOP "Bazinga" 
  
  !     the laser part
  
  
  IF(tfs > 0D0 .AND. e0>0.d0) THEN
    !ALLOCATE(vlaser(kdfull2))
    !CALL laserp(vlaser,rho)
    !ALLOCATE(Vproj(kdfull2))
    !CALL projectp(Vproj)
     STOP "laser not implemented for cuda..."
  END IF
  
  
  !       the sum
  
  !DO ind=1,nxyz
    IF(nion2 /= 0) THEN
  !     add=chpcoul(ind)-potion(ind)
  call cugetalocsubstractpotion<<<ceiling(real(kdfull2)/256),256>>>(aloc_d,chpcoul_d,chpdft_d,potion_d)
    ELSE
  !     add=chpcoul(ind)
  call cugetaloc<<<ceiling(real(kdfull2)/256),256>>>(aloc_d,chpcoul_d,chpdft_d)
    END IF
  !  IF(tfs > 0D0) THEN
  !     add = add + vlaser(ind) + Vproj(ind)
  !  END IF
  
  !  aloc(ind)=chpdft(ind)+add
  !  aloc(ind+nxyz)=chpdft(ind+nxyz)+add
  !call cugetaloc<<<ceiling(real(kdfull2)/256),256>>>(aloc_d,chpcoul_d,chpdft_d)
  
  !END DO
  
  
  !IF(tfs > 0D0)  THEN
     !DEALLOCATE(vlaser)
     !DEALLOCATE(Vproj) 
     !STOP "laser not implemented for cuda..."
  !ENDIF
  
  
  !      optionally static dipole potential
  
  IF(tdipolxyz) THEN 
    !ind=0
    !DO jz=1,nz2
    !  addz = (jz-nz)*dz*dpolz
    !  DO jy=1,ny2
    !    addy = (jy-ny)*dy*dpoly
    !    DO jx=1,nx2
    !      addx = (jx-nx)*dx*dpolx
    !      add = addx+addy+addz
    !      ind = ind + 1
    !      aloc(ind) = aloc(ind) + add
    !      aloc(ind+nxyz) = aloc(ind+nxyz) + add
    !    END DO
    !  END DO
    !END DO
    STOP "static dipole potential not implemented for cuda..."
  END IF
  
  !usew1 = .false.
  !      usew4 = .false.
  !DEALLOCATE(rhon)
  !DEALLOCATE(chpdft)
  !DEALLOCATE(vlaser)
  
  
  !IF (izforcecorr == 1) CALL zeroforce(aloc,rho)
  
  
  end subroutine cucalclocal
  
  
  
  
  !-----tstep---------------------------------------------------------
  
  SUBROUTINE cutstep(it)
  
  USE params
  USE cukinetic
  USE kinetic
  USE cutools
  #if(twostsic)
  USE twost, ONLY:tnearest
  #endif
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)                      :: it
  COMPLEX(DP),DEVICE,DIMENSION(:,:),ALLOCATABLE :: q1_d,q2_d
  COMPLEX(DP),DIMENSION(:,:),ALLOCATABLE :: q0,actpsi !test
  COMPLEX(DP),DIMENSION(:,:),ALLOCATABLE :: q1,q2     !test
  REAL(DP) , DIMENSION(:), ALLOCATABLE :: aloc,rho    !test
  REAL(DP) , DIMENSION(:), ALLOCATABLE :: rhocu,aloccu!test
  LOGICAL :: tenerg
  REAL(DP) :: tin,tout
  INTEGER :: ind, ishift, ithr, itsub,  nb, nlocact, nup, i
  INTEGER :: ncount_init, ncount_rate, ncount_max, ncount_syst, ncount_fin
  REAL(DP) :: ri, dt, pr, time_init, time_fin, time_cpu
  
  #if(parayes)
  myn = 0
  #endif
  ALLOCATE(q1(2*kdfull2,0:nthr)) !!
  ALLOCATE(q2(kdfull2,0:nthr)) !!
  ALLOCATE(q0(kdfull2,nstate)) !! 
  ALLOCATE(actpsi(kdfull2,nstate)) !! 
  ALLOCATE(aloc(2*kdfull2),rho(2*kdfull2)) !! 
  allocate(rhocu(2*kdfull2))
  allocate(aloccu(2*kdfull2))

  q0   = psi_d !!
  aloc = aloc_d !!
  rho  = rho_d !!
  CALL cpu_time(time_init)
  CALL system_clock(ncount_init,ncount_rate,ncount_max)
  !IF (ntref > 0 .AND. it > ntref) nabsorb = 0           ! is that the correct place?
  
  !itsub = MOD(it,ipasinf) + 1
  
  ri = -dt1*0.5D0
  dt = dt1*0.5D0
  nlocact = numspin*nxyz
  !CALL tstep(q0,aloc,rho,it)
  
  !     half time step in coordinate space
  !     local phase field on workspace 'psi_d'
  aloc(:) = 0.d0
  !write(*,*) "call cucoordstep"
  IF (numspin==2) THEN
    call cucoordstep<<<kdfull2/256,256>>>(psi_d,aloc_d,dt)
  ELSE
    call cucoordstepsp1<<<kdfull2/256,256>>>(psi_d,aloc_d,dt)
  ENDIF
  
  DO ind=1,nlocact
     pr=-dt*aloc(ind)
     q1(ind,0)=CMPLX(COS(pr),SIN(pr),DP)
  END DO
  DO nb=1,nstate
    ishift = (ispin(nrel2abs(nb))-1)*nxyz
    q0(:,nb) = q1(ishift+1:ishift+kdfull2,0)*q0(:,nb)
  END DO
  
  !     half non-local step
  
  IF(ipsptyp == 1 .AND. tnonlocany) THEN
  !CALL cpu_time(tin)
      allocate(q1_d(kdfull2,nstate),q2_d(kdfull2,nstate))
      tenerg = MOD(it, jinfo) == 0 .OR. MOD(it, jstinf) == 0 .OR. &
      MOD(it, jenergy) == 0
      CALL cunonlocstep(psi_d,q1_d,q2_d,dt,tenerg,6)   ! 4
  !CALL cpu_time(tout)
      !write(*,*) "time for cunonlocstep",tout-tin
      !DO nb=1,nstate
        !  tenerg = itsub == ipasinf
        !  CALL nonlocstep(q0(1,nb),q1,q2,dt,tenerg,nb,6)   ! 4
      !END DO
  END IF
  
  
  !       one full time step for the kinetic energy
  
  ithr=0
  !DO nb=1,nstate
  !CALL cpu_time(tin)
  !CALL cukinprop()
  DO nb=1,nstate
        CALL kinprop(q0(1,nb))
  END DO
      !write(*,*) "call cukinprop"
  !CALL cpu_time(tout)
   !   write(*,*) "time step for the kinetic energy",tout-tin

  
  
  psi_d = q0
  #if(findiff|numerov)
      STOP "No for cuda cutstep"
  #endif
  !END DO
  
  CALL flush(7)
  
  !     half non-local step
  
  IF(ipsptyp == 1 .AND. tnonlocany) THEN
  !CALL cpu_time(tin)
        tenerg = MOD(it, jinfo) == 0 .OR. MOD(it, jstinf) == 0 .OR. &
      MOD(it, jenergy) == 0
    CALL cunonlocstep(psi_d,q1_d,q2_d,dt,tenerg,6)
    deallocate(q1_d,q2_d)
  !CALL cpu_time(tout)
     ! write(*,*) "time for cunonlocstep",tout-tin
      !DO nb=1,nstate
        !  tenerg = itsub == ipasinf
        !  CALL nonlocstep(q0(1,nb),q1,q2,dt,tenerg,nb,6)   ! 4
      !END DO
  END IF
  
  
  CALL dyn_mfield(rho,aloc,q0,dt,it)
  !write(*,*) "call cudyn_mfield"
  !CALL cpu_time(tin)
  CALL cudyn_mfield()


  rhocu = rho_d
  aloccu = aloc_d
  !custat  = cudaDeviceSynchronize()
  !i=1
  !do while (i<=size(rho,1))
  !if (ABS(rhocu(i)-rho(i))>1.d-8) then
  !	write(*,'(a,I12,2x,3F12.8)') 'error rho ',i,rhocu(i),rho(i),ABS(rhocu(i)-rho(i))
  !end if
  !i = i+1
  !end do
  i=1
  do while (i<=size(aloc,1))
  if (ABS(aloccu(i)-aloc(i))>1.d-8) then
    write(*,'(a,I12,2x,3F12.8)') 'error ',i,aloccu(i),aloc(i),ABS(aloccu(i)-aloc(i))
  end if
  i = i+1
  end do
  !STOP "Bazinga !"










  !CALL cpu_time(tout)
   !   write(*,*) "time for cudyn_mfield",tout-tin
  !     half time step in coordinate space:
  
  !write(*,*) "call cucoordstep"
  IF (numspin==2) THEN
    call cucoordstep<<<kdfull2/256,256>>>(psi_d,aloc_d,dt)
  ELSE
    call cucoordstepsp1<<<kdfull2/256,256>>>(psi_d,aloc_d,dt)
  ENDIF
  
  
          !actpsi = psi_d
    !DO ind = 1,kdfull2
    !DO nb = 1 , nstate
    !IF (ABS(actpsi(ind,nb)-q0(ind,nb))/ABS(q0(ind,nb)) >= 1.d-8) &
    !write(*,'(a,3E18.5)') "error"  ,ABS(actpsi(ind,nb)),ABS(q0(ind,nb)),ABS(actpsi(ind,nb)-q0(ind,nb))
    !END DO
    !END DO
    !STOP
  CALL cpu_time(time_fin)
  time_cpu = time_fin-time_init
  CALL system_clock(ncount_fin,ncount_rate,ncount_max)
  ncount_syst=ncount_fin-ncount_init
  IF(myn == 0)THEN
    WRITE(6,'(a,2(1pg13.5))') ' CPU time in TSTEP',time_cpu,ncount_syst*1D-4
    WRITE(7,'(a,2(1pg13.5))') ' CPU time in TSTEP',time_cpu,ncount_syst*1D-4
    CALL FLUSH(6)
    CALL FLUSH(7)
  END IF
  
  #if(extended)
  IF (izforcecorr /= -1) THEN
    STOP "izforcecorr for CUDA"
  END IF

  IF ((jescmask > 0 .AND. MOD(it,jescmask) == 0) .OR. &
      (jescmaskorb > 0 .AND. MOD(it,jescmaskorb) == 0)  ) STOP "no mask for cuda"!CALL  escmask(it)
  #endif

  CALL flush(6)
  CALL flush(7)
  
  end subroutine cutstep
  
  attributes(host) subroutine cunonlocstep(inpsi_d,q1_d,q2_d,ri,tenerg,norder)
  USE params
  USE cutools
  IMPLICIT NONE
  COMPLEX(DP) , DEVICE , INTENT(IN OUT)  :: inpsi_d(kdfull2,nstate)
  COMPLEX(DP) , DEVICE , INTENT(IN OUT)  :: q1_d(kdfull2,nstate)
  COMPLEX(DP) , DEVICE , INTENT(IN OUT)  :: q2_d(kdfull2,nstate)
  REAL(DP), INTENT(IN)                         :: ri
  LOGICAL, INTENT(IN)                      :: tenerg
  INTEGER, INTENT(IN)                      :: norder
  
  REAL(DP) ::  sumadd
  COMPLEX(DP) :: rii,cfac
  
  !write(*,*) "cutotnonlocalc"; call flush(6);
  CALL cutotnonlocalc(inpsi_d,q1_d,0)
  IF(tenerg) THEN !  add nonloc.pot energy
    !!sumadd = 0D0
    !write(*,*) "call cugetnonlocenergies"; call flush(6)
    !call cugetnonlocenergies<<<nstate,256,256*8>>>(inpsi_d,q1_d,enonloc_d)
    !!DO  i=1,nxyz
    !!  sumadd  = REAL(qact(i),DP)*REAL(q1(i),DP) +AIMAG(qact(i))*AIMAG(q1(i))  + sumadd
    !!END DO
     !write(*,*) "sizes ", size(enonlo,1) , size(enonloc_d,1)
     !enonlo(:nstate) = enonloc_d(:nstate)
     !epotsp(:nstate) = enonlo(:nstate) + epotsp(:nstate)
    
  END IF
  
   cfac=-ri*eye
  rii=cfac
  !write(*,*) "addenscontribtaylor"; call flush(6);
  call addenscontribtaylor<<<kdfull2/256,256>>>(inpsi_d,q1_d,rii)
  
  !write(*,*) "cutotnonlocalc"; call flush(6);
  CALL cutotnonlocalc(q1_d,q2_d,0)
  rii=rii*cfac/2D0
  call addenscontribtaylor<<<kdfull2/256,256>>>(inpsi_d,q2_d,rii)
  IF(norder <= 2) RETURN
  
  CALL cutotnonlocalc(q2_d,q1_d,0)
  rii=rii*cfac/3D0
  call addenscontribtaylor<<<kdfull2/256,256>>>(inpsi_d,q1_d,rii)
  IF(norder <= 3) RETURN
  
  CALL cutotnonlocalc(q1_d,q2_d,0)
  rii=rii*cfac/4D0
  call addenscontribtaylor<<<kdfull2/256,256>>>(inpsi_d,q2_d,rii)
  IF(norder <= 4) RETURN
  
  CALL cutotnonlocalc(q2_d,q1_d,0)
  rii=rii*cfac/5D0
  call addenscontribtaylor<<<kdfull2/256,256>>>(inpsi_d,q1_d,rii)
  IF(norder <= 5) RETURN
  
  CALL cutotnonlocalc(q1_d,q2_d,0)
  rii=rii*cfac/6D0
  call addenscontribtaylor<<<kdfull2/256,256>>>(inpsi_d,q2_d,rii)
  
  
  end subroutine cunonlocstep
  
  attributes(host) subroutine cutstep_exp(it,timagtime)
  
  USE params
  USE cukinetic
  USE cutools
  USE util, ONLY:wfovlp
  IMPLICIT NONE
  INTEGER, INTENT(IN)                      :: it
  LOGICAL, INTENT(IN)			 :: timagtime
  INTEGER(DP) :: nb,i
  COMPLEX(DP) :: cdtact
  LOGICAL,PARAMETER :: tnorotate=.true.
  !COMPLEX(DP) :: qwork(kdfull2,nstate), test(kdfull2,nstate)
  !COMPLEX(DP) :: q1(kdfull2)
  !REAL(DP)    :: aloc(2*kdfull2),aloctest(2*kdfull2),aloccpu(2*kdfull2),rhocpu(2*kdfull2)
  !---------------------------------------
  
  IF(ifsicp==5) STOP "ifsicp=5 not for cuda ..."
  
  !aloc = aloc_d
  kpsi_d	   = psi_d
  
  IF(.NOT.timagtime) THEN
    cdtact = CMPLX(dt1/2D0,0D0,DP)
    IF(tnorotate .OR. ifsicp .NE. 8) THEN
      DO nb=1,nstate
         !qwork(:,nb) = psi_d(:,nb)!q0(:,nb)
        !CALL exp_evol(qwork(:,nb),aloc,nb,4,cdtact,q1)
         CALL cuexp_evol(kpsi_d,4,cdtact,nb)
      END DO
    ELSE
  #if(twostsic)
       STOP "twostsic not for cuda ..."
      !qwork = q0
      !CALL exp_evolp(qwork,aloc,4,cdtact,q1,q0)
  #else
      STOP " IFSICP==8 requires compilation with option twostsic"
  #endif
    END IF
  
  #if(twostsic)
     STOP "twostsic not for cuda ..."
    !IF(tnearest .AND. ifsicp==8) CALL eval_unitrot(qwork,q0)
  #endif
    !test = kpsi_d
     !DO nb = 1 ,nstate 
    !DO i = 1,kdfull2
    !if(abs(test(i,nb)-qwork(i,nb))>=1.d-4) then
    !write(*,*) "error",abs(test(i,nb)-qwork(i,nb))
    !end if
    !END DO
    !END DO
  
     !CALL dyn_mfield(rhocpu,aloccpu,test,dt1*0.5D0,it)
     CALL cudyn_mfield(kpsi_d)
       !aloctest = aloc_d
    !DO i = 1,kdfull2
    !if(abs(aloctest(i)-aloccpu(i))>=1.d-4) then
    !write(*,*) "error aloc",abs(aloctest(i)-aloccpu(i))
    !end if
    !END DO
  
       !aloctest = rho_d
    !DO i = 1,kdfull2
    !if(abs(aloctest(i)-rhocpu(i))>=1.d-4) then
    !write(*,*) "error rho",abs(aloctest(i)-rhocpu(i))
    !end if
    !END DO
  
  END IF
  
  !     full time step to next wavefunctions
  !     use exponential evolution to fourth order
  
  IF(timagtime) THEN
    cdtact = CMPLX(0D0,-dt1,DP)
  ELSE
    cdtact = CMPLX(dt1,0D0,DP)
  END IF
  
  !nterms = 4
  IF(tnorotate .OR. ifsicp .NE. 8) THEN
    DO nb=1,nstate
      !CALL exp_evol(q0(:,nb),aloc,nb,nterms,cdtact,q1)
       CALL cuexp_evol(psi_d,4,cdtact,nb)
    END DO
  ELSE
    STOP "how to do expevol ? (cuda) ..."
  #if(twostsic)
    !qwork = q0
    !CALL exp_evolp(q0,aloc,nterms,cdtact,q1,qwork)
  #endif
  END IF
  
  #if(twostsic)
  !IF(tnearest .AND. ifsicp==8 .AND. .NOT.timagtime) CALL eval_unitrot(q0,qwork)
  STOP "twostsic not for cuda ..."
  #endif
  #if(extended)
  IF (ntref > 0 .AND. it > ntref) nabsorb = 0
  #endif
  !     compute mean field at new time
  
  !CALL dyn_mfield(rho,aloc,q0,dt1,it)
  CALL cudyn_mfield()
  return
  end subroutine cutstep_exp
  
  attributes(host) subroutine cuexp_evol(actpsi_d,norder,dtact,nbe)
  use params, only:ispin,nrel2abs,nxyz
  use cutools
  IMPLICIT NONE
  !  ---------------------------------------
  COMPLEX(DP), DEVICE, INTENT(IN OUT)      :: actpsi_d(kdfull2,nstate)
  INTEGER, INTENT(IN)                      :: norder
  COMPLEX(DP), INTENT(IN)                  :: dtact
  INTEGER(DP), INTENT(IN)			 :: nbe
  !COMPLEX(DP) :: qwork(kdfull2),test(kdfull2),qact(kdfull2) !test
  !REAL(DP) :: aloc(2*kdfull2) !test
  INTEGER :: i,nb !test
  !LOGICAL :: iferror !test
  !  ---------------------------------------
  INTEGER :: ilocbas, isig, nterm
  COMPLEX(DP) :: dti,cfac
  !  -------------------
  
  IF (ispin(nrel2abs(nbe)) == 1) THEN
    ilocbas = 1
  ELSE IF (ispin(nrel2abs(nbe)) == 2) THEN
    ilocbas = nxyz+1
  ELSE
    STOP " EXPEVOL: spin index must be 1 or 2"
  END IF
  IF(ABS(IMAG(dtact))>1D-10) THEN
    isig = -1
  ELSE
    isig = 1
  END IF
  !isig=-1
  dti = dtact*CMPLX(0D0,1D0,DP)
   cfac = CMPLX(1D0,0D0,DP)
  !DO  i=1,nxyz
  !  qwork(i) = qact(i)
  !END DO
  call affectpsiwork<<<kdfull2/256,256>>>(actpsi_d,psiwork_d,nbe)
  !call affectpsiwork<<<kdfull2/MAX_THREADS_PER_BLOCKS,MAX_THREADS_PER_BLOCKS>>>(actpsi_d,psiwork_d,nbe)
  !aloc = aloc_d
  !qact = psiwork_d
  !qwork(:) = qact(:)
  DO nterm=1,norder
    !CALL hpsi(qwork,aloc(ilocbas),nbe,isig*nterm)
    call cuhpsi(nbe,isig*nterm)
    cfac = -dti/nterm*cfac
      !DO  i=1,nxyz
      !qact(i) = qact(i) + cfac*qwork(i)
      !END DO
     call addcontribtaylor<<<kdfull2/256,256>>>(actpsi_d,psiwork_d,nbe,cfac)
  
  END DO
          !iferror = .FALSE.
    !test = actpsi_d(:,nbe)
    !DO i = 1,kdfull2
    !if(abs(test(i)-qact(i))>=1.d-8) then
    !write(*,*) "error",abs(test(i)-qact(i))
          !iferror = .TRUE.
    !end if
    !END DO
         ! IF(iferror) write(*,*) "State :",nbe,"  ERROR"
         ! IF(.NOT.iferror) write(*,*) "State :",nbe," no error"
  RETURN
  
  end subroutine cuexp_evol
  
  subroutine cuhpsi(nbe,itpri)
  USE params
  USE cutools
  USE util, ONLY:wfovlp
  USE cukinetic
  implicit none
  !  ---------------------------------------
  INTEGER(DP), INTENT(IN)                      :: nbe
  INTEGER		,INTENT(IN)		     :: itpri
  REAL(DP) :: tnorm
  LOGICAL :: tpri
  COMPLEX(DP) , DEVICE, ALLOCATABLE  :: q1_d(:),q2_d(:)
  !COMPLEX(DP) ,  ALLOCATABLE  :: q1(:),work(:),q2(:)
  real(DP) , DEVICE  :: spvariance_d,spamoy_d
  real(DP) :: locspvar, locenergy
  LOGICAL,PARAMETER :: ttest=.FALSE.
  INTEGER :: is
  INTEGER :: shft
  !  ---------------------------------------
  allocate(q1_d(kdfull2),q2_d(kdfull2))
  !allocate(q1(kdfull2),q2(kdfull2),work(kdfull2))
  shft = ispin(nrel2abs(nbe))-1
  !  -------------------
  tpri =  ABS(itpri)==1
  tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(kxbox*kybox*kybox,DP))
  call cucopy1dto3dworking<<<gridfft3d,tBlockfft3d>>>(psiwork_d,ffta_d)
  call cufftExecZ2Z(plan3d, ffta_d(:,:,:), ffta_d(:,:,:), CUFFT_FORWARD)
  call cucopy3dto1dworking<<<gridfft3d,tBlockfft3d>>>(q1_d,ffta_d,tnorm) 
  call multpsiwork<<<kdfull2/256,256>>>(akv_d,q1_d)
  tnorm =SQRT(8D0*pi*pi*pi)/SQRT(REAL(kxbox*kybox*kybox,DP))
  call cucopy1dto3dworking<<<gridfft3d,tBlockfft3d>>>(q1_d,ffta_d)
  call cufftExecZ2Z(plan3d, ffta_d(:,:,:), ffta_d(:,:,:), CUFFT_INVERSE)
  call cucopy3dto1dworking<<<gridfft3d,tBlockfft3d>>>(q2_d,ffta_d,tnorm) 
  !  -------------------
  IF(ipsptyp == 1) THEN
  !work = psiwork_d
  call cunonlocalc(q1_d,0)
  !CALL nonlocalc(work,q1,0)
  IF(tpri) then 
    call cuwfovlp<<<kdfull2/256,256,256*16>>>(psiwork_d,q1_d,enonloc_d(nbe))
    !q2 = q1_d
    !do is = 1,kdfull2; if(abs(q2(is))>1.d-1) write(*,*) q2(is); enddo
    !STOP
    !do is = 1,kdfull2
  !	if (abs(q2(is)-q1(is))>=1.d-6) write(*,*) abs(q2(is)),abs(q1(is))
   ! end do 
    enonlo(nbe) = enonloc_d(nbe)
    !write(*,*) "nbe ",nbe,"enonlo ",enonlo(nbe)
    !enonlo(nbe) = wfovlp(work,q1)
    !write(*,*) "nbe ",nbe,"enonlo ",enonlo(nbe)
    !if (nbe==20) STOP "end"
  ENDIF
  call addmultpsiwork<<<kdfull2/256,256>>>(aloc_d(shft*nxyz+1:),q1_d,psiwork_d)
  ELSE
  call affectmultpsiwork<<<kdfull2/256,256>>>(aloc_d(shft*nxyz+1:),q1_d,psiwork_d)
  END IF
  
  IF(ifsicp==5) THEN
    STOP "ipsptyp == 5 not implemented for cuda"
  END IF
  
  #if(twostsic)
    STOP "twostsic not implemented for cuda"
  #endif
  IF(tpri) THEN
    call cuwfovlp<<<kdfull2/256,256,256*16>>>(psiwork_d,q1_d,epot_d(nbe))
    epotsp(nbe) = epot_d(nbe)
    amoy(nbe) = ekinsp(nbe)+epotsp(nbe)
    !q2 = q1+q2
    call addpsiwork<<<kdfull2/256,256>>>(q2_d,q1_d)
    call cuwfovlp<<<kdfull2/256,256,256*16>>>(q2_d,q2_d,spvariance_d)
    call cuwfovlp<<<kdfull2/256,256,256*16>>>(psiwork_d,q2_d,spamoy_d)
    locspvar  = spvariance_d
    locenergy = spamoy_d
    spvariance(nbe) = SQRT(MAX(locspvar-locenergy**2,1D-99))
    is=ispin(nrel2abs(nbe))
    IF(ttest) WRITE(*,'(a,2i4,5(1pg13.5))') &
     ' HPSI: nbe,is,esp,var=',nbe,is,amoy(nbe),spvariance(nbe), &
        ekinsp(nbe),epotsp(nbe),locspvar
    CALL flush(6)
  ELSE
    call addpsiwork<<<kdfull2/256,256>>>(q2_d,q1_d)
  END IF
  
  IF(itpri<0) THEN
    call addvaluedpsiwork<<<kdfull2/256,256>>>(-amoy(nbe),psiwork_d,q2_d)
  ELSE
    psiwork_d = q2_d
  END IF
  
  deallocate(q1_d,q2_d)
  !deallocate(q1,work,q2)
  end subroutine cuhpsi
  
  
  subroutine cutotnonlocalc(inpsi_d,aux_d,ionact)
  use params
  use cutools
  IMPLICIT NONE
  
  COMPLEX(DP), device, INTENT(IN) :: inpsi_d(kdfull2,nstate)
  COMPLEX(DP), device, INTENT(IN OUT) :: aux_d(kdfull2,nstate)
  INTEGER , INTENT(IN)	:: ionact
  INTEGER :: i, ion!,ind
  INTEGER :: maxion, minion
  !LOGICAL :: error
  !COMPLEX(DP) :: psi(kdfull2),aux(kdfull2),actsum
  !COMPLEX(DP) ,DEVICE :: actsum_d(nstate)
  
  !DO i=1,kdfull2
  !  aux(i)=CMPLX(0D0,0D0,DP)
  !END DO
  !psi = inpsi_d(:,1)
  !aux = (0.d0,0.d0)
  !write(*,*) "reset2d"
  call reset2d<<<dim3(kdfull2/256,nstate,1),dim3(256,1,1)>>>(aux_d)
  IF(ionact == 0) THEN
    minion = 1
    maxion = nion
  ELSE IF(ionact <= nion .AND. ionact > 0) THEN
    minion = ionact
    maxion = ionact
  ELSE
    STOP ' NONLOCALC: ion out of range'
  END IF
  DO ion=minion,maxion
    ! write(*,*) "launch cutoteffnl4c" ,ifin(ion),nstate
    IF(tnonloc(ion)) call cutoteffnl4c<<<nstate,128,131*16>>>(inpsi_d,aux_d,ion) !CALL effnl4c(psi,aux,ion)
  !DO i = 1 , nstate
    !write(*,*) ""
    !psi = inpsi_d(:,i)
    !aux = (0.d0,0.d0)
    !IF(tnonloc(ion)) call effnl4c(psi,aux,ion)
    !psi = aux_d(:,i)
    !error = .false.
    !DO ind = 1,kdfull2
    !IF (ABS(aux(ind)-psi(ind))/ABS(aux(ind)) >= 1.d-8) then
    !error = .true.
    !write(*,'(a,3ES14.7)') "error"  ,ABS(aux(ind)),ABS(psi(ind)),ABS(aux(ind)-psi(ind))
    !endif
    !END DO
          !actsum = actsum_d(i)
    !write(*,*) "sumr01_d ",actsum
    !write(*,*) "nb ",i,"error ",error
  !ENDDO 
  END DO
  RETURN
  
  
  end subroutine cutotnonlocalc
  
  subroutine cunonlocalc(aux_d,ionact)
  use params
  use cutools
  IMPLICIT NONE
  
  COMPLEX(DP), device, INTENT(IN OUT) :: aux_d(kdfull2)
  INTEGER , INTENT(IN)	:: ionact
  INTEGER :: i, ion
  INTEGER :: maxion, minion
  
  !DO i=1,kdfull2
  !  aux(i)=CMPLX(0D0,0D0,DP)
  !END DO
  call reset1d<<<kdfull2/256,256>>>(aux_d)
  IF(ionact == 0) THEN
    minion = 1
    maxion = nion
  ELSE IF(ionact <= nion .AND. ionact > 0) THEN
    minion = ionact
    maxion = ionact
  ELSE
    STOP ' NONLOCALC: ion out of range'
  END IF
  DO ion=minion,maxion
     !write(*,*) "launch cueffnl4c" ,ifin(ion)
    IF(tnonloc(ion)) call cueffnl4c<<<1,256,256*16>>>(psiwork_d,aux_d,ion) !CALL effnl4c(psi,aux,ion)
  END DO
  RETURN
  
  
  end subroutine cunonlocalc
  
  
  
  subroutine cukinprop()
  !       propagation with exp(-i*dt*e_kin)
  
  USE params
  USE cutools
  USE FFTW
  USE cukinetic
  IMPLICIT NONE
  !#if(fftw_cpu && !oldkinprop)
  !INTEGER :: ithr
  !#endif
  REAL(DP):: facnr
  !COMPLEX(DP) ,allocatable :: q1(:,:),actffta(:,:,:,:),fftacpu(:,:,:,:)
  !integer xx,yy,zz
  integer :: ii
  !allocate(q1(kdfull2,nstate),actffta(nx2,ny2,nz2,nstate),fftacpu(nx2,ny2,nz2,nstate))
  
  !q1 = psi_d
  
  
  !tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))
  !  here version using 3D FFTW
  #if(fftw_cpu && !oldkinprop)
  !  ithr = 0
  facnr = 1D0/(nx2*ny2*nz2)
  DO ii = 1 ,nstate
  !write(*,*) "call cucopy1dto3d"; call flush(6)
  call cucopy1dto3d<<<gridfft3d,tBlockfft3d>>>(psi_d,ffta_d,ii)
  !DO ii = 1 ,nstate
  !CALL copy1dto3d(q1(:,ii),fftacpu(:,:,:,ii),nx2,ny2,nz2)
  !END DO
  
  !fftacpu = ffta_d
  !CALL fftw_execute_dft(pforw(0),fftacpu(:,:,:,ii),fftacpu(:,:,:,ii))
  call cufftExecZ2Z(plan3d, ffta_d(:,:,:), ffta_d(:,:,:), CUFFT_FORWARD)
  !write(*,*) "call cufftExecZ2Z"; call flush(6)
  !call cufftExecZ2Z(planmany, ffta_d, ffta_d, CUFFT_FORWARD)
  !ffta_d = fftacpu
  !DO ii = 1 ,nstate
  !fftacpu(:,:,:,ii) = akprop*fftacpu(:,:,:,ii)
  !END DO
  !write(*,*) "call cumultakprop"; call flush(6)
  call cumultakprop<<<gridfft3d,tBlockfft3d>>>(ffta_d,akprop_d)
  
  !fftacpu = ffta_d
  
  !CALL fftw_execute_dft(pback(0),fftacpu(:,:,:,ii),fftacpu(:,:,:,ii))
  call cufftExecZ2Z(plan3d, ffta_d(:,:,:), ffta_d(:,:,:), CUFFT_INVERSE)
  
  !write(*,*) "call cufftExecZ2Z"; call flush(6)
  !call cufftExecZ2Z(planmany, ffta_d, ffta_d, CUFFT_INVERSE)
  
  !ffta_d = fftacpu
  !DO ii = 1 ,nstate
  !CALL secopy3dto1d(fftacpu(:,:,:,ii),q1(:,ii),facnr,nx2,ny2,nz2)
  !END DO
  !psi_d = q1
  !write(*,*) "call cucopy3dto1d"; call flush(6)
  call cucopy3dto1d<<<gridfft3d,tBlockfft3d>>>(psi_d,ffta_d,facnr,ii) 
  END DO
  !actffta= ffta_d
  !DO xx = 1,nx2
  !DO yy = 1,ny2
  !DO zz = 1,nz2
  !DO ii = 1,nstate
  ! if (abs(actffta(xx,yy,zz,ii)-fftacpu(xx,yy,zz,ii))/abs(fftacpu(xx,yy,zz,ii))>1.d-3 &
  !	.AND. abs(fftacpu(xx,yy,zz,ii)) >= 1.d-8) &
  !	write(*,*) "error",abs(actffta(xx,yy,zz,ii)),abs(fftacpu(xx,yy,zz,ii))
  !END DO
  !END DO
  !END DO
  !END DO
  !STOP
  
  
  
  #else
   STOP "A such kinprop for CUDA is not implemented"
  #endif
  
  end subroutine cukinprop
  
  END MODULE cudadyn
  