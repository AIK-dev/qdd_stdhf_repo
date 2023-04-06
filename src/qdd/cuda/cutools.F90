
MODULE cutools
USE params, only: DP
USE cuparams
USE cudafor
CONTAINS

attributes(global) subroutine cucalcrhospindep(inpsi,occ,outrho) 
  use params
	implicit none
	complex(DP)  :: inpsi(:,:)
	real(DP)     :: occ(:)
	real(DP)     :: outrho(:)

	integer  :: i,j

	 i = (blockIdx%x-1)*blockDim%x + threadIdx%x ! MAX kdfull2
	 !j = (blockIdx%y-1)*blockDim%y + threadIdx%y ! MAX nstate/2

	if (i<=size(inpsi,1)) then
		outrho(i)		= 0.d0
		outrho(i+size(inpsi,1)) = 0.d0
	    j=1
	    do while (j<= size(occ,1)/2) !warning numsp==2 ?
	    outrho(i+size(inpsi,1))	= outrho(i+size(inpsi,1)) +&
	      occ(j+size(occ,1)/2)*CONJG(inpsi(i,j+size(occ,1)/2))*inpsi(i,j+size(occ,1)/2)
	    outrho(i)			= outrho(i)+ occ(j)*CONJG(inpsi(i,j))*inpsi(i,j)
	    j =j+1
	    end do
	    outrho(i)			= outrho(i)+outrho(i+size(outrho,1)/2)
	    outrho(i+size(outrho,1)/2) 	= outrho(i)-2.0*outrho(i+size(outrho,1)/2)
	    outrho(i+size(outrho,1)/2)	= outrho(i+size(outrho,1)/2)/MAX(outrho(i),1.D-8)
	end if

end subroutine cucalcrhospindep

attributes(global) subroutine cucalcrhotot(inpsi,occ,outrho) 
  use params
	implicit none
	complex(DP)  :: inpsi(:,:)
	real(DP)     :: occ(:)
	real(DP)     :: outrho(:,:)

	integer  :: i,j,spin

	 i = (blockIdx%x-1)*blockDim%x + threadIdx%x ! MAX kdfull2
	 !j = (blockIdx%y-1)*blockDim%y + threadIdx%y ! MAX nstate/2

	if (i<=size(inpsi,1)) then
		outrho(i,:)		= 0.d0
	    j=1
	    do while (j<= size(occ,1))
	    if (j<=size(occ,1)/2) then
		spin=1
	    else
		spin=2
	    endif
	    outrho(i,spin)	= outrho(i,spin)+ occ(j)*CONJG(inpsi(i,j))*inpsi(i,j)
	    j =j+1
	    end do
	end if

end subroutine cucalcrhotot

attributes(global) subroutine cukerneletdhf(inpsi,inoccup,spenrg,&
						actoccup,avtrans,&
						itrans,xtrans,&
						jjumpaverwidth,jjumpVtrans,strength,tau)
	implicit none
	complex(DP)  	:: inpsi(:,:)
	real(DP)	:: inoccup(:)
	real(DP)	:: spenrg(:)
	real(DP)	:: actoccup(:,:)
	real(DP)	:: avtrans(:)
	integer		:: itrans(:)
	real(DP)	:: xtrans(:)
	real(DP), value	:: jjumpaverwidth,jjumpVtrans,strength,tau

	integer :: i,nst,ind
	integer :: i1,i2,i3,i4
	real(DP):: edif,xdum
	complex(DP) :: trprb

	i = blockDim%y*(threadIdx%y-1)+threadIdx%x
	i1 = threadIdx%x
	i2 = threadIdx%y
	nst = size(inoccup,1)
    do while (i1 <= nst)
	i2 = threadIdx%y
        do while (i2 <= nst)
            if (i2/=i1) then
                do i3 = 1, nst
                    if (i3/=i1 .and. i3/=i2) then
                        do i4 = 1, nst
                            if (i4/=i1 .and. i4/=i2 .and. i4/=i3 .and. cugettotspin(i1,i2,i3,i4)==0) then

                                edif = spenrg(i1)+spenrg(i2) - spenrg(i3)-spenrg(i4)
                                if (abs(edif)<jjumpaverwidth) then
                                    trprb=cmplx(0.d0,0.d0)
                                    do ind = 1, size(inpsi,1)
                                        trprb=trprb+ &
                                            conjg(inpsi(ind,i1)*inpsi(ind,i2))*inpsi(ind,i3)*inpsi(ind,i4)
                                    end do
                                    trprb=trprb*dvol_d*jjumpVtrans
                                    xdum=strength*REAL(conjg(trprb)*trprb,8)*( &
                                        (inoccup(i3)*inoccup(i4)*(1-inoccup(i1))*(1-inoccup(i2))) &
                                        -(inoccup(i1)*inoccup(i2)*(1-inoccup(i3))*(1-inoccup(i4))))
                                    actoccup(i,i1)=actoccup(i,i1)+xdum
                                    if (xdum>1.d-8) itrans(i) = itrans(i) + 1
                                    avtrans(i) = avtrans(i) + REAL(CONJG(trprb)*trprb,8)
                                    xtrans(i)=xtrans(i)+(strength/tau)*REAL(conjg(trprb)*trprb,8)*( &
                                        (inoccup(i3)*inoccup(i4)*(1-inoccup(i1))*(1-inoccup(i2))) &
                                        +(inoccup(i1)*inoccup(i2)*(1-inoccup(i3))*(1-inoccup(i4))))
                                end if
                            end if
                        end do
                    end if
                end do
            end if
	    i2 = i2 + blockDim%y
        end do
	i1 = i1 + blockDim%x
    end do
	
	call syncthreads()
	ind = blockDim%x*blockDim%y/2

	do while (ind>= 1)
		if (i <= ind) &
		actoccup(i,:) = actoccup(i,:) + actoccup(i+ind,:)
		avtrans(i)    = avtrans(i)    + avtrans(i+ind)
		itrans(i)     = itrans(i)     + itrans(i+ind)
		xtrans(i)     = xtrans(i)     + xtrans(i+ind)
	ind = ind/2
	call syncthreads()
	end do

end subroutine cukerneletdhf

attributes(device) integer function cugettotspin(a,b,c,d)
	implicit none
	integer	:: a,b,c,d
	cugettotspin= ispin_d(a)+ispin_d(b)-ispin_d(c)-ispin_d(d)
end  function cugettotspin
attributes(global) subroutine cugradkspacestep(psik,akv,Apsi,epswf,e0dmp)
	implicit none
	complex(DP)  	:: psik(:)
	real(DP)  	:: akv(:)
	complex(DP)  	:: Apsi(:)
	real(DP),value	:: epswf,e0dmp

	integer :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x

	if (i<=size(psik,1)) then
	psik(i)=psik(i)-epswf/(e0dmp+akv(i))*Apsi(i)
	end if

end subroutine cugradkspacestep

attributes(global) subroutine caddALpsi(outALpsi,invpsi,injgradpsi,inhbarm)
	implicit none
	complex(DP)  :: outALpsi(:)
	complex(DP)  :: invpsi(:)
	complex(DP)  :: injgradpsi(:)
	real(DP),value:: inhbarm

	integer :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x

	if (i<=size(outALpsi,1)) then
	!ALpsi=vpsi+ALpsi-hbarm*jgradpsi
	outALpsi(i) = outALpsi(i) + invpsi(i) - inhbarm*injgradpsi(i)
	end if
end subroutine caddALpsi 

attributes(global) subroutine cucompjgradpsi(outjgradpsi,ingradpsi,ingradApsi,inlambdaj,curr0,curr1,muj)
	implicit none
	complex(DP)  :: outjgradpsi(:)
	complex(DP)  :: ingradpsi(:,:)
	complex(DP)  :: ingradApsi(:,:)
	real(DP)     :: inlambdaj(:,:)
	real(DP)     :: curr0(:,:)
	real(DP)     :: curr1(:,:)
	real(DP),value:: muj

	integer :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x

	if (i<=size(outjgradpsi,1)) then
        !jgradpsi=         A(:)*grad(:)+gradApsi
        !jgradpsi=jgradpsi+A(:)*grad(:)+gradApsi
        !jgradpsi=jgradpsi+A(:)*grad(:)+gradApsi
  outjgradpsi(i) = cmplx(0.d0,1.d0)*(inlambdaj(i,1)+2*muj*(curr1(i,1)-curr0(i,1)))*ingradpsi(i,1)+ingradApsi(i,1)&
		  +cmplx(0.d0,1.d0)*(inlambdaj(i,2)+2*muj*(curr1(i,2)-curr0(i,2)))*ingradpsi(i,2)+ingradApsi(i,2)&
		  +cmplx(0.d0,1.d0)*(inlambdaj(i,3)+2*muj*(curr1(i,3)-curr0(i,3)))*ingradpsi(i,3)+ingradApsi(i,3)
	end if
end subroutine cucompjgradpsi   

!attributes(device)  exec. GPU called GPU 
!attributes(host)  exec. CPU called by CPU
!attributes(global) exec. GPU called by CPU

attributes(global) subroutine cucompcurrpenality(outApsi,inpsi,inlambdaj,curr0,curr1,muj)
	implicit none
	complex(DP)  :: outApsi(:,:)
	complex(DP)  :: inpsi(:)
	real(DP)     :: inlambdaj(:,:)
	real(DP)     :: curr0(:,:)
	real(DP)     :: curr1(:,:)
	real(DP),value:: muj

	integer :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x

	if (i<=size(outApsi,1)) then
     ! A=eye*(2*muj* (curr1(:,1)-curr0(:,1))+lambdaj(:,1) )
     ! Apsi=A(:)*psi1(:,N)
     ! A=eye*(2*muj* (curr1(:,2)-curr0(:,2))+lambdaj(:,2) )
     ! Apsi=A(:)*psi1(:,N)
     !  A=eye*(2*muj* (curr1(:,3)-curr0(:,3))+lambdaj(:,3) )
     !  Apsi=A(:)*psi1(:,N)
	  outApsi(i,1) = cmplx(0.d0,1.d0)*(inlambdaj(i,1)+2*muj*(curr1(i,1)-curr0(i,1)))*inpsi(i)
	  outApsi(i,2) = cmplx(0.d0,1.d0)*(inlambdaj(i,2)+2*muj*(curr1(i,2)-curr0(i,2)))*inpsi(i)
	  outApsi(i,2) = cmplx(0.d0,1.d0)*(inlambdaj(i,2)+2*muj*(curr1(i,3)-curr0(i,3)))*inpsi(i)
	end if
end subroutine cucompcurrpenality


attributes(global) subroutine cucomprhopenality(outALpsi,inpsi,inlambda,rhotot0,rhotot1,mu,ispin)
	implicit none
	complex(DP)  :: outALpsi(:)
	complex(DP)  :: inpsi(:)
	real(DP)     :: inlambda(:,:)
	real(DP)     :: rhotot0(:,:)
	real(DP)     :: rhotot1(:,:)
	real(DP),value:: mu
	integer	,value:: ispin

	integer :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x

	if (i<=size(outALpsi,1)) then
     ! compute lambda*psi
     !      ALpsi=lambda(:,ispin(N))*psi1(:,N)!real space
     ! compute rho*psi
     !      ALpsi=ALpsi+mu*(rhototloc( :,ispin(N) )-rhotot0( :,ispin(N) ) )*psi1(:,N)!real space
		outALpsi(i) = (inlambda(i,ispin)+mu*(rhotot1(i,ispin)-rhotot0(i,ispin)))*inpsi(i)
	end if
end subroutine cucomprhopenality

attributes(global) subroutine cupdatelangmult(lbd,lbdj,rhotot0,rhotot1,curr0,curr1,mu,muj)
	implicit none
	real(DP)     :: lbd(:,:)
	real(DP)     :: lbdj(:,:)
	real(DP)     :: rhotot0(:,:)
	real(DP)     :: rhotot1(:,:)
	real(DP)     :: curr0(:,:)
	real(DP)     :: curr1(:,:)
	real(DP),value :: mu,muj

	integer :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x

	if (i<=size(lbd,1)) then
	!lambda=lambda+2.d0*(rhotot1-rhotot0)*mu
  	!lambdaj=lambdaj+2.d0*(curr1-curr0)*muj
	lbd(i,1)=lbd(i,1)+2.d0*(rhotot1(i,1)-rhotot0(i,1))*mu
	lbd(i,2)=lbd(i,2)+2.d0*(rhotot1(i,2)-rhotot0(i,2))*mu
	lbdj(i,1)=lbdj(i,1)+2.d0*(curr1(i,1)-curr0(i,1))*muj
	lbdj(i,2)=lbdj(i,2)+2.d0*(curr1(i,2)-curr0(i,2))*muj
	lbdj(i,3)=lbdj(i,3)+2.d0*(curr1(i,3)-curr0(i,3))*muj
	end if

end subroutine cupdatelangmult

attributes(global) subroutine cusumequili(rhotot0,rhotot1,curr0,curr1,err,errj,ma1,sumcurr0,sumcurr1)
	implicit none
	REAL(DP) 	:: rhotot0(:,:),rhotot1(:,:),curr0(:,:),curr1(:,:)
	REAL(DP)	:: err,errj,ma1,sumcurr0,sumcurr1
	REAL(DP),shared :: temp(blockDim%x)

	integer :: i,j,inext

	i =threadIdx%x
	j=i
	temp(j) = 0.D0
	do while (i<=size(rhotot0,1))
	if (blockIdx%x==1) then
		temp(j) = temp(j) + abs(rhotot1(i,1)-rhotot0(i,1))&
				  + abs(rhotot1(i,2)-rhotot0(i,2))
	elseif (blockIdx%x==2) then
		temp(j) = temp(j) + abs(curr1(i,1)-curr0(i,1))&
				  + abs(curr1(i,2)-curr0(i,2))&
				  + abs(curr1(i,3)-curr0(i,3))
	elseif (blockIdx%x==3) then
		temp(j) = temp(j) + abs(rhotot1(i,1))&
				  + abs(rhotot1(i,2))
	elseif (blockIdx%x==4) then
		temp(j) = temp(j) + abs(curr0(i,1))&
				  + abs(curr0(i,2))&
				  + abs(curr0(i,3))
	else
		temp(j) = temp(j) + abs(curr1(i,1))&
				  + abs(curr1(i,2))&
				  + abs(curr1(i,3))
	endif
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
if (j==1) then
	if (blockIdx%x==1) then
		err = temp(1)
	elseif (blockIdx%x==2) then
		errj= temp(1)
	elseif (blockIdx%x==3) then
		ma1= temp(1)*dvol_d
	elseif (blockIdx%x==4) then
		sumcurr0= temp(1)*dvol_d
	else
		sumcurr1= temp(1)*dvol_d
	endif
endif
end subroutine cusumequili


attributes(global) subroutine cucalcrhonospin(inpsi,occ,outrho) 
	implicit none
	complex(DP)  :: inpsi(:,:)
	real(DP)     :: occ(:)
	real(DP)     :: outrho(:)

	integer  :: i,j

	 i = (blockIdx%x-1)*blockDim%x + threadIdx%x ! MAX kdfull2
	 !j = (blockIdx%y-1)*blockDim%y + threadIdx%y ! MAX nstate/2

	if (i<=size(inpsi,1)) then
		outrho(i)		= 0.d0
		outrho(i+size(inpsi,1)) = 0.d0
	    j=1
	    do while (j<= size(occ,1))
	    outrho(i)			= outrho(i)+ occ(j)*CONJG(inpsi(i,j))*inpsi(i,j)
	    j =j+1
	    end do
	    outrho(i+size(outrho,1)/2) 	= 0.d0
	end if

end subroutine cucalcrhonospin


attributes(global) subroutine affectpsiwork(inpsiarray,outpsiwork,innb)
	implicit none
	complex(DP)	:: inpsiarray(:,:)   
	complex(DP)     :: outpsiwork(:)
	integer(DP)	,value	:: innb 
	!----------------------------
	
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(inpsiarray,1)) then
		outpsiwork(i) = inpsiarray(i,innb) 
	end if
end subroutine affectpsiwork

attributes(global) subroutine affectpsi(inpsiarray,outpsiwork)
	implicit none
	complex(DP)	:: inpsiarray(:)
	complex(DP)     :: outpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(inpsiarray,1)) then
		outpsiwork(i) = inpsiarray(i) 
	end if
end subroutine affectpsi

attributes(global) subroutine addcontribtaylor(outpsiarray,inpsiwork,innb,cfactor)
	implicit none
	complex(DP)	:: outpsiarray(:,:)
	complex(DP)     :: inpsiwork(:)
	integer(DP)	,value	:: innb 
	complex(DP)	,value  :: cfactor
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiarray,1)) then
		outpsiarray(i,innb) = outpsiarray(i,innb)++cfactor*inpsiwork(i)
	end if
end subroutine addcontribtaylor

attributes(global) subroutine addenscontribtaylor(outpsiarray,inpsiwork,cfactor)
	implicit none
	complex(DP)	:: outpsiarray(:,:)
	complex(DP)     :: inpsiwork(:,:)
	complex(DP)	,value  :: cfactor
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiarray,1)) then
		outpsiarray(i,:) = outpsiarray(i,:)+cfactor*inpsiwork(i,:)
	end if
end subroutine addenscontribtaylor

attributes(global) subroutine multpsiwork(infield,outpsiwork)
	implicit none
	real(DP)	:: infield(:)
	complex(DP)     :: outpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = infield(i)*outpsiwork(i)
	end if
end subroutine multpsiwork

attributes(global) subroutine normpsiwork(cnorm,outpsiwork)
	implicit none
	complex(DP)	:: cnorm
	complex(DP)     :: outpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = outpsiwork(i)/sqrt(real(cnorm))
	end if
end subroutine normpsiwork

attributes(global) subroutine affectmultpsiwork(infield,outpsiwork,inpsiwork)
	implicit none
	real(DP)	:: infield(:)
	complex(DP)     :: outpsiwork(:)
	complex(DP)     :: inpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = infield(i)*inpsiwork(i)
	end if
end subroutine affectmultpsiwork

attributes(global) subroutine affectmultpsiworkcst(infield,outpsiwork,inpsiwork,cst)
	implicit none
	real(DP)	:: infield(:)
	complex(DP)     :: outpsiwork(:)
	complex(DP)     :: inpsiwork(:)
	real(DP) , value :: cst
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = cst*infield(i)*inpsiwork(i)
	end if
end subroutine affectmultpsiworkcst

attributes(global) subroutine addmultpsiwork(infield,outpsiwork,inpsiwork)
	implicit none
	real(DP)	:: infield(:)
	complex(DP)     :: outpsiwork(:)
	complex(DP)     :: inpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = outpsiwork(i)+infield(i)*inpsiwork(i)
	end if
end subroutine addmultpsiwork


attributes(global) subroutine addweightedpsiwork(inconst,outpsiwork,inpsiwork)
	implicit none
	real(DP) 	:: inconst
	complex(DP)     :: outpsiwork(:)
	complex(DP)     :: inpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = outpsiwork(i)-inconst*inpsiwork(i)
	end if
end subroutine addweightedpsiwork

attributes(global) subroutine addvaluedpsiwork(inconst,outpsiwork,inpsiwork)
	implicit none
	real(DP) ,value	:: inconst
	complex(DP)     :: outpsiwork(:)
	complex(DP)     :: inpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = outpsiwork(i)-inconst*inpsiwork(i)
	end if
end subroutine addvaluedpsiwork

attributes(global) subroutine addcomplexpsiwork(outpsiwork,inpsiwork,inconst)
	implicit none
	complex(DP) 	:: inconst
	complex(DP)     :: outpsiwork(:)
	complex(DP)     :: inpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = outpsiwork(i)-inconst*inpsiwork(i)
	end if
end subroutine addcomplexpsiwork


attributes(global) subroutine addpsiwork(outpsiwork,inpsiwork)
	implicit none
	complex(DP)     :: outpsiwork(:)
	complex(DP)     :: inpsiwork(:)
	integer :: i

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(outpsiwork,1)) then
		outpsiwork(i) = outpsiwork(i)+inpsiwork(i)
	end if
end subroutine addpsiwork

attributes(global) subroutine cugetecback(rhojel,chpcoul,ecbackout)
	implicit none
	REAL(DP) 	:: rhojel(:)
	REAL(DP)	:: chpcoul(:)
	REAL(DP)	:: ecbackout
	REAL(DP),shared :: temp(blockDim%x)

	integer :: i,j,inext

	i =threadIdx%x
	j=i
	temp(j) = 0.D0
	do while (i<=size(rhojel,1))
	temp(j) = temp(j) - chpcoul(i)*rhojel(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
if (j==1) ecbackout=temp(1)
end subroutine cugetecback

attributes(global) subroutine cugetecbackion(rhojel,chpcoul,ecbackout)
	implicit none
	REAL(DP) 	:: rhojel(:)
	REAL(DP)	:: chpcoul(:)
	REAL(DP)	:: ecbackout
	REAL(DP),shared :: temp(blockDim%x)

	integer :: i,j,inext

	i =threadIdx%x
	j=i
	temp(j) = 0.D0
	do while (i<=size(rhojel,1)/2)
	temp(j) = temp(j) - chpcoul(i)*rhojel(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
if (j==1) ecbackout=temp(1)
end subroutine cugetecbackion

attributes(global) subroutine cugetecrhoion(inrho,chpcoul,potion,ecbackout)
	implicit none
	REAL(DP) 	:: inrho(:)
	REAL(DP)	:: chpcoul(:)
	REAL(DP)	:: potion(:)
	REAL(DP)	:: ecbackout
	REAL(DP),shared :: temp(blockDim%x)

	integer :: i,j,inext

	i =threadIdx%x
	j=i
	temp(j) = 0.D0
	do while (i<=size(inrho,1)/2)
	temp(j) = temp(j) + chpcoul(i)*inrho(i) - potion(i)*inrho(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
if (j==1) ecbackout=temp(1)
end subroutine cugetecrhoion

attributes(global) subroutine cugetecrho(inrho,chpcoul,ecbackout)
	implicit none
	REAL(DP) 	:: inrho(:)
	REAL(DP)	:: chpcoul(:)
	REAL(DP)	:: ecbackout
	REAL(DP),shared :: temp(blockDim%x)

	integer :: i,j,inext

	i =threadIdx%x
	j=i
	temp(j) = 0.D0
	do while (i<=size(inrho,1)/2)
	temp(j) = temp(j) + chpcoul(i)*inrho(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
if (j==1) ecbackout=temp(1)
end subroutine cugetecrho

attributes(global) subroutine cukineticenerg(inkpsi,inakv,outekin,sum0ex)
	implicit none
	COMPLEX(DP)	:: inkpsi(:,:)
	REAL(DP)	:: inakv(:)
	REAL(DP)	:: outekin(:)
	REAL(DP),value 	:: sum0ex
	REAL(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext

	i =threadIdx%x
	j=i
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + (REAL(inkpsi(i,blockIdx%x),DP)*REAL(inkpsi(i,blockIdx%x),DP) +&
			    AIMAG(inkpsi(i,blockIdx%x))*AIMAG(inkpsi(i,blockIdx%x)))*inakv(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outekin(blockIdx%x)=temp(1)/sum0ex
	
end subroutine cukineticenerg

attributes(global) subroutine cukineticvar(inkpsi,inakv,outekin,sum0ex)
	implicit none
	COMPLEX(DP)	:: inkpsi(:,:)
	REAL(DP)	:: inakv(:)
	REAL(DP)	:: outekin(:)
	REAL(DP),value 	:: sum0ex
	REAL(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext

	i =threadIdx%x
	j=i
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + (REAL(inkpsi(i,blockIdx%x),DP)*REAL(inkpsi(i,blockIdx%x),DP) +&
			    AIMAG(inkpsi(i,blockIdx%x))*AIMAG(inkpsi(i,blockIdx%x)))*inakv(i)*inakv(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outekin(blockIdx%x)=temp(1)/sum0ex
	
end subroutine cukineticvar

!warning nstateup = nstatedown = nstate/2
attributes(global) subroutine cucalc_epotloc(inkpsi,inpot,outepot)
	implicit none
	COMPLEX(DP)	:: inkpsi(:,:)
	REAL(DP)	:: inpot(:)
	REAL(DP)	:: outepot(:)

	COMPLEX(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext,nb,shft

	nb = blockIdx%x
	i  = threadIdx%x
	j  = i 
	shft = (ispin_d(nb)-1)*size(inkpsi,1)

	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + CONJG(inkpsi(i,nb))*inpot(i+shft)*inkpsi(i,nb)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outepot(nb)=REAL(temp(1),DP)*dvol_d
	
end subroutine cucalc_epotloc

attributes(global) subroutine cucalc_epotnonloc(inpsi,inkpsi,inpot,outepot)
	implicit none
	COMPLEX(DP)	:: inpsi(:,:)
	COMPLEX(DP)	:: inkpsi(:,:)
	REAL(DP)	:: inpot(:)
	REAL(DP)	:: outepot(:)
	
	COMPLEX(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext,nb,shft

	nb = blockIdx%x
	i  = threadIdx%x
	j  = i 

    !  ishift = (ispin(nrel2abs(nb))-1)*nxyz
	!if (nb > size(outepot,1)/2) then 
		shft = (ispin_d(nb)-1)*size(inkpsi,1)
	!else
	!	shft = 0
	!end if
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + CONJG(inpsi(i,nb))*(inkpsi(i,nb)+inpot(i+shft)*inpsi(i,nb))
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outepot(nb)=REAL(temp(1),DP)*dvol_d
	
end subroutine cucalc_epotnonloc


attributes(global) subroutine cucalc_epotsp1nonloc(inpsi,inkpsi,inpot,outepot)
	implicit none
	COMPLEX(DP)	:: inpsi(:,:)
	COMPLEX(DP)	:: inkpsi(:,:)
	REAL(DP)	:: inpot(:)
	REAL(DP)	:: outepot(:)

	COMPLEX(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext,nb

	nb = blockIdx%x
	i  = threadIdx%x
	j  = i 
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + CONJG(inpsi(i,nb))*(inkpsi(i,nb)+inpot(i)*inpsi(i,nb))
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outepot(nb)=REAL(temp(1),DP)*dvol_d
	
end subroutine cucalc_epotsp1nonloc

attributes(global) subroutine cuwfovlp(inkpsi,inpsi,outovp)
	implicit none
	COMPLEX(DP)	:: inkpsi(:)
	COMPLEX(DP)	:: inpsi(:)
	REAL(DP)	:: outovp
	COMPLEX(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext

	i  = threadIdx%x
	j  = i 
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + CONJG(inkpsi(i))*inpsi(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outovp=REAL(temp(1),DP)*dvol_d
	
end subroutine cuwfovlp

attributes(global) subroutine cuwfovlptot(inkpsi,inpsi,outovp)
	implicit none
	COMPLEX(DP)	:: inkpsi(:,:)
	COMPLEX(DP)	:: inpsi(:,:)
	REAL(DP)	:: outovp(:)
	COMPLEX(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext,n

	i  = threadIdx%x
	n  = blockIdx%x
	j  = i 
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + CONJG(inkpsi(i,n))*inpsi(i,n)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outovp(n)=REAL(temp(1),DP)*dvol_d
	
end subroutine cuwfovlptot

attributes(global) subroutine cucomplexwfovlp(inkpsi,inpsi,outovp)
	implicit none
	COMPLEX(DP)	:: inkpsi(:)
	COMPLEX(DP)	:: inpsi(:)
	COMPLEX(DP)	:: outovp
	COMPLEX(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext

	i  = threadIdx%x
	j  = i 
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + CONJG(inkpsi(i))*inpsi(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outovp=temp(1)*dvol_d
	
end subroutine cucomplexwfovlp

attributes(global) subroutine cutwobodykernel(inkpsi,inpsi,outovp,xval,yval,zval)
	implicit none
	COMPLEX(DP)	:: inkpsi(:)
	COMPLEX(DP)	:: inpsi(:)
	COMPLEX(DP)	:: outovp(:)
	COMPLEX(DP),shared :: temp(blockDim%x)
	REAL(DP)	:: xval(:),yval(:),zval(:)
	
	integer :: i,j,inext,r
	real(DP) :: dist
	
	r  = blockIdx%x 
	i  = threadIdx%x
	j  = i 
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	dist = sqrt((xval(i)-xval(r))**2+(yval(i)-yval(r))**2+(zval(i)-zval(r))**2)
	if(r/=i) temp(j) = temp(j) + CONJG(inkpsi(i))*inpsi(i)*exp(-dist/1.3)/dist
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outovp(r)=temp(1)!*dvol_d
	
end subroutine cutwobodykernel

attributes(global) subroutine cucdot_product(inkpsi,inpsi,outovp,dvolin)
	implicit none
	COMPLEX(DP)	:: inkpsi(:)
	COMPLEX(DP)	:: inpsi(:)
	REAL(DP)	:: outovp
	REAL(DP) ,value :: dvolin
	COMPLEX(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext

	i  = threadIdx%x
	j  = i 
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	temp(j) = temp(j) + conjg(inkpsi(i))*inpsi(i)
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outovp=REAL(temp(1),DP)*dvolin
	
end subroutine cucdot_product

attributes(global) subroutine cugetnonlocenergies(inkpsi,inpsi,outnonloc)
	implicit none
	COMPLEX(DP)	:: inkpsi(:,:)
	COMPLEX(DP)	:: inpsi(:,:)
	REAL(DP)	:: outnonloc(:)
	REAL(DP),shared :: temp(blockDim%x)
	
	integer :: i,j,inext
	i  = threadIdx%x
	j  = i 
	temp(j) = 0.D0
	do while (i<=size(inkpsi,1))
	!temp(j) = temp(j) + CONJG(inkpsi(i,blockIdx%x))*inpsi(i,blockIdx%x)
	temp(j) = temp(j) + REAL(inkpsi(i,blockIdx%x),DP)*REAL(inpsi(i,blockIdx%x),DP) +&
			  + AIMAG(inkpsi(i,blockIdx%x))*AIMAG(inpsi(i,blockIdx%x))
	i = i + blockDim%x
	end do

	call syncthreads()
	inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
if (j==1) outnonloc(blockIdx%x)=temp(1)*dvol_d

end subroutine cugetnonlocenergies

attributes(global) subroutine cucoordstep(inoutpsi,inaloc,indt)
	implicit none
	COMPLEX(DP)	:: inoutpsi(:,:)
	REAL(DP)	:: inaloc(:)
	REAL(DP) ,value :: indt

	integer :: i,shft,nb

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(inoutpsi,1)) then
	nb = 1 
	do while (nb<=size(inoutpsi,2))
		if (nb > size(inoutpsi,2)/2) then 
			shft = size(inoutpsi,1)	
		else
			shft = 0
		end if
	inoutpsi(i,nb) = CMPLX(COS(-indt*inaloc(i+shft)),SIN(-indt*inaloc(i+shft)),DP)*inoutpsi(i,nb)
	nb = nb +1	
	end do
	end if
end subroutine cucoordstep

attributes(global) subroutine cucoordstepsp1(inoutpsi,inaloc,indt)
	implicit none
	COMPLEX(DP)	:: inoutpsi(:,:)
	REAL(DP)	:: inaloc(:)
	REAL(DP) ,value :: indt

	integer :: i,nb

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(inoutpsi,1)) then
	nb = 1 
	do while (nb<=size(inoutpsi,2))
	inoutpsi(i,nb) = CMPLX(COS(-indt*inaloc(i)),SIN(-indt*inaloc(i)),DP)*inoutpsi(i,nb)
	nb = nb +1	
	end do
	end if
end subroutine cucoordstepsp1

attributes(global) subroutine cuapplymask(inoutpsi,mask)
	implicit none
	COMPLEX(DP)	:: inoutpsi(:,:)
	REAL(DP)	:: mask(:)

	integer :: i,nb

	i= (blockIdx%x-1)*blockDim%x + threadIdx%x 

	if (i<=size(inoutpsi,1)) then
	nb = 1 
	do while (nb<=size(inoutpsi,2))
	inoutpsi(i,nb) = mask(i)*inoutpsi(i,nb)
	nb = nb +1	
	end do
	end if
end subroutine cuapplymask


attributes(global) subroutine curhofld(inrhoinp,outrhokr,outrhoki,nx,ny,nz)
	implicit none
	REAL(DP)		:: inrhoinp(:)
	REAL(DP)		:: outrhokr(:)
	REAL(DP)		:: outrhoki(:)
	INTEGER  , value 	:: nx,ny,nz

	integer  :: i,j,i1, i2, i3

	i1= (blockIdx%x-1)*blockDim%x + threadIdx%x !X
	i2= (blockIdx%y-1)*blockDim%y + threadIdx%y !Y
	i3= (blockIdx%z-1)*blockDim%z + threadIdx%z !Z
	if (i1<=nx .and. i2<=ny .and. i3 <= nz) then
	i = 4*(i3-1)*nx*ny+&
	    2*(i2-1)*nx + i1

	j = (i3-1)*nx*ny+&
	    (i2-1)*nx + i1
	outrhoki(i)=0.d0
	outrhokr(i)=inrhoinp(j)
	end if
end subroutine curhofld

attributes(global) subroutine reset1d(inarr)
IMPLICIT NONE
COMPLEX(DP) :: inarr(:)

integer   :: i

i = (blockIdx%x-1)*blockDim%x + threadIdx%x

if (i<=size(inarr,1)) then
inarr(i) = (0.d0,0.d0)
end if

end subroutine reset1d

attributes(global) subroutine reset2d(inarr)
IMPLICIT NONE
COMPLEX(DP) :: inarr(:,:)

integer   :: i

i = (blockIdx%x-1)*blockDim%x + threadIdx%x

if (i<=size(inarr,1)) then
inarr(i,blockIdx%y) = (0.d0,0.d0)
end if

end subroutine reset2d

attributes(global) subroutine reset2dr(inarr)
IMPLICIT NONE
REAL(DP) :: inarr(:,:)

integer   :: i

i = (blockIdx%x-1)*blockDim%x + threadIdx%x

if (i<=size(inarr,1)) then
inarr(i,blockIdx%y) = 0.d0
end if

end subroutine reset2dr


attributes(global) SUBROUTINE cutoteffnl4c(psi,aux,ion)

!     **************************

USE cuparams ,only: np_d,h0_11g_d,h0_22g_d,h1_11g_d,h0_33g_d,h1_22g_d,h2_11g_d,h0_12g_d,&
		  ifin_d,p0_1_d,p0_2_d,p1_2_d,p1_1x_d,p1_1y_d,p1_1z_d,p0_3_d,p1_2x_d,p1_2y_d,&
		  p1_2z_d,p2_xy_d,p2_xz_d,p2_yz_d,p2_xy2_d,p2_z2_d,icount_d
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)     :: psi(:,:)
COMPLEX(DP), INTENT(IN OUT) :: aux(:,:)
COMPLEX(DP) , shared 	    :: temp(blockDim%x)
INTEGER,INTENT(IN) , value  :: ion
!COMPLEX(DP) , INTENT(OUT)   :: outsum(:)

INTEGER :: i,j,ii, nfin,inext,l
REAL(DP) :: h0_11, h0_12, h0_13, h0_21, h0_22, h0_23, h0_31, h0_32, h0_33, h1_11, h1_12, h1_21, h1_22, h2_11

COMPLEX(DP) ,shared :: sumr01,sumr02,sumr03,sumr04,sumr05,sumr06
COMPLEX(DP) ,shared :: sumr07,sumr08,sumr09
COMPLEX(DP) ,shared :: sumr16,sumr17,sumr18,sumr19,sumr20
COMPLEX(DP) :: f0_1,f0_2,f0_3,f1_1x,f1_1y,f1_1z,f1_2x,f1_2y,f1_2z, &
               f2_xy,f2_xz,f2_yz,f2_xy2,f2_z2


REAL(DP),PARAMETER :: fac1_12=-0.422577127364258D0 !   -0.5D0*SQRT(5D0/7D0)
REAL(DP),PARAMETER :: fac0_23=-0.629940788348712D0 !  -0.5D0*SQRT(100D0/63D0)
REAL(DP),PARAMETER :: fac0_13= 0.243975018237133D0 !  0.5D0*SQRT(5D0/21D0)
REAL(DP),PARAMETER :: fac0_12=-0.387298334620742D0 !  -0.5D0*SQRT(3D0/5D0)
REAL(DP),PARAMETER :: pi     = 3.141592653589793D0
 

h0_11=h0_11g_d(np_d(ion))

h0_22=h0_22g_d(np_d(ion))
h1_11=h1_11g_d(np_d(ion))

h0_33=h0_33g_d(np_d(ion))
h1_22=h1_22g_d(np_d(ion))
h2_11=h2_11g_d(np_d(ion))

IF(h0_12g_d(np_d(ion)).gt.-1D10) THEN
  h0_12=h0_12g_d(np_d(ion))
ELSE
  h0_12 = fac0_12*h0_22
ENDIF
h0_21 = h0_12

h0_13 = fac0_13*h0_33
h0_31 = h0_13

h0_23 = fac0_23*h0_33
h0_32 = h0_23

h1_12 = fac1_12*h1_22
h1_21 = h1_12

i= threadIdx%x

!  initialize the result of the integration:
if (i==1) then
sumr01 = 0D0
sumr02 = 0D0
sumr03 = 0D0
sumr04 = 0D0
sumr05 = 0D0
sumr06 = 0D0
sumr07 = 0D0
sumr08 = 0D0
sumr09 = 0D0
sumr16 = 0D0
sumr17 = 0D0
sumr18 = 0D0
sumr19 = 0D0
sumr20 = 0D0
end if
call syncthreads()
nfin=ifin_d(ion)
j = i

!  compute all the integrals over r'-dependant part of the
!  wavefunction times the r'-dependent part of the pseudo:

IF(ABS(h0_11) > 1.0D-20) THEN

	temp(j) = (0D0,0D0)
	call syncthreads()
	do while (i<= nfin)
	ii = icount_d(i,ion)
	temp(j) = temp(j)+ psi(ii,blockIdx%x)*p0_1_d(i,ion)
        i = i + blockDim%x
	end do
	call syncthreads()
	inext = blockDim%x/2
	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
	if (j==1) sumr01=temp(1)*dvol_d
	call syncthreads()
END IF
IF(ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > 1.0D-20) THEN

do l = 1 , 4
	temp(j) = (0D0,0D0)
	call syncthreads()
	i= threadIdx%x
	do while (i<= nfin)
	ii = icount_d(i,ion)
	if(l==1) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p0_2_d(i,ion)
	if(l==2) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p1_1x_d(i,ion)
	if(l==3) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p1_1y_d(i,ion)
	if(l==4) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p1_1z_d(i,ion)
        i = i + blockDim%x
	end do
	call syncthreads()
	inext = blockDim%x/2
	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
	if (j==1) then
	 if (l==1) sumr02=temp(1)*dvol_d
	 if (l==2) sumr04=temp(1)*dvol_d
	 if (l==3) sumr05=temp(1)*dvol_d
	 if (l==4) sumr06=temp(1)*dvol_d
	end if
	call syncthreads()
end do

END IF
IF(ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > 1.0D-20) THEN
do l = 1 , 9 
	temp(j) = (0D0,0D0)
	call syncthreads()
	i= threadIdx%x
	do while (i<= nfin)
	ii = icount_d(i,ion)
	if(l==1) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p0_3_d(i,ion)
	if(l==2) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p1_2x_d(i,ion)
	if(l==3) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p1_2y_d(i,ion)
	if(l==4) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p1_2z_d(i,ion)
	if(l==5) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p2_xy_d(i,ion)
	if(l==6) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p2_xz_d(i,ion)
	if(l==7) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p2_yz_d(i,ion)
	if(l==8) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p2_xy2_d(i,ion)
	if(l==9) temp(j) = temp(j)+ psi(ii,blockIdx%x)*p2_z2_d(i,ion)
        i = i + blockDim%x
	end do
	call syncthreads()
	inext = blockDim%x/2
	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
	if (j==1) then
	 if (l==1) sumr03=temp(1)*dvol_d
	 if (l==2) sumr07=temp(1)*dvol_d
	 if (l==3) sumr08=temp(1)*dvol_d
	 if (l==4) sumr09=temp(1)*dvol_d
	 if (l==5) sumr16=temp(1)*dvol_d
	 if (l==6) sumr17=temp(1)*dvol_d
	 if (l==7) sumr18=temp(1)*dvol_d
	 if (l==8) sumr19=temp(1)*dvol_d
	 if (l==9) sumr20=temp(1)*dvol_d
	end if
	call syncthreads()

end do
END IF


!if(threadIdx%x==1) outsum(blockIdx%x) = sumr01
!  Now, the result of the r' integral is multiplied with the
!  r-dependent part and summed up over all l-contributions 

!ii = icount_d(i,ion)
IF(ABS(h0_11) > 1.0D-20) THEN
	i= threadIdx%x
	f0_1 = 2.0D0*( h0_11*sumr01 + h0_12*sumr02 + h0_13*sumr03 )/(4D0*pi)
	do while(i<=nfin)
        ii = icount_d(i,ion)
        aux(ii,blockIdx%x)= aux(ii,blockIdx%x) + p0_1_d(i,ion)*f0_1
        i = i + blockDim%x
        end do
END IF
IF(ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > 1.0D-20) THEN
	i= threadIdx%x
	f0_2 = 2.0D0*( h0_21*sumr01 + h0_22*sumr02 + h0_23*sumr03 )/(4D0*pi)
	f1_1x = 2.0D0*( h1_11*sumr04 + h1_12*sumr07 )*3D0/(4D0*pi)
	f1_1y = 2.0D0*( h1_11*sumr05 + h1_12*sumr08 )*3D0/(4D0*pi)
	f1_1z = 2.0D0*( h1_11*sumr06 + h1_12*sumr09 )*3D0/(4D0*pi)
	do while(i<=nfin)
        ii = icount_d(i,ion)
	aux(ii,blockIdx%x)= aux(ii,blockIdx%x) + p0_2_d(i,ion)*f0_2   &
	+(p1_1x_d(i,ion)*f1_1x + p1_1y_d(i,ion)*f1_1y + p1_1z_d(i,ion)*f1_1z) 
        i = i + blockDim%x
        end do
END IF
IF(ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > 1.0D-20) THEN
	i= threadIdx%x
	f0_3 = 2.0D0*( h0_31*sumr01 + h0_32*sumr02 + h0_33*sumr03 /4D0*pi)
	f1_2x = 2.0D0*( h1_21*sumr04 + h1_22*sumr07 )*3D0/(4D0*pi)
	f1_2y = 2.0D0*( h1_21*sumr05 + h1_22*sumr08 )*3D0/(4D0*pi)
	f1_2z = 2.0D0*( h1_21*sumr06 + h1_22*sumr09 )*3D0/(4D0*pi)
	f2_xy = 2.0D0* h2_11*sumr16*15D0/(4D0*pi)
	f2_xz = 2.0D0* h2_11*sumr17*15D0/(4D0*pi)
	f2_yz = 2.0D0* h2_11*sumr18*15D0/(4D0*pi)
	f2_xy2 = 2.0D0* h2_11*sumr19*5D0/(16D0*pi)
	f2_z2 = 2.0D0* h2_11*sumr20*5D0/(16D0*pi)
	do while(i<=nfin)
        ii = icount_d(i,ion)
	aux(ii,blockIdx%x)= aux(ii,blockIdx%x) + p0_3_d(i,ion)*f0_3  &
	+(p1_2x_d(i,ion)*f1_2x + p1_2y_d(i,ion)*f1_2y + p1_2z_d(i,ion)*f1_2z) &
	+(p2_xy_d(i,ion)*f2_xy + p2_xz_d(i,ion)*f2_xz + p2_yz_d(i,ion)*f2_yz) &
	+(p2_xy2_d(i,ion)*f2_xy2 + p2_z2_d(i,ion)*f2_z2)
        i = i + blockDim%x
        end do
END IF


RETURN


END SUBROUTINE cutoteffnl4c


attributes(global) SUBROUTINE cueffnl4c(psi,aux,ion)

!     **************************

USE cuparams ,only: np_d,h0_11g_d,h0_22g_d,h1_11g_d,h0_33g_d,h1_22g_d,h2_11g_d,h0_12g_d,&
		  ifin_d,p0_1_d,p0_2_d,p1_2_d,p1_1x_d,p1_1y_d,p1_1z_d,p0_3_d,p1_2x_d,p1_2y_d,&
		  p1_2z_d,p2_xy_d,p2_xz_d,p2_yz_d,p2_xy2_d,p2_z2_d,icount_d,dvol_d
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)     :: psi(:)
COMPLEX(DP), INTENT(IN OUT) :: aux(:)
COMPLEX(DP) , shared 	    :: temp(blockDim%x)
INTEGER,INTENT(IN) , value  :: ion

INTEGER :: i,ii, nfin,inext,j,l
REAL(DP) :: h0_11, h0_12, h0_13, h0_21, h0_22, h0_23, h0_31, h0_32, h0_33, h1_11, h1_12, h1_21, h1_22, h2_11

COMPLEX(DP) ,shared :: sumr01,sumr02,sumr03,sumr04,sumr05,sumr06
COMPLEX(DP) ,shared :: sumr07,sumr08,sumr09
COMPLEX(DP) ,shared :: sumr16,sumr17,sumr18,sumr19,sumr20
COMPLEX(DP) :: f0_1,f0_2,f0_3,f1_1x,f1_1y,f1_1z,f1_2x,f1_2y,f1_2z, &
               f2_xy,f2_xz,f2_yz,f2_xy2,f2_z2


REAL(DP),PARAMETER :: fac1_12=-0.422577127364258D0 !   -0.5D0*SQRT(5D0/7D0)
REAL(DP),PARAMETER :: fac0_23=-0.629940788348712D0 !  -0.5D0*SQRT(100D0/63D0)
REAL(DP),PARAMETER :: fac0_13= 0.243975018237133D0 !  0.5D0*SQRT(5D0/21D0)
REAL(DP),PARAMETER :: fac0_12=-0.387298334620742D0 !  -0.5D0*SQRT(3D0/5D0)
REAL(DP),PARAMETER :: pi     = 3.141592653589793D0
 

h0_11=h0_11g_d(np_d(ion))

h0_22=h0_22g_d(np_d(ion))
h1_11=h1_11g_d(np_d(ion))

h0_33=h0_33g_d(np_d(ion))
h1_22=h1_22g_d(np_d(ion))
h2_11=h2_11g_d(np_d(ion))

IF(h0_12g_d(np_d(ion)).gt.-1D10) THEN
  h0_12=h0_12g_d(np_d(ion))
ELSE
  h0_12 = fac0_12*h0_22
ENDIF
h0_21 = h0_12

h0_13 = fac0_13*h0_33
h0_31 = h0_13

h0_23 = fac0_23*h0_33
h0_32 = h0_23

h1_12 = fac1_12*h1_22
h1_21 = h1_12

i= (blockIdx%x-1)*blockDim%x + threadIdx%x

!  initialize the result of the integration:
if (i==1) then
sumr01 = 0D0
sumr02 = 0D0
sumr03 = 0D0
sumr04 = 0D0
sumr05 = 0D0
sumr06 = 0D0
sumr07 = 0D0
sumr08 = 0D0
sumr09 = 0D0
sumr16 = 0D0
sumr17 = 0D0
sumr18 = 0D0
sumr19 = 0D0
sumr20 = 0D0
end if
call syncthreads()
nfin=ifin_d(ion)


!  compute all the integrals over r'-dependant part of the
!  wavefunction times the r'-dependent part of the pseudo:

IF(ABS(h0_11) > 1.0D-20) THEN

	temp(j) = (0D0,0D0)
	call syncthreads()
	if (i<= nfin) then
	ii = icount_d(i,ion)
	temp(j) = temp(j)+ psi(ii)*p0_1_d(i,ion)
	end if  
	call syncthreads()
	inext = blockDim%x/2
	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(i+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
	if (i==1) sumr01=temp(1)*dvol_d
	call syncthreads()

ELSE IF(ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > 1.0D-20) THEN

do l = 1 , 4
	temp(j) = (0D0,0D0)
	call syncthreads()
	do while (i<= nfin)
	ii = icount_d(i,ion)
	if(l==1) temp(j) = temp(j)+ psi(ii)*p0_2_d(i,ion)
	if(l==2) temp(j) = temp(j)+ psi(ii)*p1_1x_d(i,ion)
	if(l==3) temp(j) = temp(j)+ psi(ii)*p1_1y_d(i,ion)
	if(l==4) temp(j) = temp(j)+ psi(ii)*p1_1z_d(i,ion)
        i = i + blockDim%x
	end do
	call syncthreads()
	inext = blockDim%x/2
	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
	if (j==1) then
	 if (l==1) sumr02=temp(1)*dvol_d
	 if (l==2) sumr04=temp(1)*dvol_d
	 if (l==3) sumr05=temp(1)*dvol_d
	 if (l==4) sumr06=temp(1)*dvol_d
	end if
	call syncthreads()
end do

ELSE IF(ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > 1.0D-20) THEN
do l = 1 , 9 
	temp(j) = (0D0,0D0)
	call syncthreads()
	do while (i<= nfin)
	ii = icount_d(i,ion)
	if(l==1) temp(j) = temp(j)+ psi(ii)*p0_3_d(i,ion)
	if(l==2) temp(j) = temp(j)+ psi(ii)*p1_2x_d(i,ion)
	if(l==3) temp(j) = temp(j)+ psi(ii)*p1_2y_d(i,ion)
	if(l==4) temp(j) = temp(j)+ psi(ii)*p1_2z_d(i,ion)
	if(l==5) temp(j) = temp(j)+ psi(ii)*p2_xy_d(i,ion)
	if(l==6) temp(j) = temp(j)+ psi(ii)*p2_xz_d(i,ion)
	if(l==7) temp(j) = temp(j)+ psi(ii)*p2_yz_d(i,ion)
	if(l==8) temp(j) = temp(j)+ psi(ii)*p2_xy2_d(i,ion)
	if(l==9) temp(j) = temp(j)+ psi(ii)*p2_z2_d(i,ion)
        i = i + blockDim%x
	end do
	call syncthreads()
	inext = blockDim%x/2
	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
	call syncthreads()
	if (j==1) then
	 if (l==1) sumr03=temp(1)*dvol_d
	 if (l==2) sumr07=temp(1)*dvol_d
	 if (l==3) sumr08=temp(1)*dvol_d
	 if (l==4) sumr09=temp(1)*dvol_d
	 if (l==5) sumr16=temp(1)*dvol_d
	 if (l==6) sumr17=temp(1)*dvol_d
	 if (l==7) sumr18=temp(1)*dvol_d
	 if (l==8) sumr19=temp(1)*dvol_d
	 if (l==9) sumr20=temp(1)*dvol_d
	end if
	call syncthreads()

end do
END IF



!  Now, the result of the r' integral is multiplied with the
!  r-dependent part and summed up over all l-contributions 

ii = icount_d(i,ion)

IF(ABS(h0_11) > 1.0D-20) THEN
	f0_1 = 2.0D0*( h0_11*sumr01 + h0_12*sumr02 + h0_13*sumr03 )/(4D0*pi)
	aux(ii)= aux(ii) + p0_1_d(i,ion)*f0_1
END IF
IF(ABS(h1_11) + ABS(h0_22) + ABS(h0_12) > 1.0D-20) THEN
	f0_2 = 2.0D0*( h0_21*sumr01 + h0_22*sumr02 + h0_23*sumr03 )/(4D0*pi)
	f1_1x = 2.0D0*( h1_11*sumr04 + h1_12*sumr07 )*3D0/(4D0*pi)
	f1_1y = 2.0D0*( h1_11*sumr05 + h1_12*sumr08 )*3D0/(4D0*pi)
	f1_1z = 2.0D0*( h1_11*sumr06 + h1_12*sumr09 )*3D0/(4D0*pi)
	aux(ii)= aux(ii) + p0_2_d(i,ion)*f0_2   &
	+(p1_1x_d(i,ion)*f1_1x + p1_1y_d(i,ion)*f1_1y + p1_1z_d(i,ion)*f1_1z) 
END IF
IF(ABS(h2_11) + ABS(h1_22) + ABS(h0_33) > 1.0D-20) THEN
	f0_3 = 2.0D0*( h0_31*sumr01 + h0_32*sumr02 + h0_33*sumr03 )/(4D0*pi)
	f1_2x = 2.0D0*( h1_21*sumr04 + h1_22*sumr07 )*3D0/(4D0*pi)
	f1_2y = 2.0D0*( h1_21*sumr05 + h1_22*sumr08 )*3D0/(4D0*pi)
	f1_2z = 2.0D0*( h1_21*sumr06 + h1_22*sumr09 )*3D0/(4D0*pi)
	f2_xy = 2.0D0* h2_11*sumr16*15D0/(4D0*pi)
	f2_xz = 2.0D0* h2_11*sumr17*15D0/(4D0*pi)
	f2_yz = 2.0D0* h2_11*sumr18*15D0/(4D0*pi)
	f2_xy2 = 2.0D0* h2_11*sumr19*5D0/(16D0*pi)
	f2_z2 = 2.0D0* h2_11*sumr20*5D0/(16D0*pi)
	aux(ii)= aux(ii) + p0_3_d(i,ion)*f0_3  &
	+(p1_2x_d(i,ion)*f1_2x + p1_2y_d(i,ion)*f1_2y + p1_2z_d(i,ion)*f1_2z) &
	+(p2_xy_d(i,ion)*f2_xy + p2_xz_d(i,ion)*f2_xz + p2_yz_d(i,ion)*f2_yz) &
	+(p2_xy2_d(i,ion)*f2_xy2 + p2_z2_d(i,ion)*f2_z2)
END IF


RETURN
END SUBROUTINE cueffnl4c

!~ !     *****************************

attributes(global) subroutine cucopy1dto3d(inpsi,fft3d,innb)
implicit none
	COMPLEX(DP)		:: inpsi(:,:)
	COMPLEX(DP)		:: fft3d(:,:,:)
	integer		,value	:: innb

	integer  :: i,i1, i2, i3,nb

	i1= (blockIdx%x-1)*blockDim%x + threadIdx%x !X
	i2= (blockIdx%y-1)*blockDim%y + threadIdx%y !Y
	i3= (blockIdx%z-1)*blockDim%z + threadIdx%z !Z

	i = (i3-1)*size(fft3d,1)*size(fft3d,2)+&
	    (i2-1)*size(fft3d,1) + i1

	if (i1<=size(fft3d,1) .AND. i2<=size(fft3d,2) .AND. i3<=size(fft3d,3)) then
	!nb =1
	!do while (nb<=size(inpsi,2))
	fft3d(i1,i2,i3)=inpsi(i,innb)
	!nb=nb+1
	!end do
	end if

end subroutine cucopy1dto3d

attributes(global) subroutine cucopy1dto3dworking(inpsi,fft3d)
implicit none
	COMPLEX(DP)		:: inpsi(:)
	COMPLEX(DP)		:: fft3d(:,:,:)

	integer  :: i,i1, i2, i3

	i1= (blockIdx%x-1)*blockDim%x + threadIdx%x !X
	i2= (blockIdx%y-1)*blockDim%y + threadIdx%y !Y
	i3= (blockIdx%z-1)*blockDim%z + threadIdx%z !Z

	i = (i3-1)*size(fft3d,1)*size(fft3d,2)+&
	    (i2-1)*size(fft3d,1) + i1

	if (i1<=size(fft3d,1) .AND. i2<=size(fft3d,2) .AND. i3<=size(fft3d,3)) then
	fft3d(i1,i2,i3)=inpsi(i)
	end if

end subroutine cucopy1dto3dworking

attributes(global) subroutine cucopy3dto1d(inpsi,fft3d,coef,innb)
implicit none
	COMPLEX(DP)		:: inpsi(:,:)
	COMPLEX(DP)		:: fft3d(:,:,:)
	real(DP)	,value  :: coef
	integer		,value  :: innb

	integer  :: i,i1, i2, i3!,nb

	i1= (blockIdx%x-1)*blockDim%x + threadIdx%x !X
	i2= (blockIdx%y-1)*blockDim%y + threadIdx%y !Y
	i3= (blockIdx%z-1)*blockDim%z + threadIdx%z !Z

	i = (i3-1)*size(fft3d,1)*size(fft3d,2)+&
	    (i2-1)*size(fft3d,1) + i1

	if (i1<=size(fft3d,1) .AND. i2<=size(fft3d,2) .AND. i3<=size(fft3d,3)) then
	!nb =1
	!do while (nb<=size(inpsi,2))
	!inpsi(i,nb)=coef*fft3d(i1,i2,i3,nb)
	inpsi(i,innb)=coef*fft3d(i1,i2,i3)
	!nb=nb+1
	!end do

	end if

end subroutine cucopy3dto1d

attributes(global) subroutine cucopy3dto1dworking(inpsi,fft3d,coef)
implicit none
	COMPLEX(DP)		:: inpsi(:)
	COMPLEX(DP)		:: fft3d(:,:,:)
	real(DP)	,value  :: coef
	integer		,value  :: innb

	integer  :: i,i1, i2, i3

	i1= (blockIdx%x-1)*blockDim%x + threadIdx%x !X
	i2= (blockIdx%y-1)*blockDim%y + threadIdx%y !Y
	i3= (blockIdx%z-1)*blockDim%z + threadIdx%z !Z

	i = (i3-1)*size(fft3d,1)*size(fft3d,2)+&
	    (i2-1)*size(fft3d,1) + i1

	if (i1<=size(fft3d,1) .AND. i2<=size(fft3d,2) .AND. i3<=size(fft3d,3)) then

	inpsi(i)=coef*fft3d(i1,i2,i3)

	end if

end subroutine cucopy3dto1dworking


attributes(global) subroutine cumultakprop(fft3d,akpr)
implicit none
	COMPLEX(DP)		:: fft3d(:,:,:)
	COMPLEX(DP)	        :: akpr(:,:,:)

	integer  :: i1, i2, i3

	i1= (blockIdx%x-1)*blockDim%x + threadIdx%x !X
	i2= (blockIdx%y-1)*blockDim%y + threadIdx%y !Y
	i3= (blockIdx%z-1)*blockDim%z + threadIdx%z !Z

	if (i1<=size(fft3d,1) .AND. i2<=size(fft3d,2) .AND. i3<=size(fft3d,3)) then

	fft3d(i1,i2,i3)=akpr(i1,i2,i3)*fft3d(i1,i2,i3)

	end if

end subroutine cumultakprop


!attributes(global) subroutine curhofld(inrhoinp,outrhokr,outrhoki)
!	implicit none
!	REAL(DP)		:: inrhoinp(:)
!	REAL(DP)		:: outrhokr(:)
!	REAL(DP)		:: outrhoki(:)
!	!INTEGER  , value 	:: nxi,nyi,nzi
!
!	integer  :: i,j,i1, i2, i3
!
!	i1= (blockIdx%x-1)*blockDim%x + threadIdx%x
!	i2= (blockIdx%y-1)*blockDim%y + threadIdx%y
!	i3= (blockIdx%z-1)*blockDim%z + threadIdx%z
!	i = 4*(i3-1)*(griddim%x)*(blockDim%x)*(griddim%y)*(blockDim%y)+&
!	    2*(i2-1)*(griddim%x)*(blockDim%x) + i1
!
!	j = (i3-1)*(griddim%x)*(blockDim%x)*(griddim%y)*(blockDim%y)+&
!	    (i2-1)*(griddim%x)*(blockDim%x) + i1
!	outrhoki(i)=0.d0
!	outrhokr(i)=inrhoinp(j)
!
!end subroutine curhofld

attributes(global) subroutine constructrhospdeg(outrohspu,inrho,factup)
	implicit none
	REAL(DP)		:: outrohspu(:)
	REAL(DP)		:: inrho(:)
	REAL(DP),value		:: factup

	integer  :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	

	if (i<=size(inrho,1)/2) then
	outrohspu(i)=inrho(i)*factup
	outrohspu(i+size(inrho,1)/2)=1D0
	end if
end subroutine constructrhospdeg

attributes(global) subroutine constructrhospu(outrohspu,inrho,factup)
	implicit none
	REAL(DP)		:: outrohspu(:)
	REAL(DP)		:: inrho(:)
	REAL(DP),value		:: factup

	integer  :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	

	if (i<=size(inrho,1)/2) then
	outrohspu(i)=inrho(i)*(1D0+inrho(i+size(inrho,1)/2))*factup
	outrohspu(i+size(inrho,1)/2)=1D0
	end if
end subroutine constructrhospu

attributes(global) subroutine constructrhospd(outrohspd,inrho,factdw)
	implicit none
	REAL(DP)		:: outrohspd(:)
	REAL(DP)		:: inrho(:)
	REAL(DP),value		:: factdw

	integer  :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	

	if (i<=size(inrho,1)/2) then
	outrohspd(i)=inrho(i)*(1D0-inrho(i+size(inrho,1)/2))*factdw
	outrohspd(i+size(inrho,1)/2)=-1D0
	end if
end subroutine constructrhospd

attributes(global) subroutine sicalocaverdensity(outaloc,inrhoup,inrhodw)
	implicit none
	REAL(DP)		:: outaloc(:)
	REAL(DP)		:: inrhoup(:),inrhodw(:)

	integer  :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	

	if (i<=size(outaloc,1)/2) then
	outaloc(i) = outaloc(i) - inrhoup(i)
	outaloc(i+size(outaloc,1)/2) = outaloc(i+size(outaloc,1)/2) - inrhodw(i+size(outaloc,1)/2)
	end if
end subroutine sicalocaverdensity

attributes(global) subroutine sicalocaverdensitydeg(outaloc,inrhoup)
	implicit none
	REAL(DP)		:: outaloc(:)
	REAL(DP)		:: inrhoup(:)

	integer  :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	

	if (i<=size(outaloc,1)/2) then
	outaloc(i) = outaloc(i) - inrhoup(i)
	end if
end subroutine sicalocaverdensitydeg


attributes(global) subroutine compcoulombadsic(outrho1,outrho2,inrho)
	implicit none
	REAL(DP)		:: outrho1(:)
	REAL(DP)		:: outrho2(:),inrho(:)

	integer  :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	

	if (i<=size(inrho,1)/2) then
		outrho2(i)=  inrho(i)
    		outrho1(i)=  inrho(i)*inrho(i+size(inrho,1)/2)
	end if
end subroutine compcoulombadsic

attributes(global) subroutine alocfromcoulsic(outaloc,incoulsum,incouldif,facup,facdw)
	implicit none
	real(DP)  	 :: outaloc(:)
	real(DP)	 :: incoulsum(:),incouldif(:)
	real(DP) , value :: facup,facdw


	integer  :: i

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	

	if (i<=size(outaloc,1)/2) then
	outaloc(i) = outaloc(i) - 0.5D0*(incoulsum(i)+incouldif(i))*facup 
        outaloc(i + size(outaloc,1)/2) = outaloc(i + size(outaloc,1)/2) - facdw*0.5D0*(incoulsum(i)-incouldif(i))
	end if
end subroutine alocfromcoulsic


attributes(global) subroutine alocfromcoulsicsp1(outaloc,incoul,fac)
        implicit none
        real(DP)         :: outaloc(:)
        real(DP)         :: incoul(:)
        real(DP) , value :: fac


        integer  :: i

        i = (blockIdx%x-1)*blockDim%x+threadIdx%x


        if (i<=size(outaloc,1)/2) then
        outaloc(i) = outaloc(i) - fac*incoul(i)
        end if
end subroutine alocfromcoulsicsp1


attributes(global) subroutine cugetaloc(outaloc,inchpcoul,inchpdft)
	implicit none 
	REAL(DP)		:: outaloc(:)
	REAL(DP)		:: inchpcoul(:)
	REAL(DP)		:: inchpdft(:)

	integer  :: i


i = (blockIdx%x-1)*blockDim%x+threadIdx%x

if(i<=size(outaloc,1)) then
outaloc(i) = inchpcoul(i)+inchpdft(i)
outaloc(i+size(outaloc,1)/2) = inchpcoul(i)+inchpdft(i+size(outaloc,1)/2)
end if

end subroutine cugetaloc

attributes(global) subroutine cugetalocsubstractpotion(outaloc,inchpcoul,inchpdft,inpotion)
	implicit none 
	REAL(DP)		:: outaloc(:)
	REAL(DP)		:: inchpcoul(:)
	REAL(DP)		:: inchpdft(:)
	REAL(DP)		:: inpotion(:)

	integer  :: i


i = (blockIdx%x-1)*blockDim%x+threadIdx%x

if(i<=size(outaloc,1)) then
outaloc(i) = inchpcoul(i)+inchpdft(i)-inpotion(i)
outaloc(i+size(outaloc,1)/2) = inchpcoul(i)+inchpdft(i+size(outaloc,1)/2)-inpotion(i)
end if

end subroutine cugetalocsubstractpotion

attributes(global) subroutine cucalc_lda(inrho,outchpdft,dvol,energout)
! Pade' approximant to Perdew-Wang 92 density functional
	implicit none 
	REAL(DP)		:: inrho(:)
	REAL(DP)		:: outchpdft(:)
	REAL(DP)	,value  :: dvol
	REAL(DP)		:: energout
	REAL(DP) , shared	:: temp(blockDim%x)

	INTEGER :: ii,j,inext
	REAL(DP) :: ec, rp, xi
	REAL(DP) :: a0, a1, a2, a3, da0, da1, da2, da3
	REAL(DP) :: b1, b2, b3, b4, db1, db2, db3, db4
	REAl(DP) :: t, t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t12, t13, t15, t17, &
		  & t22, t23, t24, t25, t26, t28, t29, t34, t35, t36, t37, t42, t44, t48, &
		  & t53, t58, t63, t64, t65, t68, t70, t71, t72, t77, t82,  t83, t88, t93, t98, t102, t109, t135

DATA a0  /0.458165293D0 /
DATA da0 /0.119086804D0 /

DATA a1  /2.2170586D0 /
DATA da1 /0.615740256D0 /

DATA a2  /0.740555173D0 /
DATA da2 /0.157420151D0 /

DATA a3  /0.019682278D0 /
DATA da3 /0.003532336D0 /

DATA b1  /1.000000000D0 /
DATA db1 /0.000000000D0 /

DATA b2  /4.504130959D0 /
DATA db2 /0.2361297D0 /

DATA b3  /1.1106363D0 /
DATA db3 /0.205200460D0 /

DATA b4  /0.023592917D0/
DATA db4 /0.004200005D0 /
ii = (blockIdx%x-1)*blockDim%x+threadIdx%x
j=ii
temp(j) = 0.D0
do while (ii<=size(outchpdft,1)/2)
  rp     = MAX(inrho(ii),1D-16)
  xi     = inrho(ii+size(inrho,1)/2)
  t1 = xi*rp
  t2 = rp
  t3 = 1D0/t2
  t4 = xi
  t6 = (1D0+t4)**(1D0/3D0)
  t7 = t6**2
  t8 = t7**2
  t10 = (1D0-t4)**(1D0/3D0)
  t11 = t10**2
  t12 = t11**2
  t13 = t8+t12-2
  t15 = 2D0**(1D0/3D0)
  t17 = 1D0/(2D0*t15-2D0)
  t22 = 3D0**(1D0/3D0)
  t23 = (a1+da1*t13*t17)*t22
  t24 = 4D0**(1D0/3D0)
  t25 = t24**2
  t26 = 1/3.141592653589793D0
  t28 = (t26*t3)**(1D0/3D0)
  t29 = t25*t28
  t34 = t22**2
  t35 = (a2+da2*t13*t17)*t34
  t36 = t28**2
  t37 = t24*t36
  t42 = (a3+da3*t13*t17)*t26
  t44 = a0+da0*t13*t17+t23*t29/4.0D0+t35*t37/4D0+3D0/4D0*t42*t3
  t48 = (b1+db1*t13*t17)*t22
  t53 = (b2+db2*t13*t17)*t34
  t58 = (b3+db3*t13*t17)*t26
  t63 = (b4+db4*t13*t17)*t22
  t64 = t36**2
  t = t48*t29/4D0+t53*t37/4D0+3D0/4D0*t58*t3+3D0/16D0*t63*t25*t64
  t68 = 1D0/t
  t70 = t2**2
  t71 = 1D0/t70
  t72 = t1*t71
  t77 = 4D0/3D0*t6*(t3-t72)+4D0/3D0*t10*(-t3+t72)
  t82 = t22*t25
  t83 = t82*t28
  t88 = 1D0/t36*t26*t71
  t93 = t34*t24*t36
  t98 = 1D0/t28*t26*t71
  t102 = t17*t26*t3
  t109 = t**2
  t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4D0-t23*t25*t88/12D0+da2  &
      *t77*t17*t93/4D0-t35*t24*t98/6D0+3D0/4D0*da3*t77*t102-3D0/4D0*t42  &
      *t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4-t48*t25*t88/12D0+db2*t77*t17  &
      *t93/4D0-t53*t24*t98/6D0+3D0/4D0*db3*t77*t102-3D0/4D0*t58*t71+3D0  &
      /16D0*db4*t77*t17*t82*t64-t63*t25*t28*t26*t71/4D0)

	outchpdft(ii)      = -t135  * 2.D0!e2

  t77 = 4D0/3D0*t6*(-t3-t72)+4D0/3D0*t10*(t3+t72)
  t82 = t22*t25
  t83 = t82*t28
  t88 = 1D0/t36*t26*t71
  t93 = t34*t24*t36
  t98 = 1D0/t28*t26*t71
  t102 = t17*t26*t3
  t109 = t**2
  t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4D0-t23*t25*t88/12+da2  &
      *t77*t17*t93/4D0-t35*t24*t98/6D0+3D0/4D0*da3*t77*t102-3D0/4D0*t42  &
      *t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4D0-t48*t25*t88/12D0+db2*t77*t17 &
      *t93/4D0-t53*t24*t98/6D0+3D0/4D0*db3*t77*t102-3D0/4D0*t58*t71+3D0  &
      /16D0*db4*t77*t17*t82*t64-t63*t25*t28*t26*t71/4D0)

       outchpdft(ii+size(outchpdft,1)/2) = -t135  * 2.D0!e2

  
  
  
  t1=rp
  t4 = xi
  t6 = (1D0+t4)**(1D0/3D0)
  t7 = t6**2
  t8 = t7**2
  t10 = (1D0-t4)**(1D0/3D0)
  t11 = t10**2
  t12 = t11**2
  t13 = t8+t12-2D0
  t15 = 2D0**(1D0/3D0)
  t17 = 1D0/(2D0*t15-2D0)
  t22 = 3D0**(1D0/3D0)
  t24 = 4D0**(1D0/3D0)
  t25 = t24**2
  t26 = 1D0/3.1415926
  t28 = (t26*t3)**(1D0/3D0)
  t29 = t25*t28
  t34 = t22**2
  t36 = t28**2
  t37 = t24*t36
  t65 = t36**2
  t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*t22*t29/4D0+(a2+da2*t13*t17  &
      )*t34*t37/4D0+3D0/4D0*(a3+da3*t13*t17)*t26*t3)/((b1+db1*t13*t17)*  &
      t22*t29/4D0+(b2+db2*t13*t17)*t34*t37/4D0+3D0/4D0*(b3+db3*t13*t17)*t26 &
      *t3+3D0/16D0*(b4+db4*t13*t17)*t22*t25*t65)
  
!    xc energy-density is now:   -e2*t70/rp
  
!    next step is to compose rearrangement energy
  
  t1 = xi*rp
  t2 = rp
  t3=ABS((t1+t2)/2D0)   !  *e2
  t4=ABS((t1-t2)/2D0)   !  *e2
t5= outchpdft(ii)*t3+outchpdft(ii+size(outchpdft,1)/2)*t4
       temp(j) = (-t70*2.D0 - 0.5D0*t5) + temp(j)
ii = ii + blockDim%x
end do

call syncthreads()
inext = blockDim%x/2

	do while (inext>= 1)
		if (j <= inext) &
		temp(j) = temp(j) + temp(j+inext)
	inext = inext/2
	call syncthreads()
	end do
if (j==1) energout=temp(1)*dvol

end subroutine cucalc_lda

attributes(global) subroutine curesult(inrhoinp,outrhokr,e2,nx,ny,nz)
	implicit none
	REAL(DP)		:: inrhoinp(:)
	REAL(DP)		:: outrhokr(:)
	REAL(DP) ,value		:: e2
	INTEGER  ,value		:: nx,ny,nz

	integer  :: i,j,i1, i2, i3
	
	i1  = (blockIdx%x-1)*blockDim%x + threadIdx%x
	i2  = (blockIdx%y-1)*blockDim%y + threadIdx%y
	i3  = (blockIdx%z-1)*blockDim%z + threadIdx%z

      IF(i3 <= nz .AND. i2 <= ny .AND. i1 <= nx) THEN
	i   = 4*nx*ny*(i3-1)+2*nx*(i2-1)+i1 

	j = 4*nx*ny*(min(i3,nz)-1)/4+2*nx*(min(i2,ny)-1)/2+min(i1,nx)
        inrhoinp(j) = e2*outrhokr(i)
      END IF

end subroutine curesult

attributes(global) subroutine cumultarrays(inrhokr,inrhoki,mult)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	real(DP) ,value	:: mult

	integer :: i 
	i = (blockIdx%x-1)*blockDim%x + threadIdx%x

	if (i<=size(inrhokr,1)) then
		inrhokr(i) = mult*inrhokr(i) 
		inrhoki(i) = mult*inrhoki(i)
	end if
end subroutine cumultarrays


attributes(global) subroutine cunettocharg(outrhon,inelecrho,inionrho)
        implicit none
        real(DP)        :: outrhon(:)
        real(DP)        :: inelecrho(:)
        real(DP)        :: inionrho(:)

        integer :: i
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (i<=size(outrhon,1)) then
                outrhon(i) = inelecrho(i)-inionrho(i)
        end if
end subroutine cunettocharg

attributes(global) subroutine cuaffectcharg(outrhon,inelecrho)
        implicit none
        real(DP)        :: outrhon(:)
        real(DP)        :: inelecrho(:)

        integer :: i
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (i<=size(outrhon,1)) then
                outrhon(i) = inelecrho(i)
        end if
end subroutine cuaffectcharg

attributes(global) subroutine cusubstractpotion(outchpcoul,inpotion)
        implicit none
        real(DP)        :: outchpcoul(:)
        real(DP)        :: inpotion(:)

        integer :: i
        i = (blockIdx%x-1)*blockDim%x + threadIdx%x

        if (i<=size(outchpcoul,1)) then
                outchpcoul(i) = outchpcoul(i)-inpotion(i)
        end if
end subroutine cusubstractpotion

attributes(global) subroutine cukineticfactor(inrhokr,inrhoki,kkr,kki)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	real(DP) 	:: kkr(:)
	real(DP)	:: kki(:)

	real(DP) :: save2
	integer :: i 	
	i = (blockIdx%x-1)*blockDim%x + threadIdx%x

	if (i<=size(inrhokr,1)) then
	 save2     = kkr(i)*inrhokr(i)+kki(i)*inrhoki(i)
	  inrhoki(i) = kkr(i)*inrhoki(i)+kki(i)*inrhokr(i)
	  inrhokr(i) = save2
	end if
end subroutine cukineticfactor


attributes(global) subroutine compXYZgrad(inpsi,outfftax,outfftay,outfftaz)
	implicit none
	complex(DP)	:: inpsi(:)
	complex(DP)	:: outfftax(:)
	complex(DP)	:: outfftay(:)
	complex(DP)	:: outfftaz(:)

	integer :: i1,i2,i3,ind,off

	i1 = threadIdx%x
	i2 = blockIdx%y
	i3 = blockIdx%x

	ind = (i3-1)*blockDim%x*gridDim%y+(i2-1)*blockDim%x+i1
	off = (i3-1)*blockDim%x*gridDim%y+(i2-1)*blockDim%x
	outfftax(mod(i1+blockDim%x/2,blockDim%x)+1+off)   = inpsi(ind)
	off = (i3-1)*gridDim%y*blockDim%x+(i1-1)*gridDim%y
	outfftay(mod(i2+gridDim%y/2,gridDim%y)+1  +off)   = inpsi(ind)
	off = (i2-1)*gridDim%x*blockDim%x+(i1-1)*gridDim%x
	outfftaz(mod(i3+gridDim%x/2,gridDim%x)+1  +off)   = inpsi(ind)
end subroutine compXYZgrad

attributes(global) subroutine comp3grad(inpsi,outfftax,outfftay,outfftaz)
	implicit none
	complex(DP)	:: inpsi(:,:)
	complex(DP)	:: outfftax(:)
	complex(DP)	:: outfftay(:)
	complex(DP)	:: outfftaz(:)

	integer :: i1,i2,i3,ind,off

	i1 = threadIdx%x
	i2 = blockIdx%y
	i3 = blockIdx%x

	ind = (i3-1)*blockDim%x*gridDim%y+(i2-1)*blockDim%x+i1
	off = (i3-1)*blockDim%x*gridDim%y+(i2-1)*blockDim%x
	outfftax(mod(i1+blockDim%x/2,blockDim%x)+1+off)   = inpsi(ind,1)
	off = (i3-1)*gridDim%y*blockDim%x+(i1-1)*gridDim%y
	outfftay(mod(i2+gridDim%y/2,gridDim%y)+1  +off)   = inpsi(ind,2)
	off = (i2-1)*gridDim%x*blockDim%x+(i1-1)*gridDim%x
	outfftaz(mod(i3+gridDim%x/2,gridDim%x)+1  +off)   = inpsi(ind,3)
end subroutine comp3grad


attributes(global) subroutine cuapplygrad(infftax,infftay,infftaz,dkx,dky,dkz)
	implicit none
	complex(DP)	:: infftax(:)
	complex(DP)	:: infftay(:)
	complex(DP)	:: infftaz(:)
	real(DP) ,value :: dkx,dky,dkz

	integer :: i1,i2,i3,off

	i1 = threadIdx%x
	i2 = blockIdx%y
	i3 = blockIdx%x

	off = (i3-1)*blockDim%x*gridDim%y+(i2-1)*blockDim%x
	IF(i1 == (blockDim%x/2+1)) THEN
		infftax(i1+off) = (0D0,0D0)
	ELSE IF(i1 > (blockDim%x/2+1)) THEN
		infftax(i1+off) = infftax(i1+off)*cmplx(0D0,(i1-blockDim%x-1)*dkx)
	ELSE
		infftax(i1+off) = infftax(i1+off)*cmplx(0D0,(i1-1)*dkx)
	END IF

	off = (i3-1)*gridDim%y*blockDim%x+(i1-1)*gridDim%y
	IF(i2 == (gridDim%y/2+1)) THEN
		infftay(i2+off) = (0D0,0D0)
	ELSE IF(i2 > (gridDim%y/2+1)) THEN
		infftay(i2+off) = infftay(i2+off)*cmplx(0D0,(i2-gridDim%y-1)*dky)
	ELSE
		infftay(i2+off) = infftay(i2+off)*cmplx(0D0,(i2-1)*dky)
	END IF

	off = (i2-1)*gridDim%x*blockDim%x+(i1-1)*gridDim%x
	IF(i3 == (gridDim%x/2+1)) THEN
		infftaz(i3+off) = (0D0,0D0)
	ELSE IF(i3 > (gridDim%x/2+1)) THEN
		infftaz(i3+off) = infftaz(i3+off)*cmplx(0D0,(i3-gridDim%x-1)*dkz)
	ELSE
		infftaz(i3+off) = infftaz(i3+off)*cmplx(0D0,(i3-1)*dkz)
	END IF
end subroutine cuapplygrad

attributes(global) subroutine decompXYZgrad(outgrad,infftax,infftay,infftaz)
	implicit none
	complex(DP)	:: outgrad(:,:)
	complex(DP)	:: infftax(:)
	complex(DP)	:: infftay(:)
	complex(DP)	:: infftaz(:)

	integer :: i1,i2,i3,ind,off

	i1 = threadIdx%x
	i2 = blockIdx%y
	i3 = blockIdx%x

	ind = (i3-1)*blockDim%x*gridDim%y+(i2-1)*blockDim%x+i1
	off = (i3-1)*blockDim%x*gridDim%y+(i2-1)*blockDim%x
	outgrad(ind,1) = infftax(mod(i1+blockDim%x/2,blockDim%x)+1+off)/blockDim%x	
	off = (i3-1)*gridDim%y*blockDim%x+(i1-1)*gridDim%y
	outgrad(ind,2) = infftay(mod(i2+gridDim%y/2,gridDim%y)+1  +off)/gridDim%y
	off = (i2-1)*gridDim%x*blockDim%x+(i1-1)*gridDim%x
	outgrad(ind,3) = infftaz(mod(i3+gridDim%x/2,gridDim%x)+1  +off)/gridDim%x
end subroutine decompXYZgrad

attributes(global) subroutine cuaddgradtocurr(outcurr,ingrad,inpsi,occupin)
	implicit none
	real(DP)	:: outcurr(:,:)
	complex(DP)	:: ingrad(:,:)
	complex(DP)	:: inpsi(:)
	real(DP)	:: occupin

	integer :: i,j

	i = (blockIdx%x-1)*blockDim%x+threadIdx%x
	j = (blockIdx%y-1)*blockDim%y+threadIdx%y
	if (i<=size(outcurr,1) .and. j<=size(outcurr,2)) then
	outcurr(i,j) = outcurr(i,j) + occupin*AIMAG(CONJG(inpsi(i))*ingrad(i,j))
	end if
end subroutine cuaddgradtocurr

attributes(global) subroutine compXFFT(inrhokr,outfftax,nx,nxy1)
	implicit none
	real(DP) 	:: inrhokr(:)
	complex(DP)	:: outfftax(:)
	integer , value :: nx,nxy1

        integer          :: i0,i1
	integer 	:: i ,j

	i1  = ((blockIdx%x-1)*(gridDim%y)+(blockIdx%y-1))*2*nx
        i0 = (blockIdx%x-1)*nxy1 + (blockIdx%y-1)*2*nx
	i = threadIdx%x
	if (i<=size(outfftax,1)) then
		if (i<=nx+1) then
			j = i0 + nx -1+i
		else
			j = i0 - nx -1+i ! a verif
		end if
		outfftax(i1+i) = CMPLX(inrhokr(j),0D0,DP)
	end if
end subroutine compXFFT

attributes(global) subroutine compXBFT(inrhokr,inrhoki,outfftax,nx,nxy1)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftax(:)
	integer , value :: nx,nxy1

        integer          :: i0,i1
        integer         :: i ,j

        i1  = ((blockIdx%x-1)*(gridDim%y)+(blockIdx%y-1))*2*nx
        i0 = (blockIdx%x-1)*nxy1 + (blockIdx%y-1)*2*nx
        i = threadIdx%x

	if (i<=size(outfftax,1)) then
		if (i<=nx+1) then
			j = i0 + nx -1+i
		outfftax(i+i1) = CMPLX(inrhokr(j),inrhoki(j),DP)
		else
			j = i0 +3*nx-i+1 ! a verif
		outfftax(i+i1) = CMPLX(inrhokr(j),-inrhoki(j),DP)
		end if
	end if
end subroutine compXBFT


attributes(global) subroutine compYFFT(inrhokr,inrhoki,outfftay,ny,nxi,nxy1)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftay(:)
	integer , value :: ny,nxi,nxy1

	integer 	:: i ,j,i1,i0

	i  = threadIdx%x
	i1  = ((blockIdx%x-1)*(gridDim%y)+(blockIdx%y-1))*2*ny
	i0 = (blockIdx%x-1)*nxy1 + (blockIdx%y-1+nxi/2)
	if (i<=size(outfftay,1)) then
		if (i<=ny+1) then
			j = i0+nxi*(ny-2)+i*nxi
		else
			j = i0 - nxi +(i-ny-1)*nxi
		end if
		outfftay(i+i1) = CMPLX(inrhokr(j),inrhoki(j),DP)
	end if
end subroutine compYFFT

attributes(global) subroutine compYBFT(inrhokr,inrhoki,outfftay,ny,nxi,nxy1)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftay(:)
	integer , value :: ny,nxi,nxy1

	integer 	:: i ,j,i1,i0

	i  = threadIdx%x
	i1  = ((blockIdx%x-1)*(gridDim%y)+(blockIdx%y-1))*2*ny
	i0 = (blockIdx%x-1)*nxy1 + (blockIdx%y-1+nxi/2)
	if (i<=size(outfftay,1)) then
		if (i<=ny+1) then
			j = i0+nxi*(ny-2)+i*nxi
		else
			j = i0  +(i-ny-2)*nxi
		end if
		outfftay(i+i1) = CMPLX(inrhokr(j),inrhoki(j),DP)
	end if
end subroutine compYBFT


attributes(global) subroutine compZFFT(inrhokr,inrhoki,outfftaz,nxyf,nyf,nzh,kfftz)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftaz(:)
	integer , value :: nxyf,nyf,nzh,kfftz

        integer         :: i1,i2,offset
	integer 	:: i ,j
	
        i1 = blockIdx%x+gridDim%x-2
        i2 = blockIdx%y
        offset = ((i2-1)*gridDim%x+(i1-gridDim%x+1))*blockDim%x
	i = MOD((threadIdx%x)+nzh,kfftz)+1
	if (i+offset<=size(outfftaz,1)) then
	!      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
		j = (threadIdx%x-1)*nxyf+(i2-1)*nyf+i1
		outfftaz(i+offset) = CMPLX(inrhokr(j),inrhoki(j),DP)
	end if
end subroutine compZFFT


attributes(global) subroutine compZBFT(inrhokr,inrhoki,outfftaz,nxyf,nyf,nzh,kfftz)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftaz(:)
	integer , value :: nxyf,nyf,nzh,kfftz

        integer         ::i1,i2,offset
	integer 	:: i ,j
	i1 = blockIdx%x+gridDim%x-2
        i2 = blockIdx%y
        !offset = ((i2-1)*gridDim%y+i2)*blockDim%x
	offset = ((i2-1)*gridDim%x+(i1-1))*blockDim%x
	i = MOD(( threadIdx%x)+nzh,kfftz)+1
	if (i<=size(outfftaz,1)) then
	!      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
		j = (threadIdx%x-1)*nxyf+(i2-1)*nyf+i1
		outfftaz(i+offset) = CMPLX(inrhokr(j),inrhoki(j),DP)
	end if
end subroutine compZBFT

attributes(global) subroutine decompXFFT(inrhokr,inrhoki,outfftax,nx,nxy1)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftax(:)
	integer , value :: nx,nxy1

        integer          :: i0,i1
	integer 	:: i ,j

	i = threadIdx%x
	i1  = ((blockIdx%x-1)*(gridDim%y)+(blockIdx%y-1))*2*nx
        i0 = (blockIdx%x-1)*nxy1 + (blockIdx%y-1)*2*nx

	if (i<=size(outfftax,1)) then
		if (i<=nx+1) then
			j = i0 + nx -1+i
		else
			j = i0 - nx -1+i
		end if
		inrhokr(j) = REAL(outfftax(i+i1),DP)
		inrhoki(j) = aimag(outfftax(i+i1)) 
	end if
end subroutine decompXFFT

attributes(global) subroutine decompXBFT(inrhokr,inrhoki,outfftax,nx,i0)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftax(:)
	integer , value :: nx,i0

	integer 	:: i ,j

	i = (blockIdx%x-1)*blockDim%x + threadIdx%x

	if (i<=size(outfftax,1)) then
		if (i<=nx+1) then
			j = i0 + nx -1+i
		else
			j = i0 - nx -1+i
		end if
		inrhokr(j) = REAL(outfftax(i),DP)
		inrhoki(j) = aimag(outfftax(i)) 
	end if
end subroutine decompXBFT

attributes(global) subroutine decompYFFT(inrhokr,inrhoki,outfftay,ny,nxi,nxy1)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP) 	:: inrhoki(:)
	complex(DP)	:: outfftay(:)
	integer , value :: ny,nxi,nxy1

	integer 	:: i ,j,i1,i0

	i  = threadIdx%x
	i1  = ((blockIdx%x-1)*(gridDim%y)+(blockIdx%y-1))*2*ny
	i0 = (blockIdx%x-1)*nxy1 + (blockIdx%y-1+nxi/2)
	if (i<=size(outfftay,1)) then
		if (i<=ny+1) then
			j = i0+nxi*(ny-2)+i*nxi
		else
			j = i0 - nxi +(i-ny-1)*nxi
		end if
		inrhokr(j) = REAL(outfftay(i+i1),DP)
		inrhoki(j) = aimag(outfftay(i+i1)) 
	end if
end subroutine decompYFFT

attributes(global) subroutine decompYBFT(inrhokr,inrhoki,outfftay,ny,nxi,nxy1)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP) 	:: inrhoki(:)
	complex(DP)	:: outfftay(:)
	integer , value :: ny,nxi,nxy1

	integer 	:: i ,j,i1,i0

	i  = threadIdx%x
	i1  = ((blockIdx%x-1)*(gridDim%y)+(blockIdx%y-1))*2*ny
	i0 = (blockIdx%x-1)*nxy1 + (blockIdx%y-1+nxi/2)
	if (i<=size(outfftay,1)) then
		if (i<=ny+1) then
			j = i0+nxi*(ny-2)+i*nxi
		else
			j = i0  +(i-ny-2)*nxi
		end if
		inrhokr(j) = REAL(outfftay(i+i1),DP)
		inrhoki(j) = aimag(outfftay(i+i1)) 
	end if
end subroutine decompYBFT

attributes(global) subroutine decompZFFT(inrhokr,inrhoki,outfftaz,nxyf,nyf,nzh,kfftz)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftaz(:)
	integer , value :: nxyf,nyf,nzh,kfftz

        integer         :: i1,i2,offset
	integer 	:: i ,j
        
        i1 = blockIdx%x+gridDim%x-2
        i2 = blockIdx%y
        offset = ((i2-1)*gridDim%x+(i1-gridDim%x+1))*blockDim%x
	i = MOD(( threadIdx%x)+nzh,kfftz)+1
	if (i+offset<=size(outfftaz,1)) then
		j = ( threadIdx%x-1)*nxyf+(i2-1)*nyf+i1
		inrhokr(j) = REAL(outfftaz(i+offset),DP)
		inrhoki(j) = aimag(outfftaz(i+offset)) 
	end if
end subroutine decompZFFT

attributes(global) subroutine decompZBFT(inrhokr,inrhoki,outfftaz,nxyf,nyf,nzh,kfftz)
	implicit none
	real(DP) 	:: inrhokr(:)
	real(DP)	:: inrhoki(:)
	complex(DP)	:: outfftaz(:)
	integer , value :: nxyf,nyf,nzh,kfftz

        integer         :: i1,i2
	integer 	:: i ,j

        i1= blockIdx%x
        i2= blockIdx%y
	i = MOD((threadIdx%x)+nzh,kfftz)+1
	if (i<=size(outfftaz,1)) then
		j = ( threadIdx%x-1)*nxyf+(i2-1)*nyf+i1
		inrhokr(j) = REAL(outfftaz(i),DP)
		inrhoki(j) = aimag(outfftaz(i)) 
	end if
end subroutine decompZBFT


END MODULE cutools
