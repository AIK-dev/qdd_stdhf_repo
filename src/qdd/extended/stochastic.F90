MODULE sto
use params
contains
subroutine stochastic(psi,rho,aloc,eout,it)

implicit none
	COMPLEX(DP), INTENT(IN OUT)         :: psi(kdfull2,kstate)
	INTEGER    , INTENT(IN)             :: it
	COMPLEX(DP), INTENT(IN OUT)         :: Eout(kstate,kstate)

	REAL(DP) , INTENT(IN OUT)           :: aloc(2*kdfull2)
	REAL(DP) , INTENT(IN OUT)           :: rho(2*kdfull2)
	REAL(DP)                :: energref,spenergref
	COMPLEX(DP)::MoP(kstate,kstate)
	!--------------------------
	INTEGER       :: i,nb,ordering(kstate)
	REAL(DP)      :: sumEout
   !energy has just been computed, use params
  WRITE(*,*) 'START STOCHASTIC'
  energref=energy;
  spenergref=0.d0
  do i=1,nstate; 
  	energref=energref+eout(i,i)*occup(i);
  	spenergref=spenergref+(amoy(i)+eout(i,i))*occup(i)
  enddo
  write(*,*)'enegref, spenergeref', energref,spenergref
  !reorder, the occupied first.
  write(*,*) 'before hdiag'
  call hdiag_part_hole(psi,aloc,MoP,Eout)

end subroutine stochastic
SUBROUTINE hdiag_part_hole(psi,aloc,actMoP,Eout)!mv
	USE params 
	USE util, ONLY: wfovlp
	IMPLICIT NONE
	!--------------------------
	COMPLEX(DP), INTENT(IN OUT)         :: psi(kdfull2,kstate)
	REAL(DP), INTENT(IN)                :: aloc(2*kdfull2)
	COMPLEX(DP), INTENT(IN OUT)         :: actMoP(kstate,kstate)
	COMPLEX(DP), INTENT(IN)             :: Eout(kstate,kstate)
	!--------------------------
	COMPLEX(DP),allocatable   :: hmatr(:,:),vect(:,:)
	COMPLEX(DP),allocatable   :: eoutact(:,:),psiaux(:)
	 REAL(DP),allocatable      :: eigen(:)
	 !--------------------------
	INTEGER       :: nb , st
	INTEGER       :: nelecup,nelecdown
	allocate(hmatr(kstate,kstate),vect(kstate,kstate),eoutact(kstate,kstate),eigen(kstate))
	write(*,*) 'before sortmop'
	call sort_MoP(psi,.true.,actMoP)
	 eoutact=matmul(transpose(actMoP),Eout)
     eoutact=matmul(eoutact,actMoP)
	 hmatr(:,:) = cmplx(0.d0,0.d0,8)
	 vect(:,:) = cmplx(0.d0,0.d0,8)
	 allocate(psiaux(kdfull2))
	 do nb = 1, nstate
	 	psiaux=psi(:,nb)
	   	call hpsi(psiaux,aloc,nb,0)
	 	do st = 1, nstate
	 		if(ispin(nb).eq.ispin(st)) then 
       			hmatr(nb,st) = wfovlp(psi(:,st),psiaux)
     		endif
	 	end do
	 	hmatr(nb,nb)=hmatr(nb,nb)+eoutact(nb,nb)
	 end do
	 deallocate (psiaux)
	 write(*,'(20f5.2)') ((abs(hmatr(nb,st)),nb=1,nstate),st=1,nstate)

	! call cdiag(hmatr(1:nelect,1:nelect)&
	!           ,eigen(1:nelect)&
    !           ,vect(1:nelect,1:nelect)&
    !           ,nelect,nelect)
    !       write(*,*) 'eigen 1'
    ! write(*,'(10 f8.4)')eigen

    ! call cdiag(hmatr(nelect+1:nstate,nelect+1:nstate)&
	!           ,eigen(nelect+1:nstate)&
    !           ,vect(nelect+1:nstate,nelect+1:nstate)&
    !           ,nstate-nelect,nstate-nelect)
	
  
    
    !   write(*,*) 'vect'
	! write(*,'(20f5.2)') ((abs(vect(nb,st)),nb=1,nstate),st=1,nstate)
    ! call calcispin(vect)
    ! write(*,*) 'eigen'
    ! write(*,'(10 f8.4)')eigen

  	! write(*,*) 'ispindiag'
    ! write(*,'(10 f8.4)')ispin
  	! psi =matmul(psi,vect)
    ! actMoP=matmul(actmop,vect)
end SUBROUTINE hdiag_part_hole	

SUBROUTINE get_spenerg(psi,aloc,spnerg)

	USE params
	USE util, ONLY: wfovlp
#if(pgi)
	USE cudadyn
	USE cutools
#endif
	IMPLICIT NONE
	
	COMPLEX(DP), INTENT(IN OUT)         :: psi(kdfull2,kstate)
	REAL(DP), INTENT(IN)                :: aloc(2*kdfull2)


	REAL(DP), INTENT(OUT)               :: spnerg(kstate)


	COMPLEX(DP)  , allocatable  :: whpsi(:)
	INTEGER  :: nb
	
#if(pgi)
call cucalc_ekin()
call cucalc_epot<<<nstate,256,256*16>>>(psi_d,aloc_d,epot_d,dvol)
spnerg(:nstate) = epot_d(:nstate)
spnerg(:nstate) = spnerg(:nstate) + ekinsp(:nstate)
#else
	allocate(whpsi(kdfull2))
	  DO nb=1,nstate
	    whpsi = psi(:,nb)
	    CALL hpsi(whpsi,aloc,nb,0) ! with "0", "nb" is not usefull
	    spnerg(nb) = REAL(wfovlp(psi(:,nb),whpsi),DP)
	  END DO
     DEALLOCATE(whpsi)
#endif
RETURN

END SUBROUTINE get_spenerg
SUBROUTINE sort_MoP(psi,reverse,MoP)! psi, spnergy, occup
        USE params !occup, ispin are in params
        USE util
		use rta_module
        IMPLICIT NONE
        complex(8), intent(in out) :: psi(kdfull2,kstate)
        logical, intent(in) :: reverse
        complex(8), intent(in out) :: MoP(kstate,kstate)
		real(DP)::occloc(1:kstate)

        integer :: ordering(kstate)
        integer :: i,io=0,j
		MoP=cmplx(0.d0,0.d0,DP)
		write(*,*) 'sort Mop'
		   ! fill ordering with indices of ordered occup
		do i=1,nstate
		 	if (occup(i)>0.5) then
				io=io+1
				ordering(io)=i
			endif
		enddo
		do i=1,nstate
		 	if (occup(i)<0.5) then
				io=io+1
				ordering(io)=i
			endif
		enddo
		do i=1,nstate;MoP(ordering(i),i)=1.d0;enddo
		
		write(*,'(20i3)')ordering
		write(*,'(20f4.1)')(occup(ordering(i)), i=1,nstate)
		write(*,*)'mop'
		write(*,'(20f4.1)')((abs(mop(i,j)),j=1,nstate),i=1,nstate)
		
		

        !call qsort(nstate,occup(1:nstate),ordering(1:nstate),reverse)
		! put the other arrays in the right order
        psi=matmul(psi(:,1:nstate),Mop(1:nstate,1:nstate))
		write(*,*) 'psi done'
        call ireorder(ispin,ordering,nstate)
        call rreorder(occup(1:nstate),ordering(1:nstate))
write(*,*) 'fin sort mop'
END SUBROUTINE sort_MoP

subroutine calcispin(vect)
use params, only:ispin,DP,nstate,kstate,occup
complex(DP) ::vect(kstate,kstate),mat(kstate,kstate),mato(kstate,kstate)
real(DP):: ispindiag(kstate),occupdiag(kstate)
integer::i,j

mat=0.d0
mat0=0.d0
do i=1,nstate
    mat(i,i)=real(ispin(i),DP)
    mato(i,i)=real(occup(i),DP)
enddo
mat= matmul(transpose(conjg(vect)),matmul(mat ,vect))
mato=matmul(transpose(conjg(vect)),matmul(mato,vect))

do i=1,nstate
    ispindiag(i)=mat(i,i)
    occupdiag(i)=mato(i,i)
enddo
write(*,*) 'ispindiag'
    write(*,'(10 f8.4)')ispindiag
write(*,*) 'occupdiag'
    write(*,'(10 f8.4)')occupdiag
    occup=nint(occupdiag)
    ispin=nint(ispindiag)
end subroutine calcispin
end module sto