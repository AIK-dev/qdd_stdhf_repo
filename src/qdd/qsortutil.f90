  SUBROUTINE qsort(nb,array,indarray,reverse)

	integer, intent(in) :: nb
	real(8), intent(in out) :: array(nb)
	integer, intent(out) :: indarray(nb)
	logical, intent(in) :: reverse
	integer :: iki

	do iki = 1, nb
	    indarray(iki)=iki
	end do
	! write(*,*) ''
	! write(*,*) 'QSORT'
	! write(*,'(2x,50f5.2)') (array(:))
	call quicksort(nb,array,indarray,reverse)
	! write(*,'(2x,50f5.2)') array(:)
	! write(*,'(2x,50I5)') indarray(:)
	! write(*,*) ''

contains
recursive subroutine quicksort(nb,array,indarray,reverse)

	integer, intent(in) :: nb
	real(8), intent(in out) :: array(nb)
	integer, intent(in out) :: indarray(nb)
	logical, intent(in) :: reverse
	integer :: p
    write(*,*)'nb',nb
	if (nb>0) then
	call partition(nb,array,indarray,p,reverse)
	call quicksort((p-1),array(:(p-1)),indarray(:(p-1)),reverse)
	call quicksort((nb-p),array((p+1):nb), &
	    indarray((p+1):nb),reverse)
	end if

end subroutine quicksort

subroutine partition(nb,array,indarray,p,reverse)

	integer, intent(in) :: nb
	real(8), intent(in out) :: array(nb)
	integer, intent(in out) :: indarray(nb)
	integer, intent(out) :: p
	logical, intent(in) :: reverse
	integer :: iki
	real(8) :: pval

	p=(nb+1)/2
	pval=array(p)

	call swap(nb,array,indarray,p,nb)

	p=1

	if (reverse) then
	do iki = 1, nb-1
	    if (array(iki)>=pval) then
		call swap(nb,array,indarray,iki,p)
		p=p+1
	    end if
	end do
	else
	do iki = 1, nb-1
	    if (array(iki)<=pval) then
		call swap(nb,array,indarray,iki,p)
		p=p+1
	    end if
	end do
	end if

	call swap(nb,array,indarray,p,nb)

end subroutine partition

subroutine swap(nb,array,indarray,i1,i2)

	integer, intent(in) :: nb
	real(8), intent(in out) :: array(nb)
	integer, intent(in out) :: indarray(nb)
	integer, intent(in) :: i1,i2
	integer :: a
	real(8) :: x

	x=array(i2)
	a=indarray(i2)
	array(i2)=array(i1)
	indarray(i2)=indarray(i1)
	array(i1)=x
	indarray(i1)=a

end subroutine swap
end subroutine qsort