subroutine quadvar(quadratic,m,mquad,mnquad,quadind,nquadind)
	
	implicit none
	
	! SCALAR ARGUMENTS
	integer, intent(in)  :: m
	integer, intent(out) :: mquad,mnquad

	! ARRAY ARGUMENTS
	integer,allocatable,intent(out) :: quadind(:),nquadind(:)
	logical,intent(in) :: quadratic(m)
	
	! LOCAL SCALAR
	integer :: i,allocerr
	
	! Define variables mquad and mnquad	
	mquad  = 0
	mnquad = 0 
	do i = 1, m
		if ( quadratic(i) ) then
			mquad = mquad + 1
		else
			mnquad = mnquad + 1     
		end if
	end do
	
	allocate(quadind(mquad),nquadind(mnquad),stat=allocerr)
	if ( allocerr .ne. 0 ) then
		write(*,*) 'Allocation error in main program'
		stop
	end if	

	! Define variables quadind and nquadind	
	mquad  = 0
	mnquad = 0 
	do i = 1, m
		if ( quadratic(i) ) then
			mquad = mquad + 1
			quadind(mquad) = i
		else
			mnquad = mnquad + 1     
			nquadind(mnquad) = i        
		end if
	end do
	
end subroutine	

!***********************************************************************
!***********************************************************************

subroutine evalqstep(stp,ind,flag)

	use globals, only: n => nfinner, x => xinner, d
	use myproblem, only: quadstep

	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in) :: ind
	integer, intent(out) :: flag
	real(kind=8), intent(out) :: stp
	
	call quadstep(n,x,d,stp,ind,flag)

end subroutine evalqstep
