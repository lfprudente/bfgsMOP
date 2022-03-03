subroutine scalefactor(n,m,x,scaleF)	

	use globals, only: sF
	use myproblem
	
	implicit none
	
	! SCALAR ARGUMENTS
	integer, intent(in)  :: n,m
	logical, intent(in)  :: scaleF

	! ARRAY ARGUMENTS
	real(kind=8), intent(in)  :: x(n)
		
	! LOCAL SCALAR
	integer :: i,ind
	real(kind=8) :: eps
	
	! LOCAL ARRAY
	real(kind=8) :: g(n)
	
	eps = 1.0d-08

! Compute objective functions scaling factor

	if ( scaleF ) then
		do ind = 1,m
			call evalg(n,x,g,ind)

			sF(ind) = 1.0d0
			do i = 1,n
				sF(ind) = max( sF(ind), abs( g(i) ) )
			end do
			
			sF(ind) = max( eps, 1.0d0 / sF(ind) )
		end do
	else
		sF(1:m) = 1.0d0
	end if
	
end subroutine scalefactor
