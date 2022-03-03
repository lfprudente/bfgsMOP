subroutine BFGSupdateCautions(n,m,x,xprev,JF,JFprev,theta,B)

implicit none
	
	! SCALAR ARGUMENTS 
	integer, intent(in) :: n,m
	real(kind=8), intent(in) :: theta
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: x(n),xprev(n),JF(m,n),JFprev(m,n)
	real(kind=8), intent(inout) :: B(n,n,m)
	
	! LOCAL SCALARS
	integer :: j,ind
	real(kind=8) :: rho,Dxs,sTy,sTBs,denominator,eps
	
	! LOCAL ARRAYS
	real(kind=8) :: s(n),y(n),Bs(n),BssTB(n,n),yyT(n,n),ysTB(n,n),BsyT(n,n)
	
	eps = 1.0d-6
	
	s = x - xprev
	
	! Compute f(x,s)
	
	Dxs = maxval( matmul( JF,s ) ) 
	
	do ind = 1,m
		y = JF(ind,:) - JFprev(ind,:)
		
		! Compute <s,y>
		
		sTy = dot_product( s, y )
		
		! Compute Bs
		
		Bs = matmul( B(:,:,ind), s )
		
		! Compute Bss^TB
		
		do j = 1,n
			BssTB(:,j) = Bs(j) * Bs
		end do
		
		! Compute yy^T
		
		do j = 1,n
			yyT(:,j) = y(j) * y
		end do
		
		! Compute sTBs
		
		sTBs = dot_product( s, Bs )
		
		! Update B
		
		if ( sTy >= eps * min( abs(theta), 1.0d0 ) ) then
			B(:,:,ind) = B(:,:,ind) - BssTB / sTBs + yyT / sTy
		end if
				
		
		
	end do	

end subroutine BFGSupdateCautions
