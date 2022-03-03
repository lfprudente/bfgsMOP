subroutine BFGSupdate(n,m,x,xprev,JF,JFprev,strconvex,theta,B)

implicit none
	
	! SCALAR ARGUMENTS 
	integer,intent(in) :: n,m
	real(kind=8), intent(in) :: theta
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: x(n),xprev(n),JF(m,n),JFprev(m,n)
	real(kind=8), intent(inout) :: B(n,n,m)
	logical,intent(in) :: strconvex(m)
	
	! LOCAL SCALARS
	integer :: j,ind
	real(kind=8) :: rho,Dxs,sTy,sTBs,denominator,eps
	
	! LOCAL ARRAYS
	real(kind=8) :: s(n),y(n),Bs(n),BssTB(n,n),yyT(n,n),ysTB(n,n),BsyT(n,n)
	
	eps = 1.0d-6
	
	eps = eps * min( abs(theta), 1.0d0 )
	
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
		
		if ( strconvex(ind) .or. sTy > eps ) then
			B(:,:,ind) = B(:,:,ind) - BssTB / sTBs + yyT / sTy
		else

			! Define rho (actually, rho^(-1))
			
			rho = Dxs - dot_product( JFprev(ind,:),s ) 
			
			denominator = ( rho - sTy ) ** 2 + rho * sTBs
			
			do j = 1,n
				ysTB(:,j) = Bs(j) * y
			end do
			
			BsyT = transpose(ysTB)
			
			B(:,:,ind) = B(:,:,ind) - rho * BssTB / denominator &
				         + sTBs * yyT / denominator &
			             + ( rho - sTy ) * ( ysTB + BsyT ) / denominator
		end if

	end do	

end subroutine BFGSupdate
