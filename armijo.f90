subroutine armijo(evalf,stp,stpmin,m,f0,g0,fend,ftol,iprint,nfev,info)
					
	implicit none
	
	! SCALAR ARGUMENTS 
	integer, intent(in) :: m
	integer, intent(out) :: nfev, info
	real(kind=8), intent(in) :: stpmin, ftol
	real(kind=8), intent(inout) :: stp
	logical,intent(in) :: iprint
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: f0(m), g0(m)
	real(kind=8), intent(out) :: fend(m)

!	
!   Given continuous differentiable functions f_i:R->R for i = 1,...,m,
!	this subroutine is designed to find a positive step-size that 
!	satisfies the sufficient decrease condition
!
!		f_i(stp) <= f_i(0) + ftol * stp * f'_max(0) for all i = 1,...,m, 
!
!   where f'_max(0) = max{f'_i(0) : i in {1,...m}}.
!
!	Mandatory external subroutine: evalf
!
!	On entry the user most provide the following variables: stp,m,f0,g0,ftol.
!
!	evalf is an external subroutine (see the interface block).
!		evalf(x,f,ind) should return the evalution of the function f_ind
!	    at x.
!
!	stp is a double precision variable.
!		On entry stp is the initial trial step-size. A positive initial 
!		estimate must be provided.
!		On exit stp is the final step-size. If info = 0 then stp satisfies 
!		the sufficient decrease condition. See 'info' 
!		below to see the others stopping criteria.
!
!	stpmin is a double precision variable.
!		On entry stpmin is a positive lower bound for the step.
!		On exit stpmin is unchanged.
!
!	m is an integer variable.
!		On entry m is the number of functions.
!		On exit m is unchanged.
!
!	f0(m) is a double precision array.
!		On initial entry f0(ind) is the value of the function f_ind at 0.
!		On exit f0 is unchanged.
!
!	g0(m) is a double precision array.
!		On initial entry g0(ind) is the value of the derivative of the 
!		function f_ind at 0.
!		On exit g0 is unchanged.
!
!	fend(m) is a double precision array.
!		On exit fend(ind) contains the value of f_ind(stp)
!
!	ftol is a double precision variable.
!		On entry ftol specifies a nonnegative parameter for the
!		sufficient decrease condition.
!		On exit ftol is unchanged.
!
!	iprint is a logical variable.
!		If iprint = .true. then the output is printed on the screen,
!		otherwise no output is shown.
!
!	nfev is a integer variable.
!		On exit nfev is the number of evaluations of the functions.
!
!	info is a integer variable.
!		On exit info contains the flag of the solution.
!	     0: stp satisfies the Armijo condition
!		 1: stp = stpmin   
!		-1: an error ocurred during the execution
!
	
	! LOCAL SCALARS
	integer :: i,iA,nfevbt
	real(kind=8), parameter :: zero = 0.0d0
	real(kind=8) :: f,maxg0,ftest,fiA
	logical :: sdc
	
	! INTERFACE OF EXTERNAL SUBROUTINES 
	interface
		! SUBROUTINE evalf:
		subroutine evalf(x,f,ind)
		! SCALAR ARGUMENTS
		integer, intent(in) :: ind
		real(kind=8), intent(in) :: x
		real(kind=8), intent(out) :: f
		end subroutine evalf	
	end interface
		
	!-------------------------------------------------------------------
	!     Initialization
	!-------------------------------------------------------------------  
			
	! Print problem information
	
	if ( iprint ) then
		write(*,1000)
		write(*,1010) ftol
		write(*,1020) m, stpmin
	end if
	
	! Check for errors in input data
	
	if ( stp < stpmin ) then
		write(*,1030)
		info = - 1
		return
	end if
	
	! Counters

	nfev = 0

  ! Compute maxg0
  
	maxg0 = maxval(g0)
		
	! Check if maxg0 is negative and define ftest 

	if ( maxg0 >= zero ) then
		write(*,1040)
		info = - 1
		return
	end if
	
	ftest = ftol * maxg0	
	
	! Test the sufficient descent condition (sdc) at stp 

	sdc = .true.     
	do i = 1, m
			
		call evalf(stp,f,i)
		nfev = nfev + 1	
		
		fend(i) = f	
		
		if ( f > f0(i) + ftest * stp ) then
			sdc = .false.
			iA  = i
			fiA = f
			exit
		end if				
	end do
	
	if ( sdc ) then
		info = 0
		return
	else
		call backtrackingMO(evalf,stp,stpmin,m,f0,g0,iA,fiA,fend,ftest,nfevbt,info)
		nfev = nfev + nfevbt
	end if
	
	! Non-executable statements
	
	1000 format(/,1X,'---------------------------------------------------------',/,1X,& 
			         ' An Armijo line search algorithm for vector Optimization ',/,1X,& 
				     '---------------------------------------------------------')

    1010 format(/,1X,'Sufficient descent condition tolerance:',1X,1P,D8.2)
 
	1020 format(/,1X,'Number of functions:',1X,I3,/,1X,'stpmin =',1X,1P,D8.2)	
			
	1030 format(/,'Error: stp < stpmin.')
	
	1040 format(/,'Error: the search direction is not a descent direction.')
	
end subroutine armijo

!***********************************************************************
!***********************************************************************	

subroutine backtrackingMO(evalf,stp,stpmin,m,f0,g0,iA,fiA,fend,ftest,nfev,info)	


	implicit none
	
	! SCALAR ARGUMENTS
	integer, intent(in) :: m,iA
	integer, intent(out) :: nfev,info
	real(kind=8), intent(in) :: stpmin,ftest,fiA
	real(kind=8), intent(inout) :: stp
	!logical,intent(in) :: iprint
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: f0(m), g0(m)
	real(kind=8), intent(out) :: fend(m)
	
	! LOCAL SCALARS
	integer :: i,ind,outiter,infoBT
	real(kind=8), parameter :: two = 2.0d0, sigma1 = 0.5d-1, sigma2 = 9.5d-1
	real(kind=8) :: f,stpq,stpt
	logical :: sdc
	
		! INTERFACE OF EXTERNAL SUBROUTINES 
	interface
		! SUBROUTINE evalf:
		subroutine evalf(x,f,ind)
		! SCALAR ARGUMENTS
		integer, intent(in) :: ind
		real(kind=8), intent(in) :: x
		real(kind=8), intent(out) :: f
		end subroutine evalf
	end interface

	! Counters
	
	outiter = 0
	nfev    = 0
		
	!-------------------------------------------------------------------
	!     Main loop
	!-------------------------------------------------------------------    

	Main_loop: do
	
		! Test the vector Armijo condition
		
		if ( outiter == 0 ) then
			sdc = .false.
			ind = iA
			f   = fiA
		elseif ( infoBT == 0 ) then
			sdc = .true.     
			do i = 1, m				
			
				if ( i == ind ) cycle
			
				call evalf(stp,f,i)
				nfev = nfev + 1		

				fend(i) = f
				
				if ( f > f0(i) + ftest * stp ) then
					sdc = .false.
					ind = i
					exit
				end if	
			end do
			
		end if
		
		! Finish backtracking with the current point
		
		if ( sdc ) then
			info = 0
			return
		end if
		
		! Test if stp is too small
		
		if ( stp <= stpmin ) then
			stp = stpmin
			info = 1
			
			do i = 1, m				
				
				if ( i == ind ) cycle
				
				call evalf(stp,f,i)
				nfev = nfev + 1	
				
				fend(i) = f
			end do
			return
		end if

		outiter = outiter + 1		
		
		! Compute new trial stepsize based on f_ind
		
		Inner_loop: do
		
			! Test Armijo condition for f_ind
			
			if ( f <= f0(ind) + ftest * stp ) then
				infoBT = 0
				exit Inner_loop
			end if
			
			if ( stp <= stpmin ) exit Inner_loop
			
			stpq = ( (g0(ind) / ( (f0(ind)-f) / stp + g0(ind) ) ) / two ) * stp
			
			stpt = stpq
			
			if ( stpt >= sigma1 * stp .and. stpt <= sigma2 * stp ) then
				stp = stpt
			else
				stp = stp / two
			end	if
			
			call evalf(stp,f,ind)
			nfev = nfev + 1	
                    
		end do Inner_loop
		
		fend(ind) = f
	
	end do Main_loop

	!--------------------------------------------------------------------- 
	!     End of main loop
	!---------------------------------------------------------------------
	
end subroutine backtrackingMO
