subroutine lsvecopt(evalf,evalg,quadstep,stp,stmin,stmax,m,mquad,      &
					mnquad,quadratic,quadind,nquadind,f0,g0,fend,ftol,gtol, &
					iprint,nfev,ngev,outiter,totiniter,info,LStype,tol)
					
	use globals, only: JFin

	implicit none
	
	! SCALAR ARGUMENTS
	integer, intent(in) :: m, mquad, mnquad, LStype
	integer, intent(out) :: nfev, ngev, outiter, totiniter, info
	real(kind=8), intent(in) :: stmin, stmax, ftol, gtol
	real(kind=8), intent(in), optional :: tol
	real(kind=8), intent(inout) :: stp
	logical,intent(in) :: iprint
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: f0(m), g0(m)
	real(kind=8), intent(out) :: fend(m)
	logical, intent(in) :: quadratic(m)
	integer, intent(in) :: quadind(mquad), nquadind(mnquad)

!	Reference:
!
!   L. R. Lucambio Pérez, and L. F. Prudente, "A Wolfe line search 
!	algorithm for vector optimization", ACM Transactions on Mathematical 
!   Software 45(4), pp. 37:1-37:23, 2019.
!	
!   Given continuous differentiable functions f_i:R->R for i = 1,...,m,
!	this subroutine is designed to find a positive step-size that 
!	satisfies the sufficient decrease condition
!
!		f_i(stp) <= f_i(0) + ftol * stp * f'_max(0) for all i = 1,...,m, 
!
!	and the curvature condition according to the line search type required:
!
!   * standard Wolfe condition (LStype = 1)
!
!		max{f'_i(stp) : i in {1,...m}} >= gtol * f'_max(0),
!
!   * strong Wolfe condition (LStype = 2)
!
!		abs(max{f'_i(stp) : i in {1,...m}}) <= - gtol * f'_max(0),
!
!   * restricted Wolfe condition (LStype = 3)
!
!		gtol * f'_max(0) <= max{f'_i(stp) : i in {1,...m}} <= tol,
!
!   where f'_max(0) = max{f'_i(0) : i in {1,...m}}.
!
!	Let F:Rn->Rp be a continuous vector-valued differentiable function, 
!	and K be a closed, convex, and pointed cone with non-empty interior 
!	in Rp. Assume that the generator set C of K is finite, i.e., there 
!	exist w_1,...,w_m belonging to K^*\{0} such that C = {w_1,...,w_m} 
!	and K^* = cone(conv(C)). If d in Rn is a K-descent direction of F at
!    x, then taking
!
!   f_i(stp) := < w_i, F( x + stp * d) >, for all i = 1,...,m,
!
!	the above sufficient decrease condition and curvature condition
!	correspond to the vector-valued Wolfe conditions for the step-size stp.
!
!	Mandatory external subroutines: evalf, evalg, and quadstep.
!
!	On entry the user most provide the following variables: stp,stmin,
!	stmax,m,mquad,mnquad,quadratic,quadind,nquadind,f0,g0,ftol,gtol.
!
!	evalf is an external subroutine (see the interface block).
!		evalf(x,f,ind) should return the evalution of the function f_ind
!	    at x.
!
!	evalg is an external subroutine (see the interface block).
!		evalg(x,g,ind) should return the evalution of the derivative of 
!		f'_ind at x.
!
!	quadstep is an external subroutine (see the interface block).
!		If f_ind is a convex quadratic function, then quadstep(x,ind,flag)
!	    should return the minimizer x of f_ind. If quadstep is coded for 
!		f_ind, then the output argument flag must be set to zero, 
!		otherwise flag must be definided as a non-null value.	
!
!	stp is a double precision variable.
!		On entry stp is the initial trial step-size. A positive initial 
!		estimate must be provided.
!		On exit stp is the final step-size. If info = 0 then stp satisfies 
!		the sufficient decrease and curvature conditions. See 'info' 
!		below to see the others stopping criteria.
!
!	stmin is a double precision variable.
!		On entry stmin is a nonnegative lower bound for the step-size.
!		On exit stmin is unchanged.
!
!	stpax is a double precision variable.
!		On entry stmax is a positive upper bound for the step.
!		On exit stmax is unchanged.
!
!	m is an integer variable.
!		On entry m is the number of functions.
!		On exit m is unchanged.
!
!	mquad is an integer variable.
!		On entry mquad is the number of (convex) quadratic functions.
!		On exit mquad is unchanged.
!
!	mnquad is an integer variable.
!		On entry mnquad is the number of non-quadratic functions.
!		On exit mnquad is unchanged.
!
!	quadratic(m) is a logical array.
!		On entry: quadratic(ind) = .true., if f_ind is a quadratic function
!	    	      quadratic(ind) = .false., otherwise.
!		If you are not sure if f_ind is a quadratic function, then set
!		quadratic(ind) = .false.
!		On exit quadratic is unchanged.
!
!	quadind(mquad) is a integer array.
!		On entry quadind must contain the indices ind's for which f_ind 
!		is a (convex) quadratic function.
!		On exit quadind is unchanged.
!
!	nquadind(mnquad) is a integer array.
!		On entry nquadind must contain the indices ind's for which f_ind 
!		is a non-quadratic function.
!		On exit nquadind is unchanged.
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
!	ftol is a double precision variable.
!		On entry ftol specifies a nonnegative parameter for the
!		sufficient decrease condition. ftol must be set less then or
!		equal to 1/2.
!		On exit ftol is unchanged.
!
!	gtol is a double precision variable.
!		On entry gtol specifies a nonnegative parameter for the
!		curvature condition. gtol must be less than ftol.
!		On exit gtol is unchanged.
!
!	iprint is a logical variable.
!		If iprint = .true. then the output is printed on the screen,
!		otherwise no output is shown.
!
!	nfev is a integer variable.
!		On exit nfev is the number of evaluations of the functions.
!
!	ngev is a integer variable.
!		On exit ngev is the number of evaluations of the derivatives of 
!		the functions.
!
!	outiter is a integer variable.
!		On exit outiter is the number of outer iterations.
!
!	totiniter is a integer variable.
!		On exit totiniter is the total number of inner iterations.
!
!	info is a integer variable.
!		On exit info contains the flag of the solution.
!	   	 0: stp satisfies the strong Wolfe conditions.
!		 1: stp = stpmin   
!		 2: stp = stpmax   
!		 3: rouding erros prevent progress
!		 4: the working interval is too small
!		 5: the number of iterations is exhausted
!		-1: an error ocurred during the execution of lsvecopt
!
	
	! LOCAL SCALARS
	real(kind=8), parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0,   &
	                           smallnum = -1.0d+99, bignum = 1.0d+99
	integer :: ind, indwork, i, j, initer, maxoutiter, & 
	           flag, MTinfo
	real(kind=8) :: f,g,xtol,tolLS,ftolinner,gtolinner,maxg0,maxg,gtest,&
	                ftest,stpmin,stpmax,stp0,stpquad
	logical :: sdc, cc, brackt

	! LOCAL ARRAYS
	integer :: isave(3)
	real(kind=8) :: dsave(13)
	character(len=60) :: task
	
	! INTERFACE OF MORE AND THUENTE SUBROUTINES 
	interface
		! SUBROUTINE evalf:
		subroutine evalf(x,f,ind)
		! SCALAR ARGUMENTS
		integer, intent(in) :: ind
		real(kind=8), intent(in) :: x
		real(kind=8), intent(out) :: f
		end subroutine evalf

		! SUBROUTINE evalg:
		subroutine evalg(x,g,ind)
		! SCALAR ARGUMENTS
		integer, intent(in) :: ind
		real(kind=8), intent(in) :: x
		real(kind=8), intent(out) :: g
		end subroutine evalg

		! SUBROUTINE quadstep:
		subroutine quadstep(x,ind,flag)
		! SCALAR ARGUMENTS
		integer, intent(in) :: ind
		integer, intent(out) :: flag
		real(kind=8), intent(out) :: x
		end subroutine quadstep 
		
		subroutine dcsrch(stp,f,g,ftol,gtol,tol,xtol,task,stpmin,stpmax,isave,dsave)
		! SCALAR ARGUMENTS
		real(kind=8), intent(in) :: ftol, gtol, tol, xtol, stpmin, stpmax
		real(kind=8), intent(inout) :: f, g, stp
		character(len=60), intent(inout) :: task    
		! ARRAY ARGUMENTS    
		integer, intent(inout) :: isave(2)
		real(kind=8), intent(inout) :: dsave(13)
		end subroutine dcsrch
							
		! SUBROUTINE dcstep:
		subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)
		! SCALAR ARGUMENTS
		real(kind=8), intent(in) :: fp,dp, stpmin,stpmax
		real(kind=8), intent(inout) :: stx,fx,dx,sty,fy,dy,stp
		logical, intent(inout) :: brackt
		end subroutine dcstep
	
	end interface
		
	!-------------------------------------------------------------------
	!     Initialization
	!-------------------------------------------------------------------  
	
	xtol = 1.d-20

	ftolinner = min( 1.1d0 * ftol, 0.75d0 * ftol + 0.25d0 * gtol )
	gtolinner = max( 0.9d0 * gtol, 0.25d0 * ftol + 0.75d0 * gtol )
	
	maxoutiter = 100
	
	stpmin = stmin
	stpmax = stmax

	brackt = .false.
		
	! Print problem information
	
	if ( iprint ) then
		write(*,1000)
		write(*,1010) ftol, gtol
		write(*,1020) m, stpmin, stpmax
	end if
	
	! Check for errors in input data
	
	if ( ftol >= gtol ) then
		write(*,1021)
		info = - 1
		return
	end if
	
	if ( mquad >= 1 .and. ftol > 0.5d0 ) then
		write(*,1022)
		info = - 1
		return
	end if
	
	if ( stp < stpmin ) then
		write(*,1023)
		info = - 1
		return
	end if
	
  if ( stp > stpmax ) then 
		write(*,1024)
		info = - 1
		return
	end if
	
	if ( LStype /= 1 .and. LStype /= 2 .and. LStype /= 3 ) then
		write(*,1025)
		info = - 1
		return
	end if
	
	if ( LStype == 3 .and. .not. present(tol) ) then
		write(*,1026)
		info = - 1
		return
	end if
	
	if ( LStype == 3 .and. tol < zero ) then
		write(*,1027)
		info = - 1
		return
	end if
	
	! Counters

	outiter   = 0
	totiniter = 0
	nfev      = 0
	ngev      = 0

  ! Compute maxg0
  
	maxg0 = maxval(g0)
	
	! Check if maxg0 is negative, and define ftest and gtest

	if ( maxg0 >= zero ) then
		write(*,1028)
		info = - 1
		return
	end if
	
	ftest =  ftol * maxg0	
	gtest = -gtol * maxg0	
	
	! Define the tolerance for the line search
	
	if ( LStype == 1 )  then
		tolLS = bignum
	elseif ( LStype == 2 )  then
		tolLS = gtest
	elseif ( LStype == 3 ) then		
		tolLS = tol
	end if
	
	! Initialize fend
	
	do i = 1,m
		fend(i) = 0.0d0
	end do
	
	! If there exists a convex quadratic function, then modify the initial
	! trial stepsize by setting
	! stp <-- min{ stp, stpquad }, 
	! where stpquad is the smallest minimizer of the quadratics functions.
	! Moreover, define the index for the bracketing phase by 
	! ind = argmin{ g0(i): i \in nquadind }. 
	
	if ( mquad >= 1 ) then
		do i = 1, mquad
			call quadstep(stp0,quadind(i),flag)
			if ( flag /= 0 ) then
				write(*,1030)
				stop
			end if
			
			if ( i == 1 ) then
				stpquad = stp0
			else
				stpquad = min( stpquad, stp0 )
			end if
		end do
		
		! Check whether stpquad < stpmin. In this case, stop with 
		! stp = stpmin
		
		if ( stpquad < stpmin ) then
			stp = stpmin

			if ( iprint ) then
				write(*,1029) stp
				write(*,1090) nfev, ngev
			end if
			info = 1
			return
		end if 
				
		if ( mquad == m ) then
			stp = min( stpmax, stpquad )
		else
			stp    = min( stp, stpquad )
			stpmax = min( stpmax, stpquad )
			ind    = nquadind( minloc( g0(nquadind),dim=1 ) )
		end if
	else	
		ind = minloc(g0,dim=1)
	end if
	
	if ( iprint ) then
		write(*,1040) gtest
		write(*,1050)
	end if
	
	! Test the sufficient descent condition (sdc) and the curvature 
	! condition (cc) at stp or try to identify a function for which the 
	! interval (0,stp) brackets a desired step-size
	
	! Sufficient descent condition (sdc)

	sdc = .true.     
	do j = 1, mnquad
		i = nquadind(j)
		
		call evalf(stp,f,i)
		nfev = nfev + 1
		
		fend(i) = f
		
		if ( f > f0(i) + ftest * stp ) then
			sdc    = .false.
			brackt = .true.
			ind    = i
			stpmax = stp 
			exit
		end if				
	end do
	
	! Curvature condition (cc) 

	cc = .false.
	if ( sdc ) then
		if ( mquad >= 1 .and. stp == stpquad ) then
			maxg = zero
		else
			maxg = smallnum
		end if
		
		do j = 1, mnquad
			i = nquadind(j)
			
			call evalg(stp,g,i)
			ngev = ngev + 1
			
			maxg = max( maxg , g )	
			
			if ( g > tolLS ) then
				brackt = .true.
				ind    = i
				stpmax = stp
				exit
			end if	
		end do	
		
		if ( maxg >= -gtest .and. maxg <= tolLS ) then
			cc = .true.
		end if
	end if 	
		
	!-------------------------------------------------------------------
	!     Main loop
	!-------------------------------------------------------------------    

	Main_loop: do
		
		!---------------------------------------------------------------
		!     Print information
		!---------------------------------------------------------------
		
		if ( iprint ) then
			if ( outiter == 0 )	then
				write(*,1060) outiter, stp, stpmin, stpmax, sdc, cc, brackt
			else
				write(*,1070) outiter, indwork, stp, stpmin, stpmax, sdc,  &
							  cc, brackt, initer
			end if
		end if

		!---------------------------------------------------------------
		!     Test stopping criteria
		!---------------------------------------------------------------
		
		! Test the fulfillment of the sufficient descent condition (sdc)
		! and the curvature condition (cc)
		
		if ( sdc .and. cc ) then
			if ( iprint ) then
				write(*,1080) stp
				write(*,1090) nfev, ngev
			end if
			info = 0
			
			! Evaluate the quadratic functions at stp
			if ( mquad >= 1 ) then
				do j = 1, mquad
				
					i = quadind(j)
					
					call evalf(stp,f,i)
					nfev = nfev + 1
					
					fend(i) = f
					
					call evalg(stp,g,i)
					ngev = ngev + 1
				
				end do
			end if
			
			return
		end if

		! Test for stopping at stp = stpmin   

		if ( outiter > 0 .and. MTinfo == 1 ) then
			if ( iprint ) then
				write(*,1100) stp
				write(*,1090) nfev, ngev
			end if
			info = 1
			return
		end if
		
		! Test for stopping at stp = stpmax  

		if ( .not. brackt .and. stp == stmax ) then
			if ( iprint ) then
				write(*,1110) stp
				write(*,1090) nfev, ngev
			end if
			info = 2
			
			if ( mquad >= 1 ) then
				do j = 1, mquad
				
					i = quadind(j)
					
					call evalf(stp,f,i)
					nfev = nfev + 1
					
					fend(i) = f
					
					call evalg(stp,g,i)
					ngev = ngev + 1
				end do
			end if
			
			return
		end if
		
		! Test whether rouding erros prevent progress

		if ( outiter > 0 .and. MTinfo == 3 ) then
		    if ( iprint ) then
				write(*,1120) stp
				write(*,1090) nfev, ngev
			end if
			info = 3
			return
		end if
		
		! Test whether the working interval is too small

		if ( outiter > 0 .and. MTinfo == 4 ) then
		    if ( iprint ) then
				write(*,1130) stp
				write(*,1090) nfev, ngev
			end if
			info = 4
			return
		end if
		
		! Test whether the number of iterations is exhausted
		
		if ( outiter >= maxoutiter ) then
			if ( iprint ) then
				write(*,1140) stp
				write(*,1090) nfev, ngev
			end if
			info = 5
			return
		end if
		
		!-------------------------------------------------------------------  
		!     Iteration
		!------------------------------------------------------------------- 

		outiter = outiter + 1

		f = f0(ind)
		g = maxg0
		
		task   = 'START'
		initer = 0
	
		!---------------------------------------------------------------
		!     Call the (modified) Moré and Thuente line search
		!---------------------------------------------------------------
		 
		! J. J. Moré and D. J. Thuente, Line Search Algorithms with Guaranteed
		! Sufficient Decrease, ACM Trans. Math. Softw., 20 (1994), pp. 286–307.
		! http://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/
		
		Inner_loop: do 
		
			if ( task /= 'START' ) initer = initer + 1

			call dcsrch(stp,f,g,ftolinner,gtolinner,tolLS,xtol,task,stpmin,stpmax,isave,dsave)
			
			! Identify the flag of the Moré and Thuente algorithm
	
			if ( task == 'FG' ) then
					MTinfo = -1
			elseif ( task == 'CONVERGENCE' ) then
					MTinfo =  0
			elseif ( task == 'WARNING: STP = STPMIN' ) then
					MTinfo =  1
			elseif ( task == 'WARNING: STP = STPMAX' ) then
					MTinfo =  2
			elseif ( task == 'WARNING: ROUNDING ERRORS PREVENT PROGRESS' ) then
					MTinfo =  3
			elseif ( task == 'WARNING: XTOL TEST SATISFIED' ) then
					MTinfo =  4
			end if
			
			! 'isave(1) = 1' means that 'brackt = true'
			
			if ( isave(1) == 1 ) brackt = .true.
			
			! Test for convergence of the inner method
			
			if ( MTinfo == 0 .or. MTinfo == 1 .or. MTinfo == 2 .or.    &
					 MTinfo == 3 .or. MTinfo == 4 ) exit Inner_loop	
			
			! If 'brackt = false', then the actual step-size stp is greater 
			! than the previous step-size stpprev. It follows that
			! stp \in [ min{ delta * stpprev, stpmax}, stpmax ], where 
			! delta = 1.1. In this case, exit the inner loop.
			
			if ( initer > 0 .and. .not. brackt ) exit Inner_loop			
			
			! Evaluate f and g
			
			call evalf(stp,f,ind)
			call evalg(stp,g,ind)
			nfev = nfev + 1
			ngev = ngev + 1		
			
		end do Inner_loop
		
		totiniter = totiniter + initer	
		
		!---------------------------------------------------------------
		!     Prepare for the next iteration
		!---------------------------------------------------------------
		
		! Save the working index and the store the objective function at stp 
		
		indwork = ind
		
		fend(ind) = f

		if ( MTinfo /= 1 .and. MTinfo /= 3 .and. MTinfo /= 4 ) then
		
			! Test the sufficient descent condition (sdc) and the curvature 
			! condition (cc) at stp or try to identify a function for which 
			! the interval (0,stp) brackets a desired step-size
			
			! Sufficient descent condition (sdc)
			
			sdc = .true.     
			do j = 1, mnquad
			
!				i = nquadind(j)
!				if ( MTinfo == 0 .and. i == indwork ) cycle
				
!				if ( MTinfo == -1 .and. j == 1 ) then
!				     i = indwork
!				end if
				
				
				if ( MTinfo == -1 .and. j == 1 ) then
				     i = indwork
				else
				     i = nquadind(j)
				end if
				
				if ( MTinfo == 0 .and. i == indwork ) cycle
				
				if ( MTinfo == -1 .and. j > 1 .and. i == indwork ) i = nquadind(1)
				
				
				call evalf(stp,f,i)
				nfev = nfev + 1		
				
				fend(i) = f
				
				if ( f > f0(i) + ftest * stp ) then
					sdc    = .false.
					brackt = .true.
					ind    = i
					stpmax = stp 
					exit
				end if				
			end do
			
			! Curvature condition (cc) 
			
			cc = .false.
			if ( sdc ) then
			
				if ( mquad >= 1 .and. stp == stpquad ) then
					if ( MTinfo == 0 ) then
						maxg = max( zero, g )
					else
						maxg = zero
					end if
				elseif (  MTinfo == 0 ) then
					maxg = g
				else
					maxg = smallnum
				end if
				
				do j = 1, mnquad
			
!					i = nquadind(j)
!					if ( MTinfo == 0 .and. i == indwork ) cycle
					
					
					if ( MTinfo == -1 .and. j == 1 ) then
					    i = indwork
					else
						i = nquadind(j)
					end if
					
					if ( MTinfo == 0 .and. i == indwork ) cycle
					
					if ( MTinfo == -1 .and. j > 1 .and. i == indwork ) i = nquadind(1)
					
					call evalg(stp,g,i)
					ngev = ngev + 1
					
					maxg = max( maxg , g )	
					
					if ( g > tolLS ) then
						brackt = .true.
						ind    = i
						stpmax = stp
						exit
					end if	
								
				end do	

				if ( maxg >= -gtest .and. maxg <= tolLS ) then
					cc = .true.
				end if
			end if 
		
  	end if  	
		
	!--------------------------------------------------------------------- 
	!     Iterate
	!---------------------------------------------------------------------    

	end do Main_loop

	!--------------------------------------------------------------------- 
	!     End of main loop
	!---------------------------------------------------------------------
	
	! Non-executable statements
	
	1000 format(/,1X,'-----------------------------------------------------',/,1X,& 
					 'lsvecopt: A Wolfe line search for vector Optimization',/,1X,& 
				     '-----------------------------------------------------')

  1010 format(/,1X,'Sufficient descent condition tolerance:',1X,1P,D8.2,&
              /,1X,'Curvature condition tolerance:         '1X,1P,D8.2)
 
	1020 format(/,1X,'Number of functions:',1X,I3,/,1X,'stpmin =',1X,1P,D8.2,/,1X,'stpmax =',1X,1P,D8.2)	
	
	1021 format(/,'Error of lsvecopt: ftol >= gtol.')
	
	1022 format(/,'Error of lsvecopt: ftol > 1/2.')
	
	1023 format(/,'Error of lsvecopt: stp < stpmin.')
	
	1024 format(/,'Error of lsvecopt: stp > stpmax.')
	
	1025 format(/,'Error of lsvecopt: invalid value for LStype.')
	
	1026 format(/,'Error of lsvecopt: tol is mandatory for LStype = 2.')
	
	1027 format(/,'Error of lsvecopt: tol must be a non-negative number.')
	
	1028 format(/,'Error of lsvecopt: the search direction is not a descent direction.')
	
	1029 format(/,1X,'Flag of lsvecopt: the smallest minimizer of the quadratics functions',/,19X,&
	                 'is less that stpmin. The final step is equal to stpmin.',/,/,1X,'Final step: ',1P,D8.2)
	
	1030 format(/,'Error of lsvecopt: the minimizer of a quadratic function was not provided.')
	
	1040 format(/,1X,'Curvature condition test: -gtol * maxg0 =',1X,1P,D8.2)
	
	1050 format(/,2X,'out',1X,'ind',4X,'stp',5X,'stpmin',4X,'stpmax',3X,'SDC',3X,'CC',2X,'brackt',1X,'inner')
	
	1060 format(1X,I3,3X,'-',2X,1P,D8.2,2X,1P,D8.2,2X,1P,D8.2,3X,L1,4X,L1,5X,L1,6X,'-')
	
	1070 format(1X,I3,1X,I3,2X,1P,D8.2,2X,1P,D8.2,2X,1P,D8.2,3X,L1,4X,L1,5X,L1,4X,I3)		
	
	1080 format(/,1X,'Flag of lsvecopt: a point satisfying the sufficient descent condition',/,19X,&
	                 'and the curvature condition was found.',/,/,1X,'Final step: ',1P,D8.2)						
	
	1090 format(/,1X,'Number of functions evaluations:   ',I4,/,1X,'Number of derivatives evaluations: ',I4)
	
	1100 format(/,1X,'Flag of lsvecopt: the final step is equal to stpmin.',/,/,1X,'Final step: ',1P,D8.2)
	
	1110 format(/,1X,'Flag of lsvecopt: the final step is equal to stpmax.',/,/,1X,'Final step: ',1P,D8.2)
	
	1120 format(/,1X,'Flag of lsvecopt: rouding erros prevent progress.',/,/,1X,'Final step: ',1P,D8.2)
	
	1130 format(/,1X,'Flag of lsvecopt: the working interval is too small.',/,/,1X,'Final step: ',1P,D8.2)
	
	1140 format(/,1X,'Flag of lsvecopt: maximum of iterations reached.',/,/,1X,'Final step: ',1P,D8.2)

end subroutine lsvecopt
