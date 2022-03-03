subroutine StBFGSWolfe(n,m,x,l,u,scaleF,checkder,outiter,time,nfev,ngev,theta,inform)

	use globals , only: sF,nfinner,xinner,d,JF,JFin,B,quadratic
	use myproblem

	implicit none
	
	! SCALARS ARGUMENTS
	integer,intent(in)       :: n,m
	logical,intent(in)       :: scaleF,checkder
	integer,intent(out)      :: outiter,nfev,ngev,inform
	real(kind=4),intent(out) :: time
	real(kind=8),intent(out) :: theta
	
	! ARRAY ARGUMENTS
	real(kind=8),intent(inout) :: x(n)
	real(kind=8),intent(in)    :: l(n),u(n)
	
	! LOCAL SCALARS
	integer :: ind,j,mquad,mnquad,infoIS,infoLS,nfevLS,ngevLS,infofun,allocerr,&
	           maxoutiter,informfg,LStype,outiterLS,totiniterLS
	real(kind=8) :: stp,stpmin,stpmax,epsopt,ftol,gtol,tolLS
	real(kind=8), parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0,   &
	                           bignum = 1.0d+99
	real(kind=4) :: timeini,timefin,timeiniIS,timefinIS
	logical :: iprint,iprintLS

	! LOCAL ARRAYS
	integer,allocatable :: quadind(:),nquadind(:)
	real(kind=8) :: g(n),xprev(n),Fx(m),JFprev(m,n),phi0(m),gphi0(m),fend(m)
	
	! EXTERNAL SUBROUTINES
	external :: evalphi,evalgphi,evalqstep
	
  !
  ! inform:
  !
  !  0: Optimality satisfied
  !  1: Maximum number of iterations reached
  ! -1: An error occurred during the execution of the algorithm
  !
  ! March of 2022.
  
	interface
		subroutine quadvar(quadratic,m,mquad,mnquad,quadind,nquadind)
		! SCALAR ARGUMENTS
		integer, intent(in)  :: m
		integer, intent(out) :: mquad,mnquad
		! ARRAY ARGUMENTS
		integer,allocatable,intent(out) :: quadind(:),nquadind(:)
		logical,intent(in) :: quadratic(m)
		end subroutine quadvar
	end interface
  
  	!---------------------------------------------------------------------   
	!     Initialization
	!---------------------------------------------------------------------
		
	! Check derivatives
	
	if ( checkder ) call checkdF(n,m,x,l,u)
	
	! Start timing
		
	call cpu_time(timeini)
	
	! Set default parameters 

	epsopt = 5.0d0 * sqrt( 2.0d0 ** (-52) )
	ftol   = 1.d-04
	gtol   = 1.d-01
		
	maxoutiter = 2000
	stpmin     = 1.0d-15
	stpmax     = 1.0d+10
		
	! Allocating the global arrays
	
	allocate(d(n),xinner(n),sF(m),quadratic(m),JF(m,n),JFin(m,n),B(n,n,m),stat=allocerr)
	if ( allocerr .ne. 0 ) then
		write(*,*) 'Allocation error in main program'
		stop
	end if
	
	! Do not use quadratic informations of the objectives functions
	
	do ind = 1,m
		quadratic(ind) = .false.
	end do
	
	call quadvar(quadratic,m,mquad,mnquad,quadind,nquadind)
	
	! Saving the inner informations
	
	nfinner = n
	xinner  = x
	
	! Print problem information	
	
	iprint   = .true.	
	iprintLS = .false.
	
	if ( iprint ) then
		write(*,1000)
		write(*,1010) n, m
	end if
	
	! Scale problem
	
	call scalefactor(n,m,x,scaleF)	

	! Counters
	
	outiter = 0	
	nfev    = 0
	ngev    = 0
	
	!------------------------------------------------------------------- 
	!     Main loop
	!-------------------------------------------------------------------    

	Main_loop: do
	
		!---------------------------------------------------------------
		!     Prepare the iteration
		!---------------------------------------------------------------
		
		! Compute the Jacobian JF
		
		if ( outiter > 0 ) JFprev = JF
		
		if ( outiter > 0 .and. ( infoLS == 0 .or. infoLS == 2 ) ) then
			JF = JFin
		else
			do ind = 1,m
				call sevalg(n,x,g,ind,infofun)
				ngev = ngev + 1
				
				JF(ind,:) = g
				
				if ( infofun /= 0 ) then
					inform = -1
					
					! Stop timing
			
					call cpu_time(timefin)
					
					time = timefin - timeini
			
					deallocate(d,xinner,sF,quadratic,quadind,nquadind,JF,JFin,B)
					
					return
				end if
			end do
		
		end if
		
		! Compute the BFGS matrices
	
		if ( outiter == 0 ) then
			do ind = 1,m
				B(:,:,ind) = 0.0d0
				forall (j = 1:n) B(j,j,ind) = 1.0d0
			end do
		else
			call BFGSupdateCautions(n,m,x,xprev,JF,JFprev,theta,B)
		end if
			
		! Compute the search direction
		
		call evaldN(n,m,d,infoIS)
		
		! Compute theta	
		
		if ( infoIS /= 0 ) then
			theta = bignum
			
			! Stop timing
			
			call cpu_time(timefin)
			
			time = timefin - timeini
					
			inform = -1
			
			deallocate(d,xinner,sF,quadratic,quadind,nquadind,JF,JFin,B)
			return
		else
			theta = - bignum
			do ind = 1,m
				theta = max( theta, dot_product( JF(ind,:),d )             &
						 + 0.5d0 * dot_product( matmul( B(:,:,ind),d ),d ) )
			end do
		end if
		
		! Print information

		if ( iprint ) then
			if ( outiter == 0 ) then
				write(*,1030) epsopt
				write(*,1040)
				write(*,1050) outiter,abs(theta),infoIS,nfev,ngev
			else
				if ( mod(outiter,10) == 0 ) write(*,1040)
				write(*,1060) outiter,abs(theta),infoLS,infoIS,nfev,ngev
			end if
		end if
	
		!-------------------------------------------------------------------
		!     Test stopping criteria
		!-------------------------------------------------------------------
		
		! Test optimality	
		
		if ( abs( theta ) <= epsopt  ) then
	
			inform = 0
			
			! Stop timing
			
			call cpu_time(timefin)
			
			time = timefin - timeini
  
			if ( iprint ) write(*,1070) nfev, ngev, time
			deallocate(d,xinner,sF,quadratic,quadind,nquadind,JF,JFin,B)
			return
			
		end if		
		
		! Test whether the number of iterations is exhausted
		
		if ( outiter >= maxoutiter ) then
				inform = 1
				
				! Stop timing
			
				call cpu_time(timefin)
			
				time = timefin - timeini
					
				if ( iprint ) write(*,1080) nfev, ngev, time

				deallocate(d,xinner,sF,quadratic,quadind,nquadind,JF,JFin,B)
				return
		end if
		
		! Test for errors in the inner solver
		
		if (  infoIS /= 0 ) then
				inform = 2
				
				! Stop timing
			
				call cpu_time(timefin)
			
				time = timefin - timeini
					
				if ( iprint ) write(*,1090) nfev, ngev, time

				deallocate(d,xinner,sF,quadratic,quadind,nquadind,JF,JFin,B)
				return
		end if

		!-------------------------------------------------------------------
		!     Iteration
		!-------------------------------------------------------------------
		
		! Increment outiter
				
		outiter = outiter + 1
						
		! Compute initial step
			
		stp = one

		! Compute phi0 = [phi_1(0),...,phi_m(0)] = [f1(x),...,fm(x)]
		
		if ( outiter > 1 .and. ( infoLS == 0 .or. infoLS == 2 ) ) then
			phi0 = fend
		else
			do ind = 1,m
				call sevalf(n,x,phi0(ind),ind,informfg)
				nfev = nfev + 1   
			end do
		
		end if
		
		! Compute gphi0 = [gphi_1(0),...,gphi_m(0)]
		
		gphi0 = matmul( JF, d )			
													  
		! Compute the stepsize satisfying the Wolfe conditions
		
		LStype = 1   ! Standard Wolfe condition	
		tolLS  = bignum
		
		call lsvecopt(evalphi,evalgphi,evalqstep,stp,stpmin,stpmax,m,mquad,&
					mnquad,quadratic,quadind,nquadind,phi0,gphi0,fend,ftol, &
					gtol,iprintLS,nfevLS,ngevLS,outiterLS,totiniterLS,infoLS,LStype,tolLS)
					
		! Check for errors in lsvecopt
		
		if ( infoLS == - 1 ) then
			inform = - 1

			deallocate(d,xinner,sF,quadratic,quadind,nquadind,JF,JFin,B)
			return
		end if
				 
		nfev = nfev + nfevLS
		ngev = ngev + ngevLS
		
		! Update x	
		
		xprev = x
				
		x = x + stp * d
		
		! Saving the inner informations
		
		xinner  = x
		
	!---------------------------------------------------------------------
	!     Iterate
	!---------------------------------------------------------------------    

	end do Main_loop

	!--------------------------------------------------------------------- 
	!     End of main loop
	!---------------------------------------------------------------------
	
	! Non-executable statements
	
	1000 format(/,1X,'-------------------------------------------------------------------------------',/,1X,& 
					 '         Standard BFGS-Wolfe method for Multiobjective Optimization            ',/,1X,& 
					 '-------------------------------------------------------------------------------')
					 								 
	1010 format(/,1X,'Number of variables:',1X,I6,/,1X,'Number of functions:',1X,I6)	
	
	1030 format(/,1X,'Optimality tolerance:',1X,1P,D7.1)
	
	1040 format(/,4X,'out',2X,'|theta|',3X,'LS',1X,'IS',1X,'#evalf',2X,'#evalg')
	
	1050 format(1X,I5,3X,1P,D8.2,3X,'-',2X,I1,1X,I6,2X,I6)
	
	1060 format(1X,I5,3X,1P,D8.2,3X,I1,2X,I1,1X,I6,2X,I6)		
	
	1070 format(/,1X,'Flag of MOPsolver: solution was found',/,/,1X, &
					 'Number of functions evaluations:               ',I6,/,1X, &
					 'Number of derivatives evaluations:             ',I6/,/,1X, &
					 'Total CPU time in seconds: ',0P,F8.2)						
	
	1080 format(/,1X,'Flag of MOPsolver: maximum of iterations reached',/,/,1X,&
					 'Number of functions evaluations:               ',I6,/,1X, &
					 'Number of derivatives evaluations:             ',I6/,/,1X, &
					 'Total CPU time in seconds: ',0P,F8.2)
					 
	1090 format(/,1X,'Flag of MOPsolver: it was not possible to solve the subproblem',/,/,1X,&
					 'Number of functions evaluations:               ',I6,/,1X, &
					 'Number of derivatives evaluations:             ',I6/,/,1X, &
					 'Total CPU time in seconds: ',0P,F8.2)

end subroutine StBFGSWolfe
