module myproblem

use globals, only: problem

implicit none

! LOCAL SCALARS
integer :: dimm,dimk
logical :: penalize
real(kind=8) :: alpha
real(kind=8),parameter :: penparam = 1.0d10

! LOCAL ARRAYS
real(kind=8),pointer :: lin(:),uin(:)


  ! SUBROUTINES
  public :: inip, evalf, evalg, quadstep
  private:: penalizef,penalizeg

contains

	subroutine inip(n,m,x,l,u,strconvex,scaleF,checkder)
	
	use globals, only: seed

	implicit none

	! SCALAR ARGUMENTS
	integer, intent(out) :: n,m
	logical, intent(out) :: scaleF,checkder
	
	! ARRAY ARGUMENTS
	real(kind=8), allocatable, intent(out) :: x(:)
	real(kind=8), allocatable, target, intent(out) :: l(:),u(:)
	logical, allocatable, intent(out) :: strconvex(:)
	
	! LOCAL SCALARS
	integer :: i,allocerr
	real, parameter :: pi = 3.1415927
	real(kind=8), parameter :: zero = 0.0d0, one = 1.0d0
	real(kind=8) :: a,b
	
	! FUNCTIONS
	real(kind=8) :: drand
		
	! ----------------------------------------------------------------------

	! AP1
	! Example 1 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP1' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		strconvex(3) = .true.
		
		! Box constraints
	
		l(:) = - 1.0d1
		u(:) =   1.0d1
		
		! Initial point
		
		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------

	! AP2
	! Example 2 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP2' ) then 
	
		! Number of variables
		
		n = 1
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		
		! Box constraints
	
		l(:) = - 1.0d2
		u(:) =   1.0d2
		
		! Initial point
		
		a = - 1.0d2
		b =   1.0d2
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
	end if	
	
! ----------------------------------------------------------------------

	! AP3
	! Example 3 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP3' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints
		
		l(:) = - 1.0d2
		u(:) =   1.0d2
		
		! Initial point
		
		a = - 1.0d2
		b =   1.0d2
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	! AP4
	! Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"
	
	if ( problem == 'AP4' ) then 
	
		! Number of variables
		
		n = 3
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		strconvex(3) = .true.
		
		! Box constraints
	
		l(:) = - 1.0d1
		u(:) =   1.0d1
		
		! Initial point
		
		a = - 1.0d1
		b =   1.0d1
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------

	!  BK1
	!  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'BK1' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		
		! Box constraints
	
		l(:) = - 5.0d0
		u(:) =   1.0d1
		
		! Initial point
		
		a = - 5.0d0
		b =   1.0d1
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if

! ----------------------------------------------------------------------

	!   DD1
	!   I. Das and J. E. Dennis. Normal-boundary intersection: A new method for generating
	!		the Pareto surface in nonlinear multicriteria optimization problems. SIAM J. Optim.,
	!		8(3):631–657, 1998.
	
	if ( problem == 'DD1' ) then 

		! Number of variables

		n = 5
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .false.
		
		! Box constraints
	
		l(:) = - 2.0d1
		u(:) =   2.0d1
		
		! Initial point

		a = - 1.0d0
		b =   1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		x = 1.0d1 * drand(seed) * x / norm2(x)
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			
	
! ----------------------------------------------------------------------

	!   DGO1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'DGO1' ) then 
	
		! Number of variables

		n = 1
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints
	
		l(:) = - 1.0d1
		u(:) =   1.3d1
		
		! Initial point
		
		a = - 1.0d1
		b =   1.3d1
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	!   DGO2
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'DGO2' ) then 

		! Number of variables

		n = 1
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		
		! Box constraints
	
		l(:) = - 9.0d0
		u(:) =   9.0d0
		
		! Initial point
		
		a = - 9.0d0
		b =   9.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------


	!   DTLZ1
	!   Scalable test problems for evolutionary multi-objective optimization
	
	if ( problem == 'DTLZ1' ) then 
	
		! Number of objectives
		
		m = 3
		
		dimm = m

		! Number of variables
		
		dimk = 5

		n = dimk + m - 1
		
!		dimk = n - m + 1
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints
	
		l(:) = 0.0d0
		u(:) = 1.0d0
		
		! Initial point
		
		a = 0.0d0
		b = 1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
	
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
			
		return
	end if	
	
! ----------------------------------------------------------------------


	!   DTLZ2
	!   Scalable test problems for evolutionary multi-objective optimization
	
	if ( problem == 'DTLZ2' ) then 
	
		! Number of objectives
		
		m = 3
		
		dimm = m

		! Number of variables
		
		dimk = 5

		n = dimk + m - 1
		
!		dimk = n - m + 1
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints
	
		l(:) = 0.0d0
		u(:) = 1.0d0
		
		! Initial point
		
		a = 0.0d0
		b = 1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
	
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
			
		return
	end if
	
! ----------------------------------------------------------------------


	!   DTLZ3
	!   Scalable test problems for evolutionary multi-objective optimization
	
	if ( problem == 'DTLZ3' ) then 
	
		! Number of objectives
		
		m = 3
		
		dimm = m

		! Number of variables
		
		dimk = 5

		n = dimk + m - 1
		
!		dimk = n - m + 1
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints
	
		l(:) = 0.0d0
		u(:) = 1.0d0
		
		! Initial point
		
		a = 0.0d0
		b = 1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
	
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
			
		return
	end if	

! ----------------------------------------------------------------------


	!   DTLZ4
	!   Scalable test problems for evolutionary multi-objective optimization
	
	if ( problem == 'DTLZ4' ) then 
	
		! Number of objectives
		
		m = 3
		
		dimm = m

		! Number of variables
		
		dimk = 5

		n = dimk + m - 1
		
!		dimk = n - m + 1
		
		! Define parameter alpha
		
		alpha = 2.0d0
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints
	
		l(:) = 0.0d0
		u(:) = 1.0d0
		
		! Initial point
		
		a = 0.0d0
		b = 1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
	
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
			
		return
	end if	

! ----------------------------------------------------------------------		
	

	!   FA1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'FA1' ) then 
	
		! Number of variables

		n = 3
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		strconvex(3) = .false.
		
		! Box constraints
	
		l(:) = 0.0d0
		u(:) = 1.0d0
		
		! Initial point
		
		a = 0.0d0
		b = 1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			
	
! ----------------------------------------------------------------------

	!   Far1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'Far1' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints
	
		l(:) = - 1.0d0
		u(:) =   1.0d0
		
		! Initial point

		a = - 1.0d0
		b =   1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------
	
	! FDS
	! NEWTON’S METHOD FOR MULTIOBJECTIVE OPTIMIZATION
	
	if ( problem == 'FDS' ) then 
	
		! Number of variables
		
		n = 5
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		strconvex(3) = .true.
		
		! Box constraints
	
		l(:) = - 2.0d0
		u(:) =   2.0d0
		
		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			
	
! ----------------------------------------------------------------------
	
	! FF1
	! C. M. Fonseca and P. J. Fleming, “An overview of evolutionary algorithms in multiobjective optimization
	
	if ( problem == 'FF1' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0
		
		! Initial point
		
		a = - 1.0d0
		b =   1.0d0
		
		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
				
	! ----------------------------------------------------------------------	

	!   Hil1
	!   Generalized Homotopy Approach to Multiobjective Optimization
	
	if ( problem == 'Hil1' ) then 
	
		! Number of variables

		n = 2

		! Number of objectives

		m = 2		
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			write(*,*) 'Allocation error in main program'
			stop
		end if
		
		! Strictly strconvex?
			
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = zero
		u(:) = one

		! Initial point

		a =  0.0d0
		b =  1.0d0

		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------	

	!   IKK1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( problem == 'IKK1' ) then 
	
		! Number of variables

		n = 2

		! Number of objectives

		m = 3	
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			write(*,*) 'Allocation error in main program'
			stop
		end if
		
		! Strictly strconvex?
			
		strconvex(1) = .false.
		strconvex(2) = .false.
		strconvex(3) = .false.
		
		! Box constraints

		l(:) = - 5.0d1
		u(:) =   5.0d1

		! Initial point

		a = - 5.0d1
		b =   5.0d1

		do i = 1,n
			x(i) = a + ( b - a ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------	

	!   IM1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( problem == 'IM1' ) then 
	
		! Number of variables

		n = 2

		! Number of objectives

		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(1) = 1.0d0
		u(1) = 4.0d0
		
		l(2) = 1.0d0
		u(2) = 2.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if

! ----------------------------------------------------------------------	

	! JOS1
	! Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why Does It Work and How?
	
	if ( problem == 'JOS1' ) then 

		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = - 1.0d2
		u(:) =   1.0d2
		
		! Initial point
		
		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------	

	! JOS4
	! Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why Does It Work and How?
	! See also: NEWTON’S METHOD FOR MULTIOBJECTIVE OPTIMIZATION
	
	if ( problem == 'JOS4' ) then 

		! Number of variables

		n = 20
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 1.0d-2
		u(:) = 1.0d0
		
		! Initial point
		
		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	

! ----------------------------------------------------------------------

	! 	KW2
	!   I.Y. Kim, O.L. de Weck, Adaptive weighted-sum method for bi-objective optimization: Pareto front generation
	
	if ( problem == 'KW2' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = -3.0d0
		u(:) =  3.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------	

	!   LE1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( problem == 'LE1' ) then 
	
		! Number of variables

		n = 2

		! Number of objectives

		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = -5.0d0
		u(:) =  1.0d1

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------
	
	! Lov1
	! "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( problem == 'Lov1' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = -1.0d1
		u(:) =  1.0d1

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------
	
	! Lov2
	! "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( problem == 'Lov2' ) then 

		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = -0.75d0
		u(:) =  0.75d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------
	
	! Lov3
	! "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( problem == 'Lov3' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = -2.0d1
		u(:) =  2.0d1

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------
	
	! Lov4
	! "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( problem == 'Lov4' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = -2.0d1
		u(:) =  2.0d1

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------
	
	! Lov5
	! "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( problem == 'Lov5' ) then 
	
		! Number of variables
		
		n = 3
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = -2.0d0
		u(:) =  2.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------
	
	! Lov6
	! "Singular Continuation: Generating Piecewise Linear Approximations to Pareto Sets via Global Analysis"
	
	if ( problem == 'Lov6' ) then 
	
		! Number of variables
		
		n = 6
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(1) = 0.1d0
		u(1) = 0.425
		
		l(2:6) = - 0.16d0
		u(2:6) =   0.16d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			

! ----------------------------------------------------------------------

	!  LTDZ
	!	 Combining convergence and diversity in evolutionary multiobjective optimization
	
	if ( problem == 'LTDZ' ) then 
	
		! Number of variables

		n = 3
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		strconvex(3) = .false.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------	

	! 	MGH9 (Gaussian)

	! More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	! Softw. 7(1), 17–41 (1981)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)

	if ( problem == 'MGH9' ) then 
	
		! Number of variables

		n = 3
		
		! Number of objectives
		
		m = 15
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints

		l(:) = - 2.0d0
		u(:) =   2.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
	
! ----------------------------------------------------------------------

	! 	MGH16 (Brown and Dennis)

	! More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	! Softw. 7(1), 17–41 (1981)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)
	
	if ( problem == 'MGH16' ) then 
	
		! Number of variables

		n = 4
		
		! Number of objectives
		
		m = 5
				
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		strconvex(3) = .false.
		strconvex(4) = .false.
		strconvex(5) = .false.
		
		! Box constraints

		l(1) = - 2.5d1
		u(1) =   2.5d1
		
		l(2:3) = - 5.0d0
		u(2:3) =   5.0d0

		l(4) = - 1.0d0
		u(4) =   1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	! 	MGH26 (Trigonometric)

	! More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	! Softw. 7(1), 17–41 (1981)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)

	
	if ( problem == 'MGH26' ) then 
	
		! Number of variables

		n = 4
														
		! Number of objectives
		
		m = n
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------

	! 	MGH33 (Linear function - rank 1)

	! More, J.J., Garbow, B.S., Hillstrom, K.E.: Testing unconstrained optimization software. ACM T. Math.
	! Softw. 7(1), 17–41 (1981)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)
	
	if ( problem == 'MGH33' ) then 
	
		! Number of variables

		n = 10
		
		! Number of objectives (m>=n)
		
		m = n
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex? (Each function is only convex)
				
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------	

	!   MHHM2
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit

	
	if ( problem == 'MHHM2' ) then 
	
		! Number of variables

		n = 2

		! Number of objectives

		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		strconvex(3) = .true.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
	
! ----------------------------------------------------------------------

	!   MLF1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'MLF1' ) then 
	
		! Number of variables

		n = 1
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 2.0d1

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	!   MLF2
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'MLF2' ) then 

		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = - 1.0d2
		u(:) =   1.0d2

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------
	
	!  MMR1
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	
	
	if ( problem == 'MMR1' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(1) = 0.1d0
		u(1) = 1.0d0
		
		l(2) = 0.0d0
		u(2) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------
	
	!  MMR2
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	
	
	if ( problem == 'MMR2' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------
	
	!  MMR3
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	
	
	if ( problem == 'MMR3' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = -1.0d0
		u(:) =  1.0d0
		
		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------
	
	!  MMR4
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	
	
	if ( problem == 'MMR4' ) then 

		! Number of variables

		n = 3
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 4.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!  MOP2
	!  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit		
	
	if ( problem == 'MOP2' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!  MOP3	
	!  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit	
	
	if ( problem == 'MOP3' ) then 
	
		! Number of variables
	
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = - pi
		u(:) =   pi

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	!  MOP5
	!  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'MOP5' ) then 

		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		strconvex(3) = .false.
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!  MOP6
	!  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'MOP6' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			
	
! ----------------------------------------------------------------------

	!  MOP7
	!  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'MOP7' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		strconvex(3) = .true.
		
		! Box constraints

		l(:) = - 4.0d2
		u(:) =   4.0d2

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!  PNR
	!	 M. Preuss, B. Naujoks, and G. Rudolph, Pareto set and EMOA behaviour for simple
	!  multimodal multiobjective functions, In: Parallel Problem O Solving from Nature-PPSN IX.
	!  Springer. Berlin, 2006, pp. 513–522.
	
	if ( problem == 'PNR' ) then 

		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = - 2.0d0
		u(:) =   2.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
				
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if				
	
! ----------------------------------------------------------------------	
	
	!   QV1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'QV1' ) then 
	
		! Number of variables
		
		n = 10
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = - 5.0d0
		u(:) =   5.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------	

	! SD
	! Stadler, W., Dauer, J.: Multicriteria optimization in engineering: a tutorial and survey. In: Kamat, M.P.
    ! (ed.) Progress in Aeronautics and Astronautics: Structural Optimization: Status and Promise, vol. 150,
    ! pp. 209–249. American Institute of Aeronautics and Astronautics, Reston (1992)
 
	
	if ( problem == 'SD' ) then 
	
		! Number of variables
		
		n = 4
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		
		! Box constraints

		l(1) = 1.0d0
		u(1) = 3.0d0
		
		l(2:3) = sqrt(2.0d0)
		u(2:3) = 3.0d0
		
		l(4) = 1.0d0
		u(4) = 3.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
	end if		
	
! ----------------------------------------------------------------------	

	! SLCDT1
	! Convergence of stochastic search algorithms to finite size pareto set approximations
 
	
	if ( problem == 'SLCDT1' ) then 
	
		! Number of variables
		
		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = - 1.5d0
		u(:) =   1.5d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
	end if	
	
! ----------------------------------------------------------------------

	!  SLCDT2
	!	 Convergence of stochastic search algorithms to finite size pareto set approximations

	
	if ( problem == 'SLCDT2' ) then 
	
		! Number of variables

		n = 10
		
		! Number of objectives
		
		m = 3
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		strconvex(3) = .false.
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
	end if			
	
	
! ----------------------------------------------------------------------

	!   SP1
	!  	A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'SP1' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = - 1.0d2
		u(:) =   1.0d2

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!   SSFYY2
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'SSFYY2' ) then 
	
		! Number of variables
	
		n = 1
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = - 1.0d2
		u(:) =   1.0d2

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!   SK1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'SK1' ) then 

		! Number of variables
	
		n = 1
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = - 1.0d2
		u(:) =   1.0d2

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	!   SK2
	!  	A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'SK2' ) then 
	
		! Number of variables

		n = 4
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = - 1.0d1
		u(:) =   1.0d1

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	!   TKLY1
	!  	A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'TKLY1' ) then 
	
		! Number of variables

		n = 4
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(1) = 0.1d0
		u(1) = 1.0d0
		
		l(2:4) = 0.0d0
		u(2:4) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	! 	Toi4
	
	! Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	! Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)
	
	if ( problem == 'Toi4' ) then 
	
		! Number of variables

		n = 4
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
	
		! Box constraints

		l(:) = - 2.0d0
		u(:) =   5.0d0
		
		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	! 	Toi8 (TRIDIA)
	
	! Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	! Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)
	
	if ( problem == 'Toi8' ) then 
	
		! Number of variables

		n = 3
		
		! Number of objectives
		
		m = n
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0
		
		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
		
! ----------------------------------------------------------------------

	! 	Toi9 (Shifted TRIDIA)
	
	! Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	! Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)
	
	if ( problem == 'Toi9' ) then 
	
		! Number of variables

		n = 4
												
		! Number of objectives
		
		m = n
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints

		l(:) = - 1.0d0
		u(:) =   1.0d0
		
		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	! 	Toi10 (Rosenbrock)
	
	! Toint, P.L.: Test problems for partially separable optimization and results for the routine PSPMIN.
	! Technical Report, The University of Namur, Department of Mathematics, Belgium (1983)
	! See also: Mita,K. ,Fukuda,E.H. ,Yamashita,N.: Nonmonotone linesearches for unconstrained multiobjective 
	! optimization problems. J. Glob. Optim. 75(1), 63–90 (2019)
	
	if ( problem == 'Toi10' ) then 
	
		! Number of variables

		n = 4
										
		! Number of objectives
		
		m = n - 1
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		do i = 1,m
			strconvex(i) = .false.
		end do
		
		! Box constraints

		l(:) = - 2.0d0
		u(:) =   2.0d0
		
		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if					
	
! ----------------------------------------------------------------------

	!   VU1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'VU1' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = - 3.0d0
		u(:) =   3.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!   VU2
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'VU2' ) then 
	
		! Number of variables

		n = 2
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .true.
		
		! Box constraints

		l(:) = - 3.0d0
		u(:) =   3.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if
	
! ----------------------------------------------------------------------

	!   ZDT1
	!   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( problem == 'ZDT1' ) then 
	
		! Number of variables

		n = 30
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 1.0d-2
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 	
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if	
	
! ----------------------------------------------------------------------

	!   ZDT2
	!   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( problem == 'ZDT2' ) then 
	
		! Number of variables

		n = 30
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!   ZDT3
	!   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( problem == 'ZDT3' ) then 
	
		! Number of variables

		n = 30
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 1.0d-2
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			
	
! ----------------------------------------------------------------------

	!   ZDT4
	!   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( problem == 'ZDT4' ) then 
	
		! Number of variables

		n = 30
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(1) = 1.0d-2
		u(1) = 1.0d0
		
		l(2:n) = - 5.0d0
		u(2:n) =   5.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if		
	
! ----------------------------------------------------------------------

	!   ZDT6
	!   Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
	
	if ( problem == 'ZDT6' ) then 
	
		! Number of variables

		n = 10
		
		! Number of objectives
		
		m = 2
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .false.
		strconvex(2) = .false.
		
		! Box constraints

		l(:) = 0.0d0
		u(:) = 1.0d0

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .true.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			
	
! ----------------------------------------------------------------------

	!   ZLT1
	!   A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'ZLT1' ) then 
	
		! Number of variables

		n = 10
		
		! Number of objectives
		
		m = 5
		
		allocate(x(n),l(n),u(n),strconvex(m),stat=allocerr)
		if ( allocerr .ne. 0 ) then
			 write(*,*) 'Allocation error in main program'
			 stop
		end if
		
		! Strictly strconvex?
		
		strconvex(1) = .true.
		strconvex(2) = .true.
		strconvex(3) = .true.
		strconvex(4) = .true.
		strconvex(5) = .true.
		
		! Box constraints

		l(:) = - 1.0d3
		u(:) =   1.0d3

		! Initial point

		do i = 1,n
			x(i) = l(i) + ( u(i) - l(i) ) * drand(seed)
		end do 
		
		! Penalize the problem?
		
		penalize = .false.
		
		if ( penalize ) then
			lin => l
			uin => u
		end if 
		
		! Scale the problem?
			
		scaleF   = .true.

		! Check derivatives?

		checkder = .false.
		
		return
	end if			
		
	end subroutine inip

	!***********************************************************************
	!***********************************************************************

	subroutine evalf(n,x,f,ind)
	
	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in) :: n,ind
	real(kind=8), intent(out) :: f
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: x(n)
	
	! LOCAL SCALARS
	integer :: i,j,k,m
	real, parameter :: pi = 3.141592653589793
	real(kind=8) :: faux,A1,A2,A3,B1,B2,a,b,t,y
	
	! LOCAL ARRAYS
	real(kind=8) :: MM(3,3),p(3)
	real(kind=8),allocatable :: theta(:) 
	
! ----------------------------------------------------------------------

  ! AP1: Exemple 1 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP1' ) then 
			
		if ( ind == 1 ) then
			f = 0.25d0 * ( ( x(1) - 1.0d0 ) ** 4 + 2.0d0 * ( x(2) - 2.0d0 ) ** 4 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = exp( ( x(1) + x(2) ) / 2.0d0 ) + x(1) ** 2 + x(2) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
		
		if ( ind == 3 ) then
			f = 1.0d0/6.0d0 * ( exp( - x(1) ) + 2.0d0 * exp( - x(2) ) )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------

  ! AP2: Exemple 2 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP2' ) then 
			
		if ( ind == 1 ) then
			f = x(1) ** 2 - 4.0d0
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(1) - 1.0d0 ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------

  ! AP3: Exemple 3 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP3' ) then 
			
		if ( ind == 1 ) then
			f = 0.25d0 * ( ( x(1) - 1.0d0 ) ** 4 + 2.0d0 * ( x(2) - 2.0d0 ) ** 4 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(2) - x(1) ** 2 ) ** 2 + ( 1.0d0 - x(1) ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

 ! AP4: Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP4' ) then 
			
		if ( ind == 1 ) then
			f = 1.0d0/9.0d0 * ( ( x(1) - 1.0d0 ) ** 4 + 2.0d0 * ( x(2) - 2.0d0 ) ** 4 &
			+ 3.0d0 * ( x(3) - 3.0d0 ) ** 4 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + x(1) ** 2 + x(2) ** 2 + x(3) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
		
		if ( ind == 3 ) then
			f = 1.0d0/1.2d1 * ( 3.0d0 * exp( -x(1) ) + 4.0d0 * exp( - x(2) ) + 3.0d0 * exp( -x(3) ) )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------	

	!  BK1
	!  A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit
	
	if ( problem == 'BK1' ) then 
		
		if ( ind == 1 ) then
			f = x(1) ** 2 + x(2) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(1) - 5.0d0 ) ** 2 + ( x(2) - 5.0d0 ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
					
	end if			
	
	
! ----------------------------------------------------------------------	

!  DD1

	if ( problem == 'DD1' ) then 
		
		if ( ind == 1 ) then
			f = x(1) ** 2 + x(2) ** 2 + x(3) ** 2 + x(4) ** 2 + x(5) ** 2 
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 3.0d0 * x(1) + 2.0d0 * x(2) - x(3) / 3.0d0 + 1.0d-2 * ( x(4) - x(5) ) ** 3
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if
	
! ----------------------------------------------------------------------

	!  DGO1

	if ( problem == 'DGO1' ) then 
			
		if ( ind == 1 ) then
			f =  sin( x(1) )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f =  sin( x(1) + 0.7d0 )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
	end if		
	
! ----------------------------------------------------------------------

	!  DGO2

	if ( problem == 'DGO2' ) then 
			
		if ( ind == 1 ) then
			f = x(1) ** 2  
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f =  9.0d0 - sqrt( 8.1d1 - x(1) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
	end if		
	
		
	
! ----------------------------------------------------------------------

	!  DTLZ1

	if ( problem == 'DTLZ1' ) then 
	
		k = dimk
		m = dimm
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
		end do
		faux = 1.0d2 * ( k + faux )
		
		f = 0.5d0 * ( 1.0d0 + faux )					
		do i = 1,m-ind
			f = f * x(i)
		end do
		
		if ( ind > 1 ) f = f * ( 1.0d0 - x(m-ind+1) )
				
		if ( penalize ) call penalizef(n,x,f)
		return	
		
	end if		
	
! ----------------------------------------------------------------------

	!  DTLZ2

	if ( problem == 'DTLZ2' ) then 
	
		k = dimk
		m = dimm
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2
		end do
		
		f = 1.0d0 + faux					
		do i = 1,m-ind
			f = f * cos( x(i) * pi / 2.0d0 )
		end do
		
		if ( ind > 1 ) f = f * sin( x(m-ind+1) * pi / 2.0d0 )
				
		if ( penalize ) call penalizef(n,x,f)
		return	
		
	end if
	
! ----------------------------------------------------------------------

	!  DTLZ3

	if ( problem == 'DTLZ3' ) then 
	
		k = dimk
		m = dimm
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
		end do
		faux = 1.0d2 * ( k + faux )
		
		f = 1.0d0 + faux					
		do i = 1,m-ind
			f = f * cos( x(i) * pi / 2.0d0 )
		end do
		
		if ( ind > 1 ) f = f * sin( x(m-ind+1) * pi / 2.0d0 )
				
		if ( penalize ) call penalizef(n,x,f)
		return	
		
	end if	
	
! ----------------------------------------------------------------------

	!  DTLZ4

	if ( problem == 'DTLZ4' ) then 
	
		k = dimk
		m = dimm
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2
		end do
		
		f = 1.0d0 + faux					
		do i = 1,m-ind
			f = f * cos( x(i)**alpha * pi / 2.0d0 )
		end do
		
		if ( ind > 1 ) f = f * sin( x(m-ind+1)**alpha * pi / 2.0d0 )
				
		if ( penalize ) call penalizef(n,x,f)
		return	
		
	end if	

! ----------------------------------------------------------------------

	!  FA1

	if ( problem == 'FA1' ) then 
			
		if ( ind == 1 ) then
			f = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) )
			f = ( x(2) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(2) + 1.0d0  ) ) ** 0.5d0 )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
		
		if ( ind == 3 ) then
			faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) )
			f = ( x(3) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(3) + 1.0d0  ) ) ** 0.1d0 )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
	end if		
	
! ----------------------------------------------------------------------

	!  Far1

	if ( problem == 'Far1' ) then 
			
		if ( ind == 1 ) then
			f =  - 2.0d0 * exp( 1.5d1 * ( - ( x(1) - 0.1d0 ) ** 2 - x(2) ** 2 ) ) &
				 - exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				 + exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				 + exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ** 2 - ( x(2) + 0.6d0 ) ** 2 ) ) &
				 + exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ** 2 - ( x(2) + 0.6d0 ) ** 2 ) ) 
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f =  2.0d0 * exp( 2.0d1 * ( - x(1) ** 2 - x(2) ** 2 ) ) &
				 + exp( 2.0d1 * ( - ( x(1) - 0.4d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				 - exp( 2.0d1 * ( - ( x(1) + 0.5d0 ) ** 2 - ( x(2) - 0.7d0 ) ** 2 ) ) &
				 - exp( 2.0d1 * ( - ( x(1) - 0.5d0 ) ** 2 - ( x(2) + 0.7d0 ) ** 2 ) ) &
				 + exp( 2.0d1 * ( - ( x(1) + 0.4d0 ) ** 2 - ( x(2) + 0.8d0 ) ** 2 ) ) 
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if
	
! ----------------------------------------------------------------------

	! FDS
	! NEWTON’S METHOD FOR MULTIOBJECTIVE OPTIMIZATION
			
	if ( problem == 'FDS' ) then 
			
		if ( ind == 1 ) then
			f = 0.0d0
			do i = 1,n
				f = f + i * ( x(i) - i ) ** 4
			end do	
			f = f / ( n ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = exp( sum(x)/n ) + norm2(x) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
		
		if ( ind == 3 ) then
			f = 0.0d0
			do i = 1,n
				f = f + i * ( n - i + 1.0d0 ) * exp( - x(i) )
			end do	
			f = f / ( n * ( n + 1.0d0 ) ) 
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------
	
	! FF1 
	! C. M. Fonseca and P. J. Fleming, “An overview of evolutionary algorithms in multiobjective optimization

	if ( problem == 'FF1' ) then 

		if ( ind == 1 ) then
			f = 1.0d0 - exp( - ( x(1) - 1.0d0 ) ** 2 - ( x(2) + 1.0d0 ) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 1.0d0 - exp( - ( x(1) + 1.0d0 ) ** 2 - ( x(2) - 1.0d0 ) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------
	
	!  Hil1

	if ( problem == 'Hil1' ) then 
		a = 2.0d0 * pi / 3.6d2 * ( 4.5d1 + 4.0d1 * sin( 2.0d0 * pi * x(1) ) &
			+ 2.5d1 * sin( 2.0d0 * pi * x(2) ) )
		b = 1.0d0 + 0.5d0 * cos( 2.0d0 * pi * x(1) )
		
		if ( ind == 1 ) then
			f = cos( a ) * b
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = sin( a ) * b
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if			
	
! ----------------------------------------------------------------------
	
	!  IKK1

	if ( problem == 'IKK1' ) then 
		
		if ( ind == 1 ) then
			f = x(1) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(1) - 2.0d1 ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
		
		if ( ind == 3 ) then
			f = x(2) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if			
	
! ----------------------------------------------------------------------
	
	!  IM1

	if ( problem == 'IM1' ) then 
		
		if ( ind == 1 ) then
			f = 2.0d0 * sqrt( x(1) )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = x(1) * ( 1.0d0 - x(2) ) + 5.0d0
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if
	
! ----------------------------------------------------------------------	
	
	! JOS1 
	! Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why Does It Work and How?
	
	if ( problem == 'JOS1' ) then 

		if ( ind == 1 ) then
			f = 0.0d0
			do i = 1,n
					f = f + x(i) ** 2
			end do
			f = f / n
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 0.0d0
			do i = 1,n
					f = f + ( x(i) - 2.0d0 ) ** 2  
			end do
			f = f / n
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
	
	end if
	
! ----------------------------------------------------------------------	
	
	! JOS4
	! Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why Does It Work and How?
	
	if ( problem == 'JOS4' ) then 

		if ( ind == 1 ) then
		
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
		
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			
			f = faux * ( 1.0d0 - ( x(1) / faux ) ** 0.25d0 - ( x(1) / faux ) ** 4.0d0 )
			
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
	
	end if	
	
! ----------------------------------------------------------------------

	! 	KW2
	
	if ( problem == 'KW2' ) then 
			
		if ( ind == 1 ) then
			f = - 3.0d0 * ( 1.0d0 - x(1) )**2 * exp( -x(1)**2 - ( x(2) + 1.0d0 ) ** 2 ) &
				+ 1.0d1 * ( x(1) / 5.0d0 - x(1)**3 - x(2)**5 ) * exp( - x(1)**2 - x(2)**2 ) &
				+ 3.0d0 * exp( -( x(1) + 2.0d0 )**2 - x(2)**2 ) - 0.5d0 * ( 2.0d0 * x(1) + x(2) )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
			
		if ( ind == 2 ) then
			f = - 3.0d0 * ( 1.0d0 + x(2) )**2 * exp( -x(2)**2 - ( 1.0d0 - x(1) ) ** 2 ) &
				+ 1.0d1 * ( - x(2) / 5.0d0 + x(2)**3 + x(1)**5 ) * exp( - x(1)**2 - x(2)**2 ) &
				+ 3.0d0 * exp( -( 2.0d0 - x(2) )**2 - x(1)**2 ) 
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	! 	LE1
	
	if ( problem == 'LE1' ) then 
											
		if ( ind == 1 ) then
			f = ( x(1) ** 2 + x(2) ** 2 ) ** 0.125d0
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( ( x(1) - 0.5d0 ) ** 2 + ( x(2) - 0.5d0 ) ** 2 ) ** 0.25d0
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
			
	end if		
	
! ----------------------------------------------------------------------
	
	! Lov1  
	
	if ( problem == 'Lov1' ) then 

		if ( ind == 1 ) then
			f = - ( -1.05d0 * x(1) ** 2 - 0.98d0 * x(2) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = - ( -0.99d0 * ( x(1) - 3.0d0 ) ** 2 - 1.03d0 * ( x(2) - 2.5d0 ) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------
	
	! Lov2
	
	if ( problem == 'Lov2' ) then 

		if ( ind == 1 ) then
			f = x(2)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(2) - x(1) ** 3 ) / ( x(1) + 1.0d0 )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------
	
	! Lov3  
	
	if ( problem == 'Lov3' ) then 
			
		if ( ind == 1 ) then
				f = - ( - x(1) ** 2 - x(2) ** 2 )
				if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
				f = - ( - ( x(1) - 6.0d0 ) ** 2 + ( x(2) + 0.3d0 ) ** 2 )
				if ( penalize ) call penalizef(n,x,f)
			return
		end if	
		
	end if		

! ----------------------------------------------------------------------

	! Lov4  
	
	if ( problem == 'Lov4' ) then 
	
		if ( ind == 1 ) then
			f = - x(1) ** 2 - x(2) ** 2 - 4.0d0 * ( exp( - ( x(1) + 2.0d0 ) ** 2 - x(2) ** 2 ) + &
				exp( - ( x(1) - 2.0d0 ) ** 2 - x(2) ** 2 ) )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = - ( x(1) - 6.0d0 ) ** 2 - ( x(2) + 0.5d0 ) ** 2
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

	! Lov5
	
	if ( problem == 'Lov5' ) then 
	
		MM(1,:) = (/ -1.0d0  , -0.03d0,  0.011d0/)
		MM(2,:) = (/ -0.03d0 , -1.0d0 ,  0.07d0 /)
		MM(3,:) = (/  0.011d0,  0.07d0, -1.01d0 /)
		
		p(:) = (/ x(1)  , x(2) - 0.15d0,  x(3)/)
		a = 0.35d0
		
		A1 = sqrt( 2.0d0 * pi / a ) * exp( dot_product( p, matmul(MM,p) ) / a ** 2 )
		
		p(:) = (/ x(1)  , x(2) + 1.1d0,  0.5d0 * x(3)/)
		a = 3.0d0
		
		A2 = sqrt( 2.0d0 * pi / a ) * exp( dot_product( p, matmul(MM,p) ) / a ** 2 )
		
		faux = A1 + A2
	
		if ( ind == 1 ) then
			f = sqrt(2.0d0)/2.0d0 * ( x(1) + faux )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = sqrt(2.0d0)/2.0d0 * ( - x(1) + faux )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------

	! Lov6
	
	if ( problem == 'Lov6' ) then 
	
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
		end if
		
		if ( ind == 2 ) then
			f = 1.0d0 - sqrt( x(1) ) - x(1) * sin( 1.0d1 * pi * x(1) ) &
			+ x(2) ** 2 + x(3) ** 2 + x(4) ** 2 + x(5) ** 2 + x(6) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

	! 	LTDZ
	!	Combining convergence and diversity in evolutionary multiobjective optimization
	
	if ( problem == 'LTDZ' ) then 
											
		if ( ind == 1 ) then
			f = 3.0d0 - ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 3.0d0 - ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 3 ) then
			f = 3.0d0 - ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * sin( x(1) * pi / 2.0d0 )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
			
	end if						

! ----------------------------------------------------------------------

	! 	MGH9 
	
	if ( problem == 'MGH9' ) then 
	
		if ( ind == 1 .or. ind == 15 ) then
			y = 9.0d-4
		elseif ( ind == 2 .or. ind == 14 ) then 		
			y = 4.4d-3
		elseif ( ind == 3 .or. ind == 13 ) then 		
			y = 1.75d-2
		elseif ( ind == 4 .or. ind == 12 ) then 		
			y = 5.4d-2		
		elseif ( ind == 5 .or. ind == 11 ) then 		
			y = 1.295d-1		
		elseif ( ind == 6 .or. ind == 10 ) then 		
			y = 2.42d-1
		elseif ( ind == 7 .or. ind == 9  ) then 		
			y = 3.521d-1
		elseif ( ind == 8 ) then 		
			y = 3.989d-1
		end if
		
		t = ( 8.0d0 - ind ) / 2.0d0
				
		f = x(1) * exp( - x(2) * ( t - x(3) ) ** 2 / 2.0d0 ) - y			
		
		if ( penalize ) call penalizef(n,x,f)
		return
		
	end if	
	
! ----------------------------------------------------------------------

	! 	MGH16 
	
	if ( problem == 'MGH16' ) then 
		
		t = ind / 5.0d0
				
		f = ( x(1) + t * x(2) - exp(t) ) ** 2 + ( x(3) + x(4) * sin(t) - cos(t) ) ** 2
		
		if ( penalize ) call penalizef(n,x,f)
		return
			
	end if				
	
! ----------------------------------------------------------------------

	! 	MGH26 
	
	if ( problem == 'MGH26' ) then 
								
		t = 0.0d0
		do i = 1,n
			t = t + cos(x(i))
		end do
		
		f = ( n - t + ind * ( 1.0d0 - cos(x(ind)) ) - sin(x(ind)) ) ** 2
		
		if ( penalize ) call penalizef(n,x,f)
		return
			
	end if				
	
! ----------------------------------------------------------------------

	! 	MGH33
	
	if ( problem == 'MGH33' ) then 
								
		faux = 0.0d0					
		do i = 1,n
			faux = faux + i * x(i)
		end do
		
		f = ( ind * faux - 1.0d0 ) ** 2			
		
		if ( penalize ) call penalizef(n,x,f)
		return
			
	end if	
	
	! ----------------------------------------------------------------------

	! 	MHHM2
	
	if ( problem == 'MHHM2' ) then 
											
		if ( ind == 1 ) then
			f = ( x(1) - 0.8d0 ) ** 2 + ( x(2) - 0.6d0 ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(1) - 0.85d0 ) ** 2 + ( x(2) - 0.7d0 ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
			
		if ( ind == 3 ) then
			f = ( x(1) - 0.9d0 ) ** 2 + ( x(2) - 0.6d0 ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
				
	end if		
	
! ----------------------------------------------------------------------

	!  MLF1

	if ( problem == 'MLF1' ) then 
			
		if ( ind == 1 ) then
			f = ( 1.0d0 + x(1) / 2.0d1 ) * sin( x(1) )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( 1.0d0 + x(1) / 2.0d1 ) * cos( x(1) )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	!  MLF2

	if ( problem == 'MLF2' ) then 
			
		if ( ind == 1 ) then
			f = 5.0d0 - ( ( x(1) ** 2 + x(2) - 1.1d1 ) ** 2 &
				+ ( x(1) + x(2) ** 2 - 7.0d0 ) ** 2 ) / 2.0d2
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 5.0d0 - ( ( 4.0d0 * x(1) ** 2 + 2.0d0 * x(2) - 1.1d1 ) ** 2 &
				+ ( 2.0d0 * x(1) + 4.0d0 *  x(2) ** 2 - 7.0d0 ) ** 2 ) / 2.0d2
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	!  MMR1 		
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	

!	if ( problem == 'MMR1' ) then 
			
!		if ( ind == 1 ) then
!			f = 1.0d0 + x(1) ** 2
!			if ( penalize ) call penalizef(n,x,f)
!			return
!		end if
		
!		if ( ind == 2 ) then
!			f = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ** 2 ) -&
!			 exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ** 2 )
!			f = f / ( 1.0d0 + x(1) ** 2 )
!			if ( penalize ) call penalizef(n,x,f)
!			return	
!		end if	
			
!	end if		
	
	if ( problem == 'MMR1' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ** 2 ) -&
			 exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ** 2 )
			f = f / x(1)
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if	
			
	end if
	
! ----------------------------------------------------------------------

	!  MMR2
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	

	if ( problem == 'MMR2' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then			
			faux = x(1) / ( 1.0d0 + 1.0d1 * x(2) )
			f = 1.0d0 - faux ** 2 - faux * sin( 8.0d0 * pi * x(1) )
			f = f * ( 1.0d0 + 1.0d1 * x(2) )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if	
			
	end if		
		
! ----------------------------------------------------------------------

	!  MMR3
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	

	if ( problem == 'MMR3' ) then 
			
		if ( ind == 1 ) then
			f = x(1) ** 3
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(2) - x(1) ) ** 3
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if	
			
	end if		
	
! ----------------------------------------------------------------------

	!  MMR4
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	

	if ( problem == 'MMR4' ) then 
			
		if ( ind == 1 ) then
			f = x(1) - 2.0d0 * x(2) - x(3) - 3.6d1 / ( 2.0d0 * x(1) + x(2) + 2.0d0 * x(3) + 1.0d0 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = - 3.0d0 * x(1) + x(2) - x(3)
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if	
			
	end if						

! ----------------------------------------------------------------------

	!  MOP 2	

	if ( problem == 'MOP2' ) then 
			
		if ( ind == 1 ) then
			f = 0.0d0
			do i = 1,n
				f = f + ( x(i) - 1.0d0 / ( sqrt( real(n) ) ) ) ** 2
			end do
			
			f = 1.0d0 - exp ( - f )	
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 0.0d0
			do i = 1,n
				f = f + ( x(i) + 1.0d0 / ( sqrt( real (n) ) ) ) ** 2
			end do
			
			f = 1.0d0 - exp ( - f )	
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

	!  MOP 3	

	if ( problem == 'MOP3' ) then 
			
		if ( ind == 1 ) then
			A1 = 0.5d0 * sin(1.0d0) - 2.0d0 * cos(1.0d0) + sin(2.0d0) - 1.5d0 * cos(2.0d0) 
			A2 = 1.5d0 * sin(1.0d0) - cos(1.0d0) + 2.0d0 * sin(2.0d0) - 0.5d0 * cos(2.0d0)
			B1 = 0.5d0 * sin(x(1)) - 2.0d0 * cos(x(1)) + sin(x(2)) - 1.5d0 * cos(x(2)) 
			B2 = 1.5d0 * sin(x(1)) - cos(x(1)) + 2.0d0 * sin(x(2)) - 0.5d0 * cos(x(2))	
			f = - 1.0d0 - ( A1 - B1 ) ** 2 - ( A2 - B2 ) ** 2
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = - ( x(1) + 3.0d0 ) ** 2 - ( x(2) + 1.0d0 ) ** 2
			f = - f	
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if	
			
	end if			
	
! ----------------------------------------------------------------------

	!  MOP 5	

	if ( problem == 'MOP5' ) then 
			
		if ( ind == 1 ) then
			f = 0.5d0 * ( x(1) ** 2 + x(2) ** 2 ) + sin( x(1) ** 2 + x(2) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( 3.0d0 * x(1) - 2.0d0 * x(2) + 4.0d0 ) ** 2 / 8.0d0 + &
				( x(1) - x(2) + 1.0d0 ) ** 2 / 2.7d1 + 1.5d1
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
		
		if ( ind == 3 ) then
			f = 1.0d0 / ( x(1) ** 2 + x(2) ** 2 + 1.0d0 ) - 1.1d0 * exp( - x(1) ** 2 - x(2) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

	!  MOP6

	if ( problem == 'MOP6' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			a = 1.0d0 + 1.0d1 * x(2) 
			t = x(1) / a
			f = a * ( 1.0d0 - t ** 2 - t * sin( 8.0d0 * pi * x(1) ) )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if	
	

! ----------------------------------------------------------------------

	!  MOP 7

	if ( problem == 'MOP7' ) then 
			
		if ( ind == 1 ) then
			f =  ( x(1) - 2.0d0 ) ** 2 / 2.0d0 &
			+ ( x(2) + 1.0d0 ) ** 2 / 1.3d1 + 3.0d0
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f =  ( x(1) + x(2) - 3.0d0 ) ** 2 / 3.6d1 &
			+ ( - x(1) + x(2) + 2.0d0 ) ** 2 / 8.0d0 - 1.7d1
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
		
		if ( ind == 3 ) then
			f =  ( x(1) + 2.0d0 * x(2) - 1.0d0 ) ** 2 / 1.75d2 &
			+ ( - x(1) + 2.0d0 * x(2) ) ** 2 / 1.7d1 - 1.3d1
			if ( penalize ) call penalizef(n,x,f)
			return
		end if	
				
	end if			
	
! ----------------------------------------------------------------------

	! 	PNR
	
	if ( problem == 'PNR' ) then 
											
		if ( ind == 1 ) then
			f = x(1) ** 4 + x(2) ** 4 - x(1) ** 2 + x(2) ** 2 - 1.0d1 * x(1) * x(2) + 2.0d1
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = x(1) ** 2 + x(2) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
			
	end if
	
! ----------------------------------------------------------------------

	!  QV1

	if ( problem == 'QV1' ) then 
			
		if ( ind == 1 ) then
			f = 0.0d0
			do i = 1,n
				f = f + x(i) ** 2 - 1.0d1 * cos( 2.0d0 * pi * x(i) ) + 1.0d1
			end do
		
			f = ( f / n ) ** 0.25 
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 0.0d0
			do i = 1,n
				f = f + ( x(i) - 1.5d0 ) ** 2 - 1.0d1 * cos( 2.0d0 * pi * ( x(i) -1.5d0 ) ) + 1.0d1
			end do
		
			f = ( f / n ) ** 0.25
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	! SD

	if ( problem == 'SD' ) then 
			
		if ( ind == 1 ) then
			f = 2.0d0 * x(1) + sqrt(2.0d0) * ( x(2) + x(3) ) + x(4)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 2.0d0 / x(1) + 2.0d0 * sqrt(2.0d0) / x(2) + 2.0d0 * sqrt(2.0d0) / x(3) + 2.0d0 / x(4)
			if ( penalize ) call penalizef(n,x,f)
			return		
		end if	
			
	end if		
	
! ----------------------------------------------------------------------

	! SLCDT1

	if ( problem == 'SLCDT1' ) then 
			
		if ( ind == 1 ) then
			f = 0.5d0 * ( sqrt( 1.0d0 + ( x(1) + x(2) ) ** 2 ) + &
					sqrt( 1.0d0 + ( x(1) - x(2) ) ** 2 ) + x(1) - x(2) ) +  &
					0.85d0 * exp( - ( x(1) + x(2) ) ** 2 )		
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 0.5d0 * ( sqrt( 1.0d0 + ( x(1) + x(2) ) ** 2 ) + &
					sqrt( 1.0d0 + ( x(1) - x(2) ) ** 2 ) - x(1) + x(2) ) +  &
					0.85d0 * exp( - ( x(1) + x(2) ) ** 2 )	
			if ( penalize ) call penalizef(n,x,f)
			return		
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

	! 	SLCDT2
	!	Convergence of stochastic search algorithms to finite size pareto set approximations 
	
	if ( problem == 'SLCDT2' ) then 
										
		if ( ind == 1 ) then
			f = ( x(1) - 1.0d0 ) ** 4
			do i = 2,n
				f = f + ( x(i) - 1.0d0 ) ** 2
			end do
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(2) + 1.0d0 ) ** 4
			do i = 1,n
				if ( i /= 2 ) f = f + ( x(i) + 1.0d0 ) ** 2  
			end do
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 3 ) then
			f = ( x(3) - 1.0d0 ) ** 4
			do i = 1,n
				if ( i /= 3 ) f = f + ( x(i) - ( - 1.0d0 ) ** (i+1) ) ** 2  
			end do
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
			
	end if		
						
! ----------------------------------------------------------------------

	!  SP1

	if ( problem == 'SP1' ) then 
			
		if ( ind == 1 ) then
			f = ( x(1) - 1.0d0 ) ** 2 + ( x(1) - x(2) ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(2) - 3.0d0 ) ** 2 + ( x(1) - x(2) ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if				
		
! ----------------------------------------------------------------------

	!  SSFYY2

	if ( problem == 'SSFYY2' ) then 
		
		if ( ind == 1 ) then
			f = 1.0d1 + x(1) ** 2 - 1.0d1 * cos( x(1) * pi / 2.0d0 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( x(1) - 4.0d0 ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if			
				
! ----------------------------------------------------------------------

	!  SK1

	if ( problem == 'SK1' ) then 
			
		if ( ind == 1 ) then
			f = - x(1) ** 4 - 3.0d0 * x(1) ** 3 + 1.0d1 * x(1) ** 2 + 1.0d1 * x(1) + 1.0d1
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = - 0.5d0 * x(1) ** 4 + 2.0d0 * x(1) ** 3 + 1.0d1 * x(1) ** 2 - 1.0d1 * x(1) + 5.0d0
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	!  SK2

	if ( problem == 'SK2' ) then 
			
		if ( ind == 1 ) then
			f = - ( x(1) - 2.0d0 ) ** 2 - ( x(2) + 3.0d0 ) ** 2 &
			- ( x(3) - 5.0d0 ) ** 2 - ( x(4) - 4.0d0 ) ** 2 + 5.0d0
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = ( sin( x(1) ) + sin( x(2) ) + sin( x(3) ) + sin( x(4) ) ) / &
			( 1.0d0 + ( x(1) ** 2 + x(2) ** 2 + x(3) ** 2 + x(4) ** 2 ) / 1.0d2 )
			f = - f
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if		
	
	
! ----------------------------------------------------------------------

	!  TKLY1

	if ( problem == 'TKLY1' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
		
			A1 = ( 2.0d0 - exp( - ( ( x(2) - 0.1d0 ) / 4.0d-3 ) ** 2 ) &
				         - 0.8d0 * exp( - ( ( x(2) - 0.9d0 ) / 4.0d-1 ) ** 2 ) )
			A2 = ( 2.0d0 - exp( - ( ( x(3) - 0.1d0 ) / 4.0d-3 ) ** 2 ) &
				         - 0.8d0 * exp( - ( ( x(3) - 0.9d0 ) / 4.0d-1 ) ** 2 ) )
			A3 = ( 2.0d0 - exp( - ( ( x(4) - 0.1d0 ) / 4.0d-3 ) ** 2 ) &
				         - 0.8d0 * exp( - ( ( x(4) - 0.9d0 ) / 4.0d-1 ) ** 2 ) )
			f = A1 * A2 * A3 / x(1)
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if	
	
	
! ----------------------------------------------------------------------

	! 	Toi4
	
	if ( problem == 'Toi4' ) then 
			
		if ( ind == 1 ) then
			f = x(1) ** 2 + x(2) ** 2 + 1.0d0
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = 0.5d0 * ( ( x(1) - x(2) ) ** 2 + ( x(3) - x(4) ) ** 2  ) + 1.0d0
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if					
	
! ----------------------------------------------------------------------

	! 	Toi8
	
	if ( problem == 'Toi8' ) then 
			
		if ( ind == 1 ) then
			f = ( 2.0d0 * x(1) - 1.0d0 ) ** 2 
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind /= 1 ) then
			f = ind * ( 2.0d0 * x(ind-1) - x(ind) ) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	! 	Toi9
	
	if ( problem == 'Toi9' ) then 
			
		if ( ind == 1 ) then
			f = ( 2.0d0 * x(1) - 1.0d0 ) ** 2 + x(2) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind > 1 .and. ind < n ) then
			f = ind * ( 2.0d0 * x(ind-1) - x(ind) ) ** 2  &
			   - ( ind - 1.0d0 ) * x(ind-1)**2 + ind * x(ind) ** 2
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
		
		if ( ind == n ) then
			f = n * ( 2.0d0 * x(n-1) - x(n) ) ** 2 - ( n - 1.0d0 ) * x(n-1)**2
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	! 	Toi10 (Rosenbrock)
	
	if ( problem == 'Toi10' ) then 
			
		f = 1.0d2 * ( x(ind+1) - x(ind) ** 2 ) ** 2 + ( x(ind+1) - 1.0d0 ) ** 2
		if ( penalize ) call penalizef(n,x,f)
			return
			
	end if	
	
! ----------------------------------------------------------------------

	!  VU1

	if ( problem == 'VU1' ) then 
			
		if ( ind == 1 ) then
			f = 1.0d0 / ( x(1) ** 2 + x(2) ** 2 + 1.0d0 )
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = x(1) ** 2 + 3.0d0 * x(2) ** 2 + 1.0d0
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  VU2

	if ( problem == 'VU2' ) then 
			
		if ( ind == 1 ) then
			f = x(1) + x(2) + 1.0d0
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			f = x(1) ** 2 + 2.0d0 * x(2) - 1.0d0
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	!  ZDT1

	if ( problem == 'ZDT1' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			
			f = faux * ( 1.0d0 - sqrt( x(1) / faux ) )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if			
	
! ----------------------------------------------------------------------

	!  ZDT2

	if ( problem == 'ZDT2' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			
			f = faux * ( 1.0d0 - ( x(1) / faux ) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if				
	
! ----------------------------------------------------------------------

	!  ZDT3

	if ( problem == 'ZDT3' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			t = x(1) / faux
			
			f = faux * ( 1.0d0 - sqrt( t ) - t * sin( 1.0d1 * pi * x(1) ) )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if			
		
! ----------------------------------------------------------------------

	!  ZDT4

	if ( problem == 'ZDT4' ) then 
			
		if ( ind == 1 ) then
			f = x(1)
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 0.0d0
			do i = 2,n
				faux = faux + x(i) ** 2 - 1.0d1 * cos( 4.0d0 * pi * x(i) )
			end do
			faux = faux + 1.0d0 + 1.0d1 * ( n - 1 )
			t = x(1) / faux
			
			f = faux * ( 1.0d0 - sqrt( t ) )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  ZDT6

	if ( problem == 'ZDT6' ) then 
			
		if ( ind == 1 ) then
			f = 1.0d0 - exp( -4.0d0 * x(1) ) * ( sin( 6.0d0 * pi * x(1) ) ) ** 6
			if ( penalize ) call penalizef(n,x,f)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * ( sum(x(2:n)) / ( n - 1 ) ) ** 0.25d0
			
			f = faux * ( 1.0d0 - ( ( 1.0d0 - exp( -4.0d0 * x(1) ) * &
			    ( sin( 6.0d0 * pi * x(1) ) ) ** 6 ) / faux ) ** 2 )
			if ( penalize ) call penalizef(n,x,f)
			return	
		end if
			
	end if			
	
! ----------------------------------------------------------------------

	!  ZLT1

	if ( problem == 'ZLT1' ) then 
		
		f = ( x(ind) - 1.0d0 ) ** 2
		do i = 1,n
			if ( i /= ind ) f = f + x(i) ** 2
		end do
		if ( penalize ) call penalizef(n,x,f)
		return
			
	end if				
			
	end subroutine evalf	

	!***********************************************************************
	!***********************************************************************

	subroutine evalg(n,x,g,ind)
	
	implicit none
	
	! SCALAR ARGUMENTS
	integer, intent(in) :: n,ind
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in)  :: x(n)
	real(kind=8), intent(out) :: g(n)
	
	! LOCAL SCALARS
	
	integer :: i,j,k,l,m
	real, parameter :: pi = 3.141592653589793
	real(kind=8) :: faux,gaux1,gaux2,A1,A2,A3,B1,B2,a,b,t
	
	! LOCAL ARRAYS
	real(kind=8) :: MM(3,3),p(3)
	real(kind=8),allocatable :: theta(:),thetader(:)
	
! ----------------------------------------------------------------------

! AP1: Exemple 1 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP1' ) then 
			
		if ( ind == 1 ) then
			g(1) = ( x(1) - 1.0d0 ) ** 3 
			g(2) = 2.0d0 * ( x(2) - 2.0d0 ) ** 3
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) + 2.0d0 * x(1)
			g(2) = 0.5d0 * exp( ( x(1) + x(2) ) / 2.0d0 ) + 2.0d0 * x(2)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
		
		if ( ind == 3 ) then
			g(1) = - 1.0d0/6.0d0 * exp( - x(1) ) 
			g(2) = - 1.0d0/3.0d0 * exp( - x(2) ) 
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
		
	end if	
	
! ----------------------------------------------------------------------

! AP2: Exemple 2 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP2' ) then 
			
		if ( ind == 1 ) then
			g(1) = 2.0d0 * x(1)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * ( x(1) - 1.0d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if		
	
	
! ----------------------------------------------------------------------

! AP3: Exemple 3 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP3' ) then 
			
		if ( ind == 1 ) then
			g(1) = ( x(1) - 1.0d0 ) ** 3 
			g(2) = 2.0d0 * ( x(2) - 2.0d0 ) ** 3
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = - 4.0d0 * x(1) * ( x(2) - x(1) ** 2 ) - 2.0d0 * ( 1.0d0 - x(1) )
			g(2) = 2.0d0 * ( x(2) - x(1) ** 2 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

! AP4: Exemple 4 of "A modified Quasi-Newton method for vector optimization problem"

	if ( problem == 'AP4' ) then 
		
		if ( ind == 1 ) then
			g(1) = 4.0d0/9.0d0 * ( x(1) - 1.0d0 ) ** 3 
			g(2) = 8.0d0/9.0d0 * ( x(2) - 2.0d0 ) ** 3 
			g(3) = 1.2d1/9.0d0 * ( x(3) - 3.0d0 ) ** 3 
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(1) 
			g(2) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(2) 
			g(3) = 1.0d0/3.0d0 * exp( ( x(1) + x(2) + x(3) ) / 3.0d0 ) + 2.0d0 * x(3) 
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
		
		if ( ind == 3 ) then
			g(1) = - 1.0d0/4.0d0 * exp( -x(1) ) 
			g(2) = - 1.0d0/3.0d0 * exp( -x(2) ) 
			g(3) = - 1.0d0/4.0d0 * exp( -x(3) ) 
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------	

	!  BK1

	if ( problem == 'BK1' ) then 
			
		if ( ind == 1 ) then
			do i = 1,n
				g(i) = 2.0d0 * x(i)
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			do i = 1,n
				g(i) = 2.0d0 * ( x(i) - 5.0d0 ) 
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------	

	!  DD1

	if ( problem == 'DD1' ) then 
			
		if ( ind == 1 ) then
			do i = 1,n
				g(i) = 2.0d0 * x(i)
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 3.0d0
			g(2) = 2.0d0
			g(3) = - 1.0d0 / 3.0d0
			g(4) = 3.0d-2 * ( x(4) - x(5) ) ** 2
			g(5) = - 3.0d-2 * ( x(4) - x(5) ) ** 2
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  DGO1

	if ( problem == 'DGO1' ) then 
			
		if ( ind == 1 ) then
			g(1) =  cos( x(1) )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) =  cos( x(1) + 0.7d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
	end if		
	
! ----------------------------------------------------------------------

	!  DGO2

	if ( problem == 'DGO2' ) then 
			
		if ( ind == 1 ) then
			g(1) =  2.0d0 * x(1)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) =  x(1) / sqrt( 8.1d1 - x(1) ** 2 )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
	end if		
	
! ----------------------------------------------------------------------

	!  DTLZ1

	if ( problem == 'DTLZ1' ) then 
	
!		k = dimk
!		m = dimm
	
!		faux = 0.0d0
!		do i = m,n
!			faux = faux + ( x(i) - 0.5d0 ) ** 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
!		end do
!		faux = 1.0d2 * ( k + faux )
		
!		f = 0.5d0 * ( 1.0d0 + faux )			
!		do i = 1,m-ind
!			f = f * x(i)
!		end do
		
!		if ( ind > 1 ) f = f * ( 1.0d0 - x(m-ind+1) )
		
		
		k = dimk
		m = dimm
		
		do i = 1,n
			g(i) = 0.0d0
		end do
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
		end do
		faux = 1.0d2 * ( k + faux )
		
		do i = 1,m-ind
			g(i) = 0.5d0 * ( 1.0d0 + faux )
		end do 
		
		do i = 1,m-ind
			do j = 1,m-ind
				if ( j == i ) cycle
				g(i) = g(i) * x(j)
			end do
			if ( ind > 1 ) g(i) = g(i) * ( 1.0d0 - x(m-ind+1) )
		end do
		
		if ( ind > 1 ) then
			g(m-ind+1) = - 0.5d0 * ( 1.0d0 + faux )
			do j = 1,m-ind
				g(m-ind+1) = g(m-ind+1) * x(j)
			end do
		end if
				
		
		faux = 0.5d0
		do i = 1,m-ind
			faux = faux * x(i)
		end do

		if ( ind > 1 ) faux = faux * ( 1.0d0 - x(m-ind+1) )
		
		do i = m,n
			g(i) = faux * 1.0d2 * ( 2.0d0 * ( x(i) - 0.5d0 ) + 2.0d1 * pi * sin( 2.0d1 * pi * ( x(i) - 0.5d0 ) ) )
		end do
					
		if ( penalize ) call penalizeg(n,x,g)
		return

	end if
	
! ----------------------------------------------------------------------

	!  DTLZ2

	if ( problem == 'DTLZ2' ) then 
	
!		k = dimk
!		m = dimm
	
!		faux = 0.0d0
!		do i = m,n
!			faux = faux + ( x(i) - 0.5d0 ) ** 2
!		end do
		
!		f = 1.0d0 + faux					
!		do i = 1,m-ind
!			f = f * cos( x(i) * pi / 2.0d0 )
!		end do
		
!		if ( ind > 1 ) f = f * sin( x(m-ind+1) * pi / 2.0d0 )
		
		k = dimk
		m = dimm
		
		do i = 1,n
			g(i) = 0.0d0
		end do
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2
		end do
		faux = 1.0d0 + faux
		
		do i = 1,m-ind
			g(i) = faux
		end do 
		
		do i = 1,m-ind
			do j = 1,m-ind
				if ( j == i ) then
					g(i) = - g(i) * pi / 2.0d0 * sin( x(i) * pi / 2.0d0 )
				else
					g(i) = g(i) * cos( x(j) * pi / 2.0d0 )
				end if
			end do
			if ( ind > 1 ) g(i) = g(i) * sin( x(m-ind+1) * pi / 2.0d0 )
		end do
		
		if ( ind > 1 ) then
			g(m-ind+1) = faux * pi / 2.0d0 * cos( x(m-ind+1) * pi / 2.0d0 )
			do j = 1,m-ind
				g(m-ind+1) = g(m-ind+1) * cos( x(j) * pi / 2.0d0 )
			end do
		end if
				
		
		faux = 1.0d0
		do i = 1,m-ind
			faux = faux * cos( x(i) * pi / 2.0d0 )
		end do

		if ( ind > 1 ) faux = faux * sin( x(m-ind+1) * pi / 2.0d0 )
		
		do i = m,n
			g(i) = faux * 2.0d0 * ( x(i) - 0.5d0 )
		end do
					
		if ( penalize ) call penalizeg(n,x,g)
		return

	end if
	
! ----------------------------------------------------------------------

	!  DTLZ3

	if ( problem == 'DTLZ3' ) then 
	
!		k = dimk
!		m = dimm
	
!		faux = 0.0d0
!		do i = m,n
!			faux = faux + ( x(i) - 0.5d0 ) ** 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
!		end do
!		faux = 1.0d2 * ( k + faux )
		
!		f = 1.0d0 + faux					
!		do i = 1,m-ind
!			f = f * cos( x(i) * pi / 2.0d0 )
!		end do
		
!		if ( ind > 1 ) f = f * sin( x(m-ind+1) * pi / 2.0d0 )
		
		k = dimk
		m = dimm
		
		do i = 1,n
			g(i) = 0.0d0
		end do
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2 - cos( 2.0d1 * pi * ( x(i) - 0.5d0 ) )
		end do
		faux = 1.0d2 * ( k + faux )
		faux = faux + 1.0d0
		
		do i = 1,m-ind
			g(i) = faux
		end do 
		
		do i = 1,m-ind
			do j = 1,m-ind
				if ( j == i ) then
					g(i) = - g(i) * pi / 2.0d0 * sin( x(i) * pi / 2.0d0 )
				else
					g(i) = g(i) * cos( x(j) * pi / 2.0d0 )
				end if
			end do
			if ( ind > 1 ) g(i) = g(i) * sin( x(m-ind+1) * pi / 2.0d0 )
		end do
		
		if ( ind > 1 ) then
			g(m-ind+1) = faux * pi / 2.0d0 * cos( x(m-ind+1) * pi / 2.0d0 )
			do j = 1,m-ind
				g(m-ind+1) = g(m-ind+1) * cos( x(j) * pi / 2.0d0 )
			end do
		end if
				
		
		faux = 1.0d0
		do i = 1,m-ind
			faux = faux * cos( x(i) * pi / 2.0d0 )
		end do

		if ( ind > 1 ) faux = faux * sin( x(m-ind+1) * pi / 2.0d0 )
		
		do i = m,n
			g(i) = faux * 1.0d2 * ( 2.0d0 * ( x(i) - 0.5d0 ) + 2.0d1 * pi * sin( 2.0d1 * pi * ( x(i) - 0.5d0 ) ) )
		end do
					
		if ( penalize ) call penalizeg(n,x,g)
		return

	end if
	
! ----------------------------------------------------------------------

	!  DTLZ4

	if ( problem == 'DTLZ4' ) then 
	
!		k = dimk
!		m = dimm
	
!		faux = 0.0d0
!		do i = m,n
!			faux = faux + ( x(i) - 0.5d0 ) ** 2
!		end do
		
!		f = 1.0d0 + faux					
!		do i = 1,m-ind
!			f = f * cos( x(i)**alpha * pi / 2.0d0 )
!		end do
		
!		if ( ind > 1 ) f = f * sin( x(m-ind+1)**alpha * pi / 2.0d0 )
		
		k = dimk
		m = dimm
		
		do i = 1,n
			g(i) = 0.0d0
		end do
	
		faux = 0.0d0
		do i = m,n
			faux = faux + ( x(i) - 0.5d0 ) ** 2
		end do
		faux = 1.0d0 + faux
		
		do i = 1,m-ind
			g(i) = faux
		end do 
		
		do i = 1,m-ind
			do j = 1,m-ind
				if ( j == i ) then
					g(i) = - g(i) * pi / 2.0d0 * alpha * x(i)**( alpha - 1.0d0 ) * sin( x(i)**alpha * pi / 2.0d0 )
				else
					g(i) = g(i) * cos( x(j)**alpha * pi / 2.0d0 )
				end if
			end do
			if ( ind > 1 ) g(i) = g(i) * sin( x(m-ind+1)**alpha * pi / 2.0d0 )
		end do
		
		if ( ind > 1 ) then
			g(m-ind+1) = faux * pi / 2.0d0 * alpha * x(i)**( alpha - 1.0d0 ) * cos( x(m-ind+1)**alpha * pi / 2.0d0 )
			do j = 1,m-ind
				g(m-ind+1) = g(m-ind+1) * cos( x(j)**alpha * pi / 2.0d0 )
			end do
		end if
				
		
		faux = 1.0d0
		do i = 1,m-ind
			faux = faux * cos( x(i)**alpha * pi / 2.0d0 )
		end do

		if ( ind > 1 ) faux = faux * sin( x(m-ind+1)**alpha * pi / 2.0d0 )
		
		do i = m,n
			g(i) = faux * 2.0d0 * ( x(i) - 0.5d0 )
		end do
					
		if ( penalize ) call penalizeg(n,x,g)
		return

	end if
		
! ----------------------------------------------------------------------

	!  FA1

	if ( problem == 'FA1' ) then 
			
		if ( ind == 1 ) then
			g(1) = 4.0d0 * exp( -4.0d0 * x(1) )  / ( 1.0d0 - exp( -4.0d0 ) )
			g(2) = 0.0d0
			g(3) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			!faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) )
			!f = ( x(2) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(2) + 1.0d0  ) ) ** 0.5d0 )
						
			a = exp( -4.0d0 * x(1) )
			b = 1.0d0 - exp( -4.0d0 )
			t = ( 1.0d0 - a ) / ( b * ( x(2) + 1.0d0 ) )
			
			g(1) =  - 2.0d0 * a / b *  t ** (-0.5d0)
			g(2) = 1.0d0 - 0.5d0 * t ** 0.5d0 
			g(3) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
		
		if ( ind == 3 ) then
			!faux = ( 1.0d0 - exp( -4.0d0 * x(1) ) ) / ( 1.0d0 - exp( -4.0d0 ) )
			!f = ( x(3) + 1.0d0 ) * ( 1.0d0 - ( faux / ( x(3) + 1.0d0  ) ) ** 0.1d0 )

			a = exp( -4.0d0 * x(1) )
			b = 1.0d0 - exp( -4.0d0 )
			t = ( 1.0d0 - a ) / ( b * ( x(3) + 1.0d0  ) )
			
			g(1) = - 0.4d0 * t ** (-0.9d0) * a / b 
			g(2) = 0.0d0
			g(3) = 1.0d0 - 0.9d0 * t ** 0.1d0
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
	end if		
	
! ----------------------------------------------------------------------

	!  Far1

	if ( problem == 'Far1' ) then 
			
		if ( ind == 1 ) then
			g(1) = 6.0d1 * ( x(1) - 0.1d0 ) * exp( 1.5d1 * ( - ( x(1) - 0.1d0 ) ** 2 - x(2) ** 2 ) ) &
				 + 4.0d1 * ( x(1) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				 - 4.0d1 * ( x(1) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				 - 4.0d1 * ( x(1) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ** 2 - ( x(2) + 0.6d0 ) ** 2 ) ) &
				 - 4.0d1 * ( x(1) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ** 2 - ( x(2) + 0.6d0 ) ** 2 ) ) 
		
			g(2) = 6.0d1 * x(2) * exp( 1.5d1 * ( - ( x(1) - 0.1d0 ) ** 2 - x(2) ** 2 ) ) &
				 + 4.0d1 * ( x(2) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				 - 4.0d1 * ( x(2) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				 - 4.0d1 * ( x(2) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.6d0 ) ** 2 - ( x(2) + 0.6d0 ) ** 2 ) ) &
				 - 4.0d1 * ( x(2) + 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.6d0 ) ** 2 - ( x(2) + 0.6d0 ) ** 2 ) ) 		 
					 
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = - 8.0d1 * x(1) * exp( 2.0d1 * ( - x(1) ** 2 - x(2) ** 2 ) ) &
					 - 4.0d1 * ( x(1) - 0.4d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.4d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
					 + 4.0d1 * ( x(1) + 0.5d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.5d0 ) ** 2 - ( x(2) - 0.7d0 ) ** 2 ) ) &
					 + 4.0d1 * ( x(1) - 0.5d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.5d0 ) ** 2 - ( x(2) + 0.7d0 ) ** 2 ) ) &
					 - 4.0d1 * ( x(1) + 0.4d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.4d0 ) ** 2 - ( x(2) + 0.8d0 ) ** 2 ) ) 
						 
			g(2) = - 8.0d1 * x(2) * exp( 2.0d1 * ( - x(1) ** 2 - x(2) ** 2 ) ) &
				   - 4.0d1 * ( x(2) - 0.6d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.4d0 ) ** 2 - ( x(2) - 0.6d0 ) ** 2 ) ) &
				   + 4.0d1 * ( x(2) - 0.7d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.5d0 ) ** 2 - ( x(2) - 0.7d0 ) ** 2 ) ) &
				   + 4.0d1 * ( x(2) + 0.7d0 ) * exp( 2.0d1 * ( - ( x(1) - 0.5d0 ) ** 2 - ( x(2) + 0.7d0 ) ** 2 ) ) &
				   - 4.0d1 * ( x(2) + 0.8d0 ) * exp( 2.0d1 * ( - ( x(1) + 0.4d0 ) ** 2 - ( x(2) + 0.8d0 ) ** 2 ) ) 
						 
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	! FDS
	! NEWTON’S METHOD FOR MULTIOBJECTIVE OPTIMIZATION

	if ( problem == 'FDS' ) then 
			
		if ( ind == 1 ) then
			do i = 1,n
				g(i) = 4.0d0 * i * ( x(i) - i ) ** 3 / n ** 2
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			do i = 1,n
				g(i) = exp( sum(x)/n )/n + 2.0d0 * x(i)
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
		
		if ( ind == 3 ) then
			do i = 1,n
				g(i) = - i * ( n - i + 1.0d0 ) * exp( - x(i) ) / ( n * ( n + 1.0d0 ) ) 
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------

	! FF1 
	! C. M. Fonseca and P. J. Fleming, “An overview of evolutionary algorithms in multiobjective optimization

	if ( problem == 'FF1' ) then 
		
		if ( ind == 1 ) then
			g(1) = 2.0d0 * ( x(1) - 1.0d0 ) * exp( - ( x(1) - 1.0d0 ) ** 2 - ( x(2) + 1.0d0 ) ** 2 )
			g(2) = 2.0d0 * ( x(2) + 1.0d0 ) * exp( - ( x(1) - 1.0d0 ) ** 2 - ( x(2) + 1.0d0 ) ** 2 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * ( x(1) + 1.0d0 ) * exp( - ( x(1) + 1.0d0 ) ** 2 - ( x(2) - 1.0d0 ) ** 2 )
			g(2) = 2.0d0 * ( x(2) - 1.0d0 ) * exp( - ( x(1) + 1.0d0 ) ** 2 - ( x(2) - 1.0d0 ) ** 2 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if						
	
! ----------------------------------------------------------------------
	

	!  Hil1

	if ( problem == 'Hil1' ) then 
		a = 2.0d0 * pi / 3.6d2 * ( 4.5d1 + 4.0d1 * sin( 2.0d0 * pi * x(1) ) &
			+ 2.5d1 * sin( 2.0d0 * pi * x(2) ) )
		b = 1.0d0 + 0.5d0 * cos( 2.0d0 * pi * x(1) )
		
		if ( ind == 1 ) then
			g(1) = - 1.6d2 * pi ** 2 / 3.6d2 * cos( 2.0d0 * pi * x(1) ) * sin( a ) *  b &
			- pi * sin( 2.0d0 * pi * x(1) ) * cos( a ) 
			g(2) = - 1.0d2 * pi ** 2 / 3.6d2 * cos( 2.0d0 * pi * x(2) ) * sin( a ) *  b 
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 1.6d2 * pi ** 2 / 3.6d2 * cos( 2.0d0 * pi * x(1) ) * cos( a ) *  b &
			- pi * sin( 2.0d0 * pi * x(1) ) * sin( a ) 
			g(2) = 1.0d2 * pi ** 2 / 3.6d2 * cos( 2.0d0 * pi * x(2) ) * cos( a ) *  b 
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if
	
! ----------------------------------------------------------------------
	
	!  IKK1

	if ( problem == 'IKK1' ) then 
		
		if ( ind == 1 ) then
			g(1) = 2.0d0 * x(1)
			g(2) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * ( x(1) - 2.0d1 )
			g(2) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
		
		if ( ind == 3 ) then
			g(1) = 0.0d0
			g(2) = 2.0d0 * x(2)
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if				
	
! ----------------------------------------------------------------------
	
	!  IM1

	if ( problem == 'IM1' ) then 
		
		if ( ind == 1 ) then
			g(1) = 1.0d0 / sqrt( x(1) )
			g(2) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then			
			g(1) = ( 1.0d0 - x(2) )
			g(2) = - x(1)
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if
	
! ----------------------------------------------------------------------	
	
	! JOS1 
	! Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why Does It Work and How?
	
	if ( problem == 'JOS1' ) then
	
		if ( ind == 1 ) then					
			do i = 1,n
					g(i) = 2.0d0 * x(i) / n
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			do i = 1,n
					g(i) = 2.0d0 * ( x(i) - 2.0d0 ) / n
			end do	
			if ( penalize ) call penalizeg(n,x,g)
			return				
		end if
			
	end if
	
! ----------------------------------------------------------------------	
	
	! JOS4
	! Dynamic Weighted Aggregation for Evolutionary Multi-Objetive Optimization: Why Does It Work and How?
	
	if ( problem == 'JOS4' ) then 

		if ( ind == 1 ) then
			g(1) = 1.0d0
			g(2:n) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			t = x(1) / faux
			
			g(1) = - 0.25d0 * t ** (-0.75d0) - 4.0d0 * t ** 3
			do i = 2,n
				g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 - 0.75d0 * t ** 0.25d0 + 3.0d0 * t ** 4.0d0 )
			end do
			
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
	
	end if		

! ----------------------------------------------------------------------

	! 	KW2
	
	if ( problem == 'KW2' ) then 
			
		if ( ind == 1 ) then
			g(1) = 6.0d0 * ( 1.0d0 - x(1) ) * exp( -x(1)**2 - ( x(2) + 1.0d0 ) ** 2 ) &
				  + 6.0d0 * ( 1.0d0 - x(1) )**2 * exp( -x(1)**2 - ( x(2) + 1.0d0 ) ** 2 ) * x(1) &
				  + 1.0d1 * ( 1.0d0 / 5.0d0 - 3.0d0 * x(1)**2 ) * exp( - x(1)**2 - x(2)**2 ) &
				  - 2.0d1 * ( x(1) / 5.0d0 - x(1)**3 - x(2)**5 ) * exp( - x(1)**2 - x(2)**2 ) * x(1) &
				  - 6.0d0 * exp( -( x(1) + 2.0d0 )**2 - x(2)**2 ) * ( x(1) + 2.0d0 ) - 1.0d0
			g(2) = 6.0d0 * ( 1.0d0 - x(1) )**2 * exp( -x(1)**2 - ( x(2) + 1.0d0 ) ** 2 ) * ( x(2) + 1.0d0 ) &
						- 5.0d1 * x(2)**4 * exp( - x(1)**2 - x(2)**2 ) &
						- 1.0d1 * ( x(1) / 5.0d0 - x(1)**3 - x(2)**5 ) * exp( - x(1)**2 - x(2)**2 ) * 2.0d0 * x(2)&
						- 6.0d0 * exp( -( x(1) + 2.0d0 )**2 - x(2)**2 ) * x(2) - 0.5d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = - 6.0d0 * ( 1.0d0 + x(2) )**2 * exp( -x(2)**2 - ( 1.0d0 - x(1) ) ** 2 ) * ( 1.0d0 - x(1) ) &
						+ 5.0d1 * x(1)**4 * exp( - x(1)**2 - x(2)**2 ) &
						- 2.0d1 * ( - x(2) / 5.0d0 + x(2)**3 + x(1)**5 ) * exp( - x(1)**2 - x(2)**2 ) * x(1) &
						- 6.0d0 * exp( -( 2.0d0 - x(2) )**2 - x(1)**2 ) * x(1)
			
			g(2) = - 6.0d0 * ( 1.0d0 + x(2) ) * exp( -x(2)**2 - ( 1.0d0 - x(1) ) ** 2 ) &
						 + 6.0d0 * ( 1.0d0 + x(2) )**2 * exp( -x(2)**2 - ( 1.0d0 - x(1) ) ** 2 ) * x(2) &
						 + 1.0d1 * ( - 1.0d0 / 5.0d0 + 3.0d0 * x(2)**2 ) * exp( - x(1)**2 - x(2)**2 ) &
						 - 2.0d1 * ( - x(2) / 5.0d0 + x(2)**3 + x(1)**5 ) * exp( - x(1)**2 - x(2)**2 ) * x(2)&
						 + 6.0d0 * exp( -( 2.0d0 - x(2) )**2 - x(1)**2 ) * ( 2.0d0 - x(2) )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if

	
! ----------------------------------------------------------------------

	! 	LE1
	
	if ( problem == 'LE1' ) then 
											
		if ( ind == 1 ) then			
			t = 0.25d0 * ( x(1) ** 2 + x(2) ** 2 ) ** (-0.875d0)
			
			g(1) = x(1) * t
			g(2) = x(2) * t
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			t = 0.5d0 * ( ( x(1) - 0.5d0 ) ** 2 + ( x(2) - 0.5d0 ) ** 2 ) ** (-0.75d0)
			
			g(1) = ( x(1) - 0.5d0 ) * t
			g(2) = ( x(2) - 0.5d0 ) * t
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
			
	end if
	
! ----------------------------------------------------------------------
	
	! Lov1  
	
	if ( problem == 'Lov1' ) then 
			
		if ( ind == 1 ) then
			g(1) = 2.1d0 * x(1)
			g(2) = 2.0d0 * 0.98d0 * x(2)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * 0.99d0 * ( x(1) - 3.0d0 )
			g(2) = 2.0d0 * 1.03d0 * ( x(2) - 2.5d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------
	
	! Lov2
	
	if ( problem == 'Lov2' ) then 

		if ( ind == 1 ) then
			g(1) = 0.0d0
			g(2) = 1.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) =  ( - 3.0d0 * x(1) ** 2 * ( x(1) + 1.0d0 ) - ( x(2) - x(1) ** 3 ) )
			g(1) = - g(1) / ( x(1) + 1.0d0 ) ** 2
			
			g(2) = - 1.0d0 / ( x(1) + 1.0d0 ) 
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------
	
	! Lov3 
	
	if ( problem == 'Lov3' ) then 
		
		if ( ind == 1 ) then
			g(1) =  2.0d0 * x(1)
			g(2) =  2.0d0 * x(2)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) =  2.0d0 * ( x(1) - 6.0d0 )
			g(2) = - 2.0d0 * ( x(2) + 0.3d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------
	
	! Lov4 
	
	if ( problem == 'Lov4' ) then 
			
		if ( ind == 1 ) then
			g(1) =  2.0d0 * x(1) - 8.0d0 * ( ( x(1) + 2.0d0 ) * exp( - ( x(1) + 2.0d0 ) ** 2 - x(2) ** 2 ) &
			+ ( x(1) - 2.0d0 ) * exp( - ( x(1) - 2.0d0 ) ** 2 - x(2) ** 2 ) )
			g(2) =  2.0d0 * x(2) - 8.0d0 * ( x(2) * exp( - ( x(1) + 2.0d0 ) ** 2 - x(2) ** 2 ) &
			+ x(2) * exp( - ( x(1) - 2.0d0 ) ** 2 - x(2) ** 2 ) )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) =  2.0d0 * ( x(1) - 6.0d0 )
			g(2) =  2.0d0 * ( x(2) + 0.5d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if	
	
! ----------------------------------------------------------------------

	! Lov5
	
	if ( problem == 'Lov5' ) then 
	
		MM(1,:) = (/ -1.0d0  , -0.03d0,  0.011d0/)
		MM(2,:) = (/ -0.03d0 , -1.0d0 ,  0.07d0 /)
		MM(3,:) = (/  0.011d0,  0.07d0, -1.01d0 /)
		
		p = (/ x(1), x(2) - 0.15d0,  x(3)/)
		a = 0.35d0
		
		A1 = sqrt( 2.0d0 * pi / a ) * exp( dot_product( p, matmul( MM, p ) ) / a ** 2 )
		
		p(:) = (/ x(1) , x(2) + 1.1d0,  0.5d0 * x(3) /)
		a = 3.0d0
		
		A2 = sqrt( 2.0d0 * pi / a ) * exp( dot_product( p, matmul( MM, p ) ) / a ** 2 )
	
		if ( ind == 1 ) then
			g(1) = sqrt(2.0d0)/2.0d0 + sqrt(2.0d0)/2.0d0 * A1 * &
			( 2.0d0 * MM(1,1)* x(1) + 2.0d0 * MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) - 0.15d0 ) ) / 0.35d0**2 &
			+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,1) * x(1) + MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) + 1.1d0 ) ) / 3.0d0**2
			g(1) = - g(1)
			
			g(2) = sqrt(2.0d0)/2.0d0 * A1 * &
			( 2.0d0 * MM(1,2)* x(1) + 2.0d0 * MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) - 0.15d0 ) ) / 0.35d0**2 &
			+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,2) * x(1) + MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) + 1.1d0 ) ) / 3.0d0**2
			g(2) = - g(2)
			
			g(3) = sqrt(2.0d0)/2.0d0 * A1 * &
			( 2.0d0 * MM(1,3)* x(1) + 2.0d0 * MM(3,3) * x(3) + 2.0d0 * MM(2,3) * ( x(2) - 0.15d0 ) ) / 0.35d0**2 &
			+ sqrt(2.0d0)/2.0d0 * A2 * ( MM(1,3) * x(1) + MM(3,3) * x(3) / 2.0d0 + MM(2,3) * ( x(2) + 1.1d0 ) ) / 3.0d0**2
			g(3) = - g(3)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = - sqrt(2.0d0)/2.0d0 + sqrt(2.0d0)/2.0d0 * A1 * &
			( 2.0d0 * MM(1,1)* x(1) + 2.0d0 * MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) - 0.15d0 ) ) / 0.35d0**2 &
			+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,1) * x(1) + MM(1,3) * x(3) + 2.0d0 * MM(1,2) * ( x(2) + 1.1d0 ) ) / 3.0d0**2
			g(1) = - g(1)
			
			g(2) = sqrt(2.0d0)/2.0d0 * A1 * &
			( 2.0d0 * MM(1,2)* x(1) + 2.0d0 * MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) - 0.15d0 ) ) / 0.35d0**2 &
			+ sqrt(2.0d0)/2.0d0 * A2 * ( 2.0d0 * MM(1,2) * x(1) + MM(2,3) * x(3) + 2.0d0 * MM(2,2) * ( x(2) + 1.1d0 ) ) / 3.0d0**2
			g(2) = - g(2)
			
			g(3) = sqrt(2.0d0)/2.0d0 * A1 * &
			( 2.0d0 * MM(1,3)* x(1) + 2.0d0 * MM(3,3) * x(3) + 2.0d0 * MM(2,3) * ( x(2) - 0.15d0 ) ) / 0.35d0**2 &
			+ sqrt(2.0d0)/2.0d0 * A2 * ( MM(1,3) * x(1) + MM(3,3) * x(3) / 2.0d0 + MM(2,3) * ( x(2) + 1.1d0 ) ) / 3.0d0**2
			g(3) = - g(3)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------

	! Lov6
	
	if ( problem == 'Lov6' ) then 
	
		if ( ind == 1 ) then
			g(1)   = 1.0d0
			g(2:6) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = - 0.5d0 / sqrt( x(1) ) - sin( 1.0d1 * pi * x(1) ) - 1.0d1 * pi * x(1) * cos( 1.0d1 * pi * x(1) )
			do i = 2,6
				g(i) = 2.0d0 * x(i)
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
			
	end if		
	
! ----------------------------------------------------------------------	

	! 	LTDZ
	!	Combining convergence and diversity in evolutionary multiobjective optimization
	
	if ( problem == 'LTDZ' ) then 
											
		if ( ind == 1 ) then
			g(1) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * sin( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 )
			g(2) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 )
			g(3) = - cos( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 )
			g = - g
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * sin( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 )
			g(2) = - pi / 2.0d0 * ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * cos( x(2) * pi / 2.0d0 )
			g(3) = - cos( x(1) * pi / 2.0d0 ) * sin( x(2) * pi / 2.0d0 )
			g = - g
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 3 ) then
			g(1) = pi / 2.0d0 * ( 1.0d0 + x(3) ) * sin( x(1) * pi / 2.0d0 ) * sin( x(1) * pi / 2.0d0 ) &
				 - pi / 2.0d0 * ( 1.0d0 + x(3) ) * cos( x(1) * pi / 2.0d0 ) * cos( x(1) * pi / 2.0d0 )
			g(2) = 0.0d0
			g(3) = - cos( x(1) * pi / 2.0d0 ) * sin( x(1) * pi / 2.0d0 )
			g = - g
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
			
	end if		

! ----------------------------------------------------------------------

	! 	MGH9 
	
	if ( problem == 'MGH9' ) then 
				
		t = ( 8.0d0 - ind ) / 2.0d0
		
		g(1) = exp( - x(2) * ( t - x(3) ) ** 2 / 2.0d0 )
		g(2) = - x(1) * exp( - x(2) * ( t - x(3) ) ** 2 / 2.0d0 ) * ( t - x(3) ) ** 2 / 2.0d0
		g(3) = x(1) * exp( - x(2) * ( t - x(3) ) ** 2 / 2.0d0 ) * x(2) * ( t - x(3) ) 
						
		if ( penalize ) call penalizeg(n,x,g)
			return
			
	end if		
	
! ----------------------------------------------------------------------

	! 	MGH16 
	
	if ( problem == 'MGH16' ) then 
			
		t = ind / 5.0d0
		
		g(1) = 2.0d0 * ( x(1) + t * x(2) - exp(t) )
		g(2) = 2.0d0 * t * ( x(1) + t * x(2) - exp(t) )
		g(3) = 2.0d0 * ( x(3) + x(4) * sin(t) - cos(t) )
		g(4) = 2.0d0 * sin(t) * ( x(3) + x(4) * sin(t) - cos(t) )
							
		if ( penalize ) call penalizeg(n,x,g)
			return
			
	end if		
	
! ----------------------------------------------------------------------

	! 	MGH26 
	
	if ( problem == 'MGH26' ) then 
								
		t = 0.0d0
		do i = 1,n
			t = t + cos(x(i))
		end do
		
		gaux1 = 2.0d0 * ( n - t + ind * ( 1.0d0 - cos(x(ind)) ) - sin(x(ind)) )
		
		
		do i = 1,n
			g(i) = gaux1 * sin(x(i));
		end do
		
!		g(1) = gaux1 * sin(x(1))
!		g(2) = gaux1 * sin(x(2))
!		g(3) = gaux1 * sin(x(3))
!		g(4) = gaux1 * sin(x(4))
		
		g(ind) = g(ind)  + gaux1 * ( ind * sin(x(ind)) - cos(x(ind)) )
				
		if ( penalize ) call penalizeg(n,x,g)
			return
		
	end if				
	
! ----------------------------------------------------------------------

	! 	MGH33
	
	if ( problem == 'MGH33' ) then 

		faux = 0.0d0					
		do i = 1,n
			faux = faux + i * x(i)
		end do
		
		faux = 2.0d0 * ( ind * faux - 1.0d0 ) 
	
		do i = 1,n
			g(i) = faux * real(i) * real(ind)
		end do
						
		if ( penalize ) call penalizeg(n,x,g)
		return
			
	end if				
	
! ----------------------------------------------------------------------

	! 	MHHM2
	
	if ( problem == 'MHHM2' ) then 
											
		if ( ind == 1 ) then
			g(1) = 2.0d0 * ( x(1) - 0.8d0 )
			g(2) = 2.0d0 * ( x(2) - 0.6d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then			
			g(1) = 2.0d0 * ( x(1) - 0.85d0 )
			g(2) = 2.0d0 * ( x(2) - 0.7d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
			
		if ( ind == 3 ) then			
			g(1) = 2.0d0 * ( x(1) - 0.9d0 )
			g(2) = 2.0d0 * ( x(2) - 0.6d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
				
	end if	
	
! ----------------------------------------------------------------------

	!  MLF1

	if ( problem == 'MLF1' ) then 
			
		if ( ind == 1 ) then
			g(1) = sin( x(1) ) / 2.0d1 + ( 1.0d0 + x(1) / 2.0d1 ) * cos( x(1) )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = cos( x(1) ) / 2.0d1 - ( 1.0d0 + x(1) / 2.0d1 ) * sin( x(1) )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  MLF2

	if ( problem == 'MLF2' ) then 
			
		if ( ind == 1 ) then
			g(1) = ( 2.0d0 * x(1) * ( x(1) ** 2 + x(2) - 1.1d1 ) + ( x(1) + x(2) ** 2 - 7.0d0 )  ) / 1.0d2
			g(2) = ( ( x(1) ** 2 + x(2) - 1.1d1 ) + 2.0d0 * x(2) * ( x(1) + x(2) ** 2 - 7.0d0 )  ) / 1.0d2
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = ( 8.0d0 * x(1) * ( 4.0d0 * x(1) ** 2 + 2.0d0 * x(2) - 1.1d1 ) &
				 + 2.0d0 * ( 2.0d0 * x(1) + 4.0d0 *  x(2) ** 2 - 7.0d0 ) ) / 1.0d2
			g(2) = ( 2.0d0 * ( 4.0d0 * x(1) ** 2 + 2.0d0 * x(2) - 1.1d1 ) &
				 + 8.0d0 * x(2) * ( 2.0d0 * x(1) + 4.0d0 *  x(2) ** 2 - 7.0d0 ) ) / 1.0d2
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	!  MMR1			
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	
	
!	if ( problem == 'MMR1' ) then 
			
!			if ( ind == 1 ) then
!				g(1) = 2.0d0 * x(1) 
!				g(2) = 0.0d0	
!				if ( penalize ) call penalizeg(n,x,g)
!			return
!			end if
			
!			if ( ind == 2 ) then
!				g(1) = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ** 2 ) &
!					   - exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ** 2 )
!				g(1) = - 2.0d0 * x(1) * g(1) / ( 1.0d0 + x(1) ** 2 ) ** 2
				
!				g(2) = 1.0d1 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ** 2 ) * ( x(2) - 0.6d0 ) &
!				   + 1.25d3 * exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ** 2 ) * ( x(2) - 0.2d0 ) 
!				g(2) = g(2) / ( 1.0d0 + x(1) ** 2 )   
!				if ( penalize ) call penalizeg(n,x,g)
!			return	
!			end if	
			
!	end if	

	if ( problem == 'MMR1' ) then 
			
			if ( ind == 1 ) then
				g(1) = 1.0d0 
				g(2) = 0.0d0	
				if ( penalize ) call penalizeg(n,x,g)
			return
			end if
			
			if ( ind == 2 ) then
				g(1) = 2.0d0 - 0.8d0 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ** 2 ) &
					   - exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ** 2 )
				g(1) = - g(1) / x(1) ** 2 
				
				g(2) = 1.0d1 * exp( - ( ( x(2) - 0.6d0 ) / 0.4d0 ) ** 2 ) * ( x(2) - 0.6d0 ) &
				   + 1.25d3 * exp( - ( ( x(2) - 0.2d0 ) / 0.04d0 ) ** 2 ) * ( x(2) - 0.2d0 ) 
				g(2) = g(2) / x(1)   
				if ( penalize ) call penalizeg(n,x,g)
			return	
			end if	
			
	end if		
	
! ----------------------------------------------------------------------

!  MMR2
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	

	if ( problem == 'MMR2' ) then 
			
		if ( ind == 1 ) then
			g(1) = 1.0d0 
			g(2) = 0.0d0	
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
		
!			faux = x(1) / ( 1.0d0 + 1.0d1 * x(2) )
!			f = 1.0d0 - faux ** 2 - faux * sin( 8.0d0 * pi * x(1) )
!			f = f * ( 1.0d0 + 1.0d1 * x(2) )

			a = 1.0d0 + 1.0d1 * x(2)
			faux = x(1) / a
			
			g(1) = - 2.0d0 * faux / a - sin( 8.0d0 * pi * x(1) ) / a - 8.0d0 * pi * faux * cos( 8.0d0 * pi * x(1) )
			g(1) = g(1) * a
			
			g(2) = 1.0d1 * ( 1.0d0 - faux ** 2 - faux * sin( 8.0d0 * pi * x(1) ) )
			g(2) = g(2) + a * ( 2.0d1 * faux * x(1) / a**2 + 1.0d1 / a**2 * x(1) * sin( 8.0d0 * pi * x(1) ) ) 
		
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if	
			
	end if		
		
! ----------------------------------------------------------------------

	!  MMR3
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	

	if ( problem == 'MMR3' ) then 
			
		if ( ind == 1 ) then
			g(1) = 3.0d0 * x(1) ** 2
			g(2) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = - 3.0d0 * ( x(2) - x(1) ) ** 2
			g(2) = 3.0d0 * ( x(2) - x(1) ) ** 2
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if	
			
	end if			
	
! ----------------------------------------------------------------------

	!  MMR4
	!  Box-constrained multi-objective optimization: A gradient-like method without ‘‘a priori’’ scalarization	

	if ( problem == 'MMR4' ) then 
			
		if ( ind == 1 ) then	
			t = ( 2.0d0 * x(1) + x(2) + 2.0d0 * x(3) + 1.0d0 ) ** 2
			g(1) = 1.0d0 + 7.2d1 / t
			g(2) = - 2.0d0 + 3.6d1 / t
			g(3) = - 1.0d0 + 7.2d1 / t
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then			
			g(1) = - 3.0d0
			g(2) =   1.0d0
			g(3) = - 1.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if	
			
	end if					
		
! ----------------------------------------------------------------------

	!  MOP 2	

	if ( problem == 'MOP2' ) then 
			
		if ( ind == 1 ) then
			faux = 0.0d0
			do i = 1,n
				faux = faux + ( x(i) - 1.0d0 / ( sqrt(real(n)) ) ) ** 2
			end do
			
			do i = 1,n
				g(i) = 2.0d0 * ( x(i) - 1.0d0 / sqrt(real(n)) ) * exp( - faux )
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 0.0d0
			do i = 1,n
				faux = faux + ( x(i) + 1.0d0 / ( sqrt(real(n)) ) ) ** 2
			end do
			
			do i = 1,n
				g(i) = 2.0d0 * ( x(i) + 1.0d0 / sqrt(real(n)) ) * exp( - faux )
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if	
			
	end if			

! ----------------------------------------------------------------------

	!  MOP 3	

	if ( problem == 'MOP3' ) then 
			
		if ( ind == 1 ) then
			A1 = 0.5d0 * sin(1.0d0) - 2.0d0 * cos(1.0d0) + sin(2.0d0) - 1.5d0 * cos(2.0d0) 
			A2 = 1.5d0 * sin(1.0d0) - cos(1.0d0) + 2.0d0 * sin(2.0d0) - 0.5d0 * cos(2.0d0)
			B1 = 0.5d0 * sin(x(1)) - 2.0d0 * cos(x(1)) + sin(x(2)) - 1.5d0 * cos(x(2)) 
			B2 = 1.5d0 * sin(x(1)) - cos(x(1)) + 2.0d0 * sin(x(2)) - 0.5d0 * cos(x(2))	
			g(1) = 2.0d0 * ( A1 - B1 ) *  ( - 0.5d0 * cos(x(1)) - 2.0d0 * sin(x(1)) ) +&
				   2.0d0 * ( A2 - B2 ) *  ( - 1.5d0 * cos(x(1)) - sin(x(1)) )
			
			g(2) = 2.0d0 * ( A1 - B1 ) *  ( - cos(x(2)) - 1.5d0 * sin(x(2)) ) +&
				   2.0d0 * ( A2 - B2 ) *  ( - 2.0d0 * cos(x(2)) - 0.5d0 * sin(x(2)) )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) =  2.0d0 * ( x(1) + 3.0d0 )
			g(2) =  2.0d0 * ( x(2) + 1.0d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if	
			
	end if
	
! ----------------------------------------------------------------------

	!  MOP 5	

	if ( problem == 'MOP5' ) then 
			
		if ( ind == 1 ) then
			g(1) = x(1) + 2.0d0 * x(1) * cos( x(1) ** 2 + x(2) ** 2 )
			g(2) = x(2) + 2.0d0 * x(2) * cos( x(1) ** 2 + x(2) ** 2 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 3.0d0 * ( 3.0d0 * x(1) - 2.0d0 * x(2) + 4.0d0 ) / 4.0d0 &
				 + 2.0d0 * ( x(1) - x(2) + 1.0d0 ) / 2.7d1
			g(2) = - 2.0d0 * ( 3.0d0 * x(1) - 2.0d0 * x(2) + 4.0d0 ) / 4.0d0 &
				 - 2.0d0 * ( x(1) - x(2) + 1.0d0 ) / 2.7d1	 	 
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
		
		if ( ind == 3 ) then
			g(1) = - 2.0d0 * x(1) / ( x(1) ** 2 + x(2) ** 2 + 1.0d0 ) ** 2 &
				   + 2.2d0 * x(1) * exp( - x(1) ** 2 - x(2) ** 2 )
			g(2) = - 2.0d0 * x(2) / ( x(1) ** 2 + x(2) ** 2 + 1.0d0 ) ** 2 &
				   + 2.2d0 * x(2) * exp( - x(1) ** 2 - x(2) ** 2 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if	
				
	end if		
	
! ----------------------------------------------------------------------

	!  MOP6

	if ( problem == 'MOP6' ) then 
			
		if ( ind == 1 ) then
			g(1) = 1.0d0
			g(2) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			a = 1.0d0 + 1.0d1 * x(2) 
			b = sin( 8.0d0 * pi * x(1) )
			t = x(1) / a
			
			g(1) = - 2.0d0 * t - b - 8.0d0 * pi * x(1) * cos( 8.0d0 * pi * x(1) )
			g(2) = 1.0d1 * ( 1.0d0 - t ** 2 - t * b ) + 1.0d1 * x(1) / a * ( 2.0d0 * t + b )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		

! ----------------------------------------------------------------------

	!  MOP 7

	if ( problem == 'MOP7' ) then 
			
		if ( ind == 1 ) then
			g(1) =  x(1) - 2.0d0
			g(2) =  2.0d0 * ( x(2) + 1.0d0 ) / 1.3d1
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) =  ( x(1) + x(2) - 3.0d0 ) / 1.8d1 - ( - x(1) + x(2) + 2.0d0 ) / 4.0d0
			g(2) =  ( x(1) + x(2) - 3.0d0 ) / 1.8d1 + ( - x(1) + x(2) + 2.0d0 ) / 4.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
		
		if ( ind == 3 ) then
			g(1) =  2.0d0 * ( x(1) + 2.0d0 * x(2) - 1.0d0 ) / 1.75d2 &
			- 2.0d0 * ( - x(1) + 2.0d0 * x(2) ) / 1.7d1
			g(2) =  4.0d0 * ( x(1) + 2.0d0 * x(2) - 1.0d0 ) / 1.75d2 &
			+ 4.0d0 * ( - x(1) + 2.0d0 * x(2) ) / 1.7d1
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if			
			
	end if		
	
! ----------------------------------------------------------------------

	! 	PNR
	
	if ( problem == 'PNR' ) then 
										
		if ( ind == 1 ) then
			g(1) = 4.0d0 * x(1) ** 3 - 2.0d0 * x(1) - 1.0d1 * x(2)
			g(2) = 4.0d0 * x(2) ** 3 + 2.0d0 * x(2) - 1.0d1 * x(1)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * x(1)
			g(2) = 2.0d0 * x(2)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  vu

	if ( problem == 'QV1' ) then 
			
		if ( ind == 1 ) then
			faux = 0.0d0
			do i = 1,n
				faux = faux + x(i) ** 2 - 1.0d1 * cos( 2.0d0 * pi * x(i) ) + 1.0d1
			end do
			faux = 0.25 * ( faux / n ) ** (-0.75d0)
			
			do i = 1,n
				g(i) = faux * ( 2.0d0 * x(i) + 2.0d1 * pi * sin( 2.0d0 * pi * x(i) )  ) / n
			end do
			
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 0.0d0
			do i = 1,n
				faux = faux + ( x(i) - 1.5d0 ) ** 2 - 1.0d1 * cos( 2.0d0 * pi * ( x(i) -1.5d0 ) ) + 1.0d1
			end do
		
			faux = 0.25 * ( faux / n ) ** (-0.75d0)
			
			do i = 1,n
				g(i) = faux * ( 2.0d0 * ( x(i) - 1.5d0 ) + 2.0d1 * pi * sin( 2.0d0 * pi * ( x(i) - 1.5d0 ) )  ) / n
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if

! ----------------------------------------------------------------------

	! SD

	if ( problem == 'SD' ) then 
			
		if ( ind == 1 ) then
			g(1) = 2.0d0
			g(2) = sqrt(2.0d0)
			g(3) = sqrt(2.0d0)
			g(4) = 1.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = - 2.0d0 / x(1) ** 2
			g(2) = - 2.0d0 * sqrt(2.0d0) / x(2) ** 2
			g(3) = - 2.0d0 * sqrt(2.0d0) / x(3) ** 2
			g(4) = - 2.0d0 / x(4) ** 2
			if ( penalize ) call penalizeg(n,x,g)
			return		
		end if	
			
	end if		

! ----------------------------------------------------------------------

	! SLCDT1

	if ( problem == 'SLCDT1' ) then 
	
		if ( ind == 1 ) then			
			g(1) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ** 2 ) + &
						( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ** 2 ) + 1.0d0 )  &
						 - 2.0d0 * 0.85d0 *  ( x(1) + x(2) ) * exp( - ( x(1) + x(2) ) ** 2 )
			g(2) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ** 2 ) - &
						( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ** 2 ) - 1.0d0 )  &
						 - 2.0d0 * 0.85d0 * ( x(1) + x(2) ) * exp( - ( x(1) + x(2) ) ** 2 )				 
			if ( penalize ) call penalizeg(n,x,g)
			return			 
		end if
		
		if ( ind == 2 ) then
			g(1) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ** 2 ) + &
						( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ** 2 ) - 1.0d0 )  &
						 - 2.0d0 * 0.85d0 * exp( - ( x(1) + x(2) ) ** 2 ) * ( x(1) + x(2) )
			g(2) = 0.5d0 * ( ( x(1) + x(2) )/sqrt( 1.0d0 + ( x(1) + x(2) ) ** 2 ) - &
						( x(1) - x(2) )/sqrt( 1.0d0 + ( x(1) - x(2) ) ** 2 ) + 1.0d0 )  &
						 - 2.0d0 * 0.85d0 * ( x(1) + x(2) ) * exp( - ( x(1) + x(2) ) ** 2 )
			if ( penalize ) call penalizeg(n,x,g)
			return			 
		end if 	
			
	end if	
	
! ----------------------------------------------------------------------

	! 	SLCDT2
	!	Convergence of stochastic search algorithms to finite size pareto set approximations 
	
	if ( problem == 'SLCDT2' ) then 
											
		if ( ind == 1 ) then
			g(1) = 4.0d0 * ( x(1) - 1.0d0 ) ** 3
			do i = 2,n
				g(i) = 2.0d0 * ( x(i) - 1.0d0 )
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(2) = 4.0d0 * ( x(2) + 1.0d0 ) ** 3
			
			do i = 1,n
				if ( i /= 2 ) g(i) = 2.0d0 * ( x(i) + 1.0d0 )
			end do				
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
		
		if ( ind == 3 ) then
			g(3) = 4.0d0 * ( x(3) - 1.0d0 ) ** 3
			do i = 1,n
				if ( i /= 3 ) g(i) = 2.0d0 * ( x(i) - ( - 1.0d0 ) ** (i+1) ) 
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
			
	end if			
	
! ----------------------------------------------------------------------

	!  SP1

	if ( problem == 'SP1' ) then 
		
		if ( ind == 1 ) then
			g(1) = 2.0d0 * ( x(1) - 1.0d0 ) + 2.0d0 * ( x(1) - x(2) )
			g(2) = - 2.0d0 * ( x(1) - x(2) )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * ( x(1) - x(2) )
			g(2) = 2.0d0 * ( x(2) - 3.0d0 ) - 2.0d0 * ( x(1) - x(2) )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if			
	
! ----------------------------------------------------------------------

	!  SSFYY2

	if ( problem == 'SSFYY2' ) then 
			
		if ( ind == 1 ) then
			g(1) = 2.0d0 * x(1) + 5.0d0 * pi * sin( x(1) * pi / 2.0d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * ( x(1) - 4.0d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	!  SK1

	if ( problem == 'SK1' ) then 
			
		if ( ind == 1 ) then
			g(1) = 4.0d0 * x(1) ** 3 + 9.0d0 * x(1) ** 2 - 2.0d1 * x(1) - 1.0d1
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * x(1) ** 3 - 6.0d0 * x(1) ** 2 - 2.0d1 * x(1) + 1.0d1
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	!  SK2

	if ( problem == 'SK2' ) then 
			
		if ( ind == 1 ) then
			g(1) = 2.0d0 * ( x(1) - 2.0d0 )
			g(2) = 2.0d0 * ( x(2) + 3.0d0 )
			g(3) = 2.0d0 * ( x(3) - 5.0d0 )
			g(4) = 2.0d0 * ( x(4) - 4.0d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + ( x(1) ** 2 + x(2) ** 2 + x(3) ** 2 + x(4) ** 2 ) / 1.0d2
			
			g(1) = ( - cos( x(1) ) * faux + ( sin( x(1) ) + sin( x(2) ) &
			+ sin( x(3) ) + sin( x(4) ) ) * x(1) / 5.0d1 ) / faux ** 2
			g(2) = ( - cos( x(2) ) * faux + ( sin( x(1) ) + sin( x(2) ) &
			+ sin( x(3) ) + sin( x(4) ) ) * x(2) / 5.0d1 ) / faux ** 2
			g(3) = ( - cos( x(3) ) * faux + ( sin( x(1) ) + sin( x(2) ) &
			+ sin( x(3) ) + sin( x(4) ) ) * x(3) / 5.0d1 ) / faux ** 2
			g(4) = ( - cos( x(4) ) * faux + ( sin( x(1) ) + sin( x(2) ) &
			+ sin( x(3) ) + sin( x(4) ) ) * x(4) / 5.0d1 ) / faux ** 2
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if					
	
! ----------------------------------------------------------------------

	!  TKLY1

	if ( problem == 'TKLY1' ) then 
			
		if ( ind == 1 ) then
			g(1)   = 1.0d0
			g(2:4) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			A1 = ( 2.0d0 - exp( - ( ( x(2) - 0.1d0 ) / 4.0d-3 ) ** 2 ) &
				         - 0.8d0 * exp( - ( ( x(2) - 0.9d0 ) / 4.0d-1 ) ** 2 ) )
			A2 = ( 2.0d0 - exp( - ( ( x(3) - 0.1d0 ) / 4.0d-3 ) ** 2 ) &
				         - 0.8d0 * exp( - ( ( x(3) - 0.9d0 ) / 4.0d-1 ) ** 2 ) )
			A3 = ( 2.0d0 - exp( - ( ( x(4) - 0.1d0 ) / 4.0d-3 ) ** 2 ) &
				         - 0.8d0 * exp( - ( ( x(4) - 0.9d0 ) / 4.0d-1 ) ** 2 ) )
				         
				         
			
			g(1) = - A1 * A2 * A3 / x(1) ** 2
			g(2) = A2 * A3 / x(1) * ( 5.0d2 * exp( - ( ( x(2) - 0.1d0 ) / 4.0d-3 ) ** 2 ) * ( ( x(2) - 0.1d0 ) / 4.0d-3 )  &
			+ 4.0d0 * exp( - ( ( x(2) - 0.9d0 ) / 4.0d-1 ) ** 2 ) * ( ( x(2) - 0.9d0 ) / 4.0d-1 ) )
			g(3) = A1 * A3 / x(1) * ( 5.0d2 * exp( - ( ( x(3) - 0.1d0 ) / 4.0d-3 ) ** 2 ) * ( ( x(3) - 0.1d0 ) / 4.0d-3 )  &
			+ 4.0d0 * exp( - ( ( x(3) - 0.9d0 ) / 4.0d-1 ) ** 2 ) * ( ( x(3) - 0.9d0 ) / 4.0d-1 ) )
			g(4) = A1 * A2 / x(1) * ( 5.0d2 * exp( - ( ( x(4) - 0.1d0 ) / 4.0d-3 ) ** 2 ) * ( ( x(4) - 0.1d0 ) / 4.0d-3 )  &
			+ 4.0d0 * exp( - ( ( x(4) - 0.9d0 ) / 4.0d-1 ) ** 2 ) * ( ( x(4) - 0.9d0 ) / 4.0d-1 ) )
			
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	! 	Toi4
	
	if ( problem == 'Toi4' ) then 
			
		if ( ind == 1 ) then
			g(1) = 2.0d0 * x(1)
			g(2) = 2.0d0 * x(2)
			g(3) = 0.0d0
			g(4) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = x(1) - x(2)
			g(2) = - ( x(1) - x(2) )
			g(3) = x(3) - x(4)
			g(4) = - ( x(3) - x(4) )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	! 	Toi8
	
	if ( problem == 'Toi8' ) then 
	
		g(:) = 0.0d0
		
		if ( ind == 1 ) then
			g(1) = 4.0d0 * ( 2.0d0 * x(1) - 1.0d0 )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind /= 1 ) then
			g(ind-1) = 4.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) )
			g(ind)   = - 2.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) )
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if	
	
! ----------------------------------------------------------------------

	! 	Toi9
	
	if ( problem == 'Toi9' ) then 
	
		g(:) = 0.0d0
		
		if ( ind == 1 ) then
			g(1) = 4.0d0 * ( 2.0d0 * x(1) - 1.0d0 )
			g(2) = 2.0d0 * x(2)
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind > 1 .and. ind < n ) then
			g(ind-1) =  4.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) ) &
					  - 2.0d0 * ( ind - 1.0d0 ) * x(ind-1)
			g(ind)   = - 2.0d0 * ind * ( 2.0d0 * x(ind-1) - x(ind) ) + 2.0d0 * ind * x(ind)
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
		
		if ( ind == n ) then
			g(n-1) = 4.0d0 * n * ( 2.0d0 * x(n-1) - x(n) ) - 2.0d0 * ( n - 1.0d0 ) * x(n-1)
			g(n)   = - 2.0d0 * n * ( 2.0d0 * x(n-1) - x(n) )
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
			
	end if
	
! ----------------------------------------------------------------------

	! 	Toi10 (Rosenbrock)
	
	if ( problem == 'Toi10' ) then 
			
		g(:) = 0.0d0
		
		g(ind)   = - 4.0d2 * ( x(ind+1) - x(ind) ** 2 ) * x(ind)
		g(ind+1) = 2.0d2 * ( x(ind+1) - x(ind) ** 2 ) + 2.0d0 * ( x(ind+1) - 1.0d0 )
		if ( penalize ) call penalizeg(n,x,g)
			return
						
	end if			
		
! ----------------------------------------------------------------------

	!  VU1

	if ( problem == 'VU1' ) then 
			
		if ( ind == 1 ) then
			g(1) = - 2.0d0 * x(1) / ( x(1) ** 2 + x(2) ** 2 + 1.0d0 ) ** 2
			g(2) = - 2.0d0 * x(2) / ( x(1) ** 2 + x(2) ** 2 + 1.0d0 ) ** 2
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * x(1)
			g(2) = 6.0d0 * x(2)
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  VU2

	if ( problem == 'VU2' ) then 
			
		if ( ind == 1 ) then
			g(1) = 1.0d0
			g(2) = 1.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			g(1) = 2.0d0 * x(1)
			g(2) = 2.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  ZDT1

	if ( problem == 'ZDT1' ) then 
			
		if ( ind == 1 ) then
			g(1) = 1.0d0
			g(2:n) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			t = x(1) / faux
			
			g(1) = -0.5d0 * t ** (-0.5d0)
			do i = 2,n
				g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 - sqrt(t) / 2.0d0 )
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  ZDT2

	if ( problem == 'ZDT2' ) then 
			
		if ( ind == 1 ) then
			g(1) = 1.0d0
			g(2:n) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			t = x(1) / faux
			
			g(1) = -2.0d0 * t
			do i = 2,n
				g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 + t ** 2 )
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  ZDT3

	if ( problem == 'ZDT3' ) then 
			
		if ( ind == 1 ) then
			g(1) = 1.0d0
			g(2:n) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 1.0d0 + 9.0d0 * sum(x(2:n)) / ( n - 1 )
			t = x(1) / faux
			
			a = sin( 1.0d1 * pi * x(1) )
			
			g(1) = -0.5d0 * t ** (-0.5d0) - a - 1.0d1 * pi * x(1) * cos( 1.0d1 * pi * x(1) )
			do i = 2,n
				g(i) = 9.0d0 / ( n - 1 ) * ( 1.0d0 - 0.5d0 * sqrt(t) ) 
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  ZDT4

	if ( problem == 'ZDT4' ) then 
			
		if ( ind == 1 ) then
			g(1) = 1.0d0
			g(2:n) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then
			faux = 0.0d0
			do i = 2,n
				faux = faux + x(i) ** 2 - 1.0d1 * cos( 4.0d0 * pi * x(i) )
			end do
			faux = faux + 1.0d0 + 1.0d1 * ( n - 1 )
			t = x(1) / faux
			
			g(1) = -0.5d0 * t ** (-0.5d0)
			do i = 2,n
				a = sin( 4.0d0 * pi * x(i) )
				g(i) = ( 2.0d0 * x(i) + 4.0d1 * pi * a ) * ( 1.0d0 - sqrt(t) / 2.0d0 )
			end do
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  ZDT6

	if ( problem == 'ZDT6' ) then 
			
		if ( ind == 1 ) then			
			a = exp( -4.0d0 * x(1) )
			b = sin( 6.0d0 * pi * x(1) )
			
			g(1) = 4.0d0 * a * b ** 6 - 3.6d1 * pi * a * b ** 5 * cos( 6.0d0 * pi * x(1) )
			g(2:n) = 0.0d0
			if ( penalize ) call penalizeg(n,x,g)
			return
		end if
		
		if ( ind == 2 ) then			
			a = exp( -4.0d0 * x(1) )
			b = sin( 6.0d0 * pi * x(1) )
			gaux1 = 1.0d0 - a * b ** 6 
			gaux2 = 1.0d0 + 9.0d0 * ( sum(x(2:n)) / ( n - 1 ) ) ** 0.25d0
			t = gaux1 / gaux2
			
			A1 = 9.0d0 * 0.25d0 * ( sum(x(2:n)) / ( n - 1 ) ) ** (-0.75d0) / ( n - 1 )
			
			g(1) = -2.0d0 * t * ( 4.0d0 * a * b ** 6 - 3.6d1 * pi * a * b ** 5 * cos( 6.0d0 * pi * x(1) ) )
			do i = 2,n
				g(i) = A1 * ( 1.0d0 + t ** 2 )
			end do
			
			if ( penalize ) call penalizeg(n,x,g)
			return	
		end if
			
	end if		
	
! ----------------------------------------------------------------------

	!  ZLT1

	if ( problem == 'ZLT1' ) then 
		
		g(ind) = 2.0d0 * ( x(ind) - 1.0d0 )
		do i = 1,n
			if ( i /= ind ) g(i) = 2.0d0 * x(i)
		end do
			
	end if			
				

	end subroutine evalg
	
!***********************************************************************
!***********************************************************************

	subroutine quadstep(n,x,d,stp,ind,flag)
	
	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in) :: n,ind
	real(kind=8), intent(out) :: stp
	integer, intent(out) :: flag
	
	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: x(n),d(n)
		
	!     If quadratic(ind) is set equal to true in function inip, then
	!     YOU MUST to provide the minimizer x of the associated function.
	
	flag = -1
	
	end subroutine quadstep 
	

!***********************************************************************
!***********************************************************************	
	
	subroutine penalizef(n,x,f)
	
	implicit none

	! SCALAR ARGUMENTS
	integer,intent(in) :: n
	real(kind=8),intent(inout) :: f
	
	! ARRAY ARGUMENTS
	real(kind=8),intent(in) :: x(n)
	
	! LOCAL SCALARS
	integer :: i
	
	do i = 1,n
		f = f + penparam/3.0d0 * ( max( 0.0d0, x(i) - uin(i) )**3 + max( 0.0d0, - x(i) + lin(i) )**3  )
	end do

	end subroutine penalizef
	
!***********************************************************************
!***********************************************************************	
	
	subroutine penalizeg(n,x,g)
	
	implicit none

	! SCALAR ARGUMENTS
	integer,intent(in) :: n
	
	! ARRAY ARGUMENTS
	real(kind=8),intent(in) :: x(n)
	real(kind=8),intent(inout) :: g(n)
	
	! LOCAL SCALARS
	integer :: i
	
	do i = 1,n
		g(i) = g(i) + penparam * ( max( 0.0d0, x(i) - uin(i) )**2 - max( 0.0d0, - x(i) + lin(i) )**2 )
	end do

	end subroutine penalizeg

end module myproblem
