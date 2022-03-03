! ******************************************************************
! ******************************************************************

!subroutine evaldN(nvar,nf,d,multLg,inform)
subroutine evaldN(nvar,nf,d,inform)

	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in)      :: nvar,nf
	integer, intent(out)     :: inform

	! ARRAY ARGUMENTS
	real(kind=8), intent(out)  :: d(nvar)!,multLg(nf)

	! LOCAL SCALARS
	logical      :: checkder
	integer      :: n,m,allocerr,hnnzmax,jcnnzmax,nvparam,info,alinfo
	real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
	   f,nlpsupn,snorm

	! LOCAL ARRAYS
	character(len=80)     :: specfnm,outputfnm,vparam(10)
	logical               :: coded(11)
	logical,      pointer :: equatn(:),linear(:)
	real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

	! EXTERNAL SUBROUTINES
	external :: myevalfN,myevalgN,myevalhN,myevalcN,myevaljacN,myevalhcN, &
			  myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp

	! This routine calculates the steepest descend direction by solving
	!
	! Minimize     tau
	! subject to   [JF(x)(u-x)]_i <= tau, i=1,...,m.
	!              u \in C
	!            
	! inform:
	!
	! 0: Feasibility, optimality and complementarity satisfied
	! 2: Too large penalty parameter. Infeasible problem?
	! 3: Maximum number of algencan iterations reached

	! Number of variables

	n = nvar + 1

	! Set lower bounds, upper bounds, and initial guess

	allocate(x(n),l(n),u(n),stat=allocerr)
	if ( allocerr .ne. 0 ) then
		write(*,*) 'Allocation error in main program'
		stop
	end if

	
	l(:) = - 1.0d+20
	u(:) =   1.0d+20

	x(:) = 0.0d0

	! Constraints

	m = nf

	allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
	if ( allocerr .ne. 0 ) then
		write(*,*) 'Allocation error in main program'
		stop
	end if

	equatn(1:nf) = .false.
	linear(1:nf) = .false.

	lambda(1:m) = 0.0d0

	! Coded subroutines

	coded(1:6)  = .true.  ! fsub, gsub, hsub, csub, jacsub, hcsub
	coded(7:11) = .false. ! fcsub,gjacsub,gjacpsub,hlsub,hlpsub

	! Upper bounds on the number of sparse-matrices non-null elements

	jcnnzmax = n * m
	hnnzmax  = m * nvar * ( nvar + 1 ) / 2

	! Checking derivatives?

	checkder = .false.

	! Parameters setting

	epsfeas   = 1.0d-08
	epsopt    = 1.0d-08

	efstain   = sqrt( epsfeas )
	eostain   = epsopt ** 1.5d0

	efacc     = sqrt( epsfeas )
	eoacc     = sqrt( epsopt )

	outputfnm = ''
	specfnm   = ''
	
	nvparam = 2
	vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 00'
	vparam(2) = 'OUTER-ITERATIONS-LIMIT 100'

	! Call Algencan

	! E. G. Birgin and J. M. Martı́nez, Practical augmented Lagrangian 
	! methods for constrained optimization, SIAM, 2014.

	call algencan(myevalfN,myevalgN,myevalhN,myevalcN,myevaljacN,myevalhcN, &
	   myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
	   hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
	   specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
	   checkder,f,cnorm,snorm,nlpsupn,info,alinfo)
	   	   
	inform = alinfo    
		 
	d = x(1:nvar)       
!	multLg = lambda

	deallocate(x,l,u,lambda,equatn,linear,stat=allocerr)
	if ( allocerr .ne. 0 ) then
		write(*,*) 'Deallocation error in main program'
		stop
	end if

end subroutine evaldN

! ******************************************************************
! ******************************************************************

subroutine myevalfN(n,x,f,flag)

	implicit none

	! SCALAR ARGUMENTS
	integer,      intent(in)  :: n
	integer,      intent(out) :: flag
	real(kind=8), intent(out) :: f

	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: x(n)

	flag = 0	

	f = x(n)  

end subroutine myevalfN

! ******************************************************************
! ******************************************************************

subroutine myevalgN(n,x,g,flag)

	implicit none

	! SCALAR ARGUMENTS
	integer,      intent(in)  :: n
	integer,      intent(out) :: flag

	! ARRAY ARGUMENTS
	real(kind=8), intent(in)  :: x(n)
	real(kind=8), intent(out) :: g(n)

	flag = 0

	g(1:n-1) = 0.0d0

	g(n) = 1.0d0

end subroutine myevalgN

! ******************************************************************
! ******************************************************************

subroutine myevalhN(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

	implicit none

	! SCALAR ARGUMENTS
	logical,      intent(out) :: lmem
	integer,      intent(in)  :: lim,n
	integer,      intent(out) :: flag,hnnz

	! ARRAY ARGUMENTS
	integer,      intent(out) :: hcol(lim),hrow(lim)
	real(kind=8), intent(in)  :: x(n)
	real(kind=8), intent(out) :: hval(lim)

	flag = 0
	lmem = .false.

	hnnz = 0
  
end subroutine myevalhN

! ******************************************************************
! ******************************************************************

subroutine myevalcN(n,x,ind,c,flag)

	use globals, only: JF, B

	implicit none

	! SCALAR ARGUMENTS
	integer,      intent(in)  :: ind,n
	integer,      intent(out) :: flag
	real(kind=8), intent(out) :: c

	! ARRAY ARGUMENTS
	real(kind=8), intent(in)  :: x(n)
	
	! LOCAL SCALAR
	integer :: nvar
	
	flag = 0
	
	nvar = n - 1

	c = dot_product( JF(ind,:),x(1:nvar) ) &
	    + 0.5d0 * dot_product( matmul( B(:,:,ind),x(1:nvar) ),x(1:nvar) ) - x(n)

	end subroutine myevalcN

! ******************************************************************
! ******************************************************************

subroutine myevaljacN(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

	use globals, only: JF, B

	implicit none

	! SCALAR ARGUMENTS
	logical, intent(out) :: lmem
	integer, intent(in)  :: ind,lim,n
	integer, intent(out) :: flag,jcnnz

	! ARRAY ARGUMENTS
	integer,      intent(out) :: jcvar(lim)
	real(kind=8), intent(in)  :: x(n)
	real(kind=8), intent(out) :: jcval(lim)

	! LOCAL SCALAR
	integer :: i, nvar
	
	flag = 0
	
	nvar = n - 1

	lmem = .false.

	jcnnz = n

	if ( jcnnz .gt. lim ) then
		lmem = .true.
		return
	end if

	do i = 1,nvar
		jcvar(i) = i
		jcval(i) = JF(ind,i) +  dot_product( B(i,:,ind),x(1:nvar) )
	end do

	jcvar(n) = n
	jcval(n) = - 1.0d0


end subroutine myevaljacN

! ******************************************************************
! ******************************************************************

subroutine myevalhcN(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

	use globals, only: B

	implicit none

	! SCALAR ARGUMENTS
	logical, intent(out) :: lmem
	integer, intent(in)  :: ind,lim,n
	integer, intent(out) :: flag,hcnnz

	! ARRAY ARGUMENTS
	integer,      intent(out) :: hccol(lim),hcrow(lim)
	real(kind=8), intent(in)  :: x(n)
	real(kind=8), intent(out) :: hcval(lim)
	
	! LOCAL SCALAR
	integer :: i,j,k,nvar
	
	flag = 0
	
	nvar = n - 1

	lmem = .false.
	
	hcnnz = nvar * ( nvar + 1 ) / 2
	
	if ( hcnnz .gt. lim ) then
		lmem = .true.
		return
	end if
	
	k = 0
	do i = 1,nvar
		do j = 1,i
			k = k + 1
			hcrow(k) = i
			hccol(k) = j
			hcval(k) = B(i,j,ind)
		end do
	end do

  end subroutine myevalhcN

! ******************************************************************
! ******************************************************************

subroutine myevalfc(n,x,f,m,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: m,n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: c(m)

  flag = - 1

end subroutine myevalfc

! ******************************************************************
! ******************************************************************

subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,m,n
  integer,      intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: jcfun(lim),jcvar(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n),jcval(lim)

  flag = - 1

end subroutine myevalgjac

! ******************************************************************
! ******************************************************************

subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,   intent(inout) :: gotj
  integer,   intent(in)    :: m,n
  integer,   intent(out)   :: flag
  character, intent(in)    :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)    :: x(n)
  real(kind=8), intent(inout) :: p(m),q(n)
  real(kind=8), intent(out)   :: g(n)

  flag = - 1

end subroutine myevalgjacp

! ******************************************************************
! ******************************************************************

subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,m,n
  integer,      intent(out) :: flag,hlnnz
  real(kind=8), intent(in)  :: sf

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hlcol(lim),hlrow(lim)
  real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
  real(kind=8), intent(out) :: hlval(lim)

  flag = - 1

end subroutine myevalhl

! ******************************************************************
! ******************************************************************

subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(inout) :: goth
  integer,      intent(in)    :: m,n
  integer,      intent(out)   :: flag
  real(kind=8), intent(in)    :: sf

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
  real(kind=8), intent(out) :: hp(n)

  flag = - 1

end subroutine myevalhlp

!! ******************************************************************
!! ******************************************************************

!subroutine evaldSD(nvar,nf,lC,uC,d,multLg,inform)

!	use globals, only: xinner

!	implicit none

!	! SCALAR ARGUMENTS
!	integer, intent(in)      :: nvar,nf
!	integer, intent(out)     :: inform

!	! ARRAY ARGUMENTS
!	real(kind=8), intent(out) :: d(nvar),multLg(nf)
!	real(kind=8), intent(in)  :: lC(nvar),uC(nvar)

!	! LOCAL SCALARS
!	logical      :: checkder
!	integer      :: n,m,allocerr,hnnzmax,jcnnzmax,nvparam,alinfo
!	real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
!	   f,nlpsupn,snorm

!	! LOCAL ARRAYS
!	character(len=80)     :: specfnm,outputfnm,vparam(10)
!	logical               :: coded(11)
!	logical,      pointer :: equatn(:),linear(:)
!	real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

!	! EXTERNAL SUBROUTINES
!	external :: myevalfSD,myevalgSD,myevalhSD,myevalcSD,myevaljacSD,myevalhcSD, &
!			  myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp

!	! This routine calculates the steepest descend direction by solving
!	!
!	! Minimize     alpha + 1/2 ||d||^2
!	! subject to   [JF(x)d]_i <= alpha, i=1,...,m. 
!	!
!	!            
!	! inform:
!	!
!	! 0: Feasibility, optimality and complementarity satisfied
!	! 2: Too large penalty parameter. Infeasible problem?
!	! 3: Maximum number of algencan iterations reached


!	! Number of variables

!	n = nvar + 1

!	! Set lower bounds, upper bounds, and initial guess

!	allocate(x(n),l(n),u(n),stat=allocerr)
!	if ( allocerr .ne. 0 ) then
!	 write(*,*) 'Allocation error in main program'
!	 stop
!	end if
	
!	l(1:nvar) = lC(:) - xinner(:)
!	u(1:nvar) = uC(:) - xinner(:)

!	l(n) = - 1.0d+20
!	u(n) =   1.0d+20
	
!	l(:) = - 1.0d+20
!	u(:) =   1.0d+20

!	x(:) = 0.0d0

!	! Constraints

!	m = nf

!	allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
!	if ( allocerr .ne. 0 ) then
!	 write(*,*) 'Allocation error in main program'
!	 stop
!	end if

!	equatn(1:m) = .false.
!	linear(1:m) = .true.
	
!	lambda(1:m) = 0.0d0

!	! Coded subroutines

!	coded(1:6)  = .true.  ! fsub, gsub, hsub, csub, jacsub, hcsub
!	coded(7:11) = .false. ! fcsub,gjacsub,gjacpsub,hlsub,hlpsub

!	! Upper bounds on the number of sparse-matrices non-null elements

!	jcnnzmax = n * m
!	hnnzmax  = n - 1

!	! Checking derivatives?

!	checkder = .false.

!	! Parameters setting

!	epsfeas   = 1.0d-10
!	epsopt    = 1.0d-08

!	efstain   = sqrt( epsfeas )
!	eostain   = epsopt ** 1.5d0

!	efacc     = sqrt( epsfeas )
!	eoacc     = sqrt( epsopt )

!	outputfnm = ''
!	specfnm   = ''

!	nvparam = 2
!	vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 00'
!	vparam(2) = 'OUTER-ITERATIONS-LIMIT 100'

!	! Call Algencan

!	! E. G. Birgin and J. M. Martı́nez, Practical augmented Lagrangian 
!	! methods for constrained optimization, SIAM, 2014.

!	call algencan(myevalfSD,myevalgSD,myevalhSD,myevalcSD,myevaljacSD,myevalhcSD, &
!	   myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
!	   hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
!	   specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
!	   checkder,f,cnorm,snorm,nlpsupn,inform,alinfo)
	   
!	inform = alinfo     
		 
!	d = x(1:nvar)       
!	multLg = lambda

!	deallocate(x,l,u,lambda,equatn,linear,stat=allocerr)
!	if ( allocerr .ne. 0 ) then
!	 write(*,*) 'Deallocation error in main program'
!	 stop
!	end if

!end subroutine evaldSD

!! ******************************************************************
!! ******************************************************************

!subroutine myevalfSD(n,x,f,flag)

!  implicit none

!  ! SCALAR ARGUMENTS
!  integer,      intent(in)  :: n
!  integer,      intent(out) :: flag
!  real(kind=8), intent(out) :: f

!  ! ARRAY ARGUMENTS
!  real(kind=8), intent(in) :: x(n)

!  flag = 0

!  f = 0.5d0 * dot_product(x(1:n-1),x(1:n-1)) + x(n)

!end subroutine myevalfSD

!! ******************************************************************
!! ******************************************************************

!subroutine myevalgSD(n,x,g,flag)

!	implicit none

!	! SCALAR ARGUMENTS
!	integer,      intent(in)  :: n
!	integer,      intent(out) :: flag

!	! ARRAY ARGUMENTS
!	real(kind=8), intent(in)  :: x(n)
!	real(kind=8), intent(out) :: g(n)

!	! LOCAL SCALAR
!	integer :: i

!	flag = 0

!	do i = 1,n-1
!		g(i) = x(i)
!	end do

!	g(n) = 1.0d0

!end subroutine myevalgSD

!! ******************************************************************
!! ******************************************************************

!subroutine myevalhSD(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

!	implicit none

!	! SCALAR ARGUMENTS
!	logical,      intent(out) :: lmem
!	integer,      intent(in)  :: lim,n
!	integer,      intent(out) :: flag,hnnz

!	! ARRAY ARGUMENTS
!	integer,      intent(out) :: hcol(lim),hrow(lim)
!	real(kind=8), intent(in)  :: x(n)
!	real(kind=8), intent(out) :: hval(lim)

!	! LOCAL SCALAR
!	integer :: i

!	flag = 0
!	lmem = .false.

!	hnnz = n-1

!	if ( hnnz .gt. lim ) then
!		lmem = .true.
!		return 	
!	end if

!	do i = 1,n-1
!		hrow(i) = i
!		hcol(i) = i
!		hval(i) = 1.0d0
!	end do

!end subroutine myevalhSD

!! ******************************************************************
!! ******************************************************************

!subroutine myevalcSD(n,x,ind,c,flag)

!	use globals, only: JF 

!	implicit none

!	! SCALAR ARGUMENTS
!	integer,      intent(in)  :: ind,n
!	integer,      intent(out) :: flag
!	real(kind=8), intent(out) :: c

!	! ARRAY ARGUMENTS
!	real(kind=8), intent(in)  :: x(n)

!	flag = 0 
!	c = dot_product(JF(ind,:),x(1:n-1))-x(n)

!end subroutine myevalcSD

!! ******************************************************************
!! ******************************************************************

!subroutine myevaljacSD(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

!	use globals, only: JF

!	implicit none

!	! SCALAR ARGUMENTS
!	logical, intent(out) :: lmem
!	integer, intent(in)  :: ind,lim,n
!	integer, intent(out) :: flag,jcnnz

!	! ARRAY ARGUMENTS
!	integer,      intent(out) :: jcvar(lim)
!	real(kind=8), intent(in)  :: x(n)
!	real(kind=8), intent(out) :: jcval(lim)

!	! LOCAL SCALAR
!	integer :: i

!	! LOCAL ARRAY 
!	real(kind=8) :: g(n-1)

!	flag = 0
!	lmem = .false.

!	jcnnz = n

!	if ( jcnnz .gt. lim ) then
!		lmem = .true.
!		return
!	end if

!	do i = 1,n-1
!			jcvar(i) = i
!			jcval(i) = JF(ind,i)
!	end do

!	jcvar(n) = n
!	jcval(n) = - 1.0d0

!end subroutine myevaljacSD

!! ******************************************************************
!! ******************************************************************

!subroutine myevalhcSD(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

!	implicit none

!	! SCALAR ARGUMENTS
!	logical, intent(out) :: lmem
!	integer, intent(in)  :: ind,lim,n
!	integer, intent(out) :: flag,hcnnz

!	! ARRAY ARGUMENTS
!	integer,      intent(out) :: hccol(lim),hcrow(lim)
!	real(kind=8), intent(in)  :: x(n)
!	real(kind=8), intent(out) :: hcval(lim)

!	flag = 0
!	lmem = .false.

!	hcnnz = 0

!end subroutine myevalhcSD
