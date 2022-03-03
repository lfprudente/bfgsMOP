! ******************************************************************
! ******************************************************************

program algencanma

  implicit none

  ! LOCAL SCALARS
  logical      :: checkder
  integer      :: allocerr,hnnzmax,inform,jcnnzmax,m,n,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,snorm

  ! LOCAL ARRAYS
  character(len=15)     :: strtmp
  character(len=80)     :: specfnm,outputfnm,vparam(10)
  logical               :: coded(11)
  logical,      pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
              myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! Number of variables

  n = 2

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  l(1) = - 10.0d0
  u(1) =   10.0d0

  l(2) = - 1.0d+20
  u(2) =   1.0d+20

  x(1:2) = 0.0d0

  ! Constraints

  m = 2

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  equatn(1:m) = .false.
  lambda(1:m) = 0.0d0

  linear(1) = .false.
  linear(2) = .true.

  ! Coded subroutines

  coded(1:11) = .false.
  coded(7)  = .true. ! fcsub
  coded(9)  = .true. ! gjacpsub
  coded(11) = .true. ! hlpsub

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax  = 0
  hnnzmax   = 0

  ! Checking derivatives?

  checkder = .true.

  ! Parameters setting

  epsfeas   = 1.0d-08
  epsopt    = 1.0d-08

  efstain   = sqrt( epsfeas )
  eostain   = epsopt ** 1.5d0

  efacc     = sqrt( epsfeas )
  eoacc     = sqrt( epsopt )

  outputfnm = ''
  specfnm   = ''

  nvparam = 0

  ! Optimize

  call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
       checkder,f,cnorm,snorm,nlpsupn,inform)

  deallocate(x,l,u,lambda,equatn,linear,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in main program'
     stop
  end if

  stop

end program algencanma

! ******************************************************************
! ******************************************************************

subroutine myevalf(n,x,f,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  flag = - 1

end subroutine myevalf

! ******************************************************************
! ******************************************************************

subroutine myevalg(n,x,g,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: n
  integer,      intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n)

  flag = - 1

end subroutine myevalg

! ******************************************************************
! ******************************************************************

subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out) :: lmem
  integer,      intent(in)  :: lim,n
  integer,      intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  flag = - 1

end subroutine myevalh

! ******************************************************************
! ******************************************************************

subroutine myevalc(n,x,ind,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: ind,n
  integer,      intent(out) :: flag
  real(kind=8), intent(out) :: c

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)

  flag = - 1

end subroutine myevalc

! ******************************************************************
! ******************************************************************

subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in)  :: ind,lim,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: jcvar(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: jcval(lim)

  flag = - 1

end subroutine myevaljac

! ******************************************************************
! ******************************************************************

subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in)  :: ind,lim,n
  integer, intent(out) :: flag,hcnnz

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: hccol(lim),hcrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hcval(lim)

  flag = - 1

end subroutine myevalhc

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

  flag = 0

  ! Objective function

  f = x(2)

  ! Constraints

  c(1) = x(1) ** 2 + 1 - x(n)
  c(2) = 2.0d0 - x(1) - x(n)

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

  ! The meaning of argument work follows: work = 'j' or 'J' means that
  ! q is an input array argument and that p = Jacobian x q must be
  ! computed, while work = 't' or 'T' means that p is an input array
  ! argument and that q = Jacobian^t x p must be computed. Moreover, a
  ! capital letter (i.e. 'J' or 'T') means that the gradient g of the
  ! objective function must also be computed. A lower letter (i.e.
  ! 'j' or 't') means that only the product with the Jacobian is
  ! required. In the later case, input array argument g MUST NOT be
  ! modified nor referenced.

  flag = 0

  ! Gradient of the objective function

  if ( work .eq. 'J' .or. work .eq. 'T' ) then
     g(1) = 0.0d0
     g(2) = 1.0d0
  end if

  ! p = Jacobian x q
  if ( work .eq. 'j' .or. work .eq. 'J' ) then
     p(1) = 2.0d0 * x(1) * q(1) - q(2)
     p(2) = - q(1) - q(2)

     ! q = Jacobian^t x p
  else ! if ( work .eq. 't' .or. work .eq. 'T' ) then
     q(1) = 2.0d0 * x(1) * p(1) - p(2)
     q(2) = - p(1) - p(2)
  end if

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

  flag = 0

  hp(1) = 2.0d0 * sc(1) * lambda(1) * p(1)
  hp(2) = 0.0d0

end subroutine myevalhlp
