! ******************************************************************
! ******************************************************************

program algencanma

  implicit none

  ! LOCAL SCALARS
  logical :: checkder
  integer :: hnnzmax,i,inform,jcnnzmax,m,n,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,snorm

  ! LOCAL ARRAYS
  character(len=80) :: specfnm,outputfnm,vparam(10)
  logical :: coded(11)
  logical, pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! COMMON SCALARS
  integer :: ncirc

  ! COMMON BLOCKS
  common /pdata/ ncirc

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,myevalfc, &
       myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! Problem data

  ncirc = 91

  ! Number of variables

  n = 2

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n))

  l(1:n) = - 1.0d+20
  u(1:n) =   1.0d+20

  x(1:n) =   0.0d0

  ! Constraints

  m = 2 + ncirc

  allocate(equatn(m),linear(m),lambda(m))

  equatn(1:m) = .false.

  linear(1:2) = .true.
  linear(3:m) = .false.

  lambda(1:m) = 0.0d0

  ! Coded subroutines

  coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
  coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = 2 * m
  hnnzmax  = 2 * ncirc + 3 * m

  ! Checking derivatives?

  checkder = .false.

  ! Parameters setting

  epsfeas   =   1.0d-08
  epsopt    =   1.0d-08

  efstain   =   1.0d+20
  eostain   = - 1.0d+20

  efacc     = - 1.0d+20
  eoacc     = - 1.0d+20

  outputfnm = ''
  specfnm   = ''

  nvparam   = 0

  ! Optimize

  call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,  &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax,  &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,  &
       checkder,f,cnorm,snorm,nlpsupn,inform)

  deallocate(x,l,u,lambda,equatn,linear)

  stop

end program algencanma

! ******************************************************************
! ******************************************************************

subroutine myevalf(n,x,f,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! Compute objective function

  flag = 0

  f = x(1) + x(2)

end subroutine myevalf

! ******************************************************************
! ******************************************************************

subroutine myevalg(n,x,g,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: g(n)

  ! Compute gradient of the objective function

  flag = 0

  g(1:2) = 1.0d0

end subroutine myevalg

! ******************************************************************
! ******************************************************************

subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,n
  integer, intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  ! Compute (lower triangle of the) Hessian of the objective function

  flag = 0
  lmem = .false.

  hnnz = 0

end subroutine myevalh

! ******************************************************************
! ******************************************************************

subroutine myevalc(n,x,ind,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: ind,n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: c

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! COMMON SCALARS
  integer :: ncirc

  ! LOCAL SCALARS
  integer :: icirc
  real(kind=8) :: alpha,r

  ! COMMON BLOCKS
  common /pdata/ ncirc

  ! Compute ind-th constraint

  flag = 0

  if ( ind .eq. 1 ) then
     c = x(1) - 3.0d0 * x(2) - 1.0d0

  else if ( ind .eq. 2 ) then
     c = - x(1) + x(2) - 1.0d0

  else ! if ( 3 .le. ind .and. ind .le. ncirc + 2 ) then
     icirc = ind - 2

     alpha = ( icirc - 1 ) * 1.57079632679d0 / ( ncirc - 1 )

     r = 5.0d0

     c = ( x(1) - ( 1.0d0 - r ) * cos( alpha ) ) ** 2 &
       + ( x(2) - ( 1.0d0 - r ) * sin( alpha ) ) ** 2 - r ** 2
  end if

end subroutine myevalc

! ******************************************************************
! ******************************************************************

subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: jcvar(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: jcval(lim)

  ! COMMON SCALARS
  integer :: ncirc

  ! LOCAL SCALARS
  integer :: icirc
  real(kind=8) :: alpha,r

  ! COMMON BLOCKS
  common /pdata/ ncirc

  ! Compute gradient of the ind-th constraint

  flag = 0
  lmem = .false.

  if ( ind .eq. 1 ) then
     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = 1
     jcval(1) =   1.0d0
     jcvar(2) = 2
     jcval(2) = - 3.0d0

  else if ( ind .eq. 2 ) then
     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = 1
     jcval(1) = - 1.0d0
     jcvar(2) = 2
     jcval(2) =   1.0d0

  else ! if ( 3 .le. ind .and. ind .le. ncirc + 2 ) then
     icirc = ind - 2

     alpha = ( icirc - 1 ) * 1.57079632679d0 / ( ncirc - 1 )

     r = 5.0d0

     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = 1
     jcval(1) = 2.0d0 * ( x(1) - ( 1.0d0 - r ) * cos( alpha ) )
     jcvar(2) = 2
     jcval(2) = 2.0d0 * ( x(2) - ( 1.0d0 - r ) * sin( alpha ) )
  end if

end subroutine myevaljac

! ******************************************************************
! ******************************************************************

subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,hcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hccol(lim),hcrow(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hcval(lim)

  ! COMMON SCALARS
  integer :: ncirc

  ! LOCAL SCALARS
  integer :: icirc
  real(kind=8) :: alpha,r

  ! COMMON BLOCKS
  common /pdata/ ncirc

  ! Compute lower-triangle of the ind-th constraint's Hessian

  flag = 0
  lmem = .false.

  if ( ind .eq. 1 ) then
     hcnnz = 0

  else if ( ind .eq. 2 ) then
     hcnnz = 0

  else ! if ( 3 .le. ind .and. ind .le. ncirc + 2 ) then
     hcnnz = 2

     if ( hcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     hcrow(1) = 1
     hccol(1) = 1
     hcval(1) = 2.0d0
     hcrow(2) = 2
     hccol(2) = 2
     hcval(2) = 2.0d0
  end if


end subroutine myevalhc

! ******************************************************************
! ******************************************************************

subroutine myevalfc(n,x,f,m,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: c(m)

  flag = - 1

end subroutine myevalfc

! ******************************************************************
! ******************************************************************

subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,m,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: jcfun(lim),jcvar(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: g(n),jcval(lim)

  flag = - 1

end subroutine myevalgjac

! ******************************************************************
! ******************************************************************

subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: gotj
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  character, intent(in) :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(inout) :: p(m),q(n)
  real(kind=8), intent(out) :: g(n)

  flag = - 1

end subroutine myevalgjacp

! ******************************************************************
! ******************************************************************

subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,m,n
  integer, intent(out) :: flag,hlnnz
  real(kind=8), intent(in) :: sf

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hlcol(lim),hlrow(lim)
  real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
  real(kind=8), intent(out) :: hlval(lim)

  flag = - 1

end subroutine myevalhl

! ******************************************************************
! ******************************************************************

subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: goth
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  real(kind=8), intent(in) :: sf

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: lambda(m),p(n),sc(m),x(n)
  real(kind=8), intent(out) :: hp(n)

  flag = - 1

end subroutine myevalhlp
