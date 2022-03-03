! ******************************************************************
! ******************************************************************

program algencanma

  implicit none

  ! LOCAL SCALARS
  logical :: checkder
  integer :: hnnzmax,i,inform,jcnnzmax,m,n,npairs,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,seed,snorm

  ! LOCAL ARRAYS
  character(len=80) :: specfnm,outputfnm,vparam(10)
  logical :: coded(11)
  logical, pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! COMMON SCALARS
  integer :: p

  ! COMMON BLOCKS
  common /pdata/ p

  ! FUNCTIONS
  real(kind=8) :: drand

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,myevalfc, &
       myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! Problem data

  p = 3
  npairs = p * ( p - 1 ) / 2

  ! Number of variables

  n = 2 * p + 2

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n))

  l(1:2*p) = - 1.0d+20
  u(1:2*p) =   1.0d+20

  l(2*p+1:n) = 0.0d0
  u(2*p+1:n) = 1.0d+20

  seed = 654321.0d0
  do i = 1,2*p
     x(i) = - 5.0d0 + 10.0d0 * drand(seed)
  end do
  x(2*p+1:n) = 10.0d0

  ! Constraints

  m = npairs + 4 * p

  allocate(equatn(m),linear(m),lambda(m))

  equatn(1:m) = .false.

  linear(1:npairs) = .false.
  linear(npairs+1:m) = .true.

  lambda(1:m) = 0.0d0

  ! Coded subroutines

  coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
  coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = 4 * npairs + 2 * ( 4 * p )
  hnnzmax  = 1 + 6 * npairs + 10 * npairs + 3 * ( 4 * p )

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

  f = x(n-1) * x(n)

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

  g(1:n-2) = 0.0d0
  g(n-1)   = x(n)
  g(n)     = x(n-1)

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

  hnnz = 1

  if ( hnnz .gt. lim ) then
     lmem = .true.
     return
  end if

  hrow(1) = n
  hcol(1) = n - 1
  hval(1) = 1.0d0

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
  integer :: p

  ! LOCAL SCALARS
  integer :: i,j

  ! COMMON BLOCKS
  common /pdata/ p

  ! Compute ind-th constraint

  flag = 0

  if ( 1 .le. ind .and. ind .le. p * ( p - 1 ) / 2 ) then
     call pair(p,ind,i,j)
     c = ( dble(i) + dble(j) ) ** 2 - ( x(2*i-1) - x(2*j-1) ) ** 2 - ( x(2*i) - x(2*j) ) ** 2

  else if ( ind .le. p * ( p - 1 ) / 2 + p ) then
     i = ind - p * ( p - 1 ) / 2
     c = - 0.5d0 * x(n-1) + dble(i) - x(2*i-1)

  else if ( ind .le. p * ( p - 1 ) / 2 + 2 * p ) then
     i = ind - p * ( p - 1 ) / 2 - p
     c = - 0.5d0 * x(n-1) + dble(i) + x(2*i-1)

  else if ( ind .le. p * ( p - 1 ) / 2 + 3 * p ) then
     i = ind - p * ( p - 1 ) / 2 - 2 * p
     c = - 0.5d0 * x(n) + dble(i) - x(2*i)

  else
     i = ind - p * ( p - 1 ) / 2 - 3 * p
     c = - 0.5d0 * x(n) + dble(i) + x(2*i)
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
  integer :: p

  ! LOCAL SCALARS
  integer :: i,j

  ! COMMON BLOCKS
  common /pdata/ p

  ! Compute gradient of the ind-th constraint

  flag = 0
  lmem = .false.

  if ( 1 .le. ind .and. ind .le. p * ( p - 1 ) / 2 ) then
     call pair(p,ind,i,j)

     jcnnz = 4

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = 2 * i - 1
     jcval(1) = - 2.0d0 * ( x(2*i-1) - x(2*j-1) )
     jcvar(2) = 2 * j - 1
     jcval(2) =   2.0d0 * ( x(2*i-1) - x(2*j-1) )
     jcvar(3) = 2 * i
     jcval(3) = - 2.0d0 * ( x(2*i) - x(2*j) )
     jcvar(4) = 2 * j
     jcval(4) =   2.0d0 * ( x(2*i) - x(2*j) )

  else if ( ind .le. p * ( p - 1 ) / 2 + p ) then
     i = ind - p * ( p - 1 ) / 2

     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = n - 1
     jcval(1) = - 0.5d0
     jcvar(2) = 2 * i - 1
     jcval(2) = - 1.0d0

  else if ( ind .le. p * ( p - 1 ) / 2 + 2 * p ) then
     i = ind - p * ( p - 1 ) / 2 - p

     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = n - 1
     jcval(1) = - 0.5d0
     jcvar(2) = 2 * i - 1
     jcval(2) =   1.0d0

  else if ( ind .le. p * ( p - 1 ) / 2 + 3 * p ) then
     i = ind - p * ( p - 1 ) / 2 - 2 * p

     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if
     jcvar(1) = n
     jcval(1) = - 0.5d0
     jcvar(2) = 2 * i
     jcval(2) = - 1.0d0

  else
     i = ind - p * ( p - 1 ) / 2 - 3 * p

     jcnnz = 2

     if ( jcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     jcvar(1) = n
     jcval(1) = - 0.5d0
     jcvar(2) = 2 * i
     jcval(2) =   1.0d0
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
  integer :: p

  ! LOCAL SCALARS
  integer :: i,j

  ! COMMON BLOCKS
  common /pdata/ p

  ! Compute lower-triangle of the ind-th constraint's Hessian

  flag = 0
  lmem = .false.

  if ( 1 .le. ind .and. ind .le. p * ( p - 1 ) / 2 ) then
     call pair(p,ind,i,j)

     hcnnz = 6

     if ( hcnnz .gt. lim ) then
        lmem = .true.
        return
     end if

     hcrow(1) = 2 * i - 1
     hccol(1) = 2 * i - 1
     hcval(1) = - 2.0d0

     hcrow(2) = 2 * i
     hccol(2) = 2 * i
     hcval(2) = - 2.0d0

     hcrow(3) = 2 * j - 1
     hccol(3) = 2 * j - 1
     hcval(3) = - 2.0d0

     hcrow(4) = 2 * j
     hccol(4) = 2 * j
     hcval(4) = - 2.0d0

     hcrow(5) = 2 * j - 1
     hccol(5) = 2 * i - 1
     hcval(5) =   2.0d0

     hcrow(6) = 2 * j
     hccol(6) = 2 * i
     hcval(6) =   2.0d0

  else
     hcnnz = 0
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

! ******************************************************************
! ******************************************************************

subroutine pair(p,ind,i,j)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: p,ind
  integer, intent(out) :: i,j

  ! LOCAL SCALARS
  integer :: k

  i = 1
  k = ind
10 if ( k .gt. p - i ) then
     k = k - ( p - i )
     i = i + 1
     go to 10
  end if
  j = i + k

end subroutine pair

! ******************************************************************
! ******************************************************************

subroutine drawsol(n,x)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! COMMON SCALARS
  integer :: p

  ! LOCAL SCALARS
  integer :: i
  real(kind=8) :: scale

  ! COMMON BLOCKS
  common /pdata/ p

  open(unit=10,file='solution.mp')

  ! SCALING
  scale = 10.0d0 / ( max( x(n-1), x(n) ) )
  write(10,10) scale

  ! RECTANGULAR CONTAINER
  write(10,20) -0.5d0*x(n-1),-0.5d0*x(n), 0.5d0*x(n-1),-0.5d0*x(n), &
       0.5d0*x(n-1), 0.5d0*x(n),-0.5d0*x(n-1), 0.5d0*x(n)

  ! CIRCULAR ITEMS
  do i = 1,p
     write(10,30) 2.0d0*dble(i),2.0d0*dble(i),x(2*i-1),x(2*i)
  end do

  write(10,40)

  close(10)

  ! NON-EXECUTABLE STATEMENTS

10 format('beginfig(1);'/,'u = ',f20.10,' cm;') 
20 format('draw (',f20.10,'u,',f20.10,'u)--', &
          '     (',f20.10,'u,',f20.10,'u)--', &
          '     (',f20.10,'u,',f20.10,'u)--', &
          '     (',f20.10,'u,',f20.10,'u)--cycle;')
30 format('draw fullcircle', &
          /,'     xscaled  ',f20.10,'u', &
          /,'     yscaled  ',f20.10,'u', &
          /,'     shifted (',f20.10,'u,',f20.10,'u);')
40 format('endfig;',/,'end;') 

end subroutine drawsol

