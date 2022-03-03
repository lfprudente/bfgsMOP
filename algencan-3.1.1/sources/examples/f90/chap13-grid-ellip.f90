! ******************************************************************
! ******************************************************************

program algencanma

  use modgrid

  implicit none

  ! LOCAL SCALARS
  logical :: checkder
  character :: answer
  integer :: allocerr,hnnzmax,i,inform,j,jcnnzmax,m,n,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,snorm,angle,posx,posx1,posx2,posy,posy1,posy2,a,b,step

  ! LOCAL ARRAYS
  character(len=80) :: specfnm,outputfnm,vparam(10)
  logical :: coded(11)
  logical, pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! FUNCTIONS
  integer :: II

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,myevalfc, &
       myevalgjac,myevalgjacp,myevalhl,myevalhlp

  write(*,*) 'Enter nabs and nord: '
  read(*,*) nabs,nord

  ! Number of variables

  n = 2 * nabs * nord

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  l(1:n) = - 1.0d+20
  u(1:n) =   1.0d+20

  do i = 1,nabs
     angle = ( 0.5d0 - ( dble(i) - 1 ) / ( nabs - 1 ) ) * 3.14159d0
     posx1 = cos(angle)
     posy1 = sin(angle)
     l(II(i,1,'x'))    = posx1 
     l(II(i,1,'y'))    = posy1 
     u(II(i,1,'x'))    = l(II(i,1,'x'))
     u(II(i,1,'y'))    = l(II(i,1,'y'))

     a = 6.0d0
     b = 3.0d0
     posy2 = sqrt( 1.0d0 / ( cos(angle)**2 / sin(angle)**2 / a**2 + 1.0d0/b**2 ) )
     if ( i .gt. nabs / 2 ) posy2 = - posy2
     if ( abs( posy2 ) .le. 1.0d-08 ) then
        posx2 = 6.0d0
     else
        posx2 = posy2 * cos(angle) / sin(angle)
     end if
     l(II(i,nord,'x')) = posx2
     l(II(i,nord,'y')) = posy2
     u(II(i,nord,'x')) = l(II(i,nord,'x'))
     u(II(i,nord,'y')) = l(II(i,nord,'y'))
  end do

  do j = 1,nord
     posy = 1.0d0 + 2.0d0 * ( dble(j) - 1 ) / ( nord - 1 )
     l(II(1,j,'x'))    = 0.0d0
     l(II(1,j,'y'))    = posy 
     u(II(1,j,'x'))    = l(II(1,j,'x'))
     u(II(1,j,'y'))    = l(II(1,j,'y'))

     l(II(nabs,j,'x')) = 0.0d0
     l(II(nabs,j,'y')) = - posy
     u(II(nabs,j,'x')) = l(II(nabs,j,'x'))
     u(II(nabs,j,'y')) = l(II(nabs,j,'y'))
  end do

  ! Set initial guess

  do i = 1,nabs
     x(II(i,1,'x')) = l(II(i,1,'x'))
     x(II(i,1,'y')) = l(II(i,1,'y'))

     x(II(i,nord,'x')) = l(II(i,nord,'x'))
     x(II(i,nord,'y')) = l(II(i,nord,'y'))

     do j = 2,nord - 1
        step = ( dble(j) - 1 ) / ( nord - 1 )
        x(II(i,j,'x')) = l(II(i,1,'x')) + ( l(II(i,nord,'x')) - l(II(i,1,'x')) ) * step
        x(II(i,j,'y')) = l(II(i,1,'y')) + ( l(II(i,nord,'y')) - l(II(i,1,'y')) ) * step
     end do
  end do

  ! Problem data

  write(*,*) 'Enter constants for fS and fA in the objective function: '
  read(*,*) sigmafS,sigmafA

  write(*,*) '[C]onstrained or [U]nconstrained? '
  read(*,*) answer

  if ( answer .eq. 'C' .or. answer .eq. 'c' ) then
     m = 1

     write(*,*) 'Enter constants for fS and fA in the constraint: '
     read(*,*) sigmacS,sigmacA

     write(*,*) 'Enter constant for the right-hand-side: '
     read(*,*) constc
  else
     m = 0
  end if

  ! Constraints

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  equatn(1:m) = .false.
  linear(1:m) = .false.
  lambda(1:m) = 0.0d0

  ! Coded subroutines

  coded(1:11) = .false.
  coded(7)    = .true.  ! evalfc
  coded(8)    = .true.  ! evalgjac
  coded(10)   = .true.  ! evalhl

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = n
  hnnzmax  = 6 * ( 2 * nabs * nord - nabs - nord ) + &
             36 * ( nabs - 1 ) * ( nord - 1 )

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
  vparam(1) = 'ITERATIONS-OUTPUT-DETAIL 11'
  vparam(2) = 'NEWTON-LINE-SEARCH-INNER-SOLVER'

  ! Optimize

  call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,  &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax,  &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,  &
       checkder,f,cnorm,snorm,nlpsupn,inform)

  call drawsol(n,x)

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
  integer, intent(in) :: n
  integer, intent(out) :: flag
  real(kind=8), intent(out)  :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  flag = - 1

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

  flag = - 1

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
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  flag = - 1

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

  flag = - 1

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

  flag = - 1

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

  flag = - 1

end subroutine myevalhc

! ******************************************************************
! ******************************************************************

subroutine myevalfc(n,x,f,m,c,flag)

  use modgrid

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: c(m)

  ! FUNCTIONS
  integer :: II

  ! LOCAL SCALARS
  integer :: i,j
  real(kind=8) :: area,f1,f2

  flag = 0

  f1 = 0.0d0
  do i = 1,nabs 
     do j = 1,nord
        if ( j .lt. nord ) then
           ! Distance from P(i,j) to P(i,j+1)
           f1 = f1 + ( x(II(i,j,'x')) - x(II(i,j+1,'x')) ) ** 2 &
                   + ( x(II(i,j,'y')) - x(II(i,j+1,'y')) ) ** 2
        end if
        if ( i .lt. nabs ) then
           ! Distance from P(i,j) to P(i+1,j)
           f1 = f1 + ( x(II(i,j,'x')) - x(II(i+1,j,'x')) ) ** 2 &
                   + ( x(II(i,j,'y')) - x(II(i+1,j,'y')) ) ** 2
        end if
     end do
  end do
  f1 = 0.5d0 * f1 / ( ( nabs - 1 ) * nord + ( nord - 1 ) * nabs )

  f2 = 0.0d0
  do i = 1,nabs - 1
     do j = 1,nord - 1
        area = x(II(i  ,j  ,'x')) * x(II(i+1,j  ,'y')) -  x(II(i+1,j  ,'x')) * x(II(i  ,j  ,'y')) &
             + x(II(i+1,j  ,'x')) * x(II(i+1,j+1,'y')) -  x(II(i+1,j+1,'x')) * x(II(i+1,j  ,'y')) &
             + x(II(i+1,j+1,'x')) * x(II(i  ,j+1,'y')) -  x(II(i  ,j+1,'x')) * x(II(i+1,j+1,'y')) &
             + x(II(i  ,j+1,'x')) * x(II(i  ,j  ,'y')) -  x(II(i  ,j  ,'x')) * x(II(i  ,j+1,'y'))

        f2 = f2 + area ** 2
     end do
  end do
  f2 = 0.5d0 * f2 / ( ( nabs - 1 ) * ( nord - 1 ) )

  f = sigmafS * f1 + sigmafA * f2

  if ( m .ne. 0 ) c(1) = sigmacS * f1 + sigmacA * f2 - constc

end subroutine myevalfc

! ******************************************************************
! ******************************************************************

subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

  use modgrid

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,m,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: jcfun(lim),jcvar(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: g(n),jcval(lim)

  ! FUNCTIONS
  integer :: II

  ! LOCAL SCALARS
  integer :: i,ind,j
  real(kind=8) :: area,diffx,diffy

  ! LOCAL ARRAYS
  real(kind=8) :: g1(n),g2(n)

  flag = 0
  lmem = .false.

  g1(1:n) = 0.0d0

  do i = 1,nabs 
     do j = 1,nord
        if ( j .lt. nord ) then
           ! Distance from P(i,j) to P(i,j+1)
           diffx = x(II(i,j,'x')) - x(II(i,j+1,'x'))
           g1(II(i,j,'x'))   = g1(II(i,j,'x'))   + diffx
           g1(II(i,j+1,'x')) = g1(II(i,j+1,'x')) - diffx

           diffy = x(II(i,j,'y')) - x(II(i,j+1,'y'))
           g1(II(i,j,'y'))   = g1(II(i,j,'y'))   + diffy
           g1(II(i,j+1,'y')) = g1(II(i,j+1,'y')) - diffy
        end if
        if ( i .lt. nabs ) then
           ! Distance from P(i,j) to P(i+1,j)
           diffx = x(II(i,j,'x')) - x(II(i+1,j,'x'))
           g1(II(i,j,'x'))   = g1(II(i,j,'x'))   + diffx
           g1(II(i+1,j,'x')) = g1(II(i+1,j,'x')) - diffx

           diffy = x(II(i,j,'y')) - x(II(i+1,j,'y'))
           g1(II(i,j,'y'))   = g1(II(i,j,'y'))   + diffy
           g1(II(i+1,j,'y')) = g1(II(i+1,j,'y')) - diffy
        end if
     end do
  end do

  g1(1:n) = g1(1:n) / ( ( nabs - 1 ) * nord + ( nord - 1 ) * nabs )

  g2(1:n) = 0.0d0

  do i = 1,nabs - 1
     do j = 1,nord - 1
        area = x(II(i  ,j  ,'x')) * x(II(i+1,j  ,'y')) -  x(II(i+1,j  ,'x')) * x(II(i  ,j  ,'y')) &
             + x(II(i+1,j  ,'x')) * x(II(i+1,j+1,'y')) -  x(II(i+1,j+1,'x')) * x(II(i+1,j  ,'y')) &
             + x(II(i+1,j+1,'x')) * x(II(i  ,j+1,'y')) -  x(II(i  ,j+1,'x')) * x(II(i+1,j+1,'y')) &
             + x(II(i  ,j+1,'x')) * x(II(i  ,j  ,'y')) -  x(II(i  ,j  ,'x')) * x(II(i  ,j+1,'y'))

        g2(II(i  ,j  ,'x')) = g2(II(i  ,j  ,'x')) + area * ( + x(II(i+1,j  ,'y')) - x(II(i  ,j+1,'y')) )
        g2(II(i  ,j  ,'y')) = g2(II(i  ,j  ,'y')) + area * ( - x(II(i+1,j  ,'x')) + x(II(i  ,j+1,'x')) )
        g2(II(i  ,j+1,'x')) = g2(II(i  ,j+1,'x')) + area * ( - x(II(i+1,j+1,'y')) + x(II(i  ,j  ,'y')) )
        g2(II(i  ,j+1,'y')) = g2(II(i  ,j+1,'y')) + area * ( + x(II(i+1,j+1,'x')) - x(II(i  ,j  ,'x')) )
        g2(II(i+1,j  ,'x')) = g2(II(i+1,j  ,'x')) + area * ( - x(II(i  ,j  ,'y')) + x(II(i+1,j+1,'y')) )
        g2(II(i+1,j  ,'y')) = g2(II(i+1,j  ,'y')) + area * ( + x(II(i  ,j  ,'x')) - x(II(i+1,j+1,'x')) )
        g2(II(i+1,j+1,'x')) = g2(II(i+1,j+1,'x')) + area * ( - x(II(i+1,j  ,'y')) + x(II(i  ,j+1,'y')) )
        g2(II(i+1,j+1,'y')) = g2(II(i+1,j+1,'y')) + area * ( + x(II(i+1,j  ,'x')) - x(II(i  ,j+1,'x')) )
     end do
  end do

  g2(1:n) = g2(1:n) / ( ( nabs - 1 ) * ( nord - 1 ) )

  g(1:n) = sigmafS * g1(1:n) + sigmafA * g2(1:n)

  if ( n .gt. lim ) then
     lmem = .true.
     return
  end if

  if ( m .eq. 0 ) then
     jcnnz = 0
  else
     jcnnz = n
     jcfun(1:n) = 1
     jcvar(1:n) = (/ (i, i = 1, n) /)
     jcval(1:n) = sigmacS * g1(1:n) + sigmacA * g2(1:n)
  end if

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

  use modgrid

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

  ! FUNCTIONS
  integer :: II

  ! LOCAL SCALARS
  integer :: h1nnz,i,j
  real(kind=8) :: area,dareadijx,dareadijy,dareadijpx,dareadijpy, &
       dareadipjx,dareadipjy,dareadipjpx,dareadipjpy

  flag = 0
  lmem = .false.

  if ( 6 * ( 2 * nabs * nord - nabs - nord ) &
       + 36 * ( nabs - 1 ) * ( nord - 1 ) .gt. lim ) then
     lmem = .true.
     return
  end if

  hlnnz = 0

  do i = 1,nabs 
     do j = 1,nord
        if ( j .lt. nord ) then
           ! Distance from P(i,j) to P(i,j+1)

           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j,'x')
           hlcol(hlnnz) = II(i,j,'x')
           hlval(hlnnz) = 1.0d0
 
           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j+1,'x')
           hlcol(hlnnz) = II(i,j,'x')
           hlval(hlnnz) = - 1.0d0
 
           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j+1,'x')
           hlcol(hlnnz) = II(i,j+1,'x')
           hlval(hlnnz) = 1.0d0
 
           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j,'y')
           hlcol(hlnnz) = II(i,j,'y')
           hlval(hlnnz) = 1.0d0
 
           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j+1,'y')
           hlcol(hlnnz) = II(i,j,'y')
           hlval(hlnnz) = - 1.0d0
 
           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j+1,'y')
           hlcol(hlnnz) = II(i,j+1,'y')
           hlval(hlnnz) = 1.0d0
        end if
        if ( i .lt. nabs ) then
           ! Distance from P(i,j) to P(i+1,j)

           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j,'x')
           hlcol(hlnnz) = II(i,j,'x')
           hlval(hlnnz) = 1.0d0

           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i+1,j,'x')
           hlcol(hlnnz) = II(i,j,'x')
           hlval(hlnnz) = - 1.0d0

           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i+1,j,'x')
           hlcol(hlnnz) = II(i+1,j,'x')
           hlval(hlnnz) = 1.0d0

           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i,j,'y')
           hlcol(hlnnz) = II(i,j,'y')
           hlval(hlnnz) = 1.0d0

           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i+1,j,'y')
           hlcol(hlnnz) = II(i,j,'y')
           hlval(hlnnz) = - 1.0d0

           hlnnz = hlnnz + 1
           hlrow(hlnnz) = II(i+1,j,'y')
           hlcol(hlnnz) = II(i+1,j,'y')
           hlval(hlnnz) = 1.0d0
        end if
     end do
  end do

  hlval(1:hlnnz) = hlval(1:hlnnz) / ( ( nabs - 1 ) * nord + ( nord - 1 ) * nabs )

  h1nnz = hlnnz

  do i = 1,nabs - 1
     do j = 1,nord - 1

        area = x(II(i  ,j  ,'x')) * x(II(i+1,j  ,'y')) -  x(II(i+1,j  ,'x')) * x(II(i  ,j  ,'y')) &
             + x(II(i+1,j  ,'x')) * x(II(i+1,j+1,'y')) -  x(II(i+1,j+1,'x')) * x(II(i+1,j  ,'y')) &
             + x(II(i+1,j+1,'x')) * x(II(i  ,j+1,'y')) -  x(II(i  ,j+1,'x')) * x(II(i+1,j+1,'y')) &
             + x(II(i  ,j+1,'x')) * x(II(i  ,j  ,'y')) -  x(II(i  ,j  ,'x')) * x(II(i  ,j+1,'y'))

        dareadijx   = + x(II(i+1,j  ,'y')) - x(II(i  ,j+1,'y'))
        dareadijy   = - x(II(i+1,j  ,'x')) + x(II(i  ,j+1,'x'))
        dareadijpx  = - x(II(i+1,j+1,'y')) + x(II(i  ,j  ,'y'))
        dareadijpy  = + x(II(i+1,j+1,'x')) - x(II(i  ,j  ,'x'))
        dareadipjx  = - x(II(i  ,j  ,'y')) + x(II(i+1,j+1,'y'))
        dareadipjy  = + x(II(i  ,j  ,'x')) - x(II(i+1,j+1,'x'))
        dareadipjpx = - x(II(i+1,j  ,'y')) + x(II(i  ,j+1,'y'))
        dareadipjpy = + x(II(i+1,j  ,'x')) - x(II(i  ,j+1,'x'))

!!$     g2(II(i  ,j  ,'x')) = g2(II(i  ,j  ,'x')) + area * dareadijx  

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j  ,'x')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx * dareadijx

!!$     g2(II(i  ,j  ,'y')) = g2(II(i  ,j  ,'y')) + area * dareadijy  

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j  ,'y')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx * dareadijy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j  ,'y')
        hlcol(hlnnz) = II(i  ,j  ,'y')
        hlval(hlnnz) = dareadijy * dareadijy

!!$     g2(II(i  ,j+1,'x')) = g2(II(i  ,j+1,'x')) + area * dareadijpx 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j+1,'x')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx * dareadijpx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j+1,'x')
        hlcol(hlnnz) = II(i  ,j  ,'y')
        hlval(hlnnz) = dareadijy * dareadijpx + area

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j+1,'x')
        hlcol(hlnnz) = II(i  ,j+1,'x')
        hlval(hlnnz) = dareadijpx * dareadijpx

!!$     g2(II(i  ,j+1,'y')) = g2(II(i  ,j+1,'y')) + area * dareadijpy 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j+1,'y')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx  * dareadijpy - area 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j+1,'y')
        hlcol(hlnnz) = II(i  ,j  ,'y')
        hlval(hlnnz) = dareadijy  * dareadijpy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j+1,'y')
        hlcol(hlnnz) = II(i  ,j+1,'x')
        hlval(hlnnz) = dareadijpx  * dareadijpy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i  ,j+1,'y')
        hlcol(hlnnz) = II(i  ,j+1,'y')
        hlval(hlnnz) = dareadijpy  * dareadijpy

!!$     g2(II(i+1,j  ,'x')) = g2(II(i+1,j  ,'x')) + area * dareadipjx 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'x')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx * dareadipjx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'x')
        hlcol(hlnnz) = II(i  ,j  ,'y')
        hlval(hlnnz) = dareadijy * dareadipjx - area

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'x')
        hlcol(hlnnz) = II(i  ,j+1,'x')
        hlval(hlnnz) = dareadijpx * dareadipjx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'x')
        hlcol(hlnnz) = II(i  ,j+1,'y')
        hlval(hlnnz) = dareadijpy * dareadipjx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'x')
        hlcol(hlnnz) = II(i+1,j  ,'x')
        hlval(hlnnz) = dareadipjx * dareadipjx

!!$     g2(II(i+1,j  ,'y')) = g2(II(i+1,j  ,'y')) + area * dareadipjy 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'y')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx * dareadipjy + area

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'y')
        hlcol(hlnnz) = II(i  ,j  ,'y')
        hlval(hlnnz) = dareadijy * dareadipjy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'y')
        hlcol(hlnnz) = II(i  ,j+1,'x')
        hlval(hlnnz) = dareadijpx * dareadipjy 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'y')
        hlcol(hlnnz) = II(i  ,j+1,'y')
        hlval(hlnnz) = dareadijpy * dareadipjy 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'y')
        hlcol(hlnnz) = II(i+1,j  ,'x')
        hlval(hlnnz) = dareadipjx * dareadipjy 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j  ,'y')
        hlcol(hlnnz) = II(i+1,j  ,'y')
        hlval(hlnnz) = dareadipjy * dareadipjy

!!$     g2(II(i+1,j+1,'x')) = g2(II(i+1,j+1,'x')) + area * dareadipjpx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'x')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx * dareadipjpx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'x')
        hlcol(hlnnz) = II(i  ,j  ,'y')
        hlval(hlnnz) = dareadijy * dareadipjpx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'x')
        hlcol(hlnnz) = II(i  ,j+1,'x')
        hlval(hlnnz) = dareadijpx * dareadipjpx 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'x')
        hlcol(hlnnz) = II(i  ,j+1,'y')
        hlval(hlnnz) = dareadijpy * dareadipjpx + area

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'x')
        hlcol(hlnnz) = II(i+1,j  ,'x')
        hlval(hlnnz) = dareadipjx * dareadipjpx

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'x')
        hlcol(hlnnz) = II(i+1,j  ,'y')
        hlval(hlnnz) = dareadipjy * dareadipjpx - area

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'x')
        hlcol(hlnnz) = II(i+1,j+1,'x')
        hlval(hlnnz) = dareadipjpx * dareadipjpx

!!$     g2(II(i+1,j+1,'y')) = g2(II(i+1,j+1,'y')) + area * dareadipjpy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i  ,j  ,'x')
        hlval(hlnnz) = dareadijx * dareadipjpy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i  ,j  ,'y')
        hlval(hlnnz) = dareadijy * dareadipjpy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i  ,j+1,'x')
        hlval(hlnnz) = dareadijpx * dareadipjpy - area

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i  ,j+1,'y')
        hlval(hlnnz) = dareadijpy * dareadipjpy 

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i+1,j  ,'x')
        hlval(hlnnz) = dareadipjx * dareadipjpy + area

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i+1,j  ,'y')
        hlval(hlnnz) = dareadipjy * dareadipjpy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i+1,j+1,'x')
        hlval(hlnnz) = dareadipjpx * dareadipjpy

        hlnnz = hlnnz + 1
        hlrow(hlnnz) = II(i+1,j+1,'y')
        hlcol(hlnnz) = II(i+1,j+1,'y')
        hlval(hlnnz) = dareadipjpy * dareadipjpy
     end do
  end do

  hlval(h1nnz+1:hlnnz) = hlval(h1nnz+1:hlnnz) / ( ( nabs - 1 ) * ( nord - 1 ) )

  if ( m .eq. 0 ) then
     hlval(1:h1nnz)       = sigmafS * sf * hlval(1:h1nnz)
     hlval(h1nnz+1:hlnnz) = sigmafA * sf * hlval(h1nnz+1:hlnnz)
  else
     hlval(1:h1nnz)       = ( sf * sigmafS + sc(1) * lambda(1) * sigmacS ) * hlval(1:h1nnz)
     hlval(h1nnz+1:hlnnz) = ( sf * sigmafA + sc(1) * lambda(1) * sigmacA ) * hlval(h1nnz+1:hlnnz)
  end if

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

integer function II(i,j,c)

  use modgrid

  ! SCALAR ARGUMENTS
  character, intent(in) :: c
  integer, intent(in) :: i,j

  if ( .not. ( 1 .le. i .and. i .le. nabs ) .or. &
       .not. ( 1 .le. j .and. j .le. nord ) .or. &
       .not. ( c .eq. 'x' .or. c .eq. 'y' ) ) then
     write(*,*) 'Wrong call to indexing function II. Check your code!'
     write(*,*) 'i = ',i,' j = ',j,' c = ',c
     stop
  end if

  if ( c .eq. 'x' ) then
     II = 2 * ( ( i - 1 ) * nord + j ) - 1
  else ! if ( c .eq. 'y' ) then
     II = 2 * ( ( i - 1 ) * nord + j )
  end if

  return

end function II

! ******************************************************************
! ******************************************************************

subroutine drawsol(n,x)

  use modgrid

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! FUNCTIONS
  integer :: II

  ! LOCAL SCALARS
  integer :: i,j
  real(kind=8) :: scale,maxx,maxy,minx,miny,xx1,xx2,yy1,yy2
  logical :: inside1,inside2

  open(unit=10,file='grid.mp')

  ! SCALING
  do i = 1,nabs
     do j = 1,nord
        minx = min( minx, x(II(i,j,'x')) )
        maxx = max( maxx, x(II(i,j,'x')) )

        miny = min( miny, x(II(i,j,'y')) )
        maxy = max( maxy, x(II(i,j,'y')) )
     end do
  end do

  scale = 4.85d0 / ( max( maxx - minx, maxy - miny ) )
  write(10,10) scale

  ! SEGMENTS
  do i = 1,nabs
     do j = 1,nord
        if ( i .lt. nabs ) then
!!$           xx1 = x(II(i,j,'x'))
!!$           yy1 = x(II(i,j,'y'))
!!$           xx2 = x(II(i+1,j,'x'))
!!$           yy2 = x(II(i+1,j,'y'))
!!$           inside1 = .false.
!!$           if ( xx1 .ge. -0.4d0 .and. xx1 .le. 1.4d0 .and. &
!!$                yy1 .ge.  0.4d0 .and. yy1 .le. 2.0d0 ) then
!!$              inside1 = .true.
!!$           end if
!!$           inside2 = .false.
!!$           if ( xx2 .ge. -0.4d0 .and. xx2 .le. 1.4d0 .and. &
!!$                yy2 .ge.  0.4d0 .and. yy2 .le. 2.0d0 ) then
!!$              inside2 = .true.
!!$           end if
!!$           if ( inside1 .or. inside2 ) then
!!$              write(10,30) xx1,yy1,xx2,yy2
!!$           end if
           write(10,30) x(II(i,j,'x')),x(II(i,j,'y')), &
                        x(II(i+1,j,'x')),x(II(i+1,j,'y'))
        end if
        if ( j .lt. nord ) then
!!$           xx1 = x(II(i,j,'x'))
!!$           yy1 = x(II(i,j,'y'))
!!$           xx2 = x(II(i,j+1,'x'))
!!$           yy2 = x(II(i,j+1,'y'))
!!$           inside1 = .false.
!!$           if ( xx1 .ge. -0.4d0 .and. xx1 .le. 1.4d0 .and. &
!!$                yy1 .ge.  0.4d0 .and. yy1 .le. 2.0d0 ) then
!!$              inside1 = .true.
!!$           end if
!!$           inside2 = .false.
!!$           if ( xx2 .ge. -0.4d0 .and. xx2 .le. 1.4d0 .and. &
!!$                yy2 .ge.  0.4d0 .and. yy2 .le. 2.0d0 ) then
!!$              inside2 = .true.
!!$           end if
!!$           if ( inside1 .or. inside2 ) then
!!$              write(10,30) xx1,yy1,xx2,yy2
!!$           end if
           write(10,30) x(II(i,j,'x')),x(II(i,j,'y')), &
                        x(II(i,j+1,'x')),x(II(i,j+1,'y'))
        end if
     end do
  end do

  ! DOTS
  do i = 1,nabs
     do j = 1,nord
!!$        xx1 = x(II(i,j,'x'))
!!$        yy1 = x(II(i,j,'y'))
!!$        if ( xx1 .ge. -0.4d0 .and. xx1 .le. 1.4d0 .and. &
!!$             yy1 .ge.  0.4d0 .and. yy1 .le. 2.0d0 ) then
!!$           write(10,20) xx1,yy1
!!$        end if
      write(10,20) x(II(i,j,'x')),x(II(i,j,'y'))
      ! write(10,21) i,j,x(II(i,j,'x')),x(II(i,j,'y'))
     end do
  end do

  write(10,40)

  close(10)

  ! NON-EXECUTABLE STATEMENTS

10 format('beginfig(1);'/,'u = ',f20.10,' cm;') 
20 format('drawdot (',f20.10,'u,',f20.10,'u) withpen pencircle scaled 1.25;')
21 format('dotlabel ( btex $(',i2,',',i2,')$ etex, (',f20.10,'u,',f20.10,'u)) withpen pencircle scaled 2.5;')
30 format('draw (',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.5 withcolor .7white;')
40 format('endfig;',/,'end;') 

end subroutine drawsol
