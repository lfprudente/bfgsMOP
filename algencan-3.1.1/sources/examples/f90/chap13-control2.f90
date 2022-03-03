! ******************************************************************
! ******************************************************************

program algencanma

  use modcontrol

  implicit none

  ! LOCAL SCALARS
  logical :: checkder
  integer :: allocerr,hnnzmax,hnnzmax1,hnnzmax2,hnnzmax3,i, &
       inform,j,jcnnzmax,m,n,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,snorm

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

  ! Problem data

  Ndisc  =   1000000
  deltat =   1.0d0 / Ndisc
  x0     = - 2.0d0
  y0     =   4.0d0
  z0     =   0.0d0

  ! Number of variables

  n = 4 * Ndisc

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  l(1:n) = - 1.0d+20
  u(1:n) =   1.0d+20

  x(1:n) =   0.0d0

  ! Constraints

  m = 3 * Ndisc

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  equatn(1:m) = .true.

  do j = 0,Ndisc - 1
     linear(3*j+1) = .true.
     linear(3*j+2) = .false.
     linear(3*j+3) = .false.
  end do

  lambda(1:m) =  0.0d0

  ! Coded subroutines

  coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
  coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = 12 * Ndisc

  hnnzmax1 = 0
  hnnzmax2 = 6 * Ndisc
  hnnzmax3 = 31 * Ndisc
  hnnzmax  = hnnzmax1 + hnnzmax2 + hnnzmax3

  ! Checking derivatives?

  checkder = .false.

  ! Parameters setting

  epsfeas   = 1.0d-08
  epsopt    = 1.0d-08

  efstain   = sqrt( epsfeas )
  eostain   = epsopt ** 1.5d0

  !efacc   = sqrt( epsfeas )
  !eoacc   = sqrt( epsopt )
  efacc     = 1.0d+20
  eoacc     = 1.0d+20

  outputfnm = ''
  specfnm   = ''

  nvparam = 0

  ! Optimize

  call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,  &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax,  &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,  &
       checkder,f,cnorm,snorm,nlpsupn,inform)

  do i = 1,Ndisc
     write(81,100) (1.0d0*(i-1))/(Ndisc-1),x(II('x',i))
     write(82,100) (1.0d0*(i-1))/(Ndisc-1),x(II('y',i))
     write(83,100) (1.0d0*(i-1))/(Ndisc-1),x(II('u',i-1))
  end do

  deallocate(x,l,u,lambda,equatn,linear,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in main program'
     stop
  end if

100 format(1P,D16.8,1X,1P,D16.8)

  stop

end program algencanma

! ******************************************************************
! ******************************************************************

subroutine myevalf(n,x,f,flag)

  use modcontrol

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag
  real(kind=8), intent(out)  :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! FUNCTIONS
  integer :: II

  flag = 0

  f = x(II('z',Ndisc))

end subroutine myevalf

! ******************************************************************
! ******************************************************************

subroutine myevalg(n,x,g,flag)

  use modcontrol

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: g(n)

  ! FUNCTIONS
  integer :: II

  flag = 0

  g(1:n) = 0.0d0

  g(II('z',Ndisc)) = 1.0d0

end subroutine myevalg

! ******************************************************************
! ******************************************************************

subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  use modcontrol

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,n
  integer, intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  flag = 0
  lmem = .false.

  hnnz = 0

end subroutine myevalh

! ******************************************************************
! ******************************************************************

subroutine myevalc(n,x,ind,c,flag)

  use modcontrol

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: ind,n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: c

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! FUNCTIONS
  integer :: II

  ! LOCAL SCALARS
  integer :: i

  flag = 0

  i = ( ind - 1 ) / 3

  if ( mod(ind,3) .eq. 1 ) then
     if ( i .eq. 0 ) then
        c = x(II('x',i+1)) - x0 - deltat * y0
     else
        c = x(II('x',i+1)) - x(II('x',i)) - deltat * x(II('y',i)) 
     end if

  else if ( mod(ind,3) .eq. 2 ) then
     if ( i .eq. 0 ) then
        c = x(II('y',i+1)) - y0 - deltat * &
             ( - x0 - ( x0 ** 2 - 1.0d0 ) * y0 + x(II('u',i)) )
     else
        c = x(II('y',i+1)) - x(II('y',i)) - deltat * &
             ( - x(II('x',i)) - ( x(II('x',i)) ** 2 - 1.0d0 ) * x(II('y',i)) + x(II('u',i)) )
     end if

  else ! if ( mod(ind,3) .eq. 3 ) then
     if ( i .eq. 0 ) then
        c = x(II('z',i+1)) - z0 - deltat * 0.5d0 * ( x0 ** 2 + y0 ** 2 + x(II('u',i)) ** 2 )
     else
        c = x(II('z',i+1)) - x(II('z',i)) - deltat * 0.5d0 * &
             ( x(II('x',i)) ** 2 + x(II('y',i)) ** 2 + x(II('u',i)) ** 2 )
     end if
  end if

end subroutine myevalc

! ******************************************************************
! ******************************************************************

subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  use modcontrol

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: jcvar(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: jcval(lim)

  ! FUNCTIONS
  integer :: II

  ! LOCAL SCALARS
  integer :: i

  flag = 0
  lmem = .false.

  i = ( ind - 1 ) / 3

  if ( mod(ind,3) .eq. 1 ) then

     if ( i .eq. 0 ) then
        jcnnz = 1

        if ( jcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        jcvar(1) = II('x',i+1)
        jcval(1) = 1.0d0
     else
        jcnnz = 3

        if ( jcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        jcvar(1) = II('x',i+1)
        jcval(1) = 1.0d0

        jcvar(2) = II('x',i)
        jcval(2) = - 1.0d0

        jcvar(3) = II('y',i)
        jcval(3) = - deltat
     end if

  else if ( mod(ind,3) .eq. 2 ) then
     if ( i .eq. 0 ) then
        jcnnz = 2

        if ( jcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        jcvar(1) = II('y',i+1)
        jcval(1) = 1.0d0

        jcvar(2) = II('u',i)
        jcval(2) = - deltat
     else
        jcnnz = 4

        if ( jcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        jcvar(1) = II('y',i+1)
        jcval(1) = 1.0d0

        jcvar(2) = II('x',i)
        jcval(2) = - deltat * ( - 1.0d0 - 2.0d0 * x(II('x',i)) * x(II('y',i)) )

        jcvar(3) = II('y',i)
        jcval(3) = - 1.0d0 + deltat * ( x(II('x',i)) ** 2 - 1.0d0 )

        jcvar(4) = II('u',i)
        jcval(4) = - deltat
     end if

  else ! if ( mod(ind,3) .eq. 3 ) then
     if ( i .eq. 0 ) then
        jcnnz = 2

        if ( jcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        jcvar(1) = II('z',i+1)
        jcval(1) = 1.0d0

        jcvar(2) = II('u',i)
        jcval(2) = - deltat * x(II('u',i))
     else
        jcnnz = 5

        if ( jcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        jcvar(1) = II('z',i+1)
        jcval(1) = 1.0d0

        jcvar(2) = II('x',i)
        jcval(2) = - deltat * x(II('x',i))

        jcvar(3) = II('y',i)
        jcval(3) = - deltat * x(II('y',i))

        jcvar(4) = II('z',i)
        jcval(4) = - 1.0d0

        jcvar(5) = II('u',i)
        jcval(5) = - deltat * x(II('u',i))
     end if
  end if

end subroutine myevaljac

! ******************************************************************
! ******************************************************************

subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  use modcontrol

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,hcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hccol(lim),hcrow(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hcval(lim)

  ! FUNCTIONS
  integer :: II

  ! LOCAL SCALARS
  integer :: i

  flag = 0
  lmem = .false.

  i = ( ind - 1 ) / 3

  if ( mod(ind,3) .eq. 1 ) then
     hcnnz = 0

  else if ( mod(ind,3) .eq. 2 ) then
     if ( i .ne. 0 ) then
        hcnnz = 2

        if ( hcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        hcrow(1) = II('x',i)
        hccol(1) = II('x',i)
        hcval(1) = - deltat * 2.0d0 * x(II('y',i))

        hcrow(2) = max( II('x',i), II('y',i) )
        hccol(2) = min( II('x',i), II('y',i) )
        hcval(2) = - deltat * 2.0d0 * x(II('x',i))
     end if

  else ! if ( mod(ind,3) .eq. 3 ) then
     if ( i .eq. 0 ) then
        hcnnz = 1

        if ( hcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        hcrow(1) = II('u',i)
        hccol(1) = II('u',i)
        hcval(1) = - deltat
     else
        hcnnz = 3

        if ( hcnnz .gt. lim ) then
           lmem = .true.
           return
        end if

        hcrow(1) = II('x',i)
        hccol(1) = II('x',i)
        hcval(1) = - deltat

        hcrow(2) = II('y',i)
        hccol(2) = II('y',i)
        hcval(2) = - deltat

        hcrow(3) = II('u',i)
        hccol(3) = II('u',i)
        hcval(3) = - deltat
     end if
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

integer function II(c,i)

  use modcontrol

  character, intent(in) :: c
  integer, intent(in) :: i

  if ( ( c .eq. 'u' .and. .not. ( 0 .le. i .and. i .le. Ndisc - 1 ) ) .or. &
       ( c .eq. 'x' .and. .not. ( 1 .le. i .and. i .le. Ndisc     ) ) .or. &
       ( c .eq. 'y' .and. .not. ( 1 .le. i .and. i .le. Ndisc     ) ) .or. &
       ( c .eq. 'z' .and. .not. ( 1 .le. i .and. i .le. Ndisc     ) ) .or. &
       ( c .ne. 'u' .and. c .ne. 'x' .and. c .ne. 'y' .and. c .ne. 'z' ) ) then
     write(*,*) 'Wrong call to indexing function II. Check your code!'
     write(*,*) 'c = ',c,' i = ',i
     stop
  end if

  if ( c .eq. 'u' ) then
     II = 4 * i + 1
  else if ( c .eq. 'x' ) then
     II = 4 * i - 2
  else if ( c .eq. 'y' ) then
     II = 4 * i - 1
  else ! if ( c .eq. 'z' ) then
     II = 4 * i
  end if

  return

end function II
