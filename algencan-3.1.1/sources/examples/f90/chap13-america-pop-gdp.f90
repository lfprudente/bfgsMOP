! ******************************************************************
! ******************************************************************

program algencanma

  use modamerica

  implicit none

  ! LOCAL SCALARS
  logical :: checkder
  integer :: allocerr,hnnzmax,hnnzmax1,hnnzmax2,hnnzmax3,i,i1,i2, &
       inform,ip,j,jcnnzmax,k,m,n,nvparam,totnbord
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,snorm,totcurval,tottarval

  ! LOCAL ARRAYS
  character(len=80) :: specfnm,outputfnm,vparam(10)
  logical :: coded(11)
  logical, pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,myevalfc, &
       myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! Compute 'centers' of coutries

  do j = 1,nstates
     w(1,j) = 0.0d0
     w(2,j) = 0.0d0
     do i = 1,nbord(j)
        k = pbord(i,j)
        w(1,j) = w(1,j) + put(1,k)
        w(2,j) = w(2,j) + put(2,k)
     end do
     w(1,j) = w(1,j) / nbord(j)
     w(2,j) = w(2,j) / nbord(j)
  end do

  ! Scale target values

  totcurval = 0.0d0
  do j = 1,nstates
     do i = 1,nbord(j)
        ip = mod(i,nbord(j)) + 1
        i1 = pbord(i ,j)
        i2 = pbord(ip,j)
        totcurval = totcurval + &
             0.5d0 * ( put(2,i2) * put(1,i1) - put(2,i1) * put(1,i2) )
     end do
  end do

  tottarval = sum( tarval(1:nstates) )

  tarval(1:nstates) = tarval(1:nstates) * ( totcurval / tottarval )

  ! Number of variables

  n = 2 * np + 5 * nstates

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  do i = 1,np
     l(2*i-1) = - 1.0d+20
     l(2*i)   = - 1.0d+20

     u(2*i-1) =   1.0d+20
     u(2*i)   =   1.0d+20

     x(2*i-1) = put(1,i)
     x(2*i)   = put(2,i)
  end do

  do j = 1,nstates
     l(2*np+5*j-4) = - 0.1d0 * abs( w(1,j) ) 
     l(2*np+5*j-3) = - 0.1d0 * abs( w(2,j) )
     l(2*np+5*j-2) =   0.0d0 
     l(2*np+5*j-1) =   0.0d0
     l(2*np+5*j)   =  -0.1d0

     u(2*np+5*j-4) =   0.1d0 * abs( w(1,j) ) 
     u(2*np+5*j-3) =   0.1d0 * abs( w(2,j) )
     u(2*np+5*j-2) =   2.0d0 
     u(2*np+5*j-1) =   2.0d0
     u(2*np+5*j)   =   0.1d0

     x(2*np+5*j-4) =   0.0d0 
     x(2*np+5*j-3) =   0.0d0
     x(2*np+5*j-2) =   1.0d0 
     x(2*np+5*j-1) =   1.0d0
     x(2*np+5*j)   =   0.0d0
  end do

  ! Constraints

  m = nstates

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  equatn(1:m) = .true.
  linear(1:m) = .false.
  lambda(1:m) =  0.0d0

  ! Coded subroutines

  coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
  coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp

  ! Upper bounds on the number of sparse-matrices non-null elements

  totnbord = sum( nbord(1:nstates) )

  jcnnzmax = 4 * totnbord

  hnnzmax1 = 24 * totnbord
  hnnzmax2 = 2 * totnbord
  hnnzmax3 = 0
  do j = 1,nstates
     hnnzmax3 = hnnzmax3 + 4 * nbord(j) * ( 4 * nbord(j) + 1 ) / 2
  end do

  hnnzmax  = hnnzmax1 + hnnzmax2 + hnnzmax3

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

  nvparam = 0

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

  use modamerica

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag
  real(kind=8), intent(out)  :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! LOCAL SCALARS
  integer :: i,j,k
  real(kind=8) :: ctheta,stheta,theta

  ! LOCAL ARRAYS
  real(kind=8) :: d(2),diff(2),putmw(2),t(2),v(2)

  flag = 0

  f = 0.0d0
  do j = 1,nstates
     t(1)  = x(2*np+5*j-4)
     t(2)  = x(2*np+5*j-3)
     d(1)  = x(2*np+5*j-2)
     d(2)  = x(2*np+5*j-1)
     theta = x(2*np+5*j)
     ctheta = cos(theta)
     stheta = sin(theta)
     do i = 1,nbord(j)
        k = pbord(i,j)

        putmw(1:2) = put(:,k) - w(:,j)

        v(1) = t(1) + w(1,j) + ctheta * d(1) * putmw(1) - stheta * d(2) * putmw(2)
        v(2) = t(2) + w(2,j) + stheta * d(1) * putmw(1) + ctheta * d(2) * putmw(2)

        diff(1:2)  = x(2*k-1:2*k) - v(1:2)

        f = f + 0.5d0 * ( diff(1) ** 2 + diff(2) ** 2 )
     end do
  end do

end subroutine myevalf

! ******************************************************************
! ******************************************************************

subroutine myevalg(n,x,g,flag)

  use modamerica

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: g(n)

  ! LOCAL SCALARS
  integer :: i,j,k
  real(kind=8) :: ctheta,stheta,theta

  ! LOCAL ARRAYS
  real(kind=8) :: d(2),diff(2),putmw(2),t(2),v(2)

  flag = 0

  g(1:n) = 0.0d0

  do j = 1,nstates
     t(1)  = x(2*np+5*j-4)
     t(2)  = x(2*np+5*j-3)
     d(1)  = x(2*np+5*j-2)
     d(2)  = x(2*np+5*j-1)
     theta = x(2*np+5*j)
     ctheta = cos(theta)
     stheta = sin(theta)
     do i = 1,nbord(j)
        k = pbord(i,j)

        putmw(1:2) = put(:,k) - w(:,j)

        v(1) = t(1) + w(1,j) + ctheta * d(1) * putmw(1) - stheta * d(2) * putmw(2)
        v(2) = t(2) + w(2,j) + stheta * d(1) * putmw(1) + ctheta * d(2) * putmw(2)

        diff(1:2)  = x(2*k-1:2*k) - v(1:2)

        g(2*k-1) = g(2*k-1) + diff(1)
        g(2*k)   = g(2*k)   + diff(2)

        g(2*np+5*j-4) = g(2*np+5*j-4) - diff(1)
        g(2*np+5*j-3) = g(2*np+5*j-3) - diff(2)
        g(2*np+5*j-2) = g(2*np+5*j-2) - ( diff(1) * ctheta + diff(2) * stheta ) * putmw(1) 
        g(2*np+5*j-1) = g(2*np+5*j-1) + ( diff(1) * stheta - diff(2) * ctheta ) * putmw(2) 
        g(2*np+5*j)   = g(2*np+5*j) - diff(1) * ( - stheta * d(1) * putmw(1) - ctheta * d(2) * putmw(2) ) &
                                    - diff(2) * (   ctheta * d(1) * putmw(1) - stheta * d(2) * putmw(2) )
      end do
  end do

end subroutine myevalg

! ******************************************************************
! ******************************************************************

subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  use modamerica

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,n
  integer, intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  ! LOCAL SCALARS
  integer :: i,j,k
  real(kind=8) :: ctheta,stheta,theta, &
       dv1dd1,dv1dd2,dv2dd1,dv2dd2,dv1dtheta,dv2dtheta

  ! LOCAL ARRAYS
  real(kind=8) :: d(2),diff(2),putmw(2),t(2),v(2)

  flag = 0
  lmem = .false.

  hnnz = 0

  do j = 1,nstates
     t(1)  = x(2*np+5*j-4)
     t(2)  = x(2*np+5*j-3)
     d(1)  = x(2*np+5*j-2)
     d(2)  = x(2*np+5*j-1)
     theta = x(2*np+5*j)
     ctheta = cos(theta)
     stheta = sin(theta)
     do i = 1,nbord(j)
        k = pbord(i,j)

        putmw(1:2) = put(:,k) - w(:,j)

        v(1) = t(1) + w(1,j) + ctheta * d(1) * putmw(1) - stheta * d(2) * putmw(2)
        v(2) = t(2) + w(2,j) + stheta * d(1) * putmw(1) + ctheta * d(2) * putmw(2)

        diff(1:2)  = x(2*k-1:2*k) - v(1:2)

        dv1dd1 =   ctheta * putmw(1)
        dv1dd2 = - stheta * putmw(2)

        dv2dd1 =   stheta * putmw(1)
        dv2dd2 =   ctheta * putmw(2)

        dv1dtheta = - stheta * d(1) * putmw(1) - ctheta * d(2) * putmw(2)
        dv2dtheta =   ctheta * d(1) * putmw(1) - stheta * d(2) * putmw(2)

        if ( hnnz + 24 .gt. lim ) then
           lmem = .true.
           return
        end if

        hnnz = hnnz + 1
        hrow(hnnz) = 2*k-1
        hcol(hnnz) = 2*k-1
        hval(hnnz) = 1.0d0

        hnnz = hnnz + 1
        hrow(hnnz) = 2*k
        hcol(hnnz) = 2*k
        hval(hnnz) = 1.0d0

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-4
        hcol(hnnz) = 2*k-1
        hval(hnnz) = - 1.0d0

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-4
        hcol(hnnz) = 2*np+5*j-4
        hval(hnnz) = 1.0d0

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-3
        hcol(hnnz) = 2*k
        hval(hnnz) = - 1.0d0

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-3
        hcol(hnnz) = 2*np+5*j-3
        hval(hnnz) = 1.0d0

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-2
        hcol(hnnz) = 2*k-1
        hval(hnnz) = - dv1dd1

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-2
        hcol(hnnz) = 2*k
        hval(hnnz) = - dv2dd1

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-2
        hcol(hnnz) = 2*np+5*j-4
        hval(hnnz) = dv1dd1 

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-2
        hcol(hnnz) = 2*np+5*j-3
        hval(hnnz) = dv2dd1

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-2
        hcol(hnnz) = 2*np+5*j-2
        hval(hnnz) = dv1dd1 ** 2 + dv1dd1 * dv2dd1

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-1
        hcol(hnnz) = 2*k-1
        hval(hnnz) = - dv1dd2

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-1
        hcol(hnnz) = 2*k
        hval(hnnz) = - dv2dd2

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-1
        hcol(hnnz) = 2*np+5*j-4
        hval(hnnz) = dv1dd2

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-1
        hcol(hnnz) = 2*np+5*j-3
        hval(hnnz) =  dv2dd2

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-1
        hcol(hnnz) = 2*np+5*j-2
        hval(hnnz) = dv1dd2

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j-1
        hcol(hnnz) = 2*np+5*j-1
        hval(hnnz) = dv2dd2 ** 2 + dv1dd2 ** 2

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j
        hcol(hnnz) = 2*k-1
        hval(hnnz) = - dv1dtheta

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j
        hcol(hnnz) = 2*k
        hval(hnnz) = - dv2dtheta

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j
        hcol(hnnz) = 2*np+5*j-4
        hval(hnnz) = dv1dtheta

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j
        hcol(hnnz) = 2*np+5*j-3
        hval(hnnz) = dv2dtheta

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j
        hcol(hnnz) = 2*np+5*j-2
        hval(hnnz) = ( dv1dtheta - diff(2) ) * dv1dd1 + ( dv2dtheta + diff(1) ) * dv2dd1

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j
        hcol(hnnz) = 2*np+5*j-1
        hval(hnnz) = ( dv1dtheta - diff(2) ) * dv1dd2 + ( dv2dtheta + diff(1) ) * dv2dd2  

        hnnz = hnnz + 1
        hrow(hnnz) = 2*np+5*j
        hcol(hnnz) = 2*np+5*j
        hval(hnnz) = ( dv1dtheta - diff(2) ) * dv1dtheta + ( dv2dtheta + diff(1) ) * dv2dtheta
      end do
  end do

end subroutine myevalh

! ******************************************************************
! ******************************************************************

subroutine myevalc(n,x,ind,c,flag)

  use modamerica

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: ind,n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: c

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! LOCAL SCALARS
  integer :: i,i1,i2,ip
  real(kind=8) :: area

  flag = 0

  area = 0.0d0

  do i = 1,nbord(ind)
     ip = mod(i,nbord(ind)) + 1
     i1 = pbord(i ,ind)
     i2 = pbord(ip,ind)
     area = area + x(2*i2) * x(2*i1-1) - x(2*i1) * x(2*i2-1)
  end do

  c = 0.5d0 * area / tarval(ind) - 1.0d0

end subroutine myevalc

! ******************************************************************
! ******************************************************************

subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  use modamerica

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: jcvar(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: jcval(lim)

  ! LOCAL SCALARS
  integer :: i,i1,i2,ip

  flag = 0
  lmem = .false.

  jcnnz = 4 * nbord(ind)

  if ( jcnnz .gt. lim ) then
     lmem = .true.
     return
  end if

  jcnnz = 0

  do i = 1,nbord(ind)
     ip = mod(i,nbord(ind)) + 1
     i1 = pbord(i ,ind)
     i2 = pbord(ip,ind)

     jcvar(jcnnz+1) = 2 * i2
     jcval(jcnnz+1) = 0.5d0 * x(2*i1-1) / tarval(ind)

     jcvar(jcnnz+2) = 2 * i1 - 1
     jcval(jcnnz+2) = 0.5d0 * x(2*i2) / tarval(ind)

     jcvar(jcnnz+3) = 2 * i1
     jcval(jcnnz+3) = - 0.5d0 * x(2*i2-1) / tarval(ind)

     jcvar(jcnnz+4) = 2 * i2 - 1
     jcval(jcnnz+4) = - 0.5d0 * x(2*i1) / tarval(ind)

     jcnnz = jcnnz + 4
  end do

!!$  jcnnz = 2 * nbord(ind)
!!$
!!$  if ( jcnnz .gt. lim ) then
!!$     lmem = .true.
!!$     return
!!$  end if
!!$
!!$  jcval(1:jcnnz) = 0.0d0
!!$
!!$  do i = 1,nbord(ind)
!!$     i1 = pbord(i ,ind)
!!$     jcvar(2*i-1) = 2 * i1 - 1
!!$     jcvar(2*i)   = 2 * i1
!!$  end do
!!$
!!$  do i = 1,nbord(ind)
!!$     ip = mod(i,nbord(ind)) + 1
!!$     i1 = pbord(i ,ind)
!!$     i2 = pbord(ip,ind)
!!$
!!$     jcval(2*ip)   = jcval(2*ip)   + 0.5d0 * x(2*i1-1) / tarval(ind)
!!$     jcval(2*i-1)  = jcval(2*i-1)  + 0.5d0 * x(2*i2)   / tarval(ind)
!!$     jcval(2*i)    = jcval(2*i)    - 0.5d0 * x(2*i2-1) / tarval(ind)
!!$     jcval(2*ip-1) = jcval(2*ip-1) - 0.5d0 * x(2*i1)   / tarval(ind)
!!$  end do

end subroutine myevaljac

! ******************************************************************
! ******************************************************************

subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  use modamerica

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,hcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hccol(lim),hcrow(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hcval(lim)

  ! LOCAL SCALARS
  integer :: i,i1,i2,ip

  flag = 0
  lmem = .false.

  if ( 2 * nbord(ind) .gt. lim ) then
     lmem = .true.
     return
  end if

  hcnnz = 0

  do i = 1,nbord(ind)
     ip = mod(i,nbord(ind)) + 1
     i1 = pbord(i ,ind)
     i2 = pbord(ip,ind)

     hcrow(hcnnz+1) = max(  2 * i1 - 1, 2 * i2 )
     hccol(hcnnz+1) = min(  2 * i1 - 1, 2 * i2 )
     hcval(hcnnz+1) = 0.5d0 / tarval(ind)

     hcrow(hcnnz+2) = max(  2 * i1, 2 * i2 - 1 )
     hccol(hcnnz+2) = min(  2 * i1, 2 * i2 - 1 )
     hcval(hcnnz+2) = - 0.5d0 / tarval(ind)

     hcnnz = hcnnz + 2
  end do

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

subroutine drawsol(n,x)

  use modamerica

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! LOCAL SCALARS
  integer :: i,j,k

  open(unit=10,file='map.mp')

  write(10,*) ' beginfig(1)'
  write(10,*) ' prologues:=3;'
  write(10,*) ' path p;'
  write(10,*) ' color coral,green,orange,purple,yellow;'
  write(10,*) ' coral  = (240/255,128/255,128/255);'
  write(10,*) ' green  = (162/255,205/255, 90/255);'
  write(10,*) ' orange = (255/255,165/255, 79/255);'
  write(10,*) ' purple = (159/255,121/255,238/255);'
  write(10,*) ' yellow = (255/255,215/255,  0/255);'
  write(10,*) ' u = 1.0cm;'
  write(10,*) ' pickup pencircle scaled 1.0pt;'

  do i = 1,nstates
     write(10,*) ' p := '

     do j = 1,nbord(i)
        k = pbord(j,i)
        write(10,100) x(2*k-1),x(2*k)
     end do

     write(10,200) color(i)
  end do

  write(10,*) ' endfig;'
  write(10,*) ' end'

  close(10)

  ! NON-EXECUTABLE STATEMENTS

100 format(' (',f10.3,'u,',f10.3,'u)--')
200 format(' cycle; fill p withcolor ',A6,'; draw p;')

end subroutine drawsol
