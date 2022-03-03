! ******************************************************************
! ******************************************************************

program algencanma

  implicit none

  ! LOCAL SCALARS
  logical :: checkder
  integer :: allocerr,hnnzmax,hnnzmax1,hnnzmax2,hnnzmax3,i,inform, &
       jcnnzmax,m,n,nvparam
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,seed,snorm

  ! LOCAL ARRAYS
  character(len=80) :: specfnm,outputfnm,vparam(10)
  logical :: coded(11)
  logical, pointer :: equatn(:),linear(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)

  ! COMMON SCALARS
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

  ! FUNCTIONS
  real(kind=8) :: drand

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,myevalfc, &
       myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! Problem data

  ndim   = 3
  nite   = 1000
  r      = 1.0d0
  objrad = 15.0d0

  ! Number of variables

  n = ndim * nite

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  l(1:n) = - 1.0d+20
  u(1:n) =   1.0d+20

  seed = 123456.0d0
  do i = 1,n
     x(i) = - objrad + 2.0d0 * objrad * drand(seed)
  end do

  ! Constraints

  m = nite

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if

  equatn(1:m) = .false.
  linear(1:m) = .false.
  lambda(1:m) = 0.0d0

  ! Coded subroutines

  coded(1:6)  = .true.  ! evalf, evalg, evalh, evalc, evaljac, evalhc
  coded(7:11) = .false. ! evalfc,evalgjac,evalgjacp,evalhl,evalhlp

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = nite * ndim
  hnnzmax1 = nite * ndim * ( nite * ndim + 1 ) / 2
  hnnzmax2 = nite * ndim
  hnnzmax3 = nite * ( ndim * ( ndim + 1 ) / 2 )
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

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! COMMON SCALARS
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! LOCAL SCALARS
  integer :: i,ind1,ind2,j

  ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

  ! Compute objective function

  flag = 0

  f = 0.0d0
  do i = 1,nite
     ind1 = ndim * ( i - 1 )
     do j = i + 1,nite
        ind2 = ndim * ( j - 1 )
        f = f + max( 0.0d0, ( 2.0d0 * r ) ** 2 - & 
             sum( ( x(ind1+1:ind1+ndim) - x(ind2+1:ind2+ndim) ) ** 2 ) ) ** 2
     end do
  end do

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

  ! COMMON SCALARS
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! LOCAL SCALARS
  integer :: i,ind1,ind2,j
  real(kind=8) :: fparc

  ! LOCAL ARRAYS
  real(kind=8) :: val(ndim),xdiff(ndim)

 ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

  ! Compute gradient of the objective function

  flag = 0

  g(1:n) = 0.0d0

  do i = 1,nite
     ind1 = ndim * ( i - 1 )
     do j = i + 1,nite
        ind2 = ndim * ( j - 1 )
        xdiff(1:ndim) = x(ind1+1:ind1+ndim) - x(ind2+1:ind2+ndim)
        fparc = ( 2.0d0 * r ) ** 2 - sum( xdiff(1:ndim) ** 2 )
        if ( fparc .gt. 0.0d0 ) then
           val(1:ndim) = 4.0d0 * fparc * xdiff(1:ndim)
           g(ind1+1:ind1+ndim) = g(ind1+1:ind1+ndim) - val(1:ndim)
           g(ind2+1:ind2+ndim) = g(ind2+1:ind2+ndim) + val(1:ndim)
        end if
     end do
  end do

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

  ! COMMON SCALARS
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! LOCAL SCALARS
  integer :: i,ind1,ind2,j,k,l
  real(kind=8) :: fparc,val

  ! LOCAL ARRAYS
  real(kind=8) :: diagblock(nite,ndim,ndim),xdiff(ndim)

  ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

  ! Compute (lower triangle of the) Hessian of the objective function

  flag = 0
  lmem = .false.

  diagblock(1:nite,1:ndim,1:ndim) = 0.0d0

  hnnz = 0
  do i = 1,nite
     ind1 = ndim * ( i - 1 )
     do j = i + 1,nite
        ind2 = ndim * ( j - 1 )
        xdiff(1:ndim) = x(ind1+1:ind1+ndim) - x(ind2+1:ind2+ndim)
        fparc = ( 2.0d0 * r ) ** 2 - sum( xdiff(1:ndim) ** 2 )
        if ( fparc .gt. 0.0d0 ) then
           do k = 1,ndim
              val = 8.0d0 * xdiff(k) ** 2 - 4.0d0 * fparc
              diagblock(i,k,k) = diagblock(i,k,k) + val
              diagblock(j,k,k) = diagblock(j,k,k) + val

              if ( hnnz + 1 .gt. lim ) then
                 lmem = .true.
                 return
              end if

              hrow(hnnz+1) = ind2 + k
              hcol(hnnz+1) = ind1 + k
              hval(hnnz+1) = - val
              hnnz = hnnz + 1

              do l = 1,k - 1
                 val = 8.0d0 * xdiff(k) * xdiff(l)
                 diagblock(i,k,l) = diagblock(i,k,l) + val
                 diagblock(j,k,l) = diagblock(j,k,l) + val

                 if ( hnnz + 1 .gt. lim ) then
                    lmem = .true.
                    return
                 end if

                 hrow(hnnz+1) = ind2 + k
                 hcol(hnnz+1) = ind1 + l
                 hval(hnnz+1) = - val
                 hnnz = hnnz + 1
              end do

              do l = k + 1,ndim
                 val = 8.0d0 * xdiff(k) * xdiff(l)

                 if ( hnnz + 1 .gt. lim ) then
                    lmem = .true.
                    return
                 end if

                 hrow(hnnz+1) = ind2 + k
                 hcol(hnnz+1) = ind1 + l
                 hval(hnnz+1) = - val
                 hnnz = hnnz + 1
              end do
           end do
        end if
     end do
  end do

  do i = 1,nite
     do k = 1,ndim
        do l = 1,k
           if ( hnnz + 1 .gt. lim ) then
              lmem = .true.
              return
           end if

           hrow(hnnz+1) = ndim * ( i - 1 ) + k
           hcol(hnnz+1) = ndim * ( i - 1 ) + l
           hval(hnnz+1) = diagblock(i,k,l)
           hnnz = hnnz + 1
        end do
     end do
  end do

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
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! LOCAL SCALARS
  integer :: ini

  ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

  ! Compute ind-th constraint

  flag = 0

  ini = ndim * ( ind - 1 )
  c = sum( x(ini+1:ini+ndim) ** 2 ) - ( objrad -r ) ** 2

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
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! LOCAL SCALARS
  integer :: i,ini

  ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

 ! Compute gradient of the ind-th constraint

  flag = 0
  lmem = .false.

  jcnnz = ndim

  if ( jcnnz .gt. lim ) then
     lmem = .true.
     return
  end if

  ini = ndim * ( ind - 1 )

  do i = 1,ndim
     jcvar(i) = ini + i
     jcval(i) = 2.0d0 * x(ini+i)
  end do

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
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! LOCAL SCALARS
  integer :: i,ini

  ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

  ! Compute gradient of the ind-th constraint

  flag = 0
  lmem = .false.

  hcnnz = ndim

  if ( hcnnz .gt. lim ) then
     lmem = .true.
     return
  end if

  ini = ndim * ( ind - 1 )

  do i = 1,ndim
     hcrow(i) = ini + i
     hccol(i) = ini + i
     hcval(i) = 2.0d0
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

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! COMMON SCALARS
  integer :: ndim,nite
  real(kind=8) :: r,objrad

  ! LOCAL SCALARS
  integer :: i,j

  ! COMMON BLOCKS
  common /pdata/ ndim,nite,r,objrad

  open(10,file='solution.xyz')

  write(10,9000) nite

  do i = 1,nite
     write(10,9010) (x((i-1)*ndim+j),j=1,ndim)
  end do

  close(10)

  ! NON-EXECUTABLE STATEMENTS

9000 format(I10,/,'Solution')
9010 format('H',3(1X,F20.10))

end subroutine drawsol
