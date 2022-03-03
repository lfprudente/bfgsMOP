! ******************************************************************
! ******************************************************************

program algencanma

  implicit none

  ! LOCAL SCALARS
  logical :: checkder
  integer :: allocerr,hnnzmax,i,inform,j,jcnnzmax,m,n,nvparam,pij,sl0norm
  real(kind=8) :: cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt, &
       f,nlpsupn,seed,sl1norm,snorm

  ! LOCAL ARRAYS
  character(len=80) :: specfnm,outputfnm,vparam(10)
  logical :: coded(11)
  logical, pointer :: equatn(:),linear(:)
  integer, pointer :: signali(:)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)
  real(kind=8), pointer :: recoveredsignal(:),s(:),s0(:),tmp(:), &
       tmp2(:),truesignal(:)

  ! COMMON SCALARS
  integer :: l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! COMMON ARRAYS
  integer, pointer :: P(:,:)
  real(kind=8), pointer :: y(:)

  ! COMMON BLOCKS
  common /pdata/ y,P,l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! FUNCTIONS
  real(kind=8) :: drand

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,myevalfc, &
       myevalgjac,myevalgjacp,myevalhl,myevalhlp

  ! Problem data

  ! Read original signal from file

  open(20,file='phantom.dat')

  read(20,*) nsignrows,nsigncols

  nsign = nsignrows * nsigncols

  l2nsign = 0
10 if ( 2 ** l2nsign .lt. nsign ) then
     l2nsign = l2nsign + 1
     go to 10
  end if

  if ( 2 ** l2nsign .ne. nsign ) then
     write(*,*)  'ERROR: nsign must be a power of two.'
     stop
  end if

  allocate(truesignal(nsign),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error'
     stop
  end if

  do i = 1,nsignrows
     read(20,*) (truesignal((i-1)*nsigncols+j),j=1,nsigncols)
  end do

  close(20)

  ! Compute known solution s for comparison purposes only (s represents
  ! the coefficients of the true signal as a linear combination of the
  ! rows of Psi)

  allocate(s(nsign),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error'
     stop
  end if

  call multpsi(0,truesignal,s)

  sl0norm = 0
  sl1norm = 0.0d0
  do i = 1,nsign
     if ( s(i) .ne. 0.0d0 ) then
        sl0norm = sl0norm + 1
     end if
     sl1norm = sl1norm + abs( s(i) )
  end do

  write(*,*) 'Known solution l0-norm =',dble(sl0norm) / nsign
  write(*,*) 'Known solution l1-norm =',sl1norm

  ! Set number of samples (msampl < nsign)

  ! msampl = floor( 1.5d0 * sl0norm )
  ! msampl = floor( 2.0d0 * sl0norm )
  ! msampl = floor( 2.5d0 * sl0norm )
  ! msampl = floor( 3.0d0 * sl0norm )
  msampl = floor( 3.5d0 * sl0norm )
  ! msampl = floor( 4.0d0 * sl0norm )
  ! msampl = floor( 4.5d0 * sl0norm )
  ! msampl = floor( 5.0d0 * sl0norm )

  ! Set parameter L for the sampling matrix Phi (lPhi < msampl)

  lPhi = 2

  ! Generate permutation matrices for sampling matrix Phi

  allocate(P(nsign,lPhi),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error'
     stop
  end if

  seed = 123456.0d0

  do j = 1,lPhi
     do i = 1,nsign
        pij    = floor( ( nsign - i + 1 ) * drand(seed) )
        P(i,j) = i + min( pij, nsign - i )
     end do
  end do

  ! Sample y from the true signal

  allocate(y(msampl),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error'
     stop
  end if

  call multphi(0,truesignal,y)

  ! Show problem information

  write(*,*) 'Problem generated using signal file: phantom.dat'
  write(*,*) 'nsignrows: ',nsignrows
  write(*,*) 'nsigncols: ',nsigncols
  write(*,*) 'nsign    : ',nsign
  write(*,*) 'msampl   : ',msampl
  write(*,*) 'lPhi     : ',lPhi

  ! Number of variables

  n = 2 * nsign

  ! Set lower bounds, upper bounds, and initial guess

  allocate(x(n),l(n),u(n))

  l(1:n) = 0.0d0
  u(1:n) = 1.0d+20

  ! Initial guess s0 = Theta^T ( Theta Theta^T )^{-1} y, where Theta = Phi Psi^T.
  ! tmp is computed by CG as the solution of ( Theta Theta^T ) tmp = y.
  ! (Since Psi is an orthogonal matrix, Theta Theta^T = Phi Phi^T.)
  ! Then s0 is computed as s0 = Theta^T tmp.

  allocate(s0(nsign),tmp(msampl),tmp2(nsign),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error'
     stop
  end if

  call cg(msampl,tmp,y)
  call multphi(1,tmp,tmp2)
  call multpsi(0,tmp2,s0)

  do i = 1,nsign
     if ( s0(i) .ge. 0.0d0 ) then
        x(i)         =   s0(i)
        x(nsign + i) =   0.0d0
     else
        x(i)         =   0.0d0
        x(nsign + i) = - s0(i)
     end if
  end do

  ! Constraints

  m = msampl

  allocate(equatn(m),linear(m),lambda(m))

  equatn(1:m) = .true.
  linear(1:m) = .true.
  lambda(1:m) =  0.0d0

  ! Coded subroutines

  coded(1:11) = .false.
  coded(7)    = .true.  ! evalfc
  coded(9)    = .true.  ! evalgjacp
  coded(11)   = .true.  ! evalhlp

  ! Upper bounds on the number of sparse-matrices non-null elements

  jcnnzmax = 0
  hnnzmax  = 0

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

  nvparam   = 1
  vparam(1) = 'OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED'

  ! Optimize

  call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,  &
       myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax,  &
       hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,  &
       checkder,f,cnorm,snorm,nlpsupn,inform)

  ! Build recovered signal

  allocate(recoveredsignal(nsign),signali(nsign),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error'
     stop
  end if

  s(1:nsign) = x(1:nsign) - x(nsign+1:2*nsign)

  sl0norm = 0
  sl1norm = 0.0d0
  do i = 1,nsign
     if ( s(i) .ne. 0.0d0 ) then
        sl0norm = sl0norm + 1
     end if
     sl1norm = sl1norm + abs( s(i) )
  end do

  write(*,*) 'Solution l0-norm =',dble(sl0norm) / nsign
  write(*,*) 'Solution l1-norm =',sl1norm

  call multpsi(1,s,recoveredsignal)

  signali(1:nsign) = nint( max( 0.0d0, min( 255.0d0 * recoveredsignal(1:nsign), 255.0d0 ) ) )

  open(10,file='recovered.pgm')

  write(10,'(A/)') 'P2'
  write(10,*) '#'
  write(10,*) nsigncols,nsignrows
  write(10,*) '255'

  do i = 1,nsignrows
     write(10,*) (signali((i-1)*nsigncols+j),j=1,nsigncols)
  end do

  close(10)

  deallocate(x,l,u,lambda,equatn,linear)

  deallocate(P,recoveredsignal,s,s0,signali,tmp,tmp2,truesignal,y)

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

  ! Compute gradient of the objective function

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
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  ! Compute (lower triangle of the) Hessian of the objective function

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

  ! Compute ind-th constraint

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

  ! Compute gradient of the ind-th constraint

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

  ! Compute lower-triangle of the ind-th constraint's Hessian

  flag = - 1

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

  ! COMMON SCALARS
  integer :: l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! COMMON ARRAYS
  integer, pointer :: P(:,:)
  real(kind=8), pointer :: y(:)

  ! COMMON BLOCKS
  common /pdata/ y,P,l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! LOCAL ARRAYS
  real(kind=8) :: tmp(nsign),uv(nsign)

  flag = 0

  f = sum( x(1:n) )

  uv(1:nsign) = x(1:nsign) - x(nsign+1:2*nsign)

  call multpsi(1,uv,tmp)
  call multphi(0,tmp,c)

  c(1:msampl) = c(1:msampl) - y(1:msampl)

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

subroutine myevalgjacp(n,x,g,m,pp,q,work,gotj,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: gotj
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  character, intent(in) :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(inout) :: pp(m),q(n)
  real(kind=8), intent(out) :: g(n)

  ! COMMON SCALARS
  integer :: l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! COMMON ARRAYS
  integer, pointer :: P(:,:)
  real(kind=8), pointer :: y(:)

  ! COMMON BLOCKS
  common /pdata/ y,P,l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! LOCAL ARRAYS
  real(kind=8) :: tmp(nsign),uv(nsign)

  flag = 0

  ! Gradient of the objective function

  if ( work .eq. 'J' .or. work .eq. 'T' ) then
     g(1:n) = 1.0d0
  end if

  if ( work .eq. 'j' .or. work .eq. 'J' ) then
     ! pp = Jacobian x q
     uv(1:nsign) = q(1:nsign) - q(nsign+1:2*nsign)
     call multpsi(1,uv,tmp)
     call multphi(0,tmp,pp)

  else ! if ( work .eq. 't' .or. work .eq. 'T' ) then
     ! q = Jacobian^t x pp
     call multphi(1,pp,tmp)
     call multpsi(0,tmp,q)
     q(nsign+1:2*nsign) = - q(1:nsign)
  end if

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

  flag = 0

  hp(1:n) = 0.0d0

end subroutine myevalhlp

! ******************************************************************
! ******************************************************************

subroutine multpsi(work,z,w)

  implicit none 

  ! SCALAR ARGUMENTS
  integer, intent(in) :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: z(*)
  real(kind=8), intent(out) :: w(*)

  ! Psi in R^{nsign x nsign} is such that the true signal can be 
  ! represented as a 'sparse' combination of its ROWS.

  ! Returns w = Psi z if work=0 and w = Psi^T z if work=1.

  ! COMMON SCALARS
  integer :: l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! COMMON ARRAYS
  integer, pointer :: P(:,:)
  real(kind=8), pointer :: y(:)

  ! COMMON BLOCKS
  common /pdata/ y,P,l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! LOCAL SCALARS
  integer :: i,j
  real(kind=8) :: cte

  ! LOCAL ARRAYS
  real(kind=8) :: tmp(nsign)

  cte = 1.0d0 / sqrt( 2.0d0 )

  w(1:nsign) = z(1:nsign)

  if ( work .eq. 0 ) then 
     do i = l2nsign,1,-1
        do j = 1,2 ** ( i - 1 )
           tmp(j)          = cte * ( w(2 * j - 1) + w(2 * j) )
           tmp(2**(i-1)+j) = cte * ( w(2 * j - 1) - w(2 * j) )
        end do

        w(1:2**i) = tmp(1:2**i)
     end do

  else
     do i = 1,l2nsign 
        do j = 1,2 ** ( i - 1 )
           tmp(2*j-1) = cte * ( w(j) + w(2 ** ( i - 1 ) + j ) )
           tmp(2*j)   = cte * ( w(j) - w(2 ** ( i - 1 ) + j ) )
        end do

        w(1:2**i) = tmp(1:2**i)
     end do
  end if

end subroutine multpsi

! ******************************************************************
! ******************************************************************

subroutine multphi(work,z,w)

  implicit none 

  ! SCALAR ARGUMENTS
  integer :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: z(*)
  real(kind=8), intent(out) :: w(*)

  ! Phi in R^{msampl x nsign} represents the sampling matrix.

  ! Returns w = Phi z if work=0 and w = Phi^T z if work=1.

  ! COMMON SCALARS
  integer :: l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! COMMON ARRAYS
  integer, pointer :: P(:,:)
  real(kind=8), pointer :: y(:)

  ! COMMON BLOCKS
  common /pdata/ y,P,l2nsign,lPhi,msampl,nsign,nsigncols,nsignrows

  ! LOCAL SCALARS
  integer :: i,ibstart,j,jstart,k,nbrows,nones,quo,quob,res,resb

  ! LOCAL ARRAYS
  integer :: pcol(nsign)
  real(kind=8) :: tmp(nsign),tmp2(nsign)

  if ( work .eq. 0 ) then

     quo = msampl / lPhi
     res = msampl - quo * lPhi

     ibstart = 1
     do k = 1,lPhi
        if ( k .le. res ) then
           nbrows = quo + 1
        else
           nbrows = quo
        end if

        pcol(1:nsign) = P(1:nsign,k)

        call perm(0,nsign,pcol,z,tmp)

        quob = nsign / nbrows
        resb = nsign - quob * nbrows

        jstart = 1
        do i = 1,nbrows
           if ( i .le. resb ) then
              nones = quob + 1
           else
              nones = quob
           end if
           w(ibstart + i - 1) = sum( tmp(jstart:jstart+nones-1) ) / sqrt( dble(nones) )
           jstart = jstart + nones
        end do

        ibstart = ibstart + nbrows
     end do

  else
     w(1:nsign) = 0.0d0

     quo = msampl / lPhi
     res = msampl - quo * lPhi

     ibstart = 1
     do k = 1,lPhi
        if ( k .le. res ) then
           nbrows = quo + 1
        else
           nbrows = quo
        end if

        quob = nsign / nbrows
        resb = nsign - quob * nbrows

        jstart = 1
        do i = 1,nbrows
           if ( i .le. resb ) then
              nones = quob + 1
           else
              nones = quob
           end if

           tmp(jstart:jstart+nones-1) = z(ibstart + i - 1) / sqrt( dble(nones) )

           jstart = jstart + nones
        end do

        pcol(1:nsign) = P(1:nsign,k)

        call perm(1,nsign,pcol,tmp,tmp2)

        w(1:nsign) = w(1:nsign) + tmp2(1:nsign)

        ibstart = ibstart + nbrows
     end do
  end if

end subroutine multphi

! ******************************************************************
! ******************************************************************

subroutine perm(work,n,pcol,z,w)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n,work

  ! ARRAY ARGUMENTS
  integer, intent(in) :: pcol(n)
  real(kind=8), intent(in) :: z(n)
  real(kind=8), intent(out) :: w(n)

  ! Returns w = P z if work = 0 and w = P^T z if work = 1.

  ! LOCAL SCALARS
  integer :: i

  w(1:n) = z(1:n)

  if ( work .eq. 0 ) then
     do i = n,1,-1
        call dswap(w(i),w(pcol(i)))
     end do

  else
     do i = 1,n
        call dswap(w(i),w(pcol(i)))
     end do
  end if

end subroutine perm

! ******************************************************************
! ******************************************************************

subroutine dswap(a,b)

  implicit none

  ! SCALAR ARGUMENTS
  real(kind=8), intent(inout) :: a,b

  ! LOCAL SCALARS
  real(kind=8) :: tmp

  tmp = a
  a   = b
  b   = tmp

end subroutine dswap

! ******************************************************************
! ******************************************************************

subroutine cg(n,x,b)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAYS ARGUMENTS
  real(kind=8), intent(in) :: b(n)
  real(kind=8), intent(out) :: x(n)

  ! LOCAL SCALARS
  integer :: iter
  real(kind=8) :: alpha,beta,rnorm2,rnorm2prev

  ! LOCAL ARRAYS
  real(kind=8) :: r(n),Ap(n),p(n),tmp(n),tmp2(n)

  x(1:n) = 0.0d0

  r(1:n) = - b(1:n)

  rnorm2 = sum( r(1:n) ** 2 )

  p(1:n) = - r(1:n)

  iter = 0

  do while ( iter .le. 100 .and. .not. rnorm2 .le. 1.0d-08 .and. &
       ( iter .ne. 0 .and. rnorm2 .le. rnorm2prev ) )

     call multphi(1,p,tmp)
     call multphi(0,tmp,Ap)

     alpha = rnorm2 / sum( p(1:n) * Ap(1:n) )

     x(1:n) = x(1:n) +  alpha * p(1:n)

     r(1:n) = r(1:n) + alpha * Ap(1:n)

     rnorm2prev = rnorm2

     rnorm2 = sum( r(1:n) ** 2 )

     beta = rnorm2 / rnorm2prev

     p(1:n) = - r(1:n) + beta * p(1:n)

     iter = iter + 1

  end do

end subroutine cg

