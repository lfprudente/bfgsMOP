! ******************************************************************
! ******************************************************************

subroutine calcal(nind,x,m,lambda,rho,equatn,linear,al,inform)

  use modrspace

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: m,nind
  integer,      intent(inout) :: inform
  real(kind=8), intent(out)   :: al

  ! ARRAY ARGUMENTS
  logical,      intent(in) :: equatn(m),linear(m)
  real(kind=8), intent(in) :: lambda(m),rho(m),x(nind)

  ! Set xfull with values in x
  xfull(ind(1:nind)) = x(1:nind)

  ! Compute augmented Lagrangian
  call sevalal(nfull,xfull,m,lambda,rho,equatn,linear,al,inform)
  if ( inform .ne. 0 ) return

end subroutine calcal

! ******************************************************************
! ******************************************************************

subroutine calcnal(nind,x,m,lambda,rho,equatn,linear,nal,inform)

  use modrspace

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,nind
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)  :: equatn(m),linear(m)
  real(kind=8), intent(in)  :: lambda(m),rho(m),x(nind)
  real(kind=8), intent(out) :: nal(nfull)

  ! Set xfull with values in x
  xfull(ind(1:nind)) = x(1:nind)

  ! Compute the gradient of the augmented Lagrangian
  call sevalnal(nfull,xfull,m,lambda,rho,equatn,linear,nal,inform)
  if ( inform .ne. 0 ) return

  ! Shrink nal to the reduced space 
  ! (nal must be the shrinked full-space gradient)
  call shrink(nind,nal)

end subroutine calcnal

! ******************************************************************
! ******************************************************************

subroutine calchal(nind,x,m,lambda,rho,equatn,linear,hrow,hcol, & 
     hval,hnnz,lim,inform)

  use modrspace

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: lim,m,nind
  integer, intent(inout) :: inform
  integer, intent(out)   :: hnnz

  ! ARRAY ARGUMENTS
  logical,      intent(in)  :: equatn(m),linear(m)
  integer,      intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in)  :: lambda(m),rho(m),x(nind)
  real(kind=8), intent(out) :: hval(lim)

  ! This subroutine computes the Hessian of the augmented Lagrangian
  ! in the reduced space.

  ! LOCAL SCALARS
  integer :: col,i,k,row

  ! LOCAL ARRAYS
  integer :: wi(nfull)

  ! Set xfull with values in x
  xfull(ind(1:nind)) = x(1:nind)

  ! Compute the Hessian of the augmented Lagrangian
  call sevalhal(nfull,xfull,m,lambda,rho,equatn,linear,hrow,hcol,hval, &
       hnnz,lim,inform)
  if ( inform .ne. 0 ) return

  ! Shrink representation of H
  wi(1:nfull) = 0
  wi(ind(1:nind)) = (/ (i, i=1,nind) /)

  k = 0
  do i = 1,hnnz
     row = wi(hrow(i))
     col = wi(hcol(i))

     if ( row .ne. 0 .and. col .ne. 0 ) then
        k = k + 1
        hrow(k) = row
        hcol(k) = col
        hval(k) = hval(i)
     end if
  end do

  hnnz = k

end subroutine calchal

! ******************************************************************
! ******************************************************************

subroutine calchalp(nind,x,m,lambda,rho,equatn,linear,p,hp,gothl,inform)

  use modrspace

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out)   :: gothl
  integer, intent(in)    :: m,nind
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)  :: equatn(m),linear(m)
  real(kind=8), intent(in)  :: lambda(m),p(nind),rho(m),x(nind)
  real(kind=8), intent(out) :: hp(nind)

  ! LOCAL ARRAYS
  real(kind=8) :: hpfull(nfull),pfull(nfull)

  ! Set xfull with values in x
  xfull(ind(1:nind)) = x(1:nind)

  ! Complete p with zeroes
  pfull(1:nfull) = 0
  pfull(ind(1:nind)) = p(1:nind)

  ! Compute the Hessian-vector product
  call sevalhalp(nfull,xfull,m,lambda,rho,equatn,linear,pfull,hpfull,gothl,inform)
  if ( inform .ne. 0 ) return

  ! Shrink hpfull to the reduced space
  hp(1:nind) = hpfull(ind(1:nind))

end subroutine calchalp

! ******************************************************************
! ******************************************************************

subroutine capplyhpre(nind,m,rho,equatn,gotp,r,z)

  use modrspace

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,nind
  logical, intent(inout) :: gotp

  ! ARRAY ARGUMENTS
  logical,      intent(in)  :: equatn(m)
  real(kind=8), intent(in)  :: r(nind),rho(m)
  real(kind=8), intent(out) :: z(nind)

  ! LOCAL ARRAYS
  real(kind=8) :: rfull(nfull),zfull(nfull)

  ! Complete r with zeroes
  rfull(1:nfull) = 0
  rfull(ind(1:nind)) = r(1:nind)

  ! Solve P zfull = rfull
  call applyhpre(nfull,m,rho,equatn,gotp,rfull,zfull)

  ! Shrink zfull to the reduced space
  z(1:nind) = zfull(ind(1:nind))

end subroutine capplyhpre

! ******************************************************************
! ******************************************************************

subroutine csetp(nind,x,inform)

  use modrspace
  use problvls, only: ssetp

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: nind
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(nind)

  ! Set xfull with values in x
  xfull(ind(1:nind)) = x(1:nind)

  ! Set point
  call ssetp(nfull,xfull,inform)
  if ( inform .ne. 0 ) return

end subroutine csetp

! ******************************************************************
! ******************************************************************

subroutine shrink(nind,v)

  use modrspace

  implicit none

  ! This subroutine shrinks vector v from the full dimension space
  ! (dimension n) to the reduced space (dimension nind).

  ! SCALAR ARGUMENTS
  integer, intent(in) :: nind

  ! ARRAY ARGUMENTS
  real(kind=8), intent(inout) :: v(nfull)

  ! LOCAL SCALARS
  integer :: i,indi
  real(kind=8) :: tmp

  do i = 1,nind
     indi = ind(i)
     if ( i .ne. indi ) then
        tmp     = v(indi)
        v(indi) = v(i)
        v(i)    = tmp
     end if
  end do

end subroutine shrink

! ******************************************************************
! ******************************************************************

subroutine expand(nind,v)

  use modrspace

  implicit none

  ! This subroutine expands vector v from the reduced space
  ! (dimension nind) to the full space (dimension n).

  ! SCALAR ARGUMENTS
  integer, intent(in) :: nind

  ! ARRAY ARGUMENTS
  real(kind=8), intent(inout) :: v(nfull)

  ! LOCAL SCALARS
  integer      :: i,indi
  real(kind=8) :: tmp

  do i = nind,1,- 1
     indi = ind(i)
     if ( i .ne. indi ) then
        tmp     = v(indi)
        v(indi) = v(i)
        v(i)    = tmp
     end if
  end do

end subroutine expand
