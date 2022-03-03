! ******************************************************************
! ******************************************************************

subroutine dogleg(n,g,hnnz,hlin,hcol,hval,l,pd,delta,p,dlinfo)

  use modmachconst
  use modouttyp

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(in)    :: pd
  integer,      intent(in)    :: hnnz,n
  integer,      intent(inout) :: dlinfo
  real(kind=8), intent(inout) :: delta,l

  ! ARRAY ARGUMENTS
  integer,      intent(in)    :: hcol(hnnz),hlin(hnnz)
  real(kind=8), intent(in)    :: g(n),hval(hnnz)
  real(kind=8), intent(inout) :: p(n)

  ! Compute approximate minimizer of the quadratic model using Dogleg
  ! method when the Hessian is positive definite and Cauchy point
  ! otherwise.

  ! dlinfo:
  !
  ! 0: successfull exit. Both H and g are null;
  ! 1: successfull exit.

  ! LOCAL SCALARS
  integer      :: col,i,lin
  real(kind=8) :: a,b,c,d,delta2,geucn,geucn2,gthg,pbeucn2, &
                  pueucn,pueucn2,putpb,r

  ! LOCAL ARRAYS
  real(kind=8) :: pu(n)

  ! Presentation

  if ( iprintinn .ge. 5 ) then
     write(* ,1000)
     write(10,1000)
  end if

  ! Initialization

  dlinfo = 1
  delta2 = delta**2

  ! If H is not positive definite, compute Cauchy point

  if ( .not. pd ) then
     go to 100
  end if

  pbeucn2 = 0.0d0
  do i = 1,n
     pbeucn2 = pbeucn2 + p(i)**2
  end do

  ! If Newton step is inside the trust region, this step is taken

  if ( pbeucn2 .le. delta2 ) then
     go to 500
  end if

  ! If Newton step is outside the trust region, compute the
  ! unconstrained minimizer pu of the quadratic function

  ! Compute g^T H g

  do i = 1,n
     pu(i) = 0.0d0
  end do

  do i = 1,hnnz
     lin = hlin(i)
     col = hcol(i)

     !     pu(lin) = pu(lin) + hval(i) * g(col)
     !     if ( lin .ne. col ) then
     !         pu(col) = pu(col) + hval(i) * g(lin)
     !     end if

     if ( lin .eq. col ) then
        pu(lin) = pu(lin) + ( hval(i) + l ) * g(col)
     else
        pu(lin) = pu(lin) + hval(i) * g(col)
        pu(col) = pu(col) + hval(i) * g(lin)
     end if
  end do

  gthg = 0.0d0
  do i = 1,n
     gthg = gthg + g(i) * pu(i)
  end do

  ! Compute g^T g.

  geucn2 = 0.0d0
  do i = 1,n
     geucn2 = geucn2 + g(i)**2
  end do

  ! Compute pu

  do i = 1,n
     pu(i) = - geucn2 * g(i) / gthg
  end do

  ! If uncontrained minimizer is outside the trust region, it is
  ! truncated at the border

  pueucn2 = 0.0d0
  do i = 1,n
     pueucn2 = pueucn2 + pu(i)**2
  end do
  pueucn = sqrt( pueucn2 )

  if ( pueucn2 .ge. delta2 ) then
     r = delta / pueucn
     do i = 1,n
        p(i) = r * pu(i)
     end do
     go to 500
  end if

  ! Compute step length in directions pu and pb. Direction p is a
  ! linear combination of pu and pb. To compute the step length in
  ! each direction, we have to solve (in r):
  ! \| pu + (r - 1) (pb - pu) \|^2 = delta^2.

  putpb = 0.0d0
  do i = 1,n
     putpb = putpb + pu(i) * p(i)
  end do

  a = pbeucn2 + pueucn2 - 2.0d0 * putpb
  b = - 2.0d0 * pbeucn2 - 4.0d0 * pueucn2 + 6.0d0 * putpb
  c = pbeucn2 + 4.0d0 * pueucn2 - 4.0d0 * putpb - (delta**2)
  d = b**2 - 4.0d0 * a * c

  if ( ( 2.0d0 * abs( a ) .lt. macheps23 ) .or. &
       ( d .lt. 0.0d0 ) ) then
     go to 200
  end if

  r = (- b + sqrt( d )) / (2.0d0 * a)

  if ( ( r .lt. 1.0d0 ) .or. ( r .gt. 2.0d0 ) ) then
     r = (- b - sqrt( d )) / (2.0d0 * a)
  end if
  r = r - 1.0d0

  do i = 1,n
     p(i) = r * p(i) + (1.0d0 - r) * pu(i)
  end do

  go to 500

  ! Compute Cauchy point, when H is not positive definite

100 continue

  ! Compute g^T H g

  do i = 1,n
     pu(i) = 0.0d0
  end do

  do i = 1,hnnz
     lin = hlin(i)
     col = hcol(i)

     !     pu(lin) = pu(lin) + hval(i) * g(col)
     !     if ( lin .ne. col ) then
     !         pu(col) = pu(col) + hval(i) * g(lin)
     !     end if

     if ( lin .eq. col ) then
        pu(lin) = pu(lin) + ( hval(i) + l ) * g(col)
     else
        pu(lin) = pu(lin) + hval(i) * g(col)
        pu(col) = pu(col) + hval(i) * g(lin)
     end if
  end do

  gthg = 0.0d0
  do i = 1,n
     gthg = gthg + g(i) * pu(i)
  end do

  ! Compute g^T g

  geucn2 = 0.0d0
  do i = 1,n
     geucn2 = geucn2 + g(i)**2
  end do
  geucn = sqrt( geucn2 )

  ! Compute step length

200 if ( abs( gthg ) .le. macheps23 .and. geucn .le. macheps23 ) then

     dlinfo = 0

     if ( iprintinn .ge. 5 ) then
        write(* ,1020)
        write(10,1020)
     end if

     return
  end if

  ! Compute p = coef * g

  if ( gthg .le. 0.0d0 .or. geucn2 * geucn .lt. delta * gthg ) then
     do i = 1,n
        p(i) = - delta * g(i) / geucn
     end do
  else
     do i = 1,n
        p(i) = - geucn2 * g(i) / gthg
     end do
  end if

  ! Termination

500 continue

  if ( iprintinn .ge. 5 ) then
     write(* ,1010)
     write(10,1010)
  end if

  if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
     write(*, 1030) min0(n,nprint),(p(i),i=1,min0(n,nprint))
     write(10,1030) min0(n,nprint),(p(i),i=1,min0(n,nprint))
  end if

 1000 format(/,5X,'Computation of Dogleg direction.')
 1010 format(  5X,'Dogleg computed successfully.')
 1020 format(  5X,'Null direction was computed.')
 1030 format(/,5X,'Dogleg direction (first ',I7,' components): ', &
             /,1(5X,6(1X,1P,D11.4)))

end subroutine dogleg

! ******************************************************************
! ******************************************************************

! Compute approximate minimizer of the quadratic model using Dogleg
! method when the Hessian is positive definite and Cauchy point
! otherwise.

! On Entry
!
! n        integer
!          dimension
!
! g        double precision g(n)
!          vector used to define the quadratic function
!
! hnnz     integer
!          number of nonzero elements of H
!
! hlin     integer hlin(hnnz)
!          row indices of nonzero elements of H
!
! hcol     integer hcol(hnnz)
!          column indices of nonzero elements of H
!
! hval     double precision hval(hnnz)
!          nonzero elements of H, that is,
!          H(hlin(i),hcol(i)) = hval(i). Since H is symmetric, just
!          one element H(i,j) or H(j,i) must appear in its sparse
!          representation. If more than one pair corresponding to
!          the same position of H appears in the sparse
!          representation, the multiple entries will be summed.
!          If any hlin(i) or hcol(i) is out of range, the entry will
!          be ignored
!
! pd       logical
!          indicates if the last Cholesky decomposition of moresor
!          was successfull. That is, if the last matrix used by
!          moresor was positive definite
!
! delta    double precision
!          trust-region radius
!
! On Return
!
! p        double precision p(n)
!          solution to problem
!          minimize     psi(w)
!          subjected to ||w|| <= delta
!
! dlinfo   integer
!          This output parameter tells what happened in this
!          subroutine, according to the following conventions:
!
!          0 = successfull exit. Both H and g are null;
!
!          1 = successfull exit.
