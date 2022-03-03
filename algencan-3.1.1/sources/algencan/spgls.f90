! ******************************************************************
! ******************************************************************

subroutine spgls(n,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg, &
     xp,fp,alpha,d,evalaldim,setpdim,lsinfo,inform)

  use modouttyp
  use problvlv, only: fcnt

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)       :: m,n
  integer, intent(inout)    :: inform
  integer, intent(out)      :: lsinfo
  real(kind=8), intent(in)  :: f,lamspg
  real(kind=8), intent(out) :: alpha,fp

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: g(n),l(n),lambda(m),rho(m),u(n)
  real(kind=8), intent(inout) :: x(n)
  real(kind=8), intent(out)   :: d(n),xp(n)

  ! SUBROUTINE ARGUMENTS
  external :: evalaldim,setpdim

  ! This subroutine computes a line search in the Spectral Projected
  ! Gradient direction.
  ! 
  ! lsinfo:
  ! 
  ! 0: Armijo satisfied
  ! 1: Small step with functional value similar to the current one
  ! 2: Unbounded objective function?
  ! 3: Too small backtracking step. Wrong gradient?

  ! LOCAL SCALARS
  integer      :: i
  real(kind=8) :: dsupn,gtd,xsupn

  ! ------------------------------------------------------------------
  ! Compute search direction, directional derivative, dsupn, xsupn and
  ! first trial
  ! ------------------------------------------------------------------

  gtd   = 0.0d0
  dsupn = 0.0d0
  xsupn = 0.0d0
  do i = 1,n
     d(i)  = - lamspg * g(i)
     xp(i) = x(i) + d(i)
     if ( xp(i) .lt. l(i) .or. xp(i) .gt. u(i) ) then
        xp(i) = max( l(i), min( xp(i), u(i) ) )
        d(i)  = xp(i) - x(i)
     end if
     gtd   = gtd + g(i) * d(i)
     dsupn = max( dsupn, abs( d(i) ) )
     xsupn = max( xsupn, abs( x(i) ) )
  end do

  if ( iprintinn .ge. 6 ) then
     write(* ,100) xsupn,lamspg,dsupn
     write(10,100) xsupn,lamspg,dsupn
  end if

  call setpdim(n,xp,inform)
  if ( inform .ne. 0 ) return

  call evalaldim(n,xp,m,lambda,rho,equatn,linear,fp,inform)
  if ( inform .ne. 0 ) return

  alpha = 1.0d0

  if ( iprintinn .ge. 6 ) then
     write(*, 110) alpha,fp,fcnt
     write(10,110) alpha,fp,fcnt
  end if

  ! ==================================================================
  ! Backtracking
  ! ==================================================================

  call backtracking(n,x,m,lambda,rho,equatn,linear,f,d,gtd,alpha,fp, &
       xp,evalaldim,setpdim,lsinfo,inform)
  if ( inform .ne. 0 ) return

  ! ==================================================================
  ! End of backtracking
  ! ==================================================================

  ! NON-EXECUTABLE STATEMENTS

 100  format(/,5X,'SPG Line search (xsupn = ',1P,D7.1,1X,'SPGstep= ', &
                   1P,D7.1,1X,'dsupn = ',1P,D7.1,')')
 110  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)

end subroutine spgls
