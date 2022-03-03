! ******************************************************************
! ******************************************************************

subroutine backtracking(dim,x,m,lambda,rho,equatn,linear,f,d,gtd, &
     alpha,fp,xp,evalaldim,setpdim,btinfo,inform)

  use modmachconst
  use modalgconst
  use modouttyp
  use problvlv, only: fcnt

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: dim,m
  integer,      intent(inout) :: inform
  integer,      intent(out)   :: btinfo
  real(kind=8), intent(in)    :: f,gtd
  real(kind=8), intent(inout) :: alpha,fp

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: d(dim),lambda(m),rho(m),x(dim)
  real(kind=8), intent(inout) :: xp(dim)

  ! SUBROUTINE ARGUMENTS
  external evalaldim,setpdim

  ! Backtracking with quadratic interpolation.
  !
  ! btinfo:
  !
  ! 0: Armijo satisfied
  ! 2: Unbounded objective function?
  ! 3: Too small backtracking step. Wrong gradient?

  ! LOCAL SCALARS
  logical      :: samep
  integer      :: i,interp
  real(kind=8) :: atmp

  interp = 0

2010 continue

  ! Test Armijo condition

  if ( fp .le. f + alpha * gamma * gtd ) then

     ! Finish backtracking with the current point

     btinfo = 0

     if ( iprintinn .ge. 6 ) then
        write(*, 900)
        write(10,900)
     end if

     return

  end if

  ! Test f going to -inf

  if ( fp .le. fmin ) then

     ! Finish backtracking with the current point

     btinfo = 2

     if ( iprintinn .ge. 6 ) then
        write(*, 920)
        write(10,920)
     end if

     return

  end if

  ! Test if we obtained a functional value similar to the current one
  ! associated to a very small step

  samep = .true.
  do i = 1,dim
     if ( xp(i) .gt. x(i) + macheps23 * max(1.0d0,abs(x(i))) .or. &
          xp(i) .lt. x(i) - macheps23 * max(1.0d0,abs(x(i))) ) then
        samep = .false.
     end if
  end do

  if ( samep .and. fp .le. f + macheps23 * abs(f) ) then

     ! Finish backtracking with the current point

     btinfo = 3

     if ( iprintinn .ge. 6 ) then
        write(*, 930)
        write(10,930)
     end if

     return

  end if

  ! Compute new step

  interp = interp + 1

  atmp = ( - gtd * alpha ** 2 ) / ( 2.0d0 * ( fp - f - alpha * gtd ) )

  if ( atmp .ge. sigma1 * alpha .and. &
       atmp .le. sigma2 * alpha ) then
     alpha = atmp
  else
     alpha = alpha / etaint
  end if

  ! Compute new trial point

  do i = 1,dim
     xp(i) = x(i) + alpha * d(i)
  end do

  call setpdim(dim,xp,inform)
  if ( inform .ne. 0 ) return

  call evalaldim(dim,xp,m,lambda,rho,equatn,linear,fp,inform)
  if ( inform .ne. 0 ) return

  ! Print information of this iteration

  if ( iprintinn .ge. 6 ) then
     write(*, 110) alpha,fp,fcnt
     write(10,110) alpha,fp,fcnt
  end if

  go to 2010
  
  ! NON-EXECUTABLE STATEMENTS
110 format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
900 format(  5X,'Flag of backtracking: Armijo condition holds.')
920 format(  5X,'Flag of backtracking: Unbounded objective function?')
930 format(  5X,'Flag of backtracking: Too small backtracking step.')

end subroutine backtracking
