! *****************************************************************
! *****************************************************************

subroutine extrapolation(n,nind,x,l,u,m,lambda,rho,equatn,linear,g, &
     xp,fp,gp,d,alpha,amax,rbdnnz,rbdind,rbdtype,fmin,beta,etaext,  &
     maxextrap,extinfo,inform)

  use modmachconst
  use modouttyp
  use problvlv, only: fcnt

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: m,maxextrap,n,nind,rbdnnz
  integer,      intent(inout) :: inform
  integer,      intent(out)   :: extinfo
  real(kind=8), intent(in)    :: amax,beta,etaext,fmin
  real(kind=8), intent(inout) :: alpha,fp

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  character,    intent(in)    :: rbdtype(nind)
  integer,      intent(in)    :: rbdind(nind)
  real(kind=8), intent(in)    :: d(nind),g(nind),l(nind),lambda(m), &
                                 rho(m),u(nind),x(nind)
  real(kind=8), intent(inout) :: xp(nind)
  real(kind=8), intent(out)   :: gp(n) ! On return: shrinked full-space gradient

  ! Performs extrapolation.
  !
  ! extinfo:
  !
  ! 0: Success
  ! 2: Unbounded objective function?
  ! 4: beta-condition holds. No extrapolation is done
  ! 6: Maximum number of extrapolations reached
  ! 7: Similar consecutive projected points
  ! 8: Not-well-defined objective function
  ! 9: Functional value increases

  ! LOCAL SCALARS
  logical      :: projected,samep
  integer      :: extrap,i
  real(kind=8) :: atmp,fbext,ftmp,gptd,gtd

  ! LOCAL ARRAYS
  real(kind=8) :: xbext(n),xtmp(n)

  extinfo = 0

  extrap  = 0

  ! Compute directional derivative

  gtd = 0.0d0
  do i = 1,nind
     gtd = gtd + g(i) * d(i)
  end do

  call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
  if ( inform .ne. 0 ) return

  gptd = 0.0d0
  do i = 1,nind
     gptd = gptd + gp(i) * d(i)
  end do

  ! Beta condition holds. No extrapolation is done.

  if ( gptd .ge. beta * gtd ) then

     if ( iprintinn .ge. 6 ) then
        write(*, 110)
        write(10,110)
     end if

     extinfo = 4

     return

  end if

  ! Beta condition does not holds. We will extrapolate.

  if ( iprintinn .ge. 6 ) then
     write(* ,120)
     write(10,120)
  end if

  ! Save f and x before extrapolation to return in case of a
  ! not-well-defined objective function at an extrapolated point

  fbext = fp

  do i = 1,nind
     xbext(i) = xp(i)
  end do

1010 continue

  ! Test f going to -inf

  if ( fp .le. fmin ) then

     ! Finish the extrapolation with the current point

     if ( extrap .ne. 0 ) then
        call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
        if ( inform .ne. 0 ) return
     end if

     extinfo = 2

     if ( iprintinn .ge. 6 ) then
        write(*, 910)
        write(10,910)
     end if

     return

  end if

  ! Test if the maximum number of extrapolations was exceeded

  if ( extrap .ge. maxextrap ) then

     ! Finish the extrapolation with the current point

     if ( extrap .ne. 0 ) then
        call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
        if ( inform .ne. 0 ) return
     end if

     extinfo = 6

     if ( iprintinn .ge. 6 ) then
        write(*, 930)
        write(10,930)
     end if

     return

  end if

  ! Chose new step

  extrap = extrap + 1

  if ( alpha .lt. amax .and. etaext * alpha .gt. amax ) then
     atmp = amax
  else
     atmp = etaext * alpha
  end if

  ! Compute new trial point

  do i = 1,nind
     xtmp(i) = x(i) + atmp * d(i)
  end do

  if ( atmp .eq. amax ) then
     do i = 1,rbdnnz
        if ( rbdtype(i) .eq. 'L' ) then
           xtmp(rbdind(i)) = l(rbdind(i))
        else if ( rbdtype(i) .eq. 'U' ) then
           xtmp(rbdind(i)) = u(rbdind(i))
        end if
     end do
  end if

  ! Project

  if ( atmp .gt. amax ) then

     projected = .false.
     do i = 1,nind
        if ( xtmp(i) .lt. l(i) .or. xtmp(i) .gt. u(i) ) then
           projected = .true.
           xtmp(i) = max( l(i), min( xtmp(i), u(i) ) )
        end if
     end do

     ! Test if this is not the same point as the previous one. This test
     ! is performed only when xtmp is in fact a projected point.

     if ( projected ) then

        samep = .true.
        do i = 1,nind
           if ( xtmp(i) .gt. &
                xp(i) + macheps23 * max( 1.0d0, abs( xp(i) ) ) .or. &
                xtmp(i) .lt. &
                xp(i) - macheps23 * max( 1.0d0, abs( xp(i) ) ) ) then
              samep = .false.
           end if
        end do

        if ( samep ) then

           ! Finish the extrapolation with the current point

           call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
           if ( inform .ne. 0 ) return

           extinfo = 7

           if ( iprintinn .ge. 6 ) then
              write(*, 940)
              write(10,940)
           end if

           return

        end if

     end if

  end if
  
  ! Evaluate function

  call csetp(nind,xtmp,inform)
  if ( inform .ne. 0 ) return

  call calcal(nind,xtmp,m,lambda,rho,equatn,linear,ftmp,inform)

  ! If the objective function is not well defined in an extrapolated
  ! point, we discard all the extrapolated points and return to a
  ! safe region (to the last point before the extrapolation)

  if ( inform .ne. 0 ) then

     fp = fbext

     do i = 1,nind
        xp(i) = xbext(i)
     end do

     call csetp(nind,xp,inform)
     if ( inform .ne. 0 ) return

     call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
     if ( inform .ne. 0 ) return

     extinfo = 8

     if ( iprintinn .ge. 6 ) then
        write(*, 950)
        write(10,950)
     end if

     return

  end if

  ! Print information of this iteration

  if ( iprintinn .ge. 6 ) then
     write(*, 100) atmp,ftmp,fcnt
     write(10,100) atmp,ftmp,fcnt
  end if

  ! If the functional value decreases then set the current point and
  ! continue the extrapolation

  if ( ftmp .lt. fp ) then

     alpha = atmp

     fp = ftmp

     do i = 1,nind
        xp(i) = xtmp(i)
     end do

     go to 1010

  end if

  ! If the functional value does not decrease then discard the last
  ! trial and finish the extrapolation with the previous point

  call csetp(nind,xp,inform)
  if ( inform .ne. 0 ) return

  call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
  if ( inform .ne. 0 ) return

  extinfo = 9

  if ( iprintinn .ge. 6 ) then
     write(*, 960)
     write(10,960)
  end if

! NON-EXECUTABLE STATEMENTS

 100  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 110  format(  5X,'Beta condition also holds. ', &
                  'No extrapolation is done.')
 120  format(  5X,'Beta condition does not hold. We will extrapolate.')

 910  format(  5X,'Flag of Extrapolation: Unbounded objective ', &
                  'function?')
 930  format(  5X,'Flag of Extrapolation: Maximum of consecutive ', &
                  'extrapolations reached.')
 940  format(  5X,'Flag of Extrapolation: Very similar consecutive ', &
                  'projected points.')
 950  format(  5X,'Flag of Extrapolation: Not-well-defined objective ', &
                  'function in an extrapolated point.')
 960  format(  5X,'Flag of Extrapolation: Functional value increased ', &
                  'when extrapolating.')

end subroutine extrapolation
