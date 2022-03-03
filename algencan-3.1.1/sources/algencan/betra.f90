! *****************************************************************
! *****************************************************************

subroutine betra(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,    &
     trdelta,newdelta,mslamb,epsg,xeucn,d,xp,fp,gp,triter,chcnt, &
     memfail,betinfo,inform)

  use probgiven, only: hnnzlim
  use modmachconst
  use modalgconst
  use moditetyp
  use modouttyp
  use problvlv, only: fcnt

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out)   :: memfail
  integer,      intent(in)    :: m,n,nind
  integer,      intent(out)   :: betinfo,chcnt,triter
  integer,      intent(inout) :: inform
  real(kind=8), intent(in)    :: epsg,f,xeucn
  real(kind=8), intent(inout) :: mslamb,newdelta,trdelta
  real(kind=8), intent(out)   :: fp

  ! ARRAY ARGUMENTS
  logical,      intent(in)  :: equatn(m),linear(m)
  real(kind=8), intent(in)  :: g(n),l(nind),lambda(m),rho(m),u(nind),x(nind)
  real(kind=8), intent(out) :: d(nind),gp(n),xp(nind) 

  ! (On return: gp must be the shrinked full-space gradient)

  ! Solves the "unconstrained" inner problem that arises from the 
  ! active-set method to solve box-constrained minimization problem
  !
  !        Minimize f(x) subject to l <= x <= u
  ! 
  ! described in 
  !
  ! M. Andretta, E. G. Birgin e J. M. Martínez. "Practical active-set 
  ! Euclidian trust-region method with spectral projected gradients 
  ! for bound-constrained minimization". Optimization 54, pp. 
  ! 305-325, 2005.

  ! betinfo:
  !
  ! 0: Sufficient decrease of the function
  ! 1: x hits the boundary
  ! 2: Unbounded objective function?
  ! 3: Either the search direction or the step lenght is too small
  ! 4: Trust-region radius too small
  ! 5: Undefined direction
  ! 6: x is a first-order stationary point close to the boundary
  ! 7: x is a second-order stationary point
  ! 8: x is too close to the boundary

  ! LOCAL SCALARS
  logical      :: frstit,pd,samep
  integer      :: dinfo,extinfo,hnnz,i,col,row,index,rbdnnz
  real(kind=8) :: amax,ared,dbound,deucn,gsupn,phi,pred,stpl,tmp

  ! LOCAL ARRAYS
  character    :: rbdtype(nind)
  integer      :: hdiag(n),hcol(hnnzlim+n),hrow(hnnzlim+n),rbdind(nind)
  real(kind=8) :: hval(hnnzlim+n)

  ! Print presentation information

  if ( iprintinn .ge. 5 ) then
     write(*, 1000)
     write(10,1000)
  end if

  ! Initialization

  if ( .not. sameface ) then
     trdelta = max( trdelta, newdelta )
  end if

  newdelta = 0.0d0

  memfail  = .false.

  frstit   = .true.

  ! ==================================================================
  ! Compute distance to the boundary
  ! ==================================================================

  ! step 1: calculate the distance between x and the boundary.
  !         dbound is the largest positive delta such that the ball 
  !         centered at x with radius delta is still inside the box 
  !         (set of constraints).

  dbound = bignum

  do i = 1,nind
     dbound = min( dbound, x(i) - l(i) )
  end do

  do i = 1,nind
     dbound = min( dbound, u(i) - x(i) )
  end do

  ! Calculate infinite-norm of gradient.

  gsupn = 0.0d0
  do i = 1,nind
     gsupn = max( gsupn, abs( g(i) ) )
  end do

  ! ==================================================================
  ! Close to the boundary: perform inner SPG iteration
  ! ==================================================================

  ! step 2: close to the boundary, stop.

  if ( dbound .lt. 2.0d0 * trdelmin ) then

     xp(1:nind) = x(1:nind)
     gp(1:n) = g(1:n)

     fp  = f

     betinfo = 8

     if ( iprintinn .ge. 5 ) then
        write(*, 9010) 
        write(10,9010) 
     end if

     return

  end if

  ! ==================================================================
  ! Far from the boundary: perform trust-region iteration
  ! ==================================================================

  triter = triter + 1

  ! step 3: far from the boundary, solve trust-region subproblem. 
  !         Evaluate function Hessian at x.

  call calchal(nind,x,m,lambda,rho,equatn,linear,hrow,hcol,hval,hnnz,hnnzlim,inform)
  if ( inform .ne. 0 ) return

  do i = 1,nind
     hdiag(i) = 0
  end do

  do i = 1,hnnz
     row = hrow(i)
     col = hcol(i)

     if ( row .eq. col ) then
        if ( hdiag(row) .eq. 0 ) then
           hdiag(row) = i
        else
           hval(hdiag(row)) = hval(hdiag(row)) + hval(i)
           hval(i) = 0.0d0
        end if
     end if
  end do

  do i = 1,nind
     if ( hdiag(i) .eq. 0 ) then
        hnnz       = hnnz + 1      
        hrow(hnnz) = i
        hcol(hnnz) = i
        hval(hnnz) = 0.0d0
        hdiag(i)   = hnnz
     end if
  end do

  ! step 4: solve the trust-region subproblem using More-Sorensen's 
  !         algorithm to minimize "exactly" quadratics subjected to 
  !         balls. 

  ! If trust-region radius is too small, the inner algorithm stops.

100 continue

  if ( trdelta .lt. macheps * max( 1.0d0, xeucn ) ) then

     xp(1:nind) = x(1:nind)
     gp(1:n) = g(1:n)

     fp  = f

     betinfo = 4

     if ( iprintinn .ge. 5 ) then
        write(*, 9040) 
        write(10,9040) 
     end if

     return
  end if

  call moresor(nind,g,hnnz,hrow,hcol,hval,hdiag,trdelta,mssig, &
       0.0d0,mseps,msmaxit,mslamb,pd,d,chcnt,memfail,dinfo)

  if ( memfail ) then
     return
  end if

  ! If maximum allowed number of MEQB iterations is achieved, another
  ! direction d is calculated.

  if ( dinfo .eq. 5 ) then 

     if ( iprintinn .ge. 5 ) then
        write(*, 2000)
        write(10,2000)
     end if

     call dogleg(nind,g,hnnz,hrow,hcol,hval,mslamb,pd,trdelta,d,dinfo)

  end if

  ! If both internal gradient and Hessian matrix are null, subroutines 
  ! MEQB and dogleg stop with dinfo = 0 and then the inner algorithm 
  ! stops declaring "second-order stationary point".

  if ( dinfo .eq. 0 ) then

     xp(1:nind) = x(1:nind)
     gp(1:n) = g(1:n)

     fp  = f

     betinfo = 7

     if ( iprintinn .ge. 5 ) then
        write(*, 9020) 
        write(10,9020) 
     end if

     return
  end if

  ! Print direction

  if ( iprintinn .ge. 5 ) then
     write(*, 1020) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
     write(10,1020) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
  end if

  ! Evaluate the quadratic model of the objective function at d.

  call squad(nind,d,g,hnnz,hrow,hcol,hval,phi)

  ! If the value of the quadratic model at d is 0 it means that x is a 
  ! second-order stationary point. In this case, inner algorithm stops
  ! declaring this.

  if ( ( abs( phi ) .le. phieps ) .and. ( gsupn .le. epsg ) ) then

     xp(1:nind) = x(1:nind)
     gp(1:n) = g(1:n)

     fp  = f

     betinfo = 7

     if ( iprintinn .ge. 5 ) then
        write(*, 9020) 
        write(10,9020) 
     end if

     return
  end if

  ! Calculate predicted decrease of objective function

  pred = abs( phi )

  ! Calculate d Euclidian-norm

  deucn = 0.0d0
  do i = 1,nind
     deucn = deucn + d(i)**2
  end do
  deucn = sqrt( deucn )

  ! To avoid NaN and Inf directions

  if ( .not. ( deucn .le. bignum ) ) then

     trdelta = 0.25d0 * trdelta

     if ( iprintinn .ge. 5 ) then
        write(*, 2060) 
        write(10,2060) 
     end if

     if ( trdelta .lt. macheps * max( 1.0d0, xeucn ) ) then

        xp(1:nind) = x(1:nind)
        gp(1:n) = g(1:n)

        fp  = f

        betinfo = 5

        if ( iprintinn .ge. 5 ) then
           write(*, 9070) 
           write(10,9070) 
        end if

        return
     end if

     triter  = triter + 1
     frstit  = .false.

     go to 100

  end if

  ! Calculate point xp = x + d.

  xp(1:nind) = x(1:nind) + d(1:nind)

  stpl = 1.0d0

  ! Verify if xp is inside de box. If not, interior is set to
  ! false and amax is set to the biggest positive scalar such that
  ! x + amax*d is inside the box.

  call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)

  ! ==================================================================
  ! Point on the boundary
  ! ==================================================================      

  ! If xp is not interior to the box, xp = x + d is replaced
  ! by xp = x + amax*d. Now xp is definitely interior. Actually,
  ! it is in the boundary. If the objective function decreases in
  ! xp, the inner algorithm stops with xp as a solution.
  ! Otherwise, a new trust-region radius trdelta is chosen (smaller than
  ! dbound) and a new quadratic model minimizer is calculated (which
  ! is necessarily interior because of the choice of trdelta).

  if ( amax .le. 1.0d0 ) then 

     if ( iprintinn .ge. 5 ) then
        write(*, 2010) 
        write(10,2010) 
     end if

     xp(1:nind) = x(1:nind) + amax * d(1:nind)

     stpl = amax

     ! Set x(i) to l(i) or u(i) for the indices i that got to the 
     ! boundary (to prevent errors).

     do i = 1,rbdnnz
        index = rbdind(i)
        if ( rbdtype(i) .eq. 'L' ) then
           xp(index) = l(index)
        elseif ( rbdtype(i) .eq. 'U' ) then
           xp(index) = u(index)
        end if
     end do

     call csetp(nind,xp,inform)
     if ( inform .ne. 0 ) return

     call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
     if ( inform .ne. 0 ) return

     ! Print functional value

     if ( iprintinn .ge. 5 ) then
        write(*, 1030) trdelta,fp,fcnt
        write(10,1030) trdelta,fp,fcnt
     end if

     ! Test whether f is very small

     if ( fp .le. fmin ) then

        call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
        if ( inform .ne. 0 ) return

        betinfo = 2

        if ( iprintinn .ge. 5 ) then
           write(*, 9050) 
           write(10,9050) 
        end if

        return
     end if

     ! If the new point x + d is too close to the previous point x, 
     ! inner algorithm stops

     samep = .true.
     do i = 1,nind
        if ( xp(i) .gt. x(i) +  macheps * max( 1.0d0, abs( x(i) ) ) .or. &
             xp(i) .lt. x(i) -  macheps * max( 1.0d0, abs( x(i) ) ) ) then
           samep = .false.
        end if
     end do

     if ( samep .and. fp .le. f + macheps23 * abs(f) ) then

        betinfo = 3

        if ( iprintinn .ge. 5 ) then
           write(*, 9060) 
           write(10,9060) 
        end if

        return
     end if

     ! Test if function value decreases at xp

     if ( ( fp .le. f ) .or.  ( ( deucn .le. macheps23 * xeucn ) .and. &
          ( fp .le. f + macheps23 * abs( f ) ) ) ) then

        if ( iprintinn .ge. 5 ) then
           write(*, 2030) 
           write(10,2030) 
        end if

        if ( extrp4 ) then

           call extrapolation(n,nind,x,l,u,m,lambda,rho,equatn, &
                linear,g,xp,fp,gp,d,stpl,amax,rbdnnz, &
                rbdind,rbdtype,fmin,beta,etaext,maxextrap,extinfo, &
                inform)

           if ( inform .ne. 0 ) return

           if ( extinfo .eq. 2 ) then
              betinfo = 2
           else
              betinfo = 0
           end if

           ! Update the trust-region radius (which may or may not be used)

           if ( frstit ) then
              trdelta = 2.0d0 * max( trdelta, stpl * deucn )
           end if

           if ( betinfo .eq. 2 ) then
              if ( iprintinn .ge. 5 ) then
                 write(*, 9050) 
                 write(10,9050) 
              end if

              return
           end if

           betinfo = 1

           ! Calculate actual reduction of objective function

           ared = f - fp

           go to 200

        else

           ! Update the trust-region radius (which may or may not be used)

           trdelta = 2.0d0 * trdelta

           ! Calculate actual reduction of objective function

           ared = f - fp

           ! Compute gp.

           call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
           if ( inform .ne. 0 ) return

           betinfo = 1

           go to 200
        end if

     else
        tmp      = trdelmin + msrho * ( ( dbound / (1.0d0 + mssig ) ) - trdelmin )
        newdelta = trdelta
        trdelta  = max( trdelmin, tmp )

        triter   = triter + 1
        frstit   = .false.

        if ( iprintinn .ge. 5 ) then
           write(*, 2040) 
           write(10,2040) 
        end if

        go to 100
     end if
  end if

  ! ==================================================================
  ! Point interior to the box
  ! ==================================================================

  ! step 5: in this case xp is inside the box. Acceptance or
  !         rejection of the trust-region subproblem solution.

  call csetp(nind,xp,inform)
  if ( inform .ne. 0 ) return

  call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
  if ( inform .ne. 0 ) return

  ! Print functional value

  if ( iprintinn .ge. 5 ) then
     write(*, 1030) trdelta,fp,fcnt
     write(10,1030) trdelta,fp,fcnt
  end if

  ! Test whether f is very small

  if ( fp .le. fmin ) then

     call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
     if ( inform .ne. 0 ) return

     betinfo = 2

     if ( iprintinn .ge. 5 ) then
        write(*, 9050) 
        write(10,9050) 
     end if

     return
  end if

  ! If the new point x + d is too close to the previous point x, inner
  ! algorithm stops

  samep = .true.
  do i = 1,nind
     ! if ( abs( x(i) - xp(i) ) .gt. macheps * max( abs( xp(i) ), 1.0d0 ) ) then
     if ( xp(i) .gt. x(i) + macheps * max(1.0d0,abs(x(i))) .or. &
          xp(i) .lt. x(i) - macheps * max(1.0d0,abs(x(i))) )  &
          then
        samep = .false.
     end if
  end do

  ! if ( samep .and. fp - f .le. macheps23 * abs(f) ) then
  if ( samep .and. fp .le. f + macheps23 * abs(f) ) then

     betinfo = 3

     if ( iprintinn .ge. 5 ) then
        write(*, 9060) 
        write(10,9060) 
     end if

     return
  end if

  ! Calculate actual reduction of objective function

  ared = f - fp

  ! If there is not sufficient decrease of the function, the 
  ! trust-region radius is decreased and the new quadratic model 
  ! minimizer will be calculated. 

  if ( iprintinn .ge. 5 ) then
     write(*, 1010) deucn,pred,ared
     write(10,1010) deucn,pred,ared
  end if

  samep = .true.
  do i = 1,nind
     if ( xp(i) .gt. x(i)+macheps23*max(1.0d0,abs(x(i))) .or. &
          xp(i) .lt. x(i)-macheps23*max(1.0d0,abs(x(i))) ) then
        samep = .false.
     end if
  end do

  ! if ( deucn .le. macheps23 * xeucn .and. ared  .le. macheps23 * abs( f ) ) then
  ! if ( samep .and. ared .le. macheps23 * abs(f) ) then
  if ( samep .and. fp .ge. f - macheps23 * abs(f) ) then

     call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
     if ( inform .ne. 0 ) return

     go to 200

  end if

  if ( ( pred .le. phieps .or. ared .ge. tralpha*pred ) .and. &
       ( pred .gt. phieps .or. f .ge. fp ) ) then

     ! If extrapolation at step 5 is not to be performed, point xp is accepted.

     if ( extrp5 ) then

        call extrapolation(n,nind,x,l,u,m,lambda,rho,equatn,linear, &
             g,xp,fp,gp,d,stpl,amax,rbdnnz,rbdind,rbdtype, &
             fmin,beta,etaext,maxextrap,extinfo,inform)

        if ( inform .ne. 0 ) return

        if ( extinfo .eq. 2 ) then
           betinfo = 2
        else
           betinfo = 0
        end if

        if ( frstit ) then
           trdelta = max( trdelta, stpl * deucn )
        end if

        if ( betinfo .eq. 2 ) then
           if ( iprintinn .ge. 5 ) then
              write(*, 9050) 
              write(10,9050) 
           end if

           return
        end if

        ! Update actual reduction of objective function

        ared = f - fp

     else

        call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
        if ( inform .ne. 0 ) return

     end if

  else

     ! trdelta = 0.25d0 * deucn
     ! due to numerical issues, we may have deucn > trdelta
     trdelta = 0.25d0 * min( trdelta, deucn )
     triter  = triter + 1
     frstit  = .false.

     if ( iprintinn .ge. 5 ) then
        write(*, 2050) 
        write(10,2050) 
     end if

     go to 100

  end if

  ! ==================================================================
  ! Prepare for next call to this routine
  ! ==================================================================

  ! Update the trust-region radius (which may or may not be used). 
  ! This update can only be done when the current iteration was a 
  ! trust-region iteration (and not an inner SPG one).

200 continue

  if ( ared .lt. 0.25d0 * pred ) then
     trdelta = max( 0.25d0 * deucn, trdelmin )
  else
     if ( ( ared .ge. 0.5d0 * pred ) .and. &
          ( deucn .ge. trdelta - macheps23 * & 
          max( trdelta, 1.0d0 ) ) ) then
        trdelta = max( 2.0d0 * trdelta, trdelmin )
     else
        trdelta = max( trdelta, trdelmin )
     end if
  end if

  ! If new point x is in the boundary, the inner algorithm stops
  ! and returns x as solution.

  if ( stpl .ge. amax ) then

     betinfo = 1

     if ( iprintinn .ge. 5 ) then
        write(*, 9030) 
        write(10,9030) 
     end if

  else

     betinfo = 0

     if ( iprintinn .ge. 5 ) then
        write(*, 9000) 
        write(10,9000) 
     end if
  end if

! Non-executable statements
      
 1000 format(/,5X,'Trust-region iteration')
 1010 format(  5X,'deucn = ',1P,D7.1,' pred = ',1P,D11.4, &
                  ' ared = ',1P,D11.4)
 1020 format(/,5X,'Trust-region direction (first ',I7, &
                  ' components): ',/,1(5X,6(1X,1P,D11.4)))
 1030 format(/,5X,'delta = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)

 2000 format(/,5X,'Since the direction More-Sorensen calculation ', &
             /,5X,'failed, dogleg direction will be computed.')
 2010 format(/,5X,'x+d is not interior to the box.')
 2030 format(  5X,'f(x+d) < f(x), x+d is accepted.')
 2040 format(  5X,'f(x+d) >= f(x), a new direction d will be computed.')
 2050 format(  5X,'x+d did not obtain suficcient functional reduction. ' &
             /,5X,'A new direction d will be computed.')
 2060 format(/,5X,'Direction is undefined. ', &
                  'A new direction d will be computed.')
 2070 format(/,5X,'Point close to the boundary. ', &
                  'An inner SPG iteration will be used.')

 9000 format(  5X,'Flag of TR: Sufficient decrease of function.')
 9010 format(  5X,'Flag of TR: ', &
                  'First-order stationary point close to boundary.')
 9020 format(/,5X,'Flag of TR: Second-order stationary point.')
 9030 format(  5X,'Flag of TR: Point on the boundary.')      
 9040 format(  5X,'Flag of TR: Trust-region radius too small.')
 9050 format(  5X,'Flag of TR: Unbounded objective function?')
 9060 format(  5X,'Flag of TR: Very similar consecutive points.')
 9070 format(  5X,'Flag of TR: Undefined direction.')

end subroutine betra

! ******************************************************************
! ******************************************************************

subroutine squad(nred,x,g,hnnz,hrow,hcol,hval,phi)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: hnnz,nred
  real(kind=8), intent(out) :: phi

  ! ARRAY ARGUMENTS
  integer,      intent(in) :: hcol(hnnz),hrow(hnnz)
  real(kind=8), intent(in) :: g(nred),hval(hnnz),x(nred)

  ! Evaluates the quadratic model phi(x) = 1/2 x^T H x + g^T x.

  ! LOCAL SCALARS
  integer :: col,i,row

  ! LOCAL ARRAYS
  real(kind=8) :: wd(nred)

  do i = 1,nred
     wd(i) = 0.0d0
  end do

  do i = 1,hnnz
     row = hrow(i)
     col = hcol(i)

     wd(row) = wd(row) + hval(i) * x(col)
     if ( row .ne. col ) then
        wd(col) = wd(col) + hval(i) * x(row)
     end if
  end do

  phi = 0.0d0
  do i = 1,nred
     phi = phi + wd(i) * x(i)
  end do

  phi = phi * 0.5d0

  do i = 1,nred
     phi = phi + g(i) * x(i)
  end do

end subroutine squad

! *****************************************************************
! *****************************************************************

! Algorithm that finds a unconstrained minimizer of objective 
! function inside the box of constraints, hits the boundary 
! (obtaining function decrease), or finds an interior point where 
! the objective function has sufficient decrease (compared to its 
! value at x). Extrapolation may be done.
!
! When the current point x is "close to" the boundary, a Spectral 
! Projected Gradient (SPG) iteration is used to calculate the new 
! point. If this new point is at the boundary, the algorithm stops. 
! Otherwise, a new iteration begins.
!
! When x is "far from" the boundary, trust-region radius is 
! determined and d is calculated using More-Sorensen algorithm to 
! solve the trust-region subproblem (which is to find a minimizer a 
! to a function quadratic model provided that the minimizer's 
! Euclidian-norm is smaller than a given delta). The new point is
! xp = x + d. 
!
! If xp lies outside the box of constraints, it is truncated on 
! the boundary. This new y on the boundary will be candidate to be a 
! solution. If function value at new xp is smaller than function 
! value at x, inner algorithm stops with xp. Otherwise, the 
! trust-region radius is decreased so that the new solution d' to 
! the trust-region subproblem makes x + d' be interior to the box.
! More-Sorensen algorithm is used to calculate d' too.
!
! If xp lies inside the box, sufficient decrease of objective 
! function is tested. If it is true, xp is accepted as a solution
! candidate. If xp in on the boundary, inner algorithm stops and 
! if it is interior, a new iteration begins. If sufficient decrease is 
! not obtained, trust-region radius is decreased and a new quadratic 
! model minimizer is calculated (as in a classical trust-region 
! algorithm for unconstrained minimization).
!
! If the user wants, after calculating the candidate solution xp, 
! extrapolation may be performed. For this, set extrpi to true, 
! where i is the step of the inner algorithm that can call 
! extrapolation procedure. If extrp4 is true, extrapolation will be
! tried after xp hit the boundary. And if extrp5 is true,
! extrapolation will be tried after xp is calculated by
! trust-region algorithm, when xp is interior and provides
! sufficient decrease of objective function.
!
! If gradient at current point is null, inner algorithm stops 
! declaring "first-order stationary point". If quadratic model 
! minimum is 0, inner algorithm stops declaring "second-order 
! stationary point".
!
! M. Andretta, E. G. Birgin and J. M. Martinez, ''Practical active-set
! Euclidian trust-region method with spectral projected gradients for
! bound-constrained minimization'', Optimization 54, pp. 305-325, 2005.
!
! On Entry
!  
! n        integer
!          dimension of full space
!
! nind     integer
!          dimension of reduced space
!
! x        double precision x(n)
!          initial point, interior to the current face
!
! l        double precision l(n)
!          lower bounds on x
!
! u        double precision u(n)
!          upper bounds on x
!
! m        integer
! lambda   double precision lambda(m)
! rho      double precision rho(m)
! equatn   logical equatn(m)
! linear   logical linear(m)
!          These five parameters are not used nor modified by 
!          BETRA and they are passed as arguments to the user-
!          defined subroutines evalal and evalnal to compute the 
!          objective function and its gradient, respectively. 
!          Clearly, in an Augmented Lagrangian context, if BETRA is 
!          being used to solve the bound-constrained subproblems, m 
!          would be the number of constraints, lambda the Lagrange 
!          multipliers approximation and rho the penalty parameters.
!          equatn is logical :: array that, for each constraint, 
!          indicates whether the constraint is an equality constraint
!          (.true.) or an inequality constraint (.false.). Finally,
!          linear is logical :: array that, for each constraint, 
!          indicates whether the constraint is a linear constraint
!          (.true.) or a nonlinear constraint (.false.)
!
! f        double precision
!          objective function value at x
!
! g        double precision g(n)
!          gradient at x
!
! trdelta  double precision 
!          trust-region radius
!
! newdelta double precision
!          trust-region radius is set to the maximum between newdelta 
!          and trdelta
! 
! mslamb   double precision
!          value that More-Sorensen algorithm calculates to find 
!          the trust-region subproblem solution (MEQB)
!
! epsg     double precision
!          allowed error for projected gradient norm
!
! xeucn    double precision
!          x Euclidian norm
!
! xsupn    double precision
!          x sup-norm
!
! gpeucn   double precision
!          projected-gradient Euclidian norm
!
! On Return
!  
! trdelta  double precision
!          updated trut-region radius
!
! newdelta double precision
!          when the trust-region radius trdelta is decreased so that 
!          the point x + d fit the current face, newdelta is set to 
!          the previous value of trdelta. Otherwise, it is set to 0
!
! mslamb   double precision
!          updated value for next iteration (see entry parameter)
!
! d        double precision d(nind)
!          direction computed such that xp = x + d
!
! rbdind   integer rbdind(n)
!          indices of variables that reached their bounds
! 
! rbdtype  character rbdtype(n)
!          if variable rbdind(i) reached its lower bound, 
!          rbdtype(i) = 'L'. If variable rbdind(i) reached its upper
!          bound, rbdtype(i) = 'U'. 
!
! xp       double precision xp(n)
!          solution candidate, with inform as described bellow
!
! fp       double precision
!          objective function value at xp
!
! gp       double precision gp(n)
!          gradient at xp
!
! triter   integer
!          number of trust-region iterations
!
! chcnt    integer
!          number of Cholesky decompositions
!
! memfail  logical
!          true iff linear solver failed because of lack of memory
!
! betinfo  integer
!          This output parameter tells what happened in this 
!          subroutine, according to the following conventions:
!  
!          0 = Sufficient decrease of the function;
!
!          1 = x hits the boundary;
!
!          2 = Unbounded objective function?
!
!          3 = Either the search direction or the step lenght is too
!              small;
!
!          4 = Trust-region radius too small;
!
!          5 = Undefined direction;
!
!          6 = x is a first-order stationary point close to the
!              boundary;
!
!          7 = x is a second-order stationary point;
!
! inform   0 = no error occurred;
!         <0 = error in some function, gradient or Hessian routine.
