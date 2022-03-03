! ******************************************************************
! ******************************************************************

subroutine cgm(n,nind,x,m,lambda,rho,equatn,linear,l,u,g,delta,eps, &
     maxit,d,iter,cginfo,inform)

  use modmachconst
  use modalgconst
  use modalgparam, only: setalgparam,precond,hptype
  use modouttyp

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: m,maxit,n,nind
  integer,      intent(inout) :: cginfo,inform,iter
  real(kind=8), intent(in)    :: delta,eps

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: g(nind),l(nind),lambda(m),rho(m), &
                                 u(nind),x(nind)
  real(kind=8), intent(out)   :: d(nind)

  ! This subroutine implements the Conjugate Gradients method for
  ! minimizing the quadratic approximation q(d) of L(x,lambda,rho)
  ! at x
  !
  ! q(d) = 1/2 d^T H d + g^T d,
  !
  ! where H is an approximation of the Hessian matrix of the
  ! Augmented Lagrangian and g is its gradient vector,
  !
  ! subject to || d || <= delta and l <= x + d <= u.
  !
  ! In the constraint ''|| d || <= delta'', the norm will be the
  ! Euclidian-norm if the input parameter trtype is equal to 0, and
  ! it will be the sup-norm if trtype is equal to 1.
  !
  ! The method returns an approximation d of the solution such that
  !
  ! (a) ||H d + g||_2 <= eps * ||g||_2,
  !
  ! (b) ||d|| = delta or x + d is in the boundary of the box, or
  !
  ! (c) ( p such that p^t H p = 0 ) and ( d = - amax g if such p was
  !     found during the first CG iteration or the current point d
  !     of CG if such p was found in any other iteration ).
  !
  ! cginfo integer
  !
  ! termination parameter:
  !
  ! 0 = convergence with ||H d + g||_2 <= eps * ||g||_2;
  !
  ! 1 = convergence to the boundary of ||d|| <= delta;
  !
  ! 2 = convergence to the boundary of l <= x + d <= u;
  !
  ! 3 = stopping with d = dk  such that <gk,dk> <= - theta
  !     ||gk||_2 ||dk||_2 and <gk,d_{k+1}> > - theta
  !     ||gk||_2 ||d_{k+1}||_2;
  !
  ! 4 = not enough progress of the quadratic model during
  !     maxitnqmp iterations, i.e., during maxitnqmp
  !     iterations | q - qprev | <= macheps * max( | q |, 1 )
  !
  ! 6 = very similar consecutive iterates, for two
  !     consecutive iterates x1 and x2 we have that
  !
  !     | x2(i) - x1(i) | <= macheps * max ( | x1(i) |, 1 )
  !
  !     for all i.
  !
  ! 7 = stopping with p such that p^T H p = 0 and g^T p = 0;
  !
  ! 8 = too many iterations;

  ! LOCAL SCALARS
  character(len=4) :: prectmp
  logical          :: goth,gotp,negcur,restarted,samep
  integer          :: i,itertmp,itnqmp
  real(kind=8)     :: alpha,abox,amax,bestprog,bbeta,currprog,   &
                      dnorm2,gnorm2,gtd,gtp,pnorm2,ptd,pthp,ptr, &
                      q,qprev,rnorm2,ztrprev,ztr,znorm2

  ! LOCAL ARRAYS
  real(kind=8) :: p(n),hp(n),r(n),z(n)

  ! ==================================================================
  ! Initialization
  ! ==================================================================

  restarted = .false.

001 continue

  goth = .false.
  gotp = .false.

  ! gnorm2 = norm2s(nind,g)
  gnorm2 = 0.0d0
  do i = 1,nind
     gnorm2 = gnorm2 + g(i) ** 2
  end do

  iter     =      0
  itnqmp   =      0
  qprev    = bignum
  bestprog =  0.0d0

  do i = 1,nind
     d(i) = 0.0d0
     r(i) =  g(i)
  end do

  q        =  0.0d0
  gtd      =  0.0d0
  dnorm2   =  0.0d0
  rnorm2   = gnorm2

  ztr      =  0.0d0

  ! ==================================================================
  ! Print initial information
  ! ==================================================================

  if ( iprintinn .ge. 5 ) then
     write(*, 980) maxit,eps,delta,precond,hptype
     write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q

     write(10,980) maxit,eps,delta,precond,hptype
     write(10,984) iter,sqrt(rnorm2),sqrt(dnorm2),q
  end if

  ! ==================================================================
  ! Main loop
  ! ==================================================================

100 continue

  ! ==================================================================
  ! Test stopping criteria
  ! ==================================================================

  ! if ||r||_2 = ||H d + g||_2 <= eps * ||g||_2 then stop

  if ( iter .ne. 0 .and. &
       ( ( rnorm2 .le. eps ** 2 * gnorm2 .and. iter .ge. 4 ) .or. &
       ( rnorm2 .le. macheps ) ) ) then

     cginfo = 0

     if ( iprintinn .ge. 5 ) then
        write(*, 990)
        write(10,990)
     end if

     go to 500

  end if

  ! if the maximum number of iterations was achieved then stop

  if ( iter .ge. max(4, maxit) ) then

     cginfo = 8

     if ( iprintinn .ge. 5 ) then
        write(*, 998)
        write(10,998)
     end if

     go to 500

  end if

  ! ==================================================================
  ! Preconditioner
  ! ==================================================================

  if ( precond .eq. 'NONE' ) then

     do i = 1,nind
        z(i) = r(i)
     end do

     ztrprev = ztr
     ztr     = rnorm2
     znorm2  = rnorm2

  else if ( precond .eq. 'QNGN' ) then

     call capplyhpre(nind,m,rho,equatn,gotp,r,z)

     ztrprev = ztr

     ztr = 0.0d0
     do i = 1,nind
        ztr = ztr + z(i) * r(i)
     end do

     ! znorm2 = norm2s(nind,z)
     znorm2 = 0.0d0
     do i = 1,nind
        znorm2 = znorm2 + z(i) ** 2
     end do

  end if

  ! ==================================================================
  ! Compute direction
  ! ==================================================================

  if ( iter .eq. 0 ) then

     do i = 1,nind
        p(i) = - z(i)
     end do

     ptr    = - ztr
     pnorm2 =   znorm2

  else

     bbeta = ztr / ztrprev

     do i = 1,nind
        p(i) = - z(i) + bbeta * p(i)
     end do

     if ( precond .eq. 'NONE' ) then

        pnorm2 = rnorm2 - 2.0d0 * bbeta * ( ptr + alpha * pthp ) &
             + bbeta ** 2 * pnorm2
        ptr = - rnorm2 + bbeta * ( ptr + alpha * pthp )

     else if ( precond .eq. 'QNGN' ) then

        ptr = 0.0d0
        pnorm2 = 0.0d0
        do i = 1,nind
           ptr = ptr + p(i) * r(i)
           pnorm2 = pnorm2 + p(i) ** 2
        end do

     end if

  end if

  ! Force p to be a descent direction of q(d), i.e.,
  ! <\nabla q(d), p> = <H d + g, p> = <r, p> \le 0.

  if ( ptr .gt. 0.0d0 ) then

     do i = 1,nind
        p(i) = - p(i)
     end do

     ptr = - ptr

  end if

  ! ==================================================================
  ! Compute p^T H p
  ! ==================================================================

  ! hp = H p

  call calchalp(nind,x,m,lambda,rho,equatn,linear,p,hp,goth,inform)
  if ( inform .ne. 0 ) return

  ! Compute p^T hp

  pthp = 0.0d0
  do i = 1,nind
     pthp = pthp + p(i) * hp(i)
  end do

  ! ==================================================================
  ! Compute maximum step
  ! ==================================================================

  amax = bignum

  do i = 1,nind
     if ( p(i) .gt. 0.0d0 ) then
        amax = min( amax,  (   delta - d(i) ) / p(i) )
     else if ( p(i) .lt. 0.0d0 ) then
        amax = min( amax,  ( - delta - d(i) ) / p(i) )
     end if
  end do

  ! ==================================================================
  ! Compute the step
  ! ==================================================================

  negcur = .false.

  ! If p^T H p > 0 then take the conjugate gradients step

  if ( pthp .gt. 0.0d0 ) then

     alpha = min( amax, ztr / pthp )

  ! Else, if we are at iteration zero then take the maximum
  ! positive step in the minus gradient direction

  else if ( iter .eq. 0 ) then

     alpha = amax

     negcur = .true.

  ! Otherwise, stop at the current iterate

  else

     cginfo = 7

     if ( iprintinn .ge. 5 ) then
        write(*, 997)
        write(10,997)
     end if

     go to 500

  end if

  ! ==================================================================
  ! Test the angle condition
  ! ==================================================================

  ptd = 0.0d0
  gtp = 0.0d0
  do i = 1,nind
     ptd = ptd + p(i) * d(i)
     gtp = gtp + g(i) * p(i)
  end do

  ! These are gtd and dnorm2 for the new direction d which was not
  ! computed yet.

  gtd = gtd + alpha * gtp
  dnorm2 = dnorm2 + alpha ** 2 * pnorm2 + 2.0d0 * alpha * ptd

  if ( gtd .gt. 0.0d0 .or. &
       gtd ** 2 .lt. theta ** 2 * gnorm2 * dnorm2 ) then

     if ( precond .ne. 'NONE' .and. iter .eq. 0 ) then

        if ( iprintinn .ge. 5 ) then
           write(*, 986)
           write(10,986)
        end if

        restarted = .true.
        itertmp   = iter
        prectmp   = precond
        call setalgparam(val_precond   = 'NONE')
        go to 001

     end if

     cginfo = 3

     if ( iprintinn .ge. 5 ) then
        write(*, 993)
        write(10,993)
     end if

     go to 500

  end if

  ! ==================================================================
  ! Compute the quadratic model functional value at the new point
  ! ==================================================================

  qprev = q

  q = q + 0.5d0 * alpha ** 2 * pthp + alpha * ptr

  ! ==================================================================
  ! Compute new d
  ! ==================================================================

  do i = 1,nind
     d(i) = d(i) + alpha * p(i)
  end do

  ! ==================================================================
  ! Compute the residual r = H d + g
  ! ==================================================================

  do i = 1,nind
     r(i) = r(i) + alpha * hp(i)
  end do

  ! rnorm2 = norm2s(nind,r)
  rnorm2 = 0.0d0
  do i = 1,nind
     rnorm2 = rnorm2 + r(i) ** 2
  end do

  ! ==================================================================
  ! Increment number of iterations
  ! ==================================================================

  iter = iter + 1

  ! ==================================================================
  ! Print information of this iteration
  ! ==================================================================

  if ( iprintinn .ge. 5 ) then
     write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q
     write(10,984) iter,sqrt(rnorm2),sqrt(dnorm2),q
  end if

  ! ==================================================================
  ! Test other stopping criteria
  ! ==================================================================

  ! Boundary of the "trust region"

  if ( alpha .eq. amax ) then

     cginfo = 1

     if ( iprintinn .ge. 5 ) then
        if ( negcur ) then
           write(*, 987)
           write(10,987)
        end if

        write(*, 991)
        write(10,991)
     end if

     go to 500

  end if

  ! Small useful proportion

  abox = bignum

  do i = 1,nind
     if ( d(i) .gt. 0.0d0 ) then
        abox = min( abox, ( u(i) - x(i) ) / d(i) )
     else if ( d(i) .lt. 0.0d0 ) then
        abox = min( abox, ( l(i) - x(i) ) / d(i) )
     end if
  end do

  if ( abox .le. 0.1d0 ) then

     cginfo = 5

     if ( iprintinn .ge. 5 ) then
        write(* ,995)
        write(10,995)
     end if

     go to 500

  end if

  ! Two consecutive iterates are too much close

  samep = .true.
  do i = 1,nind
     if ( abs( alpha * p(i) ) .gt. &
          macheps * max( 1.0d0, abs( d(i) ) ) ) then
        samep = .false.
     end if
  end do

  if ( samep ) then

     cginfo = 6

     if ( iprintinn .ge. 5 ) then
        write(*, 996)
        write(10,996)
     end if

     go to 500

  end if

  ! Many iterations without good progress of the quadratic model

  currprog = qprev - q
  bestprog = max( currprog, bestprog )

  if ( currprog .le. epsnqmp * bestprog ) then

     itnqmp = itnqmp + 1

     if ( itnqmp .ge. maxcgitnp ) then
        cginfo = 4

        if ( iprintinn .ge. 5 ) then
           write(*, 994)
           write(10,994)
        end if

        go to 500
     endif

  else
     itnqmp = 0
  endif

  ! ==================================================================
  ! Iterate
  ! ==================================================================

  go to 100

  ! ==================================================================
  ! End of main loop
  ! ==================================================================

  ! ==================================================================
  ! Return
  ! ==================================================================

500 continue

  ! Print final information

  if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
     write(*, 985) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
     write(10,985) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
  end if

  if ( restarted ) then
     iter = iter + itertmp
     call setalgparam(val_precond = prectmp)
  end if

  return

! Non-executable statements

 980  format(/,5X,'Conjugate Gradients (maxit = ',I7,',',1X,'eps = ', &
               1P,D7.1,',',1X,'delta = ',1P,D7.1,')',                 &
             /,5X,'(Preconditioner: ',A4,',',1X,                      &
                  'Hessian-vector product type: ',A6,')')
 984  format(  5X,'CG iter = ',I7,' rnorm = ',1P,D10.4,' dnorm = ', &
                  1P,D10.4,' q = ',1P,D11.4)
 985  format(/,5X,'Truncated Newton direction (first ',I7, &
                  ' components): ',/,1(5X,6(1X,1P,D11.4)))
 986  format(  5X,'The first CG-PREC iterate did not satisfy the ', &
                  'angle condition.',                               &
             /,5X,'CG will be restarted without preconditioner)')
 987  format(  5X,'p such that p^T H p = 0 was found. ', &
                  'Maximum step was taken.')

 990  format(  5X,'Flag of CG: Convergence with small residual.')
 991  format(  5X,'Flag of CG: Convergence to the trust region ', &
                  'boundary.')
!992  format(  5X,'Flag of CG: Convergence to the box boundary.')
 993  format(  5X,'Flag of CG: The next CG iterate will not satisfy ', &
                  'the angle condition.')
 994  format(  5X,'Flag of CG: Not enough progress in the quadratic ', &
                  'model.')
 995  format(  5X,'Flag of CG: The maximum step to remain within the ', &
                  'box is smaller than 0.1.')
 996  format(  5X,'Flag of CG: Very near consecutive iterates.')
 997  format(  5X,'Flag of CG: p such that p^T H p = 0 was found.')
 998  format(  5X,'Flag of CG: Maximum number of CG iterations ', &
                  'reached.')

end subroutine cgm
