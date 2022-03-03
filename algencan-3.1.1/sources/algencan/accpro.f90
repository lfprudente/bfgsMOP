! **************************************************************************************************
! **************************************************************************************************

subroutine seco(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,efstain,eostain,perturbx0, &
     nobacktrack,f,cnorm,bdsnorm,nlnorm,fu,cunorm,outiter,totiter,best_f,best_cnorm,best_bdsnorm, &
     best_nlnorm,best_fu,best_cunorm,best_outit,best_x,best_lambda,secoinfo,inform)

  use modmachconst, only: macheps12
  use modouttyp, only: iprintout,iprintctl,nprint
  use problvls, only: ssetp,sevalg,sevalgjac,sevaljac,sevalobjc
  use probgiven, only: gjaccoded,jcnnzlim
    
  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(in) :: nobacktrack,perturbx0
  integer, intent(in) :: m,n
  integer, intent(out) :: best_outit,outiter,secoinfo,totiter
  integer, intent(inout) :: inform
  real(kind=8), intent(in) :: efstain,eostain,epsfeas,epsopt
  real(kind=8), intent(out) :: bdsnorm,cnorm,cunorm,f,fu,nlnorm,best_f,best_fu,best_cnorm, &
       best_bdsnorm,best_cunorm,best_nlnorm

  ! ARRAY ARGUMENTS
  logical, intent(in) :: equatn(m),linear(m)
  real(kind=8), intent(inout) :: l(n),u(n)
  real(kind=8), intent(inout) :: lambda(m),x(n)
  real(kind=8), intent(out) :: best_lambda(m),best_x(n)

  ! LOCAL SCALARS
  integer :: accinfo,i,iter,j,jcnnz,nlind,nuind,bstk_iter
  real(kind=8) :: cmplnrm,cmplnrmprev,cnorm2,cunorm2,epsfek,epsopk,fpenbds,rho,bstk_f,bstk_cnorm, &
       bstk_bdsnorm,bstk_nlnorm,bstk_fu,bstk_cunorm,seed
  
  ! LOCAL ARRAYS
  logical :: ind(m)
  integer :: jcfun(jcnnzlim),jclen(m),jcsta(m),jcvar(jcnnzlim),lind(n),uind(n)
  real(kind=8) :: c(m),cu(m),g(n),gpenbds(n),hpenbds(n),jcval(jcnnzlim),lmu(n),lmubar(n), &
       nl(n),tmp(n),umu(n),umubar(n),bstk_x(n),bstk_lambda(m)

  ! EXTERNAL FUNCTIONS
  real(kind=8) :: drand

  ! LOCALS THAT MAY BE PARAMETERS
  integer :: maxoutit
  real(kind=8) :: r,tau,gamma,mumax

  ! ================================================================================================
  ! Presentation
  ! ================================================================================================

  if ( iprintout .ge. 1 ) then
     write(* ,8000)
     write(10,8000)
  end if

  ! ================================================================================================
  ! Add a small (1%) perturbation to x0
  ! ================================================================================================

  if ( perturbx0 ) then
     seed = 123456.0d0
     do i = 1,n
        x(i) = x(i) + 0.01d0 * ( 2.0d0 * drand(seed) - 1.0d0 ) * max( abs( x(i) ), macheps12 )
     end do
  end if

  ! ================================================================================================
  ! Set bounds' indices arrays
  ! ================================================================================================

  ! "Ifs" below may be commented or not depending on whether
  ! artificial bound-constraints are being considered.
  
  nlind = 0
  do i = 1,n
     if ( l(i) .gt. - 1.0d+20 ) then
        nlind = nlind + 1
        lind(nlind) = i
     end if
  end do
  
  nuind = 0
  do i = 1,n
     if ( u(i) .lt. 1.0d+20 ) then
        nuind = nuind + 1
        uind(nuind) = i
     end if
  end do

  ! ================================================================================================
  ! If there are no bounds, call the subproblems' solver only once
  ! ================================================================================================

  if ( nlind .eq. 0 .and. nuind .eq. 0 ) then
     
     call accpro(n,x,l,u,m,lambda,equatn,linear,nlind,lind,lmubar,     &
          nuind,uind,umubar,rho,epsfeas,epsopt,epsfeas,epsopt,efstain, &
          eostain,nobacktrack,f,cnorm,bdsnorm,nlnorm,fu,cunorm,iter,   &
          bstk_f,bstk_cnorm,bstk_bdsnorm,bstk_nlnorm,bstk_fu,          &
          bstk_cunorm,bstk_iter,bstk_x,bstk_lambda,accinfo,inform)

     best_outit   = bstk_iter

     best_f       = bstk_f
     best_cnorm   = bstk_cnorm
     best_bdsnorm = bstk_bdsnorm
     best_nlnorm  = bstk_nlnorm
     best_fu      = bstk_fu
     best_cunorm  = bstk_cunorm
     
     best_x(1:n)      = bstk_x(1:n)
     best_lambda(1:m) = bstk_lambda(1:m)

     secoinfo = accinfo
     
     outiter  = 0
     totiter  = iter
     
     return

  end if

  ! ================================================================================================
  ! Initialization
  ! ================================================================================================

  ! Constants that might be parameters
  
  r = 0.5d0
  tau = 0.5d0
  gamma = 10.0d0
  mumax = 1.0d+20
  maxoutit = 50

  ! Iterations counters
  
  outiter = 0
  totiter = 0
  
  ! Set bounds' multipliers initial value
  
  lmu(1:nlind) = 0.0d0
  umu(1:nuind) = 0.0d0

  where ( lmu(1:nlind) .le. mumax )
     lmubar(1:nlind) = lmu(1:nlind)
  else where
     lmubar(1:nlind) = 0.0d0
  end where

  where ( umu(1:nuind) .le. mumax )
     umubar(1:nuind) = umu(1:nuind)
  else where
     umubar(1:nuind) = 0.0d0
  end where
  
  ! ================================================================================================
  ! Main loop
  ! ================================================================================================

10 continue

  ! ================================================================================================
  ! Compute objective function and constraints
  ! ================================================================================================

  call ssetp(n,x,inform)
  if ( inform .ne. 0 ) return

  call sevalobjc(n,x,f,fu,m,c,cu,inform,ignscl=.false.)
  if ( inform .ne. 0 ) return

  if ( m .eq. 0 ) then
     cnorm   = 0.0d0
     cunorm  = 0.0d0
     cnorm2  = 0.0d0
     cunorm2 = 0.0d0
  else
     cnorm   = maxval( abs( c(1:m) ) )
     cunorm  = maxval( abs( cu(1:m) ) )
     cnorm2  = 0.5d0 * sum( c(1:m)  ** 2 )
     cunorm2 = 0.5d0 * sum( cu(1:m) ** 2 )
  end if

  bdsnorm = 0.0d0
  bdsnorm = max( bdsnorm, maxval( l(lind(1:nlind)) - x(lind(1:nlind)) ) )
  bdsnorm = max( bdsnorm, maxval( x(uind(1:nuind)) - u(uind(1:nuind)) ) )

  ! ================================================================================================
  ! Gradient of the objective function and Jacobian of the constraints
  ! ================================================================================================

  if ( gjaccoded ) then

     call sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform,ignscl=.false.)
     if ( inform .ne. 0 ) return

  else ! if ( gcoded .and. ( jaccoded .or. m .eq. 0 ) ) then

     call sevalg(n,x,g,inform,ignscl=.false.)
     if ( inform .ne. 0 ) return

     if ( m .eq. 0 ) then
        jcnnz = 0
        
     else
        ind(1:m) = .true.
        call sevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,jcnnzlim,inform,ignscl=.false.)
        if ( inform .ne. 0 ) return

        jcnnz = sum( jclen(1:m) )
     
        do j = 1,m
           jcfun(jcsta(j):jcsta(j)+jclen(j)-1) = j
        end do
     end if
     
  end if

  ! ================================================================================================
  ! Compute gradient of the Lagrangian
  ! ================================================================================================

  nl(1:n) = g(1:n)

  do i = 1,jcnnz
     nl(jcvar(i)) = nl(jcvar(i)) + lambda(jcfun(i)) * jcval(i)
  end do

  ! Here we compute the norm of the gradient of the Lagrangian
  ! considering the most convenient choice of multipliers (for the
  ! bound constraints) that are non-negative with zero-tolerance and
  ! guaranteedly satisfy the complementarity constraint (as measured
  ! with the min function) with tolerance epsfeas.
  
  nlnorm = 0.0d0
  do i = 1,n
     if ( l(i) + epsfeas .lt. x(i) .and. x(i) .lt. u(i) - epsfeas ) then
        ! This is x(i) - l(i) > epsfeas and u(i) - x(i) > epsfeas,
        ! i.e. x(i) is in the interior of the box, written in a way
        ! that aims to avoid catastrophic cancelations.

        ! If we think in the bound constraints as gl(x) = l - x <= 0
        ! and gu(x) = x - u <= 0, then in this case we have -
        ! [gl(x)]_i > epsfeas and - [gu(x)]_i > epsfeas.

        ! If we aim to satisfy the complementarity constraints given
        ! by min{ - [gl(x)]_i, lmu_i } <= epsfeas and min{ -
        ! [gu(x)]_i, umu_i } <= epsfeas, for non-negative multipliers
        ! lmu_i and umu_i, then we must have 0 <= lmu_i <= epsfeas and
        ! 0 <= umu_i <= epsfeas.

        ! O.K., this is "using tolerances as much as we can", and we
        ! will make use of this chance whenever possible. Note that,
        ! at least, we are considering multipliers that satisfy
        ! non-negativity with zero-tolerance.

        ! Now, we have freedom to choose 0 <= lmu_i <= epsfeas and 0
        ! <= umu_i <= epsfeas in order to have the quantity nl(i) -
        ! lmu_i + umu_i as small as possible.

        ! If nl(i) = 0, we choose lmu_i = umu_i = 0. If nl(i) > 0, we
        ! choose lmu_i = min{ nl(i), epsfeas } and umu_i = 0. If nl(i)
        ! < 0, we choose lmu_i = 0 and umu_i = min{ - nl(i), epsfeas
        ! }.

        ! All choices above, simply translate into computing the
        ! contribution of the i-th component of the gradient to its
        ! norm as follows:

        if ( nl(i) .gt. 0.0d0 ) then
           nlnorm = max( nlnorm,   max( 0.0d0, nl(i) - epsfeas ) )
        else if ( nlnorm .lt. 0.0d0 ) then
           nlnorm = max( nlnorm, - min( 0.0d0, nl(i) + epsfeas ) )
        end if

     else if ( l(i) + epsfeas .ge. x(i) .and. x(i) .lt. u(i) - epsfeas ) then

        ! This is x near the lower and far from the upper bound. In
        ! this case we have - [gl(x)]_i <= epsfeas and - [gu(x)]_i >
        ! epsfeas. Therefore, we are free to choose any lmu_i >= 0 and
        ! 0 <= umu_i <= epsfeas.

        ! If nl(i) = 0, we choose lmu_i = umu_i = 0. If nl(i) > 0, we
        ! choose lmu_i = nl(i) and umu_i = 0. If nl(i) < 0, we choose
        ! lmu_i = 0 and umu_i = min{ - nl(i), epsfeas }.
        
        ! All choices above, simply translate into computing the
        ! contribution of the i-th component of the gradient to its
        ! norm as follows:

        if ( nl(i) .lt. 0.0d0 ) then
           nlnorm = max( nlnorm, - min( 0.0d0, nl(i) + epsfeas ) )
        end if

     else if ( l(i) + epsfeas .lt. x(i) .and. x(i) .ge. u(i) - epsfeas ) then

        ! This is x near the upper and far from the lower bound. In
        ! this case we have - [gl(x)]_i > epsfeas and - [gu(x)]_i <=
        ! epsfeas. Therefore, we are free to choose any 0 <= lmu_i <=
        ! epsfeas and umu_i >= 0.

        ! If nl(i) = 0, we choose lmu_i = umu_i = 0. If nl(i) > 0, we
        ! choose lmu_i = min{ nl(i), epsfeas } and umu_i = 0. If nl(i)
        ! < 0, we choose lmu_i = 0 and umu_i = - nl(i).
        
        ! All choices above, simply translate into computing the
        ! contribution of the i-th component of the gradient to its
        ! norm as follows:
        
        if ( nl(i) .gt. 0.0d0 ) then
           nlnorm = max( nlnorm, max( 0.0d0, nl(i) - epsfeas ) )
        end if
        
     else ! if ( l(i) + epsfeas .ge. x(i) .and. x(i) .ge. u(i) - epsfeas ) then

        ! This is x near both, the lower and the upper bound. In this
        ! case we have - [gl(x)]_i <= epsfeas and - [gu(x)]_i <=
        ! epsfeas. Therefore, we are free to choose any lmu_i >= 0 and
        ! umu_i >= 0.

        ! If nl(i) = 0, we choose lmu_i = umu_i = 0. If nl(i) > 0, we
        ! choose lmu_i = nl(i) and umu_i = 0. If nl(i) < 0, we choose
        ! lmu_i = 0 and umu_i = - nl(i).
        
        ! All choices above, simply translate into computing the
        ! contribution of the i-th component of the gradient to its
        ! norm as follows:

        ! There is no contribution in this case.
     end if
  end do
  
!!$  tmp(1:n) = nl(1:n)
!!$  tmp(lind(1:nlind)) = tmp(lind(1:nlind)) - lmu(1:nlind)
!!$  tmp(uind(1:nuind)) = tmp(uind(1:nuind)) + umu(1:nuind)
!!$
!!$  nlnorm = maxval( abs( tmp(1:n) ) )

!!$  nlnorm = maxval( abs( min( max( l(1:n), x(1:n) - nl(1:n) ), u(1:n) ) - x(1:n) ) )

  ! ================================================================================================
  ! Save best point information
  ! ================================================================================================

  if ( outiter .eq. 0 ) then
     best_outit   = outiter

     best_f       = f
     best_fu      = fu
     best_cnorm   = cnorm
     best_cunorm  = cunorm
     best_bdsnorm = bdsnorm
     best_nlnorm  = nlnorm

     best_x(1:n)      = x(1:n)
     best_lambda(1:m) = lambda(1:m)

  else

     if ( ( max( best_bdsnorm, best_cunorm ) .gt. epsfeas .and. &
            ( max( bstk_bdsnorm, bstk_cunorm ) .lt. max( best_bdsnorm, best_cunorm ) .or. &
            ( max( bstk_bdsnorm, bstk_cunorm ) .le. max( best_bdsnorm, best_cunorm ) .and. &
              bstk_fu .lt. best_fu ) .or. &
            ( max( bstk_bdsnorm, bstk_cunorm ) .le. max( best_bdsnorm, best_cunorm ) .and. &
              bstk_fu .le. best_fu .and. bstk_nlnorm .lt. best_nlnorm ) ) ) .or. &
          ( max( best_bdsnorm, best_cunorm ) .le. epsfeas .and. max( bstk_bdsnorm, bstk_cunorm ) .le. epsfeas .and. &
            ( bstk_fu .lt. best_fu .or. ( bstk_fu .le. best_fu .and. bstk_nlnorm .lt. best_nlnorm ) ) ) ) then
        best_outit   = outiter

        best_f       = bstk_f
        best_fu      = bstk_fu
        best_cnorm   = bstk_cnorm
        best_cunorm  = bstk_cunorm
        best_bdsnorm = bstk_bdsnorm
        best_nlnorm  = bstk_nlnorm

        best_x(1:n)      = bstk_x(1:n)
        best_lambda(1:m) = bstk_lambda(1:m)
     end if
  end if

  ! ================================================================================================
  ! Print information of current iterate
  ! ================================================================================================

  if ( iprintout .ge. 1 ) then
     if ( outiter .eq. 0 ) then
        write(* ,8010) outiter,fu,f,cunorm,cnorm,cunorm2,cnorm2,bdsnorm,nlnorm
        write(10,8010) outiter,fu,f,cunorm,cnorm,cunorm2,cnorm2,bdsnorm,nlnorm
     else
        write(* ,8011) rho,epsfek,epsopk,outiter,fu,f,cunorm,cnorm,cunorm2,cnorm2,bdsnorm,nlnorm
        write(10,8011) rho,epsfek,epsopk,outiter,fu,f,cunorm,cnorm,cunorm2,cnorm2,bdsnorm,nlnorm
     end if
     write(* ,8015) nprint,maxval(abs(x(1:n))),(x(i),i=1,nprint)
     write(10,8015) nprint,maxval(abs(x(1:n))),(x(i),i=1,nprint)
  end if
  
  if ( iprintctl(5) ) then
     open(20,file='seco-interrupted-tabline.out')
     write(20,9100) fu,cunorm,bdsnorm,nlnorm,f,cnorm,outiter,totiter,best_fu,best_cunorm, &
          best_bdsnorm,best_nlnorm,best_f,best_cnorm,best_outit,9,inform,999.99d0
     close(20)
  end if

  ! ================================================================================================
  ! Check stopping criteria
  ! ================================================================================================

  if ( max( cunorm, bdsnorm ) .le. epsfeas .and. nlnorm .le. epsopt ) then
     secoinfo = 0
     
     if ( iprintout .ge. 1 ) then
        write(* ,9000)
        write(10,9000)
     end if
        
     go to 500
  end if

  if ( outiter .ge. maxoutit ) then
     secoinfo = 1
     
     if ( iprintout .ge. 1 ) then
        write(* ,9010)
        write(10,9010)
     end if
        
     go to 500
  end if

  ! ================================================================================================
  ! Perform a new outer iteration
  ! ================================================================================================

  outiter = outiter + 1
  
  ! ================================================================================================
  ! Update (or set) penalty parameter
  ! ================================================================================================

  if ( outiter .eq. 1 ) then
     rho = 10.0d0 * min( max( 1.0d0, max( 1.0d0, abs( f ) ) / max( 1.0d0, fpenbds ) ), 1.0d+08 )

  else if ( cmplnrm .gt. tau * cmplnrmprev ) then 
     rho = gamma * rho
  end if

  ! ================================================================================================
  ! Set optimality requeriment for the subproblem
  ! ================================================================================================

  if ( outiter .eq. 1 ) then
     epsfek = sqrt( epsfeas )
     epsopk = sqrt( epsopt  )
     
  else if ( cunorm .le. sqrt( epsfeas ) .and. nlnorm .le. sqrt( epsopt ) ) then
     epsfek = max( min( r * cunorm, 0.1d0 * epsfek ), epsfeas )
     epsopk = max( min( r * nlnorm, 0.1d0 * epsopk ), epsopt  )
  end if

  ! ================================================================================================
  ! Call inner solver
  ! ================================================================================================

  call accpro(n,x,l,u,m,lambda,equatn,linear,nlind,lind,lmubar,nuind,uind,umubar,rho,   &
       epsfek,epsopk,epsfeas,epsopt,efstain,eostain,nobacktrack,f,cnorm,bdsnorm,nlnorm, &
       fu,cunorm,iter,bstk_f,bstk_cnorm,bstk_bdsnorm,bstk_nlnorm,bstk_fu,bstk_cunorm,   &
       bstk_iter,bstk_x,bstk_lambda,accinfo,inform)
       
  totiter = totiter + iter
  
  if ( inform .ne. 0 ) return

  ! ================================================================================================
  ! Compute (bounds) feasibility/complementarity violation
  ! ================================================================================================

  cmplnrmprev = cmplnrm
  
  cmplnrm = 0.0d0
  cmplnrm = max( cmplnrm, maxval( abs( min( - ( l(lind(1:nlind)) - x (lind(1:nlind)) ), lmubar(1:nlind) / rho ) ) ) )
  cmplnrm = max( cmplnrm, maxval( abs( min( - ( x(uind(1:nuind)) - u (uind(1:nuind)) ), umubar(1:nuind) / rho ) ) ) )

  ! ================================================================================================
  ! Update multipliers
  ! ================================================================================================

  lmu(1:nlind) = max( 0.0d0, lmubar(1:nlind) + rho * ( l(lind(1:nlind)) - x(lind(1:nlind)) ) )
  umu(1:nuind) = max( 0.0d0, umubar(1:nuind) + rho * ( x(uind(1:nuind)) - u(uind(1:nuind)) ) )

  where ( lmu(1:nlind) .le. mumax )
     lmubar(1:nlind) = lmu(1:nlind)
  else where
     lmubar(1:nlind) = 0.0d0
  end where

  where ( umu(1:nuind) .le. mumax )
     umubar(1:nuind) = umu(1:nuind)
  else where
     umubar(1:nuind) = 0.0d0
  end where

  ! ================================================================================================
  ! Iterate
  ! ================================================================================================
  
  go to 10
  
  ! ================================================================================================
  ! End of main loop
  ! ================================================================================================

500 continue

  ! ================================================================================================
  ! Non-executable statements
  ! ================================================================================================

 8000 format(/,' Sequentialy equality-constrained optimization (SECO) in action!')
 8010 format(/,' SECO iteration',53X,                                                   ' = ',  5X,I6, &
             /,' Objective function value                  (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' Sup-norm of constraints                   (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' 1/2 squared two norm of constraints       (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' Sup-norm of the bounds violation                                  ',1X,' = ',1PD11.4, &
             /,' Sup norm of gradient of the Lagrangian                            ',1X,' = ',1PD11.4)
 8011 format(/,' SECO iteration   (rho = ',1PD7.1,' epsfek = ',1PD7.1,' epsopk = ',1PD7.1,')',1X,' = ',5X,I6, &
             /,' Objective function value                  (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' Sup-norm of constraints                   (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' 1/2 squared two norm of constraints       (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' Sup-norm of the bounds violation                                  ',1X,' = ',1PD11.4, &
             /,' Sup norm of gradient of the Lagrangian                            ',1X,' = ',1PD11.4)
 8015 format(/,' Current point first ',I1,' components (sup-norm: ',1PD11.4,'): ',/,6(1X,1P,D11.4))
 9000 format(/,' Flag of SECO = Problem solved!')
 9010 format(/,' Flag of SECO = Maximum of iterations.')

 9100 format(1X,1P,D24.16,3(1X,1P,D7.1),1X,1P,D24.16,1X,1P,D7.1,2(1X,I12), &
             1X,1P,D24.16,3(1X,1P,D7.1),1X,1P,D24.16,1X,1P,D7.1,1X,I12,1X,I1,1X,I3,0P,F8.2)

end subroutine seco

! **************************************************************************************************
! **************************************************************************************************

subroutine accpro(n,x,l,u,m,lambda,equatn,linear,nlind,lind,lmu,nuind,uind,umu,rho,epsfek,epsopk, &
     epsfeas,epsopt,efstain,eostain,nobacktrack,f,cnorm,bdsnorm,nlnorm,fu,cunorm,iter,best_f,     &
     best_cnorm,best_bdsnorm,best_nlnorm,best_fu,best_cunorm,best_iter,best_x,best_lambda,        &
     accinfo,inform)

  use modmachconst, only: macheps,macheps12,bignum
  use modalgparam, only: lsssubACC,sclsubACC
  use modouttyp, only: iprintinn,iprintctl,nprint
  use problvls, only: scale,ssetp,sevalg,sevalgjac,sevaljac,sevalhl,sevalobjc
  use probgiven, only: gjaccoded,hnnzlim,jcnnzlim
  use lsslvr, only: lssana,lssend,lssfac,lssset,lsssol
    
  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(in) :: nobacktrack
  integer, intent(in) :: m,n,nlind,nuind
  integer, intent(out) :: accinfo,best_iter,iter
  integer, intent(inout) :: inform
  real(kind=8), intent(in) :: efstain,eostain,epsfeas,epsfek,epsopt,epsopk,rho
  real(kind=8), intent(out) :: bdsnorm,cnorm,cunorm,f,fu,nlnorm,best_f,best_fu,best_cnorm, &
       best_cunorm,best_bdsnorm,best_nlnorm

  ! ARRAY ARGUMENTS
  logical, intent(in) :: equatn(m),linear(m)
  integer, intent(in) :: lind(nlind),uind(nuind)
  real(kind=8), intent(in) :: l(n),lmu(nlind),u(n),umu(nuind)
  real(kind=8), intent(inout) :: lambda(m),x(n)
  real(kind=8), intent(out) :: best_x(n),best_lambda(m)

  ! LOCAL SCALARS
  character :: ittype,star
  logical :: avoidbcktr,avoidstepctrl,bcktrfail,tolstp,xcutted
  integer :: allocerr,hnnz,i,j,jcnnz,k,lim,lssinfo,nneigv,nrank,nsys,pind,pnopro,unnz
  real(kind=8) :: auglag,auglagt,cnorm2,cnorm2t,cnormt,cunorm2,cunormt,cunorm2t,bdsnormt,dlnorm, &
       dxnorm,epsadd1,epsadd2,fpenbds,ft,ftpenbds,fut,lagra,lagrat,lnorm,pval,step,xnorm

  ! LOCAL ARRAYS
  logical :: ind(m)
  integer :: hcol(hnnzlim+n),hrow(hnnzlim+n),jcfun(jcnnzlim),jclen(m),jcvar(jcnnzlim),jcsta(m), &
       ucol(hnnzlim+jcnnzlim+n+m+n),udiag(n+m),urow(hnnzlim+jcnnzlim+n+m+n)
  real(kind=8) :: adddiag(n+m),b(n+m),c(m),ct(m),cu(m),cut(m),g(n),gpenbds(n),gtpenbds(n), &
       hpenbds(n),htpenbds(n),hval(hnnzlim+n),jcval(jcnnzlim),nl(n),tmp(n), &
       uval(hnnzlim+jcnnzlim+n+m+n),xt(n)
  character(len=50) :: merfun

  ! ================================================================================================
  ! PRESENTATION
  ! ================================================================================================

  if ( iprintinn .ge. 1 ) then
     write(* ,8000)
     write(10,8000)
  end if

  ! ================================================================================================
  ! CHECK PROBLEM
  ! ================================================================================================

  if ( count( equatn(1:m) ) .ne. m ) then
     write(*,*) 'There are inequalities.'
     stop
  end if

  ! ================================================================================================
  ! INITIALIZE
  ! ================================================================================================

  iter =  0

  pnopro = 0

  tolstp        = .false.
  bcktrfail     = .false.
  avoidstepctrl = .false.
  avoidbcktr    = .false.

  call lssset(lsssubACC,sclsubACC,.true.,.false.)

  ! ================================================================================================
  ! COMPUTE OBJECTIVE FUNCTION AND CONSTRAINTS
  ! ================================================================================================
  
  call ssetp(n,x,inform)
  if ( inform .ne. 0 ) return

  call sevalobjc(n,x,f,fu,m,c,cu,inform,ignscl=.false.)
  if ( inform .ne. 0 ) return

  if ( m .eq. 0 ) then
     cnorm   = 0.0d0
     cunorm  = 0.0d0
     cnorm2  = 0.0d0
     cunorm2 = 0.0d0
  else
     cnorm   = maxval( abs( c(1:m) ) )
     cunorm  = maxval( abs( cu(1:m) ) )
     cnorm2  = 0.5d0 * sum( c(1:m)  ** 2 )
     cunorm2 = 0.5d0 * sum( cu(1:m) ** 2 )
  end if

  ! Compute bounds violation only to save "best" visited point
  
  bdsnorm = 0.0d0
  bdsnorm = max( bdsnorm, maxval( l(lind(1:nlind)) - x(lind(1:nlind)) ) )
  bdsnorm = max( bdsnorm, maxval( x(uind(1:nuind)) - u(uind(1:nuind)) ) )

  if ( nlind .gt. 0 .or. nuind .gt. 0 ) then
     call cfghpenbds(n,x,l,u,rho,nlind,lind,lmu,nuind,uind,umu,fpenbds,gpenbds,hpenbds)

     f  = f  + fpenbds
     fu = fu + fpenbds
  end if

  ! ================================================================================================
  ! MAIN LOOP
  ! ================================================================================================

100 continue

  ! ================================================================================================
  ! COMPUTE FIRST DERIVATIVES
  ! ================================================================================================

  ! Gradient of the objective function and Jacobian of the constraints

  if ( gjaccoded ) then

     call sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform,ignscl=.false.)
     if ( inform .ne. 0 ) return

  else ! if ( gcoded .and. ( jaccoded .or. m .eq. 0 ) ) then

     call sevalg(n,x,g,inform,ignscl=.false.)
     if ( inform .ne. 0 ) return

     if ( m .eq. 0 ) then
        jcnnz = 0
        
     else
        ind(1:m) = .true.
        call sevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,jcnnzlim,inform,ignscl=.false.)
        if ( inform .ne. 0 ) return

        jcnnz = sum( jclen(1:m) )
     
        do j = 1,m
           jcfun(jcsta(j):jcsta(j)+jclen(j)-1) = j
        end do
     end if
     
  end if

  if ( nlind .gt. 0 .or. nuind .gt. 0 ) then
     g(1:n) = g(1:n) + gpenbds(1:n)
  end if
  
  ! Compute gradient of the Lagrangian of f(x) + h(x)^T lambda

  nl(1:n) = g(1:n)

  do i = 1,jcnnz
     nl(jcvar(i)) = nl(jcvar(i)) + lambda(jcfun(i)) * jcval(i)
  end do

  nlnorm = maxval( abs( nl(1:n) ) )
  
  ! ================================================================================================
  ! SAVE BEST POINT INFORMATION
  ! ================================================================================================

  if ( iter .eq. 0 ) then
     best_iter    = iter
     best_f       = f
     best_fu      = fu
     best_cnorm   = cnorm
     best_cunorm  = cunorm
     best_bdsnorm = bdsnorm
     best_nlnorm  = nlnorm
     
     best_x(1:n)      = x(1:n)
     best_lambda(1:m) = lambda(1:m)

  else

     if ( ( max( best_bdsnorm, best_cunorm ) .gt. epsfeas .and. &
            ( max( bdsnorm, cunorm ) .lt. max( best_bdsnorm, best_cunorm ) .or. &
            ( max( bdsnorm, cunorm ) .le. max( best_bdsnorm, best_cunorm ) .and. &
              fu .lt. best_fu ) .or. &
            ( max( bdsnorm, cunorm ) .le. max( best_bdsnorm, best_cunorm ) .and. &
              fu .le. best_fu .and. nlnorm .lt. best_nlnorm ) ) ) .or. &
          ( max( best_bdsnorm, best_cunorm ) .le. epsfeas .and. max( bdsnorm, cunorm ) .le. epsfeas .and. &
            ( fu .lt. best_fu .or. ( fu .le. best_fu .and. nlnorm .lt. best_nlnorm ) ) ) ) then
        best_iter    = iter
        best_f       = f
        best_fu      = fu
        best_cnorm   = cnorm
        best_cunorm  = cunorm
        best_bdsnorm = bdsnorm
        best_nlnorm  = nlnorm

        best_x(1:n)      = x(1:n)
        best_lambda(1:m) = lambda(1:m)
     end if
  end if

  ! ================================================================================================
  ! WRITE INFORMATION OF THE CURRENT POINT
  ! ================================================================================================

  if ( iprintinn .ge. 1 ) then
     star = ' '
     if ( iter .eq. best_iter ) star = '*'
     write(* ,8010) iter,star,fu,f,cunorm,cnorm,cunorm2,cnorm2,nlnorm
     write(10,8010) iter,star,fu,f,cunorm,cnorm,cunorm2,cnorm2,nlnorm
     write(* ,8015) nprint,maxval(abs(x(1:n))),(x(i),i=1,nprint)
     write(10,8015) nprint,maxval(abs(x(1:n))),(x(i),i=1,nprint)
  end if

  if ( iprintctl(5) ) then
     open(20,file='accpro-tabline.out')
     write(20,9100) fu,cunorm,bdsnorm,nlnorm,f,cnorm,0,iter,best_fu,best_cunorm,best_bdsnorm, &
          best_nlnorm,best_f,best_cnorm,best_iter,9,inform,999.99d0
     close(20)
  end if

  ! ================================================================================================
  ! TEST STOPPING CRITERIA
  ! ================================================================================================

  ! Note that if tolerances for stopping are satisfied at least once,
  ! tolstp remains true forever.

  if ( cunorm .le. epsfek .and. nlnorm .le. epsopk ) then
     tolstp = .true.

     if ( iprintinn .ge. 1 ) then
        write(* ,8990)
        write(10,8990)
     end if
  end if

  if ( ( tolstp .and. ( .not. scale .or. pnopro .ge. 3 .or. iter - best_iter .ge. 3 ) ) .or. &
       ( cunorm .le. epsfek .and. nlnorm .le. macheps12 * epsopk ) ) then
     accinfo = 0

     if ( iprintinn .ge. 1 ) then
        write(* ,9000)
        write(10,9000)
     end if

     go to 500
  end if

  if ( pnopro .ge. 10 ) then
     accinfo = 2
     
     if ( iprintinn .ge. 1 ) then
        write(* ,9020)
        write(10,9020)
     end if
     
     go to 500
  end if

  if ( nlind .gt. 0 .or. nuind .gt. 0 ) then

     if ( cunorm .le. epsfek .and. fu .le. - 1.0d+10 ) then
        accinfo = 1
     
        if ( iprintinn .ge. 1 ) then
           write(* ,9010)
           write(10,9010)
        end if
        
        go to 500
     end if

     if ( iter .ge. 1000 ) then
        accinfo = 3
        
        if ( iprintinn .ge. 1 ) then
           write(* ,9030)
           write(10,9030)
        end if
        
        go to 500
     end if

  end if
  
  ! ================================================================================================
  ! DO AN ITERATION
  ! ================================================================================================

  iter = iter + 1

  ! ================================================================================================
  ! COMPUTE HESSIAN OF THE LAGRANGIAN
  ! ================================================================================================

  call sevalhl(n,x,m,lambda,hrow,hcol,hval,hnnz,hnnzlim,inform,ignscl=.false.)
  if ( inform .ne. 0 ) return

  if ( nlind .gt. 0 .or. nuind .gt. 0 ) then
     if ( hnnz + n .gt. hnnzlim ) then
        write(*,*) 'There was an error in accpro. hnnzlim should be increased in at least n!'
        stop
     end if
     
     hrow(hnnz+1:hnnz+n) = (/(i,i=1,n)/)
     hcol(hnnz+1:hnnz+n) = (/(i,i=1,n)/)
     hval(hnnz+1:hnnz+n) = hpenbds(1:n)
     hnnz = hnnz + n
  end if
  
  ! ================================================================================================
  ! ASSEMBLE THE JACOBIAN OF THE KKT SYSTEM
  ! ================================================================================================

  unnz = 0

  udiag(1:n+m) = 0

  do i = 1,hnnz
     if ( hval(i) .ne. 0.0d0 ) then
        unnz = unnz + 1
        urow(unnz) = hrow(i)
        ucol(unnz) = hcol(i)
        uval(unnz) = hval(i)
        if ( urow(unnz) .eq. ucol(unnz) ) then
           ! If there are duplicated diagonal entrances, a pointer to
           ! the last one is ok. The relevant thing is to add (below)
           ! entrances for the missing diagonal elements.
           udiag(urow(unnz)) = unnz
        end if
     end if
  end do

  do i = 1,jcnnz
     if ( jcval(i) .ne. 0.0d0 ) then
        unnz = unnz + 1
        urow(unnz) = jcfun(i) + n
        ucol(unnz) = jcvar(i)
        uval(unnz) = jcval(i)
     end if
  end do

  do i = 1,n + m
     if ( udiag(i) .eq. 0 ) then
        unnz = unnz + 1
        urow(unnz) = i
        ucol(unnz) = i
        uval(unnz) = 0.0d0

        udiag(i) = unnz
     end if
  end do

  ! ================================================================================================
  ! ANALYSE SPARSITY PATTERN
  ! ================================================================================================

  nsys = n + m

  call lssana(nsys,unnz,urow,ucol,uval,udiag,lssinfo)

  if ( lssinfo .eq. 6 ) then
     ! MEMORY FAILURE
     accinfo = 6

     if ( iprintinn .ge. 1 ) then
        write(* ,9060)
        write(10,9060)
     end if

     call lssend()
     go to 500
  end if

  ! ================================================================================================
  ! SOLVE THE NEWTONIAN SYSTEM
  ! ================================================================================================

  ! ================================================================================================
  ! SET FIRST REGULARIZATION
  ! ================================================================================================

  if ( iprintinn .ge. 1 ) then
     write(* ,8020)
     write(10,8020)
  end if

  if ( iter .eq. 1 ) then
     epsadd1 = 0.0d0
     epsadd2 = 0.0d0
     if ( m .gt. n ) epsadd2 = macheps12

  else
     if ( xcutted ) then
        epsadd1 = max( macheps12, 3.0d0 * epsadd1 )
     else
        epsadd1 = 1.0d-01 * epsadd1
     end if

     epsadd2 = 1.0d-01 * epsadd2
  end if

  ! ================================================================================================
  ! FACTORIZE THE JACOBIAN OF THE NEWTONIAN SYSTEM
  ! ================================================================================================

200 continue

  adddiag(1:n)     =   epsadd1
  adddiag(n+1:n+m) = - epsadd2

  if ( iprintinn .ge. 2 ) then
     write(* ,8030) epsadd1,epsadd2
     write(10,8030) epsadd1,epsadd2
  end if

  call lssfac(nsys,unnz,urow,ucol,uval,udiag,adddiag,pind,pval, &
       nneigv,nrank,lssinfo)

  if ( lssinfo .eq. 6 ) then
     ! MEMORY FAILURE
     accinfo = 6

     if ( iprintinn .ge. 1 ) then
        write(* ,9060)
        write(10,9060)
     end if

     call lssend()
     go to 500
  end if

  if ( nneigv .lt. m .or. nrank - nneigv .lt. n ) then
     if ( nneigv .lt. m ) then
        epsadd2 = max( macheps12, 3.0d0 * epsadd2 )
     end if
     if ( nrank - nneigv .lt. n ) then
        epsadd1 = max( macheps12, 3.0d0 * epsadd1 )
     end if

     if ( iprintinn .ge. 2 ) then
        write(* ,8040) n,m,nrank-nneigv,nneigv,n+m-nrank
        write(10,8040) n,m,nrank-nneigv,nneigv,n+m-nrank
     end if
     
     if ( epsadd1 .lt. bignum .and. epsadd2 .lt. bignum ) then
        go to 200
     end if

     accinfo = 5

     if ( iprintinn .ge. 2 ) then
        write(* ,9050)
        write(10,9050)
     end if

     call lssend()
     go to 500

  else
     if ( iprintinn .ge. 2 ) then
        write(* ,8050)
        write(10,8050)
     end if
  end if

  ! ================================================================================================
  ! SOLVE TRIANGULAR SYSTEMS
  ! ================================================================================================

  ! Set right-hand side

  b(1:n)     = - nl(1:n)
  b(n+1:n+m) = - c(1:m)

  ! Solve triangular systems and deallocate memory

  call lsssol(nsys,b)

  call lssend()

  ! ================================================================================================
  ! STEP CONTROL
  ! ================================================================================================

  ! It means that the backtracking failed in the previous iteration
  if ( .not. avoidstepctrl .and. ( bcktrfail .or. nobacktrack ) ) then
     avoidstepctrl = .true.
     if ( iprintinn .ge. 3 ) then
        write(* ,*) 'Step control disabled.'
        write(10,*) 'Step control disabled.'
     end if
  end if

  ! Control the size of the unitary step

  xcutted = .false.

  if ( .not. avoidstepctrl ) then

     xnorm  = maxval( abs( x(1:n) ) )
     dxnorm = maxval( abs( b(1:n) ) )

     lnorm  = max( 0.0d0, maxval( abs( lambda(1:m) ) ) )
     dlnorm = max( 0.0d0, maxval( abs( b(n+1:n+m) ) ) )

     if ( iprintinn .ge. 2 ) then
        write(* ,8060) xnorm,dxnorm,lnorm,dlnorm
        write(10,8060) xnorm,dxnorm,lnorm,dlnorm
     end if

     if ( dxnorm .gt. 1.0d+02 * max( 1.0d0, xnorm ) ) then
        xcutted = .true.
        b(1:n) = b(1:n) * 1.0d+2 * max( 1.0d0, xnorm ) / dxnorm

        if ( iprintinn .ge. 2 ) then
           write(* ,8061) 1.0d+02 * max( 1.0d0, xnorm )
           write(10,8061) 1.0d+02 * max( 1.0d0, xnorm )
        end if
     else
        if ( iprintinn .ge. 2 ) then
           write(* ,8062)
           write(10,8062)
        end if
     end if

     if ( dlnorm .gt. 1.0d+02 * max( 1.0d0, lnorm ) ) then
        b(n+1:n+m) = b(n+1:n+m) * 1.0d+02 * max( 1.0d0, lnorm ) / dlnorm

        if ( iprintinn .ge. 2 ) then
           write(* ,8063) 1.0d+02 * max( 1.0d0, lnorm )
           write(10,8063) 1.0d+02 * max( 1.0d0, lnorm )
        end if
     else
        if ( iprintinn .ge. 2 ) then
           write(* ,8064)
           write(10,8064)
        end if
     end if
  end if

 ! ================================================================================================
  ! BACKTRACKING
  ! ================================================================================================

  ! It means that the backtracking failed in the previous iteration
  if ( .not. avoidbcktr .and. ( bcktrfail .or. nobacktrack ) ) then
     avoidbcktr = .true.
     if ( iprintinn .ge. 3 ) then
        write(* ,*) 'Backtracking disabled.'
        write(10,*) 'Backtracking disabled.'
     end if
  end if

  ! Lagrangian and Augmented Lagrangian times (1/rho) (with rho = 1 /
  ! epsadd2), evaluated at (x,lambda)
  lagra  = f + dot_product( lambda(1:m), c(1:m)  )
  auglag = epsadd2 * lagra  + cnorm2

  if ( cunorm .gt. epsfeas ) then
     merfun = '(1/rho) times the Augmented Lagrangian function'
  else ! if ( cunorm .le. epsfeas ) then
     if ( epsadd2 .eq. 0.0d0 ) then
        merfun = 'Lagrangian function'
     else
        merfun = '(1/rho) times the Augmented Lagrangian function'
     end if
  end if

  if ( iprintinn .ge. 2 ) then
     write(* ,8070) lagra,auglag,merfun
     write(10,8070) lagra,auglag,merfun
  end if

  if ( iprintinn .ge. 3 ) then
     write(* ,8071)
     write(10,8071)
  end if

  step = 1.0d0
  xt(1:n) = x(1:n) + b(1:n)

300 continue

  call ssetp(n,xt,inform)
  if ( inform .ne. 0 ) return

  call sevalobjc(n,xt,ft,fut,m,ct,cut,inform,ignscl=.false.)
  if ( inform .ne. 0 ) return

  if ( m .eq. 0 ) then
     cnormt   = 0.0d0
     cunormt  = 0.0d0
     cnorm2t  = 0.0d0
     cunorm2t = 0.0d0
  else
     cnormt   = maxval( abs( ct(1:m) ) )
     cunormt  = maxval( abs( cut(1:m) ) )
     cnorm2t  = 0.5d0 * sum( ct(1:m)  ** 2 )
     cunorm2t = 0.5d0 * sum( cut(1:m) ** 2 )
  end if

  bdsnormt = 0.0d0
  bdsnormt = max( bdsnormt, maxval( l(lind(1:nlind)) - xt(lind(1:nlind)) ) )
  bdsnormt = max( bdsnormt, maxval( xt(uind(1:nuind)) - u(uind(1:nuind)) ) )

  if ( nlind .gt. 0 .or. nuind .gt. 0 ) then
     call cfghpenbds(n,xt,l,u,rho,nlind,lind,lmu,nuind,uind,umu,ftpenbds,gtpenbds,htpenbds)

     ft  = ft  + ftpenbds
     fut = fut + ftpenbds
  end if
  
  ! Lagrangian and Augmented Lagrangian times (1/rho) (with rho = 1 /
  ! epsadd2), evaluated at (xt,lambda)
  lagrat = ft + dot_product( lambda(1:m), ct(1:m) )
  auglagt = epsadd2 * lagrat + cnorm2t

  if ( iprintinn .ge. 3 ) then
     write(* ,8072) step,cnorm2t,lagrat,auglagt
     write(10,8072) step,cnorm2t,lagrat,auglagt
  end if

  ! Backtracking
  if ( .not. avoidbcktr .and. &
       !.not. ( cunorm .gt. epsfek .and. auglagt .le. auglag ) .and. &
       !.not. ( cunorm .le. epsfek .and. epsadd2 .eq. 0.0d0 .and. lagrat  .le. lagra  ) .and. &
       !.not. ( cunorm .le. epsfek .and. epsadd2 .ne. 0.0d0 .and. auglagt .le. auglag ) ) then
       .not. (         cunorm .le. epsfek .and. epsadd2 .eq. 0.0d0   .and. lagrat  .le. lagra  ) .and. &
       .not. ( .not. ( cunorm .le. epsfek .and. epsadd2 .eq. 0.0d0 ) .and. auglagt .le. auglag ) ) then
     step = 0.5d0 * step
     xt(1:n) = x(1:n) + step * b(1:n)

     go to 300
  end if

  if ( .not. avoidbcktr .and. &
       !.not. ( cunorm .gt. epsfek .and. auglagt .lt. auglag ) .and. &
       !.not. ( cunorm .le. epsfek .and. epsadd2 .eq. 0.0d0 .and. lagrat  .lt. lagra  ) .and. &
       !.not. ( cunorm .le. epsfek .and. epsadd2 .ne. 0.0d0 .and. auglagt .lt. auglag ) ) then
       .not. (         cunorm .le. epsfek .and. epsadd2 .eq. 0.0d0   .and. lagrat  .lt. lagra  ) .and. &
       .not. ( .not. ( cunorm .le. epsfek .and. epsadd2 .eq. 0.0d0 ) .and. auglagt .lt. auglag ) ) then
     bcktrfail = .true.
  end if

  if ( count( abs( xt(1:n) - x(1:n) ) .le. macheps * max( 1.0d0, abs( x(1:n) ) ) ) .eq. n ) then
     pnopro = pnopro + 1
  else
     pnopro = 0
  end if

  x(1:n)      = xt(1:n)
  lambda(1:m) = lambda(1:m) + step * b(n+1:n+m)

  f       = ft
  fu      = fut
  c(1:m)  = ct(1:m)
  cu(1:m) = cut(1:m)
  cnorm   = cnormt
  cunorm  = cunormt
  cnorm2  = cnorm2t
  cunorm2 = cunorm2t
  bdsnorm = bdsnormt
  
  fpenbds = ftpenbds
  gpenbds(1:n) = gtpenbds(1:n)
  hpenbds(1:n) = htpenbds(1:n)
  
  go to 100

  ! ================================================================================================
  ! END OF MAIN LOOP
  ! ================================================================================================

500 continue

! ================================================================================================
! NON-EXECUTABLE STATEMENTS
! ================================================================================================

 8000 format(/,' Acceleration Process in action!')
 8010 format(/,' Acceleration Process Iteration',33X,' = ',5X,I6,1A, &
             /,' Objective function value              (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' Sup-norm of constraints               (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' 1/2 squared two norm of constraints   (unscaled = ',1PD11.4,')',1X,' = ',1PD11.4, &
             /,' Sup-norm of the gradient of the Lagrangian                    ',1X,' = ',1PD11.4)
 8015 format(/,' Current point first ',I1,' components (sup-norm: ',1PD11.4,'): ',/,6(1X,1P,D11.4))
 8020 format(/,' Solving Newtonian system.')
 8030 format(  ' Value added to the NW = ',1PD12.4,/,' Value subtracted from to the SE = ',1PD12.4)
 8040 format(  " Wrong Jacobian's inertia. ", &
             /,' Desired POS = ',I7,' NEG = ',I7, &
             /,' Actual  POS = ',I7,' NEG = ',I7,' NULL = ',I7)
 8050 format(  ' The desired inertia was obtained.') 
 8060 format(/,' Step control. ', &
             /,' Sup-norm of x      = ',1PD11.4,' Sup-norm of dx      = ',1PD11.4, &
             /,' Sup-norm of lambda = ',1PD11.4,' Sup-norm of dlambda = ',1PD11.4)
 8061 format(  ' The sup-norm of dx was reduced to ',1PD11.4)
 8062 format(  ' The sup-norm of dx was NOT modified.')
 8063 format(  ' The sup-norm of dlambda was reduced to ',1PD11.4)
 8064 format(  ' The sup-norm of dlambda was NOT modified.')
 8070 format(/,' Backtracking (information on the current point and the search direction):', &
             /,' Value of the Lagrangian at the current point (x,lambda)',8X,' = ',1PD11.4, &
             /,' Value of (1/rho) times the Augmented Lagrangian at (x,lambda)', &
             /,' (with rho = 1 over the value subtracted from the SE)',10X,' = ',1PD11.4, &
             /,' Merit function: ',A50)
 8071 format(/,' Backtracking procedure:', &
             /,' cnorm2t : 1/2 squared two norm of constraints at trial.', &
             /,' lagrat  : Lagrangian function at trial.', &
             /,' auglagt : (1/rho) times Augmented Lagrangian function at trial.')
 8072 format(  ' Step = ',1PD11.4,' cnorm2t = ',1PD11.4,' lagrat = ',1PD11.4,' auglagt = ',1PD11.4)
 8990 format(/,' Required tolerances in the KKT system were satisfied at the last iterate.', &
             /,' However, if progress in the objective function appears to be being done, ', &
             /,' due to scaling reasons, optimization will continue.')
 9000 format(/,' Flag of Acceleration Process = KKT system solved!')
 9010 format(/,' Flag of Acceleration Process = The objective functions appears to be going to minus infinity.')
 9020 format(/,' Flag of Acceleration Process = Current point did not change ', &
           'substantially or best point was not improved during ten consecutive iterations.')
 9030 format(/,' Flag of Acceleration Process = Maximum of iterations.')
 9050 format(/," Flag of Acceleration Process = Jacobian's wrong inertia.")
 9060 format(/,' Flag of Acceleration Process = Memory failure.')

 9100 format(1X,1P,D24.16,3(1X,1P,D7.1),1X,1P,D24.16,1X,1P,D7.1,2(1X,I12), &
             1X,1P,D24.16,3(1X,1P,D7.1),1X,1P,D24.16,1X,1P,D7.1,1X,I12,1X,I1,1X,I3,0P,F8.2)

end subroutine accpro

! **************************************************************************************************
! **************************************************************************************************

subroutine cfghpenbds(n,x,l,u,rho,nlind,lind,lmu,nuind,uind,umu,fpenbds,gpenbds,hpenbds)

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n,nlind,nuind
  real(kind=8), intent(in) :: rho
  real(kind=8), intent(out) :: fpenbds

  ! ARRAY ARGUMENTS
  integer, intent(in) :: lind(nlind),uind(nuind)
  real(kind=8), intent(in) :: l(n),lmu(nlind),u(n),umu(nuind),x(n)
  real(kind=8), intent(out) :: gpenbds(n),hpenbds(n)

  ! Compute fpenbds as the penalization (augmented Lagrangian) of
  ! bounds, as well as first and second-order derivatives

  ! LOCAL SCALARS
  integer :: i
  real(kind=8) :: dpdy,dpdy2,p,y

  fpenbds = 0.0d0
  gpenbds(1:n) = 0.0d0
  hpenbds(1:n) = 0.0d0

  do i = 1,nlind
     y = l(lind(i)) - x(lind(i))
     
     if ( lmu(i) + rho * y .ge. 0.0d0 ) then
        p     = y * ( lmu(i) + 0.5d0 * rho * y )
        dpdy  = lmu(i) + rho * y
        dpdy2 = rho
     else
        p     = - 0.5d0 * lmu(i) ** 2 / rho
        dpdy  = 0.0d0
        dpdy2 = 0.0d0
     end if

     fpenbds = fpenbds + p
     gpenbds(lind(i)) = gpenbds(lind(i)) - dpdy
     hpenbds(lind(i)) = hpenbds(lind(i)) + dpdy2
  end do
  
  do i = 1,nuind
     y = x(uind(i)) - u(uind(i))
     
     if ( umu(i) + rho * y .ge. 0.0d0 ) then
        p     = y * ( umu(i) + 0.5d0 * rho * y )
        dpdy  = umu(i) + rho * y
        dpdy2 = rho
     else
        p     = - 0.5d0 * umu(i) ** 2 / rho
        dpdy  = 0.0d0
        dpdy2 = 0.0d0
     end if
     
     fpenbds = fpenbds + p
     gpenbds(uind(i)) = gpenbds(uind(i)) + dpdy
     hpenbds(uind(i)) = hpenbds(uind(i)) + dpdy2
  end do
  
end subroutine cfghpenbds

! **************************************************************************************************
! **************************************************************************************************


subroutine newtonkkt2(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,f,cnorm,nlnorm,iter,msqiter,accinfo,inform)

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: m,n
  integer,      intent(inout) :: accinfo,inform,iter,msqiter
  real(kind=8), intent(in)    :: epsfeas,epsopt
  real(kind=8), intent(inout) :: cnorm,f,nlnorm

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(inout) :: l(n),lambda(m),u(n),x(n)

  ! LOCAL SCALARS
  integer :: outiter,totiter,best_outit,secoinfo
  real(kind=8) :: bdsnorm,cunorm,efstain,eostain,fu,best_f,best_cnorm,best_bdsnorm, &
       best_nlnorm,best_fu,best_cunorm

  ! LOCAL ARRAYS
  real(kind=8) :: best_x(n),best_lambda(m)
  
  call seco(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,efstain,eostain,.false.,.false., &
       f,cnorm,bdsnorm,nlnorm,fu,cunorm,outiter,totiter,best_f,best_cnorm,best_bdsnorm, &
       best_nlnorm,best_fu,best_cunorm,best_outit,best_x,best_lambda,secoinfo,inform)

  x(1:n)      = best_x(1:n)
  lambda(1:m) = best_lambda(1:m)

  f      = best_fu
  cnorm  = max( best_cunorm, best_bdsnorm )
  nlnorm = best_nlnorm
  
  iter    = totiter
  msqiter = 0
  accinfo = secoinfo

end subroutine newtonkkt2

