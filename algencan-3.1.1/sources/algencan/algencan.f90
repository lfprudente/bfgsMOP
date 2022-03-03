! ******************************************************************
! ******************************************************************
subroutine algencan(fsub,gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub, &
     gjacpsub,hlsub,hlpsub,jcnnzmax,hnnzmax,epsfeas,epsopt,efstain, &
     eostain,efacc,eoacc,outputfnm,specfnm,nvparam,vparam,norig, &
     xorig,lorig,uorig,m,lambda,equatn,linear,coded,checkder,fu, &
     cunorm,snorm,nlnorm,inform,solinfo)

  use modmachconst, only: inimachconst,bignum,macheps,macheps12,macheps13,macheps23
  use modouttyp,    only: setouttyp,iniouttyp,iprint,iprintctl,iprintinn,iprintout,mprint,ncomp,nprint
  use modalgparam,  only: setalgparam,hptype,precond,innslvr,innercall,ignoref,lsssubACC, &
       lsssubNW,lsssubTR,sclsubACC,sclsubNW,sclsubTR,skipacc,useustp,maxoutitgiven,&
       maxinnitgiven,maxaccitgiven,rhoinigiven,rhomaxgiven
  use modalgconst,  only: n0,n1,n2

  use lssma57
  use lssma86
  use lssma97
  use problvls, only: setproblvls,sinip,sendp,sf,sevalobjc,sevalgjac,ssetp
  use problvlt, only: setproblvlt,tinip,tendp
  use problvlu, only: setproblvlu,uinip,uendp
  use problvlv, only: setproblvlv,vinip,vendp,reperr,fcnt
  use probgiven, only: iniprobgiven,firstde,seconde,truehpr,gjacpcoded,jcnnzlim,hnnzlim

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(in)    :: checkder
  integer,      intent(in)    :: hnnzmax,jcnnzmax,m,norig,nvparam
  integer,      intent(out)   :: inform,solinfo
  real(kind=8), intent(inout) :: cunorm,efacc,efstain,eoacc,eostain, &
                                 epsfeas,epsopt,fu,nlnorm,snorm

  ! ARRAY ARGUMENTS
  character(len=80), intent(in)    :: specfnm,outputfnm,vparam(nvparam)
  logical,           intent(in)    :: coded(11),linear(m)
  logical,           intent(inout) :: equatn(m)
  real(kind=8),      intent(in)    :: lorig(norig),uorig(norig)
  real(kind=8),      intent(inout) :: lambda(m),xorig(norig)

  ! EXTERNAL SUBROUTINES ARGUMENTS
  external :: fsub,gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub,gjacpsub, &
              hlsub,hlpsub

  ! LOCAL SCALARS
  logical      :: innfail,lss,scl
  integer      :: alinfo,allocerr,geninfo,iter,maxit,n,nwcalls,nwtotit, &
                  outiter,totiter,msqcalls,msqtotit,i,best_outit, &
                  secoinfo
  real(kind=8) :: bdsnorm,cnorm,f,ncnorm,rsupn,best_f,best_cnorm, &
                  best_bdsnorm,best_nlnorm,best_fu,best_cunorm
  real(kind=4) :: time

  ! LOCAL ARRAYS
  character(len=15)     :: solver,subname
  real(kind=8)          :: rho(m)
  real(kind=8), pointer :: best_lambda(:),best_x(:),l(:),u(:),x(:)
  real(kind=4)          :: tdum(2)

  ! DATA STATEMENTS
  data tdum/0.0,0.0/

  ! EXTERNAL FUNCTIONS
  real(kind=4) :: dtime

  ! ==================================================================
  ! Start timing
  ! ==================================================================

  ! CPU time for sequential applications

  time = dtime(tdum)
  
  ! ==================================================================
  ! Error tracker
  ! ==================================================================

  inform = 0

  ! ==================================================================
  ! Initialization
  ! ==================================================================

  ! Set machine-dependent constants (module machconst)
  call inimachconst()

  ! Output control (module outtyp)
  call iniouttyp()

  ! ==================================================================
  ! Set default values in 'prob' modules (probgiven, problvl[s|t|u|v])
  ! ==================================================================

  ! Initiate module probgiven
  call iniprobgiven(fsub,gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub, &
       gjacpsub,hlsub,hlpsub,coded,m,jcnnzmax,hnnzmax,inform)
  if ( inform .ne. 0 ) return

  ! Set default values for variables in module problvlv
  !call setproblvlv(val_safemode = .false.)
  call setproblvlv(val_safemode = .true., val_nsfmax = 10 )

  ! Set default values for variables in module problvlu
  call setproblvlu(val_rmfixv = .true.)

  ! Set default values for variables in module problvlt
  call setproblvlt(val_slacks = .false.)

  ! Set default values for variables in module problvls
  call setproblvls(val_scale = .true.)

  ! ==================================================================
  ! Set default values for algoritmic parameters (module algparam)
  ! ==================================================================

  ! Hessian-vector product strategy: HAPPRO, INCQUO or TRUEHP. TRUEHP
  ! is the default option. If the proper subroutines were not coded by
  ! the user, then HAPPRO is used instead. HAPPRO reduces to INCQUO in
  ! the unconstrained and bound-constrained cases and switches between
  ! INCQUO and the product by an approximation of the Hessian in the
  ! constrained case.

  if ( truehpr ) then
     call setalgparam(val_hptype = 'TRUEHP')
  else if ( .not. gjacpcoded ) then
     call setalgparam(val_hptype = 'HAPPRO')
  else
     call setalgparam(val_hptype = 'INCQUO')
  end if

  ! Conjugate Gradients preconditioner

  if ( .not. gjacpcoded ) then
     call setalgparam(val_precond = 'QNGN')
  else
     call setalgparam(val_precond = 'NONE')
  end if

  ! Subroutines for Linear systems

  if ( lss_ma57() ) then
     call setalgparam(val_lsssubTR = 'MA57')
     call setalgparam(val_sclsubTR = 'NONE')
  else
     call setalgparam(val_lsssubTR = 'NONE')
     call setalgparam(val_sclsubTR = 'NONE')
  end if

  if ( lss_ma57() ) then
     call setalgparam(val_lsssubNW = 'MA57')
     call setalgparam(val_sclsubNW = 'NONE')
  else if ( lss_ma97() ) then
     call setalgparam(val_lsssubNW = 'MA97')
     call setalgparam(val_sclsubNW = 'NONE')
  else if ( lss_ma86() ) then
     call setalgparam(val_lsssubNW = 'MA86')
     call setalgparam(val_sclsubNW = 'NONE')
  else
     call setalgparam(val_lsssubNW = 'NONE')
     call setalgparam(val_sclsubNW = 'NONE')
  end if

  call setalgparam(val_lsssubACC = lsssubNW)
  call setalgparam(val_sclsubACC = sclsubNW)

  ! Inner solver method: TR (trust regions, i.e. BETRA), NW (Newtonian
  ! system with the use of a direct solver, i.e. MA57AD), TN + TRUEHP
  ! (Truncated Newton and true Hessian-vector product) or TN + HAPPRO
  ! (Truncated Newton and product with Hessian approximation). If the
  ! proper subroutines were not provided by the user then TN + HAPPRO
  ! is used.

  if ( seconde .and. lsssubTR .ne. 'NONE' .and. norig .le. n0 ) then
     call setalgparam(val_innslvr = 'TR')
  else if ( seconde .and. lsssubNW .ne. 'NONE' .and. norig .le. n1 ) then
     call setalgparam(val_innslvr = 'NW')
  else
     call setalgparam(val_innslvr = 'TN')
!!$     if ( norig .gt. n2 ) then
!!$        call setalgparam(val_hptype  = 'HAPPRO')
!!$     end if
  end if

  ! Ignore objective function (to only find a feasible point by
  ! minimizing 1/2 of the squared infeasibility)
  call setalgparam(val_ignoref = .false.)

  ! Acceleration step

  if ( seconde .and. lsssubACC .ne. 'NONE' ) then
     call setalgparam(val_skipacc = .false.)
  else
     call setalgparam(val_skipacc = .true.)
  end if

  ! Maximum number of outer iterations, inner iterations, and initial
  ! penalty parameter

  call setalgparam(val_maxaccitgiven = - 1)
  call setalgparam(val_maxinnitgiven = - 1)
  call setalgparam(val_maxoutitgiven = - 1)
  call setalgparam(val_rhoinigiven   = - 1.0d0)
  call setalgparam(val_rhomaxgiven   = - 1.0d0)

  ! Flags used to make inner calls to Gencan

  call setalgparam(val_innercall = .false.)
  call setalgparam(val_useustp   = .false.)

  ! ==================================================================
  ! Open output file
  ! ==================================================================

  if ( outputfnm .ne. '' ) then
     open(unit=10,file=outputfnm,status='replace')
  else
     open(unit=10,               status='scratch')
  end if

  ! ==================================================================
  ! Modifiy default values using specification file
  ! ==================================================================

  call procpar(epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
       specfnm,nvparam,vparam)

  ! ==================================================================
  ! Initialize problem data structures
  ! ==================================================================

  allocate(x(norig),l(norig),u(norig),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     inform = - 93
     subname = 'ALGENCAN'
     call reperr(inform,subname)
     return
  end if

  n = norig

  x(1:norig) = xorig(1:norig)
  l(1:norig) = lorig(1:norig)
  u(1:norig) = uorig(1:norig)

  ! Initiate module problvlv
  call vinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
  if ( inform .ne. 0 ) return

  ! Initiate module problvlu
  call uinip(n,x,l,u,m,lambda,equatn,linear,coded,inform)
  if ( inform .ne. 0 ) return

  ! Initiate module problvlt
  call tinip(n,x,l,u,m,lambda,equatn,linear,coded,inform)
  if ( inform .ne. 0 ) return

  ! Initiate module problvls
  call sinip(n,x,l,u,m,lambda,equatn,linear,coded,inform)
  if ( inform .ne. 0 ) return

  call setouttyp(val_nprint = min( n, ncomp ), val_mprint = min( m, ncomp ))

  ! ==================================================================
  ! Allocate arrays modules svdgrad, svdhess, and minsq
  ! ==================================================================

  allocate(best_x(n),best_lambda(m),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     inform = - 93
     subname = 'ALGENCAN'
     call reperr(inform,subname)
     return
  end if

  call allocatearrays(n,m,jcnnzlim,hnnzlim,allocerr)

  if ( allocerr .ne. 0 ) then
     inform = - 93
     subname = 'ALGENCAN'
     call reperr(inform,subname)
     return
  end if

  ! ==================================================================
  ! Call the solver
  ! ==================================================================

  if ( .false. ) then

     solver = 'SECO'
     
     call seco(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,efstain,eostain, &
          .false.,.false.,f,cnorm,bdsnorm,nlnorm,fu,cunorm,outiter,totiter,   &
          best_f,best_cnorm,best_bdsnorm,best_nlnorm,best_fu,best_cunorm,     &
          best_outit,best_x,best_lambda,secoinfo,inform)

  else if ( .not. ignoref .and. m .gt. 0 ) then

     ! ALGENCAN for PNL problems

     solver = 'AUGLAG'
     
     call auglag(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,     &
          efstain,eostain,efacc,eoacc,f,cnorm,snorm,nlnorm,fu,      &
          cunorm,best_fu,best_cunorm,best_f,best_cnorm,best_nlnorm, &
          ncnorm,rsupn,outiter,totiter,nwcalls,nwtotit,msqcalls,    &
          msqtotit,innfail,alinfo,inform)

     solinfo = alinfo
     if ( inform .ne. 0 ) solinfo = 9

  else
     ! GENCAN for box-constrained problems and feasibility problems

     solver = 'GENCAN'
     
     maxit = 999999999

     if ( maxinnitgiven .gt. 0 ) then
        maxit = maxinnitgiven
     end if

     ! Used in feasibility problems (ignoref=true). With lambda=0 and
     ! rho=1, to minimize 1/2 of the squared infeasibility coincides
     ! with minimizing the augmented Lagrangian.

     lambda(1:m) = 0.0d0
     rho(1:m) = 1.0d0

     call gencan(n,x,l,u,m,lambda,equatn,linear,rho,epsfeas,epsopt, &
          efstain,eostain,maxit,iter,f,nlnorm,cnorm,cunorm,geninfo, &
          inform)

     solinfo  = geninfo
     if ( inform .ne. 0 ) solinfo = 9

     ncnorm   = 0.0d0
     rsupn    = 0.0d0

     outiter  = 0
     totiter  = iter
     nwcalls  = 0
     nwtotit  = 0
     msqcalls = 0
     msqtotit = 0
     innfail  = .false.

     if ( ignoref ) then
        f  = 0.0d0
        fu = 0.0d0
     else
        fu = f
     end if

     best_f       = f
     best_fu      = fu
     best_cnorm   = cnorm
     best_cunorm  = cunorm
     best_nlnorm  = nlnorm
  end if

  if ( inform .lt. 0 ) return

  ! ==================================================================
  ! End problem data structures
  ! ==================================================================

  ! Unscale Lagrange multipliers
  call sendp(n,x,l,u,m,lambda,equatn,linear,inform)
  if ( inform .lt. 0 ) return

  ! Remove slack variables  (x, l and u may be resized)
  call tendp(n,x,l,u,m,lambda,equatn,linear,inform)
  if ( inform .lt. 0 ) return

  ! Restore fixed variables (x, l and u may be resized)
  call uendp(n,x,l,u,m,lambda,equatn,linear,inform)
  if ( inform .lt. 0 ) return

  ! Save solution (primal and dual variables)
  call vendp(n,x,l,u,m,lambda,equatn,linear,inform)
  if ( inform .lt. 0 ) return

  xorig(1:n) = x(1:n)

  ! ==================================================================
  ! Stop timing
  ! ==================================================================

  ! CPU time for sequential applications

  time = dtime(tdum)
  time = tdum(1)

  if ( iprintctl(4) ) then
     write(* ,9000) time
     write(10,9000) time
  end if
  
  ! ==================================================================
  ! Close output file
  ! ==================================================================

  close(10)

  ! ==================================================================
  ! Write statistics file with table line
  ! ==================================================================

  if ( iprintctl(5) ) then
     if ( solver .eq. 'SECO' ) then
        open(20,file='seco-tabline.out')
        write(20,9100) fu,cunorm,bdsnorm,nlnorm,f,cnorm,outiter, &
                       totiter,best_fu,best_cunorm,best_bdsnorm, &
                       best_nlnorm,best_f,best_cnorm,best_outit, &
                       secoinfo,inform,time
        close(20)
        
     else
        open(20,file='algencan-tabline.out')
        write(20,9200) fu,cunorm,f,cnorm,nlnorm,best_fu,best_cunorm,    &
                       best_f,best_cnorm,best_nlnorm,ncnorm,rsupn,      &
                       inform,solinfo,innfail,n,m,outiter,totiter,fcnt, &
                       nwcalls,nwtotit,msqcalls,msqtotit,time
        close(20)
     end if
  end if

  ! ==================================================================
  ! Deallocate global arrays in modules
  ! ==================================================================

  deallocate(x,l,u,best_x,best_lambda,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     inform = - 94
     subname = 'ALGENCAN'
     call reperr(inform,subname)
     return
  end if

  call deallocatearrays(allocerr)
  if ( allocerr .ne. 0 ) then
     inform = - 94
     subname = 'ALGENCAN'
     call reperr(inform,subname)
     return
  end if

  ! ==================================================================
  ! NON-EXECUTABLE STATEMENTS
  ! ==================================================================

 9000 format(/,1X,'Total CPU time in seconds = ',0P,F8.2)

 9100 format(1X,1P,D24.16,3(1X,1P,D7.1),1X,1P,D24.16,1X,1P,D7.1,2(1X,I12),&
           1X,1P,D24.16,3(1X,1P,D7.1),1X,1P,D24.16,1X,1P,D7.1,1X,I12,1X,  &
           I1,1X,I3,0P,F8.2)

 9200 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D8.1, &
           1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D8.1,   &
           1X,1P,D7.1,1X,1P,D7.1,1X,I3,1X,I1,1X,L1,1X,I6,1X,I6,1X,I3,    &
           1X,I7,1X,I7,1X,I2,1X,I7,1X,I7,1X,I7,0P,F8.2)

end subroutine algencan
