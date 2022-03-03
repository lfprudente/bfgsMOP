! ******************************************************************
! ******************************************************************

subroutine newtonkkt(n,xo,l,u,m,lambdao,equatn,linear,epsfeas, &
     epsopt,f,cnorm,nlnorm,iter,msqiter,accinfo,inform)

  use lsslvr, only: lssana,lssend,lssfac,lssset,lsssol
  use modmachconst
  use modalgconst
  use modalgparam, only: setalgparam,lsssubACC,sclsubACC,maxaccitgiven
  use modouttyp
  use problvls, only: ssetp,sevalobjc,sevalg,sevalgjac,sevaljac,sevalhl
  use probgiven, only: gjaccoded,hnnzlim,jcnnzlim

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: m,n
  integer,      intent(inout) :: accinfo,inform,iter,msqiter
  real(kind=8), intent(in)    :: epsfeas,epsopt
  real(kind=8), intent(inout) :: cnorm,f,nlnorm

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(inout) :: l(n),lambdao(m),u(n),xo(n)

  ! accinfo:
  !
  ! 0: KKT system solved.
  ! 1: Ignored constraints were violated.
  ! 2: After correcting the Lagrange multipliers signal, optimality
  !    was lost.
  ! 3: Maximum number of iterations reached.
  ! 4: Newton seems to be diverging.
  ! 5: Singular Jacobian.
  ! 6: Insufficient space to store the KKT linear system.
  ! 7: Insufficient double precision working space for linear solver.
  ! 8: Insufficient integer working space for linear solver.

  ! COMMON SCALARS
  real(kind=8) :: epsfeascopy,epsoptcopy

  ! LOCAL SCALARS
  integer      :: cind,col,divit,fun,hlnnz,i,itmp,j,jcnnz,k,lim, &
                  lssinfo,maxit,minp,mnop,mtot,nbds,nineq,ninn, &
                  nneigv,nnon,nrank,nsys,ntot,pind,row,sind,unnz,var,vind
  real(kind=8) :: cnorm0,cnormprev,cnormref,epsact,epsadd,fu,nlnorm0, &
                  nlnormprev,nlnormref,pval,sgnvio,val

  ! LOCAL ARRAYS
  logical      :: alreadyin(m+2*n),ind(m)
  character    :: constt(2*n),status(n)
  integer      :: consti(2*n),inn(m+3*n),inp(m+2*n), &
                  jcfun(jcnnzlim+4*n+m),jclen(m+2*n), &
                  jcvar(jcnnzlim+4*n+m),jcsta(m+2*n), &
                  slaind(m+2*n),udiag(2*m+5*n), &
                  ucol(hnnzlim+jcnnzlim+4*m+11*n), &
                  urow(hnnzlim+jcnnzlim+4*m+11*n)
  real(kind=8) :: adddiag(2*m+5*n),b(2*m+5*n),c(m+2*n),cu(m),g(n), &
                  jcval(jcnnzlim+4*n+m),lambda(m+2*n),nl(m+3*n), &
                  x(m+3*n),uval(hnnzlim+jcnnzlim+4*m+11*n),csalvo(m)

  ! COMMON BLOCKS
  common /minsqstopdata/ epsfeascopy,epsoptcopy
  save   /minsqstopdata/

  ! Bound constraints are seen as inequalities. A squared slack
  ! variable is added to each inequality constraint (including te ones
  ! given by the bound constraints). It means that the total number of
  ! variables goes up to 
  !
  !     ntot = 3 n + m 
  !
  ! and the total number of constraints goes up to
  !
  !     mtot = m + 2 n. 
  !
  ! The upper limit for the number of non-null elements is:
  ! (a) (Modified) Hessian of the Lagrangian : hnnzlim + 2n + m,
  ! (b) (Modified) Jacobian of the constraints : jcnnzlim + 4n + m,
  ! (c) Jacobian of the KKT system : (a) + (b) + 5n + 2m.

  epsfeascopy = epsfeas
  epsoptcopy  = epsopt

  ! ==================================================================
  ! PRESENTATION
  ! ==================================================================

  if ( iprintout .ge. 1 ) then
     write(* ,8000)
     write(10,8000)
  end if

  ! ==================================================================
  ! INITIALIZE
  ! ==================================================================

  iter    =  0
  divit   =  0
  msqiter =  0

  if ( maxaccitgiven .gt. 0 ) then
     maxit = maxaccitgiven
  else
     maxit = maxaccit
  end if

  epsact = sqrt( epsfeas )
  epsadd = macheps12

  call lssset(lsssubACC,sclsubACC,.true.,.false.)

  ! ==================================================================
  ! SET INITIAL POINT
  ! ==================================================================

  if ( minval( u(1:n) - l(1:n) ) .le. 0.0d0 ) then
     write(*,*) 'There is a fixed variable. This is a simple incovenience, '
     write(*,*) 'but the Acceleration Process is not ready to handle it.'
     stop
  end if

  x(1:n) = xo(1:n)
  lambda(1:m) = lambdao(1:m)

  ! ==================================================================
  ! STRUCTURES FOR NEW VARIABLES AND CONSTRAINTS
  ! ==================================================================

  ! Relate constraints to slacks

  nineq = 0
  do j = 1,m
     if ( .not. equatn(j) ) then
        nineq     = nineq + 1
        sind      = n + nineq
        slaind(j) = sind
     else
        slaind(j) = 0
     end if
  end do

  nbds = 0
  do i = 1,n
     if ( l(i) .gt. - 1.0d+20 ) then
        nbds = nbds + 1
        constt(nbds) = 'L'
        consti(nbds) =  i
        cind         = m + nbds
        sind         = n + nineq + nbds
        slaind(cind) = sind
     end if

     if ( u(i) .lt. 1.0d+20 ) then
        nbds = nbds + 1
        constt(nbds) = 'U'
        consti(nbds) =  i
        cind         = m + nbds
        sind         = n + nineq + nbds
        slaind(cind) = sind
     end if
  end do

  mtot = m + nbds
  ntot = n + nineq + nbds

  alreadyin(1:mtot) = .false.

  ! ==================================================================
  ! MAIN LOOP
  ! ==================================================================

100 continue

  ! ==================================================================
  ! SET FREE AND FIXED VARIABLES (FOR CURRENT ITERATION)
  ! ==================================================================
  
  ! Note that variables are being projected onto the box. Therefore,
  ! bounds are satisfied with zero tolerance.
  
  minp = 0
  mnop = 0

  ninn = 0
  nnon = 0

  ! Set free (F) and fixed (L or U) variables (a fixed variable will
  ! be fixed forever).

  do i = 1,n
     if ( x(i) .le. l(i) ) then
        status(i) = 'L'
        x(i) = l(i)

        nnon   = nnon + 1
        inn(i) = ntot + 1 - nnon

     else if ( x(i) .ge. u(i) ) then

        status(i) = 'U'
        x(i) = u(i)

        nnon   = nnon + 1
        inn(i) = ntot + 1 - nnon

     else
        status(i) = 'F'
        ninn   = ninn + 1
        inn(i) = ninn
     end if
  end do

  ! ==================================================================
  ! COMPUTE OBJECTIVE FUNCTION AND CONSTRAINTS
  ! ==================================================================
  
  ! Objective function and (regular) constraints
  
  call ssetp(n,x,inform)
  if ( inform .ne. 0 ) return

  call sevalobjc(n,x,f,fu,m,c,cu,inform)
  if ( inform .ne. 0 ) return

  csalvo(1:m) = c(1:m)

  ! Compute constraints based on bound constraints

  do i = 1,nbds
     cind = m + i
     vind = consti(i)

     if ( constt(i) .eq. 'L' ) then
        c(cind) = l(vind) - x(vind)
     else ! if ( constt(i) .eq. 'U' ) then
        c(cind) = x(vind) - u(vind)
     end if
  end do

  ! Set slacks's values

  do j = 1,mtot
     sind = slaind(j)
     if ( sind .ne. 0 ) then
        x(sind) = sqrt( 2.0d0 * max( 0.0d0, - c(j) ) )
     end if
  end do

  ! Force complementarity with respect to the original (inequality)
  ! constraints (it requires looking at the inequalities without the
  ! slacks' effect). Values for the multipliers associated to the
  ! bound constraints will be computed after computing the gradient of
  ! the Lagrangian.

  do j = 1,m
     if ( .not. equatn(j) ) then
        if ( c(j) .lt. - epsfeas ) then
           lambda(j) = 0.0d0
        end if
     end if
  end do

  ! ==================================================================
  ! COMPUTE FEASIBILITY AND SAVE NORMS TO CHECK IMPROVEMENT
  ! ==================================================================
  
  if ( iter .eq. 0 ) then
     nlnormprev = bignum
     cnormprev  = bignum
  else
     nlnormprev = nlnorm
     cnormprev  = cnorm
  end if

  ! Note that, although some inequality constraints are being ignored
  ! to compute the Newton step, feasibility considering the whole set
  ! of constraints is being computed.

  ! Violation of original constraints (bound constraints are satisfied
  ! by the projection onto the box, see above). If original
  ! constraints are satisfied, by the way the slacks' values are
  ! computed, modified constraints are satisfied too.

  cnorm = 0.0d0
  do j = 1,m
     if ( equatn(j) ) then
        cnorm = max( cnorm, abs( c(j) ) )
     else
        cnorm = max( cnorm,      c(j)   )
     end if
  end do

  if ( iter .eq. 0 ) then
     cnorm0 = cnorm
  end if

  ! ==================================================================
  ! SET ACTIVE AND INACTIVE CONSTRAINTS (FOR CURRENT ITERATION)
  ! ==================================================================
  
  ! Regular constraints and their slacks
  
  do j = 1,m
     if ( equatn(j) ) then
        ! Active equality constraint
        minp   = minp + 1
        inp(j) = minp
     else
        sind = slaind(j)

        if ( alreadyin(j) .or. c(j) .ge. - epsact ) then
           ! Active inequality constraint and its slack
           alreadyin(j) = .true.

           minp      = minp + 1
           inp(j)    = minp

           ninn      = ninn + 1
           inn(sind) = ninn
        else
           ! Inactive inequality constraint and its slack
           mnop      = mnop + 1
           inp(j)    = mtot + 1 - mnop

           nnon      = nnon + 1
           inn(sind) = ntot + 1 - nnon
        end if
     end if
  end do

  ! Bound constraints and their slacks

  do i = 1,nbds
     cind = m + i
     vind = consti(i)

     sind = slaind(cind)

     if ( status(vind) .eq. 'F' .and. ( alreadyin(cind) .or. c(cind) .ge. - epsact ) ) then
        ! Active bound constraint and its slack
        alreadyin(cind) = .true.

        minp      = minp + 1
        inp(cind) = minp

        ninn      = ninn + 1
        inn(sind) = ninn
     else
        ! Inactive bound constraint and its slack
        mnop      = mnop + 1
        inp(cind) = mtot + 1 - mnop

        nnon      = nnon + 1
        inn(sind) = ntot + 1 - nnon
     end if
  end do

  nsys = ninn + minp

  ! ==================================================================
  ! ADD SLACKS' EFFECT TO THE CONSTRAINTS
  ! ==================================================================

  ! Slacks' effect

  do j = 1,mtot
     sind = slaind(j)
     if ( sind .ne. 0 ) then
        c(j) = c(j) + 0.5d0 * x(sind) ** 2
     end if
  end do

  ! ==================================================================
  ! COMPUTE FIRST DERIVATIVES
  ! ==================================================================

  ! Gradient of the objective function and Jacobian of the constraints

  if ( gjaccoded ) then

     call sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform)
     if ( inform .ne. 0 ) return

  else

     ! Gradient of the objective function

     call sevalg(n,x,g,inform)
     if ( inform .ne. 0 ) return

     if ( m .eq. 0 ) then
        jcnnz = 0
        
     else
        ind(1:m) = .true.
        call sevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,jcnnzlim,inform)
        if ( inform .ne. 0 ) return

        jcnnz = sum( jclen(1:m) )

        do j = 1,m
           jcfun(jcsta(j):jcsta(j)+jclen(j)-1) = j
        end do
     end if

  end if

  ! First derivative related to the slacks of the regular constraints

  do j = 1,m
     sind = slaind(j)
     
     if ( sind .ne. 0 ) then
        jcnnz = jcnnz + 1
        
        jcfun(jcnnz) = j
        jcvar(jcnnz) = sind
        jcval(jcnnz) = x(sind)
     end if
  end do
  
  call coo2csr(m,jcnnz,jcfun,jcvar,jcval,jclen,jcsta)

  ! Bound constraints with the slacks effect

  do i = 1,nbds
     cind = m + i
     sind = slaind(cind)
     vind = consti(i)

     jcsta(cind) = jcnnz + 1
     jclen(cind) = 2

     jcvar(jcnnz+1) = vind
     if ( constt(i) .eq. 'L' ) then
        jcval(jcnnz+1) = - 1.0d0
     else
        jcval(jcnnz+1) =   1.0d0
     end if

     jcvar(jcnnz+2) = sind
     jcval(jcnnz+2) = x(sind)

     jcnnz = jcnnz + jclen(cind)
  end do

  ! Compute gradient of the Lagrangian

  nl(1:n) = g(1:n)
  nl(n+1:ntot) = 0.0d0

  ! Add effect of the Jacobian of regular constraints

  do j = 1,m
     if ( lambda(j) .ne. 0.0d0 ) Then
        do i = jcsta(j),jcsta(j) + jclen(j) - 1
           nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
        end do
     end if
  end do

  ! Set Lagrange multipliers related to bound constraints (forcing
  ! complementarity)

  do i = 1,nbds
     cind = m + i
     vind = consti(i)
     if ( status(vind) .ne. 'F' ) then
        if ( constt(i) .eq. status(vind) ) then
           if ( constt(i) .eq. 'L' ) then
              lambda(cind) =   nl(vind)
           else
              lambda(cind) = - nl(vind)
           end if
        else
           lambda(cind) = 0.0d0
        end if
     else
        lambda(cind) = 0.0d0
     end if
  end do

  ! Effect of Jacobian of bound constraints

  do j = m + 1,m + nbds
     if ( lambda(j) .ne. 0.0d0 ) Then
        do i = jcsta(j),jcsta(j) + jclen(j) - 1
           nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
        end do
     end if
  end do

  ! Gradient of the Lagrangian norm

  nlnorm = maxval( abs( nl(1:ntot) ) )

  if ( iter .eq. 0 ) then
     nlnorm0 = nlnorm
  end if

  ! ==================================================================
  ! WRITE INFORMATION OF THE CURRENT POINT
  ! ==================================================================

  if ( iprintout .ge. 1 ) then
     write(* ,8010) iter,f,cnorm,nlnorm
     write(10,8010) iter,f,cnorm,nlnorm
     ! write(* ,8010) iter,f,cnorm,nlnormtmp,nlnorm
     ! write(10,8010) iter,f,cnorm,nlnormtmp,nlnorm
  end if

  ! ==================================================================
  ! TEST STOPPING CRITERIA
  ! ==================================================================

  if ( cnorm .le. epsfeas .and. nlnorm .le. epsopt ) then
     ! THE POINT SATISFIES FEASIBILITY AND OPTIMALITY

     ! (Lagrange multipliers signal must be checked)

     go to 400
  end if

  if ( iter .ge. maxit ) then
     ! MAXIMUM NUMBER OF ITERATIONS EXCEEDED

     accinfo = 3

     if ( iprintout .ge. 1 ) then
        write(* ,9030)
        write(10,9030)
     end if

     go to 500
  end if

  cnormref  = max(cnorm0 ,max(cnormprev ,max(epsfeas, 1.0d+01)))
  nlnormref = max(nlnorm0,max(nlnormprev,max(epsopt,  1.0d+01)))

  if ( cnorm .gt. cnormref  .or. nlnorm .gt. nlnormref .or. &
       cnorm .eq. cnormprev .or. nlnorm .eq. nlnormprev ) then
     ! IT SEEMS TO BE DIVERGING

     divit = divit + 1

     if ( divit .ge. 3 ) then

        accinfo = 4

        if ( iprintout .ge. 1 ) then
           write(* ,9040)
           write(10,9040)
        end if

        go to 500

     end if

  else
     divit = 0
  end if

  ! ==================================================================
  ! DO AN ITERATION
  ! ==================================================================

  iter = iter + 1

  ! ==================================================================
  ! COMPUTE SECOND DERIVATIVES
  ! ==================================================================

  ! Hessian of the Lagrangian

  call sevalhl(n,x,m,lambda,urow,ucol,uval,hlnnz,hnnzlim,inform)
  if ( inform .ne. 0 ) return

  ! Second derivatives related to slacks of the regular constraints
  ! and bound constraints

  do j = 1,mtot
     sind = slaind(j)

     if ( sind .ne. 0 ) then
        hlnnz = hlnnz + 1

        urow(hlnnz) = sind
        ucol(hlnnz) = sind
        uval(hlnnz) = lambda(j)
     end if
  end do

  ! ==================================================================
  ! ASSEMBLE THE JACOBIAN OF THE KKT SYSTEM
  ! ==================================================================

  unnz = 0

  udiag(1:nsys) = 0

  ! Hessian of the Lagrangian

  do i = 1,hlnnz
     if ( urow(i) .ge. ucol(i) ) then

        row = inn(urow(i))
        col = inn(ucol(i))
        val = uval(i)

        if ( row .le. ninn .and. col .le. ninn ) then
           if ( val .ne. 0.0d0 ) then
              ! A(row,col) = A(row,col) + val
              unnz = unnz + 1
              urow(unnz) = row
              ucol(unnz) = col
              uval(unnz) = val
              if ( row .eq. col ) udiag(row) = unnz
           end if
        end if

     end if
  end do

  ! Jacobian of the constraints

  do j = 1,mtot
     do i = jcsta(j),jcsta(j) + jclen(j) - 1

        fun = inp(j)
        var = inn(jcvar(i))
        val = jcval(i)

        if ( var .le. ninn .and. fun .le. minp ) then
           if ( val .ne. 0.0d0 ) then
              ! A(fun+ninn,var) = A(fun+ninn,var) + val
              unnz = unnz + 1
              urow(unnz) = fun + ninn
              ucol(unnz) = var
              uval(unnz) = val
           end if
        end if

     end do
  end do

  do i = 1,nsys
     if ( udiag(i) .eq. 0 ) then
        unnz = unnz + 1
        urow(unnz) = i
        ucol(unnz) = i
        uval(unnz) = 0.0d0

        udiag(i) = unnz
     end if
  end do

  ! ==================================================================
  ! ANALYSE SPARSITY PATTERN
  ! ==================================================================

  call lssana(nsys,unnz,urow,ucol,uval,udiag,lssinfo)

  if ( lssinfo .eq. 6 ) then
     ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

     accinfo = 6
     go to 500

  end if

  ! ==================================================================
  ! SOLVE THE NEWTONIAN SYSTEM
  ! ==================================================================

  ! ==================================================================
  ! COMPUTE REGULARIZATION
  ! ==================================================================

200 continue

  adddiag(1:ninn) = epsadd
  adddiag(ninn+1:nsys) = - macheps12

  ! ==================================================================
  ! FACTORIZE THE JACOBIAN OF THE NEWTONIAN SYSTEM
  ! ==================================================================

  call lssfac(nsys,unnz,urow,ucol,uval,udiag,adddiag,pind,pval, &
       nneigv,nrank,lssinfo)

  if ( lssinfo .eq. 0 .or. lssinfo .eq. 1 ) then

     if ( nneigv .ne. minp ) then
        ! WRONG INERTIA (SEE NOCEDAL AND WRIGHT)

        ! Lemma 16.3 [pg. 447]: Suppose the the Jacobian of the
        ! constraints has full rank and that the reduced Hessian
        ! Z^T H Z is positive definite. Then the Jacobian of the
        ! KKT system has ninn positive eigenvalues, minp negative
        ! eigenvalues, and no zero eigenvalues.

        ! Note that at this point we know that the matrix has no
        ! zero eigenvalues. nneigv gives the number of negative
        ! eigenvalues.

        epsadd = max( macheps12, epsadd * 10.0d0 )

        if ( iprintout .ge. 1 ) then
           itmp = ninn+minp-nneigv
           write(* ,8090) ninn,minp,itmp,nneigv,epsadd
           write(10,8090) ninn,minp,itmp,nneigv,epsadd
        end if

        go to 200
     end if

  else if ( lssinfo .eq. 2 ) then
     ! SINGULAR JACOBIAN

     epsadd = max( macheps12, epsadd * 10.0d0 )

     if ( iprintout .ge. 1 ) then
        write(* ,8080) epsadd
        write(10,8080) epsadd
     end if

     if ( epsadd .le. 1.0d+20 ) then
        go to 200
     end if

     accinfo = 5

     if ( iprintout .ge. 1 ) then
        write(* ,9050)
        write(10,9050)
     end if

     go to 500

  else if ( lssinfo .eq. 6 ) then
     ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

     accinfo = 6
     go to 500

  else if ( lssinfo .eq. 7 ) then
     ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

     accinfo = 7
     go to 500

  else ! if ( lssinfo .eq. 8 ) then
     ! INSUFFICIENT INTEGER WORKING SPACE

     accinfo = 8
     go to 500

  end if

  ! ==================================================================
  ! SOLVE TRIANGULAR SYSTEMS
  ! ==================================================================

  ! SET RHS

  do i = 1,ntot
     if ( inn(i) .le. ninn ) then
        b(inn(i)) = - nl(i)
     end if
  end do

  do j = 1,mtot
     if ( inp(j) .le. minp ) then
        b(ninn+inp(j)) = - c(j)
     end if
  end do

  ! SOLVE THE EQUATIONS AND DEALLOCATE MEMORY

  call lsssol(nsys,b)

  call lssend()

  ! ==================================================================
  ! UPDATE x AND lambda
  ! ==================================================================

  do i = 1,ntot
     if ( inn(i) .le. ninn ) then
        x(i) = x(i) + b(inn(i))
     end if
  end do

  do j = 1,mtot
     if ( inp(j) .le. minp ) then
        lambda(j) = lambda(j) + b(ninn+inp(j))
     end if
  end do

  ! ==================================================================
  ! ITERATE
  ! ==================================================================

  go to 100

  ! ==================================================================
  ! END OF MAIN LOOP
  ! ==================================================================

  ! ==================================================================
  ! CHECK LAGRANGE MULTIPLIERS SIGNAL RELATED TO INEQUALITIES AND
  ! BOUND CONSTRAINTS
  ! ==================================================================

400 continue

  sgnvio = min( minval( lambda(1:m), .not. equatn(1:m) ), minval( lambda(m+1:mtot) ) )

  if ( sgnvio .ge. 0.0d0 ) then
     accinfo = 0

     if ( iprintout .ge. 1 ) then
        write(* ,9000)
        write(10,9000)
     end if

     go to 500
  end if

  ! TRY TO CORRECT MULTIPLIERS SIGNAL

  if ( iprintout .ge. 1 ) then
     write(* ,8040)
     write(10,8040)
  end if

  call minsq(n,g,m,equatn,csalvo,nbds,lambda,constt,consti,status, &
       jcsta,jclen,jcvar,jcval,epsopt,msqiter,inform)
  if ( inform .ne. 0 ) return

  ! Compute signal violation

  sgnvio = min( minval( lambda(1:m), .not. equatn(1:m) ), minval( lambda(m+1:mtot) ) )

  ! Compute optimality

  nl(1:n) = g(1:n)
  nl(n+1:ntot) = 0.0d0

  do j = 1,mtot
     if ( lambda(j) .ne. 0.0d0 ) then
        do i = jcsta(j),jcsta(j) + jclen(j) - 1
           nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
        end do
     end if
  end do

  nlnorm = maxval( abs( nl(1:ntot) ) )

  if ( iprintout .ge. 1 ) then
     write(* ,8050) sgnvio,nlnorm
     write(10,8050) sgnvio,nlnorm
  end if

  if ( nlnorm .le. epsopt ) then
     accinfo = 0

     if ( iprintout .ge. 1 ) then
        write(* ,9000)
        write(10,9000)
     end if

     go to 500

  else
     accinfo = 2

     if ( iprintout .ge. 1 ) then
        write(* ,9020)
        write(10,9020)
     end if

     go to 500
  end if

  ! ==================================================================
  ! SET SOLUTION
  ! ==================================================================

500 continue

  xo(1:n) = x(1:n)
  lambdao(1:m) = lambda(1:m)

!!$  write(*,*) 'solucion:'
!!$  write(*,*) 'x:'
!!$  do i = 1,n
!!$     write(*,*) i,x(i)
!!$  end do
!!$  write(*,*) 'lambda:'
!!$  do i = 1,m
!!$     write(*,*) i,lambda(i)
!!$  end do

! ==================================================================
! NON-EXECUTABLE STATEMENTS
! ==================================================================

 8000 format(/,' NEWTON-KKT scheme in action!')
 8010 format(/,' NEWTON-KKT Iteration',42X,' = ',5X,I6,                 &
             /,' Objective function value',38X,' = ',1PD11.4,           &
!            /,' Maximal violation of selected nearly-active',          &
!              ' constraints',7X,' = ',1PD11.4,                         &
             /,' Maximal violation of constraints',30X,' = ',1PD11.4,   &
             /,' Sup-norm of the gradient of the Lagrangian',20X,' = ', &
                 1PD11.4)
!    /,' Sup-norm of the gradient of the Lagrangian',8X,1PD11.4,&
!        1X,' = ',1PD11.4)
 8020 format(/,' Maximal violation of constraints',30X,' = ',1PD11.4)
 8030 format(/,' Maximal violation of Lagrange multipliers', &
               ' non-negativity',6X,' = ',1PD11.4)
 8040 format(/,' GENCAN is being called to find the right', &
               ' multipliers.')
 8050 format(/,' Maximal violation of Lagrange multipliers',            &
               ' non-negativity',6X,' = ',1PD11.4,                      &
             /,' Sup-norm of the gradient of the Lagrangian',20X,' = ', &
                 1PD11.4)

 8080 format(/,' Singular Jacobian.', &
               ' epsadd was increased to ',1PD11.4)
 8090 format(/,' Wrong Jacobian inertia. ',       &
             /,' Desired POS = ',I7,' NEG = ',I7, &
             /,' Actual  POS = ',I7,' NEG = ',I7, &
             /,' epsadd was increased to ',1PD11.4)

 9000 format(/,' Flag of NEWTON-KKT = KKT system solved!')
!9010 format(/,' Flag of NEWTON-KKT = Ignored constraints were', &
!              ' violated.')
!9011 format(/,' Flag of NEWTON-KKT = Ignored bound constraints were', &
!              ' violated.')
 9020 format(/,' Flag of NEWTON-KKT = After correcting the Lagrange', &
               ' multipliers signal,',/,' optimality was lost.')
 9030 format(/,' Flag of NEWTON-KKT = Maximum of iterations reached.')
 9040 format(/,' Flag of NEWTON-KKT = Newton can not make further', &
               ' progress.')
 9050 format(/,' Flag of NEWTON-KKT = Singular Jacobian.')

end subroutine newtonkkt

! ******************************************************************
! ******************************************************************

subroutine minsq(n,g,m,equatn,csalvo,nbds,lambda,constt,consti, &
     status,jcsta,jclen,jcvar,jcval,epsopt,iter,inform)

  use modalgparam
  use problvlu, only: setproblvlu,rmfixv
  use modminsq

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,n,nbds
  integer, intent(inout) :: inform,iter
  real(kind=8), intent(in) :: epsopt

  ! ARRAY ARGUMENTS
  character,    intent(in)    :: constt(nbds),status(n)
  logical,      intent(in)    :: equatn(m)
  integer,      intent(in)    :: consti(nbds),jcsta(m+nbds),jclen(m+nbds),jcvar(*)
  real(kind=8), intent(in)    :: csalvo(m),g(n),jcval(*)
  real(kind=8), intent(inout) :: lambda(m+nbds)

  ! The problem to be solved is
  !
  !     min 0.5 || b + A^T lambda ||_2^2
  !
  !     subject to
  !
  !     lambda_i >= 0, i = m+1,...,m+nabds.
  !
  ! Columns a_i of A^T are multiplied by
  !
  ! an_i = 1 / max( 1, ||a_i||_infty ).
  !
  ! So, if we define D = diag(an_1,...,an_{m+nabds}), the problem
  ! can be rewrittem as
  !
  !     min 0.5 || b + A^T D D^{-1} lambda ||_2^2
  !
  !     subject to
  !
  !     an_i lambda_i >= 0, i = m+1,...,m+nabds.
  !
  ! Subtituting an_i lambda_i by lambda_i the problem to be solved
  ! becomes
  !
  !     min 0.5 || b + A^T D lambda ||_2^2
  !
  !     subject to
  !
  !     lambda_i >= 0, i = m+1,...,m+nabds.

  ! COMMON SCALARS
  real(kind=8) :: epsfeascopy,epsoptcopy

  ! LOCAL SCALARS
  character(len=2) :: innslvrtmp
  logical          :: rmfixvtmp,useustptmp
  integer          :: cind,geninfo,i,j,k,maxit,vind
  real(kind=8)     :: dum5,dum6,eps,msqf,msqnlpsupn

  ! LOCAL ARRAYS
  logical      :: dum2(1),dum3(1)
  integer      :: ind(m+nbds)
  real(kind=8) :: dum1(1),dum4(1),l(m+nbds), &
                  scol(m+nbds),u(m+nbds),x(m+nbds)

  ! COMMON BLOCKS
  common /minsqstopdata/ epsfeascopy,epsoptcopy
  save   /minsqstopdata/

  ! RHS

  b(1:n) = g(1:n)

  ! Matrix

  k = 0
  annz = 0

  do j = 1,m
     ! Only equality and nearly-active inequality constraints
     if ( equatn(j) .or. csalvo(j) .ge. - epsfeascopy ) then
        k = k + 1

        ind(k) = j

        do i = jcsta(j),jcsta(j) + jclen(j) - 1
           if ( jcvar(i) .le. n ) then
              annz = annz + 1
              acol(annz) = k
              arow(annz) = jcvar(i)
              aval(annz) = jcval(i)
           end if
        end do

        x(k) = lambda(j)

        if ( equatn(j) ) then
           l(k) = - 1.0d+20
           u(k) =   1.0d+20
        else
           l(k) =   0.0d0
           u(k) =   1.0d+20
        end if
     end if
  end do

  do j = 1,nbds
     vind = consti(j)

     ! Only active bound constraints
     if ( constt(j) .eq. status(vind) ) then
        k = k + 1

        cind = m + j

        ind(k) = cind

        do i = jcsta(cind),jcsta(cind) + jclen(cind) - 1
           if ( jcvar(i) .le. n ) then
              annz = annz + 1
              acol(annz) = k
              arow(annz) = jcvar(i)
              aval(annz) = jcval(i)
           end if
        end do

        x(k) = lambda(cind)

        l(k) = 0.0d0
        u(k) = 1.0d+20
     end if
  end do

  ! Dimensions

  nrows = n
  ncols = k

  ! Columns scaling

  scol(1:ncols) = 1.0d0

  do i = 1,annz
     scol(acol(i)) = max( scol(acol(i)), abs( aval(i) ) )
  end do

  scol(1:ncols) = 1.0d0 / scol(1:ncols)

  aval(1:annz) = aval(1:annz) * scol(acol(1:annz))

  x(1:ncols) = x(1:ncols) / scol(1:ncols)

  ! Call the solver

  call setalgparam(val_innercall = .true.)

  rmfixvtmp  = rmfixv
  innslvrtmp = innslvr
  useustptmp = useustp

  call setproblvlu(val_rmfixv  = .false.)
  call setalgparam(val_innslvr = 'TN')
  call setalgparam(val_useustp = .true.)

  maxit = 200
  eps = 1.0d-16

  call gencan(ncols,x,l,u,0,dum1,dum2,dum3,dum4,0.0d0,eps,0.0d0,0.0d0, &
       maxit,iter,msqf,msqnlpsupn,dum5,dum6,geninfo,inform)

  call setalgparam(val_innercall = .false.)

  call setproblvlu(val_rmfixv  = rmfixvtmp)
  call setalgparam(val_innslvr = innslvrtmp)
  call setalgparam(val_useustp = useustptmp)

  ! Set to zero multipliers associated with non-active constraints and
  ! copy unscaled solution

  lambda(1:m+nbds) = 0.0d0

  lambda(ind(1:ncols)) = x(1:ncols) * scol(1:ncols)

end subroutine minsq

! ******************************************************************
! ******************************************************************

subroutine minsqf(n,x,f,inform)

  use modminsq

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: n
  integer,      intent(inout) :: inform
  real(kind=8), intent(out)   :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! LOCAL SCALARS
  integer :: i

  ! LOCAL ARRAYS
  real(kind=8) :: p(nrows)

  p(1:nrows) = b(1:nrows)

  do i = 1,annz
     p(arow(i)) = p(arow(i)) + x(acol(i)) * aval(i)
  end do

  f = 0.5d0 * sum( p(1:nrows) ** 2 )

  f = 1.0d+08 * f

end subroutine minsqf

! ******************************************************************
! ******************************************************************

subroutine minsqg(n,x,g,inform)

  use modminsq

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: g(n)

  ! LOCAL SCALARS
  integer :: i

  ! LOCAL ARRAYS
  real(kind=8) :: p(nrows)

  p(1:nrows) = b(1:nrows)

  do i = 1,annz
     p(arow(i)) = p(arow(i)) + x(acol(i)) * aval(i)
  end do

  g(1:ncols) = 0.0d0

  do i = 1,annz
     g(acol(i)) = g(acol(i)) + p(arow(i)) * aval(i)
  end do

  g(1:ncols) = 1.0d+08 * g(1:ncols)

end subroutine minsqg

! ******************************************************************
! ******************************************************************

subroutine minsqhp(n,x,p,hp,goth,inform)

  use modminsq

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: goth
  integer, intent(in)    :: n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)    :: p(n),x(n)
  real(kind=8), intent(inout) :: hp(n)

  ! LOCAL SCALARS
  integer :: i

  ! LOCAL ARRAYS
  real(kind=8) :: tmp(nrows)

  tmp(1:nrows) = 0.0d0

  do i = 1,annz
     tmp(arow(i)) = tmp(arow(i)) + p(acol(i)) * aval(i)
  end do

  hp(1:ncols) = 0.0d0

  do i = 1,annz
     hp(acol(i)) = hp(acol(i)) + tmp(arow(i)) * aval(i)
  end do

  hp(1:ncols) = 1.0d+08 * hp(1:ncols)

end subroutine minsqhp

! ******************************************************************
! ******************************************************************

function minsqstop(n,x,inform)

  use modmachconst
  use modminsq

  implicit none

  ! FUNCTION TYPE
  logical :: minsqstop

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! COMMON SCALARS
  real(kind=8) :: epsfeascopy,epsoptcopy

  ! LOCAL SCALARS
  integer      :: i
  real(kind=8) :: pnorm

  ! LOCAL ARRAYS
  real(kind=8) :: p(nrows)

  ! COMMON BLOCKS
  common /minsqstopdata/ epsfeascopy,epsoptcopy
  save   /minsqstopdata/

  p(1:nrows) = b(1:nrows)

  do i = 1,annz
     p(arow(i)) = p(arow(i)) + x(acol(i)) * aval(i)
  end do

  pnorm = maxval( abs( p(1:nrows) ) )

  minsqstop = .false.
  if ( pnorm .le. epsoptcopy ) then
     minsqstop = .true.
  end if

end function minsqstop
