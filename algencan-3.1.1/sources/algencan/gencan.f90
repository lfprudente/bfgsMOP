! ******************************************************************
! ******************************************************************

subroutine gencan(n,x,l,u,m,lambda,equatn,linear,rho,epsfeas, &
     epsopt,efstain,eostain,maxit,iter,f,gpsupn,cnorm,cnormu, &
     geninfo,inform)

  use modmachconst
  use modalgconst
  use modalgparam
  use modprobdata, only: nbds
  use modouttyp
  use modrspace
  use moditetyp
  use modsydat
  use problvls, only: ssetp
  use problvlv, only: reperr,fcnt

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: m,maxit,n
  integer,      intent(inout) :: inform
  integer,      intent(out)   :: geninfo,iter
  real(kind=8), intent(in)    :: epsfeas,epsopt,efstain,eostain
  real(kind=8), intent(out)   :: cnorm,cnormu,f,gpsupn

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m),rho(m)
  real(kind=8), intent(inout) :: l(n),u(n),x(n)

  ! Solves the box-constrained minimization problem
  !
  !        Minimize f(x) subject to l <= x <= u
  !
  ! using the method described in
  !
  ! E. G. Birgin and J. M. Martinez, ''Large-scale active-set box-
  ! constrained optimization method with spectral projected
  ! gradients'', Computational Optimization and Applications 23, pp.
  ! 101-125, 2002.

  ! geninfo:
  !
  ! 0: Small continuous-projected-gradient norm
  ! 1: Maximum number of gencan iterations reached
  ! 3: Lack of progress
  ! 6: Unbounded objective function?

  ! LOCAL SCALARS
  logical           :: forceoi,memfail,samep,vustop
  integer           :: allocerr,cgcnt,cginfo,cgiter,cgmaxit,chcnt,i,   &
                       innit,itnp,itnfp,itngp,itnxp,lsinfo,nind,       &
                       nindprev,nwcnt,lvfit,rbdnnz,trcnt,trinfo,xind
  real(kind=8)      :: acgeps,adsupn,amax,bcgeps,cgeps,cgdel,dsupn,    &
                       fbest,fplus,gpeucnb,gieucn,gisupn,gpeucn,gpi,   &
                       gpsupnb,lamspg,mslamb,newtrdel,ssupn,tmp,trdel, &
                       tsmall,ysupn,xeucn,xsupn,maxelem

  ! LOCAL ARRAYS
  integer           :: rbdind(n)
  character         :: rbdtype(n)
  character(len=15) :: subname
  character(len=2)  :: inniter
  real(kind=8)      :: d(n),g(n),gplus(n),mindiag(n),xplus(n)

  ! DATA BLOCKS
  character(len=21) :: ittext(0:4)
  integer           :: FSTORDPNT,SECORDPNT,SMALLSTEP,TCLSBNDRY, &
                       TOOSMLTRR,UNBOUNDED,UNDEFSTEP

  data ittext(0) /'(Initial point)     '/
  data ittext(1) /'(SPG step)          '/
  data ittext(2) /'(Truncated Newton)  '/
  data ittext(3) /'(Newton w/MS-system)'/
  data ittext(4) /'(Trust region)      '/

  data UNBOUNDED  /2/
  data SMALLSTEP  /3/
  data TOOSMLTRR  /4/
  data UNDEFSTEP  /5/
  data FSTORDPNT  /6/
  data SECORDPNT  /7/
  data TCLSBNDRY  /8/

  ! EXTERNAL FUNCTIONS
  logical :: sstop

  ! EXTERNAL SUBROUTINES
  external :: sevalal

  ! ==================================================================
  ! Allocate global arrays in modules
  ! ==================================================================

  call allocategencanarrays(n,allocerr)

  if ( allocerr .ne. 0 ) then
     inform = - 93
     subname = 'GENCAN'
     call reperr(inform,subname)
     return
  end if

  ! ==================================================================
  ! Initialization
  ! ==================================================================

  ! Set some initial values:

  ! best functional and gradient-norm values
  fbest   = bignum
  gpsupnb = bignum
  gpeucnb = bignum

  ! to record a memory failure in the direct solver
  memfail = .false.

  ! for the first inner-to-the-face minimization algorithm
  inniter = 'TN'

  ! for testing lack of progress checking f, g and x
  itnfp = 0
  itngp = 0
  itnxp = 0

  ! to force a leaving-face iteration when lack of progress is
  ! detected within a face
  forceoi = .false.

  ! for calculating More-Sorensen's direction
  mslamb =  0.0d0

  ! for counting number of iterations as well as inner-to-the-face
  ! and leaving-face iterations
  iter  = 0
  cgcnt = 0
  nwcnt = 0
  trcnt = 0
  chcnt = 0
  lvfit = 0
  innit = 0

  ! just to print "Initial point" in the first ouput
  ittype = 0

  ! Print problem information

  if ( iprintinn .ge. 1 ) then
     write(* ,1000) n,epsopt
     write(10,1000) n,epsopt
  end if

  if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
     write(* ,1010) nprint,(l(i),i=1,nprint)
     write(* ,1020) nprint,(u(i),i=1,nprint)
     write(* ,1030) nprint,(x(i),i=1,nprint)

     write(10,1010) nprint,(l(i),i=1,nprint)
     write(10,1020) nprint,(u(i),i=1,nprint)
     write(10,1030) nprint,(x(i),i=1,nprint)
  end if

  ! Project initial guess. If the initial guess is infeasible,
  ! projection puts it into the box.

  do i = 1,n
     x(i) = max( l(i), min( x(i), u(i) ) )
  end do

  ! Compute function and gradient at the initial point

  call ssetp(n,x,inform)
  if ( inform .ne. 0 ) go to 500

  call sevalal(n,x,m,lambda,rho,equatn,linear,f,inform)
  if ( inform .ne. 0 ) go to 500

  call sevalnal(n,x,m,lambda,rho,equatn,linear,g,inform)
  if ( inform .ne. 0 ) go to 500

  ! For feasibility problems, compute scaled and unscaled feasibilities

  if ( ignoref ) then
     call sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)
     if ( inform .ne. 0 ) go to 500

  else
     cnorm  = 0.0d0
     cnormu = 0.0d0
  end if

  ! Compute x Euclidian and sup norms, continuous-project-gradient
  ! Euclidian and sup norms, internal gradient Euclidian norm. Set
  ! nind as the number of free variables and save in array ind their
  ! identifiers. Also set variable sameface (same face) indicating
  ! that the "previous iterate" does not belong to the current face.

  sameface = .false.

  nind   = 0
  xsupn  = 0.0d0
  xeucn  = 0.0d0
  gpsupn = 0.0d0
  gpeucn = 0.0d0
  gisupn = 0.0d0
  gieucn = 0.0d0
  do i = 1,n
     if ( abs( x(i) ) .ge. xsupn ) then
        xsupn = abs( x(i) )
        xind = i
     end if
     xeucn = xeucn + x(i) ** 2
     gpi = x(i) - g(i)
     if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
        gpi = - g(i)
     else
        gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
     end if
     gpsupn = max( gpsupn, abs( gpi ) )
     gpeucn = gpeucn + gpi ** 2
     if ( x(i) .gt. l(i) + macheps23 * max(1.0d0,abs(l(i))) .and. &
          x(i) .lt. u(i) - macheps23 * max(1.0d0,abs(u(i))) ) then
        gisupn    = max( gisupn, abs( gpi ) )
        gieucn    = gieucn + gpi ** 2
        nind      = nind + 1
        ind(nind) = i
     end if
  end do
  xeucn  = sqrt( xeucn  )
  gpeucn = sqrt( gpeucn )
  gieucn = sqrt( gieucn )

  ! Compute trut-region radius, in case betra is used as inner solver

  trdel    = max( trdelmin, trdelini * max( 1.0d0, xeucn ) )
  newtrdel = 0.0d0

  ! Initial spectral steplength

  ! Compute a small step and set the point at which the auxiliary
  ! gradient will be computed

  if ( gpsupn .ne. 0.0d0 ) then
     tsmall = macheps12 * max( 1.0d0, xsupn / gpsupn )
  else
     tsmall = 0.0d0
  end if

  do i = 1,n
     gpi = x(i) - g(i)
     if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
        gpi = - g(i)
     else
        gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
     end if
     s(i) = x(i) + tsmall * gpi
  end do

  ! Compute the gradient at the auxiliary point

  call ssetp(n,s,inform)
  if ( inform .ne. 0 ) go to 500

  call ievalnal(n,s,m,lambda,rho,equatn,linear,.false.,y,inform)
  if ( inform .ne. 0 ) go to 500

  call ssetp(n,x,inform)
  if ( inform .ne. 0 ) go to 500

  ! Compute s = x_{1/2} - x_0 and y = g_{1/2} - g_0

  sts = 0.0d0
  sty = 0.0d0
  ssupn = 0.0d0
  ysupn = 0.0d0
  yeucn = 0.0d0
  do i = 1,n
     s(i) = s(i) - x(i)
     y(i) = y(i) - g(i)
     sts  = sts + s(i) ** 2
     sty  = sty + s(i) * y(i)
     ssupn = max( ssupn, abs( s(i) ) )
     ysupn = max( ysupn, abs( y(i) ) )
     yeucn = yeucn + y(i) ** 2
  end do
  seucn = sqrt( sts )
  yeucn = sqrt( yeucn )

  ! Compute a linear relation between gpsupn and cgeps, i.e.,
  ! scalars a and b such that
  !
  !     a * log10( ||g_P(x_ini)|| ) + b = log10(cgeps_ini) and
  !
  !     a * log10( ||g_P(x_fin)|| ) + b = log10(cgeps_fin),
  !
  ! where cgeps_ini and cgeps_fin are provided. Note that if
  ! cgeps_ini is equal to cgeps_fin then cgeps will be always
  ! equal to cgeps_ini and cgeps_fin.

  if ( gpsupn .ne. 0.0d0 ) then
     acgeps = log10( cgepsf / cgepsi ) / log10( cggpnf / gpsupn )
     bcgeps = log10( cgepsi ) - acgeps * log10( gpsupn )
  else
     acgeps = 0.0d0
     bcgeps = cgepsf
  end if
  
  ! Print initial information

  if ( iprintinn .eq. 1 ) then
     if ( mod(iter,10) .eq. 0 ) then
        write(* ,1040)
        write(10,1040)
     end if
     write(* ,1050) iter,f,gpsupn,xsupn
     write(10,1050) iter,f,gpsupn,xsupn

  else if ( iprintinn .ge. 2 ) then
     write(* ,1070) iter,ittext(ittype)
     write(10,1070) iter,ittext(ittype)

     if ( iprintinn .ge. 3 ) then
        write(* ,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,chcnt
        write(10,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,chcnt
     end if

     write(* ,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn, &
                    ssupn,seucn
     write(10,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn, &
                    ssupn,seucn

     if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
        write(* ,1100) min(nprint,nind),nind, &
                       (ind(i),i=1,min(nprint,nind))
        write(* ,1110) nprint,(x(i),i=1,nprint)
        write(* ,1130) nprint,(g(i),i=1,nprint)
        write(* ,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i), &
                               i=1,nprint)
        write(10,1100) min(nprint,nind),nind, &
                       (ind(i),i=1,min(nprint,nind))
        write(10,1110) nprint,(x(i),i=1,nprint)
        write(10,1130) nprint,(g(i),i=1,nprint)
        write(10,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i), &
                               i=1,nprint)
     end if
  end if

  ! Save intermediate data for crash report

  if ( iprintctl(5) ) then
     open(20,file='gencan-tabline.out')
     write(20,1400) f,0.0d0,f,0.0d0,gpsupn,f,0.0d0,f,0.0d0,gpsupn, &
                    0.0d0,0.0d0,inform,9,.false.,n,m,0,iter,fcnt,0, &
                    0,0,0,999.9d0
     close(20)
  end if

  ! ==================================================================
  ! Main loop
  ! ==================================================================

  100 continue

  ! ==================================================================
  ! Test stopping criteria
  ! ==================================================================

  ! Test user-provided stopping criterion

  if ( useustp ) then

     vustop = sstop(n,x,m,lambda,rho,equatn,linear,inform)
     if ( inform .ne. 0 ) go to 500

     if ( vustop ) then
        geninfo = 8

        if ( iprintinn .ge. 1 ) then
           write(*, 1280)
           write(10,1280)
        end if

        go to 500
     end if

  end if

  ! Test whether the continuous-projected-gradient sup-norm
  ! is small enough to declare convergence. In the case of
  ! feasibility problems, check feasibility.

  if ( ( .not. ignoref .and. gpsupn .le. epsopt  ) .or. &
       ( ignoref .and. ( cnormu .le. epsfeas .or. &
       ( cnorm .gt. efstain .and. gpsupn .le. eostain ) ) ) ) then

     geninfo = 0

     if ( iprintinn .ge. 1 ) then
        write(*, 1200)
        write(10,1200)
     end if

     go to 500
  end if

  ! Test whether the number of iterations is exhausted

  if ( iter .ge. maxit ) then
     geninfo = 1

     if ( iprintinn .ge. 1 ) then
        write(*, 1210)
        write(10,1210)
     end if

     go to 500
  end if

  ! ==================================================================
  ! Test stopping criteria related to lack of progress
  ! ==================================================================

  if ( iter .gt. 0 ) then

     ! Test whether we have performed many iterations without
     ! moving from the current point, by checking the functional value

     if ( f .ge. fbest - macheps23 * abs( fbest ) )  then
        itnfp = itnfp + 1
     else
        itnfp = 0
     end if

     ! Test whether we have performed many iterations without
     ! moving from the current point, by checking their gradients

     if ( gpsupn .ge. gpsupnb - macheps23 * gpsupnb .or.  &
          gpeucn .ge. gpeucnb - macheps23 * gpeucnb ) then
        itngp = itngp + 1
     else
        itngp = 0
     end if

     ! Test whether we have performed many iterations without
     ! moving from the current point, by checking the step norm

     if ( ssupn .le. max( macheps23, macheps * xsupn ) .or. &
          seucn .le. max( macheps23, macheps * xeucn ) .or. &
          samep ) then
        itnxp = itnxp + 1
     else
        itnxp = 0
     end if

     itnp = 0
     if ( itnfp .ge. maxinnitnp ) itnp = itnp + 1
     if ( itngp .ge. maxinnitnp ) itnp = itnp + 1
     if ( itnxp .ge. maxinnitnp ) itnp = itnp + 1

     if ( itnp .ge. itnplevel ) then
        geninfo = 3

        if ( iprintinn .ge. 1 ) then
           write(*, 1230)
           write(10,1230)
        end if

        go to 500
     end if

  end if

  ! ==================================================================
  ! Iteration
  ! ==================================================================

  iter = iter + 1

  ! ==================================================================
  ! Compute new iterate
  ! ==================================================================

  ! We abandon the current face if the norm of the internal gradient
  ! (here, internal components of the continuous projected gradient)
  ! is smaller than eta times the norm of the continuous
  ! projected gradient. Using eta = 0.1 is a rather conservative
  ! strategy in the sense that internal iterations are preferred over
  ! SPG iterations. Replace eta = 0.1 by other tolerance in (0,1) if
  ! you find it convenient.

  if ( gieucn .le. eta * gpeucn .or. &
       gisupn .le. eta * gpsupn .or. forceoi ) then

     ! ==============================================================
     ! Some constraints should be abandoned. Compute the new iterate
     ! using an SPG iteration
     ! ==============================================================

     forceoi = .false.
     lvfit = lvfit + 1

     !     Compute spectral steplength

     if ( sty .le. 0.0d0 ) then
        lamspg = max( 1.0d0, xsupn / gpsupn )
     else
        lamspg = sts / sty
     end if
     lamspg = min( lspgma, max( lspgmi, lamspg ) )

     ! Perform safeguarded quadratic interpolation along the
     ! spectral continuous projected gradient

     call spgls(n,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg, &
          xplus,fplus,tmp,d,sevalal,ssetp,lsinfo,inform)
     if ( inform .ne. 0 ) go to 500

     call sevalnal(n,xplus,m,lambda,rho,equatn,linear,gplus,inform)
     if ( inform .ne. 0 ) go to 500

     ! Set iteration type

     ittype = 1

  else

     ! ==============================================================
     ! The new iterate will belong to the closure of the current face
     ! ==============================================================

     innit = innit + 1

     ! Save values of fixed variables for further evaluations

     nfull = n
     xfull(1:n) = x(1:n)

     ! Shrink the point, its gradient and the bounds

     call shrink(nind,x)
     call shrink(nind,g)
     call shrink(nind,l)
     call shrink(nind,u)

     ! Compute the new point using a trust-region approach

     if ( inniter .eq. 'TR' ) then

        call betra(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,trdel, &
             newtrdel,mslamb,epsopt,xeucn,d,xplus,fplus,gplus,trcnt,  &
             chcnt,memfail,trinfo,inform)

        if ( inform .ne. 0 ) go to 500

        ! Set iteration type

        ittype = 4

     end if

     ! Compute the descent direction solving the newtonian system

     if ( inniter .eq. 'NW' ) then

        ! Solve the Newtonian system
        !
        !     ( H + rho A^T A ) x = b
        !
        ! by solving the Martinez-Santos system
        !
        !     H x + A^t y = b
        !     A x - y/rho = 0
        !
        ! with a direct solver.

        call newtd(n,nind,x,l,u,g,m,lambda,rho,equatn,d,adsupn, &
             maxelem,memfail,mindiag,inform)
        if ( inform .ne. 0 ) go to 500

        nwcnt = nwcnt + 1

        if ( .not. memfail ) then

           ! Compute maximum step

           call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)

           ! Perform line search

           call tnls(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,      &
                amax,d,rbdnnz,rbdind,rbdtype,xplus,fplus,gplus,lsinfo, &
                inform)
           if ( inform .ne. 0 ) go to 500

        end if

        ! Set iteration type

        ittype = 3

     end if

     if (   inniter .eq. 'TN' .or.                               &
          ( inniter .eq. 'NW' .and. memfail ) .or.               &
          ( inniter .eq. 'TR' .and. memfail ) .or.               &
          ( inniter .eq. 'TR' .and. trinfo .eq. UNDEFSTEP ) .or. &
          ( inniter .eq. 'TR' .and. trinfo .eq. TOOSMLTRR ) .or. &
          ( inniter .eq. 'TR' .and. trinfo .eq. FSTORDPNT ) .or. &
          ( inniter .eq. 'TR' .and. trinfo .eq. SECORDPNT ) .or. &
          ( inniter .eq. 'TR' .and. trinfo .eq. TCLSBNDRY ) ) then

        ! Compute "trust-region radius" for CG

        if ( iter .eq. 1 ) then
           cgdel = max( 100.0d0, 100.0d0 * xsupn )
        else
           cgdel = max( delmin, 10.0d0 * ssupn )
        end if

        ! Set conjugate gradient stopping criteria.

        cgmaxit = min( nind, 10000 )
        cgeps   = 10.0d0 ** ( acgeps * log10( gpsupn ) + bcgeps )
        cgeps   = max( cgepsf, min( cgepsi, cgeps ) )

        ! Call Conjugate Gradients to solve the Newtonian system

        call cgm(n,nind,x,m,lambda,rho,equatn,linear,l,u,g,cgdel, &
             cgeps,cgmaxit,d,cgiter,cginfo,inform)

        cgcnt = cgcnt + cgiter

        if ( inform .ne. 0 ) go to 500

        ! Compute maximum step

        call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)

        ! Perform line search

        call tnls(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,amax, &
             d,rbdnnz,rbdind,rbdtype,xplus,fplus,gplus,lsinfo,inform)
        if ( inform .ne. 0 ) go to 500

        ! Set iteration type

        ittype = 2

     end if

     ! Set iteration type for the next iteration

     dsupn = 0.0d0
     do i = 1,nind
        dsupn = max( dsupn, abs( d(i) ) )
     end do

     if (   memfail .or.  &
          ( ittype .eq. 3 .and. adsupn .gt. maxelem ) .or.      &
          ( ittype .eq. 3 .and. amax .le. macheps12 * dsupn ) ) then
        inniter = 'TN'
     else
        inniter = innslvr
     end if

     ! Test whether to force a leaving-face iteration is recommended

     if ( nbds .gt. 0 .and.                                   &
          ( ittype .eq. 2 .and. lsinfo .eq. SMALLSTEP ) .or.  &
          ( ittype .eq. 4 .and. trinfo .eq. SMALLSTEP ) ) then
        forceoi = .true.
     end if

     ! Expand the point, its gradient and the bounds. Also expand
     ! xplus and gplus.

     call expand(nind,x)
     call expand(nind,g)
     call expand(nind,l)
     call expand(nind,u)

     xfull(ind(1:nind)) = xplus(1:nind)
     xplus(1:n) = xfull(1:n)

     call expand(nind,gplus)
  end if

  ! For feasibility problems, compute scaled and unscaled feasibilities

  if ( ignoref ) then
     call sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)
     if ( inform .ne. 0 ) go to 500

  else
     cnorm  = 0.0d0
     cnormu = 0.0d0
  end if

  ! Compute s = xplus - x and y = gplus - g

  if ( f      .lt. fbest   ) fbest   = f
  if ( gpsupn .lt. gpsupnb ) gpsupnb = gpsupn
  if ( gpeucn .lt. gpeucnb ) gpeucnb = gpeucn

  sts = 0.0d0
  sty = 0.0d0
  ssupn = 0.0d0
  ysupn = 0.0d0
  yeucn = 0.0d0
  do i = 1,n
     s(i) = xplus(i) - x(i)
     y(i) = gplus(i) - g(i)
     sts  = sts + s(i) ** 2
     sty  = sty + s(i) * y(i)
     ssupn = max( ssupn, abs( s(i) ) )
     ysupn = max( ysupn, abs( y(i) ) )
     yeucn = yeucn + y(i) ** 2
  end do
  seucn = sqrt( sts )
  yeucn = sqrt( yeucn )

  ! Compute samep as it is used within betra and backtrack

  samep = .true.
  do i = 1,n
     if ( xplus(i) .gt. x(i) + macheps23 * max(1.0d0,abs(x(i))) .or. &
          xplus(i) .lt. x(i) - macheps23 * max(1.0d0,abs(x(i))) ) then
        samep = .false.
     end if
  end do

  ! Set new point

  f = fplus

  do i = 1,n
     x(i) = xplus(i)
     g(i) = gplus(i)
  end do

  ! ==================================================================
  ! Prepare for the next iteration
  ! ==================================================================

  ! Compute continuous-project-gradient Euclidian and Sup norms,
  ! internal gradient Euclidian norm, and store in nind the number of
  ! free variables and in array ind their identifiers. Verify whether
  ! x and xprev belong to the same face.

  sameface = .true.
  nindprev = nind

  nind   = 0
  xsupn  = 0.0d0
  xeucn  = 0.0d0
  gpsupn = 0.0d0
  gpeucn = 0.0d0
  gisupn = 0.0d0
  gieucn = 0.0d0
  do i = 1,n
     if ( abs( x(i) ) .ge. xsupn ) then
        xsupn = abs( x(i) )
        xind = i
     end if
     xeucn = xeucn + x(i) ** 2
     gpi = x(i) - g(i)
     if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
        gpi = - g(i)
     else
        gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
     end if
     gpsupn = max( gpsupn, abs( gpi ) )
     gpeucn = gpeucn + gpi ** 2
     if ( x(i) .gt. l(i) + macheps23 * max(1.0d0,abs(l(i))) .and. &
          x(i) .lt. u(i) - macheps23 * max(1.0d0,abs(u(i))) ) then
        gisupn    = max( gisupn, abs( gpi ) )
        gieucn    = gieucn + gpi ** 2
        nind      = nind + 1
        if ( nind .gt. nindprev .or. ind(nind) .ne. i ) then
           sameface = .false.
        end if
        ind(nind) = i
     end if
  end do
  xeucn  = sqrt( xeucn  )
  gpeucn = sqrt( gpeucn )
  gieucn = sqrt( gieucn )

  ! Print information of this iteration

  if ( iprintinn .eq. 1 ) then
     if ( mod(iter,10) .eq. 0 ) then
        write(* ,1040)
        write(10,1040)
     end if
     write(* ,1060) iter,f,gpsupn,xsupn,ssupn,nind,lvfit,innit,fcnt
     write(10,1060) iter,f,gpsupn,xsupn,ssupn,nind,lvfit,innit,fcnt

  else if ( iprintinn .ge. 2 ) then
     write(* ,1070) iter,ittext(ittype)
     write(10,1070) iter,ittext(ittype)

     if ( iprintinn .ge. 3 ) then
        write(* ,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,chcnt
        write(10,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,chcnt
     end if

     write(* ,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn, &
                    ssupn,seucn
     write(10,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn, &
                    ssupn,seucn

     if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
        write(* ,1100) min(nprint,nind),nind, &
                       (ind(i),i=1,min(nprint,nind))
        write(* ,1110) nprint,(x(i),i=1,nprint)
        write(* ,1130) nprint,(g(i),i=1,nprint)
        write(* ,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i), &
                               i=1,nprint)
        write(10,1100) min(nprint,nind),nind, &
                      (ind(i),i=1,min(nprint,nind))
        write(10,1110) nprint,(x(i),i=1,nprint)
        write(10,1130) nprint,(g(i),i=1,nprint)
        write(10,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i), &
                               i=1,nprint)
     end if
  end if

  ! Save intermediate data for crash report

  if ( iprintctl(5) ) then
     open(20,file='gencan-tabline.out')
     write(20,1400) f,0.0d0,f,0.0d0,gpsupn,f,0.0d0,f,0.0d0,gpsupn, &
          0.0d0,0.0d0,inform,9,.false.,n,m,0,iter,fcnt,0,          &
          0,0,0,999.9d0
     close(20)
  end if

  ! ==================================================================
  ! Test line-search and trust-region related stopping criteria
  ! ==================================================================

  ! Test whether the functional value is unbounded

  if ( ( ittype .le. 3 .and. lsinfo .eq. UNBOUNDED ) .or. &
       ( ittype .eq. 4 .and. trinfo .eq. UNBOUNDED ) ) then
     geninfo = 6

     if ( iprintinn .ge. 1 ) then
        write(*, 1260)
        write(10,1260)
     end if

     go to 500
  end if

  ! ==================================================================
  ! Iterate
  ! ==================================================================

  go to 100

  ! ==================================================================
  ! End of main loop
  ! ==================================================================

  500 continue

  ! ==================================================================
  ! Deallocate global arrays in modules
  ! ==================================================================

  call deallocategencanarrays(allocerr)

  if ( allocerr .ne. 0 ) then
     inform = - 94
     subname = 'GENCAN'
     call reperr(inform,subname)
     return
  end if

  ! ==================================================================
  ! NON-EXECUTABLE STATEMENTS
  ! ==================================================================

 1000 format(/,5X,'Entry to GENCAN.',         &
             /,5X,'Number of variables: ',I7, &
             /,5X,'Required optimality tolerance: ',1P,D7.1)

 1010 format(/,5X,'Lower bounds  (first ',I7,' components):', &
             /,(5X, 6(1X,1P,D11.4)))
 1020 format(/,5X,'Upper bounds  (first ',I7,' components):', &
             /,(5X, 6(1X,1P,D11.4)))
 1030 format(/,5X,'Initial point (first ',I7,' components):', &
             /,(5X, 6(1X,1P,D11.4)))

 1040 format(/,5X,'It',9X,'function',2X,'gpsupn',3X,'xsupn',1X, &
               'stesupn',3X,'nfree',3X,'SPGit',1X,'in-face',4X,'fcnt')
 1050 format(     I7,1X,1P,D16.8,2(1X,1P,D7.1))
 1060 format(     I7,1X,1P,D16.8,3(1X,1P,D7.1),4(1X,I7))

 1070 format(/,5X,'GENCAN ITERATION                  ',12X,'= ',I7,1X, &
               A21)
 1080 format(/,5X,'Leaving-face iterations           ',12X,'= ',I7, &
             /,5X,'Inner-to-the-face iterations      ',12X,'= ',I7, &
             /,5X,'Functional evaluations            ',12X,'= ',I7, &
             /,5X,'Conjugate gradient iterations     ',12X,'= ',I7, &
             /,5X,'Newtonian system factorizations   ',12X,'= ',I7, &
             /,5X,'Trust-region iterations           ',12X,'= ',I7, &
             /,5X,'Trust-region matrix factorizations',12X,'= ',I7)

 1090 format(/,5X,'Functional value                              = ', &
                   1P,D24.16,                                         &
             /,5X,'Sup-norm of the continuous projected gradient = ', &
                   1P,D7.1,' 2-norm = ',1P,D7.1,                      &
             /,5X,'Sup-norm of the internal projection of gp     = ', &
                   1P,D7.1,' 2-norm = ',1P,D7.1,                      &
             /,5X,'Sup-norm of x (attained at x_{', I7,'})       = ', &
                   1P,D7.1,' 2-norm = ',1P,D7.1,                      &
             /,5X,'Sup-norm of x - x_{prev}                      = ', &
                   1P,D7.1,' 2-norm = ',1P,D7.1)
 1100 format(/,5X,'Current free variables (first ',I7,', total ', &
                  'number ',I7,'): ',                             &
             /,(5X, 6(1X,I7)))
 1110 format(/,5X,'Current point (first ',I7, ' components): ', &
             /,(5X, 6(1X,1P,D11.4)))
 1130 format(/,5X,'Current gradient (first ',I7,' components): ', &
             /,(5X, 6(1X,1P,D11.4)))
 1140 format(/,5X,'Current continuous projected gradient (first ',I7, &
               5X,'components): ',                                    &
             /,(5X, 6(1X,1P,D11.4)))

 1200 format(/,5X,'Flag of GENCAN: Solution was found.',/)
 1210 format(/,5X,'Flag of GENCAN: Maximum of iterations reached.',/)
 1230 format(/,5X,'Flag of GENCAN: Lack of progress in the functional ', &
                  'value, its gradient and',/,5X,'the current point. ',  &
                  'Probably, an exaggerated small norm of the ',         &
                  'continuous',/,5X,'projected gradient is being ',      &
                  'required to declare convergence.',/)
 1260 format(/,5X,'Flag of GENCAN: Objective function seems to be ', &
                  'unbounded.',/)
 1280 format(/,5X,'Flag of GENCAN: User-provided stopping criterion ', &
                  'satisfied.',/)

 1400 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D8.1, &
             1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D8.1, &
             1X,1P,D7.1,1X,1P,D7.1,1X,I3,1X,I1,1X,L1,1X,I6,1X,I6,1X,I3,  &
             1X,I7,1X,I7,1X,I2,1X,I7,1X,I7,1X,I7,0P,F8.2)

end subroutine gencan

! ******************************************************************
! ******************************************************************

subroutine compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)
  
  use modmachconst

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)  :: nind
  integer,      intent(out) :: rbdnnz
  real(kind=8), intent(out) :: amax

  ! ARRAY ARGUMENTS
  integer,      intent(out) :: rbdind(nind)
  character,    intent(out) :: rbdtype(nind)
  real(kind=8), intent(in)  :: d(nind),l(nind),u(nind),x(nind)

  ! Compute maximum step amax > 0 such that l <= x + amax d <= u.

  ! LOCAL SCALARS
  integer      :: i
  real(kind=8) :: amaxi

  rbdnnz = 0

  amax = bignum

  do i = 1,nind
     if ( d(i) .gt. 0.0d0 ) then

        amaxi = ( u(i) - x(i) ) / d(i)

        if ( amaxi .lt. amax ) then
           amax       = amaxi
           rbdnnz     = 1
           rbdind(1)  = i
           rbdtype(1) = 'U'
        else if ( amaxi .eq. amax ) then
           rbdnnz          = rbdnnz + 1
           rbdind(rbdnnz)  = i
           rbdtype(rbdnnz) = 'U'
        end if

     else if ( d(i) .lt. 0.0d0 ) then

        amaxi = ( l(i) - x(i) ) / d(i)

        if ( amaxi .lt. amax ) then
           amax      = amaxi
           rbdnnz     = 1
           rbdind(1)  = i
           rbdtype(1) = 'L'
        else if ( amaxi .eq. amax ) then
           rbdnnz          = rbdnnz + 1
           rbdind(rbdnnz)  = i
           rbdtype(rbdnnz) = 'L'
        end if

     end if
  end do

end subroutine compamax
