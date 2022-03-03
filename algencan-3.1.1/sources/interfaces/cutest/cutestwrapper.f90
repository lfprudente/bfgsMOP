! *****************************************************************
! *****************************************************************

module cutestwrapper

  integer :: mc,nc

  integer, allocatable :: ccor(:),cmap(:),slaind(:)
  real(kind=8), allocatable :: ca(:),cb(:)

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalf(n,x,f,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    flag = - 1

  end subroutine evalf

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalg(n,x,g,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)

    flag = - 1

  end subroutine evalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalh(n,x,hlin,hcol,hval,hnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: lim,n
    integer, intent(out) :: flag,hnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: hcol(lim),hlin(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    flag = - 1

  end subroutine evalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalc(n,x,ind,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ind,n
    integer, intent(out) :: flag
    real(kind=8), intent(out) :: c

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    flag = - 1

  end subroutine evalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: ind,lim,n
    integer, intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    flag = - 1

  end subroutine evaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalhc(n,x,hclin,hccol,hcval,hcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: lim,n
    integer, intent(out) :: flag,hcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: hccol(lim),hclin(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: hcval(lim)

    flag = - 1

  end subroutine evalhc

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalfc(n,x,f,m,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n
    integer, intent(out) :: flag
    real(kind=8), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    integer :: i,myid,sind,status

    ! EXTERNAL FUNCTIONS
    integer, external :: omp_get_thread_num

    flag = 0

    myid = omp_get_thread_num()

    call CUTEST_cfn_threaded(status,nc,mc,x,f,c,myid+1)

    if ( status .ne. 0 ) then
       flag = - 1
       return
    end if

    do i = m,1,-1
       c(i) = ca(i) * c(cmap(i)) + cb(i)
    end do

    do i = 1,m
       sind = slaind(i)
       if ( sind .ne. - 1 ) then
          c(i) = c(i) - x(sind)
       end if
    end do

  end subroutine evalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: lim,m,n
    integer, intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer, intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    ! LOCAL SCALARS
    integer :: fu2,fun,i,jcnnztmp,myid,sind,status,var

    ! LOCAL ARRAYS
    real(kind=8) :: dum3(1)

    ! EXTERNAL FUNCTIONS
    integer, external :: omp_get_thread_num

    flag = 0

    lmem = .false.

    myid = omp_get_thread_num()

    call CUTEST_csgr_threaded(status,nc,mc,x,dum3,.false.,jcnnz,lim,jcval,jcvar,jcfun,myid+1)

    if ( status .ne. 0 ) then
       flag = - 1
       return
    end if

    ! Remove gradient from the sparse structure

    g(1:n) = 0.0d0

    i = 1

    do while ( i .le. jcnnz )
       fun = jcfun(i)
       var = jcvar(i)

       if ( fun .eq. 0 ) then
          g(var) = g(var) + jcval(i)

          if ( i .ne. jcnnz ) then
             call intswap(jcfun(jcnnz),jcfun(i))
             call intswap(jcvar(jcnnz),jcvar(i))
             call dblswap(jcval(jcnnz),jcval(i))
          end if

          jcnnz = jcnnz - 1

       else
          i = i + 1
       end if
    end do

    ! Duplicate entrances corresponding to splitted constraints

    jcnnztmp = jcnnz

    do i = 1,jcnnztmp
       fun = jcfun(i)
       fu2 = ccor(fun)

       if ( fu2 .ne. 0 ) then
          jcnnz = jcnnz + 1
          jcfun(jcnnz) = fu2
          jcvar(jcnnz) = jcvar(i)
          jcval(jcnnz) = ca(fu2) * jcval(i)
       end if

       jcval(i) = ca(fun) * jcval(i)
    end do

    ! Add effect of slack variables

    do i = 1,m
       sind = slaind(i)
       if ( sind .ne. - 1 ) then
          jcnnz = jcnnz + 1
          jcfun(jcnnz) = i
          jcvar(jcnnz) = sind
          jcval(jcnnz) = - 1.0d0
       end if
    end do

  end subroutine evalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalgjacp(n,x,g,m,p,q,work,gotj,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(inout) :: gotj
    integer, intent(in) :: m,n
    integer, intent(out) :: flag
    character, intent(in) ::  work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    real(kind=8), intent(inout) :: p(m),q(n)

    ! The usage of this subroutine for solving problems in the CUTEst
    ! collection (in replacement of evalgjac)  was never evaluated.

    flag = - 1

  end subroutine evalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(out) :: lmem
    integer, intent(in) :: lim,m,n
    integer, intent(out) :: flag,hnnz
    real(kind=8), intent(in) :: sf

    ! ARRAY ARGUMENTS
    integer, intent(out) :: hlin(lim),hcol(lim)
    real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    integer :: i,myid,status

    ! LOCAL ARRAYS
    real(kind=8) :: v(m)

    ! EXTERNAL FUNCTIONS
    integer, external :: omp_get_thread_num

    flag = 0

    lmem = .false.

    myid = omp_get_thread_num()

    if ( sf .eq. 0.0d0 ) then
       v(1:mc) = 0.0d0

       do i = 1,m
          v(cmap(i)) = v(cmap(i)) + ca(i) * lambda(i) * sc(i)
       end do

       call CUTEST_cshc_threaded(status,nc,mc,x,v,hnnz,lim,hval,hlin,hcol,myid+1)

       ! Interchange row and column indices

       do i = 1,hnnz
          call intswap(hlin(i),hcol(i))
       end do

       return
    end if

    v(1:mc) = 0.0d0

    do i = 1,m
       v(cmap(i)) = v(cmap(i)) + ca(i) * lambda(i) * sc(i) / sf
    end do

    call CUTEST_csh_threaded(status,nc,mc,x,v,hnnz,lim,hval,hlin,hcol,myid+1)

    if ( status .ne. 0 ) then
       flag = - 1
       return
    end if

    ! Interchange row and column indices

    do i = 1,hnnz
       call intswap(hlin(i),hcol(i))
       hval(i) = hval(i) * sf
    end do

  end subroutine evalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(inout) :: goth
    integer, intent(in) :: m,n
    integer, intent(out) :: flag
    real(kind=8), intent(in) :: sf

    ! ARRAY ARGUMENTS
    real(kind=8) :: hp(n),lambda(m),p(n),sc(m),x(n)

    flag = - 1

  end subroutine evalhlp

  ! ******************************************************************
  ! ******************************************************************

  subroutine intswap(i,j)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(inout) :: i,j

    ! LOCAL SCALARS
    integer :: tmp

    tmp = i
    i   = j
    j   = tmp

  end subroutine intswap

  ! ******************************************************************
  ! ******************************************************************

  subroutine dblswap(a,b)

    implicit none

    ! SCALAR ARGUMENTS
    real(kind=8), intent(inout) :: a,b

    ! LOCAL SCALARS
    real(kind=8) :: tmp

    tmp = a
    a   = b
    b   = tmp

  end subroutine dblswap


end module cutestwrapper

! ******************************************************************
! ******************************************************************

program algencanma

  use cutestwrapper

  implicit none

  ! LOCAL SCALARS
  logical :: checkder,useslacks
  integer :: allocstat,cind,e_order,hnnzmax,inform,j,jcnnzmax,k, &
       l_order,m,n,nranges,nvparam,sind,status,v_order
  real(kind=8) :: cnorm,dum,efacc,efstain,eoacc,epsfeas,epsopt, &
       eostain,f,nlpsupn,snorm

  ! LOCAL ARRAYS
  logical :: coded(11)
  logical, allocatable :: equatn(:),equatnc(:),linear(:),linearc(:)
  real(kind=8), allocatable :: c(:),cl(:),cu(:),l(:),lc(:),lambda(:), &
       lambdac(:),u(:),uc(:),x(:),xc(:)
  character(len=10) :: pname
  character(len=80) :: specfnm,outputfnm,vparam(10)
  character(len=10), allocatable :: cnames(:),xnames(:)

  ! Set problem data

  open(10,file='OUTSDIF.d',form='formatted',status='old')

  rewind 10

  call CUTEST_cdimen(status,10,nc,mc)

  if ( status .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, ', &
          'status from subroutine CUTEST_cdimen is not equal to zero: ',status
     stop
  end if

  allocate(xc(nc),lc(nc),uc(nc),lambdac(mc),cl(mc),cu(mc),equatnc(mc), &
       linearc(mc),stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, allocation error.'
     stop
  end if

  e_order = 1
  l_order = 0
  v_order = 0

  call CUTEST_csetup_threaded(status,10,20,4,11,nc,mc,xc,lc,uc,lambdac,cl,cu,equatnc, &
       linearc,e_order,l_order,v_order)

  if ( status .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, ', &
          'status from subroutine CUTEST_csetup is not equal to zero: ',status
     stop
  end if

  close(10)

  allocate(xnames(nc),cnames(mc),stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, allocation error.'
     stop
  end if

  call CUTEST_cnames(status,nc,mc,pname,xnames,cnames)

  if ( status .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, ', &
          'status from subroutine CUTEST_cnames is not equal to zero: ',status
     stop
  end if

  ! There are two options for two-sides inequality constraints:
  !
  ! (a) Add a slack variable and transform it into an equality
  !     constraints plus a bound constraint for the slack.
  !
  ! (b) Split it into two one-side inequality constraints.
  !
  ! There is a local variable of this subroutine called useslacks.
  ! If useslacks is TRUE then we implement (a) else, we implement (b).

  useslacks = .false.

  nranges = count(.not. equatnc(1:mc) .and. cl(1:mc) .gt. - 1.0d+20 .and. cu(1:mc) .lt. 1.0d+20)

  if ( useslacks ) then
     n = nc + nranges
     m = mc
  else
     n = nc
     m = mc + nranges
  end if

  allocate(x(n),l(n),u(n),lambda(m),equatn(m),linear(m),ca(m),cb(m), &
       ccor(mc),cmap(m),slaind(m),stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, allocation error.'
     stop
  end if

  x(1:nc) = xc(1:nc)
  l(1:nc) = lc(1:nc)
  u(1:nc) = uc(1:nc)

  lambda(1:mc) = lambdac(1:mc)
  equatn(1:mc) = equatnc(1:mc)
  linear(1:mc) = linearc(1:mc)

  deallocate(xc,lc,uc,lambdac,equatnc,linearc,stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, allocation error.'
     stop
  end if

  ! Compute constraints to initialize slacks

  if ( useslacks ) then
     allocate(c(mc),stat=allocstat)
     if ( allocstat .ne. 0 ) then
        write(*,*) 'In the CUTEST interface main file, allocation error.'
        stop
     end if

     call CUTEST_cfn_threaded(status,nc,mc,x,dum,c,1)
     if ( status .ne. 0 ) then
        write(*,*) 'In the CUTEST interface main file, ', &
             'status from subroutine CUTEST_cfn is not equal to zero: ',status
        stop
     end if
  end if

  k = 0

  do j = 1,mc

     if ( equatn(j) ) then

        ! Equality constraint.

        slaind(j) =     - 1
        ccor(j)   =       0
        cmap(j)   =       j
        ca(j)     =   1.0d0
        cb(j)     = - cu(j)

     else if ( cl(j) .gt. - 1.0d+20 .and. cu(j) .lt. 1.0d+20 ) then

        ! Ranged inequality constraint: add slack or split it.

        k = k + 1

        if ( useslacks ) then

           ! Replace by c(x) - s = 0 and cl <= s <= cu.

           sind = nc + k

           l(sind) = cl(j)
           u(sind) = cu(j)
           x(sind) = max( cl(j), min( c(j), cu(j) ) )

           slaind(j) =  sind
           ccor(j)   =     0
           cmap(j)   =     j
           ca(j)     = 1.0d0
           cb(j)     = 0.0d0

           equatn(j) = .true.

        else

           ! Split into c(x) - cu <= 0 and cl - c(x) <= 0.

           cind = mc + k

           equatn(cind) = equatn(j)
           linear(cind) = linear(j)

           slaind(cind) =     - 1
           cmap(cind)   =       j
           ca(cind)     =   1.0d0
           cb(cind)     = - cu(j)

           slaind(j)    =     - 1
           ccor(j)      =    cind
           cmap(j)      =       j
           ca(j)        = - 1.0d0
           cb(j)        =   cl(j)
        end if

     else if ( cu(j) .lt. 1.0d+20 ) then

        ! Inequality constraint of type c(x) <= cu.

        slaind(j) =     - 1
        ccor(j)   =       0
        cmap(j)   =       j
        ca(j)     =   1.0d0
        cb(j)     = - cu(j)

     else if ( cl(j) .gt. - 1.0d+20 ) then

        ! Inequality constraint of type cl <= c(x).

        slaind(j) =     - 1
        ccor(j)   =       0
        cmap(j)   =       j
        ca(j)     = - 1.0d0
        cb(j)     =   cl(j)

     end if

  end do

  deallocate(cl,cu,stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, allocation error.'
     stop
  end if

  if ( useslacks ) then
     deallocate(c,stat=allocstat)
     if ( allocstat .ne. 0 ) then
        write(*,*) 'In the CUTEst interface main file, allocation error.'
        stop
     end if
  end if

  ! Lagrange multipliers approximation

  lambda(1:m) = 0.0d0

  ! In this CUTEst interface, subroutines fcsub, gjacsub, and hlsub
  ! are present. Subroutines fsub, gsub, hsub, csub, jacsub,
  ! gjacpsub, hlpsub, and hcsub are not.

  coded( 1) = .false. ! fsub
  coded( 2) = .false. ! gsub
  coded( 3) = .false. ! hsub
  coded( 4) = .false. ! csub
  coded( 5) = .false. ! jacsub
  coded( 6) = .false. ! hcsub
  coded( 7) = .true.  ! fcsub
  coded( 8) = .true.  ! gjacsub
  coded( 9) = .false. ! gjacpsub
  coded(10) = .true.  ! hlsub
  coded(11) = .false. ! hlpsub

  checkder = .false.

  call CUTEST_cdimsj(status,jcnnzmax)

  if ( status .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, ', &
          'status from subroutine CUTEST_cdimsj is not equal to zero: ',status
     stop
  end if

  if ( useslacks ) then
     jcnnzmax = jcnnzmax + nranges
  else
     jcnnzmax = 2 * jcnnzmax
  end if

  call CUTEST_cdimsh(status,hnnzmax)

  if ( status .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, ', &
          'status from subroutine CUTEST_cdimsh is not equal to zero: ',status
     stop
  end if

  ! hnnzmax must be an upper bound on the number of elements of the
  ! Hessian of the Lagrangian plus the new elements that appear when
  ! adding (to the Hessian of the Lagrangian) matrix rho \sum_j \nabla
  ! cj(x) \nabla cj(x)^t to obtain the Hessian of the Augmented
  ! Lagrangian. But this additional space is only need when Algencan
  ! uses an Euclidian trust-region approach that is recommended only
  ! for problems with no more than 500 variables. Therefore, since we
  ! are not able to estimate this quantity here (when using CUTEst),
  ! we arbitrarily add to hnnzmax the quantity 1,000,000 >= 500^2.

  hnnzmax = hnnzmax + 1000000

  epsfeas   = 1.0d-08
  epsopt    = 1.0d-08

  efstain   = sqrt( epsfeas )
  eostain   = epsopt ** 1.5d0

  efacc     = sqrt( epsfeas )
  eoacc     = sqrt( epsopt )

  outputfnm = ''
  specfnm   = 'algencan.dat'

  nvparam = 0

  write(*,*) 'Solving problem: ',pname

  call algencan(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc, &
       evalgjac,evalgjacp,evalhl,evalhlp,jcnnzmax,hnnzmax,epsfeas,epsopt, &
       efstain,eostain,efacc,eoacc,outputfnm,specfnm,nvparam,vparam,n,x, &
       l,u,m,lambda,equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn, &
       inform)

  deallocate(x,l,u,lambda,equatn,linear,ca,cb,ccor,cmap,slaind,stat=allocstat)
  if ( allocstat .ne. 0 ) then
     write(*,*) 'In the CUTEST interface main file, allocation error.'
     stop
  end if

  stop

end program algencanma

