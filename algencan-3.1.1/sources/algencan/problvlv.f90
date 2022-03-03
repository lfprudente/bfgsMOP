! CALL USER-PROVIDED SUBROUTINES

module problvlv

  implicit none

  ! SCALARS
  logical, protected :: safemode
  integer, protected :: efcnt,efccnt,egcnt,egjccnt,egjcpcnt, &
       ehcnt,ehlcnt,ehlpcnt,fcnt,nsf,nsfmax

  ! ARRAYS
  integer, allocatable, protected :: eccnt(:),ehccnt(:),ejccnt(:)

  ! SUBROUTINES
  public :: vinip,vendp,vevalf,vevalg,vevalh,vevalc,vevaljac, &
       vevalhc,vevalhl,vevalhlp,vevalfc,vevalgjac,vevalgjacp

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine setproblvlv(val_safemode,val_nsfmax)

    implicit none

    ! SCALARS
    logical, intent(in), optional :: val_safemode
    integer, intent(in), optional :: val_nsfmax

    safemode = .false.
    if ( present( val_safemode ) ) safemode = val_safemode

    nsfmax = 0
    if ( present( val_nsfmax ) ) nsfmax = val_nsfmax

    nsf = 0

  end subroutine setproblvlv

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine vinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)

    use modouttyp, only: iprintctl
    use modprobdata, only: meq,mineq,nbds,nfix

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in)    :: checkder
    integer, intent(in)    :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical,      intent(in)    :: coded(11),equatn(m),linear(m)
    real(kind=8), intent(in)    :: lambda(m)
    real(kind=8), intent(inout) :: l(n),u(n),x(n)

    ! LOCAL SCALARS
    integer :: allocerr

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! Set counters

    allocate(eccnt(m),ehccnt(m),ejccnt(m),stat=allocerr)

    if ( allocerr .ne. 0 ) then
       inform = - 93
       subname = 'VINIP'
       call reperr(inform,subname)
       return
    end if

    fcnt        = 0
    efcnt       = 0
    efccnt      = 0
    egcnt       = 0
    egjccnt     = 0
    egjcpcnt    = 0
    ehcnt       = 0
    ehlcnt      = 0
    ehlpcnt     = 0
    eccnt(1:m)  = 0
    ejccnt(1:m) = 0
    ehccnt(1:m) = 0
    
    ! Check bounds feasibility

    if ( minval( u(1:n) - l(1:n) ) .lt. 0.0d0 ) then
       inform = - 95
       subname = 'VINIP'
       call reperr(inform,subname)
       return
    end if

    ! Avoid huge bounds

    l(1:n) = max( l(1:n), - 1.0d+20 )
    u(1:n) = min( u(1:n),   1.0d+20 )

    ! Project initial guess

    x(1:n) = max( l(1:n), min( x(1:n), u(1:n) ) )

    ! Count bound constraints, equality constraints, and inequality constraints

    meq = count( equatn(1:m) )
    mineq = m - meq
    nbds = count( l(1:n) .gt. - 1.0d+20 ) + count( u(1:n) .lt. 1.0d+20 )
    nfix = count( l(1:n) .eq. u(1:n) )

    if ( iprintctl(2) ) then
       write(* ,100) n,meq,mineq,nbds,nfix
       write(10,100) n,meq,mineq,nbds,nfix
    end if

    ! Write classification line of original model

    if ( iprintctl(5) ) then
       open(50,file='class-tabline.out')
       write(50,400) n,meq,mineq,nbds,nfix
       close(50)
    end if

    ! Check derivatives

    if ( checkder ) then
       call checkd(n,x,l,u,m,lambda,inform)
       if ( inform .ne. 0 ) return
    end if

    ! Non-executable statements

100 format(/,1X,'Number of variables               : ',I7, &
           /,1X,'Number of equality constraints    : ',I7, &
           /,1X,'Number of inequality constraints  : ',I7, &
           /,1X,'Number of bound constraints       : ',I7, &
           /,1X,'Number of fixed variables         : ',I7)

400 format(  1X,I6,1X,I6,1X,I6,1X,I6,1X,I6)

  end subroutine vinip

  ! ******************************************************************
  ! ******************************************************************

  subroutine vendp(n,x,l,u,m,lambda,equatn,linear,inform)

    use modouttyp, only: iprintctl,solfnm
    use probgiven, only: fcoded,gcoded,hcoded,ccoded,jaccoded, &
         hccoded,fccoded,gjaccoded,gjacpcoded,hlcoded,hlpcoded

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical,      intent(in) :: equatn(m),linear(m)
    real(kind=8), intent(in) :: l(n),lambda(m),u(n),x(n)

    ! LOCAL SCALARS
    integer :: allocerr,aveeccnt,aveehccnt,aveejccnt,i,j,toteccnt, &
         totejccnt,totehccnt
         
    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! Report counters

    if ( iprintctl(4) ) then
       toteccnt  = 0
       totejccnt = 0
       totehccnt = 0

       do j = 1,m
          toteccnt  = toteccnt  + eccnt(j)
          totejccnt = totejccnt + ejccnt(j)
          totehccnt = totehccnt + ehccnt(j)
       end do

       if ( m .gt. 0 ) then
          aveeccnt  = toteccnt  / m
          aveejccnt = totejccnt / m
          aveehccnt = totehccnt / m
       else
          aveeccnt  = 0
          aveejccnt = 0
          aveehccnt = 0
       end if

       write(* ,100) fcoded,efcnt,gcoded,egcnt,hcoded,ehcnt,ccoded,  &
                     toteccnt,aveeccnt,jaccoded,totejccnt,aveejccnt, &
                     hccoded,totehccnt,aveehccnt,fccoded,efccnt,     &
                     gjaccoded,egjccnt,gjacpcoded,egjcpcnt,hlcoded,  &
                     ehlcnt,hlpcoded,ehlpcnt
       write(10,100) fcoded,efcnt,gcoded,egcnt,hcoded,ehcnt,ccoded,  &
                     toteccnt,aveeccnt,jaccoded,totejccnt,aveejccnt, &
                     hccoded,totehccnt,aveehccnt,fccoded,efccnt,     &
                     gjaccoded,egjccnt,gjacpcoded,egjcpcnt,hlcoded,  &
                     ehlcnt,hlpcoded,ehlpcnt
    end if

    deallocate(eccnt,ehccnt,ejccnt,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 94
       subname = 'VENDP'
       call reperr(inform,subname)
       return
    end if

  ! Save solution

    if ( solfnm .ne. '' ) then

       ! Save solution

       open(20,file=solfnm)

       ! Point

       write(20,200)
       do i = 1,n
          write(20,300) i,x(i)
       end do

       ! Lagrange multipliers

       if ( m .gt. 0 ) then
          write(20,400)
          do i = 1,m
             write(20,300) i,lambda(i)
          end do
       end if

       close(20)

    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'User-provided subroutines calls counters: ', &
          //,1X,'Subroutine fsub     (coded=',L1,'): ',I8,    &
           /,1X,'Subroutine gsub     (coded=',L1,'): ',I8,    &
           /,1X,'Subroutine hsub     (coded=',L1,'): ',I8,    &
           /,1X,'Subroutine csub     (coded=',L1,'): ',I8,    &
             1X,'(',I8,' calls per constraint in avg)',       &
           /,1X,'Subroutine jacsub   (coded=',L1,'): ',I8,    &
             1X,'(',I8,' calls per constraint in avg)',       &
           /,1X,'Subroutine hcsub    (coded=',L1,'): ',I8,    &
             1X,'(',I8,' calls per constraint in avg)',       &
           /,1X,'Subroutine fcsub    (coded=',L1,'): ',I8,    &
           /,1X,'Subroutine gjacsub  (coded=',L1,'): ',I8,    &
           /,1X,'Subroutine gjacpsub (coded=',L1,'): ',I8,    &
           /,1X,'Subroutine hlsub    (coded=',L1,'): ',I8,    &
           /,1X,'Subroutine hlpsub   (coded=',L1,'): ',I8,/)

200 format(/,'FINAL POINT:',//,2X,'INDEX',16X,'X(INDEX)')

300 format(I7,1P,D24.16)

400 format(/,'FINAL ESTIMATION OF THE LAGRANGE MULTIPLIERS: ', &
           //,2X,'INDEX',11X,'LAMBDA(INDEX)')

  end subroutine vendp

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalf(n,x,f,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: flag

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    call evalf(n,x,f,flag)

    efcnt = efcnt + 1
    fcnt  = fcnt  + 1

    if ( flag .ne. 0 ) then
       if ( iprintctl(3) ) then
          write(* ,100)
          write(10,100)
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 80
          subname = 'VEVALF'
          call reperr(inform,subname)
          return
       end if
    end if

    if ( flag .eq. 0 .and. .not. IsANumber(f) ) then
       if ( iprintctl(3) ) then
          write(* ,200)
          write(* ,300) f
          write(10,200)
          write(10,300) f
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 80
          subname = 'VEVALF'
          call reperr(inform,subname)
          return
       end if
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALF WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALF WARNING: The objective function value ',   &
         'computed by the user-supplied subroutine EVALF may be ', &
         '+Inf, -Inf or NaN.')

300 format(/,1X,'Value: ',1P,D24.16)

  end subroutine vevalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalg(n,x,g,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    integer :: flag,i

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    if ( gcoded ) then
       call evalg(n,x,g,flag)

       egcnt = egcnt + 1

       if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
             write(* ,100)
             write(10,100)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 81
             subname = 'VEVALG'
             call reperr(inform,subname)
             return
          end if
       end if

       do i = 1,n
          if ( flag .eq. 0 .and. .not. IsANumber(g(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,200)
                write(* ,300) n,i,g(i)
                write(10,200)
                write(10,300) n,i,g(i)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 81
                subname = 'VEVALG'
                call reperr(inform,subname)
                return
             end if
          end if
       end do

    else
       call ivevalg(n,x,g,inform)
       if ( inform .ne. 0 ) return
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALG WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALG WARNING: There is an element whose value may ',   &
                'be +Inf, -Inf or NaN in the gradient of the objective ', &
                'function computed by the user-supplied subroutine ',     &
                'EVALG.')

300 format(/,1X,'Dimension of the space: ',I16, &
           /,1X,'Position              : ',I16, &
           /,1X,'Value                 : ',1P,D24.16)

  end subroutine vevalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine ivevalg(n,x,g,inform)

    use modmachconst, only: macheps13

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    integer      :: i
    real(kind=8) :: fminus,fplus,step

    ! LOCAL ARRAYS
    real(kind=8) :: xtmp(n)

    xtmp(1:n) = x(1:n)

    do i = 1,n
       step = macheps13 * max( 1.0d0, abs( x(i) ) )

       xtmp(i) = x(i) + step
       call vsetp(n,xtmp,inform)
       if ( inform .ne. 0 ) return
       call vevalf(n,xtmp,fplus,inform)
       if ( inform .ne. 0 ) return

       xtmp(i) = x(i) - step
       call vsetp(n,xtmp,inform)
       if ( inform .ne. 0 ) return
       call vevalf(n,xtmp,fminus,inform)
       if ( inform .ne. 0 ) return

       g(i) = ( fplus - fminus ) / ( 2.0d0 * step )

       xtmp(i) = x(i)
    end do

  end subroutine ivevalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalh(n,x,hrow,hcol,hval,hnnz,lim,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: lim,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hcol(lim),hrow(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    logical :: lmem
    integer :: flag,i

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    call evalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

    ehcnt = ehcnt + 1

    if ( flag .ne. 0 ) then
       if ( iprintctl(3) ) then
          write(* ,100)
          write(10,100)
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 82
          subname = 'VEVALH'
          call reperr(inform,subname)
          return
       end if
    end if

    if ( lmem ) then
       inform = - 92
       subname = 'VEVALH'
       call reperr(inform,subname)
       return
    end if

    do i = 1,hnnz
       if ( hrow(i) .lt. 1 .or. hrow(i) .gt. n .or. &
            hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. &
            hcol(i) .gt. hrow(i) ) then

          if ( iprintctl(3) ) then
             write(* ,200)
             write(* ,400) n,i,hrow(i),hcol(i),hval(i)
             write(10,200)
             write(10,400) n,i,hrow(i),hcol(i),hval(i)
          end if

          hrow(i) = 1
          hcol(i) = 1
          hval(i) = 0.0d0
       end if

       if ( flag .eq. 0 .and. .not. IsANumber(hval(i)) ) then
          if ( iprintctl(3) ) then
             write(* ,300)
             write(* ,400) n,i,hrow(i),hcol(i),hval(i)
             write(10,300)
             write(10,400) n,i,hrow(i),hcol(i),hval(i)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 82
             subname = 'VEVALH'
             call reperr(inform,subname)
             return
          end if
       end if
    end do

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALH WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALH WARNING: There is an element out of range, ', &
                'or in the upper triangle, of the Hessian of the ',   &
                'objetive function computed by the user-supplied ',   &
                'subroutine EVALH. It will be ignored.')

300 format(/,1X,'VEVALH WARNING: There is an element whose value may ', &
                'be +Inf, -Inf or NaN in the Hessian of the objetive ', &
                'function computed by the user-supplied subroutine ',   &
                'EVALH.')

400 format(/,1X,'Dimension: ',I16, &
           /,1X,'Position : ',I16, &
           /,1X,'Row      : ',I16, &
           /,1X,'Column   : ',I16, &
           /,1X,'Value    : ',1P,D24.16)

  end subroutine vevalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalc(n,x,ind,c,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ind,n
    integer, intent(inout) :: inform
    real(kind=8), intent(out) :: c

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: flag

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    call evalc(n,x,ind,c,flag)

    eccnt(ind) = eccnt(ind) + 1

    if ( flag .ne. 0 ) then
       if ( iprintctl(3) ) then
          write(* ,100)
          write(10,100)
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 83
          subname = 'VEVALC'
          call reperr(inform,subname)
          return
       end if
    end if

    if ( flag .eq. 0 .and. .not. IsANumber(c) ) then
       if ( iprintctl(3) ) then
          write(* ,200) ind
          write(* ,300) c
          write(10,200) ind
          write(10,300) c
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 83
          subname = 'VEVALC'
          call reperr(inform,subname)
          return
       end if
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALC WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALC WARNING: The value of constraint ',I16,' ', &
                'computed by the user-supplied subroutine EVALC ',  &
                'may be +Inf, -Inf or NaN.')

300 format(/,1X,'Value: ',1P,D24.16)

  end subroutine vevalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: ind,lim,n
    integer, intent(inout) :: inform
    integer, intent(out) :: jcnnz
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    ! LOCAL SCALARS
    logical :: lmem
    integer :: flag,i

     ! LOCAL ARRAYS
    character(len=15) :: subname

   ! FUNCTIONS
    logical :: IsANumber

    if ( jaccoded ) then

       call evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

       ejccnt(ind) = ejccnt(ind) + 1

       if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
             write(* ,100)
             write(10,100)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 84
             subname = 'VEVALJAC'
             call reperr(inform,subname)
             return
          end if
       end if
       
       if ( lmem ) then
          inform = - 92
          subname = 'VEVALJAC'
          call reperr(inform,subname)
          return
       end if

       do i = 1,jcnnz
          if ( jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then
             if ( iprintctl(3) ) then
                write(* ,200) ind
                write(* ,400) n,i,jcvar(i),jcval(i)
                write(10,200) ind
                write(10,400) n,i,jcvar(i),jcval(i)
             end if

             jcvar(i) = 1
             jcval(i) = 0.0d0
          end if

          if ( flag .eq. 0 .and. .not. IsANumber(jcval(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,300) ind
                write(* ,400) n,i,jcvar(i),jcval(i)
                write(10,300) ind
                write(10,400) n,i,jcvar(i),jcval(i)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 84
                subname = 'VEVALJAC'
                call reperr(inform,subname)
                return
             end if
          end if
       end do

    else
       call ivevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,inform)
       if ( inform .ne. 0 ) return
    end if
 
    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALJAC WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALJAC WARNING: There is an element out of ',      &
                'range in the gradient of constraint ',I16,' ',       &
                'computed by the user-supplied subroutine EVALJAC. ', &
                'It will be ignored.')

300 format(/,1X,'VEVALJAC WARNING: There is an element whose value ',      &
                'may be +Inf, -Inf or NaN in the gradient of constraint ', &
                I16,'computed by the user-supplied subroutine ',           &
                'EVALJAC.')

400 format(/,1X,'Dimension: ',I16, &
           /,1X,'Position : ',I16, &
           /,1X,'Variable : ',I16, &
           /,1X,'Value    : ',1P,D24.16)

  end subroutine vevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine ivevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,inform)

    use modmachconst, only: macheps13

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: ind,lim,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcvar(n)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: jcval(n)

    ! LOCAL SCALARS
    integer      :: i
    real(kind=8) :: cminus,cplus,step,val

    ! LOCAL ARRAYS
    character(len=15) :: subname
    real(kind=8)      :: xtmp(n)

    xtmp(1:n) = x(1:n)

    jcnnz = 0

    do i = 1,n
       step = macheps13 * max( 1.0d0, abs( x(i) ) )

       xtmp(i) = x(i) + step
       call vsetp(n,xtmp,inform)
       if ( inform .ne. 0 ) return
       call vevalc(n,xtmp,ind,cplus,inform)
       if ( inform .ne. 0 ) return

       xtmp(i) = x(i) - step
       call vsetp(n,xtmp,inform)
       if ( inform .ne. 0 ) return
       call vevalc(n,xtmp,ind,cminus,inform)
       if ( inform .ne. 0 ) return

       val = ( cplus - cminus ) / ( 2.0d0 * step )

       if ( val .ne. 0.0d0 ) then
          if ( jcnnz + 1 .gt. lim ) then
             inform = - 92
             subname = 'IVEVALJAC'
             call reperr(inform,subname)
             return
          end if

          jcnnz = jcnnz + 1
          jcvar(jcnnz) = i
          jcval(jcnnz) = val
       end if

       xtmp(i) = x(i)
    end do

  end subroutine ivevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalhc(n,x,ind,hrow,hcol,hval,hnnz,lim,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: ind,lim,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hcol(lim),hrow(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    logical :: lmem
    integer :: flag,i

     ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    call evalhc(n,x,ind,hrow,hcol,hval,hnnz,lim,lmem,flag)

    ehccnt(ind) = ehccnt(ind) + 1

    if ( flag .ne. 0 ) then
       if ( iprintctl(3) ) then
          write(* ,100)
          write(10,100)
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 85
          subname = 'VEVALHC'
          call reperr(inform,subname)
          return
       end if
    end if

    if ( lmem ) then
       inform = - 92
       subname = 'VEVALHC'
       call reperr(inform,subname)
       return
    end if

    do i = 1,hnnz
       if ( hrow(i) .lt. 1 .or. hrow(i) .gt. n .or. &
            hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. &
            hcol(i) .gt. hrow(i) ) then

          if ( iprintctl(3) ) then
             write(* ,200) ind
             write(* ,400) n,i,hrow(i),hcol(i),hval(i)
             write(10,200) ind
             write(10,400) n,i,hrow(i),hcol(i),hval(i)
          end if

          hrow(i) = 1
          hcol(i) = 1
          hval(i) = 0.0d0
       end if

       if ( flag .eq. 0 .and. .not. IsANumber(hval(i)) ) then
          if ( iprintctl(3) ) then
             write(* ,300) ind
             write(* ,400) n,i,hrow(i),hcol(i),hval(i)
             write(10,300) ind
             write(10,400) n,i,hrow(i),hcol(i),hval(i)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 85
             subname = 'VEVALHC'
             call reperr(inform,subname)
             return
          end if
       end if
    end do

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALHC WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALHC WARNING: There is an element out of range ', &
                'or in the upper triangle of the Hessian of ',        &
                'constraint ',I16,' computed by the user-supplied ',  &
                'subroutine EVALHC. It will be ignored.')

300 format(/,1X,'VEVALHC WARNING: There is an element whose value ',      &
                'may be +Inf, -Inf or NaN in the Hessian of constraint ', &
                I16,' computed by the user-supplied subroutine ',         &
                'EVALHC.')

400 format(/,1X,'Dimension: ',I16, &
           /,1X,'Position : ',I16, &
           /,1X,'Row      : ',I16, &
           /,1X,'Column   : ',I16, &
           /,1X,'Value    : ',1P,D24.16)

  end subroutine vevalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalfc(n,x,f,m,c,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    integer :: flag,i

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    call evalfc(n,x,f,m,c,flag)

    efccnt = efccnt + 1
    fcnt   = fcnt   + 1

    if ( flag .ne. 0 ) then
       if ( iprintctl(3) ) then
          write(* ,100)
          write(10,100)
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 86
          subname = 'VEVALFC'
          call reperr(inform,subname)
          return
       end if
    end if

    if ( flag .eq. 0 .and. .not. IsANumber(f) ) then
       if ( iprintctl(3) ) then
          write(* ,200)
          write(* ,400) f
          write(10,200)
          write(10,400) f
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 86
          subname = 'VEVALFC'
          call reperr(inform,subname)
          return
       end if
    end if

    do i = 1,m
       if ( flag .eq. 0 .and. .not. IsANumber(c(i)) ) then
          if ( iprintctl(3) ) then
             write(* ,300)
             write(* ,500) n,m,i,c(i)
             write(10,300)
             write(10,500) n,m,i,c(i)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 86
             subname = 'VEVALFC'
             call reperr(inform,subname)
             return
          end if
       end if
    end do

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALFC WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALFC WARNING: The objective function value ',   &
                'computed by the user-supplied',/,1X,'subroutine ', &
                'EVALFC may be +Inf, -Inf or NaN.')

300 format(/,1X,'VEVALFC WARNING: The value of a constraint ',      &
                'computed by the user-supplied',/,1X,'subroutine ', &
                'EVALFC may be +Inf, -Inf or NaN.')

400 format(/,1X,'Value: ',1P,D24.16)

500 format(/,1X,'Dimension of the space: ',I16, &
           /,1X,'Number of constraints : ',I16, &
           /,1X,'Constraint            : ',I16, &
           /,1X,'Value                 : ',1P,D24.16)

  end subroutine vevalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: lim,m,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    ! LOCAL SCALARS
    logical :: lmem
    integer :: flag,i

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    if ( gjaccoded ) then

       call evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

       egjccnt = egjccnt + 1

       if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
             write(* ,100)
             write(10,100)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 87
             subname = 'VEVALGJAC'
             call reperr(inform,subname)
             return
          end if
       end if

       if ( lmem ) then
          inform = - 92
          subname = 'VEVALGJAC'
          call reperr(inform,subname)
          return
       end if

       do i = 1,n
          if ( flag .eq. 0 .and. .not. IsANumber(g(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,300)
                write(* ,400) n,i,g(i)
                write(10,300)
                write(10,400) n,i,g(i)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 87
                subname = 'VEVALGJAC'
                call reperr(inform,subname)
                return
             end if
          end if
       end do

       do i = 1,jcnnz
          if ( jcfun(i) .lt. 1 .or. jcfun(i) .gt. m .or. &
               jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then

             if ( iprintctl(3) ) then
                write(* ,200)
                write(* ,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
                write(10,200)
                write(10,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
             end if

             jcfun(i) = 1
             jcvar(i) = 1
             jcval(i) = 0.0d0
          end if

          if ( flag .eq. 0 .and. .not. IsANumber(jcval(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,300)
                write(* ,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
                write(10,300)
                write(10,500) n,m,i,jcfun(i),jcvar(i),jcval(i)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 87
                subname = 'VEVALGJAC'
                call reperr(inform,subname)
                return
             end if
          end if
       end do

    else
       call ivevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)
       if ( inform .ne. 0 ) return
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALGJAC WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALGJAC WARNING: There is an element out of ',      &
                'range in the gradient of the objective function or ', &
                'in the Jacobian of the constraints computed by the ', &
                'user-supplied subroutine EVALGJAC. It will be ',      &
                'ignored.')

300 format(/,1X,'VEVALGJAC WARNING: There is an element whose value ', &
                'may be +Inf, -Inf or NaN in the',/,1X,'gradient of ', &
                'the objective function or in the Jacobian of the ',   &
                'constraints',/,1X,'computed by the user-supplied ',   &
                'subroutine EVALGJAC.')

400 format(/,1X,'Dimension of the space: ',I16, &
           /,1X,'Position              : ',I16, &
           /,1X,'Value                 : ',1P,D24.16)

500 format(/,1X,'Dimension of the space: ',I16, &
           /,1X,'Number of constraints : ',I16, &
           /,1X,'Position              : ',I16, &
           /,1X,'Constraint            : ',I16, &
           /,1X,'Variable              : ',I16, &
           /,1X,'Value                 : ',1P,D24.16)

  end subroutine vevalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine ivevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)

    use modmachconst, only: macheps13

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: lim,m,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    ! LOCAL SCALARS
    integer :: i,j
    real(kind=8) :: fminus,fplus,step,val

    ! LOCAL ARRAYS
    character(len=15) :: subname
    real(kind=8)      :: cminus(m),cplus(m),xtmp(n)

    xtmp(1:n) = x(1:n)

    jcnnz = 0

    do i = 1,n
       step = macheps13 * max( 1.0d0, abs( x(i) ) )

       xtmp(i) = x(i) + step
       call vsetp(n,xtmp,inform)
       if ( inform .ne. 0 ) return
       call vevalfc(n,xtmp,fplus,m,cplus,inform)
       if ( inform .ne. 0 ) return

       xtmp(i) = x(i) - step
       call vsetp(n,xtmp,inform)
       if ( inform .ne. 0 ) return
       call vevalfc(n,xtmp,fminus,m,cminus,inform)
       if ( inform .ne. 0 ) return

       do j = 1,m
          val = ( cplus(j) - cminus(j) ) / ( 2.0d0 * step )

          if ( val .ne. 0.0d0 ) then
             if ( jcnnz + 1 .gt. lim ) then
                inform = - 92
                subname = 'IVEVALGJAC'
                call reperr(inform,subname)
                return
             end if

             jcnnz = jcnnz + 1
             jcfun(jcnnz) = j
             jcvar(jcnnz) = i
             jcval(jcnnz) = val
          end if
       end do

       g(i) = ( fplus - fminus ) / ( 2.0d0 * step )

       xtmp(i) = x(i)
    end do

  end subroutine ivevalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalgjacp(n,x,g,m,p,q,work,gotj,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    logical,   intent(inout) :: gotj
    integer,   intent(in)    :: m,n
    integer,   intent(inout) :: inform
    character, intent(in)    :: work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)    :: x(n)
    real(kind=8), intent(inout) :: p(m),q(n)
    real(kind=8), intent(out)   :: g(n)

    ! The meaning of argument work follows: work = 'j' or 'J' means that 
    ! q is an input array argument and that p = Jacobian x q must be 
    ! computed, while work = 't' or 'T' means that p is an input array 
    ! argument and that q = Jacobian^t x p must be computed. Moreover, a 
    ! capital letter (i.e. 'J' or 'T') means that the gradient of the 
    ! objective function g must also be computed. A lower letter (i.e. 
    ! 'j' or 't') means that only the product with the Jacobian is 
    ! required. In the later case, input array argument g MUST NOT be 
    ! modified nor referenced.

    ! Nov, 19, 2015: In the cases work = 'j' or 'J', a new meaning of
    ! the output array parameter p(1:m) was added. On input, each
    ! position j of p may have a null or a non-null value. Only
    ! positions j that on input, has a non-null value need to be
    ! computed.
    
    ! LOCAL SCALARS
    integer :: flag,i,j

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    call evalgjacp(n,x,g,m,p,q,work,gotj,flag)

    egjcpcnt = egjcpcnt + 1

    if ( flag .ne. 0 ) then
       if ( iprintctl(3) ) then
          write(* ,100)
          write(10,100)
       end if

       nsf = nsf + 1

       if ( safemode .and. nsf .gt. nsfmax ) then
          inform = 88
          subname = 'VEVALGJACP'
          call reperr(inform,subname)
          return
       end if
    end if

    if ( work .eq. 'J' .or. work .eq. 'T' ) then
       do i = 1,n
          if ( flag .eq. 0 .and. .not. IsANumber(g(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,200)
                write(* ,400) n,i,g(i)
                write(10,200)
                write(10,400) n,i,g(i)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 88
                subname = 'VEVALGJACP'
                call reperr(inform,subname)
                return
             end if
          end if
       end do
    end if

    if ( work .eq. 'j' .or. work .eq. 'J' ) then
       do j = 1,m
          if ( flag .eq. 0 .and. .not. IsANumber(p(j)) ) then
             if ( iprintctl(3) ) then
                write(* ,300)
                write(* ,400) m,j,p(j)
                write(10,300)
                write(10,400) m,j,p(j)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 88
                subname = 'VEVALGJACP'
                call reperr(inform,subname)
                return
             end if
          end if
       end do

    else ! if ( work .eq. 't' .or. work .eq. 'T' ) then
       do i = 1,n
          if ( flag .eq. 0 .and. .not. IsANumber(q(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,300)
                write(* ,400) n,i,q(i)
                write(10,300)
                write(10,400) n,i,q(i)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 88
                subname = 'VEVALGJACP'
                call reperr(inform,subname)
                return
             end if
          end if
       end do
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALGJACP WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALGJACP WARNING: There is an element whose ',        &
                'value may be +Inf, -Inf or NaN in the gradient of ',    &
                'the objective function computed by the user-supplied ', &
                'subroutine EVALGJACP.')

300 format(/,1X,'VEVALGJACP WARNING: There is an element whose ',      &
                'value may be +Inf, -Inf or NaN in the product of ',   &
                'the Jacobian (or transpose Jacobian) times a given ', &
                'vector, computed by the user-supplied subroutine ',   &
                'EVALGJACP.')

400 format(/,1X,'Dimension of computed array: ',I16, &
           /,1X,'Position                   : ',I16, &
           /,1X,'Value                      : ',1P,D24.16)

  end subroutine vevalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalhl(n,x,m,lambda,sf,sc,hrow,hcol,hval,hnnz,lim,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: lim,m,n
    integer,      intent(inout) :: inform
    integer,      intent(out)   :: hnnz
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hrow(lim),hcol(lim)
    real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    logical :: lmem
    integer :: flag,i

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    if ( hlcoded ) then

       call evalhl(n,x,m,lambda,sf,sc,hrow,hcol,hval,hnnz,lim,lmem,flag)

       ehlcnt = ehlcnt + 1

       if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
             write(* ,100)
             write(10,100)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 89
             subname = 'VEVALHL'
             call reperr(inform,subname)
             return
          end if
       end if

       if ( lmem ) then
          inform = - 92
          subname = 'VEVALHL'
          call reperr(inform,subname)
          return
       end if

       do i = 1,hnnz
          if ( hrow(i) .lt. 1 .or. hrow(i) .gt. n .or. &
               hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. &
               hcol(i) .gt. hrow(i) ) then

             if ( iprintctl(3) ) then
                write(* ,200)
                write(* ,400) n,i,hrow(i),hcol(i),hval(i)
                write(10,200)
                write(10,400) n,i,hrow(i),hcol(i),hval(i)
             end if

             hrow(i) = 1
             hcol(i) = 1
             hval(i) = 0.0d0
          end if

          if ( flag .eq. 0 .and. .not. IsANumber(hval(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,300)
                write(* ,400) n,i,hrow(i),hcol(i),hval(i)
                write(10,300)
                write(10,400) n,i,hrow(i),hcol(i),hval(i)
             end if

             nsf = nsf + 1
             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 89
                subname = 'VEVALHL'
                call reperr(inform,subname)
                return
             end if
          end if
       end do

    else if ( hcoded .and. ( hccoded .or. m .eq. 0 ) ) then
       call ivevalhl(n,x,m,lambda,sf,sc,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALHL WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALHL WARNING: There is an element out of range, ', &
                'or in the upper triangle, of the',/,1X,'Hessian of ', &
                'the Lagrangian computed by the user-supplied ',       &
                'subroutine EVALHL. It',/,1X,'will be ignored.')

300 format(/,1X,'VEVALHL WARNING: There is an element whose value ',  &
                'may be +Inf, -Inf or NaN in the',/,1X,'Hessian of ', &
                'the Lagrangian computed by the user-supplied ',      &
                'subroutine EVALHL.')

400 format(/,1X,'Dimension: ',I16, &
           /,1X,'Position : ',I16, &
           /,1X,'Row      : ',I16, &
           /,1X,'Column   : ',I16, &
           /,1X,'Value    : ',1P,D24.16)

  end subroutine vevalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine ivevalhl(n,x,m,lambda,sf,sc,hrow,hcol,hval,hnnz,lim,inform)

    use modalgparam, only: ignoref

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: lim,m,n
    integer,      intent(inout) :: inform
    integer,      intent(out)   :: hnnz
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hrow(lim),hcol(lim)
    real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    integer :: col,con,hnnztmp,i,ind,itmp,j,row,nextj,rnnz
    real(kind=8) :: val

    ! LOCAL ARRAYS
    integer :: hcon(lim),pos(n),rind(n),strow(n)
    real(kind=8) :: rval(n)

    ! ==================================================================
    ! COMPUTE HESSIANS
    ! ==================================================================

    ! COMPUTE HESSIAN OF THE OBJECTIVE FUNCTION

    if ( ignoref ) then
       hnnz = 0

    else
       ! Compute Hessian of the objective function
       call vevalh(n,x,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

       ! For each element of the Hessian of the objective function,
       ! set constraint index as zero
       hval(1:hnnz) = hval(1:hnnz) * sf
       hcon(1:hnnz) = 0
    end if

    if ( m .eq. 0 ) return

    ! COMPUTE HESSIANS OF THE CONSTRAINTS

    ind = 0

    do j = 1,m
       ! The 'if' statement below was added on January 4, 2013 to
       ! avoid computing the Hessian of a constraint with a null 
       ! Lagrange multiplier.
       if ( lambda(j) .ne. 0.0d0 ) then
          ! Compute Hessian of constraint j
          call vevalhc(n,x,j,hrow(hnnz+ind+1:lim),hcol(hnnz+ind+1:lim), &
               hval(hnnz+ind+1:lim),hnnztmp,lim-hnnz-ind,inform)
          if ( inform .ne. 0 ) return

          ! For each element of the Hessian, set constraint as j
          hval(hnnz+ind+1:hnnz+ind+hnnztmp) = hval(hnnz+ind+1:hnnz+ind+hnnztmp) * sc(j)
          hcon(hnnz+ind+1:hnnz+ind+hnnztmp) = j

          ind = ind + hnnztmp
       end if
    end do

    if ( ind .eq. 0 ) return

    hnnz = hnnz + ind

    ! ==================================================================
    ! SET ROW LINKED LISTS
    ! ==================================================================

    ! Initialize pointers to the first element of each row
    strow(1:n) = 0

    ! Set row linked lists
    do i = 1,hnnz
       row        = hrow(i)
       itmp       = strow(row)
       strow(row) = i
       hrow(i)    = itmp
    end do

    ! ==================================================================
    ! BUILD HESSIAN OF THE LAGRANGIAN ROW BY ROW
    ! ==================================================================

    ! Initialize array pos
    pos(1:n) = 0

    do i = 1,n
       ! Initialize the i-th row of the Hessian of the Lagrangian
       rnnz = 0

       ! Process the i-th row of all the Hessians
       j = strow(i)

10     if ( j .ne. 0 ) then

          ! Process element (i,hcol(j)) of the Hessian of constraint
          ! hcon(j) (Hessian of the objective function if hcon(j)=0)

          col = hcol(j)
          con = hcon(j)
          if ( con .eq. 0 ) then
             val = hval(j)
          else
             val = hval(j) * lambda(con)
          end if

          if ( pos(col) .ne. 0 ) then
             rval(pos(col)) = rval(pos(col)) + val

          else
             rnnz           = rnnz + 1
             pos(col)       = rnnz
             rind(pos(col)) = col
             rval(pos(col)) = val
          end if

          ! Get next element in the i-th row linked list
          j = hrow(j)
          go to 10
       end if

       ! Clean array pos
       pos(rind(1:rnnz)) = 0

       ! Set i-th row of hl (over the i-th rows of the Hessians)
       ! and mark remaining elements to be deleted
       j = strow(i)

20     if ( j .ne. 0 ) then
          nextj = hrow(j)

          if ( rnnz .ne. 0 ) then
             hrow(j) = i
             hcol(j) = rind(rnnz)
             hval(j) = rval(rnnz)
             rnnz    = rnnz - 1
          else
             hrow(j) = 0
          end if

          j = nextj
          go to 20
       end if

    end do

    ! Eliminate remaining elements (marked with hrow(j)=0)
    j = 1

30  if ( j .le. hnnz ) then
       if ( hrow(j) .eq. 0 ) then
          if ( j .ne. hnnz ) then
             hrow(j) = hrow(hnnz)
             hcol(j) = hcol(hnnz)
             hval(j) = hval(hnnz)
          end if
          hnnz = hnnz - 1
       else
          j = j + 1
       end if

       go to 30
    end if

  end subroutine ivevalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine vevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

    use probgiven
    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(inout) :: gothl
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
    real(kind=8), intent(out) :: hp(n)

    ! LOCAL SCALARS
    integer :: flag,i

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! FUNCTIONS
    logical :: IsANumber

    if ( hlpcoded ) then
       call evalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,flag)

       ehlpcnt = ehlpcnt + 1

       if ( flag .ne. 0 ) then
          if ( iprintctl(3) ) then
             write(* ,100)
             write(10,100)
          end if

          nsf = nsf + 1

          if ( safemode .and. nsf .gt. nsfmax ) then
             inform = 90
             subname = 'VEVALHLP'
             call reperr(inform,subname)
             return
          end if
       end if

       do i = 1,n
          if ( flag .eq. 0 .and. .not. IsANumber(hp(i)) ) then
             if ( iprintctl(3) ) then
                write(* ,200)
                write(* ,300) n,i,hp(i)
                write(10,200)
                write(10,300) n,i,hp(i)
             end if

             nsf = nsf + 1

             if ( safemode .and. nsf .gt. nsfmax ) then
                inform = 90
                subname = 'VEVALHLP'
                call reperr(inform,subname)
                return
             end if
          end if
       end do

    else if ( seconde ) then
       call ivevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
       if ( inform .ne. 0 ) return
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'VEVALHLP WARNING: A non-null flag was returned.',/)

200 format(/,1X,'VEVALHLP WARNING: There is an element in the ',     &
                'product of the Hessian of the Lagrangian by a ',    &
                'computed by the user-supplied subroutine EVALHLP ', &
                'whose value may be +Inf, -Inf or NaN.')

300 format(/,1X,'Dimension of the space: ',I16, &
           /,1X,'Position              : ',I16, &
           /,1X,'Value                 : ',1P,D24.16)

  end subroutine vevalhlp

  ! ******************************************************************
  ! ******************************************************************

  subroutine ivevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

    use probgiven, only: hnnzlim
    use modsvdhess

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(inout) :: gothl
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
    real(kind=8), intent(out) :: hp(n)

    ! LOCAL SCALARS
    integer :: col,i,row
    real(kind=8) :: val

    if ( .not. gothl ) then
       gothl = .true.

       call vevalhl(n,x,m,lambda,sf,sc,svdhrow,svdhcol,svdhval, &
            svdhnnz,hnnzlim,inform)
       if ( inform .ne. 0 ) return
    end if

    hp(1:n) = 0.0d0

    do i = 1,svdhnnz
       row = svdhrow(i)
       col = svdhcol(i)
       val = svdhval(i)

       hp(row) = hp(row) + p(col) * val

       if ( row .ne. col ) then
          hp(col) = hp(col) + p(row) * val
       end if
    end do

  end subroutine ivevalhlp

  ! ******************************************************************
  ! ******************************************************************

  subroutine reperr(inform,subname)

    use modouttyp, only: iprintctl

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: inform

    ! ARRAY ARGUMENTS
    character(len=15), intent(in) :: subname

    if ( iprintctl(3) ) then
       write(* ,100) inform,subname
       write(10,100) inform,subname
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'*** Error code = ',I3,' detected in subroutine ',A15,15X,                '***',/, &
           /,1X,'*** Codes meaning:                                                        ***',/, &
           /,1X,'***  80) Non-null flag (or NaN/Inf) in user-supplied evalf subroutine     ***', &
           /,1X,'***  81) Non-null flag (or NaN/Inf) in user-supplied evalg subroutine     ***', &
           /,1X,'***  82) Non-null flag (or NaN/Inf) in user-supplied evalh subroutine     ***', &
           /,1X,'***  83) Non-null flag (or NaN/Inf) in user-supplied evalc subroutine     ***', &
           /,1X,'***  84) Non-null flag (or NaN/Inf) in user-supplied evaljac subroutine   ***', &
           /,1X,'***  85) Non-null flag (or NaN/Inf) in user-supplied evalhc subroutine    ***', &
           /,1X,'***  86) Non-null flag (or NaN/Inf) in user-supplied evalfc subroutine    ***', &
           /,1X,'***  87) Non-null flag (or NaN/Inf) in user-supplied evalgjac subroutine  ***', &
           /,1X,'***  88) Non-null flag (or NaN/Inf) in user-supplied evalgjacp subroutine ***', &
           /,1X,'***  89) Non-null flag (or NaN/Inf) in user-supplied evalhl subroutine    ***', &
           /,1X,'***  90) Non-null flag (or NaN/Inf) in user-supplied evalhlp subroutine   ***', &
           /,1X,'*** -91) Mandatory subroutines were not properly provided                 ***', &
           /,1X,'*** -92) Lack of memory (increase jcnnzmax or hnnzmax)                    ***', &
           /,1X,'*** -93) Memory allocation error                                          ***', &
           /,1X,'*** -94) Memory deallocation error                                        ***', &
           /,1X,'*** -95) Bound constraints are infeasible                                 ***',/, &
           /,1X,'*** Terminating execution. ***',/)

  end subroutine reperr

  ! ******************************************************************
  ! ******************************************************************

  subroutine vsetp(n,x,inform)

    use modsvdgrad, only: svdgotc

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    svdgotc = .false.

    call setp(n,x)

  end subroutine vsetp

  ! ******************************************************************
  ! ******************************************************************

  subroutine vunsetp()

    use modsvdgrad, only: svdgotc

    implicit none

    svdgotc = .false.

    call unsetp()

  end subroutine vunsetp

  ! ******************************************************************
  ! ******************************************************************

  subroutine checkd(n,xini,l,u,m,lambdaini,inform)

    use modmachconst
    use probgiven, only: gcoded,gjaccoded,gjacpcoded,hcoded,hccoded, &
         hlcoded,hlpcoded,jaccoded

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: l(n),u(n),xini(n),lambdaini(m)

    ! This subrotutine checks the user supplied first and second
    ! derivatives subroutines (evalg, evalh, evaljac and evalhc) for
    ! computing the objective function gradient and Hessian and the
    ! constraints gradients and Hessians, respectively.

    ! LOCAL SCALARS
    character    :: answer
    integer      :: i,j
    real(kind=8) :: fakesf,seed

    ! LOCAL ARRAYS
    real(kind=8) :: x(n),lambda(m),fakesc(m)

    ! FUNCTIONS
    real(kind=8) :: drand

    ! SET A PERTURBATIONS OF X AND LAMBDA (INITIAL GUESS) AND FAKE SCALING FACTORS

    seed = 123456.0d0
    do i = 1,n
       if ( l(i) .lt. xini(i) .and. xini(i) .lt. u(i) ) then
          x(i) = xini(i) + macheps12 * ( 2.0d0 * drand(seed) - 1.0d0 ) * &
               max( 1.0d0, abs( xini(i) ) )
       else if ( xini(i) .eq. l(i) ) then
          x(i) = xini(i) + macheps12 * drand(seed) * max( 1.0d0, abs( xini(i) ) )
       else
          x(i) = xini(i) - macheps12 * drand(seed) * max( 1.0d0, abs( xini(i) ) )
       end if
       x(i) = max( l(i), min( x(i), u(i) ) )
    end do

    do j = 1,m
       lambda(j) = lambdaini(j) + macheps12 * ( 2.0d0 * drand(seed) - 1.0d0 ) * &
            max( 1.0d0, abs( lambdaini(j) ) )
    end do

    fakesf = drand(seed)
    do j = 1,m
       fakesc(j) = drand(seed)
    end do

    write(* ,100)
    write(10,100)

    do i = 1,n
       write(* ,110) i,x(i)
       write(10,110) i,x(i)
    end do

    if ( m .gt. 0 ) then

       write(* ,101)
       write(10,101)
       
       do j = 1,m
          write(* ,111) j,lambda(j)
          write(10,111) j,lambda(j)
       end do

    end if

    ! CHECK OBJECTIVE FUNCTION GRADIENT

    if ( .not. gcoded ) then
       write(* ,120) 'evalg'
       write(10,120) 'evalg'

       go to 1000
    end if

    write(* ,130) 'evalg'
    write(10,130) 'evalg'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1000

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
       call checkg(n,x,inform)
       if ( inform .ne. 0 ) return
    end if

    ! CHECK JACOBIAN OF CONSTRAINTS

1000 continue

    if ( .not. jaccoded ) then
       write(* ,120) 'evaljac'
       write(10,120) 'evaljac'

       go to 1020
    end if

    write(* ,130) 'evaljac'
    write(10,130) 'evaljac'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1020

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then

       j = 1

1010   if ( j .le. m ) then

          write(* ,220) j
          write(10,220) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
             return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
             go to 1020

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
             call checkjac(n,x,j,inform)
             if ( inform .ne. 0 ) return
          end if

          j = j + 1

          go to 1010

       end if

    end if

    ! CHECK HESSIAN OF THE OBJECTIVE FUNCTION

1020 continue

    if ( .not. hcoded ) then
       write(* ,120) 'evalh'
       write(10,120) 'evalh'

       go to 1030
    end if

    write(* ,130) 'evalh'
    write(10,130) 'evalh'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1030

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
       call checkh(n,x,inform)
       if ( inform .ne. 0 ) return
    end if

    ! CHECK HESSIANS OF THE CONSTRAINTS

1030 continue

    if ( .not. hccoded ) then
       write(* ,120) 'evalhc'
       write(10,120) 'evalhc'

       go to 1050
    end if

    write(* ,130) 'evalhc'
    write(10,130) 'evalhc'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1050

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then

       j = 1

1040   if ( j .le. m ) then

          write(* ,230) j
          write(10,230) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
             return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
             go to 1050

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
             call checkhc(n,x,j,inform)
             if ( inform .ne. 0 ) return
          end if

          j = j + 1

          go to 1040

       end if

    end if

    ! CHECK GRADIENT OF OBJECTIVE FUNCTION PLUS JACOBIAN OF CONSTRAINTS

1050 continue

    if ( .not. gjaccoded ) then
       write(* ,120) 'evalgjac'
       write(10,120) 'evalgjac'

       go to 1060
    end if

    write(* ,130) 'evalgjac'
    write(10,130) 'evalgjac'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1060

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
       call checkgjac(n,m,x,inform)
       if ( inform .ne. 0 ) return
    end if

    ! CHECK GRADIENT OF OBJECTIVE FUNCTION PLUS PRODUCT OF 
    ! JACOBIAN OF CONSTRAINTS TIMES A GIVEN VECTOR

1060 continue

    if ( .not. gjacpcoded ) then
       write(* ,120) 'evalgjacp'
       write(10,120) 'evalgjacp'

       go to 1070
    end if

    write(* ,130) 'evalgjacp'
    write(10,130) 'evalgjacp'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1070

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
       write(* ,*) 'Test of evalgjacp not implemented yet!'
       write(10,*) 'Test of evalgjacp not implemented yet!'
    end if

    ! CHECK HESSIAN OF THE LAGRANGIAN

1070 continue

    if ( .not. hlcoded ) then
       write(* ,120) 'evalhl'
       write(10,120) 'evalhl'

       go to 1080
    end if

    write(* ,130) 'evalhl'
    write(10,130) 'evalhl'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1080

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
       call checkhl(n,m,x,lambda,fakesf,fakesc,inform)
       if ( inform .ne. 0 ) return
       ! write(* ,*) 'Test of evalhl not implemented yet!'
       ! write(10,*) 'Test of evalhl not implemented yet!'
    end if

    ! CHECK HESSIAN OF THE LAGRANGIAN TIMES A VECTOR

1080 continue

    if ( .not. hlpcoded ) then
       write(* ,120) 'evalhlp'
       write(10,120) 'evalhlp'

       go to 1090
    end if

    write(* ,130) 'evalhlp'
    write(10,130) 'evalhlp'

    read(*,*) answer

    if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
       return

    else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
       go to 1090

    else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
       write(* ,*) 'Test of evalhlp not implemented yet!'
       write(10,*) 'Test of evalhlp not implemented yet!'
    end if

    ! NON-EXECUTABLE STATEMENTS

1090 continue

100 format(/,1X,'Derivatives will be tested at the perturbed ', &
         'initial guess: ')
101 format(/,1X,'Perturbed multipliers: ')
110 format(  1X,'x(',I6,') = ',1P,D15.8)
111 format(  1X,'lambda(',I6,') = ',1P,D15.8)
120 format(/,1X,'Skipping checking of uncoded ',A9,' subroutine.')
130 format(/,1X,'Would you like to check subroutine ',A9,'?', &
         /,1X,'Type Y(es), N(o), A(bort checking) or ',     &
         'S(kip checking this subroutine): ')

220 format(/,1X,'Check gradient of constraint ',I5,'?',   &
         /,1X,'Type Y(es), N(o), A(bort checking) or ', &
         'S(kip gradients of constraints): ')
230 format(/,1X,'Check Hessian matrix of constraint ',I5,'?', &
         /,1X,'Type Y(es), N(o), A(bort checking) or ',     &
         'S(kip Hessians of constraints): ')

  end subroutine checkd

  ! *****************************************************************
  ! *****************************************************************

  subroutine checkg(n,x,inform)

    use modmachconst, only: macheps13

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)

    ! This subrotutine checks the user supplied subroutine evalg for
    ! computing the gradient of the objective function using central
    ! finite differences with two different discretization steps.

    ! LOCAL SCALARS
    integer      :: i
    real(kind=8) :: fminus,fplus,gdiff1,gdiff2,maxerr,step1,step2,tmp

    ! LOCAL ARRAYS
    real(kind=8) :: g(n)

    call vsetp(n,x,inform)
    if ( inform .ne. 0 ) return
    call vevalg(n,x,g,inform)
    if ( inform .ne. 0 ) return

    write(* ,100)
    write(10,100)

    maxerr = 0.0d0

    do i = 1,n
       tmp  = x(i)

       step1 = macheps13 * max( abs( tmp ), 1.0d0 )

       x(i) = tmp + step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalf(n,x,fplus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp - step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalf(n,x,fminus,inform)
       if ( inform .ne. 0 ) return

       gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

       step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

       x(i) = tmp + step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalf(n,x,fplus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp - step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalf(n,x,fminus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp

       gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

       tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

       write(* ,110) i,g(i),gdiff1,gdiff2,tmp
       write(10,110) i,g(i),gdiff1,gdiff2,tmp

       maxerr = max( maxerr, tmp )

    end do

    write(* ,120) maxerr
    write(10,120) maxerr

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'Gradient vector of the objective function.',          &
         /,1X,'Index',13X,'evalg',2X,'Central diff (two different ', &
         'steps)',4X,'Absolute error')
110 format(  1X,I5,4(3X,1P,D15.8))
120 format(  1X,'Maximum absolute error = ',1P,D15.8)

  end subroutine checkg

  ! *****************************************************************
  ! *****************************************************************

  subroutine checkh(n,x,inform)

    use modmachconst, only: macheps12
    use probgiven, only: hnnzlim
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)

    ! This subrotutine checks the user supplied subroutine evalh for
    ! computing the Hessian of the objective function using central
    ! finite differences with two different discretization steps.

    ! LOCAL SCALARS
    logical :: nullcol
    integer :: allocerr,i,j,hnnz
    real(kind=8) :: elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

    ! LOCAL ARRAYS
    character(len=15) :: subname
    integer :: hrow(hnnzlim),hcol(hnnzlim)
    real(kind=8) :: g(n),gplus1(n),gplus2(n),hval(hnnzlim),maxcoe(n)

    ! LOCAL POINTERS
    real(kind=8), pointer :: H(:,:)
    
    ! Check viability of the test

    if ( n .gt. 1000 ) then
       write(* ,100) n
       write(10,100) n
       return
    end if

    ! Get memory
    
    allocate(H(n,n),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 94
       subname = 'CHECKH'
       call reperr(inform,subname)
       return
    end if

    ! Compute the gradient of the objective function at x

    call vsetp(n,x,inform)
    if ( inform .ne. 0 ) return
    call vevalg(n,x,g,inform)
    if ( inform .ne. 0 ) return

    ! Compute the Hessian of the objective function at x and save in a
    ! dense matrix

    call vevalh(n,x,hrow,hcol,hval,hnnz,hnnzlim,inform)
    if ( inform .ne. 0 ) return

    H(1:n,1:n) = 0.0d0

    do i = 1,hnnz
       H(hrow(i),hcol(i)) = H(hrow(i),hcol(i)) + hval(i)
    end do

    ! Test column by column

    write(* ,200)
    write(10,200)

    maxerr = 0.0d0

    do j = 1,n

       tmp  = x(j)

       step1 = macheps12 * max( abs( tmp ), 1.0d0 )

       x(j) = tmp + step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalg(n,x,gplus1,inform)
       if ( inform .ne. 0 ) return

       step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

       x(j) = tmp + step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalg(n,x,gplus2,inform)
       if ( inform .ne. 0 ) return

       x(j) = tmp

       write(* ,210) j
       write(10,210) j

       maxcoe(j) = 0.0d0

       nullcol = .true.

       do i = 1,n
          if ( i .ge. j ) then
             elem = H(i,j)
          else
             elem = H(j,i)
          end if
          hdiff1 = ( gplus1(i) - g(i) ) / step1
          hdiff2 = ( gplus2(i) - g(i) ) / step2
          tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
          if ( elem   .ne. 0.0d0 .or. &
               hdiff1 .ne. 0.0d0 .or. &
               hdiff2 .ne. 0.0d0 ) then
             if ( nullcol ) then
                nullcol = .false.
                write(* ,220)
                write(10,220)
             end if
             write(* ,230) i,elem,hdiff1,hdiff2,tmp
             write(10,230) i,elem,hdiff1,hdiff2,tmp
          end if
          maxcoe(j) = max( maxcoe(j), tmp )
       end do

       maxerr = max( maxerr, maxcoe(j) )

       if ( nullcol ) then
          write(* ,240)
          write(10,240)
       else
          write(* ,250) maxcoe(j)
          write(10,250) maxcoe(j)
       end if

    end do

    write(* ,*)
    write(10,*)

    do j = 1,n
       write(* ,260) j,maxcoe(j)
       write(10,260) j,maxcoe(j)
    end do

    write(* ,270) maxerr
    write(10,270) maxerr

    ! Free memory
    
    deallocate(H,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 94
       subname = 'CHECKH'
       call reperr(inform,subname)
       return
    end if
    
    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'The current implementation of subroutine CHECKH requires a dense ', &
         'n x n',/,1X,'matrix ( n = ',I6,') and allocating this matrix is undesired when ', &
         'n is',/,1X,'greater than 1,000. Therefore, the test will be skipped.')
200 format(/,1X,'Hessian matrix of the objective function column by ', &
         'column.')
210 format(/,1X,'Column:  ',I6)
220 format(/,1X,'Index',13X,'evalh',3X,'Incr. Quoc. (two different ', &
         'steps)',4X,'Absolute error')
230 format(  1X,I5,4(3X,1P,D15.8))
240 format(  1X,'All the elements of this column are null.')
250 format(  1X,'Maximum absolute error = ',1P,D15.8)
260 format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
270 format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

  end subroutine checkh

  ! *****************************************************************
  ! *****************************************************************

  subroutine checkjac(n,x,ind,inform)

    use modmachconst, only: macheps13

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: ind,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)

    ! This subrotutine checks the user supplied subroutine evaljac for
    ! computing the gradients of the constraints using central finite
    ! differences with two different discretization steps.

    ! LOCAL SCALARS
    logical      :: nullcol
    integer      :: i,jcnnz
    real(kind=8) :: cminus,cplus,jacdiff1,jacdiff2,maxerr,step1,step2,tmp

    ! LOCAL ARRAYS
    integer      :: jcvar(n)
    real(kind=8) :: g(n),jcval(n)

    ! COMPUTE THE GRADIENT OF THE CONSTRAINT AND SAVE IT INTO A DENSE
    ! VECTOR

    call vsetp(n,x,inform)
    if ( inform .ne. 0 ) return
    call vevaljac(n,x,ind,jcvar,jcval,jcnnz,n,inform)
    if ( inform .ne. 0 ) return

    do i = 1,n
       g(i) = 0.0d0
    end do

    do i = 1,jcnnz
       g(jcvar(i)) = g(jcvar(i)) + jcval(i)
    end do

    ! COMPARE WITH CENTRAL FINITE DIFFERENCES

    write(* ,100) ind
    write(10,100) ind

    maxerr = 0.0d0

    nullcol = .true.

    do i = 1,n
       tmp  = x(i)

       step1 = macheps13 * max( abs( tmp ), 1.0d0 )

       x(i) = tmp + step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalc(n,x,ind,cplus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp - step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalc(n,x,ind,cminus,inform)
       if ( inform .ne. 0 ) return

       jacdiff1 = ( cplus - cminus ) / ( 2.0d0 * step1 )

       step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

       x(i) = tmp + step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalc(n,x,ind,cplus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp - step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalc(n,x,ind,cminus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp

       jacdiff2 = ( cplus - cminus ) / ( 2.0d0 * step2 )

       tmp = min( abs( g(i) - jacdiff1 ), abs( g(i) - jacdiff2 ) )

       if ( g(i)     .ne. 0.0d0 .or. &
            jacdiff1 .ne. 0.0d0 .or. &
            jacdiff2 .ne. 0.0d0 ) then
          if ( nullcol ) then
             nullcol = .false.
             write(* ,110)
             write(10,110)
          end if
          write(* ,120) i,g(i),jacdiff1,jacdiff2,tmp
          write(10,120) i,g(i),jacdiff1,jacdiff2,tmp
       end if

       maxerr = max( maxerr, tmp )
    end do

    if ( nullcol ) then
       write(* ,130)
       write(10,130)
    else
       write(* ,140) maxerr
       write(10,140) maxerr
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'Gradient vector of constraints ',I5,'.')
110 format(/,1X,'Index',11X,'evaljac',2X,'Central diff (two ', &
         'different steps)',4X,'Absolute error')
120 format(  1X,I5,4(3X,1P,D15.8))
130 format(  1X,'All the elements of this gradient are null.')
140 format(  1X,'Maximum absolute error = ',1P,D15.8)

  end subroutine checkjac

  ! *****************************************************************
  ! *****************************************************************
  subroutine checkhc(n,x,ind,inform)

    use modmachconst, only: macheps12
    use probgiven, only: hnnzlim
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: ind,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)

    ! This subrotutine checks the user supplied subroutine evalhc for
    ! computing the Hessians of the constraints using finite
    ! differences.

    ! LOCAL SCALARS
    logical :: nullcol
    integer :: allocerr,i,j,hnnz,jcnnz
    real(kind=8) :: elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

    ! LOCAL ARRAYS
    character(len=15) :: subname
    integer :: hrow(hnnzlim),hcol(hnnzlim),jcvar(n)
    real(kind=8) :: g(n),gplus1(n),gplus2(n),hval(hnnzlim),jcval(n),maxcoe(n)

    ! LOCAL POINTERS
    real(kind=8), pointer :: H(:,:)
    
    ! Check viability of the test

    if ( n .gt. 1000 ) then
       write(* ,100) n
       write(10,100) n
       return
    end if

    ! Get memory
    
    allocate(H(n,n),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 94
       subname = 'CHECKHC'
       call reperr(inform,subname)
       return
    end if

    ! Compute the gradient of constraint ind at x and save it in a
    ! dense vector

    call vsetp(n,x,inform)
    if ( inform .ne. 0 ) return
    call vevaljac(n,x,ind,jcvar,jcval,jcnnz,n,inform)
    if ( inform .ne. 0 ) return

    do i = 1,n
       g(i) = 0.0d0
    end do

    do i = 1,jcnnz
       g(jcvar(i)) = g(jcvar(i)) + jcval(i)
    end do

    ! Compute the Hessian of constraint ind at x and save it in a
    ! dense matrix

    call vevalhc(n,x,ind,hrow,hcol,hval,hnnz,hnnzlim,inform)
    if ( inform .ne. 0 ) return

    H(1:n,1:n) = 0.0d0

    do i = 1,hnnz
       H(hrow(i),hcol(i)) = H(hrow(i),hcol(i)) + hval(i)
    end do

    write(* ,200) ind
    write(10,200) ind

    maxerr = 0.0d0

    do j = 1,n

       tmp  = x(j)

       ! Compute the gradient of constraint ind at xplus1 and
       ! save in a dense vector

       step1 = macheps12 * max( abs( tmp ), 1.0d0 )

       x(j) = tmp + step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevaljac(n,x,ind,jcvar,jcval,jcnnz,n,inform)
       if ( inform .ne. 0 ) return

       gplus1(1:n) = 0.0d0

       do i = 1,jcnnz
          gplus1(jcvar(i)) = gplus1(jcvar(i)) + jcval(i)
       end do

       ! Compute the gradient of constraint ind at xplus2 and
       ! save in a dense vector

       step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

       x(j) = tmp + step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevaljac(n,x,ind,jcvar,jcval,jcnnz,n,inform)
       if ( inform .ne. 0 ) return

       gplus2(1:n) = 0.0d0

       do i = 1,jcnnz
          gplus2(jcvar(i)) = gplus2(jcvar(i)) + jcval(i)
       end do

       x(j) = tmp

       write(* ,210) j
       write(10,210) j

       maxcoe(j) = 0.0d0

       nullcol = .true.

       do i = 1,n
          if ( i .ge. j ) then
             elem = H(i,j)
          else
             elem = H(j,i)
          end if
          hdiff1 = ( gplus1(i) - g(i) ) / step1
          hdiff2 = ( gplus2(i) - g(i) ) / step2
          tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )

          if ( elem   .ne. 0.0d0 .or. &
               hdiff1 .ne. 0.0d0 .or. &
               hdiff2 .ne. 0.0d0 ) then
             if ( nullcol ) then
                nullcol = .false.
                write(* ,220)
                write(10,220)
             end if
             write(* ,230) i,elem,hdiff1,hdiff2,tmp
             write(10,230) i,elem,hdiff1,hdiff2,tmp
          end if

          maxcoe(j) = max( maxcoe(j), tmp )
       end do

       maxerr = max( maxerr, maxcoe(j) )

       if ( nullcol ) then
          write(* ,240)
          write(10,240)
       else
          write(* ,250) maxcoe(j)
          write(10,250) maxcoe(j)
       end if

    end do

    write(* ,*)
    write(10,*)

    do j = 1,n
       write(* ,260) j,maxcoe(j)
       write(10,260) j,maxcoe(j)
    end do

    write(* ,270) maxerr
    write(10,270) maxerr

    ! Free memory
    
    deallocate(H,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 94
       subname = 'CHECKHC'
       call reperr(inform,subname)
       return
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'The current implementation of subroutine CHECKHC requires a dense ', &
         'n x n',/,1X,'matrix ( n = ',I6,') and allocating this matrix is undesired when ', &
         'n is',/,1X,'greater than 1,000. Therefore, the test will be skipped.')
200 format(/,1X,'Hessian matrix of constraint ',I5,' column by ', &
         'column.')
210 format(/,1X,'Column:  ',I6)
220 format(/,1X,'Index',12X,'evalhc',3X,'Incr. Quoc. (two different ', &
         'steps)',4X,'Absolute error')
230 format(  1X,I5,4(3X,1P,D15.8))
240 format(  1X,'All the elements of this column are null.')
250 format(  1X,'Maximum absolute error = ',1P,D15.8)
260 format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
270 format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

  end subroutine checkhc

  ! *****************************************************************
  ! *****************************************************************

  subroutine checkgjac(n,m,x,inform)

    use probgiven, only: jcnnzlim
    use modmachconst, only: macheps13

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)

    ! This subrotutine checks the user supplied subroutine evalgjac for
    ! computing the gradient of the objective function plus the Jacobian
    ! of the constraints. It uses central finite differences with two 
    ! different discretization steps.

    ! LOCAL SCALARS
    logical      :: nullcol
    integer      :: i,ind,jcnnz
    real(kind=8) :: fminus,fplus,gdiff1,gdiff2,jacdiff1,jacdiff2, &
         maxerr,step1,step2,tmp

    ! LOCAL ARRAYS
    integer      :: jcfun(jcnnzlim),jcvar(jcnnzlim)
    real(kind=8) :: cminus(m),cplus(m),g(n),jcval(jcnnzlim)

    ! COMPUTE GRADIENT OF OBJECTIVE FUNCTION AND JACOBIAN OF 
    ! CONSTRAINTS AT THE CURRENT POINT

    call vsetp(n,x,inform)
    if ( inform .ne. 0 ) return
    call vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform)
    if ( inform .ne. 0 ) return

    ! CHECK GRADIENT OF THE OBJECTIVE FUNCTION

    write(* ,100)
    write(10,100)

    maxerr = 0.0d0

    do i = 1,n
       tmp  = x(i)

       step1 = macheps13 * max( abs( tmp ), 1.0d0 )

       x(i) = tmp + step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalfc(n,x,fplus,m,cplus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp - step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalfc(n,x,fminus,m,cminus,inform)
       if ( inform .ne. 0 ) return

       gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

       step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

       x(i) = tmp + step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalfc(n,x,fplus,m,cplus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp - step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalfc(n,x,fminus,m,cminus,inform)
       if ( inform .ne. 0 ) return

       x(i) = tmp

       gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

       tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

       write(* ,110) i,g(i),gdiff1,gdiff2,tmp
       write(10,110) i,g(i),gdiff1,gdiff2,tmp

       maxerr = max( maxerr, tmp )

    end do

    write(* ,120) maxerr
    write(10,120) maxerr

    ! CHECK JACOBIAN OF CONSTRAINTS

    do ind = 1,m

       ! COPY GRADIENT OF IND-TH CONSTRAINT INTO A DENSE VECTOR

       do i = 1,n
          g(i) = 0.0d0
       end do

       do i = 1,jcnnz
          if ( jcfun(i) .eq. ind ) then
             g(jcvar(i)) = g(jcvar(i)) + jcval(i)
          end if
       end do

       ! COMPARE WITH CENTRAL FINITE DIFFERENCES

       write(* ,200) ind
       write(10,200) ind

       maxerr = 0.0d0

       nullcol = .true.

       do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vsetp(n,x,inform)
          if ( inform .ne. 0 ) return
          call vevalfc(n,x,fplus,m,cplus,inform)
          if ( inform .ne. 0 ) return

          x(i) = tmp - step1
          call vsetp(n,x,inform)
          if ( inform .ne. 0 ) return
          call vevalfc(n,x,fminus,m,cminus,inform)
          if ( inform .ne. 0 ) return

          jacdiff1 = ( cplus(ind) - cminus(ind) ) / ( 2.0d0 * step1)

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vsetp(n,x,inform)
          if ( inform .ne. 0 ) return
          call vevalfc(n,x,fplus,m,cplus,inform)
          if ( inform .ne. 0 ) return

          x(i) = tmp - step2
          call vsetp(n,x,inform)
          if ( inform .ne. 0 ) return
          call vevalfc(n,x,fminus,m,cminus,inform)
          if ( inform .ne. 0 ) return

          x(i) = tmp

          jacdiff2 = ( cplus(ind) - cminus(ind) ) / ( 2.0d0 * step2)

          tmp = min( abs(g(i) - jacdiff1), abs(g(i) - jacdiff2) )

          if ( g(i)     .ne. 0.0d0 .or. &
               jacdiff1 .ne. 0.0d0 .or. &
               jacdiff2 .ne. 0.0d0 ) then
             if ( nullcol ) then
                nullcol = .false.
                write(* ,210)
                write(10,210)
             end if
             write(* ,220) i,g(i),jacdiff1,jacdiff2,tmp
             write(10,220) i,g(i),jacdiff1,jacdiff2,tmp
          end if

          maxerr = max( maxerr, tmp )
       end do

       if ( nullcol ) then
          write(* ,230)
          write(10,230)
       else
          write(* ,120) maxerr
          write(10,120) maxerr
       end if

    end do

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'Gradient vector of the objective function.',   &
         /,1X,'Index',13X,'evalgjac',1X,'Central diff (two ', &
         'different steps)',4X,'Absolute error')
110 format(  1X,I5,4(3X,1P,D15.8))
120 format(  1X,'Maximum absolute error = ',1P,D15.8)
200 format(/,1X,'Gradient vector of constraints ',I5,'.')
210 format(/,1X,'Index',11X,'evalgjac',1X,'Central diff (two ', &
         'different steps)',4X,'Absolute error')
220 format(  1X,I5,4(3X,1P,D15.8))
230 format(  1X,'All the elements of this gradient are null.')

  end subroutine checkgjac

  ! *****************************************************************
  ! *****************************************************************

  subroutine checkhl(n,m,x,lambda,sf,sc,inform)

    use probgiven, only: jcnnzlim,hnnzlim
    use modmachconst, only: macheps12

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m,n
    real(kind=8), intent(in) :: sf
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: sc(m)
    real(kind=8), intent(inout) :: lambda(m),x(n)

    ! This subrotutine checks the user supplied subroutine evalhl for
    ! computing the Hessian of the Lagrangian using finite differences.

    ! LOCAL SCALARS
    logical :: nullcol
    integer :: allocerr,i,j,hlnnz,jcnnz
    real(kind=8) :: elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

    ! LOCAL ARRAYS
    character(len=15) :: subname
    integer :: hlrow(hnnzlim),hlcol(hnnzlim),jcfun(jcnnzlim),jcvar(jcnnzlim)
    real(kind=8) :: g(n),gplus1(n),gplus2(n),hlval(hnnzlim),jcval(jcnnzlim),maxcoe(n)

    ! LOCAL POINTERS
    real(kind=8), pointer :: H(:,:)
    
    ! Check viability of the test

    if ( n .gt. 1000 ) then
       write(* ,100) n
       write(10,100) n
       return
    end if

    ! Get memory
    
    allocate(H(n,n),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 94
       subname = 'CHECKHL'
       call reperr(inform,subname)
       return
    end if

    ! Compute the Hessian of the Lagrangian using the user suppplied
    ! subroutine and save it in a dense matrix.

    call vsetp(n,x,inform)
    if ( inform .ne. 0 ) return
    call vevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,hnnzlim,inform)
    if ( inform .ne. 0 ) return

    H(1:n,1:n) = 0.0d0
    do i = 1,hlnnz
       H(hlrow(i),hlcol(i)) = H(hlrow(i),hlcol(i)) + hlval(i)
    end do

    ! Compute the gradient of the Lagrangian at (x,lambda)

    call vsetp(n,x,inform)
    if ( inform .ne. 0 ) return
    call vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform)
    if ( inform .ne. 0 ) return

    g(1:n) = sf * g(1:n)
    do i = 1,jcnnz
       g(jcvar(i)) = g(jcvar(i)) + sc(jcfun(i)) * lambda(jcfun(i)) * jcval(i)
    end do

    write(* ,200)
    write(10,200)

    maxerr = 0.0d0

    do j = 1,n

       tmp  = x(j)

       step1 = macheps12 * max( abs( tmp ), 1.0d0 )

       x(j) = tmp + step1
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalgjac(n,x,gplus1,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform)
       if ( inform .ne. 0 ) return

       gplus1(1:n) = sf * gplus1(1:n)
       do i = 1,jcnnz
          gplus1(jcvar(i)) = gplus1(jcvar(i)) + sc(jcfun(i)) * lambda(jcfun(i)) * jcval(i)
       end do

       step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

       x(j) = tmp + step2
       call vsetp(n,x,inform)
       if ( inform .ne. 0 ) return
       call vevalgjac(n,x,gplus2,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform)
       if ( inform .ne. 0 ) return

       gplus2(1:n) = sf * gplus2(1:n)
       do i = 1,jcnnz
          gplus2(jcvar(i)) = gplus2(jcvar(i)) + sc(jcfun(i)) * lambda(jcfun(i)) * jcval(i)
       end do

       x(j) = tmp

       write(* ,210) j
       write(10,210) j

       maxcoe(j) = 0.0d0

       nullcol = .true.

       do i = 1,n
          if ( i .ge. j ) then
             elem = H(i,j)
          else
             elem = H(j,i)
          end if
          hdiff1 = ( gplus1(i) - g(i) ) / step1
          hdiff2 = ( gplus2(i) - g(i) ) / step2
          tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )

          if ( elem   .ne. 0.0d0 .or. &
               hdiff1 .ne. 0.0d0 .or. &
               hdiff2 .ne. 0.0d0 ) then
             if ( nullcol ) then
                nullcol = .false.
                write(* ,220)
                write(10,220)
             end if
             write(* ,230) i,elem,hdiff1,hdiff2,tmp
             write(10,230) i,elem,hdiff1,hdiff2,tmp
          end if

          maxcoe(j) = max( maxcoe(j), tmp )
       end do

       maxerr = max( maxerr, maxcoe(j) )

       if ( nullcol ) then
          write(* ,240)
          write(10,240)
       else
          write(* ,250) maxcoe(j)
          write(10,250) maxcoe(j)
       end if

    end do

    write(* ,*)
    write(10,*)

    do j = 1,n
       write(* ,260) j,maxcoe(j)
       write(10,260) j,maxcoe(j)
    end do

    write(* ,270) maxerr
    write(10,270) maxerr

    ! Free memory
    
    deallocate(H,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 94
       subname = 'CHECKHL'
       call reperr(inform,subname)
       return
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'The current implementation of subroutine CHECKHL requires a dense ', &
         'n x n',/,1X,'matrix ( n = ',I6,') and allocating this matrix is undesired when ', &
         'n is',/,1X,'greater than 1,000. Therefore, the test will be skipped.')
200 format(/,1X,'Hessian of the Lagrangian matrix column by column.')
210 format(/,1X,'Column:  ',I6)
220 format(/,1X,'Index',12X,'evalhl',3X,'Incr. Quoc. (two different ', &
         'steps)',4X,'Absolute error')
230 format(  1X,I5,4(3X,1P,D15.8))
240 format(  1X,'All the elements of this column are null.')
250 format(  1X,'Maximum absolute error = ',1P,D15.8)
260 format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
270 format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

  end subroutine checkhl

  ! *****************************************************************
  ! *****************************************************************

end module problvlv
