! DEAL WITH SLACK VARIABLES

module problvlt

  implicit none

  ! SCALARS
  logical, protected  :: slacks
  integer, private :: nws

  ! ARRAYS
  integer, allocatable, private :: slaind(:)

  ! SUBROUTINES
  public :: setproblvlt,tinip,tendp,tevalf,tevalg,tevalh,tevalc, &
       tevaljac,tevalhc,tevalhl,tevalhlp,tevalfc,tevalgjac,tevalgjacp

contains
  
  ! ******************************************************************
  ! ******************************************************************

  subroutine setproblvlt(val_slacks)

    implicit none

    ! SCALARS
    logical, intent(in), optional :: val_slacks

    slacks = .false.
    if ( present( val_slacks ) ) slacks = val_slacks

  end subroutine setproblvlt

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine tinip(n,x,l,u,m,lambda,equatn,linear,coded,inform)

    use modouttyp, only: iprintctl
    use problvlu
    use problvlv, only: reperr
    use probgiven, only: fccoded,jcnnzlim

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m
    integer, intent(inout) :: inform,n

    ! ARRAY ARGUMENTS
    logical,      intent(in)    :: coded(11),linear(m)
    logical,      intent(inout) :: equatn(m)
    real(kind=8), intent(in)    :: lambda(m)
    real(kind=8), pointer       :: l(:),u(:),x(:)

    ! LOCAL SCALARS
    integer      :: allocerr,j,nslacks
    real(kind=8) :: dum

    ! LOCAL ARRAYS
    character(len=15) :: subname
    logical :: ind(m)
    real(kind=8) :: c(m),ltmp(n),utmp(n),xtmp(n)

    if ( slacks .and. all(equatn(1:m)) ) then
       slacks = .false.
       if ( iprintctl(2) ) then
          write(* ,200)
          write(10,200)
       end if
       return
    end if

    if ( slacks ) then

       ! Allocate global arrays in module slacks

       allocate(slaind(m),stat=allocerr)

       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'TINIP'
          call reperr(inform,subname)
          return
       end if

       ! Add slacks

       nws = n
       nslacks = count( .not. equatn(1:m) )

       ! Reset jcnnzlim (maximum number of elements in the Jacobian of the constraints)

       jcnnzlim = jcnnzlim + nslacks

       ! Resize arrays
       
       xtmp(1:n) = x(1:n)
       ltmp(1:n) = l(1:n)
       utmp(1:n) = u(1:n)

       deallocate(x,l,u,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'TINIP'
          call reperr(inform,subname)
          return
       end if

       allocate(x(n+nslacks),l(n+nslacks),u(n+nslacks),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'TINIP'
          call reperr(inform,subname)
          return
       end if

       x(1:n) = xtmp(1:n)
       l(1:n) = ltmp(1:n)
       u(1:n) = utmp(1:n)

       call usetp(n,x,inform)
       if ( inform .ne. 0 ) return

       if ( fccoded ) then
          call uevalfc(n,x,dum,m,c,inform)
          if ( inform .ne. 0 ) return

       else
          ind(1:m) = .not. equatn(1:m)
          call uevalc(n,x,m,ind,c,inform)
          if ( inform .ne. 0 ) return
       end if

       do j = 1,m
          if ( equatn(j) ) then
             slaind(j) = - 1
          else
             equatn(j) = .true.

             n = n + 1
             slaind(j) = n

             l(n) = - 1.0d+20
             u(n) =   0.0d0
             x(n) = max( l(n), min( c(j), u(n) ) )
          end if
       end do

       if ( iprintctl(2) ) then
          write(* ,100) nslacks
          write(10,100) nslacks
       end if
    end if

    ! NON-EXECUTABLE STATEMENTS

100 format(/,1X,'Number of added slack variables   : ',I7)

200 format(/,1X,'All constraints are equality constraints. ', &
                'There are no slacks to be added.')

  end subroutine tinip

  ! ******************************************************************
  ! ******************************************************************

  subroutine tendp(n,x,l,u,m,lambda,equatn,linear,inform)

    use problvlv, only: reperr

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m
    integer, intent(inout) :: inform,n

    ! ARRAY ARGUMENTS
    logical,      intent(in)    :: linear(m)
    logical,      intent(inout) :: equatn(m)
    real(kind=8), intent(in)    :: lambda(m)
    real(kind=8), pointer       :: l(:),u(:),x(:)

    ! LOCAL SCALARS
    integer :: allocerr,j

    ! LOCAL ARRAYS
    character(len=15) :: subname
    real(kind=8) :: ltmp(nws),utmp(nws),xtmp(nws)

    if  ( slacks ) then

       ! Remove slacks

       xtmp(1:nws) = x(1:nws)
       ltmp(1:nws) = l(1:nws)
       utmp(1:nws) = u(1:nws)

       deallocate(x,l,u,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'TENDP'
          call reperr(inform,subname)
          return
       end if

       allocate(x(nws),l(nws),u(nws),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'TENDP'
          call reperr(inform,subname)
          return
       end if

       x(1:nws) = xtmp(1:nws)
       l(1:nws) = ltmp(1:nws)
       u(1:nws) = utmp(1:nws)

       n = nws

       do j = 1,m
          if ( slaind(j) .ne. - 1 ) then
             equatn(j) = .false.
          end if
       end do

       slacks = .false.

       ! Deallocate global arrays in module slacks

       deallocate(slaind,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'TENDP'
          call reperr(inform,subname)
          return
       end if

    end if

  end subroutine tendp

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalf(n,x,f,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    if ( .not. slacks ) then
       call uevalf(n,x,f,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalf(nws,x,f,inform)
       if ( inform .ne. 0 ) return
    end if

  end subroutine tevalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalg(n,x,g,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n)

    if ( .not. slacks ) then
       call uevalg(n,x,g,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalg(nws,x,g,inform)
       if ( inform .ne. 0 ) return

       g(nws+1:n) = 0.0d0
    end if

  end subroutine tevalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalh(n,x,hlin,hcol,hval,hnnz,lim,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: lim,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hcol(lim),hlin(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    if ( .not. slacks ) then
       call uevalh(n,x,hlin,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalh(nws,x,hlin,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return
    end if

  end subroutine tevalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalc(n,x,m,ind,c,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    integer :: j,sind

    if ( .not. slacks ) then
       call uevalc(n,x,m,ind,c,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalc(nws,x,m,ind,c,inform)
       if ( inform .ne. 0 ) return

       where ( ind(1:m) .and. slaind(1:m) .ne. -1 ) c(1:m) = c(1:m) - x(slaind(1:m))
    end if

  end subroutine tevalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,lim,inform)

    use problvlu
    use problvlv, only: reperr

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m)
    integer, intent(out) :: jclen(m),jcsta(m),jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    ! LOCAL SCALARS
    integer :: accum,j,nslacks,sind

    ! LOCAL ARRAYS
    character(len=15) :: subname
    integer :: origjclen(m),origjcsta(m),origjcvar(lim)
    real(kind=8) :: origjcval(lim)

    if ( .not. slacks ) then
       call uevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,lim,inform)
       if ( inform .ne. 0 ) return

    else
       call uevaljac(nws,x,m,ind,jcsta,jclen,jcvar,jcval,lim,inform)
       if ( inform .ne. 0 ) return

       ! Adding the derivatives with respect to the slacks variables
       ! is something that can be done 'in place' (without using
       ! additional temporary space), but coding it is very boring.
       
       nslacks = count( ind(1:m) .and. slaind(1:m) .ne. - 1 )

       if ( nslacks .gt. 0 ) then
          
          if ( sum( jclen(1:m), ind(1:m) ) + nslacks .gt. lim ) then
             inform = - 92
             subname = 'TEVALJAC'
             call reperr(inform,subname)
             return
          end if

          origjcsta(1:m) = jcsta(1:m)
          origjclen(1:m) = jclen(1:m)
          do j = 1,m
             origjcvar(jcsta(j):jcsta(j)+jclen(j)-1) = jcvar(jcsta(j):jcsta(j)+jclen(j)-1)
             origjcval(jcsta(j):jcsta(j)+jclen(j)-1) = jcval(jcsta(j):jcsta(j)+jclen(j)-1)
          end do
          
          accum = 0
          do j = 1,m
             if ( ind(j) ) then
                jcsta(j) = accum + 1
                jcvar(jcsta(j):jcsta(j)+origjclen(j)-1) = origjcvar(origjcsta(j):origjcsta(j)+origjclen(j)-1)
                jcval(jcsta(j):jcsta(j)+origjclen(j)-1) = origjcval(origjcsta(j):origjcsta(j)+origjclen(j)-1)
                jclen(j) = origjclen(j)

                sind = slaind(j)
                if ( sind .ne. - 1 ) then
                   jcvar(jcsta(j)+origjclen(j)) = sind
                   jcval(jcsta(j)+origjclen(j)) = - 1.0d0
                   jclen(j) = jclen(j) + 1
                end if

                accum = accum + jclen(j)
             end if
          end do
       end if
    end if

  end subroutine tevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalhc(n,x,ind,hlin,hcol,hval,hnnz,lim,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: ind,lim,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hcol(lim),hlin(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    if ( .not. slacks ) then
       call uevalhc(n,x,ind,hlin,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalhc(nws,x,ind,hlin,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return
    end if

  end subroutine tevalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,lim,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: lim,m,n
    integer,      intent(inout) :: inform
    integer,      intent(out)   :: hnnz
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hlin(lim),hcol(lim)
    real(kind=8), intent(in)  :: lambda(m),sc(m),x(n)
    real(kind=8), intent(out) :: hval(lim)

    if ( .not. slacks ) then
       call uevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalhl(nws,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return
    end if

  end subroutine tevalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(inout) :: gothl
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: lambda(m),p(n),sc(m),x(n)
    real(kind=8), intent(out) :: hp(n)

    if ( .not. slacks ) then
       call uevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalhlp(nws,x,m,lambda,sf,sc,p,hp,gothl,inform)
       if ( inform .ne. 0 ) return

       hp(nws+1:n) = 0.0d0
    end if

  end subroutine tevalhlp

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalfc(n,x,f,m,c,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    integer j,sind

    if ( .not. slacks ) then
       call uevalfc(n,x,f,m,c,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalfc(nws,x,f,m,c,inform)
       if ( inform .ne. 0 ) return

       do j = 1,m
          sind = slaind(j)
          if ( sind .ne. - 1 ) c(j) = c(j) - x(sind)
       end do
    end if

  end subroutine tevalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)

    use problvlu
    use problvlv, only: reperr

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
    integer :: j,sind

    ! LOCAL ARRAYS
    character(len=15) :: subname

    if ( .not. slacks ) then
       call uevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalgjac(nws,x(1:nws),g(1:nws),m,jcfun,jcvar,jcval,jcnnz,lim,inform)
       if ( inform .ne. 0 ) return

       g(nws+1:n) = 0.0d0

       do j = 1,m
          sind = slaind(j)
          if ( sind .ne. - 1 ) then
             if ( jcnnz + 1 .gt. lim ) then
                inform = - 92
                subname = 'TEVALGJAC'
                call reperr(inform,subname)
                return
             end if

             jcnnz = jcnnz +  1
             jcfun(jcnnz) = j
             jcvar(jcnnz) = sind
             jcval(jcnnz) = - 1.0d0
          end if
       end do
    end if

  end subroutine tevalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine tevalgjacp(n,x,g,m,p,q,work,gotj,inform)

    use problvlu

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

    ! LOCAL SCALARS
    integer :: j,sind

    if ( .not. slacks ) then
       call uevalgjacp(n,x,g,m,p,q,work,gotj,inform)
       if ( inform .ne. 0 ) return

    else
       call uevalgjacp(nws,x,g,m,p,q,work,gotj,inform)
       if ( inform .ne. 0 ) return

       if ( work .eq. 'J' .or. work .eq. 'T' ) then
          g(nws+1:n) = 0.0d0
       end if

       if ( work .eq. 'j' .or. work .eq. 'J' ) then
          do j = 1,m
             sind = slaind(j)
             if ( sind .ne. - 1 ) then
                p(j) = p(j) - q(sind)
             end if
          end do

       else ! if ( work .eq. 't' .or. work .eq. 'T' ) then
          q(nws+1:n) = 0.0d0

          do j = 1,m
             sind = slaind(j)
             if ( sind .ne. - 1 ) then
                q(sind) = q(sind) - p(j)
             end if
          end do
       end if
    end if

  end subroutine tevalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine tsetp(n,x,inform)

    use problvlu

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    if ( .not. slacks ) then
       call usetp(n,x,inform)
       if ( inform .ne. 0 ) return

    else
       call usetp(nws,x,inform)
       if ( inform .ne. 0 ) return
    end if

  end subroutine tsetp

  ! ******************************************************************
  ! ******************************************************************

  subroutine tunsetp()

    use problvlu

    implicit none

    call uunsetp()

  end subroutine tunsetp

  ! ******************************************************************
  ! ******************************************************************

end module problvlt
