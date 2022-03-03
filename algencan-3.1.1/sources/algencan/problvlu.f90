! DEAL WITH FIXED VARIABLES

module problvlu

  implicit none

  ! SCALARS
  logical, protected :: rmfixv

  ! ARRAYS
  integer, allocatable, private :: ycor(:), yind(:)
  real(kind=8), allocatable, private :: y(:)

  ! SUBROUTINES
  public :: setproblvlu,uevalf,uevalg,uevalh,uevalc,uevaljac, &
       uevalhc,uevalhl,uevalhlp,uevalfc,uevalgjac,uevalgjacp

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine setproblvlu(val_rmfixv)

    implicit none

    ! SCALARS
    logical, intent(in), optional :: val_rmfixv

    rmfixv = .true.
    if ( present( val_rmfixv ) ) rmfixv = val_rmfixv

  end subroutine setproblvlu

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine uinip(n,x,l,u,m,lambda,equatn,linear,coded,inform)

    use modprobdata, only: nbds
    use modouttyp, only: iprintctl
    use problvlv, only: reperr

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m
    integer, intent(inout) :: inform,n

    ! ARRAY ARGUMENTS
    logical,      intent(in) :: coded(11),equatn(m),linear(m)
    real(kind=8), intent(in) :: lambda(m)
    real(kind=8), pointer    :: l(:),u(:),x(:)

    ! LOCAL SCALARS
    integer :: allocerr,i
    real(kind=8) :: ltmp(n),utmp(n),xtmp(n)

    ! LOCAL ARRAYS
    character(len=15) :: subname

    if ( rmfixv .and. minval( u(1:n) - l(1:n) ) .gt. 0.0d0 ) then
       rmfixv = .false.
       if ( iprintctl(2) ) then
          write(* ,200)
          write(10,200)
       end if
       return
    end if

    if ( rmfixv ) then

       ! Allocate global arrays in module fixvar

       allocate(ycor(n),yind(0:n),y(n),stat=allocerr)

       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UINIP'
          call reperr(inform,subname)
          return
       end if

       ! Eliminate fixed variables (l=u) and save their values on y

       yind(0) = n

       n = 0
       do i = 1,yind(0)
          if ( l(i) .lt. u(i) ) then
             n = n + 1
             yind(n) = i
             ycor(i) = n
          else
             y(i) = l(i)
             ycor(i) = 0
          end if
       end do

       xtmp(1:n) = x(yind(1:n))
       ltmp(1:n) = l(yind(1:n))
       utmp(1:n) = u(yind(1:n))

       deallocate(x,l,u,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UINIP'
          call reperr(inform,subname)
          return
       end if

       allocate(x(n),l(n),u(n),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UINIP'
          call reperr(inform,subname)
          return
       end if

       x(1:n) = xtmp(1:n)
       l(1:n) = ltmp(1:n)
       u(1:n) = utmp(1:n)

       if ( iprintctl(2) ) then
          write(* ,100) yind(0) - n
          write(10,100) yind(0) - n
       end if

       nbds = nbds - 2 * ( yind(0) - n )

    end if

    ! Non-executable statements

100 format(/,1X,'Number of removed fixed variables : ',I7)

200 format(/,1X,'There are no fixed variables to be removed.')

  end subroutine uinip

  ! ******************************************************************
  ! ******************************************************************

  subroutine uendp(n,x,l,u,m,lambda,equatn,linear,inform)

    use problvlv, only: reperr

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m
    integer, intent(inout) :: inform,n

    ! ARRAY ARGUMENTS
    logical,      intent(in)    :: equatn(m),linear(m)
    real(kind=8), intent(in)    :: lambda(m)
    real(kind=8), pointer       :: l(:),u(:),x(:)

    ! LOCAL SCALARS
    integer :: allocerr,i,ind

    ! LOCAL ARRAYS
    character(len=15) :: subname
    real(kind=8) :: ltmp(n),utmp(n),xtmp(n)

    if ( rmfixv ) then

       ! Restore original x, l, u and n

       xtmp(1:n) = x(1:n)
       ltmp(1:n) = l(1:n)
       utmp(1:n) = u(1:n)

       deallocate(x,l,u,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UENDP'
          call reperr(inform,subname)
          return
       end if

       allocate(x(yind(0)),l(yind(0)),u(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UENDP'
          call reperr(inform,subname)
          return
       end if

       do i = yind(0),1,-1
          ind = ycor(i)
          if ( ind .ne. 0 ) then
             l(i) = ltmp(ind)
             u(i) = utmp(ind)
             x(i) = xtmp(ind)
          else
             l(i) = y(i)
             u(i) = y(i)
             x(i) = y(i)
          end if
       end do

       n = yind(0)

       rmfixv = .false.

       ! Deallocate global arrays in module fixvar

       deallocate(ycor,yind,y,stat=allocerr)

       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UENDP'
          call reperr(inform,subname)
          return
       end if

    end if

  end subroutine uendp

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalf(n,x,f,inform)

    use problvlv

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: allocerr

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       call vevalf(n,x,f,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALF'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalf(yind(0),ytmp,f,inform)
       if ( inform .ne. 0 ) return

       deallocate(ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALF'
          call reperr(inform,subname)
          return
       end if
    end if

  end subroutine uevalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalg(n,x,g,inform)

    use problvlv

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    integer :: allocerr

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:),gtmp(:)

    if ( .not. rmfixv ) then
       call vevalg(n,x,g,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(ytmp(yind(0)),gtmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALG'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalg(yind(0),ytmp,gtmp,inform)
       if ( inform .ne. 0 ) return

       g(1:n) = gtmp(yind(1:n))

       deallocate(ytmp,gtmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALG'
          call reperr(inform,subname)
          return
       end if
    end if

  end subroutine uevalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalh(n,x,hrow,hcol,hval,hnnz,lim,inform)

    use problvlv

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
    integer :: allocerr,col,i,j,row

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       call vevalh(n,x,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALH'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalh(yind(0),ytmp,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

       deallocate(ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALH'
          call reperr(inform,subname)
          return
       end if

       j = 0
       do i = 1,hnnz
          row = ycor(hrow(i))
          col = ycor(hcol(i))
          if ( row .ne. 0 .and. col .ne. 0 ) then
             j = j + 1
             hrow(j) = row
             hcol(j) = col
             hval(j) = hval(i)
          end if
       end do

       hnnz = j
    end if

  end subroutine uevalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalc(n,x,m,ind,c,inform)

    use problvlv

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    integer :: allocerr,j

    ! LOCAL ARRAYS
    character(len=15) :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       do j = 1,m
          if ( ind(j) ) then
             call vevalc(n,x,j,c(j),inform)
             if ( inform .ne. 0 ) return
          end if
       end do
       
    else
       allocate(ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALC'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       do j = 1,m
          if ( ind(j) ) then
             call vevalc(yind(0),ytmp,j,c(j),inform)
             if ( inform .ne. 0 ) return
          end if
       end do
       
       deallocate(ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALC'
          call reperr(inform,subname)
          return
       end if
    end if

  end subroutine uevalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,lim,inform)

    use problvlv

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
    integer :: allocerr,i,j,jcelemj,k,var

    ! LOCAL ARRAYS
    character(len=15) :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       k = 0
       do j = 1,m
          if ( ind(j) ) then
             jcsta(j) = k + 1
             call vevaljac(n,x,j,jcvar(jcsta(j)),jcval(jcsta(j)),jclen(j),lim-k,inform)
             if ( inform .ne. 0 ) return
             k = k + jclen(j)
          end if
       end do
       
    else
       allocate(ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALJAC'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       k = 0
       do j = 1,m
          if ( ind(j) ) then
             jcsta(j) = k + 1
             call vevaljac(yind(0),ytmp,j,jcvar(jcsta(j)),jcval(jcsta(j)),jclen(j),lim-k,inform)
             if ( inform .ne. 0 ) return

             jcelemj = 0
             do i = jcsta(j),jcsta(j)+jclen(j)-1
                var = ycor(jcvar(i))
                if ( var .ne. 0 ) then
                   jcvar(jcsta(j)+jcelemj) = var
                   jcval(jcsta(j)+jcelemj) = jcval(i)
                   jcelemj = jcelemj + 1
                end if
             end do
             jclen(j) = jcelemj
             k = k + jclen(j)
          end if
       end do

       deallocate(ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALJAC'
          call reperr(inform,subname)
          return
       end if
    end if

  end subroutine uevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalhc(n,x,ind,hrow,hcol,hval,hnnz,lim,inform)

    use problvlv

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
    integer :: allocerr,col,i,j,row

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       call vevalhc(n,x,ind,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALHC'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalhc(yind(0),ytmp,ind,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

       deallocate(ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALHC'
          call reperr(inform,subname)
          return
       end if

       j = 0
       do i = 1,hnnz
          row = ycor(hrow(i))
          col = ycor(hcol(i))
          if ( row .ne. 0 .and. col .ne. 0 ) then
             j = j + 1
             hrow(j) = row
             hcol(j) = col
             hval(j) = hval(i)
          end if
       end do

       hnnz = j
    end if

  end subroutine uevalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalhl(n,x,m,lambda,sf,sc,hrow,hcol,hval,hnnz,lim,inform)

    use problvlv

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
    integer :: allocerr,col,i,j,row

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       call vevalhl(n,x,m,lambda,sf,sc,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALHL'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalhl(yind(0),ytmp,m,lambda,sf,sc,hrow,hcol,hval,hnnz,lim,inform)
       if ( inform .ne. 0 ) return

       deallocate(ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALHL'
          call reperr(inform,subname)
          return
       end if

       j = 0
       do i = 1,hnnz
          row = ycor(hrow(i))
          col = ycor(hcol(i))
          if ( row .ne. 0 .and. col .ne. 0 ) then
             j = j + 1
             hrow(j) = row
             hcol(j) = col
             hval(j) = hval(i)
          end if
       end do

       hnnz = j
    end if

  end subroutine uevalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

    use problvlv

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
    integer :: allocerr

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: hptmp(:),w(:),ytmp(:)

    if ( .not. rmfixv ) then
       call vevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(hptmp(yind(0)),w(yind(0)),ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALHLP'
          call reperr(inform,subname)
          return
       end if

       w(1:yind(0)) = 0.0d0
       w(yind(1:n)) = p(1:n)

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalhlp(yind(0),ytmp,m,lambda,sf,sc,w,hptmp,gothl,inform)
       if ( inform .ne. 0 ) return

       hp(1:n) = hptmp(yind(1:n))

       deallocate(hptmp,w,ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALHLP'
          call reperr(inform,subname)
          return
       end if

    end if

  end subroutine uevalhlp

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalfc(n,x,f,m,c,inform)

    use problvlv

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    integer :: allocerr

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       call vevalfc(n,x,f,m,c,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALFC'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalfc(yind(0),ytmp,f,m,c,inform)
       if ( inform .ne. 0 ) return

       deallocate(ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALFC'
          call reperr(inform,subname)
          return
       end if
    end if

  end subroutine uevalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)

    use problvlv

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
    integer :: i,j,allocerr,var

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:),gtmp(:)

    if ( .not. rmfixv ) then
       call vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(ytmp(yind(0)),gtmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALGJAC'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       call vevalgjac(yind(0),ytmp,gtmp,m,jcfun,jcvar,jcval,jcnnz,lim,inform)
       if ( inform .ne. 0 ) return

       g(1:n) = gtmp(yind(1:n))

       j = 0
       do i = 1,jcnnz
          var = ycor(jcvar(i))
          if ( var .ne. 0 ) then
             j = j + 1
             jcfun(j) = jcfun(i)
             jcvar(j) = var
             jcval(j) = jcval(i)
          end if
       end do

       jcnnz = j

       deallocate(ytmp,gtmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALGJAC'
          call reperr(inform,subname)
          return
       end if
    end if

  end subroutine uevalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine uevalgjacp(n,x,g,m,p,q,work,gotj,inform)

    use problvlv

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
    integer :: allocerr

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: gtmp(:),qtmp(:),w(:),ytmp(:)

    if ( .not. rmfixv ) then
       call vevalgjacp(n,x,g,m,p,q,work,gotj,inform)
       if ( inform .ne. 0 ) return

    else
       allocate(gtmp(yind(0)),qtmp(yind(0)),w(yind(0)),ytmp(yind(0)),stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALGJACP'
          call reperr(inform,subname)
          return
       end if

       ytmp(1:yind(0)) = y(1:yind(0))
       ytmp(yind(1:n)) = x(1:n)

       if ( work .eq. 'j' .or. work .eq. 'J' ) then

          w(1:yind(0)) = 0.0d0
          w(yind(1:n)) = q(1:n)

          call vevalgjacp(yind(0),ytmp,gtmp,m,p,w,work,gotj,inform)
          if ( inform .ne. 0 ) return

       else ! if ( work .eq. 't' .or. work .eq. 'T' ) then
          call vevalgjacp(yind(0),ytmp,gtmp,m,p,qtmp,work,gotj,inform)
          if ( inform .ne. 0 ) return

          q(1:n) = qtmp(yind(1:n))
       end if

       if ( work .eq. 'J' .or. work .eq. 'T' ) then
          g(1:n) = gtmp(yind(1:n))
       end if

       deallocate(gtmp,qtmp,w,ytmp,stat=allocerr)
       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'UEVALGJACP'
          call reperr(inform,subname)
          return
       end if

    end if

  end subroutine uevalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine usetp(n,x,inform)

    use problvlv

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer,   intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: allocerr,i

    ! LOCAL ARRAYS
    character(len=15)     :: subname
    real(kind=8), pointer :: ytmp(:)

    if ( .not. rmfixv ) then
       call vsetp(n,x,inform)
       return
    end if

    allocate(ytmp(yind(0)),stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 93
       subname = 'USETP'
       call reperr(inform,subname)
       return
    end if

    ytmp(1:yind(0)) = y(1:yind(0))
    ytmp(yind(1:n)) = x(1:n)

    call vsetp(yind(0),ytmp,inform)
    if ( inform .ne. 0 ) return

    deallocate(ytmp,stat=allocerr)
    if ( allocerr .ne. 0 ) then
       inform = - 93
       subname = 'USETP'
       call reperr(inform,subname)
       return
    end if

  end subroutine usetp

  ! ******************************************************************
  ! ******************************************************************

  subroutine uunsetp()

    use problvlv

    implicit none

    call vunsetp()

  end subroutine uunsetp

  ! *****************************************************************
  ! *****************************************************************

end module problvlu
