! SCALE OBJECTIVE FUNCTION AND CONSTRAINTS

module problvls

  implicit none

  ! SCALARS
  logical, protected :: scale
  real(kind=8), protected :: sf

  ! ARRAYS
  real(kind=8), allocatable, protected :: sc(:)

  ! SUBROUTINES
  public :: sinip,sendp,sevalf,sevalg,sevalh,sevalc,sevaljac, &
       sevalhc,sevalhl,sevalhlp,sevalfc,sevalgjac,sevalgjacp, &
       sevalobjc

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine setproblvls(val_scale)

    implicit none

    ! SCALARS
    logical, intent(in), optional :: val_scale

    scale = .true.
    if ( present( val_scale ) ) scale = val_scale

  end subroutine setproblvls

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine sinip(n,x,l,u,m,lambda,equatn,linear,coded,inform)

    use modouttyp, only: iprintctl
    use modmachconst, only: bignum,macheps12
    use problvlt, only: tsetp,tevalg,tevalgjac,tevalgjacp,tevaljac
    use problvlv, only: reperr
    use probgiven, only: fcsubt,gjaccoded,gjacpcoded,jcnnzlim

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical,      intent(in)    :: coded(11),equatn(m),linear(m)
    real(kind=8), intent(in)    :: l(n),u(n),x(n)
    real(kind=8), intent(inout) :: lambda(m)

    ! LOCAL SCALARS
    logical      :: gotj
    integer      :: allocerr,i,fun,j,jcnnz
    real(kind=8) :: scmin

    ! LOCAL ARRAYS
    character(len=15) :: subname
    logical :: ind(m)
    integer :: jcfun(jcnnzlim),jclen(m),jcsta(m),jcvar(jcnnzlim)
    real(kind=8) :: dum(1),g(n),jcval(jcnnzlim),p(m),q(n)

    if ( scale .and. m .eq. 0 ) then
       scale = .false.
       if ( iprintctl(2) ) then
          write(* ,200)
          write(10,200)
       end if
       return
    end if

    if ( scale ) then

       ! Allocate global arrays in module scaling

       allocate(sc(m),stat=allocerr)

       if ( allocerr .ne. 0 ) then
          inform = - 93
          subname = 'SINIP'
          call reperr(inform,subname)
          return
       end if

       ! Scaling

       call tsetp(n,x,inform)
       if ( inform .ne. 0 ) return

     ! if ( gjacpcoded ) then
       if ( fcsubt .eq. 2 .and. gjacpcoded ) then

          ! Scaling the constraints may be avoided if computing all
          ! individual gradients of constraints (even only once) is
          ! a very expensive task. In this case, the code below may
          ! be replaced by:
          !
          ! do j = 1,m
          ! sc(j) = 1.0d0
          ! end do
          !
          ! and a single call to tevalgjacp with p=0 and work equal
          ! to 'J' or 'T' may be preserved to compute the gradient 
          ! of the objective function g.

          ! Scale constraints

          p(1:m) = 0.0d0

          gotj = .false.

          do j = 1,m
             p(j) = 1.0d0
             if ( j .eq. 1 ) then
                call tevalgjacp(n,x,g,m,p,q,'T',gotj,inform)
                if ( inform .ne. 0 ) return
             else
                call tevalgjacp(n,x,dum,m,p,q,'t',gotj,inform)
                if ( inform .ne. 0 ) return
             end if
             p(j) = 0.0d0

             sc(j) = 1.0d0
             do i = 1,n
                sc(j) = max( sc(j), abs( q(i) ) )
             end do
          end do

          do j = 1,m
             sc(j) = 1.0d0 / sc(j)
          end do

     ! else if ( gjaccoded ) then
       else if ( fcsubt .eq. 2 ) then

          call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,jcnnzlim,inform)
          if ( inform .ne. 0 ) return

          ! Scale constraints

          do j = 1,m
             sc(j) = 1.0d0
          end do

          do i = 1,jcnnz
             fun = jcfun(i)
             sc(fun) = max( sc(fun), abs( jcval(i) ) )
          end do

          do j = 1,m
             sc(j) = 1.0d0 / sc(j)
          end do

       else

          call tevalg(n,x,g,inform)
          if ( inform .ne. 0 ) return

          ! Scale constraints

          ind(1:m) = .true.
          call tevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,jcnnzlim,inform)
          if ( inform .ne. 0 ) return

          do j = 1,m
             sc(j) = 1.0d0
             do i = jcsta(j),jcsta(j)+jclen(j)-1
                sc(j) = max( sc(j), abs( jcval(i) ) )
             end do
             sc(j) = 1.0d0 / sc(j)
          end do

       end if

       sc(1:m) = max( macheps12, sc(1:m) )

       ! Compute objective function scaling factor

       sf = 1.0d0
       do i = 1,n
          sf = max( sf, abs( g(i) ) )
       end do

       sf = max( macheps12, 1.0d0 / sf )

       ! Scale Lagrange multipliers initial guess

       lambda(1:m) = lambda(1:m) / sc(1:m) * sf

       ! Report scaling factors

       scmin = bignum
       do j = 1,m
          scmin = min( scmin, sc(j) )
       end do

       if ( iprintctl(2) ) then
          write(* ,300) sf,scmin
          write(10,300) sf,scmin
       end if

       ! Inhibit scaling feature if it will have no effect

       if ( sf .eq. 1.0d0 .and. ( m .eq. 0 .or. scmin .eq. 1.0d0 ) ) then
          scale = .false.

          if ( iprintctl(2) ) then
             write(* ,400)
             write(10,400)
          end if
       end if

    end if

    ! NON-EXECUTABLE STATEMENTS

200 format(/,1X,'The scaling feature was mainly developed for ',       &
                'constrained problems. For',/,1X,'unconstrained and ', &
                'bound-constrained problem, please, set the ',         &
                'optimality',/,1X,'tolerance (related to the ',        &
                'sup-norm of the projected gradient of the',/,1X,      &
                'objective function) with a convenient value. ',       &
                'The scaling feature is being',/,1X,'disabled.')

300 format(/,1X,'Objective function scale factor   : ',1P,D7.1, &
         /,1X,'Smallest constraints scale factor : ',1P,D7.1)
    
400 format(/,1X,'Since the scaling factor(s) is(are) equal to one, scaling was inhibited.')

  end subroutine sinip

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine sendp(n,x,l,u,m,lambda,equatn,linear,inform)

    use problvlv, only: reperr

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical,      intent(in)    :: equatn(m),linear(m)
    real(kind=8), intent(in)    :: l(n),u(n),x(n)
    real(kind=8), intent(inout) :: lambda(m)

    ! LOCAL SCALARS

    ! LOCAL ARRAYS
    integer :: deallocerr
    character(len=15) :: subname

    if ( scale ) then

       ! Undo scaling

       lambda(1:m) = lambda(1:m) * sc(1:m) / sf

    end if

    ! Deallocate global arrays in module scaling

    if ( allocated(sc) ) then
       
       deallocate(sc,stat=deallocerr)

       if ( deallocerr .ne. 0 ) then
          inform = - 94
          subname = 'SENDP'
          call reperr(inform,subname)
          return
       end if

    end if
    
  end subroutine sendp

  ! *****************************************************************
  ! *****************************************************************

  subroutine sevalobjc(n,x,f,fu,m,c,cu,inform,ignscl)

    use modalgparam, only: ignoref
    use probgiven,   only: ccoded,fcoded,fccoded
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(in), optional :: ignscl
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f,fu

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m),cu(m)

    ! LOCAL SCALARS
    logical :: lignscl
    integer :: j

    ! LOCAL ARRAYS
    logical :: ind(m)
    
    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

       if ( ignoref ) then
          fu = 0.0d0
       else
          call tevalf(n,x,fu,inform)
          if ( inform .ne. 0 ) return
       end if

       ind(1:m) = .true.
       call tevalc(n,x,m,ind,cu,inform)
       if ( inform .ne. 0 ) return

    else ! if ( fccoded ) then
       call tevalfc(n,x,fu,m,cu,inform)
       if ( inform .ne. 0 ) return

       if ( ignoref ) fu = 0.0d0
    end if

    if ( scale .and. .not. lignscl ) then
       f = fu * sf
       c(1:m) = sc(1:m) * cu(1:m)
    else
       f = fu
       c(1:m) = cu(1:m)
    end if

  end subroutine sevalobjc

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalf(n,x,f,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(in), optional :: ignscl
    integer,      intent(in)    :: n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    if ( ignoref ) then
       f = 0.0d0
       return
    end if

    call tevalf(n,x,f,inform)
    if ( inform .ne. 0 ) return

    if ( scale .and. .not. lignscl ) f = f * sf

  end subroutine sevalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalg(n,x,g,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in)    :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    if ( ignoref ) then
       g(1:n) = 0.0d0
       return
    end if

    call tevalg(n,x,g,inform)
    if ( inform .ne. 0 ) return

    if ( scale .and. .not. lignscl ) g(1:n) = sf * g(1:n)

  end subroutine sevalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalh(n,x,hlin,hcol,hval,hnnz,lim,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in)    :: lim,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hcol(lim),hlin(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    if ( ignoref ) then
       hnnz = 0
       return
    end if

    call tevalh(n,x,hlin,hcol,hval,hnnz,lim,inform)
    if ( inform .ne. 0 ) return

    if ( scale .and. .not. lignscl ) hval(1:hnnz) = sf * hval(1:hnnz)

  end subroutine sevalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalc(n,x,m,ind,c,inform,ignscl)

    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in) :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    call tevalc(n,x,m,ind,c,inform)
    if ( inform .ne. 0 ) return

    if ( scale .and. .not. lignscl ) then
       where ( ind(1:m) ) c(1:m) = c(1:m) * sc(1:m)
    end if

  end subroutine sevalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,lim,inform,ignscl)

    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in) :: lim,m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m)
    integer, intent(out) :: jclen(m),jcsta(m),jcvar(lim)
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: jcval(lim)

    ! LOCAL SCALARS
    logical :: lignscl
    integer :: j
    
    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    call tevaljac(n,x,m,ind,jcsta,jclen,jcvar,jcval,lim,inform)
    if ( inform .ne. 0 ) return

    if ( scale .and. .not. lignscl ) then
       do j = 1,m
          if ( ind(j) ) then
             jcval(jcsta(j):jcsta(j)+jclen(j)-1) = &
             jcval(jcsta(j):jcsta(j)+jclen(j)-1) * sc(j)
          end if
       end do
    end if

  end subroutine sevaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalhc(n,x,ind,hlin,hcol,hval,hnnz,lim,inform,ignscl)

    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in)    :: ind,lim,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hcol(lim),hlin(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    call tevalhc(n,x,ind,hlin,hcol,hval,hnnz,lim,inform)
    if ( inform .ne. 0 ) return

    if ( scale .and. .not. lignscl ) hval(1:hnnz) = sc(ind) * hval(1:hnnz)

  end subroutine sevalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalhl(n,x,m,lambda,hlin,hcol,hval,hnnz,lim,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in)    :: lim,m,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: hnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: hlin(lim),hcol(lim)
    real(kind=8), intent(in)  :: lambda(m),x(n)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    logical :: lignscl

    ! LOCAL ARRAYS
    real(kind=8) :: w(0:m)

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    if ( scale .and. .not. lignscl ) then
       w(0)   = sf
       w(1:m) = sc(1:m)
    else
       w(0:m) = 1.0d0
    end if

    if ( ignoref ) then
       w(0) = 0.0d0
    end if

    call tevalhl(n,x,m,lambda,w(0),w(1:m),hlin,hcol,hval,hnnz,lim,inform)
    if ( inform .ne. 0 ) return

  end subroutine sevalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalhlp(n,x,m,lambda,p,hp,gothl,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    logical, intent(inout) :: gothl
    integer, intent(in)    :: m,n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: lambda(m),p(n),x(n)
    real(kind=8), intent(out) :: hp(n)

    ! LOCAL SCALARS
    logical :: lignscl

    ! LOCAL ARRAYS
    real(kind=8) :: w(0:m)

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    if ( scale .and. .not. lignscl ) then
       w(0)   = sf
       w(1:m) = sc(1:m)
    else
       w(0:m) = 1.0d0
    end if

    if ( ignoref ) then
       w(0) = 0.0d0
    end if

    call tevalhlp(n,x,m,lambda,w(0),w(1:m),p,hp,gothl,inform)
    if ( inform .ne. 0 ) return

  end subroutine sevalhlp

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalfc(n,x,f,m,c,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical,      intent(in), optional :: ignscl
    integer,      intent(in)    :: m,n
    integer,      intent(inout) :: inform
    real(kind=8), intent(out)   :: f

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: c(m)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    call tevalfc(n,x,f,m,c,inform)
    if ( inform .ne. 0 ) return

    if ( ignoref ) f = 0.0d0

    if ( scale .and. .not. lignscl ) then
       f = f * sf
       c(1:m) = c(1:m) * sc(1:m)
    end if

  end subroutine sevalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in)    :: lim,m,n
    integer, intent(inout) :: inform
    integer, intent(out)   :: jcnnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=8), intent(in)  :: x(n)
    real(kind=8), intent(out) :: g(n),jcval(lim)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,inform)
    if ( inform .ne. 0 ) return

    if ( ignoref ) then
       g(1:n) = 0.0d0
    end if

    if ( scale .and. .not. lignscl ) then
       g(1:n) = sf * g(1:n)
       jcval(1:jcnnz) = jcval(1:jcnnz) * sc(jcfun(1:jcnnz))
    end if

  end subroutine sevalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine sevalgjacp(n,x,g,m,p,q,work,gotj,inform,ignscl)

    use modalgparam, only: ignoref
    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical,   intent(in), optional :: ignscl
    logical,   intent(inout) :: gotj
    integer,   intent(in)    :: m,n
    integer,   intent(inout) :: inform
    character, intent(in)    :: work

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in)    :: x(n)
    real(kind=8), intent(inout) :: p(m),q(n)
    real(kind=8), intent(out)   :: g(n)

    ! LOCAL SCALARS
    logical :: lignscl

    lignscl = .false.
    if ( present( ignscl) ) lignscl = ignscl

    if ( scale .and. .not. lignscl ) then
       if ( work .eq. 't' .or. work .eq. 'T' ) then
          p(1:m) = p(1:m) * sc(1:m)
       end if
    end if

    call tevalgjacp(n,x,g,m,p,q,work,gotj,inform)
    if ( inform .ne. 0 ) return

    if ( ignoref ) then
       if ( work .eq. 'J' .or. work .eq. 'T' ) then
          g(1:n) = 0.0d0
       end if
    end if

    if ( scale .and. .not. lignscl ) then
       if ( work .eq. 'j' .or. work .eq. 'J' ) then
          p(1:m) = p(1:m) * sc(1:m)
       end if

       if ( work .eq. 'J' .or. work .eq. 'T' ) then
          g(1:n) = sf * g(1:n)
       end if
    end if

  end subroutine sevalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine ssetp(n,x,inform,ignscl)

    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl
    integer, intent(in) :: n
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    call tsetp(n,x,inform)
    if ( inform .ne. 0 ) return

  end subroutine ssetp

  ! ******************************************************************
  ! ******************************************************************

  subroutine sunsetp(ignscl)

    use problvlt

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in), optional :: ignscl

    call tunsetp()

  end subroutine sunsetp

  ! ******************************************************************
  ! ******************************************************************

end module problvls
