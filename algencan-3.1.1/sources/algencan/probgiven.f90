module probgiven

  implicit none

  ! SCALARS

  logical, public :: ccoded,fccoded,fcoded,gcoded,gjaccoded, &
       gjacpcoded,jaccoded,jacpcoded,hcoded,hccoded,hlcoded, &
       hlpcoded,firstde,seconde,truehpr

  integer, public :: fcsubt,hnnzlim,jcnnzlim

  ! SUBROUTINES

  pointer :: evalf,evalg,evalh,evalc,evaljac,evalhc, &
       evalfc,evalgjac,evalgjacp,evalhl,evalhlp

  ! SUBROUTINES

  public :: iniprobgiven

  interface

     subroutine evalf(n,x,f,flag)
       implicit none
       ! SCALAR ARGUMENTS
       integer :: flag,n
       real(kind=8) :: f
       ! ARRAY ARGUMENTS
       real(kind=8) :: x(n)
     end subroutine evalf

     subroutine evalg(n,x,g,flag)
       implicit none
       ! SCALAR ARGUMENTS
       integer :: flag,n
       ! ARRAY ARGUMENTS
       real(kind=8) :: g(n),x(n)
     end subroutine evalg

     subroutine evalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)
       implicit none
       ! SCALAR ARGUMENTS
       logical :: lmem
       integer :: flag,lim,n,hnnz
       ! ARRAY ARGUMENTS
       integer :: hcol(*),hrow(*)
       real(kind=8) :: hval(*),x(n)
     end subroutine evalh

     subroutine evalc(n,x,ind,c,flag)
       implicit none
       ! SCALAR ARGUMENTS
       integer :: ind,flag,n
       real(kind=8) :: c
       ! ARRAY ARGUMENTS
       real(kind=8) :: x(n)
     end subroutine evalc

     subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)
       implicit none
       ! SCALAR ARGUMENTS
       logical :: lmem
       integer :: flag,ind,jcnnz,lim,n
       ! ARRAY ARGUMENTS
       integer :: jcvar(n)
       real(kind=8) :: x(n),jcval(n)
     end subroutine evaljac

     subroutine evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)
       implicit none
       ! SCALAR ARGUMENTS
       logical :: lmem
       integer :: flag,hcnnz,ind,lim,n
       ! ARRAY ARGUMENTS
       integer :: hccol(*),hcrow(*)
       real(kind=8) :: hcval(*),x(n)
     end subroutine evalhc

     subroutine evalfc(n,x,f,m,c,flag)
       implicit none
       ! SCALAR ARGUMENTS
       integer :: flag,m,n
       real(kind=8) :: f
       ! ARRAY ARGUMENTS
       real(kind=8) :: c(m),x(n)
     end subroutine evalfc

     subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)
       implicit none
       ! SCALAR ARGUMENTS
       logical :: lmem
       integer :: flag,jcnnz,lim,m,n
       ! ARRAY ARGUMENTS
       integer :: jcfun(*),jcvar(*)
       real(kind=8) :: g(n),jcval(*),x(n)
     end subroutine evalgjac

     subroutine evalgjacp(n,x,g,m,p,q,work,gotj,flag)
       implicit none
       ! SCALAR ARGUMENTS
       logical gotj
       integer :: flag,m,n
       character work
       ! ARRAY ARGUMENTS
       real(kind=8) :: g(n),p(m),q(n),x(n)
     end subroutine evalgjacp

     subroutine evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)
       implicit none
       ! SCALAR ARGUMENTS
       logical :: lmem
       integer :: flag,hlnnz,lim,m,n
       real(kind=8) :: sf
       ! ARRAY ARGUMENTS
       integer :: hlcol(*),hlrow(*)
       real(kind=8) :: hlval(*),lambda(m),sc(m),x(n)
     end subroutine evalhl

     subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)
       implicit none
       ! SCALAR ARGUMENTS
       logical goth
       integer :: flag,m,n
       real(kind=8) :: sf
       ! ARRAY ARGUMENTS
       real(kind=8) :: hp(n),lambda(m),p(n),sc(m),x(n)
     end subroutine evalhlp

  end interface

contains

  ! ******************************************************************
  ! ******************************************************************

  subroutine iniprobgiven(fsub,gsub,hsub,csub,jacsub,hcsub,fcsub, &
       gjacsub,gjacpsub,hlsub,hlpsub,coded,m,jcnnzmax,hnnzmax,inform)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: hnnzmax,jcnnzmax,m
    integer, intent(inout) :: inform

    ! ARRAY ARGUMENTS
    logical, intent(in) :: coded(11)

    ! EXTERNAL SUBROUTINES ARGUMENTS
    external :: fsub,gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub, &
         gjacpsub,hlsub,hlpsub

    ! LOCAL ARRAYS
    character(len=15) :: subname

    ! Set jcnnzlim and hnnzlim

    jcnnzlim = jcnnzmax
    hnnzlim = hnnzmax

    ! Set user-provided subroutines names
    
    evalf     => fsub
    evalg     => gsub
    evalh     => hsub
    evalc     => csub
    evaljac   => jacsub
    evalhc    => hcsub
    evalfc    => fcsub
    evalgjac  => gjacsub
    evalgjacp => gjacpsub
    evalhl    => hlsub
    evalhlp   => hlpsub

    ! Set user-provided subroutines indicators

    fcoded     = coded(1)
    gcoded     = coded(2)
    hcoded     = coded(3)
    ccoded     = coded(4)
    jaccoded   = coded(5)
    hccoded    = coded(6)
    fccoded    = coded(7)
    gjaccoded  = coded(8)
    gjacpcoded = coded(9)
    hlcoded    = coded(10)
    hlpcoded   = coded(11)
    
    ! Check whether mandatory subroutines are being properly provided.
    
    ! For unconstrained and bound-constrained problems, EVALF must be 
    ! coded by the user. For constrained problems, EVALF and EVALC, or, 
    ! alternatively, EVALFC must be coded. (Note that EVALF/EVALC and 
    ! EVALFC should not be provided concurrently.) For feasibility 
    ! problems, a constant null objective function must be coded and the 
    ! problem solved with the IGNORE-OBJECTIVE-FUNCTION keyword. Coded 
    ! subroutines must be indicated by setting the entrances of array 
    ! named CODED within subroutine INIP.
    
    ! Moreover, to avoid odd combinations, only the following choices
    ! will be considered valid:
    !
    ! If the objective function and the constraints are given by evalf 
    ! and evalc, respectively, first derivatives must be given by evalg 
    ! and evaljac, while second derivatives must be given by evalh and 
    ! evalhc.
    !
    ! If the objective function and the constraints are given by evalfc
    ! then first derivatives may be given by evalgjac or evalgjacp. If
    ! first derivatives are given by evalgjac, second derivatives must 
    ! be given by evalhl. On the other hand, if first derivatives are
    ! given by evalgjacp, second derivatives must be given by evalhlp.
    !
    ! Any other odd combination will be ignored.
    !
    ! Variables being set below have the following meaning:
    !
    ! firstde: whether, following the rules dictated above, the user-
    ! provided subroutines will allow the method to compute first 
    ! derivatives. 
    !
    ! seconde: whether, following the rules dictated above, the user-
    ! provided subroutines will allow the method to compute the Hessian 
    ! matrix of the augmented Lagrangian (to minimize the augmented 
    ! Lagrangian subproblems using a second-order method) and/or the
    ! Hessian of the Lagrangian plus the Jacobian of the constraints
    ! (in order to solve the KKT system by Newton's method). 
    !
    ! truehpr: whether, following the rules dictated above, the user-
    ! provided subroutines will allow the method to compute the 
    ! product of the true Hessian of the augmented Lagrangian times
    ! a given vector. It would allow the method to solve the augmented
    ! Lagrangian subproblems using a truncated-Newton method using
    ! Conjugate Gradients to solve the Newtonian linear systems.
    
    fcsubt = 0
 
    if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then
       fcsubt = 1
       firstde = .false.
       seconde = .false.
       truehpr = .false.
       if ( gcoded .and. ( jaccoded .or. m .eq. 0 ) ) then
          firstde = .true.
          if ( hcoded .and. ( hccoded .or. m .eq. 0 ) ) then
             seconde = .true.
             truehpr = .true.
          end if
       end if

    else if ( fccoded ) then
       fcsubt = 2
       firstde = .false.
       seconde = .false.
       truehpr = .false.
       if ( gjaccoded ) then
          firstde = .true.
          if ( hlcoded ) then
             seconde = .true.
             truehpr = .true.
          end if
       else 
          if ( gjacpcoded .and. hlpcoded ) then
             truehpr = .true.
          end if
       end if
    end if

    if ( fcsubt .eq. 0 ) then
       write(* ,1000)
       write(10,1000)
       stop
    end if

1000 format(/,1X,'*** Mandatory subroutines are not being ',    &
          'provided properly ***',/,1X,'For unconstrained and ',&
          'bound-constrained problems, EVALF must be coded by ',&
          'the',/,1X,'user. For constrained problems, EVALF ',  &
          'and EVALC, or, alternatively, EVALFC',/,1X,'must ',  &
          'be coded. (Note that EVALF/EVALC and EVALFC should ',&
          'not be provided',/,1X,'concurrently.) For ',         &
          'feasibility problems, a constant ',                  &
          'null objective function',/,1X,'must be coded and ',  &
          'the problem solved with the ',                       &
          'IGNORE-OBJECTIVE-FUNCTION',/,1X,'keyword. Coded ',   &
          'subroutines must be indicated by setting ',          &
          'the entrances of array',/,1X,'named CODED.')

  end subroutine iniprobgiven

end module probgiven
