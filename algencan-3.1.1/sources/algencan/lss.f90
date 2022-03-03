! lssinfo:
!
! 0: Success.
! 1: Matrix not positive definite.
! 2: Rank deficient matrix.
! 6: Insufficient space to store the linear system.
! 7: Insufficient double precision working space.
! 8: Insufficient integer :: working space.

! ******************************************************************
! ******************************************************************

module lsslvr

  use lssma57
  use lssma86
  use lssma97

  implicit none

  save

  abstract interface
     subroutine ilssafsol(nsys,hnnz,hrow,hcol,hval,hdiag,d,sol,lssinfo)
       ! SCALAR ARGUMENTS
       integer, intent(in)  :: nsys,hnnz
       integer, intent(out) :: lssinfo
       ! ARRAY ARGUMENTS
       integer,      intent(in)    :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
       real(kind=8), intent(in)    :: d(nsys),hval(hnnz)
       real(kind=8), intent(inout) :: sol(nsys)
     end subroutine ilssafsol

     subroutine ilssana(nsys,hnnz,hrow,hcol,hval,hdiag,lssinfo)
       ! SCALAR ARGUMENTS
       integer, intent(in)    :: nsys,hnnz
       integer, intent(inout) :: lssinfo
       ! ARRAY ARGUMENTS
       integer,      intent(in) :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
       real(kind=8), intent(in) :: hval(hnnz)
     end subroutine ilssana

     subroutine ilssend()
     end subroutine ilssend

     subroutine ilssfac(nsys,hnnz,hrow,hcol,hval,hdiag,d,pind,pval,nneigv,nrank,lssinfo)
       ! SCALAR ARGUMENTS
       integer,      intent(in)  :: hnnz,nsys
       integer,      intent(out) :: pind,nneigv,nrank,lssinfo
       real(kind=8), intent(out) :: pval
       ! ARRAY ARGUMENTS
       integer,      intent(in)    :: hdiag(nsys)
       integer,      intent(inout) :: hrow(hnnz),hcol(hnnz)
       real(kind=8), intent(in)    :: d(nsys)
       real(kind=8), intent(inout) :: hval(hnnz)
     end subroutine ilssfac

     function ilssgetd(j)
       ! FUNCTION TYPE
       real(kind=8)        :: ilssgetd
       ! SCALAR ARGUMENT
       integer, intent(in) :: j
     end function ilssgetd

     subroutine ilssgetrow(nsys,idx,rownnz,rowind,rowval)
       ! SCALAR ARGUMENTS
       integer, intent(in)  :: nsys,idx
       integer, intent(out) :: rownnz
       ! ARRAY ARGUMENTS
       integer,      intent(out) :: rowind(nsys)
       real(kind=8), intent(out) :: rowval(nsys)
     end subroutine ilssgetrow

     subroutine ilssini(sclsub,acneig,usefac)
       ! SCALAR ARGUMENTS
       character(len=4), intent(in) :: sclsub
       logical,          intent(in) :: acneig,usefac
     end subroutine ilssini

     subroutine ilsspermind(hnnz,hrow,hcol)
       ! SCALAR ARGUMENTS
       integer, intent(in) :: hnnz
       ! ARRAY ARGUMENTS
       integer, intent(inout) :: hrow(hnnz),hcol(hnnz)
     end subroutine ilsspermind

     subroutine ilsspermvec(nsys,v)
       ! SCALAR ARGUMENTS
       integer, intent(in) :: nsys
       ! ARRAY ARGUMENTS
       integer, intent(inout) :: v(nsys)
     end subroutine ilsspermvec

     subroutine ilsssetfactors(nsys,lssinfo)
       ! SCALAR ARGUMENTS
       integer, intent(in)  :: nsys
       integer, intent(out) :: lssinfo
     end subroutine ilsssetfactors

     subroutine ilsssetrow(nsys)
       ! SCALAR ARGUMENTS
       integer, intent(in) :: nsys
     end subroutine ilsssetrow

     subroutine ilsssol(nsys,sol)
       ! SCALAR ARGUMENTS
       integer, intent(in) :: nsys
       ! ARRAY ARGUMENTS
       real(kind=8), intent(inout) :: sol(nsys)
     end subroutine ilsssol

     subroutine ilsssoltr(job,nsys,sol)
       ! SCALAR ARGUMENTS
       character(len=1), intent(in) :: job
       integer,          intent(in) :: nsys
       ! ARRAY ARGUMENTS
       real(kind=8), intent(inout) :: sol(nsys)
     end subroutine ilsssoltr

     subroutine ilssunpermind(nsys,hnnz,hrow,hcol)
       ! SCALAR ARGUMENTS
       integer, intent(in) :: hnnz,nsys
       ! ARRAY ARGUMENTS
       integer, intent(inout) :: hrow(hnnz),hcol(hnnz)
     end subroutine ilssunpermind

     subroutine ilssunpermvec(nsys,v)
       ! SCALAR ARGUMENTS
       integer, intent(in) :: nsys
       ! ARRAY ARGUMENTS
       integer, intent(inout) :: v(nsys)
     end subroutine ilssunpermvec
  end interface

  ! The pointers to the routines
  procedure(ilssafsol),      pointer, protected :: lssafsol     
  procedure(ilssana),        pointer, protected :: lssana       
  procedure(ilssend),        pointer, protected :: lssend       
  procedure(ilssfac),        pointer, protected :: lssfac       
  procedure(ilssgetd),       pointer, protected :: lssgetd      
  procedure(ilssgetrow),     pointer, protected :: lssgetrow    
  procedure(ilssini),        pointer, protected :: lssini       
  procedure(ilsspermind),    pointer, protected :: lsspermind   
  procedure(ilsspermvec),    pointer, protected :: lsspermvec   
  procedure(ilsssetfactors), pointer, protected :: lsssetfactors
  procedure(ilsssetrow),     pointer, protected :: lsssetrow    
  procedure(ilsssol),        pointer, protected :: lsssol       
  procedure(ilsssoltr),      pointer, protected :: lsssoltr     
  procedure(ilssunpermind),  pointer, protected :: lssunpermind 
  procedure(ilssunpermvec),  pointer, protected :: lssunpermvec 

contains

  subroutine lssset(lsssub,sclsub,acneig,usefac)

    implicit none

    ! SCALAR ARGUMENTS
    character(len=4), intent(in) :: lsssub,sclsub
    logical,          intent(in) :: acneig,usefac

    ! ----------------------------------------------------------------
    ! This function sets the HSL subroutine to be used to solve linear
    ! systems. If the subroutine is available, lssr will hold its
    ! name; otherwise, the truncated Newton method will be used, and
    ! lssr will hold the emty string ''
    ! ----------------------------------------------------------------

    select case (lsssub)

    case ('MA57')
       ! Assign routine pointers
       lssafsol      => lssafsol_ma57
       lssana        => lssana_ma57
       lssend        => lssend_ma57
       lssfac        => lssfac_ma57
       lssgetd       => lssgetd_ma57
       lssgetrow     => lssgetrow_ma57
       lssini        => lssini_ma57
       lsspermind    => lsspermind_ma57
       lsspermvec    => lsspermvec_ma57
       lsssetfactors => lsssetfactors_ma57
       lsssetrow     => lsssetrow_ma57
       lsssol        => lsssol_ma57
       lsssoltr      => lsssoltr_ma57
       lssunpermind  => lssunpermind_ma57
       lssunpermvec  => lssunpermvec_ma57

    case ('MA86')
       ! Assign routine pointers
       lssafsol      => lssafsol_ma86
       lssana        => lssana_ma86
       lssend        => lssend_ma86
       lssfac        => lssfac_ma86
       lssgetd       => lssgetd_ma86
       lssgetrow     => lssgetrow_ma86
       lssini        => lssini_ma86
       lsspermind    => lsspermind_ma86
       lsspermvec    => lsspermvec_ma86
       lsssetfactors => lsssetfactors_ma86
       lsssetrow     => lsssetrow_ma86
       lsssol        => lsssol_ma86
       lsssoltr      => lsssoltr_ma86
       lssunpermind  => lssunpermind_ma86
       lssunpermvec  => lssunpermvec_ma86

    case ('MA97')
       ! Assign routine pointers
       lssafsol      => lssafsol_ma97
       lssana        => lssana_ma97
       lssend        => lssend_ma97
       lssfac        => lssfac_ma97
       lssgetd       => lssgetd_ma97
       lssgetrow     => lssgetrow_ma97
       lssini        => lssini_ma97
       lsspermind    => lsspermind_ma97
       lsspermvec    => lsspermvec_ma97
       lsssetfactors => lsssetfactors_ma97
       lsssetrow     => lsssetrow_ma97
       lsssol        => lsssol_ma97
       lsssoltr      => lsssoltr_ma97
       lssunpermind  => lssunpermind_ma97
       lssunpermvec  => lssunpermvec_ma97
    end select

    call lssini(sclsub,acneig,usefac)

  end subroutine lssset

end module lsslvr
