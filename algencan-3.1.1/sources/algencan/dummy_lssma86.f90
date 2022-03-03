module lssma86

  implicit none

contains

  ! ******************************************************************
  ! ******************************************************************

  function lss_ma86()

    implicit none

    ! FUNCTION TYPE
    logical :: lss_ma86

    lss_ma86 = .false.

  end function lss_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssend_ma86()

    implicit none

  end subroutine lssend_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssini_ma86(sclsub,acneig,usefac)

    implicit none

    ! SCALAR ARGUMENTS
    character(len=4), intent(in) :: sclsub
    logical, intent(in) :: acneig,usefac

  end subroutine lssini_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssana_ma86(nsys,hnnz,hrow,hcol,hval,hdiag,lssinfo)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: nsys,hnnz
    integer, intent(inout) :: lssinfo

    ! ARRAY ARGUMENTS
    integer,      intent(in) :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
    real(kind=8), intent(in) :: hval(hnnz)

  end subroutine lssana_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssfac_ma86(nsys,hnnz,hrow,hcol,hval,hdiag,d,pind,pval, &
       nneigv,nrank,lssinfo)

    implicit none

    ! SCALAR ARGUMENTS
    integer,      intent(in)  :: hnnz,nsys
    integer,      intent(out) :: pind,nneigv,nrank,lssinfo
    real(kind=8), intent(out) :: pval

    ! ARRAY ARGUMENTS
    integer,      intent(in)    :: hdiag(nsys)
    integer,      intent(inout) :: hrow(hnnz),hcol(hnnz)
    real(kind=8), intent(in)    :: d(nsys)
    real(kind=8), intent(inout) :: hval(hnnz)

  end subroutine lssfac_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssol_ma86(nsys,sol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: sol(nsys)

  end subroutine lsssol_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssoltr_ma86(job,nsys,sol)

    implicit none

    ! SCALAR ARGUMENTS
    character(len=1), intent(in) :: job
    integer,          intent(in) :: nsys

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: sol(nsys)

  end subroutine lsssoltr_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsspermvec_ma86(nsys,v)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: v(nsys)

  end subroutine lsspermvec_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssunpermvec_ma86(nsys,v)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: v(nsys)

  end subroutine lssunpermvec_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsspermind_ma86(hnnz,hrow,hcol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: hnnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hrow(hnnz),hcol(hnnz)

  end subroutine lsspermind_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssunpermind_ma86(nsys,hnnz,hrow,hcol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys,hnnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hrow(hnnz),hcol(hnnz)

  end subroutine lssunpermind_ma86

  ! ******************************************************************
  ! ******************************************************************

  function lssgetd_ma86(j)

    implicit none

    ! FUNCTION TYPE
    real(kind=8) :: lssgetd_ma86

    ! SCALAR ARGUMENTS
    integer, intent(in) :: j

  end function lssgetd_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssetfactors_ma86(nsys,lssinfo)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys
    integer, intent(out) :: lssinfo

  end subroutine lsssetfactors_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssetrow_ma86(nsys)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

  end subroutine lsssetrow_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssgetrow_ma86(nsys,idx,rownnz,rowind,rowval)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys,idx
    integer, intent(out) :: rownnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: rowind(nsys)
    real(kind=8), intent(out) :: rowval(nsys)

  end subroutine lssgetrow_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssafsol_ma86(nsys,hnnz,hrow,hcol,hval,hdiag,d,sol,lssinfo)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys,hnnz
    integer, intent(out) :: lssinfo

    ! ARRAY ARGUMENTS
    integer,      intent(in)    :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
    real(kind=8), intent(in)    :: d(nsys),hval(hnnz)
    real(kind=8), intent(inout) :: sol(nsys)

  end subroutine lssafsol_ma86

end module lssma86
