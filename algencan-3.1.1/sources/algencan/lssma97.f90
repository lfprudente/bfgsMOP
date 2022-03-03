! lssinfo:
!
! 0: Success.
! 1: Matrix not positive definite.
! 2: Rank deficient matrix.
! 6: Insufficient space to store the linear system.
! 7: Insufficient real working space.
! 8: Insufficient integer working space.

! ******************************************************************
! ******************************************************************

module lssma97

  use hsl_ma97_double
  
  implicit none

  save

  ! SCALARS
  character(len=4),   private :: lsclsub
  type(ma97_akeep),   private :: akeep
  type(ma97_fkeep),   private :: fkeep
  type(ma97_control), private :: cntl
  type(ma97_info),    private :: info

  ! ARRAYS
  integer,      allocatable, private :: rptr(:),rcol(:)
  real(kind=8), allocatable, private :: rval(:)

  !$omp threadprivate(lsclsub,akeep,fkeep,cntl,info,rptr,rcol,rval)

contains

  function lss_ma97()

    implicit none

    ! FUNCTION TYPE
    logical :: lss_ma97

    lss_ma97 = .true.

  end function lss_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssend_ma97()

    implicit none

    ! LOCAL SCALAR
    integer :: allocstat

    ! Finalize MA97 structures
    call ma97_finalise(akeep,fkeep)

  end subroutine lssend_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssini_ma97(sclsub,acneig,usefac)

    use modmachconst

    implicit none

    ! SCALAR ARGUMENTS
    character(len=4), intent(in) :: sclsub
    logical,          intent(in) :: acneig,usefac

    lsclsub = sclsub

    ! Supress monitoring, warning and error messages
    cntl%print_level = - 1

    if ( lsclsub .eq. 'MC64' ) then
       cntl%scaling = 1
    else if ( lsclsub .eq. 'MC77' ) then
       cntl%scaling = 2
    else if ( lsclsub .eq. 'MC30' ) then
       cntl%scaling = 4
    else ! if ( lsclsub .eq. 'NONE' ) then
       cntl%scaling = 0
    end if

  end subroutine lssini_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssana_ma97(nsys,hnnz,hrow,hcol,hval,hdiag,lssinfo)

    use modouttyp

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)    :: nsys,hnnz
    integer, intent(inout) :: lssinfo

    ! ARRAY ARGUMENTS
    integer,      intent(in) :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
    real(kind=8), intent(in) :: hval(hnnz)

    ! LOCAL SCALAR
    integer :: allocstat

    call ma97_analyse_coord(nsys,hnnz,hrow,hcol,akeep,cntl,info)
    
    ! SUCCESS
    if ( info%flag .ge. 0 ) then
       lssinfo = info%flag
    else
       ! UNHANDLED ERROR
       if ( iprintctl(3) ) then
          write(* ,9030) info%flag
          write(10,9030) info%flag
       end if
    end if

    ! NON-EXECUTABLE STATEMENTS

9000 format(/,1X,'LSSANA-MA97 WARNING: Insufficient space to store ', &
                 'linear system. Increase',                           &
            /,1X,'parameter nsysmax from ',I16,' to at least ',I16,   &
            /,1X,'if you would like to try a direct linear solver ',  &
                 'again.')

9030 format(/,1X,'LSSANA-MA97 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')

  end subroutine lssana_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssfac_ma97(nsys,hnnz,hrow,hcol,hval,hdiag,d,pind,pval, &
       nneigv,nrank,lssinfo)

    use modouttyp

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

    ! LOCAL SCALARS
    integer      :: i,idiag
    real(kind=8) :: d1

    ! LOCAL ARRAYS
    real(kind=8) :: w(nsys)

    w(1:nsys) = hval(hdiag(1:nsys))
    hval(hdiag(1:nsys)) = hval(hdiag(1:nsys)) + d(1:nsys)

    call ma97_factor(4,hval,akeep,fkeep,cntl,info)

    hval(hdiag(1:nsys)) = w(1:nsys)

    if ( info%flag .ge. 0 ) then
       ! SUCCESS
       lssinfo = 0

    else
       ! UNHANDLED ERROR
       if ( iprintctl(3) ) then
          write(* ,9000) info%flag
          write(10,9000) info%flag
       end if

       stop

    end if

    ! NUMBER OF NEGATIVE EIGENVALUES
    nneigv = info%num_neg
    nrank  = info%matrix_rank

    ! NON-EXECUTABLE STATEMENTS

9000 format(/,1X,'LSSFAC-MA97 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')

  end subroutine lssfac_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssol_ma97(nsys,sol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: sol(nsys)

    call ma97_solve(sol,akeep,fkeep,cntl,info)

  end subroutine lsssol_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssoltr_ma97(job,nsys,sol)

    implicit none

    ! SCALAR ARGUMENTS
    character(len=1), intent(in) :: job
    integer,          intent(in) :: nsys

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: sol(nsys)

  end subroutine lsssoltr_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsspermvec_ma97(nsys,v)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: v(nsys)

  end subroutine lsspermvec_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssunpermvec_ma97(nsys,v)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: v(nsys)

  end subroutine lssunpermvec_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsspermind_ma97(hnnz,hrow,hcol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: hnnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hrow(hnnz),hcol(hnnz)

  end subroutine lsspermind_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssunpermind_ma97(nsys,hnnz,hrow,hcol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys,hnnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hrow(hnnz),hcol(hnnz)

  end subroutine lssunpermind_ma97

  ! ******************************************************************
  ! ******************************************************************

  function lssgetd_ma97(j)

    implicit none

    ! FUNCTION TYPE
    real(kind=8) :: lssgetd_ma97

    ! SCALAR ARGUMENTS
    integer, intent(in) :: j

  end function lssgetd_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssetfactors_ma97(nsys,lssinfo)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys
    integer, intent(out) :: lssinfo

  end subroutine lsssetfactors_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssetrow_ma97(nsys)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

  end subroutine lsssetrow_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssgetrow_ma97(nsys,idx,rownnz,rowind,rowval)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys,idx
    integer, intent(out) :: rownnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: rowind(nsys)
    real(kind=8), intent(out) :: rowval(nsys)

  end subroutine lssgetrow_ma97

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssafsol_ma97(nsys,hnnz,hrow,hcol,hval,hdiag,d,sol,lssinfo)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys,hnnz
    integer, intent(out) :: lssinfo

    ! ARRAY ARGUMENTS
    integer,      intent(in)    :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
    real(kind=8), intent(in)    :: d(nsys),hval(hnnz)
    real(kind=8), intent(inout) :: sol(nsys)

  end subroutine lssafsol_ma97

end module lssma97
