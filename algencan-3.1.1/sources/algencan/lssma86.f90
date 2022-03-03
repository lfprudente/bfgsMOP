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

module lssma86

  use hsl_ma86_double
  use hsl_mc68_integer
  use hsl_mc69_double
  
  implicit none

  save

  ! SCALARS
  character(len=4),   private :: lsclsub
  integer,            private :: lmap
  type(mc68_control), private :: cntl68
  type(mc68_info),    private :: info68
  type(ma86_control), private :: cntl
  type(ma86_info),    private :: info
  type(ma86_keep),    private :: keep

  ! ARRAYS
  integer,      allocatable, private :: map(:),mrow(:),order(:),ptr(:)
  real(kind=8), allocatable, private :: mval(:)

!$omp threadprivate(lsclsub,lmap,cntl68,info68,cntl,info,keep,map, &
!$omp               mrow,order,ptr,mval)

contains

  function lss_ma86()

    implicit none

    ! FUNCTION TYPE
    logical :: lss_ma86

    lss_ma86 = .true.

  end function lss_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssend_ma86()

    implicit none

    ! LOCAL SCALAR
    integer :: allocstat

    ! deallocate(fact,ifact,invp,iwork,keep,posfac,s,sdiag,w,work)
    deallocate(map,mrow,mval,order,ptr)

    ! Finalize MA86 structures
    call ma86_finalise(keep,cntl)

  end subroutine lssend_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssini_ma86(sclsub,acneig,usefac)

    use modmachconst

    implicit none

    ! SCALAR ARGUMENTS
    character(len=4), intent(in) :: sclsub
    logical, intent(in) :: acneig,usefac

    lsclsub = sclsub

    ! Supress monitoring, warning and error messages
    cntl%diagnostics_level = - 1

    if ( lsclsub .eq. 'MC64' ) then
       cntl%scaling = 1
    else if ( lsclsub .eq. 'MC77' ) then
       cntl%scaling = 2
    else ! if ( lsclsub .eq. 'NONE' ) then
       cntl%scaling = 0
    end if

  end subroutine lssini_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssana_ma86(nsys,hnnz,hrow,hcol,hval,hdiag,lssinfo)

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

    ! Allocates necessary memory
    allocate(mval(hnnz),order(nsys),ptr(nsys+1),stat=allocstat)
    if(allocstat .ne. 0) then
       lssinfo = 6
       return
    end if

    ! Convert the coordinate sparse format to HSL sparse format
    call mc69_coord_convert(HSL_MATRIX_REAL_SYM_INDEF,nsys,nsys,hnnz,   &
         hrow,hcol,ptr,mrow,lssinfo,val_in=hval,val_out=mval,lmap=lmap, &
         map=map)

    if( lssinfo .ge. 0 ) then
       ! Computes elimination orderings
       call mc68_order(1,nsys,ptr,mrow,order,cntl68,info68)

       call ma86_analyse(nsys,ptr,mrow,order,keep,cntl,info)

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
    else
       ! UNHANDLED ERROR
       if ( iprintctl(3) ) then
          write(* ,9040) lssinfo
          write(10,9040) lssinfo
       end if
    end if

    ! NON-EXECUTABLE STATEMENTS

9000 format(/,1X,'LSSANA-MA86 WARNING: Insufficient space to store ', &
                 'linear system. Increase',                           &
            /,1X,'parameter nsysmax from ',I16,' to at least ',I16,   &
            /,1X,'if you would like to try a direct linear solver ',  &
                 'again.')

9030 format(/,1X,'LSSANA-MA86 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')

9040 format(/,1X,'LSSANA-MC69 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')

  end subroutine lssana_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssfac_ma86(nsys,hnnz,hrow,hcol,hval,hdiag,d,pind,pval, &
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

    call mc69_set_values(HSL_MATRIX_REAL_SYM_INDEF,lmap,map,hval, &
         ptr(nsys+1)-1,mval)

    hval(hdiag(1:nsys)) = w(1:nsys)

    call ma86_factor(nsys,ptr,mrow,mval,order,keep,cntl,info)

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

9000 format(/,1X,'LSSFAC-MA86 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')

  end subroutine lssfac_ma86

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssol_ma86(nsys,sol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: sol(nsys)

    call ma86_solve(sol,order,keep,cntl,info)

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
