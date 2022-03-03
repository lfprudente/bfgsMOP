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

module lssma57

  use hsl_ma57_double

  implicit none

  save

  ! MA57 STRUCTURES
  type(zd11_type), private :: msys
  type(ma57_control), private :: cntl
  type(ma57_factors), private :: fact
  type(ma57_ainfo), private :: ainfo
  type(ma57_finfo), private :: finfo
  type(ma57_sinfo), private :: sinfo

  ! SCALARS
  character(len=4), private :: lsclsub
  integer, private :: dnnz,rnnz
  logical, private :: lacneig,lusefac

  ! ARRAYS
  integer,  allocatable, private :: dptr(:),drow(:),invperm(:), &
                                    perm(:),rptr(:),rcol(:)
  real(kind=8), allocatable, private :: dval(:),rval(:),s(:)

  !$omp threadprivate(msys,cntl,fact,ainfo,finfo,sinfo,lsclsub, &
  !$omp               dnnz,rnnz,lacneig,lusefac,dptr,drow,invperm, &
  !$omp               perm,rptr,rcol,dval,rval,s)

contains

  ! ******************************************************************
  ! ******************************************************************

  function lss_ma57()

    implicit none

    ! FUNCTION TYPE
    logical :: lss_ma57

    lss_ma57 = .true.

  end function lss_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssend_ma57()

    implicit none

    ! LOCAL SCALAR
    integer :: info

    ! Finalize MA57 structures
    call ma57_finalize(fact,cntl,info)

    ! Deallocate memory
    deallocate(msys%col,msys%row,msys%val,stat=info)
    deallocate(dptr,drow,dval,perm,invperm,s,rptr,rcol,rval,stat=info)

  end subroutine lssend_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssini_ma57(sclsub,acneig,usefac)

    use modmachconst, only: macheps23

    implicit none

    ! SCALAR ARGUMENTS
    character(len=4), intent(in) :: sclsub
    logical,          intent(in) :: acneig,usefac

    lacneig = acneig
    lusefac = usefac
    lsclsub = sclsub

    call ma57_initialize(control=cntl)

    ! Suppress monitoring, warning and error messages
    cntl%ldiag = 0

    if ( lsclsub .ne. 'NONE' ) then
       cntl%scaling = 1
    else
       cntl%scaling = 0
    end if

    if ( .not. lacneig ) then
       cntl%pivoting  = 2
       cntl%tolerance = macheps23
    end if

  end subroutine lssini_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssana_ma57(nsys,hnnz,hrow,hcol,hval,hdiag,lssinfo)

    use modouttyp

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys,hnnz
    integer, intent(inout) :: lssinfo

    ! ARRAY ARGUMENTS
    integer, intent(in) :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
    real(kind=8), intent(in) :: hval(hnnz)

    ! LOCAL SCALAR
    integer :: allocstat,i

    ! FUNCTIONS
    integer, external :: OMP_GET_THREAD_NUM

    ! Allocates and set HSL matrix representation
    allocate(msys%col(hnnz),msys%row(hnnz),msys%val(hnnz),stat=allocstat)
    if ( allocstat .ne. 0 ) then
       lssinfo = 6
       return
    end if

    ! Copies the matrix of the system to the HSL structure
    msys%n  = nsys
    msys%ne = hnnz

    msys%row(1:hnnz) = hrow(1:hnnz)
    msys%col(1:hnnz) = hcol(1:hnnz)
    msys%val(1:hnnz) = hval(1:hnnz)

!!$    write(*,*) 'lssana_ma57 says: '
!!$    write(*,*) 'msys%n : ',msys%n
!!$    write(*,*) 'msys%ne: ',msys%ne
!!$    write(*,*) 'maxval : ',maxval( abs( msys%val(1:hnnz) ) )

    ! Initialize factors structure
    call ma57_initialize(factors=fact)

    call ma57_analyse(msys,fact,cntl,ainfo)

!!$    write(*,*) 'ainfo%flag: ',ainfo%flag
!!$    write(*,*) 'ainfo%oor : ',ainfo%oor
!!$    write(*,*) 'ainfo%dup : ',ainfo%dup

    if ( ainfo%flag .ge. 0 ) then
       ! SUCCESS
       lssinfo = 0

       ! This test is to check if lssana is being called from moresor
       ! subroutine. In this case, we need to compute the permutation
       ! and inverse permutation array.
       if( .not. lacneig ) then
          if( .not. allocated(perm) .or. size(perm) .ne. nsys ) then
             deallocate(perm,invperm,stat=allocstat)
             
             allocate(perm(nsys),invperm(nsys),stat=allocstat)
             if(allocstat .ne. 0) then
                lssinfo = 6
                return
             end if
          end if

          call ma57_enquire(fact,perm=perm)
          invperm(perm(1:nsys)) = (/ (i, i = 1, nsys) /)
       end if

       return
    end if

    ! UNHANDLED ERROR
    if ( iprintctl(3) ) then
       write(* ,9000) ainfo%flag
       write(10,9000) ainfo%flag
    end if

    ! NON-EXECUTABLE STATEMENTS

9000 format(/,1X,'LSSANA-MA57 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')

  end subroutine lssana_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssfac_ma57(nsys,hnnz,hrow,hcol,hval,hdiag,d,pind,pval, &
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
    integer      :: allocstat,i
    real(kind=8) :: d1,maxmax,maxmax2

!!$    write(*,*) 'lssfac_ma57 says: '
!!$    write(*,*) 'nnz: ',hnnz
!!$    write(*,*) 'maxval: ',maxval( abs( hval(1:hnnz) ) )

    msys%val(hdiag(1:nsys)) = hval(hdiag(1:nsys)) + d(1:nsys)

!!$    write(*,*) 'after adding diagonal: '
!!$    write(*,*) 'maxval: ',maxval( abs( msys%val(1:hnnz) ) )

!!$    write(*,*) 'matrix: '
!!$    write(*,*) 'msys%n : ',msys%n
!!$    write(*,*) 'msys%ne: ',msys%ne
!!$    do i = 1,msys%ne
!!$       write(*,*) msys%row(i),msys%col(i),msys%val(i)
!!$    end do

!!$    write(*,*) 'lssana_ma57 says: '
!!$    write(*,*) 'maxval : ',maxval( abs( msys%val(1:hnnz) ) )

    call ma57_factorize(msys,fact,cntl,finfo)

!!$    allocate(s(nsys),stat=allocstat)
!!$    call ma57_enquire(fact,scaling=s)
!!$    write(*,*) 'scaling: '
!!$    do i = 1,nsys
!!$       write(*,*) i,s(i),maxval( abs( msys%val(1:hnnz) ), msys%row(1:hnnz) .eq. i )
!!$    end do

!!$    maxmax = 0.0d0
!!$    maxmax2 = 0.0d0
!!$    do i = 1,msys%ne
!!$       maxmax  = max( maxmax,  abs( msys%val(i) ) )
!!$       maxmax2 = max( maxmax2, abs( msys%val(i) ) * s(msys%row(i)) * s(msys%col(i)) )
!!$    end do
!!$    write(*,*) 'maxmax: ',maxmax,' maxmax2: ',maxmax2

    ! This test is to check if lssfac is being called from moresor
    ! subroutine. In this case, we need to compute the permutation
    ! and inverse permutation array.
    if( .not. lacneig .and. lsclsub .ne. 'NONE' ) then
       if( .not. allocated(s) .or. size(s) .ne. nsys ) then
          deallocate(s,stat=allocstat)

          allocate(s(nsys),stat=allocstat)
          if(allocstat .ne. 0) then
             lssinfo = 6
             return
          end if

          call ma57_enquire(fact,scaling=s)
       end if
    end if

    if ( finfo%flag .eq. 0 ) then

       d1 = hval(hdiag(1)) + d(1)

       if ( d1 .gt. 0.0d0 ) then
          ! SUCCESS
          lssinfo = 0
       else
          ! MATRIX IS NEGATIVE DEFINITE

          pind    = 1
          pval    = abs( d1 )

          lssinfo = 1
       end if

    else if ( finfo%flag .eq. -6 ) then
       ! MATRIX NOT POSITIVE DEFINITE

       pind    = finfo%more
       pval    = abs( finfo%pivot )

       lssinfo = 1

    else if ( finfo%flag .eq. 4 .or. finfo%flag .eq. -5 ) then
       ! RANK DEFICIENT MATRIX

       pind    = finfo%more
       pval    = abs( finfo%pivot )

       lssinfo = 2

    else if ( finfo%flag .eq. -3 ) then
       ! ALLOCATION MEMORY ERROR

       if ( iprintctl(3) ) then
          write(* ,9010)
          write(10,9010)
       end if

       lssinfo = 6

    else
       ! UNHANDLED ERROR

       if ( iprintctl(3) ) then
          write(* ,9000) finfo%flag
          write(10,9000) finfo%flag
       end if

       stop

    end if

    ! NUMBER OF NEGATIVE EIGENVALUES
    nneigv = finfo%neig
    nrank  = finfo%rank

!!$    write(*,*) 'lssfac_ma57 says: '
!!$    write(*,*) 'finfo % flag: ',finfo%flag
!!$    write(*,*) 'finfo % rank: ',finfo%rank
!!$    write(*,*) 'finfo % neig: ',finfo%neig
!!$    write(*,*) 'lssinfo     : ',lssinfo

    ! NON-EXECUTABLE STATEMENTS

9000 format(/,1X,'LSSFAC-MA57 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')
9010 format(/,1X,'LSSFAC-MA57 ERROR: Unable to allocate memory.')

  end subroutine lssfac_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssol_ma57(nsys,sol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: sol(nsys)

    call ma57_solve(msys,fact,sol,cntl,sinfo)

  end subroutine lsssol_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssoltr_ma57(job,nsys,sol)

    implicit none

    ! SCALAR ARGUMENTS
    character(len=1), intent(in) :: job
    integer,          intent(in) :: nsys

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: sol(nsys)

    ! LOCAL SCALARS
    integer :: info

    ! LOCAL ARRAYS
    real(kind=8) :: work(nsys)

    if ( job .eq. 'T' .or. job .eq. 't' ) then
       
       call ma57_part_solve(fact,cntl,'L',sol,info)

       work(1:nsys) = sol(invperm(1:nsys)) * sqrt( 1.0d0 / dval(1:nsys) )
       sol(1:nsys)  = work(1:nsys)

    else

       work(1:nsys) = sol(invperm(1:nsys)) * sqrt( 1.0d0 / dval(perm(1:nsys)) )
       sol(1:nsys)  = work(1:nsys)

       call ma57_part_solve(fact,cntl,'U',sol,info)

    end if

  end subroutine lsssoltr_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsspermvec_ma57(nsys,v)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: v(nsys)

    ! LOCAL ARRAYS
    integer :: iwork(nsys)

    iwork(1:nsys) = v(invperm(1:nsys))
    v(1:nsys) = iwork(1:nsys)

  end subroutine lsspermvec_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssunpermvec_ma57(nsys,v)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: v(nsys)
    
    ! LOCAL ARRAYS
    integer :: iwork(nsys)

    iwork(1:nsys) = v(perm(1:nsys))
    v(1:nsys) = iwork(1:nsys)

  end subroutine lssunpermvec_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsspermind_ma57(hnnz,hrow,hcol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: hnnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hrow(hnnz),hcol(hnnz)

    hcol(1:hnnz) = perm(hcol(1:hnnz))
    hrow(1:hnnz) = perm(hrow(1:hnnz))

  end subroutine lsspermind_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssunpermind_ma57(nsys,hnnz,hrow,hcol)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: hnnz,nsys

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hrow(hnnz),hcol(hnnz)

    hcol(1:hnnz) = invperm(hcol(1:hnnz))
    hrow(1:hnnz) = invperm(hrow(1:hnnz))

  end subroutine lssunpermind_ma57

  ! ******************************************************************
  ! ******************************************************************

  function lssgetd_ma57(j)

    implicit none

    ! FUNCTION TYPE
    real(kind=8) :: lssgetd_ma57

    ! SCALAR ARGUMENTS
    integer, intent(in) :: j

    lssgetd_ma57 = sqrt( dval(j) )

    return

  end function lssgetd_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssetfactors_ma57(nsys,lssinfo)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys
    integer, intent(out) :: lssinfo

    ! LOCAL SCALARS
    integer :: allocstat,i,rownnz

    ! TEMP
    real(kind=8) :: loc_sdiag(nsys)

    dnnz = max(2*finfo%ntwo + nsys, nsys)
    rnnz = finfo%nebdu

    ! Allocates memory
    if( .not. allocated(rptr) .or. size(rptr) .lt. nsys+1 .or. &
         size(rval) .lt. rnnz ) then
       deallocate(rptr,rcol,rval,stat=allocstat)

       allocate(rptr(nsys+1),rcol(rnnz),rval(rnnz),stat=allocstat)
       if(allocstat .ne. 0) then
          lssinfo = 6
          return
       end if
    end if

    if( .not. allocated(dptr) .or. size(dptr) .lt. nsys+1 .or. &
         size(dval) .lt. dnnz ) then
       deallocate(dptr,drow,dval,stat=allocstat)

       allocate(dptr(nsys+1),drow(dnnz),dval(dnnz),stat=allocstat)
       if(allocstat .ne. 0) then
          lssinfo = 6
          return
       end if
    end if

    if( .not. allocated(perm) .or. size(perm) .lt. nsys ) then
       deallocate(perm,invperm,stat=allocstat)
       
       allocate(perm(nsys),invperm(nsys),stat=allocstat)
       if(allocstat .ne. 0) then
          lssinfo = 6
          return
       end if
    end if

    if( .not. allocated(s) .or. size(s) .lt. nsys ) then
       deallocate(s,stat=allocstat)

       allocate(s(nsys),stat=allocstat)
       if(allocstat .ne. 0) then
          lssinfo = 6
          return
       end if
    end if

    ! Get factors D, S, and L
    call ma57_get_factors(fact,cntl,rnnz,rptr,rcol,rval, &
         dnnz,dptr,drow,dval,invperm,perm,s,sinfo)

    ! Unscale factors, if necessary
    if( lsclsub .ne. 'NONE' ) then
       rval(1:rnnz) = rval(1:rnnz) / s(rcol(1:rnnz))
    end if

    ! loc_sdiag(1:nsys) = 1.0d0 / dval(perm(1:nsys))

    ! Computes R = D^0.5 L^T
    do i = 1,nsys
       ! rval(rptr(i):rptr(i+1)-1) = rval(rptr(i):rptr(i+1)-1) * sqrt( dval(i) )
       rval(rptr(i):rptr(i+1)-1) = rval(rptr(i):rptr(i+1)-1) / sqrt( 1.0d0 / dval(i) )
       ! rval(rptr(i):rptr(i+1)-1) = rval(rptr(i):rptr(i+1)-1) / sqrt( loc_sdiag(invperm(i)) )
    end do

    lssinfo = 0

  end subroutine lsssetfactors_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lsssetrow_ma57(nsys)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nsys

  end subroutine lsssetrow_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssgetrow_ma57(nsys,idx,rownnz,rowind,rowval)

    use modmachconst, only: macheps

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys,idx
    integer, intent(out) :: rownnz

    ! ARRAY ARGUMENTS
    integer,      intent(out) :: rowind(nsys)
    real(kind=8), intent(out) :: rowval(nsys)

    ! LOCAL SCALARS
    integer :: i

    rownnz = 0
    do i = rptr(idx), rptr(idx+1)-1
       if( abs(rval(i)) .gt. macheps ) then
          rownnz = rownnz + 1
          rowind(rownnz) = rcol(i)
          rowval(rownnz) = rval(i)
       end if
    end do
    
  end subroutine lssgetrow_ma57

  ! ******************************************************************
  ! ******************************************************************

  subroutine lssafsol_ma57(nsys,hnnz,hrow,hcol,hval,hdiag,d,sol,lssinfo)

    use modmachconst, only: macheps23
    use modouttyp

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in)  :: nsys,hnnz
    integer, intent(out) :: lssinfo

    ! ARRAY ARGUMENTS
    integer,      intent(in)    :: hrow(hnnz),hcol(hnnz),hdiag(nsys)
    real(kind=8), intent(in)    :: d(nsys),hval(hnnz)
    real(kind=8), intent(inout) :: sol(nsys)

    ! LOCAL STRUCTURES
    type(zd11_type)    :: lmsys
    type(ma57_control) :: lcntl
    type(ma57_factors) :: lfact
    type(ma57_ainfo)   :: lainfo
    type(ma57_finfo)   :: lfinfo
    type(ma57_sinfo)   :: lsinfo

    ! LOCAL SCALARS
    integer :: allocstat,info
    
    allocate(lmsys%col(hnnz),lmsys%row(hnnz),lmsys%val(hnnz), &
         stat=allocstat)
    if( allocstat .ne. 0 ) then
       lssinfo = 6
       return
    end if

    ! Copies the matrix of the system to be solved
    lmsys%n  = nsys
    lmsys%ne = hnnz

    lmsys%row(1:hnnz) = hrow(1:hnnz)
    lmsys%col(1:hnnz) = hcol(1:hnnz)
    lmsys%val(1:hnnz) = hval(1:hnnz)

    ! Initialize MA57 structures
    call ma57_initialize(factors=lfact,control=lcntl)

    ! Suppress monitoring, warning and error messages
    lcntl%ldiag = 0

    ! Supress scaling
    lcntl%scaling = 0

    if ( .not. lacneig ) then
       lcntl%pivoting  = 2
       lcntl%tolerance = macheps23
    end if

    ! Analysis
    call ma57_analyse(lmsys,lfact,lcntl,lainfo)

    ! Analysis error handling
    if( lainfo%flag .ge. 0 ) then
       ! SUCESS
       lssinfo = 0
    else
       ! UNHANDLED ERROR
       if( iprintctl(3) ) then
          write( *,9000) lainfo%flag
          write(10,9000) lainfo%flag
       end if

       stop
    end if

    lmsys%val(hdiag(1:nsys)) = lmsys%val(hdiag(1:nsys)) + d(1:nsys)

    if( lsclsub .ne. 'NONE') then
       lmsys%val(1:hnnz) = lmsys%val(1:hnnz) * &
            s(lmsys%row(1:hnnz)) * &
            s(lmsys%col(1:hnnz))
    end if

    ! Factorization
    call ma57_factorize(lmsys,lfact,lcntl,lfinfo)

    ! Factorization error handling
    if ( lfinfo%flag .eq. 0 .or. lfinfo%flag .eq. 1 ) then
       ! SUCESS
       lssinfo = 0

    else if ( lfinfo%flag .eq. 5 .or. lfinfo%flag .eq. -6 ) then
       ! MATRIX NOT POSITIVE DEFINITE
       lssinfo = 1

    else if ( lfinfo%flag .eq. 4 .or. lfinfo%flag .eq. -5 ) then
       ! RANK DEFICIENT MATRIX
       lssinfo = 2

    else
       ! UNHANDLED ERROR
       if ( iprintctl(3) ) then
          write(* ,9000) lfinfo%flag
          write(10,9000) lfinfo%flag
       end if

       stop
    end if

    ! Scaling
    if ( lsclsub .ne. 'NONE' ) then
       sol(1:nsys) = sol(1:nsys) * s(1:nsys)
    end if

    ! Solve
    call ma57_solve(lmsys,lfact,sol,lcntl,lsinfo)

    ! Scaling
    if ( lsclsub .ne. 'NONE' ) then
       sol(1:nsys) = sol(1:nsys) * s(1:nsys)
    end if

    ! Finalize structures
    call ma57_finalize(lfact,lcntl,info)

    deallocate(lmsys%col,lmsys%row,lmsys%val)

    ! NON-EXECUTABLE STATEMENTS

9000 format(/,1X,'LSSAFSOL-MA57 ERROR: Unhandled error ',I16,'.', &
            /,1X,'See documentation for details.')

  end subroutine lssafsol_ma57

end module lssma57
