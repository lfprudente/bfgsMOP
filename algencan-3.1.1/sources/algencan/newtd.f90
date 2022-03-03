! ******************************************************************
! ******************************************************************

subroutine newtd(n,nind,x,l,u,g,m,lambda,rho,equatn,d,adsupn, &
     maxelem,memfail,mindiag,inform)

  use lsslvr, only: lssana,lssend,lssfac,lssset,lsssol
  use probgiven, only: hnnzlim,jcnnzlim
  use modmachconst
  use modalgparam
  use modouttyp
  use moditetyp

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out)   :: memfail
  integer,      intent(in)    :: m,n,nind
  integer,      intent(inout) :: inform
  real(kind=8), intent(inout) :: adsupn,maxelem

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m)
  real(kind=8), intent(in)    :: g(nind),l(nind),lambda(m),rho(m), &
                                 u(nind),x(nind)
  real(kind=8), intent(inout) :: mindiag(nind)
  real(kind=8), intent(out)   :: d(nind)

  ! This subroutine solves the Newtonian system
  !
  !         ( H + rho A^T A ) x = b
  !
  ! by solving the Martinez-Santos system
  !
  !         H x + A^t y = b
  !         A x - y/rho = 0

  ! LOCAL SCALARS
  logical      :: adddiff
  integer      :: dim,unnz,i,iter,lssinfo,nneigv,nrank,pind
  real(kind=8) :: diff,pval

  ! LOCAL ARRAYS
  integer      :: udiag(nind+m),urow(hnnzlim+jcnnzlim+nind+m), &
                  ucol(hnnzlim+jcnnzlim+nind+m)
  real(kind=8) :: adddiag(nind+m),adddiagprev(nind+m), &
                  uval(hnnzlim+jcnnzlim+nind+m),sol(nind+m)

  ! ------------------------------------------------------------------
  ! Presentation
  ! ------------------------------------------------------------------

  if ( iprintinn .ge. 5 ) then
     write(* ,1000)
     write(10,1000)
  end if

  ! ------------------------------------------------------------------
  ! Initialization
  ! ------------------------------------------------------------------

  iter    = 0
  memfail = .false.

  call lssset(lsssubNW,sclsubNW,.true.,.false.)

  ! ------------------------------------------------------------------
  ! Compute ML matrix
  ! ------------------------------------------------------------------

  call mlsyst(nind,x,g,m,lambda,rho,equatn,urow,ucol,uval,unnz, &
       hnnzlim+jcnnzlim+m+nind,udiag,sol,dim,inform)

  if ( inform .ne. 0 ) return

  maxelem = 0.0d0
  do i = 1,unnz
     maxelem = max( maxelem, abs( uval(i) ) )
  end do

  ! ------------------------------------------------------------------
  ! Analyse sparsity pattern
  ! ------------------------------------------------------------------

  call lssana(dim,unnz,urow,ucol,uval,udiag,lssinfo)

  if ( lssinfo .eq. 6 ) then
     ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

     memfail = .true.
     return

  end if

  ! ------------------------------------------------------------------
  ! Main loop
  ! ------------------------------------------------------------------

100 continue

  iter = iter + 1

  ! ------------------------------------------------------------------
  ! Compute regularization
  ! ------------------------------------------------------------------

  if ( iter .eq. 1 ) then

     if ( sameface .and. ittype .eq. 3 ) then

        do i = 1,nind
           mindiag(i) = 0.1d0 * mindiag(i)
        end do

     else

        do i = 1,nind
           if ( g(i) .eq. 0.0d0 ) then
              mindiag(i) = 0.0d0
           else
              if ( g(i) .gt. 0.0d0 ) then
                 diff = x(i) - l(i)
              else if ( g(i) .lt. 0.0d0 ) then
                 diff = u(i) - x(i)
              end if
              mindiag(i) = abs( g(i) / diff )
           end if
        end do

     end if

     do i = 1,nind
        adddiag(i) = max( macheps23, mindiag(i) - uval(udiag(i)) )
     end do

  else

     do i = 1,nind
        adddiagprev(i) = adddiag(i)
     end do

110  continue

     do i = 1,nind
        if ( mindiag(i) .eq. 0.0d0 ) then
           mindiag(i) = macheps23
        else
           mindiag(i) = 10.0d0 * mindiag(i)
        end if
     end do

     do i = 1,nind
        adddiag(i) = max( macheps23, mindiag(i) - uval(udiag(i)) )
     end do

     adddiff = .false.
     do i = 1,nind
        if ( adddiag(i) .gt. adddiagprev(i) ) then
           adddiff = .true.
        end if
     end do

     if ( .not. adddiff ) then
        go to 110
     end if

  end if

  do i = nind + 1,dim
     adddiag(i) = 0.0d0
  end do

  adsupn = 0.0d0
  do i = 1,dim
     adsupn = max( adsupn, adddiag(i) )
  end do

  if ( iprintinn .ge. 5 ) then
     write(* ,1010) adsupn
     write(10,1010) adsupn
  end if

  ! ------------------------------------------------------------------
  ! Factorize matrix
  ! ------------------------------------------------------------------

  call lssfac(dim,unnz,urow,ucol,uval,udiag,adddiag,pind,pval,nneigv,nrank,lssinfo)

  if ( lssinfo .eq. 0 .or. lssinfo .eq. 1 ) then

     if ( nneigv .ne. dim - nind ) then
        ! WRONG INERTIA (SEE NOCEDAL AND WRIGHT)

        ! Lemma 16.3 [pg. 447]: Assume that the Jacobian of the
        ! constraints has full rank and that the reduced Hessian
        ! Z^T H Z is positive definite. Then the Jacobian of the
        ! KKT system has n positive eigenvalues, m negative
        ! eigenvalues, and no zero eigenvalues.

        ! Note that at this point we know that the matrix has no
        ! zero eigenvalues. nneigv gives the number of negative
        ! eigenvalues.

        if ( iprintinn .ge. 5 ) then
           write(* ,1020) nneigv,dim - nind
           write(10,1020) nneigv,dim - nind
        end if

        go to 100

     else

        if ( iprintinn .ge. 5 ) then
           write(* ,1030)
           write(10,1030)
        end if

     end if

  else if ( lssinfo .eq. 2 ) then
     ! SINGULAR JACOBIAN

     if ( iprintinn .ge. 5 ) then
        write(* ,1040)
        write(10,1040)
     end if

     go to 100

  else if ( lssinfo .eq. 6 ) then
     ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

     memfail = .true.
     return

  else if ( lssinfo .eq. 7 ) then
     ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

     memfail = .true.
     return

  else ! if ( lssinfo .eq. 8 ) then
     ! INSUFFICIENT INTEGER WORKING SPACE

     memfail = .true.
     return

  end if

  ! ------------------------------------------------------------------
  ! Solve
  ! ------------------------------------------------------------------

  call lsssol(dim,sol)

  call lssend()

  d(1:nind) = sol(1:nind)

  if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
     write(*, 1050) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
     write(10,1050) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
  end if

! NON-EXECUTABLE STATEMENTS

 1000 format(/,5X,'Sparse factorization of the ML system.')
 1010 format(  5X,'Maximum value added to the diagonal: ',1P,D24.16)
 1020 format(  5X,'ML-matrix with wrong inertia.', &
             /,5X,'Actual number of negative eigenvalues  = ',I16,'.', &
             /,5X,'Desired number of negative eigenvalues = ',I16,'.')
 1030 format(  5X,'Direct solver finished successfully.')
 1040 format(  5X,'ML-matrix numerically singular.')
 1050 format(/,5X,'Newton direction (first ',I7,' components): ', &
             /,1(5X,6(1X,1P,D11.4)))

end subroutine newtd

! ******************************************************************
! ******************************************************************

subroutine mlsyst(nind,x,nal,m,lambda,rho,equatn,urow,ucol,uval, &
     unnz,ulim,udiag,b,dim,inform)
  
  use modsvdgrad
  use modrspace
  use problvls, only: sevalhl

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,nind,ulim
  integer, intent(inout) :: inform
  integer, intent(out)   :: dim,unnz

  ! ARRAY ARGUMENTS
  logical,      intent(in)  :: equatn(m)
  integer,      intent(out) :: ucol(ulim),udiag(nind+m),urow(ulim)
  real(kind=8), intent(in)  :: lambda(m),nal(nind),rho(m),x(nind)
  real(kind=8), intent(out) :: b(nind+m),uval(ulim)

  ! LOCAL SCALARS
  integer :: col,i,j,k,row,var

  ! LOCAL ARRAYS
  integer :: wi(nfull)

  ! This subrotuine is called from the reduced space.

  ! MATRIX

  ! Set xfull with values in x
  xfull(ind(1:nind)) = x(1:nind)

  ! Compute Hessian of the Lagrangian
  call sevalhl(nfull,xfull,m,svddpdc,urow,ucol,uval,unnz,ulim,inform)
  if ( inform .ne. 0 ) return

  ! Shrink representation of Hessian of the Lagrangian and 
  ! set diagonal-elements indices

  wi(1:nfull) = 0
  wi(ind(1:nind)) = (/ (i, i=1,nind) /)

  udiag(1:nind) = 0

  k = 0
  do i = 1,unnz
     row = wi(urow(i))
     col = wi(ucol(i))

     if ( row .ne. 0 .and. col .ne. 0 ) then
        k = k + 1
        urow(k) = row
        ucol(k) = col
        uval(k) = uval(i)

        if ( row .eq. col ) then
           udiag(row) = k
        end if
     end if
  end do

  do i = 1,nind
     if ( udiag(i) .eq. 0 ) then
        k = k + 1
        urow(k) = i
        ucol(k) = i
        uval(k) = 0.0d0

        udiag(i) = k
     end if
  end do

  ! Shrink Jacobian and add diagonal matrix - 1.0 / rho

  dim = nind

  do j = 1,m
     if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then

        dim = dim + 1

        do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
           var = wi(svdjcvar(i))

           if ( var .ne. 0 ) then
              k = k + 1
              urow(k) = dim
              ucol(k) = var
              uval(k) = svdjcval(i)
           end if
        end do

        k = k + 1
        urow(k) = dim
        ucol(k) = dim
        uval(k) = - 1.0d0 / rho(j)

        udiag(dim) = k
     end if
  end do

  unnz = k

  ! RHS

  b(1:nind)     = - nal(1:nind)
  b(nind+1:dim) = 0.0d0

end subroutine mlsyst
