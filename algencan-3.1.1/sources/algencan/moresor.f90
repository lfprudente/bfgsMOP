! ******************************************************************
! ******************************************************************

subroutine moresor(n,g,bnnz,blin,bcol,bval,bdiag,delta,sigma1, &
     sigma2,eps,maxit,l,pd,p,chcnt,memfail,msinfo)
  
  use lsslvr
  use modmachconst
  use modalgparam, only: lsssubTR,sclsubTR
  use modouttyp

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(inout) :: memfail,pd
  integer,      intent(in)    :: maxit,n
  integer,      intent(inout) :: bnnz,chcnt,msinfo
  real(kind=8), intent(in)    :: delta,eps,sigma1,sigma2
  real(kind=8), intent(inout) :: l

  ! ARRAY ARGUMENTS
  integer,      intent(inout) :: bcol(bnnz),bdiag(n),blin(bnnz)
  real(kind=8), intent(in)    :: g(n)
  real(kind=8), intent(inout) :: bval(bnnz),p(n)

  ! Solves the problem
  !
  ! minimize    psi(w) = 1/2 w^TBw + g^Tw
  ! subject to  ||w|| <= delta
  !
  ! Using the method described in "Computing a trust region step",
  ! by More and Sorensen.

  ! msinfo:
  !
  ! 0: both g and B are null
  ! 1: first convergence criterion is satisfied
  ! 2: second convergence criterion is satisfied
  ! 3: third convergence criterion is satisfied
  ! 5: maximum allowed number of iterations is achieved

  ! LOCAL SCALARS
  integer      :: col,dum,i,idx,iter,lin,lssinfo
  real(kind=8) :: b1n,d,delta2,geucn,ll,llant,ls,lsant,lu,luant,   &
                  peucn,peucn2,ptz,rpeucn2,rzeucn2,sgn,tau,teucn2, &
                  tmp,ueucn2

  ! LOCAL ARRAYS
  integer      :: wi(n)
  real(kind=8) :: t(n),wd(n),z(n)

  delta2 = delta**2

  msinfo  = 0

  memfail = .false.

  ! print *,'CALL LSSSET'
  call lssset(lsssubTR,'NONE',.false.,.true.)

  ! step 1: initialize ls (lower bound on l) with max{-bii}, where bii
  ! are the elements of the diagonal of B

  ls = -bval(bdiag(1))

  do i = 2,n
     lin = bdiag(i)
     ls  = max( ls, -bval(lin) )
  end do

  ! Calculate ||B||1, B sparse

  do i = 1,n
     wd(i) = 0.0d0
  end do

  do i = 1,bnnz
     lin = blin(i)
     col = bcol(i)

     if ( ( lin .le. n ) .and. ( col .le. n ) ) then
        wd(col) = wd(col) + abs( bval(i) )
        if ( lin .ne. col ) then
           wd(lin) = wd(lin) + abs( bval(i) )
        end if
     end if
  end do

  b1n = wd(1)
  do i = 2,n
     b1n = max( b1n, wd(i) )
  end do

  ! step 2: initialize ll (lower bound on l) with
  !         max{0, ls, ||g||/delta - ||B||1}, where ||B||1 is the
  !         1-norm of the matrix B

  geucn = 0.0d0
  do i = 1,n
     geucn = geucn + g(i)**2
  end do
  geucn = sqrt( geucn )

  ll = (geucn / delta) - b1n
  ll = max( 0.0d0, ll )
  ll = max( ls, ll )

  ! step 3: initialize lu (upper bound on l) with ||g||/delta + ||B||1

  lu = (geucn / delta) + b1n

  ! If the matrix is null, there is nothing to be done

  if ( ( abs( ll ) .le. macheps23 ) .and. &
       ( abs( lu ) .le. macheps23 ) .and. &
       ( abs( ls ) .le. macheps23 ) ) then
     msinfo = 0
     go to 21
  end if

  ! step 4: initialize iteration counter

  iter = 1

  ! print *,'CALL LSSANA'
  call lssana(n,bnnz,blin,bcol,bval,bdiag,lssinfo)

  if ( lssinfo .eq. 6 ) then
     ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

     memfail = .true.
     return

  else if ( lssinfo .eq. 7 ) then
     ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

     memfail = .true.
     return

  else if ( lssinfo .eq. 8 ) then
     ! INSUFFICIENT INTEGER WORKING SPACE

     memfail = .true.
     return

  end if

  ! step 5: safeguard of l (ensures that l is bigger than ll)

5 continue

  l = max( l, ll )

  ! step 6: safeguard of l (ensures that l is smaller than lu)

  l = min( l, lu )

  ! step 7: safeguard of l

  if ( ( l .le. ls + macheps23 * max( abs( ls ), 1.0d0 ) ) &
       .and. ( iter .ne. 1 ) ) then
     l = max( 1.0d-3 * lu, sqrt( ll*lu ) )
  end if

  ! step 8: try to use the Cholesky decomposition: (B +lI) = R^TR.
  !         If the decomposition is successfull, R is stored in the
  !         upper triangular portion of B (including the diagonal) and
  !         pd is set to true.
  !         If the decomposition fails, d and idx are set as explained
  !         before, pd is set to false, and the Euclidian-norm of u is
  !         calculated (see explanation of variable ueucn2)

  do i = 1,n
     wd(i) = l
  end do

  call lssfac(n,bnnz,blin,bcol,bval,bdiag,wd,idx,d,dum,dum,lssinfo)
  chcnt = chcnt + 1

  if ( lssinfo .eq. 6 ) then
     ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

     memfail = .true.
     return

  else if ( lssinfo .eq. 7 ) then
     ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

     memfail = .true.
     return

  else if ( lssinfo .eq. 8 ) then
     ! INSUFFICIENT INTEGER WORKING SPACE

     memfail = .true.
     return

  end if

  if ( ( lssinfo .eq. 1 ) .or. ( lssinfo .eq. 2 ) ) then
     pd = .false.
  else ! lssinfo .eq. 0
     pd = .true.
     ! print *,'CALL LSSSETFACTORS'
     call lsssetfactors(n,lssinfo)
  end if

  ! In this case (B + lI) is not positive definite, and d and idx are
  ! calculated. Because p cannot be calculated (it is not possible to
  ! solve the system using Cholesky factorization), the values of l,
  ! ll and ls are updated, the iteration counter is increased and a
  ! new iteration is started

  if ( .not. pd ) then

     !     Print information (current iteration)

     if ( iprintinn .ge. 5 ) then
        write(*, 1000) iter
        write(*, 1010) ls,ll,lu,l
        write(*, 1070)

        write(10,1000) iter
        write(10,1010) ls,ll,lu,l
        write(10,1070)
     end if

     llant = ll
     luant = lu
     lsant = ls

     if ( lu-l .le. macheps23 * max( abs( lu ),1.0d0 ) ) then
        lu = lu + macheps23 * max( abs( lu ),1.0d0 )
        l  = lu
     else
        call scalcu(n,bnnz,blin,bcol,bval,bdiag,l,idx,p,ueucn2,wd, &
             memfail)
        chcnt = chcnt + 1

        if ( memfail ) then
           go to 22
        end if

        ll = max( l, ll )
        ls = max( l + (d / ueucn2), ls )
        ll = max( ll, ls )
        l  = ls
     end if

     iter = iter + 1

     !     Test whether the number of iterations is exhausted

     if ( iter .gt. maxit ) then
        msinfo = 5

        if ( iprintinn .ge. 5 ) then
           write(*, 1090)
           write(10,1090)
        end if

        go to 22
     end if

     go to 5
  end if

  ! step 9: solve R^TRp = -g for p and calculate the squared
  !         Euclidian-norm of p

  do i = 1,n
     p(i) = -g(i)
  end do

  ! print *,'CALL LSSSOL'
  call lsssol(n,p)

  ! Euclidian-norm of Rp = p^T R^TRp = p^T (-g) = - p^Tg

  rpeucn2 = 0.0d0
  do i = 1,n
     rpeucn2 = rpeucn2 - p(i) * g(i)
  end do

  peucn2 = 0.0d0
  do i = 1,n
     peucn2 = peucn2 + p(i)**2
  end do
  peucn = sqrt( peucn2 )

  ! step 10: calculate z and tau, where tau * z is the approximation
  !          of the eigenvector associated with the smallest
  !          eigenvalue of B

  if ( peucn .lt. delta ) then

     !    Calculate z

     call scalcz(n,wi,wd,z)

     !     Calculate z Euclidian-norm

     tmp = 0.0d0
     do i = 1,n
        tmp = tmp + z(i)**2
     end do
     tmp = sqrt( tmp )

     !     Divide z by its norm

     do i = 1,n
        z(i) = z(i) / tmp
     end do

     !     Calculate the squared Euclidian-norm of the product Rz.
     !     Note that z^T R^T Rz = z^T (B + lI) z

     do i = 1,n
        wd(i) = 0.0d0
     end do

     do i = 1,bnnz
        lin = blin(i)
        col = bcol(i)

        if ( lin .eq. col ) then
           wd(lin) = wd(lin) + (bval(i) + l) * z(col)
        else
           wd(lin) = wd(lin) + bval(i) * z(col)
           wd(col) = wd(col) + bval(i) * z(lin)
        end if
     end do

     rzeucn2 = 0.0d0
     do i = 1,n
        rzeucn2 = rzeucn2 + z(i) * wd(i)
     end do

     !     Calculate tau

     ptz = 0.0d0
     do i = 1,n
        ptz = ptz + p(i) * z(i)
     end do

     if ( ptz .lt. 0.0d0 ) then
        sgn = -1.0d0
     else
        sgn =  1.0d0
     end if

     tmp = delta2 - peucn2

     tau = (ptz**2) + tmp
     tau = tmp / (ptz + sgn * sqrt( tau ))

  end if

  ! Print informations (current iteration)

  if ( iprintinn .ge. 5 ) then
     write(*, 1000) iter
     write(*, 1010) ls,ll,lu,l
     write(*, 1020) peucn

     write(10,1000) iter
     write(10,1010) ls,ll,lu,l
     write(10,1020) peucn

     if ( peucn .lt. delta ) then
        write(*, 1030) sqrt( peucn2 + 2 * tau * ptz + tau ** 2 )
        write(10,1030) sqrt( peucn2 + 2 * tau * ptz + tau ** 2 )
     else
        write(*, 1030) peucn
        write(10,1030) peucn
     end if

     write(*, 1050) delta
     write(10,1050) delta
  end if

  ! steps 11 and 12: update ll, lu and ls

  llant = ll
  luant = lu
  lsant = ls

  if ( peucn .lt. delta ) then
     lu = min( l, lu )
     ls = max( l - rzeucn2, ls )
  else
     ll = max( l, ll )
  end if

  ! step 13: update ls when B + lI is not positive definite.
  !          This was done right after the Cholesky decomposition
  !          failure

  ! step 14: update ll

  ll = max( ll, ls )

  ! step 15: convergence test

  if ( ( abs( l ) .le. eps ) .and. ( peucn .le. delta ) ) then
     msinfo = 1
     go to 21
  end if

  ! step 16: second convergence test

  if ( abs( delta - peucn ) .le. sigma1*delta ) then
     msinfo = 2
  end if

  ! step 17: convergence test for the hard case

  tmp = rpeucn2 + l * delta2
  tmp = max( sigma2, tmp )
  tmp = tmp * sigma1 * (2 - sigma1)

  if ( ( peucn .lt. delta ) .and. ( (rzeucn2*(tau**2) .le. tmp ) &
       .or. (lu - ll .le. eps) ) ) then

     msinfo = msinfo + 3
     go to 20
  end if

  if ( msinfo .eq. 2 ) then
     go to 21
  end if

  ! step 21: Calculate l to be used in the next iteration

  if ( ( abs( geucn ) .gt. eps ) .or. &
       ( ( l .le. ll + macheps23 * max( abs( ll ), 1.0d0 ) ) &
       .and. ( ls .le. ll ) ) ) then

     !     Solve R^T t = p for t and calculate t squared Euclidian-norm

     do i = 1,n
        t(i) = p(i)
     end do

     ! print *,'CALL LSSSOLTR'
     call lsssoltr('T',n,t)

     teucn2 = 0.0d0
     do i = 1,n
        teucn2 = teucn2 + t(i)**2
     end do

     !     Update l using Newton's method update

     l = l + (peucn2 / teucn2) * ((peucn - delta) / delta)
  else
     l = ls
  end if

  ! step 22: update iteration counter

  iter = iter + 1

  ! Test whether the number of iterations is exhausted

  if ( iter .gt. maxit ) then
     msinfo = 5

     if ( iprintinn .ge. 5 ) then
        write(*, 1090)
        write(10,1090)
     end if

     go to 22
  end if

  ! step 23: start a new iteration

  if ( ( abs( llant-ll ) .le.                      &
       macheps23 * max( abs( ll ), 1.0d0 ) ) .and. &
       ( abs( luant-lu ) .le.                      &
       macheps23 * max( abs( lu ), 1.0d0 ) ) .and. &
       ( abs( lsant-ls ) .le.                      &
       macheps23 * max( abs( ls ), 1.0d0 ) ) ) then

     ll = ll + macheps23 * max( abs( ll ), 1.0d0 )
     lu = lu - macheps23 * max( abs( lu ), 1.0d0 )
  end if

  go to 5

  ! steps 18, 19 and 20:

  ! The solution is given by p in 3 cases:
  ! - if the first convergence criterion is satisfied;
  ! - if only the second convergence criterion is satisfied;
  ! - if both the second and the third convergence criteria are
  !   satisfied, but the squared Euclidian-norm of R(tau*z) is
  !   strictly bigger than l*(delta2 - peucn2)

  ! The solution is given by p + tau*z when:
  ! - just the third convergence criterion is satisfied;
  ! - both the second and the third convergence criteria are
  !   satisfied, but the squared Euclidian-norm of R(tau*z) is smaller
  !   or equal to l*(delta2 - peucn2)

20 tmp = (rzeucn2 * (tau**2)) - l * (delta2 - peucn2)

  if ( ( msinfo .eq. 3 ) .or. &
       ( ( msinfo .eq. 5 ) .and. ( tmp .le. 0.0d0 ) ) ) then

     peucn2 = 0.0d0
     do i = 1,n
        p(i) = p(i) + tau * z(i)
        peucn2 = peucn2 + p(i)**2
     end do
     peucn = sqrt( peucn2 )

     msinfo = 3
  else
     msinfo = 2
  end if

  ! Print informations

21 if ( iprintinn .ge. 5 ) then
     write(*, 1060) iter
     write(*, 1010) ls,ll,lu,l
     write(*, 1080)
     write(*, 1040) peucn

     write(10,1060) iter
     write(10,1010) ls,ll,lu,l
     write(10,1080)
     write(10,1040) peucn
  end if

22 continue

  ! print *,'CALL LSSEND'
  call lssend()

! Non-executable statements

 1000 format(/,10X,'More-Sorensen iteration: ',I7)
 1010 format(  10X,'ls = ',1P,D11.4, &
               10X,'ll = ',1P,D11.4, &
               10X,'lu = ',1P,D11.4, &
               10X,'l  = ',1P,D11.4)
 1020 format(  10X,'Euclidian-norm of p: ',1P,D7.1)
 1030 format(  10X,'Euclidian-norm of p + z: ',1P,D7.1)
 1040 format(  10X,'Euclidian-norm of step: ',1P,D7.1)
 1050 format(  10X,'delta: ',1P,D7.1)
 1060 format(  10X,'Number of iterations: ',I7)
 1070 format(  10X,'Matrix is not positive definite!')

 1080 format(  10X,'Flag of More-Sorensen: ', &
                   'convergence criterion satisfied.')
 1090 format(  10X,'Flag of More-Sorensen: ', &
                   'maximum number of iterations achieved.')

end subroutine moresor

! ******************************************************************
! ******************************************************************

subroutine scalcu(n,annz,alin,acol,aval,adiag,l,idx,u,ueucn2,wd, &
     memfail)
  
  use lsslvr
  use modalgparam, only: lsssubTR

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(out)   :: memfail
  integer,      intent(in)    :: annz,idx,n
  real(kind=8), intent(in)    :: l
  real(kind=8), intent(inout) :: ueucn2

  ! ARRAY ARGUMENTS
  integer,      intent(inout) :: acol(annz),adiag(n),alin(annz)
  real(kind=8), intent(inout) :: aval(annz),u(idx),wd(idx-1)

  ! Solves a sparse linear system of the form (A + lI + ekek^Td)u = 0,
  ! where A + lI is a definite positive matrix in R^{k x k}.
  ! u is a vector of k positions with u(k) = 1.

  ! LOCAL SCALARS
  integer :: col,i,lin,lssinfo

  if ( idx .eq. 1 ) then
     u(idx) = 1.0d0
     ueucn2 = 1.0d0
     return
  end if

  do i = 1,idx
     u(i) = 0.0d0
  end do

  ! Permute columns and rows of A

  ! print *,'CALL LSSPERMIND'
  call lsspermind(annz,alin,acol)

  ! Eliminate columns that have index greater than idx and define u
  ! as idx-th column of A + lI

  do i = 1,annz
     col = acol(i)
     lin = alin(i)

     if ( ( col .eq. idx ) .and. ( lin .lt. idx ) ) then
        u(lin) = u(lin) - aval(i)
     else if ( ( lin .eq. idx ) .and. ( col .lt. idx ) ) then
        u(col) = u(col) - aval(i)
     end if

     if ( ( col .ge. idx ) .or. ( lin .ge. idx ) ) then
        acol(i) = col + n
        alin(i) = lin + n
     end if
  end do

  ! Solve system (A + lI)x = u

  ! print *,'CALL LSSSET'
  call lssset(lsssubTR,'NONE',.true.,.false.)

  ! print *,'CALL LSSPERMVEC'
  call lsspermvec(n,adiag)

  do i = 1,idx-1
     wd(i) = l
  end do

  ! print *,'CALL LSSAFSOL'
  call lssafsol(idx-1,annz,alin,acol,aval,adiag,wd,u,lssinfo)

  if ( lssinfo .eq. 6 ) then
     ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

     memfail = .true.
     return

  else if ( lssinfo .eq. 7 ) then
     ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

     memfail = .true.
     return

  else if ( lssinfo .eq. 8 ) then
     ! INSUFFICIENT INTEGER WORKING SPACE

     memfail = .true.
     return

  end if

  ! print *,'CALL LSSUNPERMVEC'
  call lssunpermvec(n,adiag)

  ! print *,'CALL LSSSET'
  call lssset(lsssubTR,'NONE',.false.,.true.)

  ! Undo modifications in acol and alin

  do i = 1,annz
     col = acol(i)

     if ( col .gt. n ) then
        acol(i) = acol(i) - n
        alin(i) = alin(i) - n
     end if
  end do

  ! print *,'CALL LSSUNPERMIND'
  call lssunpermind(n,annz,alin,acol)

  u(idx) = 1.0d0

  ueucn2 = 0.0d0
  do i = 1,idx
     ueucn2 = ueucn2 + u(i)**2
  end do

end subroutine scalcu

! ******************************************************************
! ******************************************************************

subroutine scalcz(n,rowind,rowval,z)

  use lsslvr
  use modmachconst

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  integer,      intent(inout) :: rowind(n)
  real(kind=8), intent(inout) :: rowval(n),z(n)

  ! Sparse implementation of the technique presented by Cline, Moler,
  ! Stewart and Wilkinson to estimate the condition number of a upper
  ! triangular matrix.

  ! LOCAL SCALARS
  integer      :: col,i,k,rownnz
  real(kind=8) :: acs,acz,d,ek,rki,s,sm,w,wk,wkm

  ! print *,'CALL LSSSETROW'
  call lsssetrow(n)

  do i = 1,n
     z(i) = 0.0d0
  end do

  acz = 0.0d0
  acs = 1.0d0
  ek  = 1.0d0

  do k = 1,n
     d = lssgetd(k)

     if ( abs( acs*z(k) ) .gt. macheps23 ) then
        ek = dsign(ek,-acs*z(k))
     end if

     if ( abs( ek - acs*z(k) ) .gt. abs( d ) ) then
        s   = abs( d ) / abs( ek - acs*z(k) )
        acs = acs * s
        ek  = s * ek
     end if

     wk  =  ek - acs * z(k)
     wkm = -ek - acs * z(k)
     s   = abs( wk )
     sm  = abs( wkm )

     if ( d .eq. 0.0d0 ) then
        wk  = 1.0d0
        wkm = 1.0d0
     else
        wk  = wk / d
        wkm = wkm / d
     end if

     if ( k .eq. n ) then
        go to 10
     end if

     call lssgetrow(n,k,rownnz,rowind,rowval)

     sm = sm + acs * acz

     do i = 1,rownnz

        col = rowind(i)
        if ( col .gt. k ) then
           rki    = rowval(i)
           sm     = sm - abs( acs * z(col) )
           sm     = sm + abs( acs * z(col) + wkm * rki )
           acz    = acz - abs( z(col) )
           acz    = acz + abs( z(col) + (wk/acs) * rki )
           z(col) = z(col) + (wk/acs) * rki
        end if
     end do

     s = s + acs * acz

     if ( s .lt. &
          sm - macheps23 * max( abs( sm ), 1.0d0 ) ) then
        w = wkm - wk
        wk = wkm

        do i = 1,rownnz

           col = rowind(i)
           if ( col .gt. k ) then
              rki    = rowval(i)
              acz    = acz - abs( z(col) )
              acz    = acz + abs( z(col) + (w/acs) * rki )
              z(col) = z(col) + (w/acs) * rki
           end if
        end do

     end if
     acz = acz - abs( z(k+1) )
10   z(k) = wk/acs
  end do

  ! Divide z by its 1-norm to avoid overflow

  s = 0.0d0
  do i = 1,n
     s = s + abs( z(i) )
  end do

  do i = 1,n
     z(i) = z(i) / s
  end do

  ! Solve Rz = y

  ! print *,'CALL LSSSOLTR'
  call lsssoltr(' ',n,z)

end subroutine scalcz

! ******************************************************************
! ******************************************************************

! moresor:
!
! Method to minimize sparse quadratic functions subjected to
! ||w|| <= delta
!
! minimize    psi(w) = 1/2 w^TBw + g^Tw
! subject to  ||w|| <= delta
!
! Method described in "Computing a trust region step", by More and
! Sorensen.
!
! The main ideia of this method is to find a positive scalar \mslamb
! that is a zero of the function
!
! phi(\mslamb) = 1/||p|| - 1/delta,
!
! where p is the solution of the linear system
!
! (B + \mslamb I)p = -g.
!
! Note that symmetric matrix B, vector g and positive real number
! delta are presented in the minimization problem above. I is the
! identity matrix.
!
! The method used to find the zero of that function is basically the
! Newton method to find roots.
!
! On Entry
!
! n        integer
!          dimension
!
! g        double precision g(n)
!          vector used to define the quadratic function
!
! bnnz     integer
!          number of nonzero elements of B
!
! blin     integer blin(bnnz)
!          row indices of nonzero elements of B
!
! bcol     integer bcol(bnnz)
!          column indices of nonzero elements of B
!
! bval     double precision bval(bnnz)
!          nonzero elements of B, that is,
!          B(blin(i),bcol(i)) = bval(i). Since B is symmetric, just
!          one element B(i,j) or B(j,i) must appear in its sparse
!          representation. If more than one pair corresponding to
!          the same position of B appears in the sparse
!          representation, the multiple entries will be summed.
!          If any blin(i) or bcol(i) is out of range, the entry will
!          be ignored
!
! bdiag    integer bdiag(n)
!          indices of diagonal elements of B in blin, bcol and bval
!
! delta    double precision
!          trust-region radius
!
! sigma1   double precision
!          allowed error for convergence criteria 1 and 2
!
! sigma2   double precision
!          allowed error for convergence criteria 3 (hard case)
!
! eps      double precision
!          allowed error
!
! maxit    integer
!          maximum number of allowed iterations
!
! l        double precision
!          initial value for mslamb
!
! On Return
!
! l        double precision
!          value that gives p as a solution to the minimization
!          problem, because p is also solution to
!          (B + l I)p = -g
!
! pd       logical
!          set to true if the last Cholesky decomposition is
!          successfull
!
! p        double precision p(n)
!          solution to problem
!          minimize     psi(w)
!          subjected to ||w|| <= delta
!
! chcnt    integer
!          number of Cholesky decompositions
!
! memfail  logical
!          true iff linear solver failed because of lack of memory
!
! msinfo   integer
!          stores which convergence criteria was satisfied:
!
!          0 = both g and B are null;
!
!          1 = first convergence criterion is satisfied;
!
!          2 = second convergence criterion is satisfied;
!
!          3 = third convergence criterion is satisfied
!
!          5 = maximum allowed number of iterations is achieved.

! ******************************************************************
! ******************************************************************

! scalcu:
!
! Solve a sparse linear system of the form (A + lI + ekek^Td)u = 0,
! where A + lI is a definite positive matrix in R^{k x k}.
! u is a vector of k positions with u(k) = 1.
!
! On Entry
!
! n        integer
!          dimension of A
!
! annz     integer
!          number of nonzero elements of A
!
! alin     integer alin(annz)
!          row indices of nonzero elements of A
!
! acol     integer acol(annz)
!          column indices of nonzero elements of A
!
! aval     integer aval(annz)
!          nonzero elements of A, that is,
!          A(alin(i),acol(i)) = aval(i). Since A is symmetric, just
!          one element A(i,j) or A(j,i) must appear in its sparse
!          representation. If more than one pair corresponding to
!          the same position of A appears in the sparse
!          representation, the multiple entries will be summed.
!          If any alin(i) or acol(i) is out of range, the entry will
!          be ignored
!
! adiag    integer adiag(n)
!          indices of diagonal elements of A in alin, acol and aval
!
! l        double precision
!          used to compute A + lI
!
! idx      integer
!          index k
!
! On Return
!
! u        double precision u(n)
!          system solution
!
! ueucn2   double precision
!          u squared Euclidian-norm
!
! memfail  logical
!          true iff linear solver failed because of lack of memory

! ******************************************************************
! ******************************************************************

! scalcz:
!
! Sparse implementation of the technique presented by Cline, Moler,
! Stewart and Wilkinson to estimate the condition number of a matrix.
! This technique is used by the More-Sorensen method to calculate an
! approximation to the eigenvector associated to the smallest
! eigenvalue of a matrix B (\mslamb_1).
! In this technique, when \mslamb approaches -\mslamb_1, \|Rz\|
! approaches 0. This insures that z is an approximation to the
! wanted eigenvector.
! Basically, it solves R^Ty = e, choosing e(k) as 1 or -1 (whatever
! gives maximum local growth of y). Then, it solves Rz = y.
! Note that R is a sparse matrix given by P^T D^0.5 L^T P (obtained
! applying subroutine MA27BD).
!
! On Entry
!
! n        integer
!          dimension of R
!
! rowind   integer rowind(n)
! rowval   double precision rowval(n)
!          working arrays
!
! On Return
!
! z        double precision z(n)
!          approximation of the eigenvector of B associated to
!          \mslamb_1
