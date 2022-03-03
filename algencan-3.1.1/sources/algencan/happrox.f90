! *****************************************************************
! *****************************************************************

subroutine comphapp(n,m,rho,equatn)

  use modalgconst
  use modsvdgrad
  use modsydat
  use modhappdat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: m,n

  ! ARRAY ARGUMENTS
  logical,      intent(in) :: equatn(m)
  real(kind=8), intent(in) :: rho(m)

  ! This subroutine computes an approximation H of the Hessian of the
  ! Augmented Lagrangian following a very simple idea: "discard the
  ! second order terms and then correct the remaining matrix in order
  ! to satisfy a secant equation".
  !
  ! Hence, H takes the form
  !
  ! H = B + S + rho A^t A,
  !
  ! where S is the spectral correction of (rho A^t A) and B is
  ! the BFGS correction of (S + rho A^t A). More specifically,
  !
  ! S = hlspg I,
  !
  ! where
  !
  ! hlspg = max(lspgmi, min(lspgma, s^t (y - rho A^t A s) / s^t s)),
  !
  ! D = S + rho A^t A,
  !
  ! and
  !
  ! B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ].
  !
  ! Note that this subroutine does not compute matrix H explicitly,
  ! but computes some quantities that will be used later to compute
  ! the product of H by a vector p.
  !
  ! The quantities computed by this subroutine are:
  !
  ! (a) hlspg = s^t (y - rho A^t A s) / (s^t s)
  !
  ! (b) hds = D s = ( hlspg I + rho A^t A ) s, and
  !
  ! (c) hstds = <s,hds>.

  ! LOCAL SCALARS
  integer      :: i,j
  real(kind=8) :: ats

  ! ------------------------------------------------------------------
  ! Compute hds = rho A^t A s
  ! ------------------------------------------------------------------

  do i = 1,n
     hds(i) = 0.0d0
  end do

  do j = 1,m

     if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then

        ! COMPUTE THE INNER PRODUCT <a,s>
        ats = 0.0d0
        do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
           ats = ats + svdjcval(i) * s(svdjcvar(i))
        end do

        ats = ats * rho(j)

        ! ADD rho * ats * a
        do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
           hds(svdjcvar(i)) = hds(svdjcvar(i)) + ats * svdjcval(i)
        end do

     end if

  end do

  hstds = 0.0d0
  do i = 1,n
     hstds = hstds + s(i) * hds(i)
  end do

  ! ------------------------------------------------------------------
  ! Compute hlspg = s^t (y - rho A^t A s) / (s^t s)
  ! ------------------------------------------------------------------

  if ( sty - hstds .le. 0.0d0 ) then
     hlspg = lspgmi
  else
     hlspg = max( lspgmi, min( (sty - hstds) / sts, lspgma ) )
  end if

  do i = 1,n
     hds(i) = hds(i) + hlspg * s(i)
  end do

  hstds = hstds + hlspg * sts

end subroutine comphapp

! *****************************************************************
! *****************************************************************

subroutine applyhapp(n,m,rho,equatn,goth,p,hp)

  use modmachconst
  use modsvdgrad
  use modsydat
  use moditetyp
  use modhappdat

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: goth
  integer, intent(in)    :: m,n

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m)
  real(kind=8), intent(in)    :: rho(m)
  real(kind=8), intent(inout) :: hp(n),p(n)

  ! This subroutine computes the product of the matrix computed by
  ! subroutine comphapp times vector p.

  ! LOCAL SCALARS
  integer      :: i,j
  real(kind=8) :: atp,c1,c2,ptds,pty

  ! ------------------------------------------------------------------
  ! Compute Hessian approximation
  ! ------------------------------------------------------------------

  if ( .not. goth ) then
     goth = .true.
     call comphapp(n,m,rho,equatn)
  end if

  ! ------------------------------------------------------------------
  ! Compute ( hlspg I ) p
  ! ------------------------------------------------------------------

  do i = 1,n
     hp(i) = hlspg * p(i)
  end do

  ! ------------------------------------------------------------------
  ! Add ( rho A^T A ) p
  ! ------------------------------------------------------------------

  do j = 1,m

     if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then

        ! COMPUTE THE INNER PRODUCT <a,p>
        atp = 0.0d0
        do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
           atp = atp + svdjcval(i) * p(svdjcvar(i))
        end do

        atp = atp * rho(j)

        ! ADD rho * atp * a
        do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
           hp(svdjcvar(i)) = hp(svdjcvar(i)) + atp * svdjcval(i)
        end do

     end if

  end do

  ! ------------------------------------------------------------------
  ! Add B p,
  ! where B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ]
  ! ------------------------------------------------------------------

  if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

     pty = 0.0d0
     ptds = 0.0d0
     do i = 1,n
        pty = pty + p(i) * y(i)
        ptds = ptds + p(i) * hds(i)
     end do

     c1 = pty / sty
     c2 = ptds / hstds
     do i = 1,n
        hp(i) = hp(i) + c1 * y(i) - c2 * hds(i)
     end do

  end if

end subroutine applyhapp

! *****************************************************************
! *****************************************************************

subroutine comphpre(n,m,rho,equatn)

  use modmachconst
  use modalgconst
  use modsvdgrad
  use modsydat
  use moditetyp
  use modhpredat

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: m,n

  ! ARRAY ARGUMENTS
  logical,      intent(in) :: equatn(m)
  real(kind=8), intent(in) :: rho(m)

  ! Consider the preconditioner
  !
  !     P = Q + E + diag(rho A^t A)
  !
  ! for matrix
  !
  !     H = B + S + rho A^t A,
  !
  ! where E is the spectral correction of diag(rho A^t A) and Q is
  ! the BFGS correction of (E + diag(rho A^t A)), while S and B are
  ! the spectral and BFGS corrections of matrix (rho A^t A),
  ! respectively.
  !
  ! This subroutine computes:
  !
  ! (a) pdiag = diag(rho A^t A),
  !
  ! (b) plspg such that E = plspg I,
  !
  ! (c) psmdy = s - D^-1 y, where D = E + diag(rho A^t A), and
  !
  ! (d) the inner product psmdty = <psmdy,y>.
  !
  ! These quantities will be used later, in subroutine applyp, to
  ! compute z = P^{-1} r.

  ! LOCAL SCALARS
  integer      :: i,j
  real(kind=8) :: sttmp

  ! ------------------------------------------------------------------
  ! Compute diag( rho A^t A )
  ! ------------------------------------------------------------------

  do i = 1,n
     pdiag(i) = 0.0d0
  end do

  do j = 1,m
     if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then
        do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
           pdiag(svdjcvar(i)) = pdiag(svdjcvar(i)) + rho(j) * svdjcval(i) ** 2
        end do
     end if
  end do

  ! ------------------------------------------------------------------
  ! Compute plspg = s^t (y - diag( rho A^t A ) s) / (s^t s)
  ! ------------------------------------------------------------------

  sttmp = 0.0d0
  do i = 1,n
     sttmp = sttmp + pdiag(i) * s(i) ** 2
  end do

  if ( sty - sttmp .le. 0.0d0 ) then
     plspg = lspgmi
  else
     plspg = max( lspgmi, min( ( sty - sttmp ) / sts, lspgma ) )
  end if

  ! ------------------------------------------------------------------
  ! Compute the BFGS correction Q of ( E + diag( rho A^t A ) )
  !
  ! Q = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
  !     [ <s - D^-1 y, y> s s^t ] / (s^t y)^2,
  !
  ! where D = ( E + diag( rho A^t A ) )
  ! ------------------------------------------------------------------

  if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

     psmdyty = 0.0d0
     do i = 1,n
        psmdy(i) = s(i) - y(i) / ( plspg + pdiag(i) )
        psmdyty = psmdyty + psmdy(i) * y(i)
     end do

  end if

end subroutine comphpre

! *****************************************************************
! *****************************************************************

subroutine applyhpre(n,m,rho,equatn,gotp,r,z)

  use modmachconst
  use modalgparam, only: innercall
  use modsydat
  use moditetyp
  use modhpredat

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: gotp
  integer, intent(in)    :: m,n

  ! ARRAY ARGUMENTS
  logical,      intent(in)  :: equatn(m)
  real(kind=8), intent(in)  :: rho(m),r(n)
  real(kind=8), intent(out) :: z(n)

  ! This subroutine computes the product of the inverse of the matrix
  ! computed by subroutine comphpre times vector r, i.e., z = P^{-1} r.

  ! LOCAL SCALARS
  integer      :: i
  real(kind=8) :: c1,c2,psmdytr,str

!!$  if ( innercall ) then
!!$     z = r
!!$     return
!!$  end if

  ! ------------------------------------------------------------------
  ! Compute P
  ! ------------------------------------------------------------------

  if ( .not. gotp ) then
     gotp = .true.
     call comphpre(n,m,rho,equatn)
  end if

  ! ------------------------------------------------------------------
  ! Compute ( E + diag( rho A^T A ) )^{-1} r
  ! ------------------------------------------------------------------

  do i = 1,n
     z(i) = r(i) / ( plspg + pdiag(i) )
  end do

  ! ------------------------------------------------------------------
  ! Add Q^{-1} r, where
  !
  ! Q^{-1} = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
  !          [ <s - D^-1 y, y> s s^t ] / (s^t y)^2
  !
  ! and D = ( E + diag( rho A^T A ) )
  ! ------------------------------------------------------------------

  if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

     str = 0.0d0
     psmdytr = 0.0d0
     do i = 1,n
        str = str + s(i) * r(i)
        psmdytr = psmdytr + psmdy(i) * r(i)
     end do

     c1 = str / sty
     c2 = psmdytr / sty - psmdyty * str / sty ** 2

     do i = 1,n
        z(i) = z(i) + c1 * psmdy(i) + c2 * s(i)
     end do

  end if

end subroutine applyhpre
