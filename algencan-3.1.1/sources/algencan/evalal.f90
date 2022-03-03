! *****************************************************************
! *****************************************************************

subroutine sevalal(n,x,m,lambda,rho,equatn,linear,al,inform)

  use modalgparam, only: innercall
  use modsvdgrad
  use problvls, only: sevalc,sevalf,sevalfc
  use probgiven, only: ccoded,fcoded,fccoded

  implicit none

  ! SCALAR ARGUMENTS
  real(kind=8), intent(out)   :: al
  integer,      intent(in)    :: m,n
  integer,      intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m),rho(m)
  real(kind=8), intent(inout) :: x(n)

  ! LOCAL SCALARS
  integer      :: j
  real(kind=8) :: f,p

  ! LOCAL ARRAYS
  logical :: ind(m)
  
  if ( innercall ) then
     call minsqf(n,x,al,inform)
     return
  end if

  if ( fccoded ) then

     ! COMPUTE OBJECTIVE FUNTION AND CONSTRAINTS
     call sevalfc(n,x,f,m,svdc,inform)
     if ( inform .ne. 0 ) return

     ! COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)
     al = f

     do j = 1,m
        ! ADD P(c_j, rho_j, lambda_j)
        call evalp(svdc(j),rho(j),lambda(j),equatn(j),p)
        al = al + p
     end do

  else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

     ! COMPUTE OBJECTIVE FUNCTION
     call sevalf(n,x,f,inform)
     if ( inform .ne. 0 ) return

     ! COMPUTE CONSTRAINTS
     ind(1:m) = .true.
     call sevalc(n,x,m,ind,svdc,inform)
     
     ! COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)

     al = f
     do j = 1,m
        ! ADD P(c_j, rho_j, lambda_j)
        call evalp(svdc(j),rho(j),lambda(j),equatn(j),p)
        al = al + p
     end do

  end if

  svdgotc = .true.

end subroutine sevalal

! *****************************************************************
! *****************************************************************

subroutine sevalnl(n,x,m,lambda,equatn,linear,nl,inform)

  use modsvdgrad
  use problvls, only: sevalg,sevalgjacp,sevalgjac,sevaljac
  use probgiven, only: ccoded,fcoded,fccoded,gjaccoded,gjacpcoded,jcnnzlim

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m)
  real(kind=8), intent(inout) :: nl(n),x(n)

  ! LOCAL SCALARS
  logical      :: alllin,allnol,gotj
  integer      :: i,j,jcnnz

  ! LOCAL ARRAYS
  logical      :: ind(m)
  integer      :: jcfun(jcnnzlim)
  real(kind=8) :: dum(1),p(m)

  if ( fccoded ) then

     if ( gjacpcoded ) then

        ! COMPUTE p
        svdconstrc = .false.
        do j = 1,m
           if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then
              p(j) = lambda(j)
              svdconstrc = .true.
           else
              p(j) = 0.0d0
           end if
        end do

        ! COMPUTE nl = Jacobian^T p
        gotj = .false.
        call sevalgjacp(n,x,svdg,m,p,nl,'T',gotj,inform)
        if ( inform .ne. 0 ) return

        ! COMPUTE gparc = Jacobian^t p IGNORING LINEAR CONSTRAINTS
        allnol = .true.
        alllin = .true.
        do j = 1,m
           if ( linear(j) ) then
              p(j) = 0.0d0
              allnol = .false.
           else
              alllin = .false.
           end if
        end do

        if ( allnol ) then
           do i = 1,n
              svdgparc(i) = nl(i)
           end do
        else if ( alllin ) then
           do i = 1,n
              svdgparc(i) = 0.0d0
           end do
        else
           call sevalgjacp(n,x,dum,m,p,svdgparc,'t',gotj,inform)
           if ( inform .ne. 0 ) return
        end if

        ! ADD THE GRADIENT OF THE OBJECTIVE FUNCTION
        do i = 1,n
           nl(i)       = nl(i)       + svdg(i)
           svdgparc(i) = svdgparc(i) + svdg(i)
        end do

     else ! if ( gjaccoded ) then
        ! In fact, this else is the choice for gjaccoded or 'not 
        ! gjacpcoded and not gjaccoded', in which case the calling 
        ! sequence starting with sevalgjac below will end up using 
        ! finite differences to approximate first derivatives of the 
        ! objective function and the constraints

        ! COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND 
        ! JACOBIAN OF CONSTRAINTS
        call sevalgjac(n,x,svdg,m,jcfun,svdjcvar,svdjcval,jcnnz,jcnnzlim,inform)
        if ( inform .ne. 0 ) return

        ! CONVERT JACOBIAN OF CONSTRAINTS FROM COORDINATE FORMAT TO
        ! COMPRESSED SPARSE ROW FORMAT
        call coo2csr(m,jcnnz,jcfun,svdjcvar,svdjcval,svdjclen,svdjcsta)

        ! COMPUTE \nabla L = \nabla f + \sum_j lambda_j * \nabla c_j
        svdconstrc = .false.

        nl(1:n)       = svdg(1:n)
        svdgparc(1:n) = svdg(1:n)

        do j = 1,m
           if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then
              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                 nl(svdjcvar(i)) = nl(svdjcvar(i)) + lambda(j) * svdjcval(i)
              end do

              if ( .not. linear(j) ) then
                 do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                    svdgparc(svdjcvar(i)) = svdgparc(svdjcvar(i)) + lambda(j) * svdjcval(i)
                 end do
              end if

              svdconstrc = .true.
           end if
        end do

     end if

  else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

     ! COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION
     call sevalg(n,x,svdg,inform)
     if ( inform .ne. 0 ) return

     nl(1:n)       = svdg(1:n)
     svdgparc(1:n) = svdg(1:n)

     ! ADD nl = Jacobian^T lambda

     ! SET LOGICAL ARRAY WITH REQUIRED GRADIENTS OF CONSTRAINTS
     ind(1:m) = ( equatn(1:m) .or. lambda(1:m) .gt. 0.0d0 )

     svdconstrc = any( ind(1:m) )

     ! COMPUTE REQUIRED GRADIENTS OF CONSTRAINTS
     call sevaljac(n,x,m,ind,svdjcsta,svdjclen,svdjcvar,svdjcval,jcnnzlim,inform)
     if ( inform .ne. 0 ) return
     
     do j = 1,m
        if ( ind(j) ) then
           
           ! ADD lambda_j * \nabla c_j
           do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
              nl(svdjcvar(i)) = nl(svdjcvar(i)) + lambda(j) * svdjcval(i)
           end do

           if ( .not. linear(j) ) then
              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                 svdgparc(svdjcvar(i)) = svdgparc(svdjcvar(i)) + lambda(j) * svdjcval(i)
              end do
           end if

        end if
     end do

  end if

end subroutine sevalnl

! *****************************************************************
! *****************************************************************

subroutine sevalnal(n,x,m,lambda,rho,equatn,linear,nal,inform)

  use modalgparam, only: innercall
  use modsvdgrad
  use problvls, only: sevalc,sevalfc
  use probgiven, only: ccoded,fccoded,fcoded

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m),rho(m)
  real(kind=8), intent(inout) :: nal(n),x(n)

  ! LOCAL SCALARS
  integer      :: j
  real(kind=8) :: dum

  ! LOCAL ARRAYS
  logical :: ind(m)
  
  if ( innercall ) then
     call minsqg(n,x,nal,inform)
     return
  end if

  if ( fccoded ) then
     ! COMPUTE CONSTRAINTS
     if ( .not. svdgotc .and. m .gt. 0 ) then
        call sevalfc(n,x,dum,m,svdc,inform)
        if ( inform .ne. 0 ) return
     end if

  else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then
     ! COMPUTE CONSTRAINTS
     if ( .not. svdgotc ) then
        ind(1:m) = .true.
        call sevalc(n,x,m,ind,svdc,inform)
        if ( inform .ne. 0 ) return
     end if
  end if

  svdgotc = .true.

  ! COMPUTE dP/dc
  do j = 1,m
     call evaldpdy(svdc(j),rho(j),lambda(j),equatn(j),svddpdc(j))
  end do

  ! COMPUTE GRADIENT OF THE LAGRANGIAN WITH DPDC INSTEAD OF LAMBDA
  call sevalnl(n,x,m,svddpdc,equatn,linear,nal,inform)
  if ( inform .ne. 0 ) return

end subroutine sevalnal

! *****************************************************************
! *****************************************************************

subroutine evalp(y,rho,lambda,equatn,p)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(in)  :: equatn
  real(kind=8), intent(in)  :: lambda,rho,y
  real(kind=8), intent(out) :: p

  if ( equatn ) then
     p  = y * ( lambda + 0.5d0 * rho * y )
  else
     if ( lambda + rho * y .ge. 0.0d0 ) then
        p = y * ( lambda + 0.5d0 * rho * y )
     else
        p = - 0.5d0 * lambda ** 2 / rho
     end if
  end if

end subroutine evalp

! *****************************************************************
! *****************************************************************

subroutine evaldpdy(y,rho,lambda,equatn,dpdy)

  implicit none

  ! SCALAR ARGUMENTS
  logical,      intent(in)  :: equatn
  real(kind=8), intent(in)  :: y,rho,lambda
  real(kind=8), intent(out) :: dpdy

  if ( equatn ) then
     dpdy = lambda + rho * y
  else
     dpdy = max( 0.0d0, lambda + rho * y )
  end if

end subroutine evaldpdy

! *****************************************************************
! *****************************************************************

subroutine ievalnal(n,xp,m,lambda,rho,equatn,linear,iglin,nalp,inform)

  use problvls, only: sevalc,sevalfc,sevalg,sevalgjac,sevalgjacp,sevaljac
  use probgiven, only: ccoded,fccoded,fcoded,gjaccoded,gjacpcoded,jcnnzlim
  use modalgparam, only: innercall
  use modsvdgrad

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(in)    :: iglin
  integer, intent(in)    :: m,n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m),rho(m)
  real(kind=8), intent(inout) :: xp(n)
  real(kind=8), intent(out)   :: nalp(n)

  ! This subroutine computes the gradient of the Augmented Lagrangian
  ! function at a point xp, which is near to x, taking care of the
  ! non-differentiability (i.e. considering the contribution of the 
  ! same constraints that contributed to compute the gradient of the 
  ! augmented Lagrangian at x). The augmented Lagrangian gradient must 
  ! have been previously computed at x. If iglin = true, linear 
  ! constraints are ignored.

  ! LOCAL SCALARS
  logical      :: alllin,allnol,gotj
  integer      :: i,j,jcpnnz
  real(kind=8) :: cpj,dpdcpj,dum

  ! LOCAL ARRAYS
  logical      :: ind(m)
  integer      :: jcpfun(jcnnzlim),jcplen(m),jcpsta(m),jcpvar(jcnnzlim)
  real(kind=8) :: cp(m),gp(n),jcpval(jcnnzlim),dpdcp(m)

  if ( innercall ) then
     call minsqg(n,xp,nalp,inform)
     return
  end if

  if ( fccoded ) then

     ! COMPUTE nalp = gp + Jacobian^T dpdcp IGNORING THE LINEAR 
     ! CONSTRAINTS IF iglin IS TRUE

     if ( gjacpcoded ) then

        ! COMPUTE CONSTRAINTS AT xp
        if ( m .gt. 0 ) then
           call sevalfc(n,xp,dum,m,cp,inform)
           if ( inform .ne. 0 ) return
        end if

        ! COMPUTE dpdcp
        do j = 1,m
           if ( ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) .and. &
                .not. ( iglin .and. linear(j) ) ) then
              dpdcp(j) = lambda(j) + rho(j) * cp(j)
           else
              dpdcp(j) = 0.0d0
           end if
        end do

        ! COMPUTE nalp = Jacobian^T dpdcp
        gotj = .false.
        call sevalgjacp(n,xp,gp,m,dpdcp,nalp,'T',gotj,inform)
        if ( inform .ne. 0 ) return

        ! ADD THE GRADIENT OF THE OBJECTIVE FUNCTION
        nalp(1:n) = nalp(1:n) + gp(1:n)

     else ! if ( gjaccoded ) then
        ! In fact, this else is the choice for gjaccoded or 'not 
        ! gjacpcoded and not gjaccoded', in which case the calling 
        ! sequence starting with sevalgjac below will end up using 
        ! finite differences to approximate first derivatives of the 
        ! objective function and the constraints

        ! COMPUTE CONSTRAINTS AT xp
        if ( m .gt. 0 ) then
           call sevalfc(n,xp,dum,m,cp,inform)
           if ( inform .ne. 0 ) return
        end if

        ! COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND JACOBIAN
        ! OF CONSTRAINTS AT xp
        call sevalgjac(n,xp,gp,m,jcpfun,jcpvar,jcpval,jcpnnz,jcnnzlim,inform)
        if ( inform .ne. 0 ) return

        nalp(1:n) = gp(1:n)

        ! CONVERT JACOBIAN OF CONSTRAINTS FROM COODINATE FORMAT TO
        ! COMPRESSED SPARSE ROW FORMAT
        call coo2csr(m,jcpnnz,jcpfun,jcpvar,jcpval,jcplen,jcpsta)

        ! COMPUTE dpdcp AND ADD dPdc * dcdx
        do j = 1,m
           if ( ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) .and. &
                .not. ( iglin .and. linear(j) ) ) then
              dpdcpj = lambda(j) + rho(j) * cp(j)

              do i = jcpsta(j),jcpsta(j) + jcplen(j) - 1
                 nalp(jcpvar(i)) = nalp(jcpvar(i)) + dpdcpj * jcpval(i)
              end do
           end if
        end do
     end if

  else if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then

     ! COMPUTE GRADIENT OF THE OBJECTIVE FUNCTION
     call sevalg(n,xp,gp,inform)
     if ( inform .ne. 0 ) return

     ! SET LOGICAL ARRAY WITH CONSTRAINTS (AND THEIR GRADIENTS) THAT MUST BE COMPUTED
     ind(1:m) = ( equatn(1:m) .or. svddpdc(1:m) .gt. 0.0d0 ) .and. &
          .not. ( iglin .and. linear(1:m) )
     
     ! COMPUTE SELECTED CONSTRAINTS (AND THEIR GRADIENTS)
     call sevalc(n,xp,m,ind,cp,inform)
     if ( inform .ne. 0 ) return

     call sevaljac(n,xp,m,ind,jcpsta,jcplen,jcpvar,jcpval,jcnnzlim,inform)
     if ( inform .ne. 0 ) return

     ! COMPUTE g + Jacobian^T dpdcp IGNORING THE LINEAR CONSTRAINTS 
     ! IF iglin IS TRUE

     nalp(1:n) = gp(1:n)

     do j = 1,m
        if ( ind(j) ) then
           dpdcpj = lambda(j) + rho(j) * cp(j)
           
           do i = jcpsta(j),jcpsta(j) + jcplen(j) - 1
              nalp(jcpvar(i)) = nalp(jcpvar(i)) + dpdcpj * jcpval(i)
           end do
        end if
     end do
  end if

end subroutine ievalnal

! ******************************************************************
! ******************************************************************

subroutine sevalhal2(n,x,m,lambda,rho,equatn,linear,halrow,halcol, &
     halval,halnnz,lim,inform)

  use modsvdgrad
  use problvls, only: sevalhl
  use problvlv, only: reperr

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: lim,m,n
  integer, intent(inout) :: halnnz,inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  integer,      intent(inout) :: halrow(lim),halcol(lim)
  real(kind=8), intent(in)    :: lambda(m),rho(m),x(n)
  real(kind=8), intent(inout) :: halval(lim)

  ! This subroutine computes the Hessian of the augmented Lagrangian.

  ! LOCAL SCALARS
  integer :: j,k,nn,l,var
  real(kind=8) :: val

  ! LOCAL ARRAYS
  character(len=15) :: subname
  integer :: pos(n),rvar(n)
  real(kind=8) ::rval(n)

  call sevalhl(n,x,m,svddpdc,halrow,halcol,halval,halnnz,lim,inform)
  if ( inform .ne. 0 ) return

  if ( m .eq. 0 ) return

  pos(1:n) = -1

  ! ADD \sum_j \rho_j * \nabla c_j \nabla c_j^t

  do j = 1,m
     if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then

        ! ADD \rho_j * \nabla c_j \nabla c_j^t

        nn = 0
        do k = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
           var = svdjcvar(k)
           val = svdjcval(k)
           if ( val .ne. 0.0d0 ) then
              if ( pos(var) .eq. -1 ) then
                 nn = nn + 1
                 pos(var) = nn
                 rvar(nn) = var
                 rval(nn) = val
              else
                 rval(pos(var)) = rval(pos(var)) + val
              end if
           end if
        end do

        pos(rvar(1:nn)) = -1

        do k = 1,nn
           do l = 1,nn
              if ( rvar(k) .ge. rvar(l) ) then
                 halnnz = halnnz + 1
                 if ( halnnz .gt. lim ) then
                    inform = - 92
                    subname = 'SEVALHAL'
                    call reperr(inform,subname)
                    return
                 end if
                 halrow(halnnz) = rvar(k)
                 halcol(halnnz) = rvar(l)
                 halval(halnnz) = rho(j) * rval(k) * rval(l)
              end if
           end do
        end do

     end if
  end do

end subroutine sevalhal2

! ******************************************************************
! ******************************************************************

subroutine sevalhal(n,x,m,lambda,rho,equatn,linear,halrow,halcol, &
     halval,halnnz,lim,inform)

  use modouttyp
  use modsvdgrad
  use problvls, only: sevalhl
  use problvlv, only: reperr

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: lim,m,n
  integer, intent(inout) :: halnnz,inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  integer,      intent(inout) :: halrow(lim),halcol(lim)
  real(kind=8), intent(in)    :: lambda(m),rho(m),x(n)
  real(kind=8), intent(inout) :: halval(lim)

  ! This subroutine computes the Hessian of the augmented Lagrangian.

  ! LOCAL SCALARS
  integer :: i,row,j,k,l,var

  ! LOCAL ARRAYS
  character(len=15) :: subname
  integer           :: strow(n)
  real(kind=8)      :: r(n)

  call sevalhl(n,x,m,svddpdc,halrow,halcol,halval,halnnz,lim,inform)
  if ( inform .ne. 0 ) return

  if ( m .eq. 0 ) return

  ! PUT MATRIX INTO A ROW-LINKED LIST

  r(1:n) = 0.0d0

  strow(1:n) = 0

  do i = 1,halnnz
     row = halrow(i)
     k   = strow(row)
     strow(row) = i
     halrow(i)  = k
  end do

  ! ADD \sum_j \rho_j * \nabla c_j \nabla c_j^t

  do j = 1,m

     if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then

        ! ADD \rho_j * \nabla c_j \nabla c_j^t

        do k = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1

           var = svdjcvar(k)

           ! PUT ROW svdjcvar(k) INTO A DENSE VECTOR

           row = strow(var)
           do while ( row .ne. 0 )
              r(halcol(row)) = r(halcol(row)) + halval(row)
              row = halrow(row)
           end do

           ! ADD VALUE

           do l = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
              if ( svdjcvar(l) .le. var ) then
                 r(svdjcvar(l)) = r(svdjcvar(l)) + rho(j) * svdjcval(l) * svdjcval(k)
              end if
           end do

           ! UPDATE VALUES IN HALVAL

           row = strow(var)
           do while ( row .ne. 0 )
              halval(row) = r(halcol(row))
              r(halcol(row)) = 0.0d0
              row = halrow(row)
           end do

           ! INSERT NEW ELEMENTS IN HAL REPRESENTATION

           do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
              l = svdjcvar(i)
              if ( r(l) .ne. 0.0d0 ) then

                 if ( halnnz + 1 .gt. lim ) then
                    inform = - 92
                    subname = 'SEVALHAL'
                    call reperr(inform,subname)
                    return
                 end if

                 halnnz = halnnz + 1
                 halval(halnnz) = r(l)
                 halcol(halnnz) = l
                 halrow(halnnz) = strow(var)
                 strow(var) = halnnz
                 r(l) = 0.0d0
              end if
           end do

        end do

     end if

  end do

  ! PUT MATRIX BACK INTO COORDINATE SQUEME

  do i = 1,n
     row = strow(i)
     do while ( row .ne. 0 )
        k = halrow(row)
        halrow(row) = i
        row = k
     end do
  end do

end subroutine sevalhal

! *****************************************************************
! *****************************************************************

subroutine sevalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,gothl,inform)

  use modalgparam, only: hptype,innercall
  use modsvdgrad, only: svdconstrc,svddpdc,svdjcsta,svdjclen, &
       svdjcval,svdjcvar
  use problvls, only: sevalhlp,sevalgjacp
  use probgiven, only: gjacpcoded

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out)   :: gothl
  integer, intent(in)    :: m,n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m),rho(m),x(n)
  real(kind=8), intent(inout) :: hp(n),p(n)

  ! LOCAL SCALARS
  logical      :: gotj
  integer      :: i,j
  real(kind=8) :: ajtp

  ! LOCAL ARRAYS
  real(kind=8) :: ap(m),atap(n),dum(1)

  if ( innercall ) then
     call minsqhp(n,x,p,hp,gothl,inform)
     return
  end if

  ! --------------------------------------------------------------
  ! Hessian approximation
  ! --------------------------------------------------------------

  if ( hptype .eq. 'HAPPRO' .and. svdconstrc ) then

     call applyhapp(n,m,rho,equatn,gothl,p,hp)

     ! --------------------------------------------------------------
     ! Incremental quotients
     ! --------------------------------------------------------------

  else if ( hptype .eq. 'INCQUO' .or. hptype .eq. 'HAPPRO' ) then

     call ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)
     if ( inform .ne. 0 ) return

     ! --------------------------------------------------------------
     ! True Hessian
     ! --------------------------------------------------------------

  else if ( hptype .eq. 'TRUEHP' ) then

     ! Compute Hessian of Lagrangian times p using svddpdc
     ! instead of lambda

     call sevalhlp(n,x,m,svddpdc,p,hp,gothl,inform)
     if ( inform .ne. 0 ) return

     ! Add rho A^T A p

     if ( gjacpcoded ) then
        where ( equatn(1:m) .or. svddpdc(1:m) .gt. 0.0d0 )
           ap(1:m) = 1.0d0
        elsewhere
           ap(1:m) = 0.0d0
        end where
        
        gotj = .false.
        call sevalgjacp(n,x,dum,m,ap,p,'j',gotj,inform)
        if ( inform .ne. 0 ) return

        where ( equatn(1:m) .or. svddpdc(1:m) .gt. 0.0d0 )
           ap(1:m) = ap(1:m) * rho(1:m)
        elsewhere
           ! Just for safety (since these position should have been
           ! untouched in the call to sevalgjacp above).
           ap(1:m) = 0.0d0
        end where

        call sevalgjacp(n,x,dum,m,ap,atap,'t',gotj,inform)
        if ( inform .ne. 0 ) return

        hp(1:n) = hp(1:n) + atap(1:n)
        
     else
        do j = 1,m
           if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then

              ajtp = 0.0d0
              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                 ajtp = ajtp + svdjcval(i) * p(svdjcvar(i))
              end do

              ajtp = ajtp * rho(j)

              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                 hp(svdjcvar(i)) = hp(svdjcvar(i)) + ajtp * svdjcval(i)
              end do

           end if
        end do
     end if

  end if

end subroutine sevalhalp

! *****************************************************************
! *****************************************************************

subroutine ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)

  use modmachconst
  use modalgparam
  use modsvdgrad
  use problvls, only: ssetp,sevalgjacp
  use probgiven, only: gjacpcoded

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m),rho(m),x(n)
  real(kind=8), intent(inout) :: hp(n),p(n)

  ! Computes an approximation of the product of the Hessian of the
  ! Augmented Lagrangian times a vector using incremental quotients.

  ! LOCAL SCALARS
  logical      :: gotj,iglin
  integer      :: i,j
  real(kind=8) :: ajtp,psupn,step,xsupn

  ! LOCAL ARRAYS
  real(kind=8) :: ap(m),atap(n),dum(1),gpparc(n),xp(n)

  ! ------------------------------------------------------------------
  ! Set auxiliary point
  ! ------------------------------------------------------------------

  xsupn = 0.0d0
  psupn = 0.0d0
  do i = 1,n
     xsupn = max( xsupn, abs( x(i) ) )
     psupn = max( psupn, abs( p(i) ) )
  end do

  step = macheps12 * max( xsupn / psupn, 1.0d0 )

  do i = 1,n
     xp(i) = x(i) + step * p(i)
  end do

  ! ------------------------------------------------------------------
  ! Compute gradient of the augmented Lagrangian at xp considering the
  ! same constraints considered at x and ignoring linear constraints
  ! ------------------------------------------------------------------

  call ssetp(n,xp,inform)
  if ( inform .ne. 0 ) return

  iglin = .true.
  call ievalnal(n,xp,m,lambda,rho,equatn,linear,iglin,gpparc,inform)
  if ( inform .ne. 0 ) return

  ! ------------------------------------------------------------------
  ! Compute gradients difference
  ! ------------------------------------------------------------------

  do i = 1,n
     hp(i) = ( gpparc(i) - svdgparc(i) ) / step
  end do

  ! ------------------------------------------------------------------
  ! Add contribution of linear constraints
  ! ------------------------------------------------------------------

  if ( gjacpcoded ) then

     where ( linear(1:m) )
        ap(1:m) = 1.0d0
     elsewhere
        ap(1:m) = 0.0d0
     end where
     
     gotj = .false.
     call sevalgjacp(n,x,dum,m,ap,p,'j',gotj,inform)
     if ( inform .ne. 0 ) return

     where ( linear(1:m) )
        ap(1:m) = ap(1:m) * rho(1:m)
     elsewhere
        ! Just for safety (since these position should have been
        ! untouched in the call to sevalgjacp above).
        ap(1:m) = 0.0d0
     end where
     
     call sevalgjacp(n,x,dum,m,ap,atap,'t',gotj,inform)
     if ( inform .ne. 0 ) return

     hp(1:n) = hp(1:n) + atap(1:n)

     ! Apparently, the code below was wrong (corrected on November 19, 2015)
     
!!$     do i = 1,n
!!$        if ( linear(i) ) then
!!$           ptmp(i) = p(i)
!!$        else
!!$           ptmp(i) = 0.0d0
!!$        end if
!!$     end do
!!$
!!$     gotj = .false.
!!$     call sevalgjacp(n,x,dum,m,ap,ptmp,'j',gotj,inform)
!!$     if ( inform .ne. 0 ) return
!!$
!!$     do j = 1,m
!!$        ap(j) = ap(j) * rho(j)
!!$     end do
!!$
!!$     call sevalgjacp(n,x,dum,m,ap,atap,'t',gotj,inform)
!!$     if ( inform .ne. 0 ) return
!!$
!!$     do i = 1,n
!!$        hp(i) = hp(i) + atap(i)
!!$     end do

  else
     do j = 1,m
        if ( equatn(j) .or. svddpdc(j) .gt. 0.0d0 ) then
           if ( linear(j) ) then

              ! Compute inner product <a,p>
              ajtp = 0.0d0
              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                 ajtp = ajtp + svdjcval(i) * p(svdjcvar(i))
              end do

              ajtp = ajtp * rho(j)

              ! Add rho * ajtp * a
              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                 hp(svdjcvar(i)) = hp(svdjcvar(i)) + ajtp * svdjcval(i)
              end do

           end if
        end if
     end do
  end if

end subroutine ievalhalp

! ******************************************************************
! ******************************************************************

subroutine sevalfeas(n,x,m,equatn,csupn,csupnu,inform)

  use modsvdgrad
  use problvls, only: sc,scale

  implicit none

  ! SCALAR ARGUMENTS
  real(kind=8), intent(inout) :: csupn,csupnu
  integer,      intent(in)    :: m,n
  integer,      intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m)
  real(kind=8), intent(inout) :: x(n)

  ! LOCAL SCALARS
  integer :: i

  ! Subroutine sevalnal is always called by GENCAN with the same point
  ! before calling this subroutine. Then, constraints are computed and
  ! saved.

  ! Compute infeasibility

  csupn = 0.0d0
  do i = 1,m
     if ( equatn(i) ) then
        csupn = max( csupn, abs( svdc(i) ) )
     else
        csupn = max( csupn, svdc(i) )
     end if
  end do

  if ( scale ) then
     csupnu = 0.0d0
     do i = 1,m
        if ( equatn(i) ) then
           csupnu = max( csupnu, abs( svdc(i) / sc(i) ) )
        else
           csupnu = max( csupnu, svdc(i) / sc(i) )
        end if
     end do

  else
     csupnu = csupn
  end if

end subroutine sevalfeas

! ******************************************************************
! ******************************************************************

subroutine stafeas(n,x,l,u,m,c,lambda,equatn,nc,ncsupn,inform)

  use modalgparam
  use modsvdgrad
  use problvls, only: sevalgjacp,sevaljac
  use probgiven, only: fcsubt,gcoded,gjaccoded,gjacpcoded,jaccoded,jcnnzlim

  implicit none

  ! SCALAR ARGUMENTS
  real(kind=8), intent(inout) :: ncsupn
  integer,      intent(in)    :: m,n
  integer,      intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical,      intent(in)    :: equatn(m)
  real(kind=8), intent(inout) :: c(m),l(n),lambda(m),nc(n),u(n),x(n)

  ! LOCAL SCALARS
  logical      :: gotj
  integer      :: i,j
  real(kind=8) :: nci

  ! LOCAL ARRAYS
  logical :: ind(m)
  integer :: ncjcsta(m),ncjclen(m),ncjcvar(jcnnzlim)
  real(kind=8) :: dum(1),ncval(n),p(m),ncjcval(jcnnzlim)

  if ( fcsubt .eq. 1 ) then
     ! This else is the choice when gcoded and jaccoded and also when
     ! first derivatives were not provided by the user. In such case,
     ! sevaljac below will ultimately call subroutines to compute 
     ! derivatives by finite differences.

     ! Those are the gradients for the constraints that were already
     ! computed and saved
     
     nc(1:n) = 0.0d0

!!$     ind(1:m) = .false.
!!$
!!$     do j = 1,m
!!$        if ( equatn(j) .or. c(j) .gt. 0.0d0 ) then
!!$           if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then
!!$              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
!!$                 nc(svdjcvar(i)) = nc(svdjcvar(i)) + c(j) * svdjcval(i)
!!$              end do
!!$           else
!!$              ind(j) = .true.
!!$           end if
!!$        end if
!!$     end do

     ! Those are the gradients of the constraints that are absent and,
     ! therefore, must be computed

     ind(1:m) = .true.
     
     call sevaljac(n,x,m,ind,ncjcsta,ncjclen,ncjcvar,ncjcval,jcnnzlim,inform)
     if ( inform .ne. 0 ) return
     do j = 1,m
      !  if ( ind(j) ) then
           do i = ncjcsta(j),ncjcsta(j) + ncjclen(j) - 1
              ! nc(ncjcvar(i)) = nc(ncjcvar(i)) + c(j) * ncjcval(i)
              ! The line above was replaced by the three lines below
              ! on Nov 9, 2016, since, apparently, it was wrong
              ! :). Note that it was correctly coded in in the cases
              ! above and below.
              if ( equatn(j) .or. c(j) .gt. 0.0d0 ) then
                 nc(ncjcvar(i)) = nc(ncjcvar(i)) + c(j) * ncjcval(i)
              end if
           end do
      !  end if
     end do

  else ! if ( fcsubt .eq. 2 ) then
     if ( gjacpcoded ) then

        do j = 1,m
           if ( equatn(j) .or. c(j) .gt. 0.0d0 ) then
              p(j) = c(j)
           else
              p(j) = 0.0d0
           end if
        end do

        gotj = .false.
        call sevalgjacp(n,x,dum,m,p,nc,'t',gotj,inform)
        if ( inform .ne. 0 ) return

     else
        nc(1:n) = 0.0d0

        do j = 1,m
           if ( equatn(j) .or. c(j) .gt. 0.0d0 ) then
              do i = svdjcsta(j),svdjcsta(j) + svdjclen(j) - 1
                 nc(svdjcvar(i)) = nc(svdjcvar(i)) + c(j) * svdjcval(i)
              end do
           end if
        end do
     end if
  end if

  ncsupn = 0.0d0
  do i = 1,n
     nci = x(i) - nc(i)
     if ( l(i) .le. nci .and. nci .le. u(i) ) then
        nci = - nc(i)
     else
        nci = max( l(i), min( nci, u(i) ) ) - x(i)
     end if

     ncsupn = max( ncsupn, abs( nci ) )
  end do

end subroutine stafeas

! ******************************************************************
! ******************************************************************

subroutine coo2csr(m,nnz,arow,acol,aval,alen,asta)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: m,nnz

  ! ARRAY ARGUMENTS
  integer,      intent(inout) :: acol(nnz),alen(m),arow(nnz),asta(m)
  real(kind=8), intent(inout) :: aval(nnz)

  ! This subroutines converts a matrix from coordinate format to
  ! compressed sparse row format.

  ! LOCAL SCALARS
  integer      :: i,j,col,coltmp,row,rowtmp
  real(kind=8) :: val,valtmp

  ! The if below is wrong, since arrays alen and asta remain un-initialized
  ! if ( nnz .eq. 0 ) return

  alen(1:m) = 0

  do i = 1,nnz
     row = arow(i)
     alen(row) = alen(row) + 1
  end do

  asta(1) = 1
  do i = 2,m
     asta(i) = asta(i-1) + alen(i-1)
  end do

  do i = 1,nnz

     val = aval(i)
     col = acol(i)
     row = arow(i)

     arow(i) = - 1

     do while ( row .ge. 0 )
        j = asta(row)
        asta(row) = j + 1

        valtmp = aval(j)
        coltmp = acol(j)
        rowtmp = arow(j)

        aval(j) = val
        acol(j) = col
        arow(j) = - 1

        val = valtmp
        col = coltmp
        row = rowtmp
     end do

  end do

  asta(1:m) = asta(1:m) - alen(1:m)

end subroutine coo2csr

! ******************************************************************
! ******************************************************************

function sstop(n,x,m,lambda,rho,equatn,linear,inform)

  use modalgparam, only: innercall

  implicit none

  ! FUNCTION TYPE
  logical :: sstop

  ! SCALAR ARGUMENTS
  integer, intent(in)    :: m,n
  integer, intent(inout) :: inform

  ! ARRAY ARGUMENTS
  logical, intent(in)         :: equatn(m),linear(m)
  real(kind=8), intent(in)    :: lambda(m),rho(m)
  real(kind=8), intent(inout) :: x(n)

  ! EXTERNAL FUNCTIONS
  logical :: minsqstop

  if ( innercall ) then
     sstop = minsqstop(n,x,inform)
     return
  end if

end function sstop
