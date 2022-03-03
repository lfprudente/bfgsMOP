!-----------------------------------------------------------------------   
! These routines were extracted in parts from Algencan. 
!-----------------------------------------------------------------------

subroutine checkdF(n,m,xini,l,u)

	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in)    :: m,n

	! ARRAY ARGUMENTS
	real(kind=8), intent(in) :: xini(n),l(n),u(n)

	! This subrotutine checks the user supplied first derivatives 
	! subroutines for computing the objective functions gradients 

	! LOCAL SCALARS
	character    :: answer
	integer      :: i
	real(kind=8) :: seed,eps

	! LOCAL ARRAYS
	real(kind=8) :: x(n)

	! FUNCTIONS
	real(kind=8) :: drand

	! SET A PERTURBATIONS OF X
	
	eps = 1.0d-8

	seed = 123456.0d0
	do i = 1,n
		 if ( l(i) .lt. xini(i) .and. xini(i) .lt. u(i) ) then
				x(i) = xini(i) + eps * ( 2.0d0 * drand(seed) - 1.0d0 ) * &
						 max( 1.0d0, abs( xini(i) ) )
		 else if ( xini(i) == l(i) ) then
				x(i) = xini(i) + eps * drand(seed) * max( 1.0d0, abs( xini(i) ) )
		 else
				x(i) = xini(i) - eps * drand(seed) * max( 1.0d0, abs( xini(i) ) )
		 end if
		 x(i) = max( l(i), min( x(i), u(i) ) )
	end do

	write(*,1000)

	do i = 1,n
		 write(* ,1010) i,x(i)
	end do

	! CHECK OBJECTIVE FUNCTIONS GRADIENTS

	do i = 1,m

		write(*,1020) i

		read(*,*) answer

		if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
			 return

		else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
			 cycle

		else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
			 call checkgF(n,x,i)
		end if

	end do
	
	! CHECK HESSIANS OF THE OBJECTIVE FUNCTIONS
	
!	do i = 1,m

!		write(* ,1030) i
		
!		read(*,*) answer

!		if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
!			return

!		else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
!			cycle
!		else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
!			call checkhF(n,x,i)
!		end if
		
!	end do


    ! NON-EXECUTABLE STATEMENTS

1000 format(/,1X,'Derivatives will be tested at the perturbed ', &
         'initial guess: ')
1010 format(  1X,'x(',I6,') = ',1P,D15.8)

1020 format(/,1X,'Check gradient of function ',I5,'?',   &
         /,1X,'Type Y(es), A(bort checking) or ', &
         'S(kip): ')
1030 format(/,1X,'Check Hessian of function ',I5,'?',   &
         /,1X,'Type Y(es), A(bort checking) or ', &
         'S(kip): ')     

end subroutine checkdF

! *****************************************************************
! *****************************************************************

subroutine checkgF(n,x,ind)

	use myproblem

	implicit none

	! SCALAR ARGUMENTS
	integer, intent(in)    :: n,ind

	! ARRAY ARGUMENTS
	real(kind=8), intent(inout) :: x(n)

	! This subrotutine checks the user supplied subroutine evalg for
	! computing the gradient of the objective function using central
	! finite differences with two different discretization steps.

	! LOCAL SCALARS
	integer      :: i
	real(kind=8) :: eps,fminus,fplus,gdiff1,gdiff2,maxerr,step1,step2,tmp

	! LOCAL ARRAYS
	real(kind=8) :: g(n)
	
	call evalg(n,x,g,ind)

	write(*,1000)

	maxerr = 0.0d0
	
	eps = 1.0d-5

	do i = 1,n
		 tmp  = x(i)

		 step1 = eps * max( abs( tmp ), 1.0d0 )

		 x(i) = tmp + step1
		 call evalf(n,x,fplus,ind)
		 
		 x(i) = tmp - step1	 
		 call evalf(n,x,fminus,ind)

		 gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

		 step2 = eps * max( abs( tmp ), 1.0d-03 )

		 x(i) = tmp + step2
		 call evalf(n,x,fplus,ind)

		 x(i) = tmp - step2
		 call evalf(n,x,fminus,ind)
		 
		 x(i) = tmp

		 gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

		 tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

		 write(*,1010) i,g(i),gdiff1,gdiff2,tmp

		 maxerr = max( maxerr, tmp )

	end do

	write(*,1020) maxerr

	! NON-EXECUTABLE STATEMENTS

1000 format(/,1X,'Gradient vector of the objective function.',/,1X, &
	             'Index',13X,'evalg',2X,'Central diff (two different steps)',4X,&
	             'Absolute error')
1010 format(  1X,I5,4(3X,1P,D15.8))
1020 format(  1X,'Maximum absolute error = ',1P,D15.8)

end subroutine checkgF

! *****************************************************************
! *****************************************************************

!subroutine checkhF(n,x,ind)

!	use myproblem

!	implicit none

!	! SCALAR ARGUMENTS
!	integer, intent(in) :: n,ind

!	! ARRAY ARGUMENTS
!	real(kind=8), intent(inout) :: x(n)

!	! This subrotutine checks the user supplied subroutine evalh for
!	! computing the Hessian of the objective function using central
!	! finite differences with two different discretization steps.

!	! LOCAL SCALARS
!	logical :: nullcol
!	integer :: i,j
!	real(kind=8) :: elem,hdiff1,hdiff2,maxerr,step1,step2,tmp,eps

!	! LOCAL ARRAYS
!	real(kind=8) :: g(n),gplus1(n),gplus2(n),maxcoe(n)
!	real(kind=8) :: h(n,n)

!	eps = 1.0d-8

!	! Compute the gradient of the objective function at x

!	call evalg(n,x,g,ind)

!	! Compute the Hessian of the objective function at x and save in a
!	! dense matrix

!	call evalh(n,x,h,ind)

!	! Test column by column

!	write(* ,200)

!	maxerr = 0.0d0

!	do j = 1,n

!	   tmp  = x(j)

!	   step1 = eps * max( abs( tmp ), 1.0d0 )

!	   x(j) = tmp + step1       
!	   call evalg(n,x,gplus1,ind)

!	   step2 = eps * max( abs( tmp ), 1.0d-03 )

!	   x(j) = tmp + step2      
!	   call evalg(n,x,gplus2,ind)

!	   x(j) = tmp

!	   write(* ,210) j

!	   maxcoe(j) = 0.0d0

!	   nullcol = .true.

!	   do i = 1,n
!		  if ( i .ge. j ) then
!			 elem = H(i,j)
!		  else
!			 elem = H(j,i)
!		  end if
!		  hdiff1 = ( gplus1(i) - g(i) ) / step1
!		  hdiff2 = ( gplus2(i) - g(i) ) / step2
!		  tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
!		  if ( elem   .ne. 0.0d0 .or. &
!			   hdiff1 .ne. 0.0d0 .or. &
!			   hdiff2 .ne. 0.0d0 ) then
!			 if ( nullcol ) then
!				nullcol = .false.
!				write(* ,220)
!				write(10,220)
!			 end if
!			 write(* ,230) i,elem,hdiff1,hdiff2,tmp
!			 write(10,230) i,elem,hdiff1,hdiff2,tmp
!		  end if
!		  maxcoe(j) = max( maxcoe(j), tmp )
!	   end do

!	   maxerr = max( maxerr, maxcoe(j) )

!	   if ( nullcol ) then
!		  write(* ,240)
!		  write(10,240)
!	   else
!		  write(* ,250) maxcoe(j)
!		  write(10,250) maxcoe(j)
!	   end if

!	end do

!	write(* ,*)

!	do j = 1,n
!	   write(* ,260) j,maxcoe(j)
!	end do

!	write(* ,270) maxerr


!! NON-EXECUTABLE STATEMENTS

!200 format(/,1X,'Hessian matrix of the objective function column by ', &
!	 'column.')
!210 format(/,1X,'Column:  ',I6)
!220 format(/,1X,'Index',13X,'evalh',3X,'Incr. Quoc. (two different ', &
!	 'steps)',4X,'Absolute error')
!230 format(  1X,I5,4(3X,1P,D15.8))
!240 format(  1X,'All the elements of this column are null.')
!250 format(  1X,'Maximum absolute error = ',1P,D15.8)
!260 format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
!270 format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

!end subroutine checkhF


! ******************************************************************
! ******************************************************************

function drand(ix)

  implicit none

  ! This is the random number generator of Schrage:
  !
  ! L. Schrage, A more portable Fortran random number generator, ACM
  ! Transactions on Mathematical Software 5 (1979), 132-138.

  ! FUNCTION TYPE
  real(kind=8) :: drand

  ! SCALAR ARGUMENT
  real(kind=8), intent(inout) :: ix

  ! LOCAL ARRAYS
  real(kind=8) :: a,p,b15,b16,xhi,xalo,leftlo,fhi,k

  data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

  xhi= ix/b16
  xhi= xhi - dmod(xhi,1.d0)
  xalo= (ix-xhi*b16)*a
  leftlo= xalo/b16
  leftlo= leftlo - dmod(leftlo,1.d0)
  fhi= xhi*a + leftlo
  k= fhi/b15
  k= k - dmod(k,1.d0)
  ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
  if (ix.lt.0) ix= ix + p
  drand= ix*4.656612875d-10

  return

end function drand
