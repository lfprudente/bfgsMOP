! ******************************************************************
! ******************************************************************

subroutine setp(n,x)

  use iso_c_binding

  implicit none

  ! C interface
  interface
     subroutine ampl_setp(n,x) bind(C, name="ampl_setp")
       import :: c_double, c_int
       integer(kind=c_int), value :: n
       real(kind=c_double)        :: x(n)
     end subroutine ampl_setp
  end interface

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  call ampl_setp( int(n,c_int), x )

end subroutine setp

! ******************************************************************
! ******************************************************************

subroutine unsetp()

  implicit none

end subroutine unsetp
