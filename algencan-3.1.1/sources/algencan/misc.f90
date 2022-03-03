! ******************************************************************
! ******************************************************************

logical function IsANumber(x)

  use modmachconst, only: bignum

  implicit none

  ! SCALAR ARGUMENTS
  real(kind=8), intent(in) :: x

  IsANumber = .true.
  if ( .not. abs( x ) .le. bignum ) IsANumber = .false.

end function IsANumber

