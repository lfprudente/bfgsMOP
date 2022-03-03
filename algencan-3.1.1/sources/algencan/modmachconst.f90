module modmachconst

  implicit none

  ! SCALARS
  real(kind=8), protected :: bignum,macheps,macheps12,macheps13,macheps23

  ! SUBROUTINES
  public :: inimachconst

contains

  subroutine inimachconst()

    implicit none

    bignum    = 1.0d+99
    macheps   = 1.0d-16
    macheps12 = sqrt( macheps )
    macheps13 = macheps ** ( 1.0d0 / 3.0d0 )
    macheps23 = macheps ** ( 2.0d0 / 3.0d0 )

  end subroutine inimachconst

end module modmachconst
