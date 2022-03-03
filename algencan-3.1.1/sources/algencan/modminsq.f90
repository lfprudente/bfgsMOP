module modminsq

  implicit none

  save

  ! SCALARS
  integer, public :: annz,ncols,nrows

  ! ARRAYS
  integer, allocatable, public :: arow(:),acol(:)
  real(kind=8), allocatable, public :: aval(:),b(:)

end module modminsq
