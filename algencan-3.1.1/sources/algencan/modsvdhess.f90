module modsvdhess

  implicit none

  save

  ! SCALARS
  integer, public :: svdhnnz

  ! ARRAYS
  integer,      allocatable, public :: svdhcol(:),svdhrow(:)
  real(kind=8), allocatable, public :: svdhval(:)

end module modsvdhess
