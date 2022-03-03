module modhappdat

  implicit none

  save 

  ! SCALARS
  real(kind=8), public :: hlspg,hstds

  ! ARRAYS
  real(kind=8), allocatable, public :: hds(:)

end module modhappdat
