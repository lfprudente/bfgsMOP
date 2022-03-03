module modrspace

  implicit none

  save

  ! SCALARS
  integer, public :: nfull

  ! ARRAYS
  integer,      allocatable, public :: ind(:)
  real(kind=8), allocatable, public :: xfull(:)

end module modrspace
