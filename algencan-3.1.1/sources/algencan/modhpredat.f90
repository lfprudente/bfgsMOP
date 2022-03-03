module modhpredat

  implicit none

  save

  ! SCALARS
  real(kind=8), public :: plspg,psmdyty

  ! ARRAYS
  real(kind=8), allocatable, public :: pdiag(:),psmdy(:)

end module modhpredat
