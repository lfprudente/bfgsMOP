module modsydat

  implicit none

  save

  ! SCALARS
  real(kind=8), public :: seucn,sts,sty,yeucn

  ! ARRAYS
  real(kind=8), allocatable, public :: s(:),y(:)

end module modsydat
