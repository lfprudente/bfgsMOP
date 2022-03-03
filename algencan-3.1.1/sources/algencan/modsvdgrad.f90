module modsvdgrad

  implicit none

  save

  ! SCALARS
  logical, public :: svdconstrc,svdgotc

  ! ARRAYS
  integer, allocatable, public      :: svdjcvar(:),svdjcsta(:),svdjclen(:)
  real(kind=8), allocatable, public :: svdc(:),svddpdc(:),svdg(:), &
       svdgparc(:),svdjcval(:)

end module modsvdgrad
