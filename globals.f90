module globals

implicit none

character(len=10) :: problem
integer :: nfinner
real(kind=8) :: seed
real(kind=8),allocatable :: xinner(:),d(:),sF(:),JF(:,:),JFin(:,:),B(:,:,:)
logical,allocatable :: quadratic(:)

end module globals
