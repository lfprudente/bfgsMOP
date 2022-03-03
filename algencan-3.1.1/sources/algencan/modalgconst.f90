module modalgconst

  implicit none

  ! General constants

  real(kind=8), parameter :: fmin = - 1.0d+20

  ! Line search constants

  integer,      parameter :: maxextrap  =       100
  integer,      parameter :: mininterp  =         4

  real(kind=8), parameter :: gamma      =   1.0d-04
  real(kind=8), parameter :: beta       =     0.5d0
  real(kind=8), parameter :: sigma1     =     0.1d0
  real(kind=8), parameter :: sigma2     =     0.9d0
  real(kind=8), parameter :: etaint     =     2.0d0
  real(kind=8), parameter :: etaext     =     2.0d0

  ! Safeguarding spectral step constants

  real(kind=8), parameter :: lspgma     =   1.0d+10
  real(kind=8), parameter :: lspgmi     =   1.0d-10

  ! Conjugate gradients constants

  real(kind=8), parameter :: theta      =   1.0d-06
  real(kind=8), parameter :: epsnqmp    =   1.0d-08
  integer,      parameter :: maxcgitnp  =         5

  ! BETRA constants

  logical,      parameter :: extrp4     =   .true.
  logical,      parameter :: extrp5     =   .true.

  integer,      parameter :: msmaxit    =       20

  real(kind=8), parameter :: mseps      =  1.0d-08
  real(kind=8), parameter :: mssig      =  1.0d-01
  real(kind=8), parameter :: msrho      =  9.0d-01
  real(kind=8), parameter :: phieps     =  1.0d-08
  real(kind=8), parameter :: trdelini   =  1.0d+02
  real(kind=8), parameter :: trdelmin   =  1.0d-08
  real(kind=8), parameter :: tralpha    =  1.0d-01

  ! ACCELERATION PROCESS constants

  integer,      parameter :: maxaccit   =        10

  ! GENCAN constants

  integer,      parameter :: itnplevel  =         2
  integer,      parameter :: maxinnitnp =         3
  integer,      parameter :: maxinnit   =      1000
  integer,      parameter :: n0         =       500
  integer,      parameter :: n1         =     10000
  integer,      parameter :: n2         =     20000

  real(kind=8), parameter :: delmin     =   1.0d+04
  real(kind=8), parameter :: eta        =   0.1d0  
  real(kind=8), parameter :: cggpnf     =   1.0d-04
  real(kind=8), parameter :: cgepsi     =   1.0d-01
  real(kind=8), parameter :: cgepsf     =   1.0d-08

  ! ALGENCAN constants

  integer,      parameter :: maxoutitnp =        10
  integer,      parameter :: maxoutit   =        50

  real(kind=8), parameter :: lammin     = - 1.0d+20
  real(kind=8), parameter :: lammax     =   1.0d+20

  real(kind=8), parameter :: rhofrac    =   0.5d0  
  real(kind=8), parameter :: rhomult    =   1.0d+01
  real(kind=8), parameter :: rhomax     =   1.0d+20

  logical,      parameter :: rhorestart =    .true.
  real(kind=8), parameter :: rhoinimin  =   1.0d-08
  real(kind=8), parameter :: rhoinimax  =   1.0d+08

end module modalgconst
