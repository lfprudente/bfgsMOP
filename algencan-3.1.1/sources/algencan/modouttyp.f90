module modouttyp
  
  implicit none

  save

  ! SCALARS
  logical, protected :: iprintctl(5)
  integer, protected :: iprint,iprintinn,iprintout,mprint,ncomp,nprint

  ! ARRAYS
  character(len=80), protected :: solfnm

  ! SUBROUTINES
  public :: iniouttyp,setouttyp

contains

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine iniouttyp()

    iprint = 10
    iprintout = iprint / 10
    iprintinn = mod( iprint, 10 )

    ncomp  =  6
    solfnm = ''

    iprintctl(1) = .false.  ! Banner
    iprintctl(2) = .false.  ! Parameters and problem processing
    iprintctl(3) = .false.  ! Warnings and errors messages
    iprintctl(4) = .false.  ! User-provided subs calls counters and timing
    iprintctl(5) = .false. ! Statistics files with a single table line

    open(10,err=100,file='.silent',status='old')
    close(10)

    iprint = 0
    iprintout = 0
    iprintinn = 0

    iprintctl(1:4) = .false.

100 continue

  end subroutine iniouttyp

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine setouttyp(val_iprint,val_ncomp,val_solfnm,val_iprintctl5, &
       val_mprint,val_nprint)

    implicit none

    logical, intent(in), optional :: val_iprintctl5
    integer, intent(in), optional :: val_iprint,val_ncomp,val_mprint, &
         val_nprint
    character(len=80), intent(in), optional :: val_solfnm

    if ( present( val_iprint ) ) then
       iprint = val_iprint
       iprintout = iprint / 10
       iprintinn = mod( iprint, 10 )
    end if

    if ( present( val_ncomp ) ) ncomp = val_ncomp

    if ( present( val_mprint ) ) mprint = val_mprint

    if ( present( val_nprint ) ) nprint = val_nprint

    if ( present( val_solfnm ) ) solfnm = val_solfnm

    if ( present( val_iprintctl5 ) ) iprintctl(5) = val_iprintctl5

  end subroutine setouttyp

  ! ******************************************************************
  ! ******************************************************************
  
end module modouttyp
