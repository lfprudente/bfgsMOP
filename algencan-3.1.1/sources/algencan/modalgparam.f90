module modalgparam

  implicit none

  ! SCALARS
  character(len=6), protected :: hptype
  character(len=4), protected :: precond,lsssubACC,lsssubNW,lsssubTR, &
       sclsubACC,sclsubNW,sclsubTR
  character(len=2), protected :: innslvr
  logical, protected :: innercall,ignoref,skipacc,useustp
  integer, protected :: maxaccitgiven,maxinnitgiven,maxoutitgiven
  real (kind=8), protected :: rhoinigiven,rhomaxgiven

  ! SUBROUTINES
  public :: setalgparam

  contains

  ! ******************************************************************
  ! ******************************************************************
  
  subroutine setalgparam(val_hptype,val_precond,val_lsssubACC,val_lsssubNW, &
       val_lsssubTR,val_sclsubACC,val_sclsubNW,val_sclsubTR,val_innslvr, &
       val_innercall,val_ignoref,val_skipacc,val_useustp,val_maxaccitgiven, &
       val_maxinnitgiven,val_maxoutitgiven,val_rhoinigiven,val_rhomaxgiven)

    implicit none

    character(len=6), intent(in), optional :: val_hptype
    character(len=4), intent(in), optional :: val_precond,val_lsssubACC, &
         val_lsssubNW,val_lsssubTR,val_sclsubACC,val_sclsubNW,val_sclsubTR
    character(len=2), intent(in), optional :: val_innslvr
    logical, intent(in), optional :: val_innercall,val_ignoref,val_skipacc, &
         val_useustp
    integer, intent(in), optional :: val_maxaccitgiven,val_maxinnitgiven, &
         val_maxoutitgiven
    real (kind=8), intent(in), optional :: val_rhoinigiven,val_rhomaxgiven

    if ( present( val_hptype        ) ) hptype        = val_hptype       
    if ( present( val_precond       ) ) precond       = val_precond      
    if ( present( val_lsssubACC     ) ) lsssubACC     = val_lsssubACC    
    if ( present( val_lsssubNW      ) ) lsssubNW      = val_lsssubNW     
    if ( present( val_lsssubTR      ) ) lsssubTR      = val_lsssubTR     
    if ( present( val_sclsubACC     ) ) sclsubACC     = val_sclsubACC    
    if ( present( val_sclsubNW      ) ) sclsubNW      = val_sclsubNW     
    if ( present( val_sclsubTR      ) ) sclsubTR      = val_sclsubTR     
    if ( present( val_innslvr       ) ) innslvr       = val_innslvr      
    if ( present( val_innercall     ) ) innercall     = val_innercall    
    if ( present( val_ignoref       ) ) ignoref       = val_ignoref      
    if ( present( val_skipacc       ) ) skipacc       = val_skipacc      
    if ( present( val_useustp       ) ) useustp       = val_useustp      
    if ( present( val_maxaccitgiven ) ) maxaccitgiven = val_maxaccitgiven
    if ( present( val_maxinnitgiven ) ) maxinnitgiven = val_maxinnitgiven
    if ( present( val_maxoutitgiven ) ) maxoutitgiven = val_maxoutitgiven
    if ( present( val_rhoinigiven   ) ) rhoinigiven   = val_rhoinigiven  
    if ( present( val_rhomaxgiven   ) ) rhomaxgiven   = val_rhomaxgiven  

  end subroutine setalgparam

  ! ******************************************************************
  ! ******************************************************************
  
end module modalgparam
