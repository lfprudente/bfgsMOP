module c_algencan_wrapper

  use iso_c_binding, only : c_bool, c_char, c_double, c_int, c_null_char, &
                            c_ptr, c_funptr, c_f_procpointer, c_f_pointer

  implicit none

  save

  ! INTERFACES TO C USER PROVIDED ROUTINES
  abstract interface

     subroutine evalf(n,x,f,flag) bind(C)
       import :: c_double, c_int

       integer(kind=c_int), value :: n
       integer(kind=c_int)        :: flag
       real(kind=c_double)        :: f,x(n)
     end subroutine evalf

     subroutine evalg(n,x,g,flag) bind(C)
       import :: c_double, c_int

       integer(kind=c_int), value :: n
       integer(kind=c_int)        :: flag
       real(kind=c_double)        :: g(n),x(n)
     end subroutine evalg

     subroutine evalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: lim,n
       integer(kind=c_int)         :: flag,hnnz,hcol(lim),hrow(lim)
       logical(kind=c_bool)        :: lmem
       real(kind=c_double)         :: hval(lim),x(n)
     end subroutine evalh

     subroutine evalc(n,x,ind,c,flag) bind(C)
       import :: c_double, c_int

       integer(kind=c_int), value :: ind,n
       integer(kind=c_int)        :: flag
       real(kind=c_double)        :: c,x(n)
     end subroutine evalc

     subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: n,ind,lim
       integer(kind=c_int)         :: flag,jcnnz,jcvar(lim)
       logical(kind=c_bool)        :: lmem
       real(kind=c_double)         :: x(n),jcval(lim)
     end subroutine evaljac

     subroutine evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: n,ind,lim
       integer(kind=c_int)         :: flag,hcnnz,hccol(lim),hcrow(lim)
       logical(kind=c_bool)        :: lmem
       real(kind=c_double)         :: hcval(lim),x(n)
     end subroutine evalhc

     subroutine evalfc(n,x,f,m,c,flag) bind(C)
       import :: c_double, c_int

       integer(kind=c_int), value :: m,n
       integer(kind=c_int)        :: flag
       real(kind=c_double)        :: f,c(m),x(n)
     end subroutine evalfc

     subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: m,n,lim
       integer(kind=c_int)         :: flag,jcnnz,jcfun(lim),jcvar(lim)
       logical(kind=c_bool)        :: lmem
       real(kind=c_double)         :: g(n),jcval(lim),x(n)
     end subroutine evalgjac

     subroutine evalgjacp(n,x,g,m,p,q,work,gotj,flag) bind(C)
       import :: c_bool, c_char, c_double, c_int

       integer(kind=c_int),    value :: m,n
       integer(kind=c_int)           :: flag
       logical(kind=c_bool)          :: gotj
       character(kind=c_char), value :: work
       real(kind=c_double)           :: g(n),p(m),q(n),x(n)
     end subroutine evalgjacp

     subroutine evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: m,n,lim
       integer(kind=c_int)         :: flag,hlnnz,hlcol(lim),hlrow(lim)
       logical(kind=c_bool)        :: lmem
       real(kind=c_double),  value :: sf
       real(kind=c_double)         :: hlval(lim),lambda(m),sc(m),x(n)
     end subroutine evalhl

     subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag) bind(C)
       import :: c_bool, c_double, c_int

       integer(kind=c_int),  value :: m,n
       integer(kind=c_int)         :: flag
       logical(kind=c_bool)        :: goth
       real(kind=c_double),  value :: sf
       real(kind=c_double)         :: hp(n),lambda(m),p(n),sc(m),x(n)
     end subroutine evalhlp

  end interface

  ! GLOBAL PROCEDURE POINTERS TO C USER PROVIDED ROUTINES
  procedure(evalf),     protected, pointer :: c__evalf    
  procedure(evalg),     protected, pointer :: c__evalg    
  procedure(evalh),     protected, pointer :: c__evalh    
  procedure(evalc),     protected, pointer :: c__evalc    
  procedure(evaljac),   protected, pointer :: c__evaljac  
  procedure(evalhc),    protected, pointer :: c__evalhc   
  procedure(evalfc),    protected, pointer :: c__evalfc   
  procedure(evalgjac),  protected, pointer :: c__evalgjac 
  procedure(evalgjacp), protected, pointer :: c__evalgjacp
  procedure(evalhl),    protected, pointer :: c__evalhl   
  procedure(evalhlp),   protected, pointer :: c__evalhlp

contains

  subroutine c_algencan(fsub,gsub,hsub,csub,jacsub,hcsub,fcsub,gjacsub,   &
       gjacpsub,hlsub,hlpsub,jcnnzmax,hnnzmax,epsfeas,epsopt,efstain,     &
       eostain,efacc,eoacc,outputfnm,specfnm,nvparam,vparam,norig,xorig,  &
       lorig,uorig,m,lambda,equatn,linear,coded,checkder,fu,cnormu,snorm, &
       nlpsupn,inform) bind(C, name="c_algencan")

    implicit none

    ! SCALAR ARGUMENTS
    logical(kind=c_bool), value, intent(in)    :: checkder
    integer(kind=c_int),  value, intent(in)    :: hnnzmax,jcnnzmax,m,norig,nvparam
    integer(kind=c_int),         intent(out)   :: inform
    real(kind=c_double),         intent(inout) :: cnormu,efacc,efstain,eoacc,eostain, &
                                                  epsfeas,epsopt,fu,nlpsupn,snorm
    ! ARRAY ARGUMENTS
    type(c_ptr),            intent(in)    :: vparam(nvparam)
    character(kind=c_char), intent(in)    :: specfnm(81),outputfnm(81)
    logical(kind=c_bool),   intent(in)    :: coded(11),linear(m)
    logical(kind=c_bool),   intent(inout) :: equatn(m)
    real(kind=c_double),    intent(in)    :: lorig(norig),uorig(norig)
    real(kind=c_double),    intent(inout) :: lambda(m),xorig(norig)

    ! EXTERNAL SUBROUTINES ARGUMENTS
    type(c_funptr), value :: fsub,gsub,hsub,csub,jacsub,hcsub,fcsub, &
                             gjacsub,gjacpsub,hlsub,hlpsub

    ! LOCAL SCALARS
    integer :: cnt,i,len,start

    ! LOCAL ARRAYS
    character(len=80) :: outputfnm_,specfnm_,vparam_(nvparam)
    logical           :: checkder_,coded_(11),equatn_(m),linear_(m)

    ! LOCAL POINTERS
    character(len=81), pointer :: str_ptr
    !character(kind=c_char), pointer :: str_ptr

    ! Assign procedure pointers from C
    call c_f_procpointer(fsub,     c__evalf    )
    call c_f_procpointer(gsub,     c__evalg    )
    call c_f_procpointer(hsub,     c__evalh    )
    call c_f_procpointer(csub,     c__evalc    )
    call c_f_procpointer(jacsub,   c__evaljac  )
    call c_f_procpointer(hcsub,    c__evalhc   )
    call c_f_procpointer(fcsub,    c__evalfc   )
    call c_f_procpointer(gjacsub,  c__evalgjac )
    call c_f_procpointer(gjacpsub, c__evalgjacp)
    call c_f_procpointer(hlsub,    c__evalhl   )
    call c_f_procpointer(hlpsub,   c__evalhlp  )

    ! Process outputfnm string
    len = 0
    outputfnm_ = ""
    do while( outputfnm(len+1) .ne. c_null_char .and. len .lt. 80 )
       len = len + 1
       outputfnm_(len:len) = outputfnm(len)
    end do

    ! Process specfnm string
    len = 0
    specfnm_ = ""
    do while( specfnm(len+1) .ne. c_null_char .and. len .lt. 80 )
       len = len + 1
       specfnm_(len:len) = specfnm(len)
    end do

    ! Process vparam array
    if ( nvparam .gt. 0 ) then
       do cnt = 1,nvparam
          ! Retrieves the cnt-th param and stores in str_ptr
          call c_f_pointer( vparam(cnt), str_ptr )

          ! Computes the length of the string
          len = 0
          do while( str_ptr(len+1:len+1) .ne. c_null_char .and. len .lt. 80 )
             len = len + 1
          end do

          ! Copies the string to vparam_
          vparam_(cnt) = str_ptr(1:len)
       end do
    end if

    ! Convert logicals
    coded_(1:11) = logical( coded(1:11) )
    equatn_(1:m) = logical( equatn(1:m) )
    linear_(1:m) = logical( linear(1:m) )
    checkder_    = logical( checkder )

    ! Call Algencan
    call algencan(f__evalf,f__evalg,f__evalh,f__evalc,f__evaljac,f__evalhc, &
         f__evalfc,f__evalgjac,f__evalgjacp,f__evalhl,f__evalhlp,jcnnzmax,  &
         hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm_,     &
         specfnm_,nvparam,vparam_,norig,xorig,lorig,uorig,m,lambda,equatn_, &
         linear_,coded_,checkder_,fu,cnormu,snorm,nlpsupn,inform)

  end subroutine c_algencan

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalf(n,x,f,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int), intent(in)  :: n
    integer(kind=c_int), intent(out) :: flag
    real(kind=c_double), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=c_double), intent(in) :: x(n)

    call c__evalf(n,x,f,flag)

  end subroutine f__evalf

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalg(n,x,g,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int), intent(in)  :: n
    integer(kind=c_int), intent(out) :: flag

    ! ARRAY ARGUMENTS
    real(kind=c_double), intent(in)  :: x(n)
    real(kind=c_double), intent(out) :: g(n)

    call c__evalg(n,x,g,flag)
   
  end subroutine f__evalg

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,             intent(out) :: lmem
    integer(kind=c_int), intent(in)  :: lim,n
    integer(kind=c_int), intent(out) :: flag,hnnz

    ! ARRAY ARGUMENTS
    integer(kind=c_int), intent(out) :: hcol(lim),hrow(lim)
    real(kind=c_double), intent(in)  :: x(n)
    real(kind=c_double), intent(out) :: hval(lim)

    ! LOGICAL SCALARS
    logical(kind=c_bool) :: lmem_

    call c__evalh(n,x,hrow,hcol,hval,hnnz,lim,lmem_,flag)

    hrow(1:hnnz) = hrow(1:hnnz) + 1
    hcol(1:hnnz) = hcol(1:hnnz) + 1

    lmem = logical( lmem_ )

  end subroutine f__evalh

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalc(n,x,ind,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int), intent(in)  :: ind,n
    integer(kind=c_int), intent(out) :: flag
    real(kind=c_double), intent(out) :: c

    ! ARRAY ARGUMENTS
    real(kind=c_double), intent(in)  :: x(n)

    call c__evalc(n,x,ind-1,c,flag)

  end subroutine f__evalc

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,             intent(out) :: lmem
    integer(kind=c_int), intent(in)  :: ind,lim,n
    integer(kind=c_int), intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer(kind=c_int), intent(out) :: jcvar(lim)
    real(kind=c_double), intent(in)  :: x(n)
    real(kind=c_double), intent(out) :: jcval(lim)

    ! LOCAL SCALARS
    logical(kind=c_bool) :: lmem_

    call c__evaljac(n,x,ind-1,jcvar,jcval,jcnnz,lim,lmem_,flag)

    jcvar(1:jcnnz) = jcvar(1:jcnnz) + 1

    lmem = logical( lmem_ )

  end subroutine f__evaljac

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,             intent(out) :: lmem
    integer(kind=c_int), intent(in)  :: ind,lim,n
    integer(kind=c_int), intent(out) :: flag,hcnnz

    ! ARRAY ARGUMENTS
    integer(kind=c_int), intent(out) :: hccol(lim),hcrow(lim)
    real(kind=c_double), intent(in)  :: x(n)
    real(kind=c_double), intent(out) :: hcval(lim)

    ! LOCAL SCALARS
    logical(kind=c_bool) :: lmem_

    call c__evalhc(n,x,ind-1,hcrow,hccol,hcval,hcnnz,lim,lmem_,flag)

    hcrow(1:hcnnz) = hcrow(1:hcnnz) + 1
    hccol(1:hcnnz) = hccol(1:hcnnz) + 1

    lmem = logical( lmem_ )

  end subroutine f__evalhc

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalfc(n,x,f,m,c,flag)

    implicit none

    ! SCALAR ARGUMENTS
    integer(kind=c_int), intent(in)  :: m,n
    integer(kind=c_int), intent(out) :: flag
    real(kind=c_double), intent(out) :: f

    ! ARRAY ARGUMENTS
    real(kind=c_double), intent(in)  :: x(n)
    real(kind=c_double), intent(out) :: c(m)

    call c__evalfc(n,x,f,m,c,flag)

  end subroutine f__evalfc

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,             intent(out) :: lmem
    integer(kind=c_int), intent(in)  :: lim,m,n
    integer(kind=c_int), intent(out) :: flag,jcnnz

    ! ARRAY ARGUMENTS
    integer(kind=c_int), intent(out) :: jcfun(lim),jcvar(lim)
    real(kind=c_double), intent(in)  :: x(n)
    real(kind=c_double), intent(out) :: g(n),jcval(lim)

    ! LOCAL SCALARS
    logical(kind=c_bool) :: lmem_

    call c__evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem_,flag)

    jcfun(1:jcnnz) = jcfun(1:jcnnz) + 1
    jcvar(1:jcnnz) = jcvar(1:jcnnz) + 1

    lmem = logical( lmem_ )
    
  end subroutine f__evalgjac

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalgjacp(n,x,g,m,p,q,work,gotj,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,                intent(inout) :: gotj
    integer(kind=c_int),    intent(in)    :: m,n
    integer(kind=c_int),    intent(out)   :: flag
    character,              intent(in)    :: work

    ! ARRAY ARGUMENTS
    real(kind=c_double), intent(in)    :: x(n)
    real(kind=c_double), intent(inout) :: p(m),q(n)
    real(kind=c_double), intent(out)   :: g(n)

    ! LOCAL SCALARS
    logical(kind=c_bool)   :: gotj_
    character(kind=c_char) :: work_

    gotj_ = logical( gotj, c_bool )
    work_ = transfer( work, work_ )

    call c__evalgjacp(n,x,g,m,p,q,work_,gotj_,flag)

    gotj = logical( gotj_ )

  end subroutine f__evalgjacp

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim, &
       lmem,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,             intent(out) :: lmem
    integer(kind=c_int), intent(in)  :: lim,m,n
    integer(kind=c_int), intent(out) :: flag,hlnnz
    real(kind=c_double), intent(in)  :: sf

    ! ARRAY ARGUMENTS
    integer(kind=c_int), intent(out) :: hlcol(lim),hlrow(lim)
    real(kind=c_double), intent(in)  :: lambda(m),sc(m),x(n)
    real(kind=c_double), intent(out) :: hlval(lim)

    ! LOCAL SCALARS
    logical(kind=c_bool) :: lmem_

    call c__evalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim, &
         lmem_,flag)

    hlrow(1:hlnnz) = hlrow(1:hlnnz) + 1
    hlcol(1:hlnnz) = hlcol(1:hlnnz) + 1

    lmem = logical( lmem_ )

  end subroutine f__evalhl

  ! ******************************************************************
  ! ******************************************************************

  subroutine f__evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

    implicit none

    ! SCALAR ARGUMENTS
    logical,             intent(inout) :: goth
    integer(kind=c_int), intent(in)    :: m,n
    integer(kind=c_int), intent(out)   :: flag
    real(kind=c_double), intent(in)    :: sf

    ! ARRAY ARGUMENTS
    real(kind=c_double), intent(in)  :: lambda(m),p(n),sc(m),x(n)
    real(kind=c_double), intent(out) :: hp(n)

    ! LOCAL SCALARS
    logical(kind=c_bool) :: goth_

    goth_ = logical( goth, c_bool )

    call c__evalhlp(n,x,m,lambda,sf,sc,p,hp,goth_,flag)
    
    goth = logical( goth_ )

  end subroutine f__evalhlp

end module c_algencan_wrapper
