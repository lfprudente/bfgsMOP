! ******************************************************************
! ******************************************************************

subroutine procpar(epsfeas,epsopt,efstain,eostain,efacc,eoacc, &
     outputfnm,specfnm,nvparam,vparam)

  use lssma57
  use lssma86
  use lssma97
  use probgiven,   only: firstde,seconde,truehpr
  use modalgparam, only: hptype,innslvr,lsssubTR,lsssubNW,lsssubACC,&
                         sclsubTR,sclsubNW,sclsubACC,skipacc
  use modouttyp,   only: iprint,iprintctl,ncomp,solfnm
  use problvlu,    only: rmfixv
  use problvlt,    only: slacks
  use problvls,    only: scale

  implicit none

  ! SCALAR ARGUMENTS
  integer,      intent(in)    :: nvparam
  real(kind=8), intent(inout) :: efacc,efstain,eoacc,eostain, &
                                 epsfeas,epsopt

  ! ARRAY ARGUMENTS
  character(len=80), intent(in) :: specfnm,outputfnm,vparam(nvparam)

  ! LOCAL SCALARS
  logical :: scl
  integer :: i,hslsubslen

  ! LOCAL ARRAYS
  character(len=15) :: hslsubs
  character(len=80) :: line

  ! ==================================================================
  ! BANNER
  ! ==================================================================

  if ( iprintctl(1) ) then
     write(*, 8000)
     write(10,8000)
  end if

  ! ==================================================================
  ! PROCESS ARRAY OF PARAMETERS
  ! ==================================================================

  if ( iprintctl(2) ) then
     if ( nvparam .eq. 0 ) then
        write(*, 9014)
        write(10,9014)
     else
        write(*, 9015) nvparam
        write(10,9015) nvparam
     end if
  end if

  do i = 1,nvparam
     call procline(vparam(i),epsfeas,epsopt,efstain,eostain,efacc,eoacc)
  end do

  ! ==================================================================
  ! PROCESS SPECIFICATION FILE
  ! ==================================================================

  if ( specfnm .eq. '' ) then
     if ( iprintctl(2) ) then
        write(*, 8900)
        write(10,8900)
     end if
     go to 500
  end if

  ! OPENING THE SPECIFICATION FILE
  open(20,err=300,file=specfnm,status='old')

  if ( iprintctl(2) ) then
     write(*, 9000)
     write(10,9000)
  end if

  ! MAIN LOOP

100 continue

  ! READING LINES
  read(20,fmt=1000,err=400,end=200) line

  ! PROCESS LINES

  call procline(line,epsfeas,epsopt,efstain,eostain,efacc,eoacc)

  ! IIERATE
  go to 100

  ! END OF LOOP

  ! TERMINATE PROCESSING THE SPECFICATION FILE

  ! CLOSING SPECIFICATION FILE
200 continue
  close(20)
  go to 500

  ! NO SPECIFICATION FILE
300 continue
  if ( iprintctl(2) ) then
     write(*, 9005)
     write(10,9005)
  end if
  go to 500

  ! ERROR READING THE SPECIFICATION FILE
400 continue
  if ( iprintctl(2) ) then
     write(*, 9010)
     write(10,9010)
  end if
  go to 500

  ! ==================================================================
  ! PRINTING PARAMETERS VALUES
  ! ==================================================================

500 continue

  if ( iprintctl(2) ) then

     hslsubs = ''
     hslsubslen = 0
     if ( lss_ma57() ) then
        hslsubs(hslsubslen+1:hslsubslen+5) = 'MA57 '
        hslsubslen = hslsubslen + 5
     end if
     if ( lss_ma86() ) then
        hslsubs(hslsubslen+1:hslsubslen+5) = 'MA86 '
        hslsubslen = hslsubslen + 5
     end if
     if ( lss_ma97() ) then
        hslsubs(hslsubslen+1:hslsubslen+5) = 'MA97 '
        hslsubslen = hslsubslen + 5
     end if
     if ( hslsubslen .eq. 0 ) then
        hslsubs = 'NONE'
     end if

     write(* ,4000) hslsubs,firstde,seconde,truehpr,hptype,        &
                    lsssubTR,sclsubTR,lsssubNW,sclsubNW,lsssubACC, &
                    sclsubACC,innslvr,.not.skipacc,rmfixv,slacks,  &
                    scale,epsfeas,epsopt,efstain,eostain,efacc,    &
                    eoacc,iprint,ncomp,                            &
                    "'"//trim(adjustl(specfnm))//"'",              &
                    "'"//trim(adjustl(outputfnm))//"'",            &
                    "'"//trim(adjustl(solfnm))//"'"
     write(10,4000) hslsubs,firstde,seconde,truehpr,hptype,        &
                    lsssubTR,sclsubTR,lsssubNW,sclsubNW,lsssubACC, &
                    sclsubACC,innslvr,.not.skipacc,rmfixv,slacks,  &
                    scale,epsfeas,epsopt,efstain,eostain,efacc,    &
                    eoacc,iprint,ncomp,                            &
                    "'"//trim(adjustl(specfnm))//"'",              &
                    "'"//trim(adjustl(outputfnm))//"'",            &
                    "'"//trim(adjustl(solfnm))//"'"
  end if

! NON-EXECUTABLE STATEMENTS

 1000 format(A80)
 4000 format(/,1X,'Available HSL subroutines = ',       A15, &
            //,1X,'ALGENCAN PARAMETERS:',/,                  &
             /,1X,'firstde                = ',       19X,L1, &
             /,1X,'seconde                = ',       19X,L1, &
             /,1X,'truehpr                = ',       19X,L1, &
             /,1X,'hptype in TN           = ',       14X,A6, &
             /,1X,'lsslvr in TR           = ',11X,A4,'/',A4, &
             /,1X,'lsslvr in NW           = ',11X,A4,'/',A4, &
             /,1X,'lsslvr in ACCPROC      = ',11X,A4,'/',A4, &
             /,1X,'innslvr                = ',       18X,A2, &
             /,1X,'accproc                = ',       19X,L1, &
             /,1X,'rmfixv                 = ',       19X,L1, &
             /,1X,'slacks                 = ',       19X,L1, &
             /,1X,'scale                  = ',       19X,L1, &
             /,1X,'epsfeas                = ',  8X,1P,D12.4, &
             /,1X,'epsopt                 = ',  8X,1P,D12.4, &
             /,1X,'efstain                = ',  8X,1P,D12.4, &
             /,1X,'eostain                = ',  8X,1P,D12.4, &
             /,1X,'efacc                  = ',  8X,1P,D12.4, &
             /,1X,'eoacc                  = ',  8X,1P,D12.4, &
             /,1X,'iprint                 = ',          I20, &
             /,1X,'ncomp                  = ',          I20, &
            //,1X,'Specification filename = ',          A20, &
             /,1X,'Output filename        = ',          A20, &
             /,1X,'Solution filename      = ',          A20)
 8000 format(/,1X,78('='),                                               &
             /,1X,'This is ALGENCAN 3.1.1.',                             &
             /,1X,'ALGENCAN, an Augmented Lagrangian method for ',       &
                  'nonlinear programming, is part of',/,1X,'the TANGO ', &
                  'Project: Trustable Algorithms for Nonlinear ',        &
                  'General Optimization.',/,1X,'See ',                   &
                  'http://www.ime.usp.br/~egbirgin/tango/ for details.', &
             /,1X,78('='))
 8900 format(/,1X,'The specification file is not being used.')
 9000 format(/,1X,'The specification file was found.',/, &
             /,1X,'Processing specification file:',/)
 9005 format(/,1X,'The specification file was not found in the current folder.')
 9010 format(/,1X,'Error reading the specification file.')
 9014 format(/,1X,'There are no strings to be processed in the array of parameters.')
 9015 format(/,1X,'Processing array of parameters with ',I3,' entrances:',/)

end subroutine procpar

! ******************************************************************
! ******************************************************************

subroutine procline(line,epsfeas,epsopt,efstain,eostain,efacc,eoacc)
  
  use lssma57
  use lssma86
  use lssma97
  use modalgparam, only: setalgparam,hptype,ignoref,innslvr,lsssubTR, &
       lsssubNW,lsssubACC,sclsubTR,sclsubNW,sclsubACC,skipacc, &
       maxoutitgiven,maxinnitgiven,maxaccitgiven,rhoinigiven, &
       rhomaxgiven
  use modouttyp,   only: setouttyp,iniouttyp,iprint,iprintctl,iprintinn,iprintout,mprint,ncomp,nprint
  use probgiven,   only: firstde,seconde,truehpr
  use problvlu,    only: setproblvlu
  use problvlt,    only: setproblvlt
  use problvls,    only: setproblvls

  implicit none

  ! SCALAR ARGUMENTS
  real(kind=8), intent(inout) :: efacc,efstain,eoacc,eostain, &
                                 epsfeas,epsopt

  ! ARRAY ARGUMENTS
  character(len=80), intent(in) :: line

  ! PARAMETERS
  integer, parameter :: nwords = 27

  ! LOCAL SCALARS
  logical      :: lss,scl
  integer      :: i,ifirst,ikey,ilast,inum,j,ff,ll
  real(kind=8) :: dnum

  ! LOCAL ARRAYS
  character(len= 4) :: lsssub,sclsub
  character(len=80) :: keyword,lline,str
  
  ! FUNCTIONS
  logical :: GetWord

  ! DATA BLOCKS
  character(len= 1) :: addinfo(nwords),lower(26),upper(26)
  character(len=45) :: dictionary(nwords)

  data lower /'a','b','c','d','e','f','g','h','i','j','k','l','m', &
       'n','o','p','q','r','s','t','u','v','w','x','y','z'/
  data upper /'A','B','C','D','E','F','G','H','I','J','K','L','M', &
       'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

  data dictionary( 1) /'SKIP-ACCELERATION-PROCESS                    '/
  data dictionary( 2) /'LINEAR-SYSTEMS-SOLVER-IN-ACCELERATION-PROCESS'/
  data dictionary( 3) /'TRUST-REGIONS-INNER-SOLVER                   '/
  data dictionary( 4) /'LINEAR-SYSTEMS-SOLVER-IN-TRUST-REGIONS       '/
  data dictionary( 5) /'NEWTON-LINE-SEARCH-INNER-SOLVER              '/
  data dictionary( 6) /'LINEAR-SYSTEMS-SOLVER-IN-NEWTON-LINE-SEARCH  '/
  data dictionary( 7) /'TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER    '/
  data dictionary( 8) /'MATRIX-VECTOR-PRODUCT-IN-TRUNCATED-NEWTON-LS '/
  data dictionary( 9) /'FIXED-VARIABLES-REMOVAL-AVOIDED              '/
  data dictionary(10) /'ADD-SLACKS                                   '/
  data dictionary(11) /'OBJECTIVE-AND-CONSTRAINTS-SCALING-AVOIDED    '/
  data dictionary(12) /'IGNORE-OBJECTIVE-FUNCTION                    '/
  data dictionary(13) /'FEASIBILITY-TOLERANCE                        '/
  data dictionary(14) /'OPTIMALITY-TOLERANCE                         '/
  data dictionary(15) /'STAINF-FEASIBILITY-TOLERANCE                 '/
  data dictionary(16) /'STAINF-OPTIMALITY-TOLERANCE                  '/
  data dictionary(17) /'ACC-FEASIBILITY-THRESHOLD                    '/
  data dictionary(18) /'ACC-OPTIMALITY-THRESHOLD                     '/
  data dictionary(19) /'ITERATIONS-OUTPUT-DETAIL                     '/
  data dictionary(20) /'NUMBER-OF-ARRAYS-COMPONENTS-IN-OUTPUT        '/
  data dictionary(21) /'SOLUTION-FILENAME                            '/
  data dictionary(22) /'SAVE-STATISTICS-TABLE-LINE                   '/
  data dictionary(23) /'ACCELERATION-PROCESS-ITERATIONS-LIMIT        '/
  data dictionary(24) /'INNER-ITERATIONS-LIMIT                       '/
  data dictionary(25) /'OUTER-ITERATIONS-LIMIT                       '/
  data dictionary(26) /'PENALTY-PARAMETER-INITIAL-VALUE              '/
  data dictionary(27) /'LARGEST-PENALTY-PARAMETER-ALLOWED            '/

  ! Upper-case letter means mandatory extra value. Lowe-case letter means optional

  data addinfo( 1) /' '/
  data addinfo( 2) /'S'/
  data addinfo( 3) /'s'/
  data addinfo( 4) /'S'/
  data addinfo( 5) /'s'/
  data addinfo( 6) /'S'/
  data addinfo( 7) /'s'/
  data addinfo( 8) /'S'/
  data addinfo( 9) /' '/
  data addinfo(10) /' '/
  data addinfo(11) /' '/
  data addinfo(12) /' '/
  data addinfo(13) /'D'/
  data addinfo(14) /'D'/
  data addinfo(15) /'D'/
  data addinfo(16) /'D'/
  data addinfo(17) /'D'/
  data addinfo(18) /'D'/
  data addinfo(19) /'I'/
  data addinfo(20) /'I'/
  data addinfo(21) /'S'/
  data addinfo(22) /' '/
  data addinfo(23) /'I'/
  data addinfo(24) /'I'/
  data addinfo(25) /'I'/
  data addinfo(26) /'D'/
  data addinfo(27) /'D'/

  ! PROCESS LINE

  lline(1:80) = line(1:80)

  ! Convert line to upper case
  do i = 1,80
     do j = 1,26
        if ( lline(i:i) .eq. lower(j) ) then
           lline(i:i) = upper(j)
        end if
     end do
  end do

  ! Get first word skipping blank lines
  if ( .not. GetWord(lline,1,ifirst,ilast) ) then
     go to 200
  end if

  ! Skip comments
  if ( lline(ifirst:ifirst) .eq. '*' .or. &
       lline(ifirst:ifirst) .eq. '#' ) then
     go to 200
  end if

  ! Set keyword
  keyword = lline(ifirst:ilast)

  ! Look up the keyword in the dictionary
  i = 1
130 if ( i .le. nwords ) then
     if ( keyword .ne. dictionary(i) ) then
       i = i + 1
       go to 130
     end if
  end if

  ! Ignore unknown keywords
  if ( i .gt. nwords ) then
     if ( iprintctl(2) ) then
        write(*, 9020) "'"//lline(ifirst:ilast)//"'"
        write(10,9020) "'"//lline(ifirst:ilast)//"'"
     end if
     go to 200
  end if

  ikey = i

  ! Read additional information if needed
  if ( addinfo(ikey) .ne. ' ' ) then

     ! Ignore keywords without the required additional information
     if ( len(trim(adjustl(lline(ilast+1:80)))) .eq. 0 .and. &
          'A' .le. addinfo(ikey) .and. addinfo(ikey) .le. 'Z' ) then
        if ( iprintctl(2) ) then
           write(*, 9030) dictionary(ikey)
           write(10,9030) dictionary(ikey)
        end if
        go to 200
     end if

     ! Convert additional information to the corresponding type
     if ( addinfo(ikey) .eq. 'I' .or. addinfo(ikey) .eq. 'i' ) then
        read(unit=lline(ilast+1:80),fmt=2000) inum

     else if ( addinfo(ikey) .eq. 'D' .or. addinfo(ikey) .eq. 'd' ) then
        read(unit=lline(ilast+1:80),fmt=3000) dnum

     else if ( addinfo(ikey) .eq. 'S' .or. addinfo(ikey) .eq. 's' ) then
!!$        str = ''
!!$        str = trim(adjustl(lline(ilast+1:80)))
        str = adjustl(lline(ilast+1:80))
     end if

  end if

  ! Process keyword
  if ( iprintctl(2) ) then
     if ( addinfo(ikey) .eq. ' ' ) then
        write(*, 9040) dictionary(ikey)
        write(10,9040) dictionary(ikey)
     else if ( addinfo(ikey) .eq. 'I' .or. addinfo(ikey) .eq. 'i' ) then
        write(*, 9041) dictionary(ikey),inum
        write(10,9041) dictionary(ikey),inum
     else if ( addinfo(ikey) .eq. 'D' .or. addinfo(ikey) .eq. 'd' ) then
        write(*, 9042) dictionary(ikey),dnum
        write(10,9042) dictionary(ikey),dnum
     else if ( addinfo(ikey) .eq. 'S' .or. addinfo(ikey) .eq. 's' ) then
        write(*, 9043) dictionary(ikey),"'"//trim(str)//"'"
        write(10,9043) dictionary(ikey),"'"//trim(str)//"'"
     end if
  end if

  ! Set the corresponding algencan argument
  if ( ikey .eq. 1 ) then
     call setalgparam(val_skipacc = .true.)

  else if ( ikey .eq. 2 ) then
     if ( seconde .and. lsssubACC .ne. 'NONE' ) then
        lsssub = ''
        sclsub = ''
        if ( GetWord(str,1,ff,ll) ) then
           lsssub = str(ff:ll)
           if ( GetWord(str,ll+1,ff,ll) ) sclsub = str(ff:ll)
        end if
        if ( lsssub .eq. 'MA57' .and. lss_ma57() ) then
           call setalgparam(val_lsssubACC = lsssub)
           if ( sclsub .eq. 'MC64' .or. sclsub .eq. 'NONE' ) then
              call setalgparam(val_sclsubACC = sclsub)
           else 
              call setalgparam(val_sclsubACC = 'MC64')
              if ( sclsub .ne. '' ) then
                 if ( iprintctl(2) ) then
                    write(* ,9130) dictionary(ikey)
                    write(10,9130) dictionary(ikey)
                 end if
              end if
           end if
        else if ( lsssub .eq. 'MA86' .and. lss_ma86() ) then
           call setalgparam(val_lsssubACC = lsssub)
           if ( sclsub .eq. 'MC64' .or. sclsub .eq.'MC77' .or. &
                sclsub .eq. 'NONE' ) then
              call setalgparam(val_sclsubACC = sclsub)
           else
              call setalgparam(val_sclsubACC = 'MC64')
              if ( sclsub .ne. '' ) then
                 if ( iprintctl(2) ) then
                    write(* ,9130) dictionary(ikey)
                    write(10,9130) dictionary(ikey)
                 end if
              end if
           end if
        else if ( lsssub .eq. 'MA97' .and. lss_ma97() ) then
           call setalgparam(val_lsssubACC = lsssub)
           if ( sclsub .eq. 'MC64' .or. sclsub .eq.'MC77' .or. &
                sclsub .eq. 'MC30' .or. sclsub .eq. 'NONE' ) then
              call setalgparam(val_sclsubACC = sclsub)
           else
              call setalgparam(val_sclsubACC = 'MC64')
              if ( sclsub .ne. '' ) then
                 if ( iprintctl(2) ) then
                    write(* ,9130) dictionary(ikey)
                    write(10,9130) dictionary(ikey)
                 end if
              end if
           end if
        else if ( lsssub .ne. '' ) then
           if ( iprintctl(2) ) then
              write(* ,9130) dictionary(ikey)
              write(10,9130) dictionary(ikey)
           end if
        end if
     else
        if ( .not. seconde ) then
           if ( iprintctl(2) ) then
              write(* ,9110) dictionary(ikey)
              write(10,9110) dictionary(ikey)
           end if
        else
           if ( iprintctl(2) ) then
              write(* ,9120) dictionary(ikey)
              write(10,9120) dictionary(ikey)
           end if
        end if
     end if

  else if ( ikey .eq. 3 ) then
     if ( seconde .and. lsssubTR .ne. 'NONE' ) then

        if ( ikey .eq. 3 ) then
           call setalgparam(val_innslvr = 'TR')
        end if

        lsssub = ''
        sclsub = ''
        if ( GetWord(str,1,ff,ll) ) then
           lsssub = str(ff:ll)
           if ( GetWord(str,ll+1,ff,ll) ) sclsub = str(ff:ll)
        end if
        if ( lsssub .eq. 'MA57' .and. lss_ma57() ) then
           call setalgparam(val_lsssubTR = lsssub)
           if ( sclsub .eq. 'MC64' .or. sclsub .eq. 'NONE' ) then
              call setalgparam(val_sclsubTR = sclsub)
           else
              call setalgparam(val_sclsubTR = 'MC64')
              if ( sclsub .ne. '' ) then
                 if ( iprintctl(2) ) then
                    write(* ,9130) dictionary(ikey)
                    write(10,9130) dictionary(ikey)
                 end if
              end if
           end if
        else if ( lsssub .ne. '' ) then
           if ( iprintctl(2) ) then
              write(* ,9130) dictionary(ikey)
              write(10,9130) dictionary(ikey)
           end if
        end if
     else
        if ( .not. seconde ) then
           if ( iprintctl(2) ) then
              write(* ,9110) dictionary(ikey)
              write(10,9110) dictionary(ikey)
           end if
        else
           if ( iprintctl(2) ) then
              write(* ,9120) dictionary(ikey)
              write(10,9120) dictionary(ikey)
           end if
        end if
     end if

  else if ( ikey .eq. 5 .or. ikey .eq. 6 ) then
     if ( seconde .and. lsssubNW .ne. 'NONE' ) then

        if ( ikey .eq. 5 ) then
           call setalgparam(val_innslvr = 'NW')
        end if

        lsssub = ''
        sclsub = ''
        if ( GetWord(str,1,ff,ll) ) then
           lsssub = str(ff:ll)
           if ( GetWord(str,ll+1,ff,ll) ) sclsub = str(ff:ll)
        end if
        if ( lsssub .eq. 'MA57' .and. lss_ma57() ) then
           call setalgparam(val_lsssubNW = lsssub)
           if ( sclsub .eq. 'MC64' .or. sclsub .eq. 'NONE' ) then
              call setalgparam(val_sclsubNW = sclsub)
           else
              call setalgparam(val_sclsubNW = 'MC64')
              if ( sclsub .ne. '' ) then
                 if ( iprintctl(2) ) then
                    write(* ,9130) dictionary(ikey)
                    write(10,9130) dictionary(ikey)
                 end if
              end if
           end if
        else if ( lsssub .eq. 'MA86' .and. lss_ma86() ) then
           call setalgparam(val_lsssubNW = lsssub)
           if ( sclsub .eq. 'MC64' .or. sclsub .eq.'MC77' .or. &
                sclsub .eq. 'NONE' ) then
              call setalgparam(val_sclsubNW = sclsub)
           else
              call setalgparam(val_sclsubNW = 'MC64')
              if ( sclsub .ne. '' ) then
                 if ( iprintctl(2) ) then
                    write(* ,9130) dictionary(ikey)
                    write(10,9130) dictionary(ikey)
                 end if
              end if
           end if
        else if ( lsssub .eq. 'MA97' .and. lss_ma97() ) then
           call setalgparam(val_lsssubNW = lsssub)
           if ( sclsub .eq. 'MC64' .or. sclsub .eq. 'MC77' .or. &
                sclsub .eq. 'MC30' .or. sclsub .eq. 'NONE' ) then
              call setalgparam(val_sclsubNW = sclsub)
           else
              call setalgparam(val_sclsubNW = 'MC64')
              if ( sclsub .ne. '' ) then
                 if ( iprintctl(2) ) then
                    write(* ,9130) dictionary(ikey)
                    write(10,9130) dictionary(ikey)
                 end if
              end if
           end if
        else if ( lsssub .ne. '' ) then
           if ( iprintctl(2) ) then
              write(* ,9130) dictionary(ikey)
              write(10,9130) dictionary(ikey)
           end if
        end if
     else
        if ( .not. seconde ) then
           if ( iprintctl(2) ) then
              write(* ,9110) dictionary(ikey)
              write(10,9110) dictionary(ikey)
           end if
        else
           if ( iprintctl(2) ) then
              write(* ,9120) dictionary(ikey)
              write(10,9120) dictionary(ikey)
           end if
        end if
     end if

  else if ( ikey .eq. 7 .or. ikey .eq. 8 ) then

     if ( ikey .eq. 7 ) then
        call setalgparam(val_innslvr = 'TN')
     end if

     if ( str .eq. 'TRUE-HESSIAN' .or. str .eq. 'TRUEHP' ) then
        if ( seconde ) then
           call setalgparam(val_hptype = 'TRUEHP')
        else
           if ( iprintctl(2) ) then
              write(* ,9100) dictionary(ikey)
              write(10,9100) dictionary(ikey)
           end if
        end if

     else if ( str .eq. 'HESSIAN-APPROXIMATION'  .or. str .eq. 'HAPPRO' ) then
        if ( firstde ) then
           call setalgparam(val_hptype = 'HAPPRO')
        else
           if ( iprintctl(2) ) then
              write(* ,9115) dictionary(ikey)
              write(10,9115) dictionary(ikey)
           end if
        end if

     else if ( str .eq. 'INCREMENTAL-QUOTIENT'  .or. str .eq. 'INCQUO' ) then
        call setalgparam(val_hptype = 'INCQUO')

     else if ( str .ne. '' ) then
        if ( iprintctl(2) ) then
           write(* ,9140) dictionary(ikey)
           write(10,9140) dictionary(ikey)
        end if
     end if

  else if ( ikey .eq.  9 ) then
     call setproblvlu(val_rmfixv = .false.)

  else if ( ikey .eq. 10 ) then
     call setproblvlt(val_slacks = .true.)

  else if ( ikey .eq. 11 ) then
     call setproblvls(val_scale = .false.)

  else if ( ikey .eq. 12 ) then
     call setalgparam(val_ignoref = .true.)

  else if ( ikey .eq. 13 ) then
     epsfeas = dnum

  else if ( ikey .eq. 14 ) then
     epsopt = dnum

  else if ( ikey .eq. 15 ) then
     efstain = dnum

  else if ( ikey .eq. 16 ) then
     eostain = dnum

  else if ( ikey .eq. 17 ) then
     efacc = dnum

  else if ( ikey .eq. 18 ) then
     eoacc = dnum

  else if ( ikey .eq. 19 ) then
     call setouttyp(val_iprint = inum)

  else if ( ikey .eq. 20 ) then
     call setouttyp(val_ncomp = inum)

  else if ( ikey .eq. 21 ) then
     call setouttyp(val_solfnm = str)

  else if ( ikey .eq. 22 ) then
     call setouttyp(val_iprintctl5 = .true.)

  else if ( ikey .eq. 23 ) then
     call setalgparam(val_maxaccitgiven = inum)

  else if ( ikey .eq. 24 ) then
     call setalgparam(val_maxinnitgiven = inum)

  else if ( ikey .eq. 25 ) then
     call setalgparam(val_maxoutitgiven = inum)

  else if ( ikey .eq. 26 ) then
     call setalgparam(val_rhoinigiven = dnum)

  else if ( ikey .eq. 27 ) then
     call setalgparam(val_rhomaxgiven = dnum)
  end if

  ! RETURN

200 continue

  ! NON-EXECUTABLE STATEMENTS

 2000 format(BN,I20)
 3000 format(BN,F24.0)
 9020 format(  1X,'Warning: Ignoring unknown keyword    ',A45)
 9030 format(  1X,'Warning: Ignoring incomplete keyword ',A45)
 9040 format(  1X,A45)
 9041 format(  1X,A45,13X,I20)
 9042 format(  1X,A45,9X,1P,D24.8)
 9043 format(  1X,A45,1X,A32)

 9100 format(/,1X,'Warning: Ignoring keyword ',A41,'.',                 &
             /,1X,'This option requires information to compute the ',   &
                  'true product of the Hessian of the Augmented ',      &
                  'Lagrangian times a given vector. Possibilities for ',&
                  'the subroutines coded by the user are: (a) ',        &
                  'EVALJAC, EVALH, and EVALHC, (b) EVALGJAC and ',      &
                  'EVALHL, or (c) EVALGJACP and EVALHLP. If you ',      &
                  'already coded (a), (b), or (c), set parameter ',     &
                  'CODED appropiately.',/)
 9110 format(/,1X,'Warning: Ignoring keyword ',A41,                     &
             /,1X,'This option requires information to compute the ',   &
                  'Hessian matrix of the Augmented Lagrangian. ',       &
                  'Possibilities for the subroutines coded by the ',    &
                  'user are: (a) EVALJAC, EVALH, and EVALHC, or (b) ',  &
                  'EVALGJAC and EVALHL. If you already coded (a) or ',  &
                  '(b), set parameter CODED appropiately.',/)
 9115 format(/,1X,'Warning: Ignoring keyword ',A41,                     &
             /,1X,'This option requires information to compute the ',   &
                  'Jacobian of the constraints. Possibilities for the ',&
                  'subroutines coded by the user are: (a) EVALJAC or ', &
                  '(b) EVALGJAC. If you already coded (a) or (b), set ',&
                  'parameter CODED appropiately.',/)
 9120 format(/,1X,'Warning: Ignoring keyword ',A41,                     &
             /,1X,'This option requires the availability of a ',        &
                  'subroutine from HSL to solve linear systems. ',      &
                  'Valid options are: (a) MA57 for the trust-region '   &
                  'inner solver, (b) MA57, MA86, or MA97 for the ',     &
                  'Newton line-search inner solver, and (c) MA57, ',    &
                  'MA86, or MA97 for the Acceleration Process. If you ',&
                  'have any of them and it does not appear in the ',    &
                  'Algencan preamble where it is written ''Available ', &
                  'HSL subroutines'', see the compilation ',            &
                  'instructions for details.',/)
 9130 format(/,1X,'Warning: Ignoring the parameter(s) of keyword ',A41, &
             /,1X,'Specifying the linear systems solver is optional. ', &
                  'Valid options for the linear system solver are: ',   &
                  '(a) MA57 for the trust-region inner solver, (b) ',   &
                  'MA57, MA86, and MA97 for the Newton line-search ',   &
                  'inner solver, and (c) MA57, MA86, or MA97 for the ', &
                  'Acceleration Process. If you have any of them and ', &
                  'it does not appear in the Algencan preamble where ', &
                  'it is written ''Available HSL subroutines'', see ',  &
                  'the compilation instructions for details. In ',      &
                  'connection with the linear systems solver, ',        &
                  'available choices for scaling are: (a) MC64 or ',    &
                  'NONE for MA57, (b) MC64, MC77, or NONE for MA86, ',  &
                  'and (c) MC64, MC77, MC30, or NONE for MA97.',/)
 9140 format(/,1X,'Warning: Ignoring the parameter(s) of keyword ',A41, &
             /,1X,'Specifying the matrix-vector product to be used ',   &
                  'within Conjugate Gradients in the truncated Newton ',&
                  'approach is optional. Valid options are: (a) ',      &
                  'TRUE-HESSIAN (or TRUEHP) in case the Hessian of ',   &
                  'the Augmented Lagrangian is available, (b) ',        &
                  'HESSIAN-APPROXIMATION (or HAPPRO) in case the ',     &
                  'Jacobian of the constraints is available, and (c) ', &
                  'INCREMENTAL-QUOTIENT (or INCQUO), that has no ',     &
                  'requirements.')

end subroutine procline

! ******************************************************************
! ******************************************************************

logical function GetWord(str,istart,ifirst,ilast)

  ! ARRAY ARGUMENTS
  character(len=*), intent(in) :: str
  integer, intent(in) :: istart
  integer, intent(out) :: ifirst,ilast

! Returns the positions of the first (ifirst) and the last (ilast)
! character of a sequence of characters without blanks (word) within
! the string str and starting from position str(istart:istart). In
! addition, it returns .false. if such a word does not exist and
! .true. otherwise.

  ! LOCAL SCALARS
  integer :: i

  ! Find first character

  i = istart
10 if ( i .le. len(str) ) then
     if ( str(i:i) .eq. ' ' ) then
        i = i + 1
        go to 10
     end if
  end if

  ! Blank string

  if ( i .gt. len(str) ) then
     GetWord = .false.
     return
  end if

  ifirst = i

  ! Find last character

  i = ifirst + 1
20 if ( i .le. len(str) ) then
     if ( str(i:i) .ne. ' ' ) then
        i = i + 1
        go to 20
     end if
  end if

  ilast = i - 1

  GetWord = .true.
  return
  
end function GetWord

