C     ******************************************************************
C     ******************************************************************

      program algencanma

      implicit none

C     PARAMETERS
      integer nmax,mmax
      parameter ( nmax = 2 )
      parameter ( mmax = 2 )

C     LOCAL SCALARS
      logical checkder
      integer hnnzmax,inform,jcnnzmax,m,n,nvparam
      double precision cnorm,efacc,efstain,eoacc,epsfeas,epsopt,eostain,
     +        f,nlpsupn,snorm

C     LOCAL ARRAYS
      character * 80 specfnm,outputfnm,vparam(10)
      logical coded(11),equatn(mmax),linear(mmax)
      double precision l(nmax),lambda(mmax),u(nmax),x(nmax)

C     EXTERNAL SUBROUTINES
      external myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,
     +         myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp

C     Number of variables

      n = 2

C     Set lower bounds, upper bounds, and initial guess

      l(1) = - 10.0d0
      u(1) =   10.0d0

      l(2) = - 1.0d+20
      u(2) =   1.0d+20

      x(1) =   0.0d0
      x(2) =   0.0d0

C     Constraints

      m = 2

      equatn(1) = .false.
      equatn(2) = .false.

      lambda(1) = 0.0d0
      lambda(2) = 0.0d0

      linear(1) = .false.
      linear(2) = .true.

C     Coded subroutines

      coded( 1) = .false. ! fsub
      coded( 2) = .false. ! gsub
      coded( 3) = .false. ! hsub
      coded( 4) = .false. ! csub
      coded( 5) = .false. ! jacsub
      coded( 6) = .false. ! hcsub
      coded( 7) = .true.  ! fcsub
      coded( 8) = .false. ! gjacsub
      coded( 9) = .true.  ! gjacpsub
      coded(10) = .false. ! hlsub
      coded(11) = .true.  ! hlpsub

C     Upper bounds on the number of sparse-matrices non-null elements

      jcnnzmax = 0
      hnnzmax  = 0

C     Checking derivatives?

      checkder = .true.

C     Parameters setting

      epsfeas   = 1.0d-08
      epsopt    = 1.0d-08

      efstain   = sqrt( epsfeas )
      eostain   = epsopt ** 1.5d0

      efacc     = sqrt( epsfeas )
      eoacc     = sqrt( epsopt )

      outputfnm = ''
      specfnm   = ''

      nvparam = 0

      call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc,
     +myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax,
     +hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm,
     +specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,
     +checkder,f,cnorm,snorm,nlpsupn,inform)

      stop

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalf(n,x,f,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalg(n,x,g,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,n,hnnz,lim

C     ARRAY ARGUMENTS
      integer hcol(lim),hrow(lim)
      double precision hval(lim),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalc(n,x,ind,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,flag,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,ind,jcnnz,lim,n

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,hcnnz,ind,lim,n

C     ARRAY ARGUMENTS
      integer hccol(lim),hcrow(lim)
      double precision hcval(lim),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalfc(n,x,f,m,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

      flag = 0

C     Objective function

      f = x(2)

C     Constraints

      c(1) = x(1) ** 2 + 1 - x(n)
      c(2) = 2.0d0 - x(1) - x(n)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,
     +flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,jcnnz,lim,m,n

C     ARRAY ARGUMENTS
      integer jcfun(lim),jcvar(lim)
      double precision g(n),jcval(lim),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical gotj
      integer flag,m,n
      character work

C     ARRAY ARGUMENTS
      double precision g(n),p(m),q(n),x(n)

C     The meaning of argument work follows: work = 'j' or 'J' means that 
C     q is an input array argument and that p = Jacobian x q must be 
C     computed, while work = 't' or 'T' means that p is an input array 
C     argument and that q = Jacobian^t x p must be computed. Moreover, a 
C     capital letter (i.e. 'J' or 'T') means that the gradient of the 
C     objective function g must also be computed. A lower letter (i.e. 
C     'j' or 't') means that only the product with the Jacobian is 
C     required. In the later case, input array argument g MUST NOT be 
C     modified nor referenced.

      flag = 0

C     Gradient of the objective function

      if ( work .eq. 'J' .or. work .eq. 'T' ) then
          g(1) = 0.0d0
          g(2) = 1.0d0
      end if

C     p = Jacobian x q
      if ( work .eq. 'j' .or. work .eq. 'J' ) then
          p(1) = 2.0d0 * x(1) * q(1) - q(2)
          p(2) = - q(1) - q(2)

C     q = Jacobian^t x p
      else ! if ( work .eq. 't' .or. work .eq. 'T' ) then
          q(1) = 2.0d0 * x(1) * p(1) - p(2)
          q(2) = - p(1) - p(2)
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,
     +lim,lmem,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical lmem
      integer flag,hlnnz,lim,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlcol(lim),hlrow(lim)
      double precision hlval(lim),lambda(m),sc(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

      flag = 0

      hp(1) = 2.0d0 * sc(1) * lambda(1) * p(1)
      hp(2) = 0.0d0

      end
