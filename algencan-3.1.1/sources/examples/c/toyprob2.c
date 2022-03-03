#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

void c_algencan(void *myevalf, void *myevalg, void *myevalh, void *myevalc,
	void *myevaljac, void *myevalhc, void *myevalfc, void *myevalgjac,
	void *myevalgjacp, void *myevalhl, void *myevalhlp, int jcnnzmax,
       	int hnnzmax,double *epsfeas, double *epsopt, double *efstin,
	double *eostin, double *efacc, double *eoacc, char *outputfnm,
	char *specfnm, int nvparam,char **vparam, int n, double *x,
	double *l, double *u, int m, double *lambda, _Bool *equatn,
	_Bool *linear, _Bool *coded, _Bool checkder, double *f,
	double *cnorm, double *snorm, double *nlpsupn,int *inform);

/* ******************************************************************
   ****************************************************************** */

void myevalf(int n, double *x, double *f, int *flag) {

   *flag = - 1;
}

/* ******************************************************************
   ****************************************************************** */

void myevalg(int n, double *x, double *g, int *flag) {

   *flag = - 1;
}

/* ******************************************************************
   ****************************************************************** */

void myevalh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz,
	     int lim, _Bool *lmem, int *flag) {
   *flag = - 1;
}

/* ******************************************************************
   ****************************************************************** */

void myevalc(int n, double *x, int ind, double *c, int *flag) {

   *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

void myevaljac(int n, double *x, int ind, int *jcvar, double *jcval,
	       int *jcnnz, int lim, _Bool *lmem, int *flag) {
   *flag = - 1;
}

/* ******************************************************************
   ****************************************************************** */

void myevalhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
	      int *hcnnz, int lim, _Bool *lmem, int *flag) {
   *flag = - 1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalfc(int n, double *x, double *f, int m, double *c, int *flag) {

  *flag = 0;

  /* Objective function */
  *f = x[1];

  /* Constraints */
  c[0] = pow( x[0], 2.0 ) + 1.0 - x[1];
  c[1] = 2.0 - x[0] - x[1];
}

/* *****************************************************************
   ***************************************************************** */

void myevalgjac(int n, double *x, double *g, int m, int *jcfun, int *jcvar,
		double *jcval, int *jcnnz, int lim, _Bool *lmem, int *flag) {

  *flag = 0;
  *lmem = 0;

  g[0] = 0.0;
  g[1] = 1.0;

  *jcnnz = 4;

  if( *jcnnz > lim ) {
    *lmem = 1;
    return;
  }

  jcfun[0] = 0;
  jcvar[0] = 0;
  jcval[0] = 2.0 * x[0];

  jcfun[1] = 0;
  jcvar[1] = 1;
  jcval[1] = - 1.0;

  jcfun[2] = 1;
  jcvar[2] = 0;
  jcval[2] = - 1.0;

  jcfun[3] = 1;
  jcvar[3] = 1;
  jcval[3] = - 1.0;
}

/* *****************************************************************
   ***************************************************************** */

void myevalgjacp(int n, double *x, double *g, int m, double *p, double *q,
		 char work, _Bool *gotj, int *flag) {

  *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalhl(int n, double *x, int m, double *lambda, double scalef,
	      double *scalec, int *hlrow, int *hlcol, double *hlval,
	      int *hlnnz, int lim, _Bool *lmem, int *flag) {

  *flag = 0;
  *lmem = 0;
  
  *hlnnz = 1;

  if( *hlnnz > lim ) {
    *lmem = 1;
    return;
  }
  
  hlrow[0] = 0;
  hlcol[0] = 0;
  hlval[0] = scalec[0] * lambda[0] * 2.0;
}

/* *****************************************************************
   ***************************************************************** */

void myevalhlp(int n, double *x, int m, double *lambda, double scalef,
	       double *scalec, double *p, double *hp, _Bool *goth, 
	       int *flag) {
  
  *flag = -1;
}

/* ******************************************************************
   ****************************************************************** */

int main() {
  _Bool  checkder;
  int    hnnzmax,hnnzmax1,hnnzmax2,hnnzmax3,i,jcnnzmax,inform,m,n,nvparam;
  double cnorm,efacc,efstin,eoacc,eostin,epsfeas,epsopt,f,nlpsupn,snorm;
  
  char   *specfnm, *outputfnm, **vparam;
  _Bool  coded[11],*equatn,*linear;
  double *l,*lambda,*u,*x;
  
  n = 2;
  m = 2;
  
  /* Memory allocation */
  x      = (double *) malloc(n * sizeof(double));
  l      = (double *) malloc(n * sizeof(double));
  u      = (double *) malloc(n * sizeof(double));
  lambda = (double *) malloc(m * sizeof(double));
  equatn = (_Bool  *) malloc(m * sizeof(_Bool ));
  linear = (_Bool  *) malloc(m * sizeof(_Bool ));
  
  if (     x == NULL ||      l == NULL ||      u == NULL ||
      lambda == NULL || equatn == NULL || linear == NULL ) {
    
    printf( "\nC INTERFACE ERROR: It was not possible to allocate memory.\n" );
    exit( 0 );
    
  }

  /* Initial point */
  for(i = 0; i < n; i++) x[i] = 0.0;
  
  /* Lower and upper bounds */
  l[0] = - 10.0;
  l[1] = - 1.0e20;
  
  u[0] = 10.0;
  u[1] = 1.0e20;
  
  /* For each constraint i, set equatn[i] = 1. if it is an equality
     constraint of the form c_i(x) = 0, and set equatn[i] = 0 if it is
     an inequality constraint of the form c_i(x) <= 0. */
  equatn[0] = 0;
  equatn[1] = 0;

  /* For each constraint i, set linear[i] = 1 if it is a linear
     constraint, otherwise set linear[i] = 0 */
  linear[0] = 0;
  linear[1] = 1;
  
  /* Lagrange multipliers approximation. */
  for( i = 0; i < m; i++ ) lambda[i] = 0.0;
  
  /* In this C interface evalf, evalg, evalh, evalc, evaljac and
     evalhc are present. evalfc, evalgjac, evalhl and evalhlp are
     not. */
  
  coded[0]  = 0; /* fsub     */
  coded[1]  = 0; /* gsub     */
  coded[2]  = 0; /* hsub     */
  coded[3]  = 0; /* csub     */
  coded[4]  = 0; /* jacsub   */
  coded[5]  = 0; /* hcsub    */
  coded[6]  = 1; /* fcsub    */
  coded[7]  = 1; /* gjacsub  */
  coded[8]  = 0; /* gjacpsub */
  coded[9]  = 1; /* hlsub    */
  coded[10] = 0; /* hlpsub   */
 
  /* Upper bounds on the number of sparse-matrices non-null
     elements */
  jcnnzmax = 4;
  hnnzmax1 = 0;
  hnnzmax2 = 1;
  hnnzmax3 = 6;
  hnnzmax  = hnnzmax1 + hnnzmax2 + hnnzmax3;

  /* Check derivatives? */
  checkder = 0;

  /* Parameters setting */
  epsfeas  = 1.0e-08;
  epsopt   = 1.0e-08;

  efstin   = sqrt( epsfeas );
  eostin   = pow( epsopt, 1.5 );

  efacc    = sqrt( epsfeas );
  eoacc    = sqrt( epsopt );

  outputfnm = "algencan.out";
  specfnm   = "algencan.dat";

  nvparam = 2;
  
  /* Allocates VPARAM array */
  vparam = ( char ** ) malloc( nvparam * sizeof( char * ) );

  /* Set algencan parameters */
  vparam[0] = "ITERATIONS-OUTPUT-DETAILS 10";
  vparam[1] = "TRUNCATED-NEWTON-LINE-SEARCH-INNER-SOLVER";

  nvparam = 0;

  /* Optimize */
  c_algencan(&myevalf,&myevalg,&myevalh,&myevalc,&myevaljac,&myevalhc,&myevalfc,
	     &myevalgjac,&myevalgjacp,&myevalhl,&myevalhlp,jcnnzmax,hnnzmax,
	     &epsfeas,&epsopt,&efstin,&eostin,&efacc,&eoacc,outputfnm,specfnm,
	     nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,checkder,
	     &f,&cnorm,&snorm,&nlpsupn,&inform);

  /* Memory deallocation */
  free(x     );
  free(l     );
  free(u     );
  free(lambda);
  free(equatn);
  free(linear);
}
