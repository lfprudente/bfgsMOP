/* =================================================================
   File: amplwrapper.c
   =================================================================

   =================================================================
   Module: AMPL Wrapper
   =================================================================

   Last update of any of the component of this module:

   November 3, 2016.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO home page:

   http://www.ime.usp.br/~egbirgin/tango/

   *****************************************************************
   *****************************************************************

   TANGO Project

   License:

   All the TANGO Project components are free software; you can
   redistribute it and/or modify it under the terms of the GNU General
   Public License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. You can also find the GPL on the GNU web site.

   Non-free versions of TANGO are available under terms different from
   those of the General Public License. Professors J. M. Mart\'{\i}nez
   (martinez@ime.unicamp.br, martinezimecc@gmail.com) or E. G. Birgin
   (egbirgin@ime.usp.br, egbirgin@gmail.com) should be contacted for
   more information related to such a license, future developments
   and/or technical support.

   In addition, we kindly ask you to acknowledge the TANGO Project and
   its authors in any program or publication in which you use of the
   TANGO Project components. For published works that use ALGENCAN we
   suggest referencing:

   R. Andreani, E. G. Birgin, J. M. Mart\'{\i}nez and M. L. Schuverdt, "On
   Augmented Lagrangian methods with general lower-level constraints",
   SIAM Journal on Optimization 18, pp. 1286-1309, 2007.

   and

   R. Andreani, E. G. Birgin, J. M. Mart\'{\i}nez and M. L. Schuverdt,
   "Augmented Lagrangian methods under the Constant Positive Linear
   Dependence constraint qualification", Mathematical Programming 111,
   pp. 5-32, 2008

   For published works that use GENCAN we suggest referencing:

   E. G. Birgin and J. M. Mart\'{\i}nez, "Large-scale active-set
   box-constrained optimization method with spectral projected
   gradients", Computational Optimization and Applications 23,
   pp. 101-125, 2002.

   For published works that use SPG we suggest referencing:

   E. G. Birgin, J. M. Mart\'{\i}nez and M. Raydan, "Nonmonotone spectral
   projected gradient methods on convex sets", SIAM Journal on
   Optimization 10, pp. 1196-1211, 2000,

   and

   E. G. Birgin, J. M. Mart\'{\i}nez and M. Raydan, "Algorithm 813: SPG -
   software for convex-constrained optimization", ACM Transactions on
   Mathematical Software 27, pp. 340-349, 2001.

   (See also other related works in the TANGO Project home page:

   http://www.ime.usp.br/~egbirgin/tango/)

   ***************************************************************** */

/* *****************************************************************
   ***************************************************************** */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "asl.h"
#include "getstub.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

ASL *asl;
int *cmap,*slaind,useslacks;
double *ca,*cb,objsign;
char *stub,message[200];
double *r, *J;

/* *****************************************************************
   ***************************************************************** */

void inip(int *n,double **x,double **l,double **u,int *m,double **lambda,
	  _Bool **equatn,_Bool **linear,int *hnnzmax,int *jcnnzmax) {

/* This subroutine set some problem data. */

   int flag,i,ind,jcnnz,k,nconstr,nslacks;
   double cl,cu;
   struct cgrad *cg;
   FILE *nl;

   strcpy( message, "" );

   /* Sparse jacobian option. */
   asl->i.congrd_mode = 1;

   /* Name of file that contains the description of the problem. */
   nl = jac0dim(stub, (fint)strlen(stub));

   /* Initial point primal and dual (if provided by ampl). */

   X0     = (real *) malloc( n_var * sizeof (real) );
   havex0 = (char *) malloc( n_var * sizeof (char) );

   pi0     = (real *) malloc( n_con * sizeof (real) );
   havepi0 = (char *) malloc( n_con * sizeof (char) );

   if ( X0 == NULL || havex0 == NULL || pi0 == NULL || havepi0 == NULL ) {

     write_sol( "\nAMPL INTERFACE ERROR: in inip, malloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }

   pfgh_read(nl,ASL_findgroups);

   /* Verifies whether the problem has more than a single objective function. */

   if ( n_obj > 1 ) {
      write_sol( "\nAlgencan can not handle more than an objective function.\n",
		 NULL, NULL, NULL );
      exit( 0 );
   }

   if ( n_obj > 0 && objtype[0] )
     objsign = - 1.0;
   else
     objsign =   1.0;

   /* Verifies whether the problem has integer variables. */

   if ( nbv > 0 || niv > 0 || nlvci > 0 || nlvoi > 0 ) {

     write_sol( "\nALGENCAN can not handle integer variables.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }

   /* Number of variables. */

   useslacks = 1;

   *n = n_var + ( useslacks ? nranges : 0 );

   /* Number of constraints (equalities plus inequalities). */

   *m = n_con + ( useslacks ? 0 : nranges );

   /* Memory allocation. */

   *x      = (double *) malloc(((*n)+(*m)) * sizeof(double));
   *l      = (double *) malloc(((*n)+(*m)) * sizeof(double));
   *u      = (double *) malloc(((*n)+(*m)) * sizeof(double));
   *lambda = (double *) malloc((*m) * sizeof(double));
   *equatn = (_Bool  *) malloc((*m) * sizeof(_Bool));
   *linear = (_Bool  *) malloc((*m) * sizeof(_Bool));

   cmap    = (int *)    malloc((*m) * sizeof(int));
   ca      = (double *) malloc((*m) * sizeof(double));
   cb      = (double *) malloc((*m) * sizeof(double));
   slaind  = (int *)    malloc((*m) * sizeof(int));

   r       = (real *)   malloc( n_con * sizeof (real) );
   J       = (real *)   malloc( nzc   * sizeof (real) );

   if (      *x == NULL ||      *l == NULL ||      *u == NULL ||
        *lambda == NULL || *equatn == NULL || *linear == NULL ||
         slaind == NULL ||    cmap == NULL ||      ca == NULL ||
             cb == NULL ||       r == NULL ||       J == NULL ) {

     write_sol( "\nAMPL INTERFACE ERROR: in inip, malloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );

   }

   /* Initial point. */

   for(i = 0; i < n_var; i++)
     if ( havex0[i] )
       (*x)[i] = X0[i];
     else
       (*x)[i] = 0.0;

   /* Lower and upper bounds. */

   for(i = 0; i < n_var; i++) {
     (*l)[i] = LUv[2 * i];
     (*u)[i] = LUv[2 * i + 1];
   }

  /* For each constraint i, set equatn(i) = 1. if it is an equality
     constraint of the form c_i(x) = 0, and set equatn(i) = 0 if it is
     an inequality constraint of the form c_i(x) <= 0. Moreover, set
     linear(i) = 1 if it is a linear constraint, otherwise set
     linear(i) = 0. */

   nconstr = 0;
   nslacks = 0;

   for(i = 0; i < n_con; i++) {

      cl = LUrhs[2 * i];
      cu = LUrhs[2 * i + 1];

      /* Equality constraint. */

      if ( cl == cu ) {
	 (*equatn)[nconstr] = 1;
	 (*linear)[nconstr] = (i < nlc ? 0 : 1);

	 slaind[nconstr]    = - 1;
         cmap[nconstr]      = i;
	 ca[nconstr]        = 1.0;
	 cb[nconstr]        = - cl;

	 if ( havepi0[i] )
	    (*lambda)[nconstr] = objsign * pi0[i];
	 else
	    (*lambda)[nconstr] = 0.0;
 
	 nconstr++;
      }

      /* Ranged inequality constraint: add slack or split it. */

      else if ( cl > negInfinity && cu < Infinity ) {

	 /* Replace by c(x) - s = 0 and l <= s <= u. */

	 if ( useslacks ) {
	    (*equatn)[nconstr] = 1;
	    (*linear)[nconstr] = (i < nlc ? 0 : 1);

	    k = n_var + nslacks;

	    slaind[nconstr] = k;
	    cmap[nconstr]   = i;
	    ca[nconstr]     = 1.0;
	    cb[nconstr]     = 0.0;

	    if ( havepi0[i] )
	       (*lambda)[nconstr] = objsign * pi0[i];
	    else
	       (*lambda)[nconstr] = 0.0;
 
	    (*l)[k] = cl;
	    (*u)[k] = cu;

	    (*x)[k] = conival(i, (real *)(*x), (fint *)&flag );

	    if ( flag != 0 ) {
	       sprintf(message,"AMPL INTERFACE ERROR: conival returned %d",(int)flag);
	       write_sol(message,NULL, NULL, NULL );
	       exit( 0 );
	    }

	    (*x)[k] = max(cl, min( (*x)[k], cu));

	    nconstr++;
	    nslacks++;
	 }

	 /* Split into c(x) - u <= 0 and l - c(x) <= 0. */

	 else {
	    double lambda1 = 0.0;
	    double lambda2 = 0.0;

	    if ( havepi0[i] ) {
	       double cval = conival(i, (real *)(*x), (fint *)&flag );

	       if ( flag != 0 ) {
		  sprintf(message,"AMPL INTERFACE ERROR: conival returned %d",(int)flag);
		  write_sol(message,NULL, NULL, NULL );
		  exit( 0 );
	       }

	       if ( fabs( cval - cu ) < fabs( cl - cval ) )
		  lambda1 = objsign * pi0[i];
	       else
		  lambda2 = - objsign * pi0[i];
	    }

	    (*equatn)[nconstr] = 0;
	    (*linear)[nconstr] = (i < nlc ? 0 : 1);

	    slaind[nconstr]    = - 1;
	    cmap[nconstr]      = i;
	    ca[nconstr]        = 1.0;
	    cb[nconstr]        = - cu;

	    (*lambda)[nconstr] = lambda1;

	    nconstr++;

	    (*equatn)[nconstr] = 0;
	    (*linear)[nconstr] = (i < nlc ? 0 : 1);

	    slaind[nconstr]    = - 1;
	    cmap[nconstr]      = i;
	    ca[nconstr]        = - 1.0;
	    cb[nconstr]        = cl;

	    (*lambda)[nconstr] = lambda2;

	    nconstr++;
	 }
      }

      /* Inequality constraint of type c(x) <= u. */

      else if ( cu < Infinity ) {
	 (*equatn)[nconstr] = 0;
	 (*linear)[nconstr] = (i < nlc ? 0 : 1);

	 slaind[nconstr]    = - 1;
	 cmap[nconstr]      = i;
	 ca[nconstr]        = 1.0;
	 cb[nconstr]        = - cu;

	 if ( havepi0[i] )
	    (*lambda)[nconstr] = objsign * pi0[i];
	 else
	    (*lambda)[nconstr] = 0.0;
 
	 nconstr++;
      }
      /* Inequality constraint of type l <= c(x). */

      else if ( cl > negInfinity ) {
	 (*equatn)[nconstr] = 0;
	 (*linear)[nconstr] = (i < nlc ? 0 : 1);

	 slaind[nconstr]    = - 1;
	 cmap[nconstr]      = i;
	 ca[nconstr]        = - 1.0;
	 cb[nconstr]        = cl;

	 if ( havepi0[i] )
	    (*lambda)[nconstr] = objsign * pi0[i];
	 else
	    (*lambda)[nconstr] = 0.0;
 
	 nconstr++;
      }
   }

   if ( *jcnnzmax == -1 ) {
     if ( useslacks )
       *jcnnzmax = nzc;
     else
       *jcnnzmax = 2 * nzc;
   }
   
   if ( *hnnzmax == -1 ) {
      if ( n_obj == 0 )
	 *hnnzmax = sphsetup(-1,0,1,1);
      else
	 *hnnzmax = sphsetup(-1,1,1,1);

      for ( ind=0; ind<*m; ind++ ) {
	 jcnnz = 0;
	 for(cg = Cgrad[cmap[ind]]; cg != NULL; cg = cg->next)
	    jcnnz++;
	 if ( slaind[ind] != -1 )
	    jcnnz++;
	 *hnnzmax += jcnnz * ( jcnnz + 1 ) / 2;
      }

      if ( *hnnzmax > *n * *n )
	 *hnnzmax = *n * *n;
   }

   free( X0 );
   free( havex0 );
   free( pi0 );
   free( havepi0 );
}

/* ******************************************************************
   ****************************************************************** */

void amplf(int n, double *x, double *f, int *flag) {

   /* This subroutine computes the objective function. */

   *flag = 0;

   if ( n_obj == 0 )
      *f = 0;

   else {
      *f = objval(0, (real *)x, (fint *)flag);

      if ( *flag != 0 )
	 sprintf( message, "AMPL INTERFACE ERROR: objval returned %d",*flag );

      *f *= objsign;
   }
}

/* ******************************************************************
   ****************************************************************** */

void amplg(int n, double *x, double *g, int *flag) {

/* This subroutine computes the gradient vector of the objective
   function. */

   int i;

   *flag = 0;

   if ( n_obj == 0 ) {
     for ( i = 0; i < n; i++ )
       g[i] = 0;
   }

   else {
     objgrd(0, (real *)x, (real *)g, (fint *)flag);

     if ( *flag != 0 )
       sprintf( message, "AMPL INTERFACE ERROR: objgrd returned %d",*flag);

     for ( i = 0; i < n_var; i++ )
	g[i] *= objsign;

     for ( i = n_var; i < n; i++ )
	g[i] = 0.0;
   }

}

/* ******************************************************************
   ****************************************************************** */

void amplh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz,
	   int lim, _Bool *lmem, int *flag) {
   *flag = - 1;
}

/* ******************************************************************
   ****************************************************************** */

void amplc(int n, double *x, int ind, double *c, int *flag) {

   /* This subroutine computes the ind-th constraint. */

   *flag = 0;

   *c = conival( cmap[ind], (real *)x, (fint *)flag );

   if ( *flag != 0 )
      sprintf( message, "AMPL INTERFACE ERROR: conival returned %d",*flag );

   if ( slaind[ind] == -1 )
      *c = *c * ca[ind] + cb[ind];

   else
      *c -= x[ slaind[ind] ];

}

/* ******************************************************************
   ****************************************************************** */

void ampljac(int n, double *x, int ind, int *jcvar, double *jcval,
	     int *jcnnz, int lim, _Bool *lmem, int *flag) {

/* This subroutine computes the gradient vector of the ind-th
   constraint. */

   int i;
   struct cgrad *cg;

   *flag = 0;
   *lmem = 0;

   *jcnnz = 0;

   congrd( cmap[ind], (real *)x, (real *)(&(jcval[*jcnnz])), (fint *)flag );

   if ( *flag != 0 )
      sprintf( message, "AMPL INTERFACE ERROR: congrd returned %d",*flag );

   for(cg = Cgrad[ cmap[ind] ]; cg != NULL; cg = cg->next) {
      jcvar[*jcnnz]  = cg->varno;
      jcval[*jcnnz] *= ca[ind];
      (*jcnnz)++;
   }

   if ( slaind[ind] != -1 ) {
      jcvar[*jcnnz] = slaind[ind];
      jcval[*jcnnz] = - 1.0;
      (*jcnnz)++;
   }
}

/* ******************************************************************
   ****************************************************************** */

void amplhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
	    int *hcnnz, int lim, _Bool *lmem, int *flag) {
   *flag = - 1;
}

/* *****************************************************************
***************************************************************** */

void amplfc(int n,double *x,double *f,int m,double *c,int *flag) {

   /* This subroutine computes the objective function and the
      constraints. */

   int ind;

   *flag = 0;

   if ( n_obj == 0 )
      *f = 0;

   else {
      *f = objval(0, (real *)x, (fint *)flag);

      if ( *flag != 0 )
	 sprintf( message, "AMPL INTERFACE ERROR: objval returned %d",*flag );

      *f *= objsign;
   }

   conval((real *)x, (real *)r, (fint *)flag);

   if ( *flag != 0 )
      sprintf( message, "AMPL INTERFACE ERROR: conval returned %d",*flag );

   for ( ind=0; ind<m; ind++ )
      if ( slaind[ind] == -1 )
	 c[ind] = r[cmap[ind]] * ca[ind] + cb[ind];

      else
	 c[ind] = r[cmap[ind]] - x[ slaind[ind] ];
}

/* *****************************************************************
   ***************************************************************** */

void amplgjac(int n,double *x,double *g,int m,int *jcfun,int *jcvar,
double *jcval,int *jcnnz, int lim, _Bool *lmem,int *flag) {

/* This subroutine computes the gradient vector of the objective
   function and the sparse jacobian of the constraints. */

   int i, ind;
   struct cgrad *cg;

   *flag = 0;
   *lmem = 0;

   if ( n_obj == 0 ) {
     for ( i = 0; i < n; i++ )
       g[i] = 0;
   }

   else {
     objgrd(0, (real *)x, (real *)g, (fint *)flag);

     if ( *flag != 0 )
       sprintf( message, "AMPL INTERFACE ERROR: objgrd returned %d",*flag);

     for ( i = 0; i < n_var; i++ )
	g[i] *= objsign;

     for ( i = n_var; i < n; i++ )
	g[i] = 0.0;
   }

   jacval((real *)x, (real *)J, (fint *)flag);

   if ( *flag != 0 )
      sprintf( message, "AMPL INTERFACE ERROR: jacval returned %d",*flag);

   *jcnnz = 0;
   for ( ind=0; ind<m; ind++ ) {
      for(cg = Cgrad[cmap[ind]]; cg != NULL; cg = cg->next) {
	 jcfun[*jcnnz] = ind;
	 jcvar[*jcnnz] = cg->varno;
	 jcval[*jcnnz] = J[cg->goff] * ca[ind];
	 (*jcnnz)++;
      }

      if ( slaind[ind] != -1 ) {
	 jcfun[*jcnnz] = ind;
	 jcvar[*jcnnz] = slaind[ind];
	 jcval[*jcnnz] = - 1.0;
	 (*jcnnz)++;
      }
   } 
}

/* *****************************************************************
   ***************************************************************** */

void amplgjacp(int n, double *x, double *g, int m, double *p, double *q,
	       char work, _Bool *gotj, int *flag) {

  *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void amplhl(int n,double *x,int m,double *lambda,double scalef,
double *scalec,int *hlrow,int *hlcol,double *hlval,int *hlnnz,
int lim, _Bool *lmem, int *flag) {

/* This subroutine computes the sparse Hessian of the Lagrangian. */

   int i,j;
   double *coeff,*OW;

   *flag = 0;
   *lmem = 0;

   OW    = (real *) calloc( n_obj, sizeof (real) );
   coeff = (real *) calloc( n_con, sizeof (real) );

   if ( OW == NULL || coeff == NULL ) {

     write_sol( "\nAMPL INTERFACE ERROR: in evalhl, calloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );

   }

   OW[0] = objsign * scalef;

   for ( i = 0; i < m; i++ )
     coeff[cmap[i]] += ca[i] * lambda[i] * scalec[i];

   if ( n_obj == 0 )
     sphes((real *)hlval,-1,NULL,(real *)coeff);
   else
     sphes((real *)hlval,-1,(real *)OW,(real *)coeff);

   *hlnnz = sputinfo->hcolstarts[n_var];

   for ( j = 0; j < n_var; j++ ) {
     for ( i = sputinfo->hcolstarts[j]; i < sputinfo->hcolstarts[j+1]; i++ ) {

       /* The AMPL subroutine returns the upper triangle (as lists of columns
          with elements from the first row to the diagonal). Invert indices to
          return the lower triangle. */

       hlrow[i] = j;
       hlcol[i] = sputinfo->hrownos[i];
     }
   }

   free( OW );
   free( coeff );
}

/* *****************************************************************
   ***************************************************************** */

void amplhlp(int n, double *x, int m, double *lambda, double scalef,
	     double *scalec, double *p, double *hp, _Bool *goth, 
	     int *flag) {
  
/* This subroutine computes the product of the Hessian of the
   Lagrangian by a given vector. */

   int j;
   double *coeff,*OW;

   *flag = 0;

   OW    = (real *) calloc( n_obj, sizeof (real) );
   coeff = (real *) calloc( n_con, sizeof (real) );

   if ( OW == NULL || coeff == NULL ) {

     write_sol( "\nAMPL INTERFACE ERROR: in evalhlp, calloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );

   }

   OW[0] = objsign * scalef;

   for ( j = 0; j < m; j++ )
     coeff[cmap[j]] += ca[j] * lambda[j] * scalec[j];

   hvcomp((real *)hp,(real *)p,n_obj,(real *)OW,(real *)coeff);

   free( OW );
   free( coeff );
}

/* *****************************************************************
   ***************************************************************** */

void endp(int n,double **x,double **l,double **u,int m,double **lambda,
	  _Bool **equatn,_Bool **linear) {

/* This subroutine can be used to do some extra job after the solver
   has found the solution. */

   int i;
   double * ampl_lambda;

   ampl_lambda = (real *) calloc( n_con, sizeof (real) );

   if ( ampl_lambda == NULL) {
     write_sol( "\nAMPL INTERFACE ERROR: in endp, calloc returned NULL.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }

   for ( i = 0; i < m; i++ )
     ampl_lambda[cmap[i]] += ca[i] * (*lambda)[i];

   write_sol( message, *x, ampl_lambda, NULL );

   free( ampl_lambda );
   
   free( *x );
   free( *l );
   free( *u );
   free( *lambda );
   free( *equatn );
   free( *linear );

   free( cmap );
   free( ca );
   free( cb );
   free( slaind );
   
   free( r );
   free( J );
}

/* *****************************************************************
   ***************************************************************** */

void ampl_setp(int n,double * x) {

   xknown((real *)x);
}

/* *****************************************************************
   ***************************************************************** */

void ampl_unsetp() {

   xunknown();
}

/* *****************************************************************
   ***************************************************************** */

void getOptions(char *argv[],int *jcnnzmax,int *hnnzmax,double *epsfeas,
		double *epsopt,double *efstain,double *eostain,double *efacc,
		double *eoacc,char **outputfnm,char **specfnm) {

  char* jcnnzmaxDescription = "INTEGER\n\n"
"Upper bound on the number of nonnull elements in the Jacobian of the           \n"
"This parameter is automatically set by AMPL and should be modified by the user \n"
"only in very specific situations.                                              \n";

  char* hnnzmaxDescription = "INTEGER\n\n"
"Upper bound on the number of nonnull elements in the Hessian of the Lagrangian.\n"
"If the trust-region method is used in Algencan, in addition, the number of     \n"
"nonnull elements in the matrix given by the transpose of the Jacobian times the\n"
"Jacobian must me added to hnnzmax. The latter value, automatically computed by \n"
"this AMPL interface of Algencan, may be overestimated. If this is the case,    \n"
"the should set it to a smaller value in order to run Algencan.                 \n";

  char* epsfeasDescription = "DOUBLE PRECISION\n\n"
"Feasibility tolerance for the sup-norm of the constraints. (Ignored in the     \n"
"unconstrained and bound-constrained cases.)                                    \n";

  char* epsoptDescription = "DOUBLE PRECISION\n\n"
"Optimality tolerance for the sup-norm of the projected gradient of the         \n"
"Lagrangian in the constrained case and the sup-norm of the projected           \n"
"gradient of the objective function in the unconstrained and the bound-         \n"
"constrained cases.                                                             \n";

  char* efstainDescription = "DOUBLE PRECISION\n\n"
"Feasibility tolerance for the stopping criterion related to convergence to an  \n"
"infeasible point.                                                              \n";

  char* eostainDescription = "DOUBLE PRECISION\n\n"
"Optimality tolerance for the stopping criterion related to convergence to an   \n"
"infeasible point.                                                              \n";

  char* efaccDescription = "DOUBLE PRECISION\n\n"
"Feasibility threshold to launch the acceleration process.                      \n";

  char* eoaccDescription = "DOUBLE PRECISION\n\n"
"Optimality threshold to launch the acceleration process.                       \n";

  char *outputfnmDescription = "STRING (eighty characters long)\n\n"
"Name of the output file                                                        \n";

  char *specfnmDescription = "STRING (eighty characters long)\n\n"
"Name of the specification file                                                 \n";

  /* Must be in alphabetical order. */
  keyword keywds[] = {
    KW("efacc",D_val,efacc,efaccDescription),

    KW("efstain",D_val,efstain,efstainDescription),

    KW("eoacc",D_val,eoacc,eoaccDescription),

    KW("eostain",D_val,eostain,eostainDescription),

    KW("epsfeas",D_val,epsfeas,epsfeasDescription),

    KW("epsopt",D_val,epsopt,epsoptDescription),

    KW("hnnzmax",I_val,hnnzmax,hnnzmaxDescription),

    KW("jcnnzmax",I_val,jcnnzmax,jcnnzmaxDescription),

    KW("outputfnm",C_val,outputfnm,outputfnmDescription),

    KW("specfnm",C_val,specfnm,specfnmDescription)
  };

   Option_Info Oinfo = {"algencan","ALGENCAN","algencan_options",
			keywds,nkeywds,1,"Algencan 3.1.1"};

   *hnnzmax  = - 1;
   *jcnnzmax = - 1;

   /* File name that contains the description of the problem. */
   stub = getstops(argv,&Oinfo);

   if ( stub == NULL ) {
     write_sol( "\nThe 'stub' file of the problem is missing.\n",
                NULL, NULL, NULL );
     exit( 0 );
   }
}

/* *****************************************************************
   ***************************************************************** */

int main(int argc,char *argv[]) {

  _Bool checkder;
  int hnnzmax,inform,jcnnzmax,m,n,nvparam;
  double f,nlpsupn,efacc,efstain,eoacc,eostain,epsfeas,epsopt,cnorm,snorm;

  char *outputfnm, *specfnm, **vparam;
  _Bool  coded[11],*equatn,*linear;
  double *l,*lambda,*u,*x;

  char emptyfnm = '\0';

  asl = ASL_alloc(ASL_read_pfgh);

  /* In this AMPL interface amplfc, amplgjac, and amplhl are coded. */

  coded[0]  = 0; /* amplf     */
  coded[1]  = 0; /* amplg     */
  coded[2]  = 0; /* amplh     */
  coded[3]  = 0; /* amplc     */
  coded[4]  = 0; /* ampljac   */
  coded[5]  = 0; /* amplhc    */
  coded[6]  = 1; /* amplfc    */
  coded[7]  = 1; /* amplgjac  */
  coded[8]  = 0; /* amplgjacp */
  coded[9]  = 1; /* amplhl    */
  coded[10] = 0; /* amplhlp   */

  /* Check derivatives? */
  checkder = 0;

  /* SET SOME SOLVER ARGUMENTS */
  epsfeas   = 1.0e-08;
  epsopt    = 1.0e-08;

  efstain   = sqrt( epsfeas );
  eostain   = pow( epsopt, 1.5 );

  efacc     = sqrt( epsfeas );
  eoacc     = sqrt( epsopt );

  specfnm   = NULL;
  outputfnm = NULL;

  nvparam   = 0;

  /* GET THE OPTIONS AND STUB'S FILE NAME. */
  getOptions(argv,&jcnnzmax,&hnnzmax,&epsfeas,&epsopt,&efstain,&eostain,
	     &efacc,&eoacc,&outputfnm,&specfnm);

  if ( specfnm == NULL ) {
    specfnm = &emptyfnm;
  }

  if ( outputfnm == NULL ) {
    outputfnm = &emptyfnm;
  }

  inip(&n,&x,&l,&u,&m,&lambda,&equatn,&linear,&hnnzmax,&jcnnzmax);

  /* CALL ALGENCAN */
  c_algencan(&amplf,&amplg,&amplh,&amplc,&ampljac,&amplhc,&amplfc,
	     &amplgjac,&amplgjacp,&amplhl,&amplhlp,jcnnzmax,hnnzmax,
	     &epsfeas,&epsopt,&efstain,&eostain,&efacc,&eoacc,outputfnm,
	     specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded,
	     checkder,&f,&cnorm,&snorm,&nlpsupn,&inform);

  /* WRITE ADDITIONAL OUTPUT INFORMATION CODED BY THE USER */
  endp(n,&x,&l,&u,m,&lambda,&equatn,&linear);

  ASL_free( &asl );

  return 0;
}
