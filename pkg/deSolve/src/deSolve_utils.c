/* Define some global variables and functions that operate on some of them */
/* Patterned on code odesolve_utils.c from package odesolve */
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "deSolve.h"

/* some functions for keeping track of how many SEXPs 
 * 	are PROTECTed, and UNPROTECTing them in the case of a fortran stop.
 */
 
long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}

/* Globals :*/
SEXP odesolve_deriv_func;
SEXP odesolve_jac_func;
SEXP odesolve_jac_vec;
SEXP odesolve_root_func;
SEXP odesolve_envir;
SEXP odesolve_gparms;

SEXP daspk_res_func;
SEXP daspk_jac_func;
SEXP daspk_psol_func;
SEXP daspk_envir;

SEXP vode_deriv_func;
SEXP vode_jac_func;
SEXP vode_envir;
SEXP de_gparms;


/* Parameter initialisation function
note: forcing initialisation ftion is in forcings.c*/

void initParms(SEXP Initfunc, SEXP Parms) {

  if (inherits(Initfunc, "NativeSymbol"))  {
    init_func *initializer;

    PROTECT(de_gparms = Parms);     incr_N_Protect();
    initializer = (init_func *) R_ExternalPtrAddr(Initfunc);
    initializer(Initdeparms);
  }

}


void Initdeparms(int *N, double *parms)
{
  int i, Nparms;

  Nparms = LENGTH(de_gparms);
  if ((*N) != Nparms)
    {
      PROBLEM "Confusion over the length of parms"
      ERROR;
    } 
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(de_gparms)[i];
    }
}
  
SEXP get_deSolve_gparms(void)
{
  return de_gparms;
}

/* extracting elements from a list */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
	 if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	 }
   return elmt;
}

/* output initialisation function */
/* The output:
    out and ipar are used to pass output variables (number set by nout)
    followed by other input (e.g. forcing functions) provided
    by R-arguments rpar, ipar
    ipar[0]: number of output variables, ipar[1]: length of rpar,
    ipar[2]: length of ipar */
void initOut(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar) {

  int j;
  nout   = INTEGER(nOut)[0];    /* number of output variables */
  if (isDll)  /* function is a dll */
  {
   if (nout > 0) isOut = 1;
   ntot  = neq + nout;          /* length of yout */
   lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
   lipar = 3 + LENGTH(Ipar);    /* length of ipar */

  } else                              /* function is not a dll */
  {
   isOut = 0;
   ntot = neq;
   lipar = 1;
   lrpar = 1;
  }

   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
      other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout; j++)        out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }

}
