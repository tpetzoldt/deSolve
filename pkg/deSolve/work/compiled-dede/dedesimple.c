/* File dedesimple.c */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

static double parms[2];
#define tau parms[0]
#define k parms[1]

/* Interface to dede utility functions in package deSolve */
SEXP getLagValue(SEXP T, SEXP nr) {
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL)
    fun =  (SEXP(*)(SEXP, SEXP))R_GetCCallable("deSolve", "getLagValue");
  return fun(T, nr);
}  

SEXP getLagDeriv(SEXP T, SEXP nr) {
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL)
    fun =  (SEXP(*)(SEXP, SEXP))R_GetCCallable("deSolve", "getLagDeriv");
  return fun(T, nr);
}  


/* Initializer  */
void initmod(void (* odeparms)(int *, double *)) {
  int N = 2;
  odeparms(&N, parms);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip) {

  SEXP R_T, R_nr;
  double  *T   = NULL;
  int     *nr  = NULL;

  double ytau = 1.0;

  if (ip[0] < 1) error("nout should be at least 1");

  PROTECT(R_T   = NEW_NUMERIC(1));
  PROTECT(R_nr  = NEW_INTEGER(1));

  T  = REAL(R_T);
  nr = INTEGER(R_nr);

  *T = *t - tau;
  *nr = 0;

  if (*t > tau) {
    ytau = *REAL(getLagValue(R_T, R_nr));
    //Rprintf("test %g %g %g \n", *t, y[0], ytau);
  }

  yout[0] = ytau;

  ydot[0] = k * ytau;

  unprotect(2);
}
