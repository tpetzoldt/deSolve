/* file lorenz.c */
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
//#include "deSolve.h" // defines a.o. histvar, histdvar, histtime


static double parms[3];
#define a parms[0]
#define b parms[1]
#define c parms[2]

static SEXP R_ret, R_T, R_nr;
static double *ret = NULL, *T = NULL;
static double tau = 1.22; // delay of DDE
static int *nr = NULL;

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


void inithist(int max, int maxlags, int solver, int nroot) {
  static void(*fun)(int, int, int, int) = NULL;
  if (fun == NULL)
    fun =  (void(*)(int, int, int, int))R_GetCCallable("deSolve", "inithist");

  return fun(max, maxlags, solver, nroot);
}  

void initglobal(int neq, int interpolMethod, int offset) {
  static void(*fun)(int, int, int) = NULL;
  if (fun == NULL)
    fun =  (void(*)(int, int, int))R_GetCCallable("deSolve", "initglobal");

  return fun(neq, interpolMethod, offset);
}  


void updatehistini(double t, double *y, double *dY, double *rwork, int *iwork) {
  static void(*fun)(double, double*, double*, double*, int*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double, double*, double*, double*, int*))R_GetCCallable("deSolve", "updatehistini");

  return fun(t, y, dY, rwork, iwork);
}


void updatehist(double t, double *y, double *dY, double *rwork, int *iwork) {
  static void(*fun)(double, double*, double*, double*, int*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double, double*, double*, double*, int*))R_GetCCallable("deSolve", "updatehist");

  return fun(t, y, dY, rwork, iwork);
}


/* initializer  */
void initmod(void (* odeparms)(int *, double *)) {
  int N = 3, histsize = 10000;
    odeparms(&N, parms);

    // neq, interpolmethod, offset
    //initglobal(3, 1, 3);
    //inithist(histsize, 10000, 1, 0);
}



/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{

    if (ip[0] < 1) error("nout should be at least 1");

    //PROTECT(R_ret = allocVector(REALSXP, 1));
    PROTECT(R_ret = NEW_NUMERIC(1));
    PROTECT(R_T   = NEW_NUMERIC(1));
    PROTECT(R_nr  = NEW_INTEGER(1)); // can also be vector
    
    T   = REAL(R_T);
    nr  = INTEGER(R_nr);

    *T = *t - tau;
    *nr = 0;

    if (*t > tau) {
      R_ret = getLagValue(R_T, R_nr);
      ret = REAL(R_ret);
      Rprintf("test %g %g %g \n", *t, y[0], *ret);
      yout[0] = *ret;
    } else {
      yout[0] = 0;
    }



    ydot[0] = a * y[0] + y[1] * y[2];
    ydot[1] = b * (y[1] - y[2]);
    ydot[2] = - y[0] * y[1] + c * y[1] - y[2];
    
    unprotect(3);
}
