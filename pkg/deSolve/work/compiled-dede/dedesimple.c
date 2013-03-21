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

void getlagvalue(double *T, int *nr, int N, double *yout) {
  static void(*fun)(double*, int*, int, double*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double*, int*, int, double*))R_GetCCallable("deSolve", "getlagvalue");
  return fun(T, nr, N, yout);
}

void getlagderiv(double *T, int *nr, int N, double *yout) {
  static void(*fun)(double*, int*, int, double*) = NULL;
  if (fun == NULL)
    fun =  (void(*)(double*, int*, int, double*))R_GetCCallable("deSolve", "getlagvalue");
  return fun(T, nr, N, yout);
}

/* Initializer  */
void initmod(void (* odeparms)(int *, double *)) {
  int N = 2;
  odeparms(&N, parms);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip) {



  if (ip[0] < 1) error("nout should be at least 1");

  double T = *t - tau;
  int nr[1] = {0};             // array
  double ytau[1] = {1.0};      // array


  if (*t > tau) {

    getlagvalue(&T, nr, 1, ytau);
    Rprintf("test %g %g %g \n", T, y[0], ytau[0]);

  }

  yout[0] = ytau[0];
  ydot[0] = k * ytau[0];

}
