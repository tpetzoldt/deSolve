#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "deSolve.h"

SEXP get_deSolve_gparms(void);

void getlagvalue(double *T, int* nr, int N, double* yout);

void getlagderiv(double *T, int* nr, int N, double* yout);

double glob_timesteps[] = {0, 0};

void R_init_deSolve(DllInfo *info) {
// R_RegisterCCallable("deSolve", "get_deSolve_gparms", (DL_FUNC) get_deSolve_gparms);

 // thpe: macro from package Matrix
 #define RREGDEF(name)  R_RegisterCCallable("deSolve", #name, (DL_FUNC) name)

  RREGDEF(get_deSolve_gparms);

  RREGDEF(getlagvalue);
  RREGDEF(getlagderiv);

  /* initialize global variables */
  timesteps = glob_timesteps;
}
