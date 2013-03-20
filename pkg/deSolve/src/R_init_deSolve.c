#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "deSolve.h"

SEXP get_deSolve_gparms(void);

SEXP getLagValue(SEXP T, SEXP nr);

SEXP getLagDeriv(SEXP T, SEXP nr);

void initglobal(int neq, int interpolMethod, int offset); 


double glob_timesteps[] = {0, 0};

void R_init_deSolve(DllInfo *info) {
  R_RegisterCCallable("deSolve", "get_deSolve_gparms", (DL_FUNC) get_deSolve_gparms);

  // thpe: macro from package Matrix
#define RREGDEF(name)  R_RegisterCCallable("deSolve", #name, (DL_FUNC) name)

  RREGDEF(inithist);          //deSolve.h
  RREGDEF(updatehistini);   //deSolve.h
  RREGDEF(updatehist);      //deSolve.h
  RREGDEF(initglobal);      // thpe: for testing. rename if really necessary
  RREGDEF(getLagValue);
  RREGDEF(getLagDeriv);

  /* initialize global variables */
  timesteps = glob_timesteps;
}
