#ifndef R_R_H
# include <R.h>
#endif


#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif


#include "deSolve.h"

// register native routines ---------------------------------------------------
//#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
//#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void unlock_solver();

/* .Call calls */
extern SEXP call_daspk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_DLL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_euler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_iteration(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_lsoda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_radau(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rk4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rkAuto(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rkFixed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rkImplicit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_zvode(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getLagDeriv(SEXP, SEXP);
extern SEXP getLagValue(SEXP, SEXP);
extern SEXP getTimestep();

/* Examples (manually added) */
//void initccl4(void (* odeparms)(int *, double *))
//void derivsccl4 (int *neq, double *t, double *y, double *ydot, double *out, int *ip)

extern void initccl4(SEXP);
extern void initparms(SEXP);
extern void initforcs(SEXP);
extern void derivsccl4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern void eventfun(SEXP, SEXP, SEXP);
extern void chemres(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern void scocpar(SEXP);
extern void scocforc(SEXP);
extern void scocder(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern void iniaqua(SEXP);
extern void initaqforc(SEXP);
extern void aquaphy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern void aquaphyforc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CMethodDef CEntries[] = {
    {"unlock_solver", (DL_FUNC) &unlock_solver, 0},
//
    {"initccl4",     (DL_FUNC) &initccl4,    1},
    {"initparms",    (DL_FUNC) &initparms,   1},
    {"initforcs",    (DL_FUNC) &initforcs,   1},
    {"eventfun",     (DL_FUNC) &eventfun,    3},
    {"derivsccl4",   (DL_FUNC) &derivsccl4,  6},
    {"chemres",      (DL_FUNC) &chemres,     8},
    {"scocpar",      (DL_FUNC) &scocpar,     1},
    {"scocforc",     (DL_FUNC) &scocforc,    1},
    {"scocder",      (DL_FUNC) &scocder,     6},
    {"iniaqua",      (DL_FUNC) &iniaqua,     1},
    {"initaqforc",   (DL_FUNC) &initaqforc,  1},
    {"aquaphy",      (DL_FUNC) &aquaphy,     6},
    {"aquaphyforc",  (DL_FUNC) &aquaphy,     6},
//
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"call_daspk",      (DL_FUNC) &call_daspk,      28},
    {"call_DLL",        (DL_FUNC) &call_DLL,        11},
    {"call_euler",      (DL_FUNC) &call_euler,      11},
    {"call_iteration",  (DL_FUNC) &call_iteration,  12},
    {"call_lsoda",      (DL_FUNC) &call_lsoda,      28},
    {"call_radau",      (DL_FUNC) &call_radau,      26},
    {"call_rk4",        (DL_FUNC) &call_rk4,        11},
    {"call_rkAuto",     (DL_FUNC) &call_rkAuto,     21},
    {"call_rkFixed",    (DL_FUNC) &call_rkFixed,    17},
    {"call_rkImplicit", (DL_FUNC) &call_rkImplicit, 17},
    {"call_zvode",      (DL_FUNC) &call_zvode,      21},
    {"getLagDeriv",     (DL_FUNC) &getLagDeriv,      2},
    {"getLagValue",     (DL_FUNC) &getLagValue,      2},
    {"getTimestep",     (DL_FUNC) &getTimestep,      0},
    {NULL, NULL, 0}
};



// C callable functions -------------------------------------------------------
SEXP get_deSolve_gparms(void);

void lagvalue(double T, int* nr, int N, double* ytau);

void lagderiv(double T, int* nr, int N, double* ytau);

double glob_timesteps[] = {0, 0};

// Initialization -------------------------------------------------------------
void R_init_deSolve(DllInfo *dll) {

  // thpe 2017-03-22, register entry points
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);

  // thpe xxxx, register C callable
  // R_RegisterCCallable("deSolve", "get_deSolve_gparms", (DL_FUNC) get_deSolve_gparms);

 // thpe: macro from package Matrix
 #define RREGDEF(name)  R_RegisterCCallable("deSolve", #name, (DL_FUNC) name)

  RREGDEF(get_deSolve_gparms);
  RREGDEF(lagvalue);
  RREGDEF(lagderiv);

  /* initialize global variables */
  timesteps = glob_timesteps;
}
