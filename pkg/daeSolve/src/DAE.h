#include <R.h>
#include <Rdefines.h>
/* global variables */
SEXP Time, Y, YPRIME , Rin;
extern SEXP de_gparms;

typedef void deriv_func(int *, double *, double *,double *,double *, int *);
deriv_func * derfun;

typedef void res_func(double *, double *, double *, double*, double *,
                      int*, double *, int*);
res_func * res_fun;

typedef void init_func (void (*)(int *, double *));

void updatedeforc(double *);

/* DAE globals */
extern SEXP DAE_res_func;
extern SEXP DAE_jac_func;
extern SEXP DAE_envir;

/* utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);


/* declarations for initialisations;*/
void initParms(SEXP Initfunc, SEXP Parms);
void Initdeparms(int *, double *);
void Initdeforc(int *, double *);
void initOut(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar);

/* use in daspk */
int n_eq;
int mu;
int ml;
int nrowpd;

/* output in DLL globals */
int nout, ntot, isOut, lrpar, lipar, *ipar;
double *out;

/* KS globals for the forcings */
int initForcings(SEXP list);

SEXP getListElement(SEXP list, const char *str);

long int nforc;  /* the number of forcings */

/* Input data. three vectors:
  tmat, fmat: time, forcing function data value
  imat: index to start of each forcing function in tmat, fmat*/
double * tvec;
double * fvec;
int    * ivec;
int    fmethod;

/* for each forcing function: index to current position in tmat, fmat,
 current value, interpolation factor, current forcing time, next forcing time,
 max time (to be removed).....
*/
int    * findex;
double * intpol;
int    * maxindex;

double * forcings;
