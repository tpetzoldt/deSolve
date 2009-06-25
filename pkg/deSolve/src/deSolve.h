#include <R.h>
#include <Rdefines.h>
/* global variables */
SEXP Time, Y, YPRIME , Rin;
extern SEXP de_gparms;

typedef void deriv_func(int *, double *, double *,double *,double *, int *);
void updatedeforc(double *);
deriv_func * derfun;

typedef void res_func(double *, double *, double *, double*, double *,
                      int*, double *, int*);
res_func * res_fun;

typedef void init_func (void (*)(int *, double *));


/* vode globals */
extern SEXP vode_deriv_func;
extern SEXP vode_jac_func;
extern SEXP vode_envir;

/* lsoda globals */
extern SEXP odesolve_deriv_func;
extern SEXP odesolve_jac_func;
extern SEXP odesolve_jac_vec;
extern SEXP odesolve_root_func;
extern SEXP odesolve_envir;

/* daspk globals */
extern SEXP daspk_res_func;
extern SEXP daspk_jac_func;
extern SEXP daspk_psol_func;
extern SEXP daspk_envir;

/* utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

/* declarations for initideparms;*/
void Initdeparms(int *, double *);
void Initdeforc(int *, double *);

/* use in daspk */
long int n_eq;
long int mu;
long int ml;
long int nrowpd;


/* KS globals for the forcings */
/* the number of forcings */
long int nforc;
/* Input data. three vectors:
  tmat, fmat: time, forcing function data value
  imat: index to start of each forcing function in tmat, fmat*/
double * tvec;
double * fvec;
int    * ivec;

/* for each forcing function: index to current position in tmat, fmat,
 current value, interpolation factor, current forcing time, next forcing time,
 max time (to be removed).....
*/
int    * findex;
double * curval;
double * intpol;
double * curtime;
double * nexttime;
int    * maxindex;

double * forcings;
