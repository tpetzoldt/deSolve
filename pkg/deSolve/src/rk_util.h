/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* Definitions and Utilities needed by Runge-Kutta Solvers                  */
/*==========================================================================*/


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

#include <R_ext/Applic.h>
#include <R_ext/Boolean.h>

#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif

#include "deSolve.h"

typedef void deriv_func(int *, double *, double *,double *, double *, int *);

typedef void init_func (void (*)(int *, double *));

void R_test_call(DllInfo *info);

void R_unload_test_call(DllInfo *info);

SEXP getvar(SEXP name, SEXP Rho);

SEXP getInputs(SEXP symbol, SEXP Rho);

SEXP getListElement(SEXP list, const char *str);

void matprod(int m, int n, int o, double* a, double* b, double* c);

double maxdiff(double *x, double *y, int n);

double maxerr(double *y1, double *y2, double* Atol, double* Rtol, int n);

void derivs(SEXP Func, double t, double* y, SEXP Parms, SEXP Rho,
	    double *ydot, double *yout, int j, int neq, int nout, int isDll);
	    
void denspar(double *FF, double *y0, double *y1, double dt, double *d,
  int neq, int stage, double *r);

void densout(double *r, double t0, double t, double dt, double* res, int neq);

void neville(double *xx, double *y, double tnew, double *ynew, int n, int ksig);

void shiftBuffer (double *x, int n, int k);

void initParms(SEXP Initfunc, SEXP Parms);

void setIstate(SEXP R_yout, SEXP R_istate, int *istate,
  int it_tot, int stage, int fsal, int qerr);
  

/* a reduced version of blas_matprod without NA checking */
static void blas_matprod1(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z) {
    const char *transa = "N", *transb = "N";
    int i;
    double one = 1.0, zero = 0.0;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
    	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}
