/* compile within R with system("R CMD SHLIB <filename>.cpp") */
/* Example adapted from lsoda documentation */

/* A typical C trick to get readable names for parameters */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>     // AS_NUMERIC ..

static double parms[7];

#define a parms[0]
#define b parms[1]
#define c parms[2]
#define d parms[3]
#define e parms[4]
#define f parms[5]
#define g parms[6]

#define S y[0]
#define P y[1]
#define K y[2]

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

void R_test_call(DllInfo *info) {
  /* Register routines, allocate resources. */
  Rprintf("test_call DLL loaded\n");
} 
        
void R_unload_test_call(DllInfo *info) {
  /* Release resources. */
  Rprintf("test_call DLL unloaded\n");
}

 
// ----- The model itself ----------------------------------------------------

void ilotka(void (* odeparms)(int *, double *)) {
    int N = 7;
    odeparms(&N, parms);
    Rprintf("lotka parameter vector initialized\n");
}

/* Derivatives */
void dlotka(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
    // sanity checks
    if (ip[0] < 2) error("nout should be at least 2");
    // an "external" signal
    double import = 0;
    if (10 <= *t && *t <=11) import = 0.2;

    // derivatives
    ydot[0] = import - b * S * P + g * K;
    ydot[1] = c * S * P  - d * K * P;
    ydot[2] = e * P * K  - f * K;


    // additional outputs
    yout[0] = a;   //inp[1];
    yout[1] = ydot[0]; //inp[2];
}
