/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2  */


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>     // AS_NUMERIC ...

#include <R_ext/Applic.h> // for dgemm
#include <R_ext/Boolean.h>

#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif

#include "deSolve.h"

/* --------- from call_lsoda.c ------------ */
/* give name to data types */
typedef void deriv_func(int *, double *, double *,double *, double *, int *);
/*
typedef void root_func (int *, double *, double *,int *, double *);
typedef void jac_func  (int *, double *, double *, int *,
		                    int *, double *, int *, double *, int *);
typedef void jac_vec   (int *, double *, double *, int *,
		                    int *, int *, double *, double *, int *);
*/
typedef void init_func (void (*)(int *, double *));

/* from deSolve_utils.c */

/*
SEXP de_gparms;

void Initdeparms(int *N, double *parms)
{
  int i, Nparms;

  Nparms = LENGTH(de_gparms);
  if ((*N) != Nparms)
    {
      PROBLEM "Confusion over the length of parms"
      ERROR;
    }
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(de_gparms)[i];
    }
}

SEXP get_deSolve_gparms(void)
{
  return de_gparms;
}
*/

/* Alternative method:
  - define only global pointer and
  - use dynamic memory allocation later
*/

//static double *inp = NULL;

/* some functions for keeping track of how many SEXPs
 * 	are PROTECTed, and UNPROTECTing them in the case of a fortran stop.
 */

/*

long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}
*/

void R_test_call(DllInfo *info) {
  /* Register routines, allocate resources. */
  Rprintf("test_call DLL loaded\n");
}

void R_unload_test_call(DllInfo *info) {
  /* Release resources. */
  Rprintf("test_call DLL unloaded\n");
}


// -------- test getvar from environment ------------------------------------
SEXP getvar(SEXP name, SEXP Rho) {
  SEXP ans;
  if(!isString(name) || length(name) != 1)
    error("name is not a single string");
  if(!isEnvironment(Rho))
    error("Rho should be an environment");
  ans = findVar(install(CHAR(STRING_ELT(name, 0))), Rho);
  return(ans);
}

SEXP getInputs(SEXP symbol, SEXP Rho) {
  if(!isEnvironment(Rho)) error("Rho should be an environment");
  return(getvar(symbol, Rho));
}

//--- test getvar from list -------------------------------------------------
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
   if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
     elmt = VECTOR_ELT(list, i);
     break;
   }
  return elmt;
}

// the following snippet is copied from array.c

/*
static void blas_matprod(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
    const char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0;
    LDOUBLE sum;
    //Rboolean have_na = FALSE;
    int have_na = FALSE;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	// Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
	// The test is only O(n) here
	//
	for (i = 0; i < nrx*ncx; i++)
	    if (ISNAN(x[i])) {have_na = TRUE; break;}
	if (!have_na)
	    for (i = 0; i < nry*ncy; i++)
		if (ISNAN(y[i])) {have_na = TRUE; break;}
	if (have_na) {
	    for (i = 0; i < nrx; i++)
		for (k = 0; k < ncy; k++) {
		    sum = 0.0;
		    for (j = 0; j < ncx; j++)
			sum += x[i + j * nrx] * y[j + k * nry];
		    z[i + k * nrx] = sum;
		}
	} else
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
    } else // zero-extent operations should return zeroes
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}
*/

// a reduced version without NA checking
static void blas_matprod1(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
    const char *transa = "N", *transb = "N";
    int i;
    double one = 1.0, zero = 0.0;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
    	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}

// Simple Matrix Multiplikation
void matprod(int m, int n, int o, double* a, double* b, double* c) {
  int i, j, k;
  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
    	c[i + m * j] = 0;
    	for (k = 0; k < n; k++) {
    	  c[i + m * j] += a[i + m * k] * b[k + n * j];
    	}
    }
  }
}

double maxdiff(double *x, double *y, int n) {
  double d = 0.0;
  for (int i = 0; i < n; i++) d = fmax(d, fabs(x[i] - y[i]));
  return(d);
}

double maxerr(double *y1, double *y2, double* Atol, double* Rtol, int n) {
  double err = 0, serr = 0, scal, delta;
  for (int i = 0; i < n; i++) {
    scal  = Atol[i] +  fmax(fabs(y1[i]), fabs(y2[i])) * Rtol[i]; // min??
    delta = fabs(y2[i] - y1[i]);
    err   = fmax(err, delta / scal); // one of these two lines
    serr  = err + pow(delta/scal, 2.0);    // one of these two
  }
  err = sqrt(serr); //euclidean norm
  return(err);
}

/*==========================================================================*/
/*   CALL TO THE MODEL FUNCTION                                             */
/*==========================================================================*/
void derivs(SEXP Func, double t, double* y, SEXP Parms, SEXP Rho,
	    double *ydot, double *yout, int j, int neq, int nout, int isDll) {
  SEXP Val, R_fcall;
  SEXP R_t;
  SEXP R_y;
  int i = 0;
  double *yy;
  double ytmp[neq];

  if (isDll)  {
    /************************************************************************/
    /*   Function is a DLL function                                         */
    /************************************************************************/
    deriv_func *cderivs;
    cderivs = (deriv_func *) R_ExternalPtrAddr(Func);
    //  ToDo: int nout --> int* iout
    //  for comparison: lsoda_derivs(N, t, y, ydot, yout, iout)
    cderivs (&neq, &t, y, ytmp, yout, &nout);
    if (j >= 0)
      for (i = 0; i < neq; i++)  ydot[i + neq * j] = ytmp[i];
  } else {
    /************************************************************************/
    /* Function is an R function                                            */
    /************************************************************************/
    PROTECT(R_t = ScalarReal(t)); incr_N_Protect();
    PROTECT(R_y = allocVector(REALSXP, neq)); incr_N_Protect();
    yy = REAL(R_y);
    for (i=0; i< neq; i++) yy[i] = y[i];

    PROTECT(R_fcall = lang4(Func, R_t, R_y, Parms)); incr_N_Protect();
    PROTECT(Val = eval(R_fcall, Rho)); incr_N_Protect();
    // extract the states of list "val"
    if (j >= 0)
      for (i = 0; i < neq; i++)  ydot[i + neq * j] = REAL(VECTOR_ELT(Val, 0))[i];
    // extract outputs from second list element
    if (j < 0)
      for (i = 0; i < nout; i++)  yout[i] = REAL(VECTOR_ELT(Val, 1))[i];
    my_unprotect(4);
  }
}

// Interpolation function for methods with "dense output"
void denspar(double *FF, double *y0, double *y1, double dt, double *d,
  int neq, int stage, double *r) {
  double ydiff, bspl;
  int i, j;
  for (i=0; i< neq; i++) {
   r[i]           = y0[i];
   ydiff          = y1[i] - y0[i];
   r[i + neq]     = ydiff;
   bspl           = dt * FF[i] - ydiff;
   r[i + 2 * neq] = bspl;
   r[i + 3 * neq] = ydiff - dt * FF[i + (stage - 1) * neq] - bspl;
   r[i + 4 * neq] = 0;
   for (j=0; j < stage; j++)
     r[i + 4 * neq] = r[i + 4 * neq] + d[j] * FF[i + j * neq];
     r[i + 4 * neq] = r[i + 4 * neq] * dt;
  }
}

void densout(double *r, double t0, double t, double dt, double* res, int neq) {
  double s  = (t - t0) / dt;
  double s1 = 1.0 - s;
  for (int i = 0; i < neq; i++)
    res[i] = r[i] + s * (r[i +     neq] + s1 * (r[i + 2 * neq]
                  + s * (r[i + 3 * neq] + s1 * (r[i + 4 * neq]))));
}

/* Polynomial interpolation
   ksig: number of signals
   n:    number of knots per signal
   x[0 .. n-1]:          vector of x values
   y[0 .. n-1, 0 .. ksig] array  of y values

   ToDo: check if ringbuffer is faster; rewrite eventually

*/
void neville(double *xx, double *y, double tnew, double *ynew, int n, int ksig) {
  int i, j, k;
  double x[n];
  double yy[n * ksig]; // temporary workspace
  double tscal = xx[n-1] - xx[0];
  double t = tnew / tscal;
  for (i = 0; i < n; i++) x[i] = xx[i]/tscal;

  for (i=0; i < n * ksig; i++) yy[i] = y[i];

  for (k = 0; k < ksig; k++) {
    for (j = 1; j < n; j++)
      for (i = n - 1; i >= j; i--) {
        yy[i + k * n] = ((t - x[i - j]) * yy[i + k * n]
          - (t - x[i]) * yy[i - 1 + k * n]) / (x[i] - x[i - j]);
      }
    ynew[k] = yy[n - 1 + k * n];
  }
}

void shiftBuffer (double *x, int n, int k) {
  // n = rows, k=columns
  for (int i = 0; i < (n - 1); i++)
    for (int j = 0; j < k; j++)
      x[i + j * n] = x[i + 1 + j * n];
}


void initParms(SEXP Initfunc, SEXP Parms) {
 if (inherits(Initfunc, "NativeSymbol"))  {
    PROTECT(de_gparms = Parms); incr_N_Protect();
    init_func *initializer;
    initializer = (init_func *) R_ExternalPtrAddr(Initfunc);
    initializer(Initdeparms);
   }
   }

 void setIstate(SEXP R_yout, SEXP R_istate, int *istate,
   int it_tot, int stage, int fsal, int qerr) {
   // note that indices are "C = R - 1"
   istate[11] = it_tot;                  // number of steps
   istate[12] = it_tot * (stage - fsal); // number of function evaluations
   istate[14] = qerr;                    // order of the method
   setAttrib(R_yout, install("istate"), R_istate);
 }


SEXP call_rkAuto(SEXP Xstart, SEXP Times, SEXP Func, SEXP Initfunc,
  SEXP Parms, SEXP Nout, SEXP Rho,
  SEXP Rtol, SEXP Atol, SEXP Tcrit, SEXP Verbose,
  SEXP Hmin, SEXP hmax, SEXP Hini, SEXP Rpar, SEXP Ipar,
  SEXP Method, SEXP Maxsteps) {

  /**  Initialization **/
  init_N_Protect();

  double *tt = NULL, *xs = NULL;

  double *y,  *f,  *Fj, *tmp, *FF, *rr;
  SEXP  R_yout;
  double *y0,  *y1,  *y2,  *dy1,  *dy2, *out, *yout;

  double err=0, dtnew=0, t, dt, t_ext, tmax;

  SEXP R_FSAL;
  int fsal=0; // assume no FSAL

  int i = 0, j=0, j1=0, k, it=0, it_tot=0, it_ext=0, nt = 0, neq=0;
  int accept = 0;
  int one=1;

  /**************************************************************************/
  /****** Processing of Arguments                                      ******/
  /**************************************************************************/
  int lAtol = LENGTH(Atol);
  double *atol = (double *) R_alloc((int) lAtol, sizeof(double));

  int lRtol = LENGTH(Rtol);
  double *rtol = (double *) R_alloc((int) lRtol, sizeof(double));

  for (j = 0; j < lRtol; j++) rtol[j] = REAL(Rtol)[j];
  for (j = 0; j < lAtol; j++) atol[j] = REAL(Atol)[j];

  double  tcrit = REAL(Tcrit)[0];
  double  hmin  = REAL(Hmin)[0];
  double  Hmax  = REAL(hmax)[0];
  double  hini  = REAL(Hini)[0];
  int  maxsteps = (int)REAL(Maxsteps)[0];
  int  nout     = (int)REAL(Nout)[0]; // number of external outputs is func is in a DLL
  int  verbose  = (int)REAL(Verbose)[0];

  int stage     = (int)REAL(getListElement(Method, "stage"))[0];

  SEXP R_A, R_B1, R_B2, R_C, R_D;
  double  *A, *bb1, *bb2=NULL, *cc=NULL, *dd=NULL;

  PROTECT(R_A = getListElement(Method, "A")); incr_N_Protect();
  A = REAL(R_A);

  PROTECT(R_B1 = getListElement(Method, "b1")); incr_N_Protect();
  bb1 = REAL(R_B1);

  PROTECT(R_B2 = getListElement(Method, "b2")); incr_N_Protect();
  if (length(R_B2)) bb2 = REAL(R_B2);

  PROTECT(R_C = getListElement(Method, "c")); incr_N_Protect();
  if (length(R_C)) cc = REAL(R_C);

  PROTECT(R_D = getListElement(Method, "d")); incr_N_Protect();
  if (length(R_D)) dd = REAL(R_D);

  double  qerr  = REAL(getListElement(Method, "Qerr"))[0];
  PROTECT(R_FSAL = getListElement(Method, "FSAL")); incr_N_Protect();
  if (length(R_FSAL)) fsal = INTEGER(R_FSAL)[0];

  PROTECT(Times = AS_NUMERIC(Times)); incr_N_Protect();
  tt = NUMERIC_POINTER(Times);
  nt = length(Times);

  PROTECT(Xstart = AS_NUMERIC(Xstart)); incr_N_Protect();
  xs  = NUMERIC_POINTER(Xstart);
  neq = length(Xstart);

  /**************************************************************************/
  /****** DLL, ipar, rpar (to be compatible with lsoda)                ******/
  /**************************************************************************/
  int isDll = 0;
  int ntot = 0;
  int isOut = 0; //?? do I need this?
  int lrpar= 0, lipar = 0;
  int *ipar;

  // testing code from lsoda
  if (inherits(Func, "NativeSymbol")) { /* function is a dll */
    isDll = 1;
    if (nout > 0) isOut = 1;
    ntot  = neq + nout;          /* length of yout */
    lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
    lipar = 3 + LENGTH(Ipar);    /* length of ipar */

  } else {                              /* function is not a dll */
    isDll = 0;
    isOut = 0;
    ntot = neq;
    lipar = 1;
    lrpar = 1;
  }
  ipar  = (int *) R_alloc(lipar, sizeof(int));

  if (isDll ==1) {
    ipar[0] = nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];
    /* first Nout elements of rpar reserved for output variables
       other elements are set in R-function lsodx via argument *rpar* */
    // for (j = 0; j < nout; j++)         out[j] = 0.;                  //???
    // for (j = 0; j < LENGTH(Rpar); j++) out[nout+j] = REAL(Rpar)[j];  //???
  }
  // end new testing code

  /**************************************************************************/
  /****** Allocation of Workspace                                      ******/
  /**************************************************************************/
  y0  =  (double *) R_alloc(neq, sizeof(double));
  y1  =  (double *) R_alloc(neq, sizeof(double));
  y2  =  (double *) R_alloc(neq, sizeof(double));
  dy1 =  (double *) R_alloc(neq, sizeof(double));
  dy2 =  (double *) R_alloc(neq, sizeof(double));
  f   =  (double *) R_alloc(neq, sizeof(double));
  y   =  (double *) R_alloc(neq, sizeof(double));
  Fj  =  (double *) R_alloc(neq, sizeof(double));
  tmp =  (double *) R_alloc(neq, sizeof(double));
  FF  =  (double *) R_alloc(neq * stage, sizeof(double));
  rr  =  (double *) R_alloc(neq * 5, sizeof(double));

  out  =  (double *) R_alloc(nout, sizeof(double));

  // matrix for polynomial interpolation
  int nknots = 4;  // 3rd order polynomials
  int iknots = 0;  // counter for knotes buffer
  double *yknots;
  yknots = (double *) R_alloc(neq * (nknots + 1), sizeof(double));

  // matrix for holding states and external outputs
  PROTECT(R_yout = allocMatrix(REALSXP, nt, neq + nout + 1)); incr_N_Protect();
  yout = REAL(R_yout);
  // initialize outputs with NA first
  for (i = 0; i < nt * (neq + nout + 1); i++) yout[i] = NA_REAL;

  // attribute that stores state information, similar to lsoda
  SEXP R_istate;
  int *istate;
  PROTECT(R_istate = allocVector(INTSXP, 22)); incr_N_Protect();
  istate = INTEGER(R_istate);
  istate[0] = 0; // assume succesful return
  for (i = 0; i < 22; i++) istate[i] = 0;

  //PROTECT(RSTATE = allocVector(REALSXP, 5));incr_N_Protect();
  //for (k = 0;k<5;k++) REAL(RSTATE)[k] = rwork[k+10];


  /**************************************************************************/
  /****** Initialization of Parameters (for DLL functions)             ******/
  /**************************************************************************/

  initParms(Initfunc, Parms);

  /**************************************************************************/
  /****** Initialization of Integration Loop                           ******/
  /**************************************************************************/

  yout[0]   = tt[0];              // initial time
  yknots[0] = tt[0];              // for polynomial interpolation
  for (i = 0; i < neq; i++) {
    y0[i]        = xs[i];         // initial values
    yout[(i + 1) * nt] = y0[i];   // output array
    yknots[iknots + nknots * (i + 1)] = xs[i]; // for polynomials
  }
  iknots++;

  t = tt[0];
  tmax = fmax(tt[nt], tcrit);
  dt = fmin(Hmax, hini);
  Hmax = fmin(Hmax, tmax - t);

 // Initialization of work arrays (to be on the safe side, remove this later)
  for (i = 0; i < neq; i++)  {
    y1[i] = 0;
    y2[i] = 0;
    Fj[i] = 0;
    for (j= 0; j < stage; j++)  {
      FF[i + j * neq] = 0;
    }
  }
  /**************************************************************************/
  /****** Main Loop                                                    ******/
  /**************************************************************************/
  it     = 1; // step counter; zero element is initial state
  it_ext = 0; // counter for external time step (dense output)
  it_tot = 0; // total number of time steps

  do { //<-------------- ??
    //Rprintf("it, t, dt, %d  %e  %e\n", it, t, dt);
    /******  save former results of last step if the method allows this
            (first same as last)                                       ******/
    if (fsal && accept){
      j1 = 1;
      for (i = 0; i < neq; i++) FF[i] = FF[i + neq * (stage - 1)];
    } else {
      j1 = 0;
    }
    /******  Prepare Coefficients from Butcher table ******/
    for (j = j1; j < stage; j++) {
      for(i = 0; i < neq; i++) Fj[i] = 0;
        k = 0;
        while(k < j) {
          for(i = 0; i < neq; i++)
            Fj[i] = Fj[i] + A[j + stage * k] * FF[i + neq * k] * dt;
          k++;
        }
        for (int i = 0; i < neq; i++) {
          tmp[i] = Fj[i] + y0[i];
        }
        /******  Compute Derivatives ******/
        // pass option to avoid unnecessary copying in derivs
        derivs(Func, t + dt * cc[j], tmp, Parms, Rho, FF, out, j, neq, nout, isDll);
    }

    /************************************************************************/
    /* Estimation of new values                                             */
    /************************************************************************/

    // -- alternative 1: hand-made
    //matprod(neq, stage, one, FF, bb1, dy1);
    //matprod(neq, stage, one, FF, bb2, dy2);

    // -- alternative 2: use BLAS
    //blas_matprod(FF, neq, stage, bb1, stage, one, dy1);
    //blas_matprod(FF, neq, stage, bb2, stage, one, dy2);

    // -- alternative 3: use BLAS with reduced error checking
    blas_matprod1(FF, neq, stage, bb1, stage, one, dy1);
    blas_matprod1(FF, neq, stage, bb2, stage, one, dy2);

    it_tot++; // count total number of time steps
    for (i = 0; i < neq; i++) {
      y1[i] = y0[i] +  dt * dy1[i];
      y2[i] = y0[i] +  dt * dy2[i];
    }

    /************************************************************************/
    /****** stepsize adjustment                                        ******/
    /************************************************************************/
    err = maxerr(y1, y2, atol, rtol, neq);

    dtnew = dt;
    accept =TRUE;
    if (err < 1.0e-20) {  // this will probably never occur
      accept = TRUE;
      dtnew = Hmax;
      //Rprintf("dtnew %e -- =0=   ", dtnew);
    } else if (err < 1.0) {
      accept = TRUE;
      dtnew = fmin(Hmax, dt * 0.9 * pow(err, -1.0/qerr)); // 1/qerr
      //Rprintf("dtnew %e  (++)   \n", dtnew);
    } else if (err > 1.0) {
      accept = FALSE;
      dtnew = dt * fmax(0.9 * pow(err, -1.0/qerr), 0.2); // 1/qerr
      //Rprintf("2  dtnew %e  (--)   \n", dtnew);
    }

    if (dtnew < hmin) {     // R: dt !!
      accept=TRUE;
      if (verbose) Rprintf("warning, h < Hmin\n"); // remove this later ...
      istate[0] = -2;
      dtnew = hmin;
    }
    /************************************************************************/
    /****** Interpolation and Data Storage                             ******/
    /************************************************************************/
    if (accept) {
      /**********************************************************************/
      /* case A) "Dense Output": built-in polynomial interpolation          */
      /* available for certain rk formulae, e.g. for rk45dp7                */
      /**********************************************************************/
      if (dd) { // i.e. if dd is not NULL
        denspar(FF, y0, y1, dt, dd, neq, stage, rr);
        t_ext = tt[it_ext];
        while (t_ext <= t + dt) { // <= ??
          densout(rr, t, t_ext, dt, tmp, neq);
          // store outputs
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            for (i = 0; i < neq; i++)
              yout[it_ext + nt * (1 + i)] = tmp[i];
          }
          if(it_ext < nt) t_ext = tt[++it_ext]; else break;
        }
      /**********************************************************************/
      /* case B) "Neville-Aitken-Interpolation" for integrators             */
      /* without dense output                                               */
      /**********************************************************************/
      } else {
        // (1) collect number "nknots" of knots in advanve
        yknots[iknots] = t + dt;   // time in first column
        for (i = 0; i < neq; i++) yknots[iknots + nknots * (1 + i)] = y2[i];
        if (iknots < (nknots - 1)) {
          iknots++;
        } else {
         // (2) do polynomial interpolation
         t_ext = tt[it_ext];
         while (t_ext <= t + dt) { // <= ??
          neville(yknots, &yknots[nknots], t_ext, tmp, nknots, neq);
          // (3) store outputs
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            for (i = 0; i < neq; i++)
              yout[it_ext + nt * (1 + i)] = tmp[i];
          }
          if(it_ext < nt) t_ext = tt[++it_ext]; else break;
         }
         shiftBuffer(yknots, nknots, neq + 1);
        }
      }
      /**********************************************************************/
      /* next time step                                                     */
      /**********************************************************************/
      t = t + dt;
      it++;
      for (i=0; i < neq; i++) y0[i] = y2[i];
    } // else rejected time step
    dt = fmin(dtnew, tmax - t);
    if (it_ext > nt) {
      Rprintf("error in rk_solvers.c - call_rkauto: output buffer overflow\n");
      break;
    }
    if (it_tot > maxsteps) {
      if (verbose) Rprintf("Max. number of steps exceeded\n");
      istate[0] = -1;
      break;
    }
  } while (t < tmax); // end of rk main loop

  /**************************************************************************/
  /* call derivs again to get external outputs                              */
  /**************************************************************************/
  // j = -1 suppresses internal copying
  for (int j = 0; j < nt; j++) {
    t = yout[j];
    for (i = 0; i < neq; i++) tmp[i] = yout[j + nt * (1 + i)];
    derivs(Func, t, tmp, Parms, Rho, FF, out, -1, neq, nout, isDll);
    for (i = 0; i < nout; i++) {
      yout[j + nt * (1 + neq + i)] = out[i];
    }
  }
  // attach essential internal information (codes are compatible to lsoda)
  // ToDo: respect function evaluations due to external outputs
  setIstate(R_yout, R_istate, istate, it_tot, stage, fsal, qerr);

  // release R resources
  if (verbose) Rprintf("Number of time steps it = %d, it_ext = %d, it_tot = %d\n",
    it, it_ext, it_tot);
  //Rprintf("Maxsteps %d\n", maxsteps);
  unprotect_all();
  //init_N_Protect();
  return(R_yout);
}

//----------------------------------------------------------------------------
SEXP call_rkFixed(SEXP Xstart, SEXP Times, SEXP Func, SEXP Initfunc,
  SEXP Parms, SEXP Nout, SEXP Rho,
  SEXP Tcrit, SEXP Verbose, SEXP Hini, SEXP Rpar, SEXP Ipar,
  SEXP Method, SEXP Maxsteps) {

  /**  Initialization **/
  init_N_Protect();

  double *tt = NULL, *xs = NULL;

  double *y,  *f,  *Fj, *tmp, *FF, *rr;
  SEXP  R_yout;
  double *y0,  *y1, *dy1, *out, *yout;

  double t, dt, t_ext, tmax;

  int i = 0, j=0, j1=0, k, it=0, it_tot=0, it_ext=0, nt = 0, neq=0;
  int one=1;

  /**************************************************************************/
  /****** Processing of Arguments                                      ******/
  /**************************************************************************/
  double  tcrit = REAL(Tcrit)[0];
  double  hini  = REAL(Hini)[0];
  int  maxsteps = (int)REAL(Maxsteps)[0];
  int  nout  = (int)REAL(Nout)[0]; // number of external outputs is Func is in a DLL
  int verbose = (int)REAL(Verbose)[0];

  int stage = (int)REAL(getListElement(Method, "stage"))[0];

  SEXP R_A, R_B1, R_C;
  double  *A, *bb1, *cc=NULL;

  PROTECT(R_A = getListElement(Method, "A")); incr_N_Protect();
  A = REAL(R_A);

  PROTECT(R_B1 = getListElement(Method, "b1")); incr_N_Protect();
  bb1 = REAL(R_B1);

  PROTECT(R_C = getListElement(Method, "c")); incr_N_Protect();
  if (length(R_C)) cc = REAL(R_C);

  PROTECT(Times = AS_NUMERIC(Times)); incr_N_Protect();
  tt = NUMERIC_POINTER(Times);
  nt = length(Times);

  PROTECT(Xstart = AS_NUMERIC(Xstart)); incr_N_Protect();
  xs  = NUMERIC_POINTER(Xstart);
  neq = length(Xstart);

  /**************************************************************************/
  /****** DLL, ipar, rpar (to be compatible with lsoda)                ******/
  /**************************************************************************/
  // provisional code; copy from call_rkauto
  int isDll = 0;
  if (inherits(Func, "NativeSymbol")) isDll = 1;

  /**************************************************************************/
  /****** Allocation of Workspace                                      ******/
  /**************************************************************************/
  y0  =  (double *) R_alloc(neq, sizeof(double));
  y1  =  (double *) R_alloc(neq, sizeof(double));
  dy1 =  (double *) R_alloc(neq, sizeof(double));
  f   =  (double *) R_alloc(neq, sizeof(double));
  y   =  (double *) R_alloc(neq, sizeof(double));
  Fj  =  (double *) R_alloc(neq, sizeof(double));
  tmp =  (double *) R_alloc(neq, sizeof(double));
  FF  =  (double *) R_alloc(neq * stage, sizeof(double));
  rr  =  (double *) R_alloc(neq * 5, sizeof(double));

  out  =  (double *) R_alloc(nout, sizeof(double));

  // matrix for polynomial interpolation
  int nknots = 4;  // 3rd order polynomials
  int iknots = 0;  // counter for knotes buffer
  double *yknots;
  yknots = (double *) R_alloc(neq * (nknots + 1), sizeof(double));


  // matrix for holding the outputs
  PROTECT(R_yout = allocMatrix(REALSXP, nt, neq + nout + 1)); incr_N_Protect();
  yout = REAL(R_yout);
  // initialize outputs with NA first
  for (i=0; i< nt*(neq+1); i++) yout[i] = NA_REAL;

  // attribute that stores state information, similar to lsoda
  SEXP R_istate;
  int *istate;
  PROTECT(R_istate = allocVector(INTSXP, 22)); incr_N_Protect();
  istate = INTEGER(R_istate);
  istate[0] = 0; // assume succesful return
  for (i = 0; i < 22; i++) istate[i] = 0;

  //PROTECT(RSTATE = allocVector(REALSXP, 5));incr_N_Protect();
  //for (k = 0;k<5;k++) REAL(RSTATE)[k] = rwork[k+10];

  /**************************************************************************/
  /****** Initialization of Parameters (for DLL functions)             ******/
  /**************************************************************************/

  initParms(Initfunc, Parms);

  /**************************************************************************/
  /****** Initialization of Integration Loop                           ******/
  /**************************************************************************/

  yout[0]   = tt[0];              // initial time
  yknots[0] = tt[0];              // for polynomial interpolation
  for (i = 0; i < neq; i++) {
    y0[i]        = xs[i];         // initial values
    yout[(i + 1) * nt] = y0[i];   // output array
    yknots[iknots + nknots * (i + 1)] = xs[i]; // for polynomials
  }
  iknots++;

  t = tt[0];                   // t    <- min(Times)
  tmax = fmax(tt[nt], tcrit);   // tmax <- max(Times, Tcrit)

  // Initialization of work arrays (to be on the safe side, remove this later)
  for (i = 0; i < neq; i++)  {
    y1[i] = 0;
    //y2[i] = 0;
    Fj[i] = 0;
    for (j= 0; j < stage; j++)  {
      FF[i + j * neq] = 0;
    }
  }
  /**************************************************************************/
  /****** Main Loop                                                    ******/
  /**************************************************************************/
  it     = 1; // step counter; zero element is initial state
  it_ext = 0; // counter for external time step (dense output)
  it_tot = 0; // total number of time steps

  do {
    /* select time step (possibly irregular) */
    if (hini > 0.0)
      dt = hini;
    else
      dt = tt[it] - tt[it-1];

    /******  Prepare Coefficients from Butcher table ******/
    for (j = j1; j < stage; j++) {
      for(i = 0; i < neq; i++) Fj[i] = 0;
        k = 0;
        while(k < j) {
          for(i = 0; i < neq; i++)
            Fj[i] = Fj[i] + A[j + stage * k] * FF[i + neq * k] * dt;
          k++;
        }
        for (int i = 0; i < neq; i++) {
          tmp[i] = Fj[i] + y0[i];
        }
        /******  Compute Derivatives ******/
        derivs(Func, t + dt * cc[j], tmp, Parms, Rho, FF, out, j, neq, nout, isDll);
    }

    /************************************************************************/
    /* Estimation of new values                                             */
    /************************************************************************/
    // use BLAS with reduced error checking
    blas_matprod1(FF, neq, stage, bb1, stage, one, dy1);

    it_tot++; // count total number of time steps
    for (i = 0; i < neq; i++) {
      y1[i] = y0[i] +  dt * dy1[i];
    }

    /************************************************************************/
    /****** Interpolation and Data Storage                             ******/
    /************************************************************************/
    // (1) collect number "nknots" of knots in advanve
    yknots[iknots] = t + dt;   // time in first column
    for (i = 0; i < neq; i++) yknots[iknots + nknots * (1 + i)] = y1[i];
    if (iknots < (nknots - 1)) {
      iknots++;
    } else {
     // (2) do polynomial interpolation
     t_ext = tt[it_ext];
     while (t_ext <= t + dt) { // <= ??
      neville(yknots, &yknots[nknots], t_ext, tmp, nknots, neq);
      // (3) store outputs
      if (it_ext < nt) {
        yout[it_ext] = t_ext;
        for (i = 0; i < neq; i++)
          yout[it_ext + nt * (1 + i)] = tmp[i];
      }
      if(it_ext < nt) t_ext = tt[++it_ext]; else break;
     }
     shiftBuffer(yknots, nknots, neq + 1);
    }
    /**********************************************************************/
    /* next time step                                                     */
    /**********************************************************************/
    t = t + dt;
    it++;
    for (i=0; i < neq; i++) y0[i] = y1[i];
    if (it_ext > nt) {
      Rprintf("error in rk_solvers.c - call_rk4auto: output buffer overflow\n");
      break;
    }
    if (it_tot > maxsteps) {
      if (verbose) Rprintf("Max. number of steps exceeded\n");
      break;
    }
  } while (t < tmax); // end of rk main loop

  /**************************************************************************/
  /* call derivs again to get external outputs                              */
  /**************************************************************************/
  // j = -1 suppresses unnecessary internal copying
  for (int j = 0; j < nt; j++) {
    t = yout[j];
    for (i = 0; i < neq; i++) tmp[i] = yout[j + nt * (1 + i)];
    derivs(Func, t, tmp, Parms, Rho, FF, out, -1, neq, nout, isDll);
    //Rprintf("%d %e %e \n", j, out[0], out[1]);
    for (i = 0; i < nout; i++) {
      yout[j + nt * (1 + neq + i)] = out[i];
    }
  }
  // attach essential internal information (codes are compatible to lsoda)
  // ToDo: respect function evaluations due to external outputs
  // ToDo: fsal; it_tot; qerr
  //setIstate(R_yout, R_istate, istate, it_tot, stage, fsal, qerr);


  // release R resources
  if (verbose) {
    Rprintf("Number of time steps it = %d, it_ext = %d, it_tot = %d\n", it, it_ext, it_tot);
    Rprintf("Maxsteps %d\n", maxsteps);
  }
  unprotect_all();
  //init_N_Protect();
  return(R_yout);
}

/*==========================================================================*/
/*  rk4 Fixed Step Integrator                                               */
/*    (special version for speed comparison with the general solution)      */
/*==========================================================================*/
SEXP call_rk4(SEXP Xstart, SEXP Times, SEXP Func, SEXP Initfunc,
	      SEXP Parms, SEXP Nout, SEXP Rho, SEXP Verbose,
	      SEXP Rpar, SEXP Ipar) {

  /**  Initialization **/
  init_N_Protect();

  double *tt = NULL, *xs = NULL;
  double *tmp, *FF, *out;

  SEXP  R_y, R_f, R_f1, R_f2, R_f3, R_f4;
  double *y,  *f,  *f1,  *f2,  *f3,  *f4;

  SEXP  R_y0, R_yout;
  double *y0,  *yout;


  double t, dt;
  int i = 0, it=0, nt = 0, neq=0;

  /**************************************************************************/
  /****** Check Arguments and Convert to C types                       ******/
  /**************************************************************************/

  PROTECT(Times = AS_NUMERIC(Times)); incr_N_Protect();
  tt = NUMERIC_POINTER(Times);
  nt = length(Times);

  PROTECT(Xstart = AS_NUMERIC(Xstart)); incr_N_Protect();
  xs  = NUMERIC_POINTER(Xstart);
  neq = length(Xstart);

  tmp =  (double *) R_alloc(neq, sizeof(double));
  FF  =  (double *) R_alloc(neq, sizeof(double));

  int  nout  = (int)REAL(Nout)[0]; // n of external outputs if func is in a DLL
  out  =  (double *) R_alloc(nout, sizeof(double));

  int verbose = (int)REAL(Verbose)[0];

  /**************************************************************************/
  /****** DLL, ipar, rpar (to be compatible with lsoda)                ******/
  /**************************************************************************/
  // provisional code; copy from call_rkauto
  int isDll = 0;
  if (inherits(Func, "NativeSymbol")) isDll = 1;

  /**************************************************************************/
  /****** Allocation of Workspace                                      ******/
  /**************************************************************************/

  PROTECT(R_y0 = allocVector(REALSXP, neq)); incr_N_Protect();
  PROTECT(R_f  = allocVector(REALSXP, neq)); incr_N_Protect();
  PROTECT(R_y  = allocVector(REALSXP, neq)); incr_N_Protect();
  PROTECT(R_f1 = allocVector(REALSXP, neq)); incr_N_Protect();
  PROTECT(R_f2 = allocVector(REALSXP, neq)); incr_N_Protect();
  PROTECT(R_f3 = allocVector(REALSXP, neq)); incr_N_Protect();
  PROTECT(R_f4 = allocVector(REALSXP, neq)); incr_N_Protect();
  y0 = REAL(R_y0);
  f  = REAL(R_f);
  y  = REAL(R_y);
  f1 = REAL(R_f1);
  f2 = REAL(R_f2);
  f3 = REAL(R_f3);
  f4 = REAL(R_f4);

  // matrix for holding the outputs
  PROTECT(R_yout = allocMatrix(REALSXP, nt, neq + nout + 1)); incr_N_Protect();
  yout = REAL(R_yout);

  // attribute that stores state information, similar to lsoda
  SEXP R_istate;
  int *istate;
  PROTECT(R_istate = allocVector(INTSXP, 22)); incr_N_Protect();
  istate = INTEGER(R_istate);
  istate[0] = 0; // assume succesful return
  for (i = 0; i < 22; i++) istate[i] = 0;

  //PROTECT(RSTATE = allocVector(REALSXP, 5));incr_N_Protect();
  //for (k = 0;k<5;k++) REAL(RSTATE)[k] = rwork[k+10];

  /**************************************************************************/
  /****** Initialization of Parameters (for DLL functions)             ******/
  /**************************************************************************/

  initParms(Initfunc, Parms);

  /**************************************************************************/
  /****** Initialization of Integration Loop                           ******/
  /**************************************************************************/

  yout[0] = tt[0]; //initial time
  for (i = 0; i < neq; i++) {
    y0[i]              = xs[i];
    yout[(i + 1) * nt] = y0[i];      // <--- check this
  }

  /**************************************************************************/
  /****** Main Loop                                                    ******/
  /**************************************************************************/
  for (it = 0; it < nt - 1; it++) {
    t = tt[it];
    dt = tt[it + 1] - t;
    if (verbose)
      Rprintf("Time steps = %d / %d time = %e\n", it + 1, nt, t);
    derivs(Func, t, y0, Parms, Rho, f1, out, 0, neq, nout, isDll);
    for (i = 0; i < neq; i++) {
      f1[i] = dt * f1[i];
      f[i]  = y0[i] + 0.5 * f1[i];
    }
    derivs(Func, t + 0.5*dt, f, Parms, Rho, f2, out, 0, neq, nout, isDll);
    for (i = 0; i < neq; i++) {
      f2[i] = dt * f2[i];
      f[i]  = y0[i] + 0.5 * f2[i];
    }
    derivs(Func, t + 0.5*dt, f, Parms, Rho, f3, out, 0, neq, nout, isDll);
    for (i = 0; i < neq; i++) {
      f3[i] = dt * f3[i];
      f[i] = y0[i] + f3[i];
    }
    derivs(Func, t + dt, f, Parms, Rho, f4, out, 0, neq, nout, isDll);
    for (i = 0; i < neq; i++) {
      f4[i] = dt * f4[i];
    }
    // Final computation of y
    for (i = 0; i < neq; i++) {
      f[i]  = (f1[i] + 2.0 * f2[i] + 2.0 * f3[i] + f4[i]) / 6.0;
      y[i]  = y0[i] + f[i];
      y0[i] = y[i]; // next time step
    }
    // store outputs
    if (it < nt) {
      yout[it + 1] = t + dt;
      for (i = 0; i < neq; i++) yout[it + 1 + nt * (1 + i)] = y[i];
    }
  } // end of rk main loop

  /**************************************************************************/
  /* call derivs again to get external outputs                              */
  /**************************************************************************/
  // j= -1 suppresses internal copying
  for (int j = 0; j < nt; j++) {
    t = yout[j];
    for (i = 0; i < neq; i++) tmp[i] = yout[j + nt * (1 + i)];
    derivs(Func, t, tmp, Parms, Rho, FF, out, -1, neq, nout, isDll);
    for (i = 0; i < nout; i++) {
      yout[j + nt * (1 + neq + i)] = out[i];
    }
  }
  // attach essential internal information (codes are compatible to lsoda)
  // ToDo: respect function evaluations due to external outputs
  // ToDo: it_tot; stage; fsal; qerr
  //setIstate(R_yout, R_istate, istate, it_tot, stage, fsal, qerr);

  // release R resources
  unprotect_all();
  //init_N_Protect();
  return(R_yout);
}
