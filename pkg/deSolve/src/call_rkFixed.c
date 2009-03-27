/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* General RK Solver for methods with fixed step size                       */
/*==========================================================================*/

#include "rk_util.h"

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
  int fsal=0; // fixed step methods have no FSAL

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

  double  qerr  = REAL(getListElement(Method, "Qerr"))[0];
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
  int *ipar = NULL;

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
    lipar = 3; // in lsoda: 1;
    lrpar = 1;
  }
  ipar  = (int *) R_alloc(lipar, sizeof(int));

//  if (isDll ==1) {
    ipar[0] = nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];
    /* first Nout elements of rpar reserved for output variables
       other elements are set in R-function lsodx via argument *rpar* */
    // for (j = 0; j < nout; j++)         out[j] = 0.;                  //???
    // for (j = 0; j < LENGTH(Rpar); j++) out[nout+j] = REAL(Rpar)[j];  //???
//  }
  // end new testing code

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
  for (i = 0; i < nt * (neq + 1); i++) yout[i] = NA_REAL;

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
        derivs(Func, t + dt * cc[j], tmp, Parms, Rho, FF, out, j, neq, ipar, isDll);
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
    derivs(Func, t, tmp, Parms, Rho, FF, out, -1, neq, ipar, isDll);
    for (i = 0; i < nout; i++) {
      yout[j + nt * (1 + neq + i)] = out[i];
    }
  }
  // attach essential internal information (codes are compatible to lsoda)
  // ToDo: respect function evaluations due to external outputs
  setIstate(R_yout, R_istate, istate, it_tot, stage, fsal, qerr);


  // release R resources
  if (verbose) {
    Rprintf("Number of time steps it = %d, it_ext = %d, it_tot = %d\n", it, it_ext, it_tot);
    Rprintf("Maxsteps %d\n", maxsteps);
  }
  unprotect_all();
  return(R_yout);
}
