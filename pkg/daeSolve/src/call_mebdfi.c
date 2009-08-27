/* Patterned on code from package odesolve */
#include <time.h>
#include <string.h>
#include "DAE.h"

/* definition of the call to the fortran function mebdfi - in file mebdfi.f*/
void F77_NAME(mebdfi)(
		     int *, double *, double *, double *, double *, double *, double *,
		     int *, int *, int*, double *, int *, int *, int *, int*, int *,
         double *,  double *,  double *, int *,

		     void (*)(double *, double *, double *, double *, double *,
                  double *, int *),             /* jac*/
		     void (*)(double *, double *, double *, double *, double *, int *,
                  double *, int *),             /* func */
         int *) ;

static void forc_DAE (double *t, double *y, double *yprime, double *cj,
                       double *delta, int *ires, double *yout, int *iout)
{
  updatedeforc(t);
  res_fun(t, y, yprime, cj, delta, ires, yout, iout);
}


/* interface between fortran function call and R function  */

static void DAE_res (double *t, double *y, double *yprime, double *cj,
                       double *delta, int *ires, double *yout, int *iout)
{                             
  int i;
  SEXP R_fcall, ans;

  REAL(Time)[0] = *t;
  for (i = 0; i < n_eq; i++)
    {
      REAL(Y)[i] = y[i];
      REAL (YPRIME)[i] = yprime[i];
    }
  PROTECT(R_fcall = lang4(DAE_res_func,Time, Y, YPRIME));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, DAE_envir));                  incr_N_Protect();

  for (i = 0; i < n_eq; i++)  	delta[i] = REAL(ans)[i];
  my_unprotect(2);
//  Rprintf(" %g", *t)
}

/* interface between fortran call to jacobian and R function */

static void DAE_jac (double *t, double *y, double *yprime,
                       double *pd,  double *cj, double *RPAR, int *IPAR)
{
  int i;
  SEXP R_fcall, ans;

  REAL(Rin)[0] = *t;
  REAL(Rin)[1] = *cj;  

  for (i = 0; i < n_eq; i++)
    {
      REAL(Y)[i] = y[i];
      REAL (YPRIME)[i] = yprime[i];      
    }
  PROTECT(R_fcall = lang4(DAE_jac_func, Rin, Y, YPRIME));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, DAE_envir));                 incr_N_Protect();
  for (i = 0; i < n_eq * nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* give name to data types */

typedef void jac_func(double *, double *, double *, double *, double *,
                      double *, int *);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_mebdfi(SEXP y, SEXP yprime, SEXP times, SEXP res, SEXP parms,
		SEXP rtol, SEXP atol, SEXP itol, SEXP rho, SEXP Tcrit, SEXP Hini,
    SEXP Maxord, SEXP maxIt, SEXP nind, SEXP jacfunc, SEXP initfunc,
    SEXP verbose, SEXP Mf, SEXP Mbnd, SEXP Liw, SEXP Lrw,
    SEXP nOut, SEXP Rpar, SEXP Ipar, SEXP flist)
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP   yout, yout2=NULL, dyout=NULL, ISTATE, RWORK;
  int    i, j, k, nt, latol, lrtol, lrw, liw, isDll;
  int    isForcing , Itol, *mbnd, mf, maxord, isOut;
  double *xytmp,  *xdytmp, *rwork, tin, tout, *Atol, *Rtol, tcrit, hini;
  double *delta=NULL, cj;
  int    idid, *iwork, mflag, ires, ierr;

  res_func  *Resfun;
  jac_func  *jac=NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    

  init_N_Protect();


  n_eq = LENGTH(y);
  nt = LENGTH(times);
  mflag = INTEGER(verbose)[0];        
  nout  = INTEGER(nOut)[0];
  ntot  = n_eq;

  mf = INTEGER(Mf)[0];
  maxord = INTEGER(Maxord)[0];

  tcrit = REAL(Tcrit)[0];
  hini = REAL(Hini)[0];
  ierr = 0;

/* function is a dll ?*/
  if (inherits(res, "NativeSymbol"))
    isDll = 1;
  else
    isDll = 0;

  isOut = 0;
  if (isDll == 0 && nout > 0) isOut =1;
  else if (isDll == 1 ) ntot = ntot+ nout;

/* initialise output vectors ... */
  if (isDll==1)  { /* function is a dll */
    lrpar = nout + LENGTH(Rpar);       /* length of rpar */
    lipar = 3    + LENGTH(Ipar);       /* length of ipar */
  } else  {                             /* function is not a dll */
    lipar = 1;
    lrpar = 1;
  }

  out   = (double *) R_alloc(lrpar, sizeof(double));
  ipar  = (int *)    R_alloc(lipar, sizeof(int));

  if (isDll ==1)  {
    ipar[0] = nout;          /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
      other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout;        j++) out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
  } else {
    for (j = 0; j < lrpar;       j++) out[j] = 0.;
    for (j = 0; j < lipar;       j++) ipar[j] = 0.;
  }

  /* copies of all variables that will be changed in the FORTRAN subroutine */
  xytmp = (double *) R_alloc(n_eq, sizeof(double));
   for (j = 0; j < n_eq; j++) xytmp[j] = REAL(y)[j];

  xdytmp = (double *) R_alloc(n_eq, sizeof(double));
   for (j = 0; j < n_eq; j++) xdytmp[j] = REAL(yprime)[j];

  /* tolerance */
  latol = LENGTH(atol);
  Atol  = (double *) R_alloc((int) latol, sizeof(double));
    for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol  = (double *) R_alloc((int) lrtol, sizeof(double));
    for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];

  Itol  = INTEGER(itol)[0];

  /* work arrays */
  liw = INTEGER(Liw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));
  for (j = 0; j<3; j++) iwork[j] = INTEGER(nind)[j];
  for (j = 3; j<liw; j++) iwork[j] = 0;
  iwork[13] = INTEGER(maxIt)[0];

  lrw =  INTEGER(Lrw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
  for (j = 0; j < lrw; j++) rwork[j] = 0.;

  mbnd  = (int *) R_alloc(4, sizeof(int));
  for (j = 0; j<4; j++) mbnd[j] = INTEGER(Mbnd)[j];

  /* initialise global variables... */

  PROTECT(Time = NEW_NUMERIC(1));                    incr_N_Protect();
  PROTECT(Rin  = NEW_NUMERIC(2));                    incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,n_eq));            incr_N_Protect();
  PROTECT(YPRIME = allocVector(REALSXP,n_eq));       incr_N_Protect();
  PROTECT(yout = allocMatrix(REALSXP,ntot+1,nt));    incr_N_Protect();
  if (isOut == 1) {
    PROTECT(dyout = allocMatrix(REALSXP,n_eq+1,nt));    incr_N_Protect();
  }
  
  /**************************************************************************/
  /****** Initialization of Parameters and Forcings (DLL functions)    ******/
  /**************************************************************************/
  initParms(initfunc, parms);
    //  error("till here");

  isForcing = initForcings(flist);

 /* pointers to functions res and jac, passed to the FORTRAN subroutine */

  if (isDll == 1)  {       /* DLL address passed to fortran */
      Resfun = (res_func *) R_ExternalPtrAddr(res);

      delta = (double *) R_alloc(n_eq, sizeof(double));
      for (j = 0; j < n_eq; j++) delta[j] = 0.;

      if(isForcing==1) {
        res_fun = (res_func *) R_ExternalPtrAddr(res);
        Resfun = (res_func *) forc_DAE;
      }

    } else {
      /* interface function between fortran and R passed to Fortran*/     
      Resfun = (res_func *) DAE_res;
      /* needed to communicate with R */      
      PROTECT(DAE_res_func = res);               incr_N_Protect();
      PROTECT(DAE_envir = rho);                  incr_N_Protect();

    }
  if (!isNull(jacfunc))
    {
      if (inherits(jacfunc,"NativeSymbol"))
	      jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
      else  {
	      jac = (jac_func *) DAE_jac;
	      PROTECT(DAE_jac_func = jacfunc);         incr_N_Protect();
	    }
    }
/*                      #### initial time step ####                           */
  idid = 1;
  REAL(yout)[0] = REAL(times)[0];
  for (j = 0; j < n_eq; j++)
      REAL(yout)[j+1] = REAL(y)[j];

  if (nout>0)  {
     tin = REAL(times)[0];

	   if (isDll == 1) {
       Resfun (&tin, xytmp, xdytmp, &cj, delta, &ires, out, ipar) ;
       for (j = 0; j < nout; j++)
	      REAL(yout)[n_eq + j + 1] = out[j];
     }
	   else for (j = 0; j < n_eq; j++)
	      REAL(dyout)[j] = xdytmp[j];
  }
/*                     ####   main time loop   ####                           */

  for (i = 0; i < nt-1; i++)
  {
      tin = REAL(times)[i];
      tout = REAL(times)[i+1];

	   //Rprintf(" tin, hini, itol, mf, lrw, mbnd %g %g %i %i %i %i %i %i %i\n", tin, hini, Itol, mf, lrw, mbnd[0], mbnd[1],mbnd[2] ,mbnd[3] );

      F77_CALL(mebdfi)(&n_eq, &tin, &hini, xytmp, xdytmp, &tout, &tcrit,
        &mf, &idid, &lrw, rwork, &liw, iwork, mbnd, &maxord,
	      &Itol, Rtol, Atol, out, ipar,  jac, Resfun, &ierr);
//	   Rprintf(" tin, xdytmp %g %g %g %g %g %g\n", tin, xdytmp[0], xdytmp[1],xdytmp[2] ,xdytmp[3], xdytmp[4] );
//	   Rprintf(" hini, xytmp %g %g %g %g %g %g\n", hini, xytmp[0], xytmp[1],xytmp[2] ,xytmp[3], xytmp[4] );
//   error("here");
     if (idid == 1) {
        idid = 0;
        F77_CALL(mebdfi)(&n_eq, &tin, &hini, xytmp, xdytmp, &tout, &tcrit,
         &mf, &idid, &lrw, rwork, &liw, iwork, mbnd, &maxord,
	       &Itol, Rtol, Atol, out, ipar,  jac, Resfun, &ierr);
      }

	    if (idid == -1)   {
	      warning("the integration failed to pass the error test, even after reducing h by factor 1e10");
	    }   else    if (idid == -2)   {
	      warning("the integration failed by repeated error test failures or by a test on rtol/atol. too much accuracy requested");
      }   else    if (idid == -3)   {
	      warning("the integration failed to achieve corrector convergence, even after reducing h by factor 1e10");
      }  else    if (idid == -4)    {
       warning("illegal values of input parameters - see printed message");
      }  else    if (idid == -5)    {
       warning("idid was -1 on input, but tout was not beyond t");
      }  else    if (idid == -6)    {
       warning("maximum number of integration steps exceeded");
      }  else    if (idid == -7)   {
       warning("stepsize is too small, < sqrt(uround)/100");
      }  else    if (idid == -11)   {
       warning("insufficient real workspace for the integration");
      }  else    if (idid == -12)   {
       warning("insufficient integer workspace for the integration");

      }
  // Rprintf(" i, tin, tout, ntot %i, %g, %g, %i\n", i, tin, tout, ntot);

 	  REAL(yout)[(i+1)*(ntot+1)] = tin;
	  for (j = 0; j < n_eq; j++)
	    REAL(yout)[(i+1)*(ntot + 1) + j + 1] = xytmp[j];

    if (nout>0) {
	    if (isDll == 1) {
        Resfun (&tin, xytmp, xdytmp, &cj, delta, &ires, out, ipar) ;
       for (j = 0; j < nout; j++)
	      REAL(yout)[(i+1)*(ntot + 1) + n_eq + j + 1] = out[j];
      }
	    else for (j = 0; j < n_eq; j++)
	      REAL(dyout)[(i+1)*(n_eq) + j + 1] = xdytmp[j];
    }
/**/
/*                    ####  an error occurred   ####                          */                     
    if (tin < tout) {
     if (idid >=0) idid = -10;
	   warning("Returning early from mebdfi  Results are accurate, as far as they go\n");

	/* redimension yout */
   	PROTECT(yout2 = allocMatrix(REALSXP,ntot+1,(i+2)));incr_N_Protect();
  	for (k = 0; k < i+2; k++)
	   for (j = 0; j < ntot+1; j++)
	    REAL(yout2)[k*(ntot+1) + j] = REAL(yout)[k*(ntot+1) + j];
	   break;
    }
  }    /* end main time loop */


/*     error("tillhere");
                ####   returning output   ####                           */

  PROTECT(ISTATE = allocVector(INTSXP, 15));incr_N_Protect();
  for (k = 0;k<14;k++) INTEGER(ISTATE)[k+1] = iwork[k];
  INTEGER(ISTATE)[0] = idid;


  PROTECT(RWORK = allocVector(REALSXP, 1));incr_N_Protect();
  REAL(RWORK)[0] = rwork[1];

  if (idid >= 0)
    {
      setAttrib(yout, install("istate"), ISTATE);
      setAttrib(yout, install("rstate"), RWORK);    
      if (isOut ==1) setAttrib(yout, install("dy"), dyout);
    }
  else
    {
      setAttrib(yout2, install("istate"), ISTATE);
      setAttrib(yout2, install("rstate"), RWORK);   
      if (isOut ==1) setAttrib(yout2, install("dy"), dyout);
    }
//

/*                       ####   termination   ####                            */       
  unprotect_all();
  if (idid >= 0)
    return(yout);
  else
    return(yout2);
}

