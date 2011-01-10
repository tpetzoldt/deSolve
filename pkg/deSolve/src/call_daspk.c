#include <time.h>
#include <string.h>
#include "deSolve.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Differential algebraic equation solver daspk.
   
   The C-wrappers that provide the interface between FORTRAN codes and R-code 
   are: C_res_func   : interface with R-code "res", passes function residuals  
        C_out        : interface with R-code "res", passes output variables  
        C_daejac_func: interface with R-code "jacres", passes jacobian
  
   C_forc_dae provides the interface between the residual function specified in
   a DLL and daspk, in case there are forcing functions.
   
  
   changes since 1.4
   karline: version 1.5: added forcing functions in DLL
   karline: version 1.6: added events
   karline: version 1.7: added time lags -> delay differential equations
            improving names

   to do: implement psolfunc
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* definition of the call to the FORTRAN function ddaspk - in file ddaspk.f*/
void F77_NAME(ddaspk)(void (*)(double *, double *, double *, double*,
                               double *, int*, double *, int*),
		     int *, double *, double *, double *, double *, 
		     int *,double *, double *,  int *,  double *,  int *, 
		     int *, int *, double *, int *,
		     void (*)/*(double *, double *, double *, double *, double *, double *, int *)*/,
		     void (*)(int *, double *, double *, double *, double *, double *, 
                  double *, double *, double *, int *, double *, double *, 
                        int *, double *, int *));   

static void C_forc_dae (double *t, double *y, double *yprime, double *cj,
                       double *delta, int *ires, double *yout, int *iout)
{
  updatedeforc(t);
  DLL_res_func(t, y, yprime, cj, delta, ires, yout, iout);
}

static void C_psol_func (int *neq, double *t, double *y, double *yprime, 
                        double *savr, double *wk, double *cj, double* wght,
                        double *wp, int *iwp, double *b, double *eplin, 
                        int *ierr, double *RPAR, int *IPAR)
{
/* not yet implemented */
}

/* interface between FORTRAN residual function call and R function  */

static void C_res_func (double *t, double *y, double *yprime, double *cj, 
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
  PROTECT(R_fcall = lang4(R_res_func,Time, Y, YPRIME));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));                  incr_N_Protect();

  for (i = 0; i < n_eq; i++)  	delta[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* deriv output function  */

static void C_out (int *nout, double *t, double *y, 
                       double *yprime, double *yout)
{
  int i;
  SEXP R_fcall, ans;

  REAL(Time)[0] = *t;
  for (i = 0; i < n_eq; i++)  
    {
      REAL(Y)[i] = y[i];
      REAL (YPRIME)[i] = yprime[i];      
    }
     
  PROTECT(R_fcall = lang4(R_res_func,Time, Y, YPRIME));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));                  incr_N_Protect();

  for (i = 0; i < *nout; i++) yout[i] = REAL(ans)[i + n_eq];

  my_unprotect(2);
}      

/* interface between FORTRAN call to jacobian and R function */

static void C_daejac_func (double *t, double *y, double *yprime, 
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
  PROTECT(R_fcall = lang4(R_daejac_func, Rin, Y, YPRIME));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));                    incr_N_Protect();
  for (i = 0; i < n_eq * nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* give name to data types */

typedef void C_daejac_func_type(double *, double *, double *, double *, double *,
                      double *, int *);
typedef void C_psol_func_type(int *, double *, double *, double *, double *,
                      double *, double *, double *, double *, int*, double *,
                      double *, int*, double *, int*);
typedef void C_kryljac_func_type(double *, int *, int *, double *, double *,
                          double *, double *, double *,
           double *, double *, double *, double *, int*, int*, double *, int*);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_daspk(SEXP y, SEXP yprime, SEXP times, SEXP resfunc, SEXP parms, 
		SEXP rtol, SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc, 
		SEXP psolfunc, SEXP verbose, SEXP info, SEXP iWork, SEXP rWork,  
    SEXP nOut, SEXP maxIt, SEXP bu, SEXP bd, SEXP nRowpd, SEXP Rpar,
    SEXP Ipar, SEXP flist, SEXP elag, SEXP eventfunc, SEXP elist)
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

  int    j, nt, ny, repcount, latol, lrtol, lrw, liw, isDll;
  int    maxit, isForcing, isEvent, islag, istate;
  double *xytmp,  *xdytmp, tin, tout, *Atol, *Rtol;
  double *delta=NULL, cj;
  int    *Info,  ninfo, idid, mflag, ires;
  int    *iwork, it, ntot= 0, nout;   
  double *rwork;
  

  /* pointers to functions passed to FORTRAN */
  C_res_func_type     *res_func;
  C_daejac_func_type  *daejac_func=NULL;
  C_psol_func_type    *psol_func=NULL;
  C_kryljac_func_type *kryljac_func=NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

  lock_solver(); /* prevent nested call of solvers that have global variables */

/*                      #### initialisation ####                              */    

  //init_N_Protect();
  long int old_N_Protect = save_N_Protected();  

  ny   = LENGTH(y);  
  n_eq = ny;                          /* n_eq is a global variable */
  nt = LENGTH(times);  
  mflag = INTEGER(verbose)[0];        

  ninfo=LENGTH(info);
  nrowpd = INTEGER(nRowpd)[0];  
  maxit = INTEGER(maxIt)[0];
  
/* function is a dll ?*/
  if (inherits(resfunc, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

  initOutC(isDll, &nout, &ntot, n_eq, nOut, Rpar, Ipar); 

  /* copies of all variables that will be changed in the FORTRAN subroutine */
  Info  = (int *) R_alloc(ninfo,sizeof(int));
   for (j = 0; j < ninfo; j++) Info[j] = INTEGER(info)[j];  
  
  xytmp = (double *) R_alloc(n_eq, sizeof(double));
   for (j = 0; j < n_eq; j++) xytmp[j] = REAL(y)[j];

  xdytmp = (double *) R_alloc(n_eq, sizeof(double));
   for (j = 0; j < n_eq; j++) xdytmp[j] = REAL(yprime)[j];

  latol = LENGTH(atol);
  Atol  = (double *) R_alloc((int) latol, sizeof(double));
    for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol  = (double *) R_alloc((int) lrtol, sizeof(double));
    for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];
  
  liw = LENGTH(iWork);
  iwork = (int *) R_alloc(liw, sizeof(int));   
    for (j = 0; j < liw; j++) iwork[j] = INTEGER(iWork)[j];  

  lrw = LENGTH(rWork);
  rwork = (double *) R_alloc(lrw, sizeof(double));
    for (j = 0; j < lrw; j++) rwork[j] = REAL(rWork)[j];
    
  //timesteps = (double *) R_alloc(2, sizeof(double));
  for (j = 0; j < 2; j++) timesteps[j] = 0.;

  /**************************************************************************/
  /****** Initialization of globals, Parameters and Forcings (DLLs)    ******/
  /**************************************************************************/
  initdaeglobals(nt, ntot);
  initParms(initfunc, parms);
  isForcing = initForcings(flist);
  isEvent = initEvents(elist, eventfunc);
  islag = initLags(elag, 0, 0);
 
 /* pointers to functions res_func, psol_func and daejac_func, 
    passed to the FORTRAN subroutine */

  if (isDll == 1)  {       /* DLL address passed to FORTRAN */
      res_func = (C_res_func_type *) R_ExternalPtrAddr(resfunc);
      /* no need to communicate with R - but output variables set here */            
      delta = (double *) R_alloc(n_eq, sizeof(double));
      for (j = 0; j < n_eq; j++) delta[j] = 0.;

      if(isForcing==1) {
        DLL_res_func = (C_res_func_type *) R_ExternalPtrAddr(resfunc);
        res_func = (C_res_func_type *) C_forc_dae;
      }

    } else {
      /* interface function between FORTRAN and R passed to FORTRAN */
      res_func = (C_res_func_type *) C_res_func;
      /* needed to communicate with R */      
      R_res_func = resfunc; 
      R_envir = rho; 

    }
  if (!isNull(jacfunc))
    {
      if (inherits(jacfunc,"NativeSymbol"))
     	{
     	if (Info[11] ==0) {        /*ordinary jac*/
	      daejac_func = (C_daejac_func_type *) R_ExternalPtrAddr(jacfunc);
	      } else {                /*krylov*/
	      kryljac_func = (C_kryljac_func_type *) R_ExternalPtrAddr(jacfunc);
	      }
	    }
      else  {
	    R_daejac_func = jacfunc;
	    daejac_func = C_daejac_func;
	    }
    }
  if (!isNull(psolfunc))
    {
      if (inherits(psolfunc,"NativeSymbol"))
     	{
	    psol_func = (C_psol_func_type *) R_ExternalPtrAddr(psolfunc);
	    }
      else  {
	    R_psol_func = psolfunc;
	    psol_func = C_psol_func;
	    }
    }

/*                      #### initial time step ####                           */    
  idid = 1;
  REAL(YOUT)[0] = REAL(times)[0];
  for (j = 0; j < n_eq; j++)
      REAL(YOUT)[j+1] = REAL(y)[j];

  if (islag == 1) updatehistini(REAL(times)[0], xytmp, xdytmp, rwork, iwork);
    
  if (nout>0)
    {
     tin = REAL(times)[0];

	   if (isDll == 1) res_func (&tin, xytmp, xdytmp, &cj, delta, &ires, out, ipar) ;
	   else C_out(&nout,&tin,xytmp,xdytmp,out);
	      for (j = 0; j < nout; j++)
	       REAL(YOUT)[j + n_eq + 1] = out[j]; 
    }
               
/*                     ####   main time loop   ####                           */    
               
  for (it = 0; it < nt-1; it++)
  {
      tin = REAL(times)[it];
      tout = REAL(times)[it+1];
    if (isEvent) {
      istate =  2;
      updateevent(&tin, xytmp, &istate);
      if (istate  == 1) Info[0] = 0;
      Info[3] = 1;
      rwork[0] = tout;
    }

     repcount = 0;
     do  /* iterations in case maxsteps>500* or in case islag */
	   {
     	if (Info[11] ==0) {        /*ordinary jac*/
	       F77_CALL(ddaspk) (res_func, &ny, &tin, xytmp, xdytmp, &tout,
			   Info, Rtol, Atol, &idid, 
			   rwork, &lrw, iwork, &liw, out, ipar, daejac_func, psol_func);

	      } else {                /* krylov - not yet used */
      	 F77_CALL(ddaspk) (res_func, &ny, &tin, xytmp, xdytmp, &tout,
			   Info, Rtol, Atol, &idid, 
			   rwork, &lrw, iwork, &liw, out, ipar, kryljac_func, psol_func);
        }
    /* in case timestep is asked for... */    
    timesteps [0] = rwork[10];
    timesteps [1] = rwork[11];
  
    if (islag == 1) updatehist(tin, xytmp, xdytmp, rwork, iwork);    
        
	  repcount ++;
	  if (idid == -1) 
      {Info[0]=1;
       } else     if (idid == -2)   {
	      warning("Excessive precision requested.  scale up `rtol' and `atol' e.g. by the factor %g\n",10.0);
       Info[0]=1;          
	      repcount=maxit+2;
	    }   else    if (idid == -3)   {
       warning("Error term became zero for some i: pure relative error control (ATOL(i)=0.0) for a variable which is now vanished");
       repcount=maxit+2;
      }   else    if (idid == -5)   {
	      warning("jacfun routine failed with the Krylov method"); 
        repcount = maxit+2;     
      }   else    if (idid == -6)   {
       warning("repeated error test failures on a step - singularity ?");
        repcount = maxit+2;     
      }  else    if (idid == -7)    {
       warning("repeated convergence test failures on a step - inaccurate Jacobian or preconditioner?");
       repcount = maxit+2; 
      }  else    if (idid == -8)    {
       warning("matrix of partial derivatives is singular with direct method-some equations redundant");
       repcount = maxit+2; 
      }  else    if (idid == -9)    {
       warning("repeated convergence test failures and error test failures ?");
       repcount = maxit+2; 
      }  else    if (idid == -10)   {
       warning("repeated convergence test failures on a step, because ires was -1");
       repcount = maxit+2; 
      }  else    if (idid == -11)   {
       warning("unrecoverable error from inside noninear solver, ires=-2 ");
       repcount = maxit+2; 
      }  else    if (idid == -12)   {
       warning("failed to compute initial y and yprime vectors");
       repcount = maxit+2; 
      }  else    if (idid == -13)   {
       warning("unrecoverable error inside the PSOL routine");
       repcount = maxit+2; 
      }  else    if (idid == -14)   {
       warning("Krylov linear system solver failed to converge");
       repcount = maxit+2; 
      }  else    if (idid == -33)   {
       warning("fatal error");
       repcount = maxit+2; 
      }

	} while (tin < tout && repcount < maxit);

 	  REAL(YOUT)[(it+1)*(ntot+1)] = tin;
	  for (j = 0; j < n_eq; j++)
	    REAL(YOUT)[(it+1)*(ntot + 1) + j + 1] = xytmp[j];

	  if (nout>0) {
	    if (isDll == 1) res_func (&tin, xytmp, xdytmp, &cj, delta, &ires, out, ipar) ;
 	    else C_out(&nout,&tin,xytmp,xdytmp,out);
      for (j = 0; j < nout; j++)
	       REAL(YOUT)[(it+1)*(ntot + 1) + j + n_eq + 1] = out[j]; 
               }
               
/*                    ####  an error occurred   ####                          */                     
    if (repcount > maxit || tin < tout || idid <= 0) {
      idid = 0;
      returnearly(1, it, ntot);
    	break;
    }
  }    /* end main time loop */

/*                   ####   returning output   ####                           */    
  terminate(idid, iwork, 23, 0, rwork, 3, 1);
  REAL(RWORK)[0] = rwork[6];
    
  //unprotect_all();
  restore_N_Protected(old_N_Protect);  
  unlock_solver();
  
  if (idid > 0)
    return(YOUT);
  else
    return(YOUT2);
}

