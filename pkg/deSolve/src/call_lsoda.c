#include <time.h>
#include <string.h>
#include "deSolve.h"

void F77_NAME(dlsoda)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *),
		     int *, double *, int *);

void F77_NAME(dlsode)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *),
		     int *, double *, int *);

void F77_NAME(dlsodes)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *, int *, int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, int *, double *, double *, int *),   /* jacvec */
		     int *, double *, int *);

void F77_NAME(dlsodar)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *), int *, 
		     void (*)(int *, double *, double *, int *, double *),  /* rootfunc */
         int *, int *, double *, int *);

static void lsoda_derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(odesolve_deriv_func,Time,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));           incr_N_Protect();

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);
}

static void lsoda_root (int *neq, double *t, double *y, int *ng, double *gout)
{
  int i;
  SEXP R_fcall, ans;
                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(odesolve_root_func,Time,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));          incr_N_Protect();

  for (i = 0; i < *ng; i++)   gout[i] = REAL(ans)[i];

  my_unprotect(2);
}

static void lsoda_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                             REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(odesolve_jac_func,Time,Y));    incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));          incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(2);
}

static void lsoda_jacvec (int *neq, double *t, double *y, int *j,
		    int *ian, int *jan, double *pdj, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans, J;
  PROTECT(J = NEW_INTEGER(1));                  incr_N_Protect();
                             INTEGER(J)[0] = *j;
                             REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang4(odesolve_jac_vec,Time,Y,J));    incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));          incr_N_Protect();

  for (i = 0; i < *neq ; i++)  pdj[i] = REAL(ans)[i];

  my_unprotect(3);
}

typedef void deriv_func(int *, double *, double *,double *, double *, int *);
typedef void root_func(int *, double *, double *,int *, double *);
typedef void jac_func(int *, double *, double *, int *,
		                  int *, double *, int *, double *, int *);
typedef void jac_vec(int *, double *, double *, int *,
		                  int *, int *, double *, double *, int *);
typedef void init_func(void (*)(int *, double *));

SEXP call_lsoda(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP iTask, SEXP rWork, SEXP iWork, SEXP jT, 
    SEXP nOut, SEXP lRw, SEXP lIw, SEXP Solver, SEXP rootfunc, 
    SEXP nRoot, SEXP Rpar, SEXP Ipar)

{
  SEXP yout, yout2, ISTATE, RWORK, IROOT;    

  int  i, j, k, nt, repcount, latol, lrtol, lrw, liw, isOut, maxit, solver;
  double *xytmp, *rwork, tin, tout, *Atol, *Rtol, *out, *dy, ss;
  int neq, itol, itask, istate, iopt, *iwork, jt, mflag, nout, ntot, is;
  int nroot, *jroot, isroot, *ipar, lrpar, lipar, isDll;
  
  deriv_func *derivs;
  jac_func *jac;
  jac_vec *jacvec;
  root_func *root;
  init_func *initializer;
    

  /* #### initialisation #### */    

  init_N_Protect();

  jt = INTEGER(jT)[0];  /* method flag */
  neq = LENGTH(y);
  nt = LENGTH(times);
  
  maxit = 10;  
  mflag = INTEGER(verbose)[0];
 
  nout   = INTEGER(nOut)[0]; 
  nroot  = INTEGER(nRoot)[0]; 
  solver = INTEGER(Solver)[0]; 
  
  /*  1= lsoda, 2=lsode: 3=lsodeS, 4=lsodar*/
  
  /*dummy parameters RPAR and IPAR in lsodx: used to pass output variables 
    note: fortran function lsodx and dependencies have been altered to also
    pass RPAR and IPAR*/
  
  if (inherits(func, "NativeSymbol"))  /* function is a dll */
  {
   isDll = 1;
   if (nout > 0) isOut = 1; 
   ntot  = neq + nout; /* length of yout */
   lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
   lipar = 3 + LENGTH(Ipar); /* length of ipar */

  } else                              /* function is not a dll */
  {
   isDll = 0;
   isOut = 0;
   ntot = neq;
   lipar = 1;
   lrpar = 1; 
  }
 
   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;
    ipar[1] = lrpar;
    ipar[2] = lipar;
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];
   
    for (j = 0; j < nout; j++) out[j] = 0.;  
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }
   
  /* copies of all variables that will be changed in the FORTRAN subroutine */

  xytmp = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < neq; j++) xytmp[j] = REAL(y)[j];
 
  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));

  liw = INTEGER (lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));
     for (j=0; j<LENGTH(iWork); j++) iwork[j] = INTEGER(iWork)[j];

  lrw = INTEGER(lRw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
     for (j=0; j<length(rWork); j++) rwork[j] = REAL(rWork)[j];

  /* initialise global variables... */

  PROTECT(Time = NEW_NUMERIC(1));                  incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,(neq)));         incr_N_Protect();
  PROTECT(yout = allocMatrix(REALSXP,ntot+1,nt));  incr_N_Protect();
  PROTECT(de_gparms = parms);                      incr_N_Protect();

  /* If there is an initializer, use it here */
  if (!isNull(initfunc))
	  {
	  initializer = (init_func *) R_ExternalPtrAddr(initfunc);
	  initializer(Initdeparms);
	  }

  /* pointers to functions derivs and jac, passed to the FORTRAN subroutine */

  if (inherits(func, "NativeSymbol")) 
    {
      derivs = (deriv_func *) R_ExternalPtrAddr(func);
      if (isOut) {dy = (double *) R_alloc(neq, sizeof(double));
                  for (j = 0; j < neq; j++) dy[j] = 0.; }
	  
    } else {
      derivs = (deriv_func *) lsoda_derivs;
/* KS: removed the PROTECT part...   
      PROTECT(odesolve_deriv_func = func);        incr_N_Protect();
      PROTECT(odesolve_envir = rho);              incr_N_Protect();  */
      
      odesolve_deriv_func = func;
      odesolve_envir = rho;

    }

  if (!isNull(jacfunc) && solver !=3)
    {
      if (inherits(jacfunc,"NativeSymbol"))
	    {
	     jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
	    } else  {
	     odesolve_jac_func = jacfunc;
	     jac = lsoda_jac;
	    }
    }

  if (!isNull(jacfunc) && solver ==3)
    {
      if (inherits(jacfunc,"NativeSymbol"))
	    {
	     jacvec = (jac_vec *) R_ExternalPtrAddr(jacfunc);
	    } else  {
	     odesolve_jac_vec = jacfunc;
	     jacvec = lsoda_jacvec;
	    }
    }

  if (solver == 4 && nroot > 0)  
  { jroot = (int *) R_alloc(nroot, sizeof(int));
     for (j=0; j<nroot; j++) jroot[j] = 0;
  
    if (inherits(rootfunc, "NativeSymbol")) 
    {
      root = (root_func *) R_ExternalPtrAddr(rootfunc);
    } else {
      root = (root_func *) lsoda_root;
/* and here      PROTECT(odesolve_root_func = rootfunc);   incr_N_Protect(); */
      odesolve_root_func = rootfunc; 
    }
  }

  if (latol == 1 && lrtol == 1 ) itol = 1;
  if (latol  > 1 && lrtol == 1 ) itol = 2;
  if (latol == 1 && lrtol  > 1 ) itol = 3;
  if (latol  > 1 && lrtol  > 1 ) itol = 4;

  itask = INTEGER(iTask)[0];   
  istate = 1;

  iopt = 0;
  ss = 0.;
  is = 0 ;
  for (i = 5; i < 8 ; i++) ss = ss+rwork[i];
  for (i = 5; i < 10; i++) is = is+iwork[i];
  if (ss >0 || is > 0) iopt = 1;


  REAL(yout)[0] = REAL(times)[0];
  for (j = 0; j < neq; j++)
    {
      REAL(yout)[j+1] = REAL(y)[j];
    }

    if (isOut == 1) {  /* output in DLL */
        derivs (&neq, &tin, xytmp, dy, out, ipar) ;
	      for (j = 0; j < nout; j++)
	       REAL(yout)[j + neq + 1] = out[j]; 
               }
/* #### main time loop #### */    
  for (i = 0; i < nt-1; i++)
    {
      tin = REAL(times)[i];
      tout = REAL(times)[i+1];
      repcount = 0;
      for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];
      for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];
      do
	{  /* error control */
	  if (istate == -2)
	    {
	      for (j = 0; j < lrtol; j++) Rtol[j] *= 10.0;
	      for (j = 0; j < latol; j++) Atol[j] *= 10.0;
	      warning("Excessive precision requested.  `rtol' and `atol' have been scaled upwards by the factor %g\n",10.0);
	      istate = 3;
	    }

    if (solver == 1)
    {	    
	  F77_CALL(dlsoda) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, NUMERIC_POINTER(rtol), NUMERIC_POINTER(atol), 
         &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar); 
    } else if (solver == 2) {
    F77_CALL(dlsode) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, NUMERIC_POINTER(rtol), NUMERIC_POINTER(atol), 
         &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar); 
    } else if (solver == 3) {
    F77_CALL(dlsodes) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, NUMERIC_POINTER(rtol), NUMERIC_POINTER(atol), 
         &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jacvec, &jt, out, ipar); 
    } else if (solver == 4) {
    F77_CALL(dlsodar) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, NUMERIC_POINTER(rtol), NUMERIC_POINTER(atol), 
         &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, root, &nroot, jroot, 
         out, ipar); 
    }
    
	  if (istate == -1) 
     {
      warning("an excessive amount of work (> maxsteps ) was done, but integration was successful - increase maxsteps");
     }
	  if (istate == 3 && solver == 4)
	    { istate = -20;  repcount = 50;  
      }

	  repcount ++;
	} while (tin < tout && istate >= 0 && repcount < maxit); 
      if (istate == -3)
	{
	  unprotect_all();
	  error("Illegal input to lsoda\n");
	}
      else
	{
	  REAL(yout)[(i+1)*(ntot+1)] = tin;
	  for (j = 0; j < neq; j++)
	    REAL(yout)[(i+1)*(ntot + 1) + j + 1] = xytmp[j];
	  if (isOut == 1) 
    {
        derivs (&neq, &tin, xytmp, dy, out, ipar) ;
	      for (j = 0; j < nout; j++)
	       REAL(yout)[(i+1)*(ntot + 1) + j + neq + 1] = out[j]; 
    }
	}
	  
 /* KS: added || tin < tout */
  if (istate < 0 || tin < tout) {
	  if (istate != -20) warning("Returning early.  Results are accurate, as far as they go\n");
	 /* need to redimension yout here, and add the attribute "istate" for */
	 /* the most recent value of `istate' from lsod */
	  PROTECT(yout2 = allocMatrix(REALSXP,ntot+1,(i+2)));incr_N_Protect();
	  for (k = 0; k < i+2; k++)
	   for (j = 0; j < ntot+1; j++)
	     REAL(yout2)[k*(ntot+1) + j] = REAL(yout)[k*(ntot+1) + j];
	  break;
      }
    }     /* end main time loop */

/* #### returning output if (istate == -20) istate = 3;#### */ 
  if (istate == -20 && nroot > 0) 
   { isroot = 1   ;
     PROTECT(IROOT = allocVector(INTSXP, nroot));incr_N_Protect();
     for (k = 0;k<nroot;k++) INTEGER(IROOT)[k] = jroot[k];
   } else isroot = 0;

  PROTECT(ISTATE = allocVector(INTSXP, 22));incr_N_Protect();
  for (k = 0;k<22;k++) INTEGER(ISTATE)[k+1] = iwork[k];

        
  PROTECT(RWORK = allocVector(REALSXP, 5));incr_N_Protect();
  for (k = 0;k<5;k++) REAL(RWORK)[k] = rwork[k+10];

  INTEGER(ISTATE)[0] = istate;  
  if (istate == -20) INTEGER(ISTATE)[0] = 3; 	  
  if (istate > 0)
    {
      setAttrib(yout, install("istate"), ISTATE);
      setAttrib(yout, install("rstate"), RWORK);    
      if (isroot==1) setAttrib(yout, install("iroot"), IROOT);    
      }
  else
    {
      setAttrib(yout2, install("istate"), ISTATE);
      setAttrib(yout2, install("rstate"), RWORK);    
      if (isroot==1) setAttrib(yout2, install("iroot"), IROOT);    
      }

  unprotect_all();
  if (istate > 0)
    return(yout);
  else
    return(yout2);
}

