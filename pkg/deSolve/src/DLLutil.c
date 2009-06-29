/* Functions to test compiled code implementation of ODE and DAE */

#include <time.h>
#include <string.h>
#include "deSolve.h"

SEXP call_DLL(SEXP y, SEXP dY, SEXP time, SEXP func, SEXP initfunc, SEXP parms,
		          SEXP nOut, SEXP Rpar, SEXP Ipar, SEXP Type,SEXP Tvec, SEXP Fvec,
              SEXP Ivec,SEXP initforc)
{
  SEXP   yout;

  double *ytmp, *dy, tin, *out, *delta, cj;
  int    ny, *ipar, ntot, i, j, nout, type, lrpar, lipar, ires;
  
  deriv_func *derivs;
  res_func *res;
  init_func  *initializer, *initforcings;

  init_N_Protect();

  ny   = LENGTH(y);
  nout = INTEGER(nOut)[0];
  ntot = ny+nout;
  type = INTEGER(Type)[0];

  lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
  lipar = 3 + LENGTH(Ipar);    /* length of ipar */
  out   = (double *) R_alloc(lrpar, sizeof(double));

  ipar  = (int *)    R_alloc(lipar, sizeof(int));
  ipar[0] = nout;
  ipar[1] = lrpar;
  ipar[2] = lipar;
  for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];

  for (j = 0; j < nout; j++) out[j] = 0.;  
  for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
  
  
  PROTECT(de_gparms = parms)                   ; incr_N_Protect();
  PROTECT(yout = allocVector(REALSXP,ntot))    ; incr_N_Protect();

 /* The initialisation routine */
  if (!isNull(initfunc))     	{
	     initializer = (init_func *) R_ExternalPtrAddr(initfunc);
	     initializer(Initdeparms);
  }

  if (!isNull(initforc))     	{
    nforc =LENGTH(Ivec)-2; /* nforc, fvec, ivec =globals */

    i = LENGTH(Fvec);
    fvec = (double *) R_alloc((int) i, sizeof(double));
    for (j = 0; j < i; j++) fvec[j] = REAL(Fvec)[j];

    tvec = (double *) R_alloc((int) i, sizeof(double));
    for (j = 0; j < i; j++) tvec[j] = REAL(Tvec)[j];

    i = LENGTH (Ivec)-1; /* last element: the interpolation method...*/
    ivec = (int *) R_alloc(i, sizeof(int));
    for (j = 0; j < i; j++) ivec[j] = INTEGER(Ivec)[j];

    fmethod =INTEGER(Ivec)[i];

	  initforcings = (init_func *) R_ExternalPtrAddr(initforc);
    initforcings(Initdeforc);
  }

    
  tin = REAL(time)[0];

  ytmp = (double *) R_alloc(ny, sizeof(double));
    for (j = 0; j < ny; j++) ytmp[j] = REAL(y)[j];

  dy   = (double *) R_alloc(ny, sizeof(double));
    for (j = 0; j < ny; j++) dy[j] = REAL(dY)[j]; 

  if(!isNull(initforc))  updatedeforc(&tin);

  if (type == 1)   {
    derivs = (deriv_func *) R_ExternalPtrAddr(func);

    derivs (&ny, &tin, ytmp, dy, out, ipar) ;
    for (j = 0; j < ny; j++)  REAL(yout)[j] = dy[j];

  } else {

    res = (res_func *) R_ExternalPtrAddr(func);
    delta = (double *) R_alloc(ny, sizeof(double));
    for (j = 0; j < ny; j++) delta[j] = 0.;

    res    (&tin, ytmp, dy, &cj, delta, &ires, out, ipar) ;
    for (j = 0; j < ny; j++)  REAL(yout)[j] = delta[j];

  }
                  
  if (nout > 0)   {

	   for (j = 0; j < nout; j++)
	       REAL(yout)[j + ny] = out[j]; 
  }

  unprotect_all();
  return(yout);
}