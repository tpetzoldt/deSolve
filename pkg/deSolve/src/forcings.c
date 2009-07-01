/* deals with forcing functions that are passed via arguments in the call
to the integration routines */

#include "deSolve.h"

int    finit = 0;

/*         -----     Check for presence of forcing functions     -----        */


int initForcings(SEXP flist) {

    SEXP Tvec, Fvec, Ivec, initforc;
    int i, j, isForcing = 0;
    init_func  *initforcings;


    initforc = getListElement(flist,"ModelForc");
    if (!isNull(initforc))
    	{
    	 Tvec = getListElement(flist,"tmat");
       Fvec = getListElement(flist,"fmat");
       Ivec = getListElement(flist,"imat");
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

       isForcing = 1;
       }
       return(isForcing);
}

/*         -----     INITIALISATION  called from compiled code   -----
   1. Check the length of forcing functions in solver call and code in DLL
   2. Initialise the forcing function vectors
   3. set pointer to DLL fortran common block or C globals /
*/

void Initdeforc(int *N, double *forc)
{
  int i, ii;

  if ((*N) != nforc) {
    PROBLEM "Confusion over the length of forc"
    ERROR;
  }

/* for each forcing function: index to current position of data,
 current value, interpolation factor, current forcing time, next forcing time,..
*/
   finit = 1;
   findex   = (int    *) R_alloc(nforc, sizeof(int));
   curval   = (double *) R_alloc(nforc, sizeof(double));
   intpol   = (double *) R_alloc(nforc, sizeof(double));
   curtime  = (double *) R_alloc(nforc, sizeof(double));
   nexttime = (double *) R_alloc(nforc, sizeof(double));
   maxindex = (int    *) R_alloc(nforc, sizeof(int));


/* Input is in three vectors:
   tvec, fvec: time and value;
   ivec : index to each forcing in tvec and fvec
*/
   for (i = 0; i<nforc; i++) {
     ii = ivec[i]-1;
     findex[i] = ii;
     maxindex[i] = ivec[i+1]-2;
     curval[i] = fvec[ii];
     curtime[i] = tvec[ii];
     nexttime[i] = tvec[ii+1];
     if (fmethod == 1) {
       intpol[i] = (fvec[ii+1]-fvec[ii])/(nexttime[i]-curtime[i]);
     } else  intpol[i] = 0;
     forc[i] = curval[i];
   }
   forcings = forc;      /* set pointer to c globals or fortran common block */

}

/*         -----     UPDATING forcing functions at each time     -----
*/

void updatedeforc(double *time)
{
  int i, ii, change, zerograd;
  double ntime;

/* check if initialised? */
   if (finit == 0)
     error ("error in forcing function: not initialised");

   for (i=0; i<nforc; i++) {
     ii = findex[i];
     change=0;
     zerograd=0;
     if (*time > nexttime[i])
     {
       ntime = nexttime[i];
       while (*time > ntime){
         if (ii+1 > maxindex[i]) {   /* this probably redundant...*/
           zerograd=1;
           break;
         }
         ii = ii+1;
         ntime = tvec[ii];
       }
       change=1;
     }
     if (*time < curtime[i])
     {
       ntime = curtime[i];
       while (*time < ntime){
         ii = ii-1;
         ntime = tvec[ii];
       }
       change=1;
     }
     if (change == 1) {
       findex[i] = ii;
       curval[i] = fvec[ii];
       curtime[i] = tvec[ii];
       if ((zerograd == 0) & (fmethod == 1)) {  /* fmethod 1=linear */
         nexttime[i] = tvec[ii+1];
         intpol[i] = (fvec[ii+1]-fvec[ii])/(nexttime[i]-curtime[i]); }
       else if (fmethod == 2) {
         nexttime[i] = tvec[ii+1];
         intpol[i] = 0;
       }
       else {
         nexttime[i] = DBL_MAX ; /* as large as possible */
         intpol[i] = 0;
       }
       
     }

     forcings[i]=curval[i]+intpol[i]*(*time-curtime[i]);

   }

}

