/* deals with forcing functions that are passed via arguments in the call
to the integration routines */

#include "deSolve.h"

int    finit = 0;


/*         -----     INITIALISATION     -----
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
     intpol[i] = (fvec[ii+1]-fvec[ii])/(nexttime[i]-curtime[i]);
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
       if (zerograd == 0) {
         nexttime[i] = tvec[ii+1];
         intpol[i] = (fvec[ii+1]-fvec[ii])/(nexttime[i]-curtime[i]); }
       else {
         nexttime[i] = 1e20 ; /* as large as possible */
         intpol[i] = 0;
       }
     }

     forcings[i]=curval[i]+intpol[i]*(*time-curtime[i]);

   }

}

