#include <R.h>
#include <Rinternals.h>

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif



/* bisection method: find index of neighboring values for xi in x */
void findindex(int* n, double* xi, double* x, int* i) {
  int lo = 0, up = *n, m;
  do {
    m = floor((lo + up)/2);
    if (*xi >= x[m])	lo = m; else up = m;
  } while ((up - lo) > 1);
  *i = lo;
}  

/* find distance of any x to the nearest of 'events' */
double fdivisor(double evt, double x1, double x2, int relmode) {
    double divisor = 1.0;
    if (relmode == 0)
	divisor = 1.0;
    else if (relmode == 1)
	divisor = evt;
    else if (relmode == 2)
	divisor = x1 - x2;
    else
	Rprintf("not definded mode\n");

    return(divisor);
}


void checkeventtimes(int* nevents, int* nx, 
                    double* events, double* x, 
                    int* relative, double* xout)  {

  int i = 0, j = 0;
  int relmode = *relative;
  double divisor, evt;

  /* assume DOUBLE_XMAX as maximal theroretical value*/
  for (i = 0; i < *nx; i++) {
    xout[i] = DOUBLE_XMAX;
  }
  /* case 1: first value */
  //if (relmode) divisor = (x[1] - x[0]); else divisor = 1.0;
  evt = events[0];
  xout[0] = fabs(x[0] - evt) / fdivisor(evt, x[1], x[0], relmode);
  /* case 2: last value  */
  //if (*relative) divisor = (x[*nx - 1] - x[*nx - 2]); else divisor = 1.0;
  evt = events[*nevents - 1];
  xout[*nx - 1] = fabs(x[*nx - 1] - evt) / fdivisor(evt, x[*nx - 1], x[*nx - 2], relmode);
  
  /* case 3: all others */
  for (i = 1; i < *nx - 1; i++) { //*nx; i++) {
    findindex(nevents, &x[i], events, &j);
    //Rprintf("%d -- %e \n", j, events[j]);
    if (j >= 0) {
        //if (relmode)  divisor = x[i] - x[i - 1]; else divisor = 1.0;
        evt = events[j];
	xout[i] = (x[i] - evt) / fdivisor(evt, x[i], x[i - 1], relmode);
    }
    if (j < *nevents) {
	//if (relmode) divisor = (x[i+1] - x[i]); else divisor = 1.0;
        evt = events[j + 1];
	xout[i] = fmin(fabs(xout[i]), fabs(x[i] - evt)) / fdivisor(evt, x[i + 1], x[i], relmode);
    }
  }
}


// assumes sorted x and events
void nearest_event(int* nevents, int* nx, 
                   double* events, double* x, double* xout)  {

  int i = 0, j = 0;

  /* initialize output with NA */
  for (i = 0; i < *nx; i++) {
    xout[i] = R_NaReal;
  }
  
  for (i = 0; i < *nx; i++) {
    findindex(nevents, &x[i], events, &j);
    if ((x[i] - events[j]) < (events[min(j + 1, *nevents - 1)] - x[i]))
  	  xout[i] = events[j];
 	  else
 	    xout[i] = events[min(j + 1, *nevents - 1)];
  }
}