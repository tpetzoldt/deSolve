#include <R.h>
#include <Rinternals.h>

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
void checkeventtimes(int* nevents, int* nx, 
                    double* events, double* x, 
                    int* relative, double* xout)  {

  int i = 0, j = 0;
  double divisor;
  for (i = 0; i < *nx; i++) {
    xout[i] = 1e99;
  }
  /* case 1: first value */
  if (*relative) divisor = (x[1] - x[0]); else divisor = 1.0;
  xout[0] = fabs(x[0] - events[0]) / divisor;
  /* case 2: last value  */
  if (*relative) divisor = (x[*nx - 1] - x[*nx - 2]); else divisor = 1.0;
  xout[*nx - 1] = fabs(x[*nx - 1] - events[*nevents - 1]) / divisor;
  
  /* case 3: all others */
  for (i = 1; i < *nx - 1; i++) { //*nx; i++) {
    findindex(nevents, &x[i], events, &j);
    //Rprintf("%d -- %e \n", j, events[j]);
    if (j >= 0) {
      if (*relative)  divisor = x[i] - x[i - 1]; else divisor = 1.0;
      xout[i] = (x[i] - events[j]) / divisor;
    }
    if (j < *nevents) {
      if (*relative) divisor = (x[i+1] - x[i]); else divisor = 1.0;
      xout[i] = fmin(fabs(xout[i]), fabs(x[i] - events[j + 1]) / divisor);
    }
  }
}
