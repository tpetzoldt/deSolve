/* suitable names for parameters and state variables */

#include <R.h>
static double parms[7];

#define eps parms[0]
#define m   parms[1]
#define L   parms[2]
#define L0  parms[3]
#define r   parms[4]
#define w   parms[5]
#define g   parms[6]

/*----------------------------------------------------------------------
 initialising the parameter common block
----------------------------------------------------------------------
*/
void carparc(void (* daeparms)(int *, double *)) {
  int N = 7;
  daeparms(&N, parms);
}
/* Compartments */

#define xl y[0]
#define yl y[1]
#define xr y[2]
#define yr y[3]
#define lam1 y[8]
#define lam2 y[9]

/*----------------------------------------------------------------------
 the residual function
----------------------------------------------------------------------
*/
void carresc (double *t, double *y, double *yprime, double *CJ,
             double *delta, int *ier, double *FOUT, int* ip) {

    double k, yb, xb, Lr, Ll;
    int i;
      
      if (ip[0] < 10) error("nout should be at least 10");

      k = m*eps*eps/2.;

      yb  = r*sin(w* *t) ;
      xb  = sqrt(L*L-yb*yb);

      for (i = 0; i <4; i++)
        FOUT[i] = y[i+4];

      Ll = sqrt(xl*xl+yl*yl) ;
      Lr = sqrt((xr-xb)*(xr-xb)+(yr-yb)*(yr-yb));

      FOUT[4]  =(L0-Ll)*xl/Ll +lam1*xb+2*lam2*(xl-xr)    ;
      FOUT[5]  =(L0-Ll)*yl/Ll +lam1*yb+2*lam2*(yl-yr)-k*g;
      FOUT[6]  =(L0-Lr)*(xr-xb)/Lr -2*lam2*(xl-xr)       ;
      FOUT[7]  =(L0-Lr)*(yr-yb)/Lr -2*lam2*(yl-yr)-k*g   ;

      FOUT[8]  = xb*xl+yb*yl;
      FOUT[9] = (xl-xr)*(xl-xr)+(yl-yr)*(yl-yr)-L*L;
      for (i = 0 ; i <4; i++)
         delta[i] =yprime[i]-FOUT[i];

      for (i = 4 ; i <8; i++)
         delta[i] = k*yprime[i]- FOUT[i];

      for (i = 8 ; i <10; i++)
         delta[i] = -FOUT[i];
}
