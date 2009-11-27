/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* General RK Solver for methods with fixed resp. adaptive step size        */
/* -- core function (main loop) --                                          */
/*==========================================================================*/


void rk_auto(
       /* integers */
       int fsal, int neq, int stage,
       int isDll, int isForcing, int verbose,
       int nknots, int interpolate, int maxsteps, int nt,
       /* int pointers */
       int* _iknots, int* _it, int* _it_ext, int* _it_tot, 
       int* istate,  int* ipar,
       /* double */
       double t, double tmax, double hmin, double hmax, 
       double alpha, double beta,
       /* double pointers */
       double* _dt, double* _errold,
       /* arrays */
       double* tt, double* y0, double* y1, double* y2, double* dy1, double* dy2,
       double* f, double* y, double* Fj, double* tmp,
       double* FF, double* rr, double* A, double* out, 
       double* bb1, double* bb2, double* cc, double* dd, 
       double* atol, double* rtol, double* yknots, double* yout,
       /* SEXPs */
       SEXP Func, SEXP Parms, SEXP Rho
 );


void rk_fixed(
       /* integers */
       int fsal, int neq, int stage,
       int isDll, int isForcing, int verbose,
       int nknots, int interpolate, int maxsteps, int nt,
       /* int pointers */
       int* _iknots, int* _it, int* _it_ext, int* _it_tot, 
       int* istate,  int* ipar,
       /* double */
        double t, double tmax, double hini,
       /* double pointers */
       double* _dt,
       /* arrays */
       double* tt, double* y0, double* y1,double* dy1, 
       double* f, double* y, double* Fj, double* tmp,
       double* FF, double* rr, double* A, double* out, 
       double* bb1, double* cc, 
       double* yknots, double* yout,
       /* SEXPs */
       SEXP Func, SEXP Parms, SEXP Rho
  );

