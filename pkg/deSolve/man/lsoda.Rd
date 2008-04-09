\name{lsoda}
\alias{lsoda}
\title{General solver for ordinary differential equations (ODE), switching automatically between stiff and non-stiff methods}
\description{Solving initial value problems for stiff or
    non-stiff systems of first-order ordinary differential equations
    (ODEs). 
    
    The \R function \code{lsoda} provides an interface to the
    Fortran ODE solver of the same name, written by Linda R. Petzold and Alan
    C. Hindmarsh.  
    
    The system of ODE's is written as an \R function (which
    may, of course, use \code{\link{.C}}, \code{\link{.Fortran}},
    \code{\link{.Call}}, etc., to call foreign code) or be defined in
    compiled code that has been dynamically loaded.  A vector of
    parameters is passed to the ODEs, so the solver may be used as part of
    a modeling package for ODEs, or for parameter estimation using any
    appropriate modeling tool for non-linear models in \R such as
    \code{\link{optim}}, \code{\link[nls]{nls}}, \code{\link{nlm}} or
    \code{\link[nlme]{nlme}}
    
  lsoda differs from the other integrators (except \code{lsodar}) in that it switches automatically between stiff and nonstiff methods.
  This means that the user does not have to determine whether the
  problem is stiff or not, and the solver will automatically choose the
  appropriate method.  It always starts with the nonstiff method.}
\usage{lsoda(y, times, func, parms, rtol=1e-6, atol=1e-6, tcrit=NULL, 
  jacfunc=NULL, jactype="fullint", verbose=FALSE,   
  dllname=NULL, initfunc=dllname, initpar=parms, 
  rpar=NULL, ipar=NULL, ynames=TRUE, nout=0, outnames=NULL, 
  hmin=0, hmax=NULL, hini=0, maxordn=12, maxords = 5, 
  bandup=NULL, banddown=NULL, maxsteps=5000, ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{times at which explicit estimates for \code{y} are desired.  The first value in \code{times} must be the initial time.}
  \item{func }{either a user-supplied function that computes the values of the
      derivatives in the ODE system (the \emph{model definition}) at time
      t, or a character string
      giving the name of a compiled function in a dynamically loaded
      shared library.
  
      If \code{func}  is a user-supplied function, it must be called as:
      \code{yprime = func(t, y, parms)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ODE system.  If the initial values \code{y} has a names
      attribute, the names will be available inside \code{func}.  \code{parms} is
      a vector of parameters (which
      may have a names attribute, desirable in a large system).
  
      The return value of \code{func} should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose next elements (possibly with a
      \code{\link{names}} attribute) are global values that are required at
      each point in \code{times}.
  
      If \code{func} is a string, then \code{dllname} must give the name
      of the shared library (without extension) which must be loaded
      before \code{lsoda()} is called.  see Details for more information.}
  \item{parms }{any parameters used in \code{func} that should be
      modifiable without rewriting the function.}
  \item{rtol }{relative error tolerance, either a scalar or an array as
      long as \code{y}.  See details. }
  \item{atol }{absolute error tolerance, either a scalar or an array as
      long as \code{y}.  See details.}
  \item{tcrit }{if not \code{NULL}, then \code{lsoda} cannot integrate past \code{tcrit}. The Fortran routine \code{lsoda} overshoots its targets
      (times points in the vector \code{times}), and interpolates values
      for the desired time points.  If there is a time beyond which
      integration should not proceed (perhaps because of a singularity),
      that should be provided in \code{tcrit}.  Note that it does not
      make sense (though it is not an error) to include times in
      \code{times} past \code{tcrit}, since the solver will stop and
      return at the last point in \code{times} that is earlier than
      \code{tcrit}.}
  \item{jacfunc }{if not \code{NULL}, an \R function that computes
      the jacobian of the system of differential equations
      dydot(i)/dy(j), or a string giving the name of a function or
      subroutine in \file{dllname} that computes the jacobian (see Details
      below for more about this option).  In some circumstances, supplying
      \code{jacfunc} can speed up
      the computations, if the system is stiff.  The \R calling sequence for
      \code{jacfunc} is identical to that of \code{func}.  
          \cr If the jacobian is a full matrix, \code{jacfunc} should return a matrix dydot/dy, where the ith
      row contains the derivative of \eqn{dy_i/dt} with respect to \eqn{y_j},
      or a vector containing the matrix elements by columns (the way \R and Fortran store matrices).
      \cr If the jacobian is banded, \code{jacfunc} should return a matrix containing only the 
      nonzero bands of the jacobian, rotated row-wise. See first example of lsode.}
  \item{jactype }{the structure of the jacobian, one of "fullint", "fullusr", "bandusr" or "bandint" - either full or banded and estimated internally or by user}
  \item{verbose }{a logical value that, when TRUE, triggers more
      verbose output from the ODE solver.  Will output the settings of vectors *istate* and *rstate* - see details}
  \item{dllname }{a string giving the name of the shared library (without
      extension) that contains all the compiled function or subroutine
      definitions refered to in \code{func} and \code{jacfunc}.}
  \item{initfunc }{if not NULL, the name of the initialisation function (which initialises values of parameters), as provided in \file{dllname}. See details.}
  \item{initpar }{only when \file{dllname} is specified and an initialisation function \code{initfunc} is in the dll: the parameters passed to the initialiser, to initialise the common blocks (fortran) or global variables (C, C++)}
  \item{rpar }{only when \file{dllname} is specified: a vector with double precision values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{ipar }{only when \file{dllname} is specified: a vector with integer values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{ynames }{if FALSE: names of state variables are not passed to function \code{func} ; this may speed up the simulation especially for large models}
  \item{nout }{only used if \code{dllname} is specified and the model is defined in compiled code: the number of output variables calculated in the compiled function \code{func}, present in the shared library. Note:
      it is not automatically checked whether this is indeed the number of output variables calculed in the dll - you have to perform this check in the code - see details}
  \item{outnames }{only used if \file{dllname} is specified and \code{nout} > 0: the names of output variables calculated in the compiled function \code{func}, present in the shared library}
  \item{hmin }{an optional minimum value of the integration
      stepsize. In special situations this parameter may speed up computations with
      the cost of precision. Don't use \code{hmin} if you don't know why!}
  \item{hmax }{an optional maximum value of the integration stepsize. If not specified, \code{hmax} is set to the largest difference in \code{times}, to avoid that the simulation possibly ignores short-term events. If 0, no maximal size is specified}
  \item{hini }{initial step size to be attempted; if 0, the initial step size is determined by the solver}
  \item{maxordn }{the maximum order to be allowed in case the method is non-stiff. Should be <=12. Reduce maxord to save storage space}
  \item{maxords }{the maximum order to be allowed in case the method is stiff. Should be <=5. Reduce maxord to save storage space}
  \item{bandup }{number of non-zero bands above the diagonal, in case the Jacobian is banded}
  \item{banddown }{number of non-zero bands below the diagonal, in case the Jacobian is banded}
  \item{maxsteps }{maximal number of steps during one call to the solver}
  \item{... }{additional arguments passed to \code{func} and \code{jacfunc} allowing this to be a generic function}
}
\value{
  A matrix with up to as many rows as elements in \code{times} and as
    many columns as elements in \code{y} plus the number of "global"
    values returned in the next elements of the return from \code{func},
    plus and additional column for the time value.  There will be a row
    for each element in \code{times} unless the Fortran routine `lsoda'
    returns with an unrecoverable error.  If \code{y} has a names
    attribute, it will be used to label the columns of the output value.
    
      The output will have the attributes \code{istate}, and \code{rstate}, two vectors with several useful elements. 
      See details.
      The first element of istate returns the conditions under which the last call to lsoda returned. Normal is \code{istate = 2}.
      If \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen
  
}
\author{R. Woodrow Setzer <setzer.woodrow@epa.gov>}
\examples{### lsexamp -- example from lsoda source code
  
  ## names makes this easier to read, but may slow down execution.
  parms <- c(k1=0.04, k2=1e4, k3=3e7)
  my.atol <- c(1e-6,  1e-10,  1e-6)
  times <- c(0,4 * 10^(-1:10))
  lsexamp <- function(t, y, p)
    {
      yd1 <- -p["k1"] * y[1] + p["k2"] * y[2]*y[3]
      yd3 <- p["k3"] * y[2]^2
      list(c(yd1,-yd1-yd3,yd3),c(massbalance=sum(y)))
    }
  exampjac <- function(t, y, p)
    {
      c(-p["k1"],	 p["k1"],  0,
    
          p["k2"]*y[3],
        - p["k2"]*y[3] - 2*p["k3"]*y[2],
                         2*p["k3"]*y[2],
    
        p["k2"]*y[2],  -p["k2"]*y[2],  0
        )
    }
    
  require(odesolve)
  ## measure speed (here and below)
  system.time( 
  out <- lsoda(c(1,0,0),times,lsexamp, parms, rtol=1e-4, atol= my.atol)
  )  
  out
  
  ## This is what the authors of lsoda got for the example:
  
  ## the output of this program (on a cdc-7600 in single precision)
  ## is as follows..
  ##
  ## at t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02
  ## at t =  4.0000e+00   y =  9.055333e-01  2.240655e-05  9.444430e-02
  ## at t =  4.0000e+01   y =  7.158403e-01  9.186334e-06  2.841505e-01
  ## at t =  4.0000e+02   y =  4.505250e-01  3.222964e-06  5.494717e-01
  ## at t =  4.0000e+03   y =  1.831975e-01  8.941774e-07  8.168016e-01
  ## at t =  4.0000e+04   y =  3.898730e-02  1.621940e-07  9.610125e-01
  ## at t =  4.0000e+05   y =  4.936363e-03  1.984221e-08  9.950636e-01
  ## at t =  4.0000e+06   y =  5.161831e-04  2.065786e-09  9.994838e-01
  ## at t =  4.0000e+07   y =  5.179817e-05  2.072032e-10  9.999482e-01
  ## at t =  4.0000e+08   y =  5.283401e-06  2.113371e-11  9.999947e-01
  ## at t =  4.0000e+09   y =  4.659031e-07  1.863613e-12  9.999995e-01
  ## at t =  4.0000e+10   y =  1.404280e-08  5.617126e-14  1.000000e+00
  
  ## Using the analytic jacobian speeds up execution a little :
  
  system.time( 
  outJ <- lsoda(c(1,0,0),times,lsexamp, parms, rtol=1e-4, atol= my.atol,
                jacfunc = exampjac)
  )
    
  all.equal(out, outJ) # TRUE
  
  ## Example for using hmax
  
  ## Parameters for steady state conditions
  parms  <- c(a=0.0, b=0.0, c=0.1, d=0.1, e=0.1, f=0.1, g=0.0)
  
  ## A simple resource limited Lotka-Volterra-Model
  ## Note passing parameters through using a closure
  lvmodel <- with(as.list(parms), function(t, x, parms) {
      import <- sigimp(t)
      ds <- import - b*x["s"]*x["p"] + g*x["k"]
      dp <- c*x["s"]*x["p"]  - d*x["k"]*x["p"]
      dk <- e*x["p"]*x["k"]  - f*x["k"]
      res<-c(ds, dp, dk)
      list(res)
    })
  
  ## vector of timesteps
  times  <- seq(0, 100, length=101)
  
  ## external signal with rectangle impulse
  signal <- as.data.frame(list(times = times,
                              import = rep(0,length(times))))
  
  signal$import[signal$times >= 10 & signal$times <=11] <- 0.2
  
  sigimp <- approxfun(signal$times, signal$import, rule=2)
  
  
  ## Start values for steady state
  y<-xstart <- c(s=1, p=1, k=1)
  
  ## LSODA (default step size)
  out2 <- as.data.frame(lsoda(xstart, times, lvmodel, parms))
  
  ## LSODA: with fixed maximum time step
  out3 <- as.data.frame(lsoda(xstart, times, lvmodel, parms, hmax=1))
  
  par(mfrow=c(2,2))
  plot (out2$time, out2$s,  type="l",   ylim=c(0,3))
  lines(out3$time, out3$s, col="green", lty="dotted")
  
  plot (out2$time, out2$p, type="l",    ylim=c(0,3))
  lines(out3$time, out3$p, col="green", lty="dotted")
  
  plot (out2$time, out2$k, type="l",    ylim=c(0,3))
  lines(out3$time, out3$k, col="green", lty="dotted")
  
  plot (out2$p, out2$k, type="l", ylim=range(out2$k,out3$k))
  lines(out3$p, out3$k, col="green", lty="dotted")}
\references{Hindmarsh, Alan C. (1983) ODEPACK, A Systematized Collection of ODE
    Solvers; in p.55--64 of  Stepleman, R.W. et al.[ed.] (1983)
    \emph{Scientific Computing}, North-Holland, Amsterdam.
    
    Petzold, Linda R. (1983)  Automatic Selection of Methods for Solving
    Stiff and Nonstiff Systems of Ordinary Differential Equations.
    \emph{Siam J. Sci. Stat. Comput.} \bold{4}, 136--148.
    
    Netlib: \url{http://www.netlib.org}}
\details{All the hard work is done by the Fortran subroutine \code{lsoda},
    whose documentation should be consulted for details (it is included as
    comments in the source file \file{src/opkdmain.f}).  The implementation is based on the
    12 November 2003 version of lsoda, from Netlib. 
  
  \code{lsoda} switches automatically between stiff and nonstiff methods.
  This means that the user does not have to determine whether the
  problem is stiff or not, and the solver will automatically choose the
  appropriate method.  It always starts with the nonstiff method.  
    
  
  The form of the \bold{jacobian} can be specified by \code{jactype} which can take the following values:
  
      \item jactype = "fullint" : a full jacobian, calculated internally by lsoda, the default
      \item jactype = "fullusr" : a full jacobian, specified by user function \code{jacfunc}
      \item jactype = "bandusr" : a banded jacobian, specified by user function \code{jacfunc}; the size of the bands specified by \code{bandup} and \code{banddown}
      \item jactype = "bandint" : a banded jacobian, calculated by lsoda; the size of the bands specified by \code{bandup} and \code{banddown}
  
  if \code{jactype}= "fullusr" or "bandusr" then the user must supply a subroutine \code{jacfunc}.
    
  The following description of \bold{error control} is adapted from that documentation (input arguments
  \code{rtol} and \code{atol}, above):
    
    The input parameters \code{rtol}, and \code{atol} determine the error
    control performed by the solver.  The solver will control the vector
    \bold{e} of estimated local errors in \bold{y}, according to an
    inequality of the form max-norm of ( \bold{e}/\bold{ewt} )
    \eqn{\leq}{<=} 1, where \bold{ewt} is a vector of positive error
    weights.  The values of \code{rtol} and \code{atol} should all be
    non-negative.
    The form of \bold{ewt} is:
    
    \deqn{\mathbf{rtol} \times \mathrm{abs}(\mathbf{y}) + \mathbf{atol}}{\bold{rtol} * abs(\bold{y}) + \bold{atol}}
    
    where multiplication of two vectors is element-by-element.
  
    If the request for precision exceeds the capabilities of the machine,
    the Fortran subroutine lsoda will return an error code; under some
    circumstances, the \R function \code{lsoda} will attempt a reasonable
    reduction of precision in order to get an answer.  It will write a
    warning if it does so.
  
  \bold{Models} may be defined in compiled C or Fortran code, as well as in R. 
  
  \item If \code{func} is a user-supplied \bold{R-function}, it must be called as: yprime = func(t, y, parms). 
  t is the current time point in the integration, y is the current estimate of the variables in the ODE system. 
  The return value of \code{func} should be a list, whose first element is a vector containing the derivatives of y with respect to time, 
  and whose second element contains output variables that are required at each point in time. 
  \cr An example is given below:
  
    \code{model<-function(t,Y,parameters)}\cr
    \code{\{}\cr
    \code{with (as.list(parameters),\{}\cr 
       \code{dy1 = -k1*Y[1] + k2*Y[2]*Y[3]}\cr
       \code{dy3 = k3*Y[2]*Y[2]}\cr
       \code{dy2 = -dy1 - dy3}\cr
      \code{list(c(dy1,dy2,dy3))}\cr         
      \code{\})}\cr
  \code{\}}\cr
  
  The jacobian associated to the above example is:
  
  \code{jac  <- function (t,Y,parameters)}\cr
  \code{\{}\cr
   \code{with (as.list(parameters),\{}\cr 
   \code{PD[1,1] = -k1}\cr
   \code{PD[1,2] = k2*Y[3]}\cr
   \code{PD[1,3] = k2*Y[2]}\cr
   \code{PD[2,1] = k1}\cr
   \code{PD[2,3] = -PD[1,3]}\cr
   \code{PD[3,2] = k3*Y[2]}\cr
   \code{PD[2,2] = -PD[1,2] - PD[3,2]}\cr
    \code{return(PD)}\cr
   \code{\})}\cr 
  \code{\}}\cr 
  
  This model can then be run as follows:
  
  \code{parms <- c(k1 = 0.04, k2 = 1e4, k3=3e7)}\cr       
  \code{Y     <- c(1.0,0.0,0.0)}\cr                       
  \code{times <- c(0,0.4*10^(0:11) )}\cr                  
  \code{PD  <- matrix(nrow=3,ncol=3,data=0)}\cr 
  \code{out <- lsoda(Y,times,model,parms=parms,}\cr   
  \code{jacfunc=jac,rtol=1e-4)}\cr 
     
  
    For \bold{code written in C}, the calling sequence for \code{func} must be as in the following
    example:
    
    \code{void myderivs(int *neq, double *t, double *y, double *ydot,}\cr
    \code{double *yout, int *ip)}\cr
    \code{\{}\cr
      \code{if (ip[0] <1) error("nout should be at least 1");}\cr
      \code{ydot[0] = -k1*y[0] + k2*y[1]*y[2];}\cr
      \code{ydot[2] = k3 * y[1]*y[1];}\cr
      \code{ydot[1] = -ydot[0]-ydot[2];}\cr
      
      \code{yout[0] = y[0]+y[1]+y[2];}\cr
      
    \code{\}}\cr
  
    where \code{*neq} is the number of equations, \code{*t} is the value
    of the independent variable, \code{y} points to a double precision
    array of length \code{*neq} that contains the current value of the
    state variables, and \code{ydot} points to an array that will contain
    the calculated derivatives.
    \cr \code{yout} points to a double precision vector whose first \code{nout} values are other output variables (different from the state variables y), 
    and the next values are double precision values as passed by parameter \code{rpar} when calling the integrator. The key to the elements of yout is set in \code{*ip}
    \cr \code{*ip} points to an integer vector whose length is at least 3; the first element (IP[0]) contains the number of output values (which should be equal or larger than \code{nout}), 
    its second element contains the length of \code{*yout}, and the third element contains the length of \code{*ip}; 
    next are integer values, as passed by parameter \code{ipar} when calling the integrator. 
    \cr Note that the we start by checking wheter enough room is allocated for the output variables, else an error is passed to R and the integration is stopped.
    
     
    Note that the presence of \code{*yout} and \code{*ip} differs from previous version of lsoda (package odesolve). 
    Here there was no provision for exchanging quantities that are not directly solutions of the odes. 
    
    In the example above, parameters are kept in a global variable in the C
    code declared as
  
    \code{static double parms[3];}
  
    \code{#define} statements are used to make the code more readable, as
    in \code{#define k1 parms[0]}
    
    This is the only way to pass parameters to a compiled C function
    from the calling R code.  Functions that use this mechanism must be
    accompanied by a function that intialises the parameter values and protects the
    memory. This function is defined by \code{initfunc} and has as its sole argument a pointer to a function (see
    declaration below) that fills a double array with double precision
    values, to copy the parameter
    values into the global variable.  In the example here, a function such as:
    
    \code{void myinit(void (* odeparms)(int *, double *))}\cr
    \code{\{}\cr
      \code{int N=3;}\cr
      \code{odeparms(&N, parms);}\cr
    \code{\}}\cr
  
  will be required to initialize the parameter vector.  Here
    \code{myinit} just calls \code{odeparms} with a pointer to a int
    that contains the dimension of the parameter vector, and a pointer to
    the array that will contain the parameter values.
    
    Models may also be defined in \bold{Fortran}.  For example:
    
    \code{subroutine myderivs (neq, t, y, ydot, yout, ip)}\cr
    \code{double precision t, y, ydot, k1, k2, k3}\cr
    \code{integer neq, ip(*)}\cr
    \code{dimension y(3), ydot(3), yout(*)}\cr
    \code{common /myparms/k1,k2,k3}\cr
      
    \code{if(ip(1) < 1) call rexit("nout should be at least 1")}\cr
  
    \code{ydot(1) = -k1*y(1) + k2*y(2)*y(3)}\cr
    \code{ydot(3) = k3*y(2)*y(2)}\cr
    \code{ydot(2) = -ydot(1) - ydot(3)}\cr
      
    \code{yout(1) = y(1) + y(2) + y(3)}\cr    
    \code{return}\cr
    \code{end}\cr
   
    In Fortran, parameters may be stored in a common block, in which
    case the file that contains the model function definition must also
    contain an initialisation subroutine, again set by \code{initfunc}
    
    \code{subroutine myinit(odeparms)}\cr
    \code{external odeparms}\cr
    \code{integer N}\cr
    \code{double precision parms(3)}\cr
    \code{common /myparms/parms}\cr
      
    \code{N = 3}\cr
    \code{call odeparms(N, parms)}\cr
    \code{return}\cr
    \code{end}\cr
  
    If it is desired to supply a jacobian to the solver, then the jacobian
    must be defined in compiled code if the ODE system is.  The C function
    call for such a function must be as in the following example:
    
    \code{void myjac(int *neq, double *t, double *y, int *ml,}\cr
    \code{           int *mu, double *pd, int *nrowpd, double *yout, int *ip)}\cr
    \code{\{}\cr
      \code{pd[0] = -k1;}\cr
      \code{pd[1] = k1;}\cr
      \code{pd[2] = 0.0;}\cr
      \code{pd[(*nrowpd)] = k2*y[2];}\cr
      \code{pd[(*nrowpd) + 1] = -k2*y[2] - 2*k3*y[1];}\cr
      \code{pd[(*nrowpd) + 2] = 2*k3*y[1];}\cr
      \code{pd[(*nrowpd)*2] = k2*y[1];}\cr
      \code{pd[2*(*nrowpd) + 1] = -k2 * y[1];}\cr
      \code{pd[2*(*nrowpd) + 2] = 0.0;}\cr
    \code{\}}\cr
  
  The corresponding subroutine in Fortran is:
        
  \code{subroutine myjac (neq, t, y, ml, mu, pd, nrowpd, yout, ip)}\cr
  \code{integer neq, ml, mu, nrowpd, ip}\cr
  \code{double precision y(*), pd(nrowpd,*), yout(*), t, k1, k2, k3}\cr
  \code{common /myparms/k1, k2, k3}\cr
  
  \code{pd(1,1) = -k1}\cr
  \code{pd(2,1) =  k1}\cr
  \code{pd(3,1) = 0.0}\cr
  \code{pd(1,2) = k2*y(3)}\cr
  \code{pd(2,2) = -k2*y(3) - 2*k3*y(2)}\cr
  \code{pd(3,2) = 2*k3*y(2)}\cr
  \code{pd(1,3) = k2*y(2)}\cr
  \code{pd(2,3) = -k2*y(2)}\cr
  \code{pd(3,3) = 0.0}\cr
  
  \code{return}\cr
  \code{end}\cr
  
  To run the model using e.g. the Fortran code, the code in mymod.f must first be compiled.
  \cr This can be done in R itself:  
  
  \code{system("R CMD SHLIB mymod.f")}
  
  which will create file mymod.dll
  
  After loading the DLL, the model can be run:
  
  \code{dyn.load("mymod.dll")}\cr 
  \code{out <- lsoda(Y,times,"myderivs",parms=parms,rtol=RTOL,}\cr
  \code{atol=ATOL,mf=MF,jacfunc="myjac",dllname="mymod",}\cr
  \code{initfunc="myinit",verbose=FALSE,ynames=FALSE,nout=3 )}
  
  Examples in both C and Fortran are in the \file{dynload} subdirectory of
  the \code{deSolve} package directory.
  
      The output will have the \bold{attributes *istate*, and *rstate*}, two vectors with several useful elements. 
      \cr if \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen. 
  
      the following elements of \bold{istate} are meaningful:    
  
      \item el 1 : returns the conditions under which the last call to lsoda returned.
      2  if LSODA was successful, -1 if excess work done, -2 means excess accuracy requested. (Tolerances too small),
  -3 means illegal input detected. (See printed message.), -4 means repeated error test failures. (Check all input),
  -5 means repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.),
  -6 means error weight became zero during problem. (Solution component i vanished, and atol or atol(i) = 0.)
      \item el 12 : The number of steps taken for the problem so far.
      \item el 13 : The number of function evaluations for the problem so far.,
      \item el 14 : The number of Jacobian evaluations and LU decompositions so far.,
      \item el 15 : The method order last used (successfully).,
      \item el 16 : The order to be attempted on the next step.,
      \item el 17 : if el 1 =-4,-5: the largest component in the error vector,
      \item el 20 : The method indicator for the last succesful step, 1=adams (nonstiff), 2= bdf (stiff),
      \item el 21 : The current method indicator to be attempted on th next step, 1=adams (nonstiff), 2= bdf (stiff),
  
      \bold{rstate} contains the following:
      \item 1: The step size in t last used (successfully).
      \item 2: The step size to be attempted on the next step.
      \item 3: A tolerance scale factor, greater than 1.0, computed when a request for too much accuracy was detected.
      \item 4: the value of t at the time of the last method switch, if any.}
\seealso{\code{\link{lsode}}, \code{\link{lsodes}}, 
    \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}}.}
\note{The \file{demo} directory contains some examples of using
    \code{\link[nlme]{gnls}} to estimate parameters in a
    dynamic model.}
\keyword{math}

