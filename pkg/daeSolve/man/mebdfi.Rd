\name{mebdfi}
\alias{mebdfi}
\title{Solver for Differential Algebraic Equations (DAE) up to index 3}
\description{
  Solves either:
  \itemize{
    \item a system of ordinary differential equations (ODE) of the form
      \deqn{y' = f(t,y,...)} or
    \item a system of differential algebraic equations (DAE) of the form
      \deqn{F(t,y,y') = 0}
  }

  using the Modified Extended Backward Differentiation formulas for stiff
  fully implicit inital value problems

  These formulas increase the absolute stability regions of
  the classical BDFs.

  The orders of the implemented formulae range from 1 to 8.
  
  The \R function \code{mebdfi} provides an interface to the Fortran DAE
  solver of the same name, written by T.J. Abdulla and J.R. Cash.
    
  The system of DE's is written as an \R function or can be defined in
  compiled code that has been dynamically loaded.
}
\usage{
mebdfi(y, times, func = NULL, parms, dy = NULL, res = NULL,
  nind=c(length(y),0,0), rtol = 1e-6, atol = 1e-8, jacfunc = NULL,
  jacres = NULL, jactype = "fullint", verbose = FALSE,
  tcrit = NULL,  hini = 0, ynames = TRUE, maxord = 7, bandup = NULL,
  banddown = NULL, maxsteps = 5000, dllname = NULL,
  initfunc = dllname, initpar = parms, rpar = NULL,
  ipar = NULL, nout = 0, outnames = NULL,
  forcings=NULL, initforc = NULL, fcontrol=NULL,  ...)
}

\arguments{
  \item{y }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set \code{times} = \code{NULL}.
  }
  \item{func }{cannot be used if the model is a DAE system. If an ODE
    system, \code{func} should be an \R-function that computes the
    values of the derivatives in the ODE system (the \emph{model
    definition}) at time t.
  
    \code{func} must be defined as: \code{func <- function(t, y,
    parms,...)}.  \cr \code{t} is the current time point in the
    integration, \code{y} is the current estimate of the variables in
    the ODE system.  If the initial values \code{y} has a \code{names}
    attribute, the names will be available inside \code{func}, unless
    \code{ynames} is FALSE. \code{parms} is a vector or list of
    parameters.  \code{...} (optional) are any other arguments passed to
    the function.
    
    The return value of \code{func} should be a list,
    whose first element is a vector containing the derivatives of
    \code{y} with respect to \code{time}, and whose next elements are
    global values that are required at each point in \code{times}.
    The derivatives should be specified in the same order as the specification
    of the state variables \code{y}.

    Note that it is not possible to define \code{func} as a compiled
    function in a dynamically loaded shared library. Use \code{res}
    instead.
  }
  \item{parms }{vector or list of parameters used in \code{func},
    \code{jacfunc}, or \code{res}
  }
  \item{dy }{the initial derivatives of the state variables of the DE
    system.  Ignored if an ODE.
  }
  \item{res }{if a DAE system: either an \R-function that computes the
    residual function F(t,y,y') of the DAE system (the model
    defininition) at time \code{t}, or a character string giving the
    name of a compiled function in a dynamically loaded shared library.
      
    If \code{res} is a user-supplied R-function, it must be defined as:
    \code{res <- function(t, y, dy, parms, ...)}.
    
    Here \code{t} is the current time point in the integration, \code{y}
    is the current estimate of the variables in the ODE system,
    \code{dy} are the corresponding derivatives.  If the initial
    \code{y} or \code{dy} have a \code{names} attribute, the names will be
    available inside \code{res}, unless \code{ynames} is \code{FALSE}.
    \code{parms} is a vector of parameters.
    
    The return value of \code{res} should be a list, whose first element
    is a vector containing the residuals of the DAE system,
    i.e. delta = F(t,y,y'), and whose next elements contain output
    variables that are required at each point in \code{times}.
  
    If \code{res} is a string, then \code{dllname} must give the name of
    the shared library (without extension) which must be loaded before
    \code{mebdfi()} is called (see package vignette \code{"compiledCode"}
    for more information).
  }
  \item{nind }{if a DAE system: a three-valued vector with the number of
    variables of index 1, 2, 3 respectively.
    The equations must be defined such that the index 1 variables precede
    the index 2 variables which in turn precede the index 3 variables.
    The sum of the variables of different index should equal N,
    the total number of variables.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function that computes the
    Jacobian of the system of differential equations. Only used in case
    the system is an ODE (y' = f(t,y)), specified by \code{func}. The \R
    calling sequence for \code{jacfunc} is identical to that of
    \code{func}.
    
    If the Jacobian is a full matrix, \code{jacfunc} should return a
    matrix dydot/dy, where the ith row contains the derivative of
    \eqn{dy_i/dt} with respect to \eqn{y_j}, or a vector containing the
    matrix elements by columns (the way \R and Fortran store matrices).
    
    If the Jacobian is banded, \code{jacfunc} should return a matrix
    containing only the nonzero bands of the Jacobian, rotated
    row-wise. See first example of lsode.
  }
  \item{jacres }{ \code{jacres} and not \code{jacfunc} should be used if
    the system is specified by the residual function F(t,y,y'),
    i.e. \code{jacres} is used in conjunction with \code{res}.
    
    If \code{jacres} is an \R-function, the calling sequence for
    \code{jacres} is identical to that of \code{res}, but with extra
    parameter \code{cj}.  Thus it should be defined as: \code{jacres <-
    function(t, y, dy, parms, cj, ...)}.  Here \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system, \eqn{y'} are the corresponding derivatives
    and \code{cj} is a scalar, which is normally proportional to
    the inverse of the stepsize. If the initial \code{y} or \code{dy}
    have a \code{names} attribute, the names will be available inside
    \code{jacres}, unless
    \code{ynames} is \code{FALSE}. \code{parms} is a vector of
    parameters (which may have a names attribute).
    
    If the Jacobian is a full matrix, \code{jacres} should return the
    matrix dG/dy + cj*dG/dyprime, where the ith row is the sum of the
    derivatives of \eqn{G_i} with respect to \eqn{y_j} and the scaled
    derivatives of \eqn{G_i} with respect to \eqn{dy_j}.
    
    If the Jacobian is banded, \code{jacres} should return only the
    nonzero bands of the Jacobian, rotated rowwise. See details for the
    calling sequence when \code{jacres} is a string.
  }
  \item{jactype }{the structure of the Jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by the user.
  }
  \item{verbose }{if TRUE: full output to the screen, e.g. will
    print the \code{diagnostiscs} of the integration - see details.
  }
  \item{tcrit }{the Fortran routine \code{mebdfi} overshoots its targets
    (times points in the vector \code{times}), and interpolates values
    for the desired time points.  If there is a time beyond which
    integration should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is set to 1e-6, but it may be better to set it equal to \code{rtol}.
    The solver is quite sensitive to values of \code{hini}; sometimes
    if it fails, it helps to decrease/increase \code{hini}
  }
  \item{ynames }{logical, if \code{FALSE} names of state variables are not
    passed to function \code{func}; this may speed up the simulation especially
    for large models.
  }
  \item{maxord }{the maximum order to be allowed, an integer between 2 and 7.
    The default is \code{maxord} = 7, but values of 4-5 may be better for
    difficult problems; hihger order methods are more efficient but less stable.
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case
    the Jacobian is banded (and \code{jactype} one of
    "bandint","bandusr")
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case
    the Jacobian is banded (and \code{jactype} one of
    "bandint","bandusr")
  }
  \item{maxsteps }{maximal number of steps per output interval taken by the
    solver; will be recalculated to be at least 500 and a multiple of
    500; the solver will give a warning if more than 500 steps are
    taken, but it will continue till \code{maxsteps} steps.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions referred to in \code{res} and
    \code{jacres}.  See package vignette \code{"compiledCode"}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"compiledCode"}.
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the dll: the
    parameters passed to the initialiser, to initialise the common
    blocks (fortran) or global variables (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the dll-functions whose names are
    specified by \code{res} and \code{jacres}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{res} and \code{jacres}.
  }
  \item{nout }{only used if \file{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{res}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculed in the dll - you have
    to perform this check in the code - See package vignette
    \code{"compiledCode"}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{res}, present in the shared library.
    These names will be used to label the output matrix.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme.

    See package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    vignette \code{compiledCode} from package \code{deSolve}.
  }
  \item{... }{additional arguments passed to \code{func},
    \code{jacfunc}, \code{res} and \code{jacres}, allowing this to be a
    generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{y} plus the number of "global" values
  returned in the next elements of the return from \code{func} or
  \code{res}, plus an additional column (the first) for the time value.
  There will be one row for each element in \code{times} unless the
  Fortran routine `mebdfi' returns with an unrecoverable error.  If
  \code{y} has a names attribute, it will be used to label the columns
  of the output value.
      
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## =======================================================================
## Coupled chemical reactions including an equilibrium
## modeled as (1) an ODE and (2) as a DAE
##
## The model describes three chemical species A,B,D:
## subjected to equilibrium reaction D <- > A + B
## D is produced at a constant rate, prod
## B is consumed at 1s-t order rate, r
## Chemical problem formulation 1: ODE
## =======================================================================

## Dissociation constant
K <- 1 

## parameters
pars <- c(
        ka   = 1e6,     # forward rate
        r    = 1,
        prod = 0.1)


Fun_ODE <- function (t, y, pars)
{
  with (as.list(c(y, pars)), {
    ra  <- ka*D        # forward rate
    rb  <- ka/K *A*B   # backward rate

    ## rates of changes
    dD  <- -ra + rb + prod
    dA  <-  ra - rb
    dB  <-  ra - rb - r*B
    return(list(dy = c(dA, dB, dD),
                CONC = A+B+D))
  })
}

## =======================================================================
## Chemical problem formulation 2: DAE
## 1. get rid of the fast reactions ra and rb by taking
## linear combinations   : dD+dA = prod (res1) and
##                         dB-dA = -r*B (res2)
## 2. In addition, the equilibrium condition (eq) reads:
## as ra = rb : ka*D = ka/K*A*B = >      K*D = A*B
## =======================================================================

Res_DAE <- function (t, y, yprime, pars)
{
  with (as.list(c(y, yprime, pars)), {
    ## residuals of lumped rates of changes
    res1 <- -dD - dA + prod
    res2 <- -dB + dA - r*B
    
    ## and the equilibrium equation
    eq   <- K*D - A*B

    return(list(c(res1, res2, eq),
                CONC = A+B+D))
  })
}

times <- seq(0, 100, by = 1)

## Initial conc; D is in equilibrium with A,B
y     <- c(A = 2, B = 3, D = 2*3/K)

## ODE model solved with mebdfi
ODE <- as.data.frame(mebdfi(y = y, times = times, func = Fun_ODE,
                     parms = pars, atol = 1e-8, rtol = 1e-8))

## Initial rate of change
dy  <- c(dA = 0, dB = 0, dD = 0) 
## DAE model solved with mebdfi
DAE <- as.data.frame(mebdfi(y = y, dy = dy, times = times,  
         res = Res_DAE, parms = pars, atol = 1e-8, rtol = 1e-8))

## ================
## plotting output
## ================
opa <- par(mfrow = c(2,2))

for (i in 2:5)
{
plot(ODE$time,ODE[,i],xlab = "time",
     ylab = "conc",main = names(ODE)[i],type = "l")
points(DAE$time,DAE[,i],col = "red")
}
legend("bottomright",lty = c(1,NA),pch = c(NA,1),
       col = c("black","red"),legend = c("ODE","DAE"))      
       
# difference between both implementations:
max(abs(ODE-DAE))

par(mfrow = opa)

## =============================================================================
## Car axis problem, index 3 DAE, 8 differential, 2 algebraic equations
## from
## F. Mazzia and C. Magherini. Test Set for Initial Value Problem Solvers,
## release 2.4. Department
## of Mathematics, University of Bari and INdAM, Research Unit of Bari,
## February 2008.
## Available at http://www.dm.uniba.it/~testset.
## =============================================================================

car <- function(t,y,dy,pars){
  with(as.list(c(pars,y)), {
      f <- rep(0,10)

      yb  <- r*sin(w*t)
      xb  <- sqrt(L*L-yb*yb)
      Ll  <- sqrt(xl^2+yl^2)
      Lr  <- sqrt((xr-xb)^2+(yr-yb)^2)

      f[1:4] <- y[5:8]
      k <- M*eps*eps/2

      f[5]  <- (L0-Ll)*xl/Ll +lam1*xb+2*lam2*(xl-xr)
      f[6]  <- (L0-Ll)*yl/Ll +lam1*yb+2*lam2*(yl-yr)-k*g
      f[7]  <- (L0-Lr)*(xr-xb)/Lr -2*lam2*(xl-xr)
      f[8]  <- (L0-Lr)*(yr-yb)/Lr -2*lam2*(yl-yr)-k*g

      f[9]  <- xb*xl+yb*yl
      f[10] <- (xl-xr)^2+(yl-yr)^2-L*L

      delt       <- dy-f
      delt[5:8]  <- k*dy[5:8]-f[5:8]
      delt[9:10] <- -f[9:10]

      list(delt=delt,f=f)
  })
}


# parameters
pars <- c(eps = 1e-2, M = 10, L = 1, L0 = 0.5,
          r   = 0.1,  w = 10, g = 1)

# initial conditions: state variables
yini <-  with (as.list(pars),
   c(xl=0, yl=L0, xr=L, yr=L0, xla=-L0/L,
     yla=0, xra=-L0/L, yra=0, lam1=0, lam2=0)
              )

# initial conditions: derivates
dyini <- rep(0,10)
FF    <- car(0,yini,dyini,pars)
dyini[1:4] <- yini[5:8]
dyini[5:8] <- 2/pars["M"]/(pars["eps"])^2*FF$f[5:8]

# check consistency of initial condition: delt should be = 0.
car(0,yini,dyini,pars)

# running the model
times <- seq(0,3,by=0.01)
nind  <- c(4,4,2)   # index 1, 2 and 3 variables
out   <- mebdfi(y=yini,dy=dyini,times,res=car,parms=pars,
                nind=nind, rtol=1e-5,atol=1e-5)

mf<-  par(mfrow=c(2,2))
plot(out[,1], out[,2], type="l", lwd=2)
plot(out[,1], out[,3], type="l", lwd=2)
plot(out[,1], out[,4], type="l", lwd=2)
plot(out[,1], out[,5], type="l", lwd=2)
mtext(outer=TRUE,side=3,line=-0.5,cex=1.5,"car axis")
par(mfrow=mf)
}
\references{
J. R. Cash, The integration of stiff initial value problems
         in O.D.E.S using modified extended backward differentiation
         formulae, Comp. and Maths. with applics., 9, 645-657, (1983).

J.R. Cash and S. Considine, an MEBDF code for stiff
         initial value problems, ACM Trans Math Software, 142-158,
         (1992).

J.R. Cash, Stable recursions with applications to the
         numerical solution of stiff systems, Academic Press,(1979).
}
\details{
  The mebdfi solver uses modified extended backward differentiation
  formulas of orders one through eight (specified with \code{maxord})
  to solve either:
  \itemize{
    \item an ODE system of the form \deqn{y' = f(t,y,...)}  for y = Y,
    or
    \item a DAE system of the form \deqn{F(t,y,y') = 0} for y = Y and y'
    = YPRIME.
  }
  

  The recommended value of \code{maxord} is eight, unless it is believed
  that there are severe stability problems in which case \code{maxord} =
  4 or 5 should be tried instead.

  ODEs are specified in \code{func}, DAEs are specified in \code{res}.
  
  If a DAE system, Values for Y \emph{and} YPRIME at the initial time
  must be given as input. Ideally,these values should be consistent,
  that is, if T, Y, YPRIME are the given initial values, they should
  satisfy F(T,Y,YPRIME) = 0.
  
  The form of the \bold{Jacobian} can be specified by
  \code{jactype}. This is one of:
  
  \describe{
    \item{jactype = "fullint":}{a full Jacobian, calculated internally
      by \code{mebdfi}, the default,
    }
    \item{jactype = "fullusr":}{a full Jacobian, specified by user
      function \code{jacfunc} or \code{jacres},
    }
    \item{jactype = "bandusr":}{a banded Jacobian, specified by user
      function \code{jacfunc} or \code{jacres}; the size of the bands
      specified by \code{bandup} and \code{banddown},
    }
    \item{jactype = "bandint":}{a banded Jacobian, calculated by
      \code{mebdfi}; the size of the bands specified by \code{bandup} and
      \code{banddown}.
    }
  }
  
  If \code{jactype} = "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc}.

  If jactype = "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc} or \code{jacres}.
  
  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver.  If the request for
  precision exceeds the capabilities of the machine, mebdfi will return
  an error code.
  
  \bold{res and jacres} may be defined in compiled C or Fortran code, as
  well as in an R-function. See deSolve's vignette \code{"compiledCode"}
  for details.  Examples
  in Fortran are in the \file{dynload} subdirectory of the
  \code{daeSolve} package directory.

  The diagnostics of the integration can be printed to screen
  by calling \code{\link{diagnostics}}. If \code{verbose} = \code{TRUE},
  the diagnostics will written to the screen at the end of the integration.

  See vignette("deSolve") for an explanation of each element in the vectors
  containing the diagnostic properties and how to directly access them.

}
\seealso{
  \itemize{
    \item \code{\link{daspk}} another DAE solver,
  }

  \code{\link{diagnostics}} to print diagnostic messages.
}
\note{

}
\keyword{math}

