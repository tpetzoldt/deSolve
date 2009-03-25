\name{diagnostics}
\alias{diagnostics}
\title{Print diagnostic characteristics of ODE Solvers.}
\description{
  Prints several diagnostics of the simulation to the screen, e.g.
  number of steps taken, the last step size,...
}
\usage{
diagnostics(obj)
}
\arguments{
  \item{obj }{is the output matrix as produced by one of the
    integration routines (except for runge-kutta methods).
  }
}
\value{
  none
}
\details{
  When the integration output is saved as a data.frame, then the required
  attributes are lost and method \code{diagnostics} will not work anymore.
}
\examples{
## The famous Lorenz equations: chaos in the earth's atmosphere
## Lorenz 1963. J. Atmos. Sci. 20, 130-141.

chaos <- function(t, state, parameters) {
  with(as.list(c(state)),{

    dx     <- -8/3*x+y*z
    dy     <- -10*(y-z)
    dz     <- -x*y+28*y-z

    list(c(dx, dy, dz))
  })
}

state <- c(x = 1, y = 1, z = 1)
times <- seq(0, 50, 0.01)
out   <- vode(state, times, chaos, 0)
pairs(out,pch=".")
diagnostics(out)
## the following is for testing only; remove it in the release version
\dontrun{
out.rk   <- rk(state, times, chaos, 0, method="ode45", maxsteps=1e5,atol=1e-4, rtol=1e-4)
pairs(out.rk, pch=".")
diagnostics(out.rk)

out.rkc   <- rk_c(state, times, chaos, 0, method="ode45", hmax=0.1, maxsteps=1e5)
pairs(out.rkc, pch=".")
diagnostics(out.rkc)
}
}
\keyword{ utilities }