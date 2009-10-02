\name{plot.daeSolve}
\alias{plot.daeSolve}
\alias{hist.daeSolve}
\alias{print.daeSolve}
\title{
 Plot and Histogram Method for daeSolve Objects
}
\description{
Plot the output of numeric integration routines.
}
\usage{
\method{plot}{daeSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
\method{hist}{daeSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
\method{print}{daeSolve}(x, \dots)
}

\arguments{
  \item{x }{an object of class \code{daeSolve}, as returned by the
    integrators, and to be plotted.}

  \item{which }{the name(s) or the index to the variables that should be
    plotted.  Default = all variables. }

  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.}
  
  \item{\dots}{additional graphics arguments passed to
    \code{\link{plot.default}} or \code{\link{hist}}}
}

\details{
  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow=c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.

  Other graphical parameters can be passed as well. Parameters
  \code{xlab} and \code{ylab} are vectorized, so it is possible to
  assign specific axis labels to individual plots.
}

\seealso{
}

\keyword{ hplot }

