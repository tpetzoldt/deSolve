\name{forcings}
\alias{forcings}
\title{
  Passing forcing functions to models written in compiled code.
}
\description{
  Models may be defined in compiled C or Fortran code, as well as
  in an R-function.

  How to write compiled code is described in package vignette
  \code{"compiledCode"}, which should be referred to for details.

  This vignette also contains examples on how to pass forcing functions.
  
  If the models are defined in compiled C or Fortran code, it is possible to
  use \code{deSolve}s \bold{forcing function} update algorithm. This is the
  compiled-code equivalent of \code{approxfun} or \code{approx}, which
  estimates the values of external variables ("forcings") at each time point
  by interpolation of a data series.
}

\details{

  The forcing function data series is
  privided by means of argument \code{forcings}, \code{initforc} is
  the name of the forcing function initialisation function, as provided in
  \file{dllname}, while \code{fcontrol} is a list used to finetune how the
  forcing update should be performed.

  The \bold{fcontrol} argument is a list that can supply any of the following
  components (conform the definitions in the \link{approx} function):
  \describe{
    \item{method }{specifies the interpolation method to be used.
      Choices are "linear" or "constant",}
    \item{rule }{an integer describing how interpolation is to take place
      outside the interval [min(times), max(times)].
      If \code{rule} is \code{1} then an error will be triggered and the
      calculation will stop if \code{times} extends the interval of the
      forcing function data set. If it is \code{2}, the *default*, the
      value at the closest data extreme is used, a warning will be printed if
      \code{verbose} is TRUE,

      Note that the default differs from the \code{approx} default}
    \item{f }{For method=\code{"constant"} a number between \code{0} and
      \code{1} inclusive, indicating a compromise between left- and
      right-continuous step functions. If \code{y0} and \code{y1} are the
      values to the left and right of the point then the value is
      \code{y0*(1-f)+y1*f} so that \code{f=0} is right-continuous and
      \code{f=1} is left-continuous,
      }
    \item{ties }{Handling of tied \code{times} values. Either a function
      with a single vector argument returning a single number result or the
      string "ordered".

      Note that the default is "ordered", hence the existence of ties will
      NOT be investigated; in the \code{C} code this will mean that -if ties
      exist, the first value will be used; if the dataset is not ordered,
      then nonsense will be produced.

      Alternative values for \code{ties} are \code{mean}, \code{min} etc
      }
   }
   The defaults are:

   \code{fcontrol=list(method="linear", rule = 2,  f = 0, ties = "ordered")}

   Note that only ONE specification is allowed, even if there is more than
   one forcing function data set.

   More information about models defined in compiled code is in the package
   vignette ("compiledCode").
}
\author{
  Karline Soetaert,
  
  Thomas Petzoldt,

  R. Woodrow Setzer (Maintainer)
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(aquaphy)

## show package vignette with tutorial about how to use compiled models
## + source code of the vignette
## + directory with C and FORTRAN sources
vignette("compiledCode")
edit(vignette("compiledCode"))
browseURL(paste(system.file(package = "deSolve"), "/doc", sep = ""))
}
}

\keyword{utilities}