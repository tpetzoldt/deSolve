\name{checkDLL}
\alias{checkDLL}

\title{
Check DLL/.so of compiled model.
}
\description{
Check DLL/.so of a compiled model and create a list of symbols.
}
\usage{
checkDLL(func, jacfunc, dllname, initfunc, verbose, nout, outnames, JT = 1)
}

\arguments{
  \item{func}{
 character: name of the derivative function.
}
  \item{jacfunc}{
%%     ~~Describe \code{jacfunc} here~~
}
  \item{dllname}{
%%     ~~Describe \code{dllname} here~~
}
  \item{initfunc}{
%%     ~~Describe \code{initfunc} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{nout}{
%%     ~~Describe \code{nout} here~~
}
  \item{outnames}{
%%     ~~Describe \code{outnames} here~~
}
  \item{JT}{
%%     ~~Describe \code{JT} here~~
}
}
\value{
List of symbols .....
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\seealso{
\code{\link{lsoda}}
}
\examples{
\dontrun{
symbols <- checkDLL(func = "derivs", jacfunc = NULL, dllname = "lorenzc",
             initfunc = "initmod", verbose = TRUE, nout = 0,
             outnames = NULL, JT = 1)
}

}
\keyword{math}
