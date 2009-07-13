### ============================================================================
### Interface to a generalized code for solving explicit variable and fixed
### step ODE solvers of the Runge-Kutta family, see helpfile for details.
### ============================================================================

rk <- function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, tcrit = NULL, hmin = 0, hmax = NULL, hini = hmax, ynames=TRUE,
  method = rkMethod("rk45dp7", ... ), maxsteps = 5000,
  dllname = NULL, initfunc=dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames=NULL, forcings=NULL,
  initforc = NULL, fcontrol=NULL, ...) {

    ## Check inputs
    hmax <- checkInput(y, times, func, rtol, atol,
        jacfunc=NULL, tcrit, hmin, hmax, hini, dllname)
    if (hmax == 0) hmax <- .Machine$double.xmax # i.e. practically unlimited

    n <- length(y)

    ## KS -> ThPe: maxsteps/tcrit checks are extra - should they be done in the other?
    if (maxsteps < 0)       stop("maxsteps must be positive")
    if (!is.finite(maxsteps)) maxsteps <- .Machine$integer.max
    if (is.character(method)) method <- rkMethod(method)
    if (is.null(tcrit)) tcrit <- max(times)

    ## Workaround for fixed step methods
    if (is.null(hini)) hini <- 0 # means that steps in "times" are used as they are

    ## Check if interpolation is switched off
    if (!is.null(method$interpolation) && !method$interpolation) {
      cat("\nMethod without or with disabled interpolation\n")
    } else {
      trange <- diff(range(times))
      ## ensure that we have at least 4..5 knots
      ## to allow 3rd order polynomial interpolation
      ## for methods without built-in dense output
      if ((is.null(method$d) &                             # has no "dense output"?
        (hmax > 0.2 * trange))) {                          # time steps too large?
        hini <- hmax <- 0.2 * trange
        if (hmin < hini) hmin <- hini
        cat("\nNote: Method ", method$ID,
            " needs intermediate steps for interpolation\n")
        cat("hmax decreased to", hmax, "\n")
      }
    }

    ## Model as shared object (DLL)?
    Ynames <- attr(y,"names")
    Initfunc <- NULL
    flist    <-list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    Nstates <- length(y) # assume length of states is correct

    if (is.character(func)) {
      DLL <- checkDLL(func, NULL, dllname,
                      initfunc, verbose, nout, outnames)

      Initfunc <- DLL$ModelInit
      Func     <- DLL$Func
      Nglobal  <- DLL$Nglobal
      Nmtot    <- DLL$Nmtot

      if (! is.null(forcings))
        flist <- checkforcings(forcings, times, dllname, initforc, verbose, fcontrol)

      rho <- NULL
      if (is.null(ipar)) ipar <- 0
      if (is.null(rpar)) rpar <- 0

    } else {
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
      rho <- environment(func)

      ## func is overruled, either including ynames, or not
      ## This allows to pass the "..." arguments and the parameters
        if(ynames) {
         Func   <- function(time,state,parms){
           attr(state, "names") <- Ynames
           func(time, state, parms, ...)}
        } else {                            # no ynames...
         Func   <- function(time,state,parms)
           func(time, state, parms, ...)
        }

      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      FF <- checkFuncEuler(Func, times, y, parms, rho, Nstates)
      Nglobal <- FF$Nglobal
      Nmtot   <- FF$Nmtot
    }

    ## handle length of atol and rtol
    if (Nstates %% length(atol))
      warning("length of atol does not match number of states")
    if (Nstates %% length(rtol))
      warning("length of rtol does not match number of states")

    atol <- rep(atol, length.out = Nstates)
    rtol <- rep(rtol, length.out = Nstates)

    ## Number of steps until the solver gives up
    nsteps  <- min(.Machine$integer.max, maxsteps * length(times))
    varstep <- method$varstep
    vrb <- FALSE # TRUE would force internal debugging output of the C code

    ## KS -> Thomas: still need to pass flist
    if (varstep) {                        # methods with variable step size
      out <- .Call("call_rkAuto", as.double(y), as.double(times),
        Func, Initfunc, parms,
        as.integer(Nglobal), rho, as.double(atol),
        as.double(rtol), as.double(tcrit), as.integer(vrb),
        as.double(hmin), as.double(hmax), as.double(hini),
        as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps), flist)
    } else {                              # fixed step methods
      #cat("hini=", hini, "\n")  # !!! temporary workaround; fix this
      #hini <- 0                 # !!! temporary workaround; fix this
      out <- .Call("call_rkFixed", as.double(y), as.double(times),
        Func, Initfunc, parms,
        as.integer(Nglobal), rho,
        as.double(tcrit), as.integer(vrb),
        as.double(hini), as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps), flist)
    }

    ## saving results
    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1,12,13,15), iout = c(1:3, 18))

    attr(out, "type") <- "rk"
    if (verbose) diagnostics(out)
    return(out)
}
