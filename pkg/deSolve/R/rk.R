
## Generalized solver for Runge-Kutta methods with variable time step
## This function is internal and not intended for the end user
rkAuto <- function(
  y, times, func, parms, rtol = 1e-6, atol = 1e-6,
	tcrit = NULL, verbose = FALSE,
  hmin = 0, hmax = NULL, hini = 0,
  method = rkMethod("rk45f", ... ), maxsteps = 5000, ...) {

  stage <- method$stage
  A     <- method$A
  bb1   <- method$b1
  bb2   <- method$b2
  cc    <- method$c
  qerr  <- 1/method$Qerr
  FSAL  <- ifelse(is.null(method$FSAL), FALSE, method$FSAL)
  
  ## track essential internal information
  ## experimental! Do it similar like lsoda
  istate <- numeric(23)

  Nstates <- length(y)
  FF      <- matrix(0, nrow = Nstates, ncol = stage)

  steps <- 0
  y0   <- y
  out  <- c(times[1], y0)
  accept <- FALSE
  
  if (verbose) {
    cat("method=", method$ID, "\n")
    cat("hini=", hini, "\n")
  }
  t    <- min(times)
  tmax <- max(times, tcrit) # NULL is handled automatically by max()
  dt   <- min(hmax, hini)
  hmax <- min(hmax, tmax - t) # limit hmax to the time range (esp. if hmax = Inf)
  
  ## function evaluations
  while ((t + dt) <= tmax) {
    ## save one step if the coefficients allow this (first same as last)
    if (FSAL & accept) {
          j1 <- 2
          FF[, 1] <- FF[, stage]
        } else {
          j1 <- 1
    }
    for (j in j1:stage) {
      Fj <- 0
      k  <- 1 
      while (k < j) {
        if (k>0) Fj <- Fj + A[j, k] * FF[ ,k]  * dt
        k <- k + 1
      }
      FF[, j] <- func(t + dt * cc[j], y0 + Fj, parms, ...)[[1]]
    }
    ## Estimation of new values
    dy1  <- FF %*% bb1
    dy2  <- FF %*% bb2
    y1   <- y0 + dt * dy1
    y2   <- y0 + dt * dy2

    ## stepsize adjustment after Press et. al (2007), Numerical Recipes in C
    yabs  <- pmax(abs(y0), abs(y2))
    scal  <- atol + yabs * rtol
    delta <- abs(y2 - y1)

    err <- delta / scal
    err  <- err[is.finite(err)] # remove Inf, in case of atol == rtol == 0
    if (length(err) > 0) {
      ## Press (2007): maximum is fine ...
      err <- max(err)
      ## ... but he takes the euklidean norm
      #err <- sqrt(sum(err)^2 / length(err))
    } else {
      err  <- 1
    }

    if (err == 0) {
       accept <- TRUE
       dtnew <- hmax # hmax must not be Inf or NULL !!!
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " +++ \n")
    } else if (err < 1) {  # accept
       accept <- TRUE
       ## safety factor = 0.9
       dtnew  <- max(hmax, dt * 0.9 * (err ^ -qerr))
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " + \n")
    } else if (err > 1){  # reject
       accept <- FALSE
       dtnew  <- dt * max(0.9 * (err ^ -qerr), 0.2)
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " - \n")
    } else { # err == 1
       accept <- TRUE
       dtnew <- dt
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " = \n")
    }

    ## Final check
    if (dt < hmin) {
      accept <- TRUE
      if (verbose) cat("h < hmin \n")
      warning("h < hmin")
      istate[1] <- -2
      dtnew <- hmin
    }
    ## data storage. Store imprecise results too, but warn if h < hmin
    if (accept) {
      t   <- t + dt
      y0  <- y2
      out <- rbind(out, c(t, y0))
    }
    steps <- steps + 1
    if (steps > maxsteps) 
      stop("
        An excessive amount of work (> maxsteps ) was done, 
        but integration was not successful -         
        increase maxsteps, increase atol/rtol, check your equations
        or select an alternative algorithm.
        ")
    dt  <- min(dtnew, tmax - t)
    if (t >= tmax) break
  }
  ## attach essential internal information
  ## experimental! Codes similar like lsoda
  istate[12] <- steps                   # number of steps
  istate[13] <- steps * stage           # number of function evaluations
  istate[15] <- method$Qerr             # order of the method
  attr(out, "istate") <- istate
  out
}

## Generalized sover for Runge-Kutta methods with fixed time step
rkFixed <- function( y, times, func, parms, tcrit = NULL,
  verbose=FALSE, hini=0, method = rkMethod("rk4", ... ), ...) {
    
  stage <- method$stage
  A     <- method$A
  bb1   <- method$b1
  #bb2   <- method$b2
  cc    <- method$c
  qerr  <- 1/method$Qerr
  
  FF      <- matrix(0, nrow=length(y), ncol=stage)

  y0   <- y
  out  <- c(times[1], y0)
  if (verbose) {
    cat("method=", method$ID, "\n")
    cat("hini=", hini, "\n")
  }
  t    <- min(times)
  tmax <- max(times, tcrit)                  # NULL is handled automatically by max
  dt <- hini
  ## derive internal (!) time step
  times <- unique(c(seq(t, tmax, dt), tmax)) # last step may possibly be shorter
  if (!is.matrix(A)) {                       # "A" coefficients given as subdiagonal
    for (i in 1:(length(times) - 1)) {
      t  <- times[i]
      for (j in 1:stage) {
        if (j == 1) Fj <- 0 else Fj <- A[j] * FF[ ,j - 1]
        FF[, j] <- dt * func(t + dt * cc[j], y0 + Fj, parms, ...)[[1]]
      }
      dy <- FF %*% bb1
      y1 <- y0 + dy
      out<- rbind(out, c(times[i + 1], y1))
      y0 <- y1
    }
  } else {                                   # "A" coefficients as matrix, not well tested !
    for (i in 1:(length(times) - 1)) {
      t  <- times[i]
      for (j in 1:stage) {
        Fj <- 0
        k  <- 1 
        while (k < j) {
          if (k>0) Fj <- Fj + A[j, k] * FF[ ,k]
          k <- k + 1
        }
        FF[, j] <- dt * func(t + dt * cc[j], y0 + Fj, parms, ...)[[1]]
      }
    }
    ## Estimation of new values
    dy1  <- FF %*% bb1
    y1   <- y0 + dy1
    ## data storage
    out <- rbind(out, c(t, y0))
  } # end if
  out
}

rk <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
	tcrit = NULL, verbose=FALSE, hmin = 0, hmax = NULL, hini = hmax,
  method = rkMethod("rk45f", ... ), maxsteps = 5000, ...) {

### check input
    if (!is.numeric(y))       stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times)&&!is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function")
    if (!is.numeric(rtol))   stop("`rtol' must be numeric")
    if (!is.numeric(atol))   stop("`atol' must be numeric")
    if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    if (!is.numeric(hmin))   stop("`hmin' must be numeric")
    if (hmin < 0) stop ("`hmin' must be a non-negative value")
    if (is.null(hmax))
       hmax <- ifelse (is.null(times), 0, max(abs(diff(times))))
    if (!is.numeric(hmax))   stop("`hmax' must be numeric")
    if (hmax < 0)            stop ("`hmax' must be a non-negative value")
    if (hini < 0)            stop("`hini' must be a non-negative value")

    if (!is.numeric(y))     stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func)) stop("`func' must be a function")
    
    if (is.character(method)) method <- rkMethod(method)

    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    rho <- environment(func)
    tmp <- eval(func(times[1], y, parms, ...), rho)
    if (!is.list(tmp)) stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 "must equal the length of the initial conditions vector (",
                 length(y),")", sep="")
      )
    Nglobal <- if (length(tmp) > 1) length(unlist(tmp[-1])) else 0
    Nmtot   <- attr(unlist(tmp[-1]), "names")

    ## -----------------------------------------------------------------------
    varstep <- method$varstep
    if (varstep) { # methods with variable step size
      out <- rkAuto(y, times, func, parms, rtol = rtol, atol = atol, tcrit = tcrit,
               verbose = verbose, hmin = hmin, hmax = hmax, hini = hini, 
               method = method, maxsteps = maxsteps, ...)
    } else {       # fixed step methods
      out <- rkFixed(y, times, func, parms, tcrit = tcrit,
         verbose = verbose, hini = hini, method = method, ...)
    }
    istate <- attr(out, "istate") # remember internal information
    ## -----------------------------------------------------------------------
    nm <- c("time",
      if (!is.null(attr(y, "names"))) names(y) else as.character(1:n)
    )
    
    ## interpolation: simple linear approximation,
    ## future versions may use dense output for some rk variants
    m   <- ncol(out)
    res <- matrix(0, nrow = length(times), ncol = m)
    res[,1] <- times
    for (i in 2:m) {
      res[,i] <- as.vector(approx(out[,1], out[,i], times)$y)
    }
    out <- res
    ## external outputs
    if (Nglobal > 0) {
      out2 <- matrix(nrow = nrow(out), ncol = Nglobal)
      for (i in 1:nrow(out2))
        out2[i,] <- func(out[i, 1], out[i, -1], parms, ...)[-1]
      out <- cbind(out, out2)
      nm  <- c(nm,
        if (!is.null(Nmtot)) Nmtot else as.character((n + 1) : (n + Nglobal))
      )
    }
    ## interpolation
    #m   <- ncol(out)
    #res <- matrix(0, nrow = length(times), ncol = m)
    #res[,1] <- times
    #for (i in 2:m) {
    #  res[,i] <- as.vector(approx(out[,1], out[,i], times)$y)
    #}
    dimnames(out) <- list(NULL, nm)
    attr(out, "istate") <- istate
    out
}
