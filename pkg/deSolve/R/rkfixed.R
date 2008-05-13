### Classical Runge-Kutta-fixed-step-integration

rkfixed <- function(y, times, func, parms, h = diff(times), ...) {
    if (!is.numeric(y))     stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func)) stop("`func' must be a function")
    
    if (!(length(h) %in% c(1, length(times) - 1)))
      stop("length of h must be either one or length of times - 1")

    if (length(h) == 1) {
      tseq <- seq(min(times), max(times) + h, h)
      h    <- diff(tseq)
    } else {
      tseq <- times[1] + c(0, cumsum(h))
    }
    
    if (max(tseq) < max(times))
      stop("sequence of h does not span complete time range")  # reformulate this

    n <- length(y)

    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    rho <- environment(func)
    tmp <- eval(func(times[1], y, parms), rho)
    if (!is.list(tmp)) stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 "must equal the length of the initial conditions vector (",
                 length(y),")", sep="")
      )
    Nglobal <- if (length(tmp) > 1) length(unlist(tmp[-1])) else 0
    Nmtot   <- attr(unlist(tmp[-1]), "names")

    y0   <- y
    out  <- c(times[1], y0)
    for (i in 1:(length(tseq) - 1)) {
      t  <- tseq[i]
      dt <- h[i]
      F1 <- dt * func(t,        y0,            parms,...)[[1]]
      F2 <- dt * func(t + dt/2, y0 + 0.5 * F1, parms,...)[[1]]
      F3 <- dt * func(t + dt/2, y0 + 0.5 * F2, parms,...)[[1]]
      F4 <- dt * func(t + dt  , y0 + F3,       parms,...)[[1]]
      dy <- (F1 + 2 * F2 + 2 * F3 + F4)/6
      y1 <- y0 + dy
      out<- rbind(out, c(tseq[i + 1], y1))
      y0 <- y1
    }
    nm <- c("time",
      if (!is.null(attr(y, "names"))) names(y) else as.character(1:n)
    )
    if (Nglobal > 0) {
      out2 <- matrix(nrow=nrow(out), ncol=Nglobal)
      for (i in 1:nrow(out2))
        out2[i,] <- func(out[i,1], out[i,-1], parms, ...)[-1]
      out <- cbind(out, out2)
      nm <- c(nm,
        if (!is.null(Nmtot)) Nmtot else as.character((n + 1) : (n + Nglobal))
      )
    }
    ## interpolation
    m   <- ncol(out)
    res <- matrix(0, nrow = length(times), ncol = m)
    res[,1] <- times
    for (i in 2:m) {
      res[,i] <- as.vector(approx(out[,1], out[,i], times)$y)
    }
    
    dimnames(res) <- list(NULL, nm)
    res
}
