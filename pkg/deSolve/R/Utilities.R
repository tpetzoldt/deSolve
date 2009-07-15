### ============================================================================
### S3 methods
### ============================================================================

print.deSolve <- function(x,...)
  print(as.data.frame(x))

### ============================================================================

plot.deSolve <- function (x, t = 1, which = 2:ncol(x), ...)
{
    var <- colnames(x)
    if (length(t) != 1)
        stop("'t' should contain one value")
    if (!is.numeric(t)) {
        Select <- which(var %in% t)
        t <- Select
    }
    else {
        if (max(t) > ncol(x))
            stop("index in 't' too large")
        if (min(t) < 1)
            stop("index in 't' should be >0")
    }
    if (!is.numeric(which)) {
        ln <- length(which)
        Select <- which(var %in% which)
        if (length(Select) != ln)
            stop("not all variables in 'which' are in 'x'")
        which <- Select
    }
    else {
        if (max(which) > ncol(x))
            stop("index in 'which' too large")
        if (min(which) < 1)
            stop("index in 'which' should be >0")
    }
    np <- length(which)
    dots <- list(...)
    nmdots <- names(dots)
    if (!"mfrow" %in% nmdots) {
        nc <- ceiling(sqrt(np))
        nr <- ceiling(np/nc)
        mfrow <- c(nr, nc)
    }
    else mfrow <- dots$mfrow
    if (!is.null(mfrow)) {
        mf <- par(mfrow = mfrow)
    }
    Main <- is.null(dots$main)
    dots$xlab <- if (is.null(dots$xlab))
        colnames(x)[t]
    else dots$xlab
    dots$ylab <- if (is.null(dots$ylab))
        ""
    else dots$ylab
    for (i in which) {
        if (Main)
            dots$main <- colnames(x)[i]
        do.call("plot", c(alist(x[, t], x[, i]), dots))
    }
}
