### ============================================================================
### S3 methods
### ============================================================================

print.deSolve <- function(x, ...)
  print(as.data.frame(x), ... )

### ============================================================================

plot.deSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {
    if (is.null(ask))
        ask <- prod(par("mfcol")) < length(which) && dev.interactive()
    t <- 1     # column with "times"
    var <- colnames(x)
    if (!is.numeric(which)) {
        ln <- length(which)
        Select <- which(var %in% which)
        if (length(Select) != ln)
            stop("not all variables in 'which' are in 'x'")
        which <- Select
    }
    else {
        which <- which + 1  # "which now refers to the column number
        if (max(which) > ncol(x))
            stop("index in 'which' too large")
        if (min(which) < 1)
            stop("index in 'which' should be > 0")
    }
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
        nc <- min(ceiling(sqrt(np)), 3)  # assume max 3 x 3 panels
        nr <- min(ceiling(np/nc), 3)
        mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else
        mfrow <- dots$mfrow

    if (!is.null(mfrow)) {
        mf <- par(mfrow = mfrow)
    }
    ## interactively wait if there are remaining figures
    if (ask) {
        oask <- devAskNewPage(TRUE)
	on.exit(devAskNewPage(oask))
    }

    Main <- is.null(dots$main)

    xxlab <- if (is.null(dots$xlab))  colnames(x)[t]  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""              else dots$ylab
    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
        if (Main)
            dots$main <- colnames(x)[i]
            dots$xlab <- xxlab[i-1]
            dots$ylab <- yylab[i-1]
        do.call("plot", c(alist(x[, t], x[, i]), dots))
    }
}
