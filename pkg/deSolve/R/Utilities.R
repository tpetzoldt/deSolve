### ============================================================================
### first some common functions
### ============================================================================

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)),3)
      nr <- min(ceiling(nv/nc),3)
      mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow)) {
      mf <- par(mfrow=mfrow)
    }

   ## interactively wait if there are remaining figures
    if (is.null(ask))
      ask <- prod(par("mfrow")) < nv && dev.interactive()

    return(ask)
}

## =============================================================================
## find a variable
## =============================================================================

selectvar <- function (which, var, NAallowed = FALSE) {

    if (!is.numeric(which)) {
        ln <- length(which)
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i]==var)
          if (length(ss) ==0 & ! NAallowed) 
            stop("variable ", which[i], " not in var")
          else if (length(ss) == 0)
            Select <- c(Select,NA)
          else
            Select <- c(Select,ss)
        }
    }
    else {
        Select <- which + 1  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
}

### ============================================================================
### S3 methods
### ============================================================================

print.deSolve <- function(x, ...)
  print(as.data.frame(x), ... )

### ============================================================================

hist.deSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {
    t     <- 1     # column with "times"
    var   <- colnames(x)
    which <- selectvar(which,var)

    np <- length(which)

    dots   <- list(...)
    nmdots <- names(dots)

    ## Set par mfrow and ask.
    ask <- setplotpar(nmdots, dots, np, ask)

    ## interactively wait if there are remaining figures
    if (is.null(ask))
        ask <- prod(par("mfrow")) < length(which) && dev.interactive()
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
        do.call("hist", c(alist(x[, i]), dots))
    }
}
### ============================================================================

image.deSolve <- function (x, which = NULL, ask = NULL,
  add.contour = FALSE, grid=NULL, method="image", ...) {

    dimens <- attributes(x)$dimens
    if (is.null(dimens))
      stop("cannot make an image from deSolve output which is 0-dimensional")
    else if (length(dimens) ==1)  # 1-D
      plot.ode1D(x, which, ask, add.contour, grid, method=method, ...)
    else if (length(dimens) ==2)  # 2-D
      plot.ode2D(x, which, ask, add.contour, grid, method=method, ...)
    else
      stop("cannot make an image from deSolve output with more than 2 dimensions")
}

### ============================================================================
### Plot utilities for the S3 plot method, 0-D, 1-D, 2-D
### ============================================================================
# karline: from version 1.8.2, also possible to plot observations.

plot.deSolve <- function (x, which = NULL, ask = NULL, x2 = NULL, obs = NULL, 
    obspar = list(), ...) {

    getname <- function (x)
      if (is.data.frame(x)) names(x) else colnames(x)

    var   <- colnames(x)
    if (! is.null(obs)) {
       obsname <- getname(obs) 
       if (! class(obs) %in% c("data.frame","matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
    }
    Which <- which 
    if (is.null(Which) & is.null(obs))   # All variables plotted
      Which <- 1:(ncol(x)-1)
      
    else if (is.null(Which)) {           # All common variables in x and obs plotted
     Which <- which(var %in% obsname)
     Which <- Which[Which != 1]          # remove first element (x-value)
     Which <- var[Which]                 # names rather than 
    } 

    # check x2: if single instance of "deSolve": put it in a list
    nother <- 0
    if ("deSolve" %in% class(x2))
      x2 <- list(x2)
    
    if (is.list(x2)) {
      nother <- length(x2)
      for ( i in 1:nother) {
        if (min(colnames(x2[[i]]) == colnames(x)) == 0)
          stop(" 'x2' and 'x' are not compatible - colnames not the same")
        if (min(dim(x2[[i]]) - dim(x) == c(0, 0)) == 0) 
          stop(" 'x2' and 'x' are not compatible - dimensions not the same")
      }
    } 
    
    # Position of variables in "x"
    t  <- 1     # column with "times" 
    xWhich   <- selectvar(Which,var)

    # Position of variables in "obs" (NA = not observed)
    if (! is.null(obs)) {
      ObsWhich <- selectvar(var[xWhich], obsname, NAallowed = TRUE)
      ObsWhich [ ObsWhich > ncol(obs)] <- NA
    } else 
      ObsWhich <- rep(NA, length(xWhich))

    np <- length(xWhich)

    dots   <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and
    # interactively wait if there are remaining figures

    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    xxlab <- if (is.null(dots$xlab))  colnames(x)[t]  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""              else dots$ylab
    Main <-  if (is.null(dots$main))  colnames(x)[xWhich] else dots$main

    ## allow individual xlab and ylab (vectorized) for each figure
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)
    Main  <- rep(Main,  length.out=np)

    isylim <- !is.null(dots$ylim)
    ylim  <- dots$ylim
    
    # dots for other series - change
    if (! is.null(x2)) {
      dots2 <- list()
      for ( i in 1:nother) {
        dd <- list(lty = i+1, col = i+1, bg = i+1)
        if ( length(dots$lty) == nother+1)
          dd$lty <- dots$lty[i+1]
        if ( length(dots$col) == nother+1)
          dd$col <- dots$col[i+1]
        if ( length(dots$bg) == nother+1)
          dd$bg <- dots$bg[i+1]
        dots2[[i]] <- dd 
      }  
    }

    for (i in 1:np) {
        ii <- xWhich[i]
        io <- ObsWhich[i]
        dots$ylim  <- ylim
        if (! isylim) {
          xs <- x[, ii]
          if (! is.null(x2)) 
           for (j in 1:nother) 
             xs <- c(xs,x2[[j]][,ii])
          if (! is.na(io)) xs <- c(xs,obs[,io])
            dots$ylim <- range(xs, na.rm = TRUE)
        }  
        dots$main <- Main[i]
        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]
        do.call("plot", c(alist(x[, t], x[, ii]), dots))
        if (! is.null(x2)) 
          for (j in 1:nother)
          do.call("lines", c(alist(x2[[j]][, t], x2[[j]][, ii]), dots2[[j]]))
        
        if (! is.na(io)) 
          do.call("points", c(alist(obs[, 1], obs[, io]), obspar))        
    }
}


### ============================================================================
select1dvar <- function (which,var) {

    if (!is.numeric(which)) {
        ln <- length(which)
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i]==var)
          if (length(ss) ==0)
            stop("variable ", which[i], " not in var")
          Select <- c(Select,ss)
        }
    }
    else {
        Select <- which  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
 }

### ============================================================================
# to drape a color over a persp plot.
### ============================================================================


drapecol <- function (A, 
          col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100), 
              NAcol = "white")
{
    nr <- nrow(A)
    nc <- ncol(A)
    ncol <- length(col)
    AA <- 0.25 * (A[1:(nr - 1), 1:(nc - 1)] + A[1:(nr - 1), 2:nc] +
        A[2:nr, 1:(nc - 1)] + A[2:nr, 2:nc])
    Ar <- range(AA, na.rm = TRUE)
    rn <- Ar[2] - Ar[1]
    ifelse(rn != 0, drape <- col[1 + trunc((AA - Ar[1])/rn *
        (ncol - 1))], drape <- rep(col[1], ncol))
    drape[is.na(drape)] <- NAcol
    return(drape)
}

### ============================================================================

plot.1D <- function (x, which=NULL, ask=NULL, grid=NULL, xyswap = FALSE, ...) {


    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)
    if (length(dimens) != 1)
      stop ("plot.1D only works for models solved with 'ode.1D'")

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec

    var <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,var)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and
    # interactively wait if there are remaining figures

    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }


    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "x"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  var  else dots$ylab

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)
    times <- x[,1]
    Main <-  if (is.null(dots$main)) paste("time",times) else rep(dots$main, length.out =length(times))

    for (j in 1:length(times)) {
      for (i in which) {
        dots$main <- Main[j]
        istart <- (i-1)*proddim + 1
        out <- x[j,(istart+1):(istart+prod(dimens))]

        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]

        if (! is.null(grid))
          Grid = grid
        else
          Grid = 1:length(out)

        if (! xyswap) {
          dots$xlab <- xxlab[i]
          dots$ylab <- yylab[i]
          do.call("plot", c(alist(Grid, out), dots))
        } else {
          if (is.null(labs)){
            dots$ylab <- xxlab[i]
            dots$xlab <- yylab[i]
          } else {
            dots$xlab <- xxlab[i]
            dots$ylab <- yylab[i]
            dots$ylim <- rev(range(Grid))    # y-axis reversed
          }
          do.call("plot", c(alist(out, Grid), dots))
        }
       }
     }
}

### ============================================================================

plot.ode1D <- function (x, which, ask, add.contour, grid, method="image", ...) {

# if x is vector, check if there are enough columns ...
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec

    var <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,var)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and
    # interactively wait if there are remaining figures

    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    Main <-  if (is.null(dots$main)) var else rep(dots$main, length.out =np)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "times"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""   else dots$ylab
    if (method=="persp")
      dotscol <- dots$col

    else if (method == "filled.contour")
    dots$color.palette <- if (is.null(dots$color.palette))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))else dots$color.palette
    else
    dots$col <- if (is.null(dots$col))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col

    dotslim <- dots$zlim

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)
    times <- x[,1]
    for (i in which) {
        dots$main <- Main[i]
        istart <- (i-1)*proddim + 1
        out <- x[,(istart+1):(istart+prod(dimens))]

        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]

        List <- alist(z=out,x=times)
        if (! is.null(grid)) List$y = grid

        if (method=="persp") {
           if(is.null(dotscol))
             dots$col <- drapecol(out,
               colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100))
           else
              dots$col<-drapecol(out,dotscol)
          if (is.null(dotslim))
            if (diff(range(out, na.rm=TRUE)) == 0) dots$zlim = c(0,1)
          else
            dots$zlim = dotslim

        }
        do.call(method, c(List, dots))
        if (add.contour) do.call("contour", c(List, add=TRUE))
    }
}

### ============================================================================

plot.ode2D <- function (x, which, ask, add.contour, grid, method="image",
   ...) {


# if x is vector, check if there are enough columns ...
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec

    var <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,var)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and
    # interactively wait if there are remaining figures

    Ask <- setplotpar(nmdots, dots, np, ask)

# here ask is always true...
    if (is.null(ask)) ask <- TRUE
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

#    Main <-  if (is.null(dots$main)) var else rep(dots$main, length.out =np)
    N <- np * nrow(x)
    Main <-  if (is.null(dots$main)) rep(var,length.out=N) else rep(dots$main,
      length.out =N)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "x"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  "y"   else dots$ylab

    if (method=="persp")
      dotscol <- dots$col

    else if (method == "filled.contour")
    dots$color.palette <- if (is.null(dots$color.palette))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))else dots$color.palette
    else
    dots$col <- if (is.null(dots$col))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col

    dotslim <- dots$zlim

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = N)
    yylab <- rep(yylab, length.out = N)

    ii <- 0
    for (nt in 1:nrow(x)) {
      for (i in which) {
        ii <- ii+1
        dots$main <- Main[ii]
        istart <- (i-1)*proddim + 1
        out <- x[nt,(istart+1):(istart+prod(dimens))]
        dim(out) <- dimens
        dots$xlab <- xxlab[ii]
        dots$ylab <- yylab[ii]

        List <- alist(z=out)
        if (! is.null(grid)) {
          List$x <- grid$x
          List$y <- grid$y
        }

        if (method=="persp") {
           if(is.null(dotscol))
             dots$col <- drapecol(out,
               colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100))
           else
              dots$col<-drapecol(out,dotscol)
          if (is.null(dotslim))
            if (diff(range(out, na.rm=TRUE)) == 0) dots$zlim = c(0,1)
          else
            dots$zlim = dotslim

        }
        do.call(method, c(List, dots))
        if (add.contour) do.call("contour", c(List, add=TRUE))
     }
     if (sum(par("mfrow") - c(1,1))==0)
       mtext(outer=TRUE, side=3,paste("time ",x[nt,1]), cex=1.5, line=-1.5)

   }
}


