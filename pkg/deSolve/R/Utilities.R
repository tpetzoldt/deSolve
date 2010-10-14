### ============================================================================
### ============================================================================
### S3 methods
### ============================================================================
### ============================================================================


### ============================================================================
### first some common functions
### ============================================================================

## =============================================================================
## function for checking and expanding arguments in dots (...) with default
## =============================================================================

expanddots <- function (dots, default, n) {
  dots <- if (is.null(dots)) default else dots
  rep(dots, length.out = n)
}

## =============================================================================
## functions for expanding arguments in dots  (...)
## =============================================================================

repdots <- function(dots, n) 
  if (is.function(dots)) dots else rep(dots, length.out = n)

setdots <- function(dots, n) lapply(dots, repdots, n)

## =============================================================================
## function for extracting element ii from dots  (...)
## =============================================================================

extractdots <- function(dots, index)   lapply(dots, "[", index)

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)), 3)
      nr <- min(ceiling(nv/nc), 3)
      mfrow <- c(nr, nc)
    } else if ("mfcol" %in% nmdots)
      mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow))  mf <- par(mfrow = mfrow)

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
  } else {
    Select <- which + 1  # "Select now refers to the column number
    if (max(Select) > length(var))
        stop("index in 'which' too large")
    if (min(Select) < 1)
        stop("index in 'which' should be > 0")
  }
  return(Select)
}

### ============================================================================
### print a deSolve object
### ============================================================================

print.deSolve <- function(x, ...)
  print(as.data.frame(x), ...)

### ============================================================================
### Create a histogram for a list of variables
### ============================================================================

hist.deSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {
   
    t     <- 1     # column with "times"
    varnames <- colnames(x)
    which    <- selectvar(which, varnames)

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
    
    Main  <- expanddots (dots$main, varnames[which], np)
    xxlab <- expanddots (dots$xlab, varnames[t],     np)
    yylab <- expanddots (dots$ylab,  ""        ,     np)
    
    Dots  <- setdots(dots, np)   # expand all dots to np values (no defaults)
    
    for (ii in 1:np) {
        i <- which[ii]
        dots <- extractdots(Dots, ii)
        dots$main <- Main[ii]
        dots$xlab <- xxlab[ii]
        dots$ylab <- yylab[ii]
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
# karline: from version 1.9, also possible to plot observations.

plot.deSolve <- function (x, ..., which = NULL, ask = NULL, obs = NULL, 
    obspar = list()) {

## check obs
    if (! is.null(obs)) {
       obsname <- colnames(obs) 
       if (! class(obs) %in% c("data.frame", "matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
    }

## variables to be plotted
    varnames <- colnames(x)
    Which   <- which 
    
    if (is.null(Which) & is.null(obs))   # All variables plotted
      Which <- 1 : (ncol(x)-1)
      
    else if (is.null(Which)) {           # All common variables in x and obs plotted
     Which <- which(varnames %in% obsname)
     Which <- Which[Which != 1]          # remove first element (x-value)
     Which <- varnames[Which]             # names rather than numbers
    } 

## Position of variables to be plotted in "x" and "x2"
    t  <- 1     # column with "times" 
    xWhich   <- selectvar(Which,varnames)
    np <- length(xWhich)

    ldots   <- list(...)
    ndots <- names(ldots)

    # number of figures in a row and
    # interactively wait if there are remaining figures

    ask <- setplotpar(ndots, ldots, np, ask)
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

## create two lists: x2: other deSolve objects, dots: plotting parameters
    x2 <- list()
    dots <- list()
    nd <- 0
    nother <- 0                          # number of other steady1D instances

    if (length(ldots) > 0) 
     for ( i in 1:length(ldots))
      if ("deSolve" %in% class(ldots[[i]])){
        x2[[nother <- nother+1]] <- ldots[[i]]  
        names(x2)[nother] <- ndots[i]
      } else if (! is.null(ldots[[i]])){
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }
    nmdots <- names(dots)
    
    if (nother > 0) {
     for ( i in 1:nother) {             # check compatibility of x2/x inputs
       if (! "deSolve" %in% class(x2[[i]]))
         stop ("elements in 'x2' should be of class 'deSolve'")
        
        # ks->thpe "times" can now be varied
        if (min(colnames(x2[[i]]) == varnames) == 0)
          stop(" 'x2' and 'x' are not compatible - colnames not the same")
      }
    } 
    
    nx <- nother + 1
    

## Position of variables in "obs" (NA = not observed)
    if (! is.null(obs)) {
      ObsWhich <- selectvar(varnames[xWhich], obsname, NAallowed = TRUE)
      ObsWhich [ ObsWhich > ncol(obs)] <- NA    # Ks->ks check why this was necessary...
    } else 
      ObsWhich <- rep(NA, length(xWhich))

## plotting parameters :
## for each plot (1:np): xlim, ylim, log, main, sub, xlab, ylab, ann,
##                      what        default            total #
    xxlab <- expanddots(dots$xlab,  varnames[t]      , np)
    yylab <- expanddots(dots$ylab,  ""               , np)
    Main  <- expanddots(dots$main,  varnames[xWhich] , np)
    Sub   <- expanddots(dots$sub ,  ""               , np)
    Log   <- expanddots(dots$log ,  ""               , np)

    # ylim and xlim can be lists
    isylim <- !is.null(dots$ylim)
    yylim   <- dots$ylim

    isxlim <- !is.null(dots$xlim)
    xxlim   <- dots$xlim

    if(!is.list(yylim)) yylim <- list(yylim)
    yylim <- rep(yylim, length.out = np)

    if(!is.list(xxlim)) xxlim <- list(xxlim)
    xxlim <- rep(xxlim, length.out = np)
    
    if (nx > 1) dotsx2 <- list()

## for each deSolve output (nx) within a plot: col, lty, lwd, type
    
    Type <- expanddots(dots$type, "l",      nx)
    Lty  <- expanddots(dots$lty, 1:nx,      nx)
    Lwd  <- expanddots(dots$lwd, par("lwd"),nx)
    Pch  <- expanddots(dots$pch, 1:nx,      nx)
    Cex  <- expanddots(dots$cex, 1:nx,      nx)
    Col  <- expanddots(dots$col, 1:nx,      nx)
    Bg   <- expanddots(dots$bg,  1:nx,      nx)
 
## ks->Thpe: also here? 
    Dots  <- setdots(dots, nx)   # expand all dots to np values (no defaults)

## for each output variable (plot)
    for (i in 1 : np) {
      ii <- xWhich[i]     #position in 'x'
      io <- ObsWhich[i]   #position in 'obs'
      
      ## plotting parameters for deSolve output 1
      dots      <- extractdots(Dots, i)
      dots$main <- Main[i]
      dots$sub  <- Sub[i]
      dots$log  <- Log[i]
      dots$xlab <- xxlab[i]
      dots$ylab <- yylab[i]

      if (! isylim) {
        yrange <- range(x[, ii])
        if (nother>0) 
         for (j in 1:nother) 
           yrange <- range(yrange, x2[[j]][,ii], na.rm = TRUE)
        if (! is.na(io)) yrange <- range(yrange, obs[,io], na.rm = TRUE)
          dots$ylim <- yrange
      } else {
        dots$ylim  <- yylim[[i]]
      } 

      if (! isxlim) {
        xrange <- range(x[, t])
        if (nother>0) 
         for (j in 1:nother) 
           xrange <- range(xrange, x2[[j]][,t], na.rm = TRUE)
        if (! is.na(io)) xrange <- range(xrange, obs[,1], na.rm = TRUE)
          dots$xlim <- xrange
      } else {
        dots$xlim  <- xxlim[[i]]
      } 

      dots$lty  <- Lty[1]
      dots$lwd  <- Lwd[1]
      dots$pch  <- Pch[1]
      dots$cex  <- Cex[1]
      dots$col  <- Col[1]
      dots$bg   <- Bg[1]
      dots$type <- Type[1]
      
      do.call("plot", c(alist(x[, t], x[, ii]), dots))
      
      # if other deSolve outputs
      if (nother>0) 
        for (j in 2:nx) {
          dotsx2$lty  <- Lty[j]
          dotsx2$lwd  <- Lwd[j]
          dotsx2$pch  <- Pch[j]
          dotsx2$cex  <- Cex[j]
          dotsx2$col  <- Col[j]
          dotsx2$bg   <- Bg[j]
          dotsx2$type <- Type[j]
        
          do.call("lines", c(alist(x2[[j-1]][, t], x2[[j-1]][, ii]), dotsx2))
        }
      # if observed variables
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

    varnames <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,varnames)

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
    xxlab <- expanddots(dots$xlab,  "x", np)
    yylab <- expanddots(dots$ylab,  varnames, np)

    ## allow individual xlab and ylab (vectorized)
    times <- x[,1]
    Main <-  if (is.null(dots$main)) paste("time",times) else rep(dots$main, length.out =length(times))

   ## ???
    Dots  <- setdots(dots, np)   # expand all dots to np values (no defaults)

    for (j in 1:length(times)) {
      for (i in which) {
        dots      <- extractdots(Dots, i)
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

plot.ode1D <- function (x, which, ask, add.contour, grid, method = "image", ...) {

# if x is vector, check if there are enough columns ...
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec

    varnames <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,varnames)

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

    Main  <-  expanddots(dots$main, varnames, np)
    xxlab <- expanddots(dots$xlab, "times", np)
    yylab <- expanddots(dots$ylab,"",       np)
    dotscol <- NULL
    if (method=="persp")
      dotscol <- dots$col

    else if (method == "filled.contour")
    dotscolorpalette <- if (is.null(dots$color.palette))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))else dots$color.palette
    else
    dotscol <- if (is.null(dots$col))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col

    dotslim <- dots$zlim
    Dots  <- setdots(dots, np)   # expand all dots to np values (no defaults)

    times <- x[,1]
    for (ii in 1:np) {
        i <- which[ii]
        dots      <- extractdots(Dots, ii)
        dots$main <- Main[ii]
        istart <- (i-1)*proddim + 1
        out <- x[,(istart+1):(istart+prod(dimens))]

        dots$xlab <- xxlab[ii]
        dots$ylab <- yylab[ii]

        dots$col <- dotscol 
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

        } else if (method == "filled.contour")
        dots$color.palette <- dotscolorpalette
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

    varnames <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,varnames)

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

#    Main <-  if (is.null(dots$main)) varnames else rep(dots$main, length.out =np)
    N <- np * nrow(x)
    labs <- (is.null(dots$xlab) && is.null(dots$ylab))

    Main <-  expanddots(dots$main, varnames, N)
    xxlab <- expanddots(dots$xlab,  "x"  , N)
    yylab <- expanddots(dots$ylab,  "y"  , N)

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

    Dots  <- setdots(dots, N)   # expand all dots to np values (no defaults)

    ii <- 0
    for (nt in 1:nrow(x)) {
      for (i in which) {
        ii <- ii+1
        dots      <- extractdots(Dots, ii)
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