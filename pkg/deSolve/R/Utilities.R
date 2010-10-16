### ============================================================================
### ============================================================================
### S3 methods
### karline+Thomas: from version 1.9, also possible to plot multiple
###                 outputs and to add observations.
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

# ks->Th: for xlim and ylim....
expanddotslist <- function (dots, n) {
  if (is.null(dots)) return(dots)
  dd <- if (!is.list(dots )) list(dots) else dots
  rep(dd, length.out = n)
}

## =============================================================================
## functions for expanding arguments in dots  (...)
## =============================================================================

repdots <- function(dots, n) 
  if (is.function(dots)) dots else rep(dots, length.out = n)

setdots <- function(dots, n) lapply(dots, repdots, n)

## =============================================================================
## function for extracting element 'index' from dots  (...)
## =============================================================================

extractdots <- function(dots, index) {
  ret <- lapply(dots, "[", index)
  ret <- lapply(ret, unlist) ## thpe: flatten list (experimental)
  return(ret)
}

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
   
    t        <- 1     # column with "times"
    varnames <- colnames(x)
    which    <- selectvar(which, varnames)

    np     <- length(which)
    dots   <- list(...)
    nmdots <- names(dots)

## Set par mfrow and ask.
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

## expand all dots to np values (no defaults        
    Dots  <- setdots(dots, np)   

    # different from default settings
    Dots$main <- expanddots (dots$main, varnames[which], np)
    Dots$xlab <- expanddots (dots$xlab, varnames[t],     np)
    Dots$xlab <- expanddots (dots$xlabb,  ""        ,     np)

## xlim and ylim are special: they are vectors or lists
    xxlim <- expanddotslist(dots$xlim, np)
    yylim <- expanddotslist(dots$ylim, np)

##    
    for (ii in 1:np) {
        i <- which[ii]
        dots <- extractdots(Dots, ii)
        if (! is.null(xxlim)) dots$xlim <- xxlim[[ii]]
        if (! is.null(yylim)) dots$ylim <- yylim[[ii]]
        do.call("hist", c(alist(x[, i]), dots))
    }
}
### ============================================================================
### Image, filled.contour and persp plots
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

plot.deSolve <- function (x, ..., which = NULL, ask = NULL, obs = NULL, 
    obspar = list()) {

## check observed data
    if (! is.null(obs)) {
       obsname <- colnames(obs) 
       if (! class(obs) %in% c("data.frame", "matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
    }

## variables to be plotted
    varnames <- colnames(x)
    Which   <- which 
    
    if (is.null(Which) & is.null(obs))  # All variables plotted
      Which <- 1 : (ncol(x)-1)
      
    else if (is.null(Which)) {          # All common variables in x and obs plotted
     Which <- which(varnames %in% obsname)
     Which <- Which[Which != 1]         # remove first element (x-value)
     Which <- varnames[Which]           # names rather than numbers
    } 

## Position of variables to be plotted in "x" 
    t       <- 1     # column with "times" 
    xWhich  <- selectvar(Which, varnames)
    np      <- length(xWhich)

    ldots   <- list(...)
    ndots   <- names(ldots)

## number of figures in a row and interactively wait if remaining figures
    ask <- setplotpar(ndots, ldots, np, ask)
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

## create two lists: x2:   other deSolve objects, 
##                   dots: remaining (plotting) parameters
    x2     <- list()
    dots   <- list()
    nd     <- 0
    nother <- 0            

    if (length(ldots) > 0) 
     for ( i in 1:length(ldots))
      if ("deSolve" %in% class(ldots[[i]])) {
        x2[[nother <- nother + 1]] <- ldots[[i]]  
        names(x2)[nother] <- ndots[i]
      } else if (! is.null(ldots[[i]])) {
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }
    nmdots <- names(dots)

## check compatibility of all deSolve objects    
    if (nother > 0) {
      for ( i in 1:nother) {            
        if (min(colnames(x2[[i]]) == varnames) == 0)
          stop("'x' is not compatible with other deSolve objects - colnames not the same")
      }
    } 

    nx <- nother + 1 # total number of deSolve objects to be plotted

## Position of variables in "obs" (NA = not observed)
    nobs <- 0
    if (! is.null(obs)) {
      ObsWhich <- selectvar(varnames[xWhich], obsname, NAallowed = TRUE)
      ObsWhich [ ObsWhich > ncol(obs)] <- NA  # Ks->ks check why this was necessary...
      nobs <- length(ObsWhich)
      ObsDots <- setdots(obspar, nobs)
    } else 
      ObsWhich <- rep(NA, np)

## plotting parameters : split in plot parameters and point parameters
    plotnames <- c("xlab","ylab","xlim","ylim","main","sub","log","asp",
                   "ann","axes","frame.plot","panel.first","panel.last",
                   "cex.lab","cex.axis","cex.main")    
    
    # plot.default parameters
    ii <- names(dots) %in% plotnames
    dotmain <- dots[ii]
    dotmain <- setdots(dotmain, np)  # expand to np for each plot

    # these are different from the default
    dotmain$xlab <- expanddots(dots$xlab, varnames[t]     , np)
    dotmain$ylab <- expanddots(dots$ylab, ""              , np)
    dotmain$main <- expanddots(dots$main, varnames[xWhich], np)

    # ylim and xlim can be lists and are at least two values
    isylim <- !is.null(dots$ylim)
    yylim  <- expanddotslist(dots$ylim, np)

    isxlim <- !is.null(dots$xlim)
    xxlim  <- expanddotslist(dots$xlim, np)

    # point parameters
    ip <- !names(dots) %in% plotnames
    dotpoints <- dots[ip]
    dotpoints <- setdots(dotpoints, nx)   # expand all dots to nx values

    # these are different from default
    dotpoints$type <- expanddots(dots$type, "l", nx)
    dotpoints$lty  <- expanddots(dots$lty, 1:nx, nx)
    dotpoints$pch  <- expanddots(dots$pch, 1:nx, nx)
    dotpoints$col  <- expanddots(dots$col, 1:nx, nx)
    dotpoints$bg   <- expanddots(dots$bg,  1:nx, nx)
 
## for each output variable (plot)
    iobs <- 0
    for (i in 1 : np) {
      ii <- xWhich[i]     # position of variable in 'x'
      io <- ObsWhich[i]   # position of variable in 'obs'
      
      # plotting parameters for deSolve output 1 (opens a plot)
      Dotmain   <- extractdots(dotmain, i)
      Dotpoints <- extractdots(dotpoints, 1)
      
      # ranges
      if (! isylim) {
        yrange <- range(x[, ii])
        if (nother>0) 
         for (j in 1:nother) 
           yrange <- range(yrange, x2[[j]][,ii], na.rm = TRUE)
        if (! is.na(io)) yrange <- range(yrange, obs[,io], na.rm = TRUE)
          Dotmain$ylim <- yrange
      } else {
        Dotmain$ylim  <- yylim[[i]]
      } 

      if (! isxlim) {
        xrange <- range(x[, t])
        if (nother>0) 
         for (j in 1:nother) 
           xrange <- range(xrange, x2[[j]][,t], na.rm = TRUE)
        if (! is.na(io)) xrange <- range(xrange, obs[,1], na.rm = TRUE)
          Dotmain$xlim <- xrange
      } else {
        Dotmain$xlim  <- xxlim[[i]]
      } 
      
      # first deSolve object plotted (new plot created)
      do.call("plot", c(alist(x[, t], x[, ii]), Dotmain, Dotpoints))
      
      # if other deSolve outputs
      if (nother > 0) 
        for (j in 2:nx) {
          Dotpoints <- extractdots(dotpoints, j)
          do.call("lines", c(alist(x2[[j-1]][, t], x2[[j-1]][, ii]), Dotpoints))
        }
      # ks -> ThPe if observed variables: select correct pars
      if (! is.na(io)) {
        iobs <- iobs + 1
        do.call("points", c(alist(obs[, 1], obs[, io]), 
                extractdots(ObsDots, iobs)))     
      }     
    }
}


### ============================================================================
select1dvar <- function (which, var) {

    if (!is.numeric(which)) {
        ln <- length(which)
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i] == var)
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

plot.1D <- function (x, which = NULL, ask = NULL, grid = NULL, xyswap = FALSE, ...) {

## Check settings of x
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)
    if (length(dimens) != 1)
      stop ("plot.1D only works for models solved with 'ode.1D'")

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec

    varnames <- if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which, varnames)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

## number of figures in a row and interactively wait if remaining figures
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    Dots <- setdots(dots, np) # expand all dots to np values (no defaults)

    # These are different from defaulst
    llab <- (!is.null(dots$xlab) |!is.null(dots$ylab)) 
    Dots$xlab <- expanddots(dots$xlab,  "x", np)
    Dots$ylab <- expanddots(dots$ylab,  varnames[which], np)

    if (xyswap & !llab) {
        xl <- Dots$ylab
        Dots$ylab <- Dots$xlab
        Dots$xlab <- Dots$ylab
    }

    # allow individual xlab and ylab (vectorized)
    times <- x[,1]
    Dots$main <- expanddots(Dots$main, paste("time", times), np)

    # xlim and ylim are special: 
    xxlim <- expanddotslist(dots$xlim, np)
    yylim <- expanddotslist(dots$ylim, np)

    ## thpe: additional check *before* the loop
    if (! is.null(grid))
      Grid <- grid
    else
      Grid <- 1:length(out)
    
    for (j in 1:length(times)) {
      for (i in which) {
        dots <- extractdots(Dots, i)
        if (! is.null(xxlim)) dots$xlim <- xxlim[[i]]
        if (! is.null(yylim)) dots$ylim <- yylim[[i]]
        if (is.null(yylim) & xyswap) 
           dots$ylim <- rev(range(Grid))    # y-axis 
        istart <- (i - 1) * proddim + 1
        out <- x[j, (istart+ 1 ) : (istart + prod(dimens))]

        if (! is.null(grid))
          Grid <- grid
        else
          Grid <- 1:length(out)

        if (! xyswap) {
          do.call("plot", c(alist(Grid, out), dots))
        } else {
          do.call("plot", c(alist(out, Grid), dots))
        }
       }
     }
}

### ============================================================================

plot.ode1D <- function (x, which, ask, add.contour, grid, method = "image", ...) {

## if x is vector, check if there are enough columns ...
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec * proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

## variables to be plotted
    if (is.null(which))
      which <- 1:nspec

    varnames <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which, varnames)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

## number of figures in a row and interactively wait if remaining figures
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    Dots  <- setdots(dots, np)   # expand dots to np values (no defaults)

    # different from the default
    Dots$main  <- expanddots(dots$main, varnames[which], np)
    Dots$xlab  <- expanddots(dots$xlab, "times",  np)
    Dots$ylab  <- expanddots(dots$ylab, "",       np)

    # colors - different if persp, image or filled.contour
    dotscol <- NULL
    if (method == "persp")
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

## xlim and ylim are special: - is there a better way????
    xxlim <- expanddotslist(dots$xlim, np)
    yylim <- expanddotslist(dots$ylim, np)

    times <- x[,1]

## for each output variable (plot)
    for (i in 1:np) {
      ii <- which[i]
      dots      <- extractdots(Dots, i)
      dots$col  <- dotscol 
      if (! is.null(xxlim)) dots$xlim <- xxlim[[i]]
      if (! is.null(yylim)) dots$ylim <- yylim[[i]]

      istart <- (ii-1)*proddim + 1
      out    <- x[ ,(istart+1):(istart+prod(dimens))]

      List <- alist(z = out, x = times)
      if (! is.null(grid)) List$y = grid

        if (method == "persp") {
           if(is.null(dotscol))
             dots$col <- drapecol(out,
               colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100))
           else
              dots$col <- drapecol(out, dotscol)
          if (is.null(dotslim))  # this to prevent error when range = 0
            if (diff(range(out, na.rm=TRUE)) == 0) 
              dots$zlim <- c(0, 1)
          else
            dots$zlim <- dotslim

        } else if (method == "filled.contour")
          dots$color.palette <- dotscolorpalette
        
        do.call(method, c(List, dots))
        if (add.contour) do.call("contour", c(List, add = TRUE))
    }
}

### ============================================================================

plot.ode2D <- function (x, which, ask, add.contour, grid, method = "image",
   ...) {

## if x is vector, check if there are enough columns ...
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

## variables to be plotted
    if (is.null(which))
      which <- 1:nspec

    varnames <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,varnames)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

## number of figures in a row and interactively wait if remaining figures
    Ask <- setplotpar(nmdots, dots, np, ask)

# here ask is always true by default...
    if (is.null(ask)) ask <- TRUE
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

#    Main <-  if (is.null(dots$main)) varnames else rep(dots$main, length.out =np)
    N <- np * nrow(x)

    Dots  <- setdots(dots, N)  # expand dots to np values (no defaults)
    # different from the default
    Dots$main <- expanddots(dots$main, varnames[which], N)
    Dots$xlab <- expanddots(dots$xlab,  "x"  , N)
    Dots$ylab <- expanddots(dots$ylab,  "y"  , N)

    if (method == "persp")
      dotscol <- dots$col

    else if (method == "filled.contour")
    dotscolor.palette <- if (is.null(dots$color.palette))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))else dots$color.palette
    else
    dotscol <- if (is.null(dots$col))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dotscol      

    dotslim <- dots$zlim
 
    i <- 0
    for (nt in 1:nrow(x)) {
      for (ii in which) {
        i       <- i+1
        dots    <- extractdots(Dots, i)
        istart  <- (ii-1)*proddim + 1
        out <- x[nt,(istart+1):(istart+prod(dimens))]
        dim(out) <- dimens

        List <- alist(z = out)
        if (! is.null(grid)) {
          List$x <- grid$x
          List$y <- grid$y
        }

        if (method == "persp") {
           if(is.null(dotscol))
             dots$col <- drapecol(out,
               colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100))
           else
              dots$col <- drapecol(out, dotscol)
          if (is.null(dotslim))
            if (diff(range(out, na.rm = TRUE)) == 0) 
              dots$zlim <- c(0, 1)
          else
            dots$zlim = dotslim

        } else if (method == "image")
          dots$col <- dotscol
        else if (method == "filled.contour")
          dots$color.palette <- dotscolor.palette        
        
        do.call(method, c(List, dots))
        if (add.contour) do.call("contour", c(List, add = TRUE))
     }
     if (sum(par("mfrow") - c(1, 1)) == 0)
       mtext(outer = TRUE, side = 3, paste("time ", x[nt, 1]), 
         cex = 1.5, line = -1.5)

   }
}