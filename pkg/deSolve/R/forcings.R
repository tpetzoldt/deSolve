### ============================================================================
### Chekc forcing function data set
### ============================================================================


checkforcings <- function (forcings,times,dllname,initforc,verbose) {

 if (is.null(initforc))
   stop(paste("initforc should be loaded if there are forcing functions ",initforc))
 if (is.loaded(initforc, PACKAGE = dllname, type = "") ||
     is.loaded(initforc, PACKAGE = dllname, type = "Fortran")) {
 ModelForc <- getNativeSymbolInfo(initforc, PACKAGE = dllname)$address
 } else
   stop(paste("initforc should be loaded if there are forcing functions ",initforc))

  if (is.data.frame(forcings)) forcings <- list(a=forcings)
  if (! is.list(forcings)) forcings <- list(a=forcings)
  nf <- length(forcings)
  #1 check if each forcing function consists of a 2-columned matrix
  for (i in 1:nf) {
    if (ncol(forcings[[i]]) != 2)
      stop("forcing function data sets should consist of two-colum matrix")
  }
  #2 time span of forcing function data sets should embrace simulation time...
  r_t <- range(times)
  for (i in 1:nf) {
    r_f <- range(forcings[[i]][,1])
    if (r_f[1] > r_t[1]) {
      mint <- c(r_t[1],forcings[[i]][1,2] )
      forcings[[i]] <- rbind(mint,forcings[[i]])
     if(verbose)
       warning(paste("extrapolating forcing function data sets to first timepoint",i))
    }
    nr   <- nrow(forcings[[i]])
    if (r_f[2] < r_t[2]) {
      maxt <- c(r_t[2],forcings[[i]][nr,2] )
      forcings[[i]] <- rbind(forcings[[i]],maxt)
     if(verbose)
       warning(paste("extrapolating forcing function data sets to last timepoint",i))
    }
#    maxt <- c(forcings[[i]][nr,1]*1e6,forcings[[i]][nr,2] )    # arbitrary
#    forcings[[i]] <- rbind(forcings[[i]],maxt)
  }
  #3 all forcings in one matrix; index to start/end
  fmat <- tmat <- NULL
  imat <- rep(1,nf+1)
  for (i in 1:nf) {
    tmat <- c(tmat, forcings[[i]][,1])
    fmat <- c(fmat, forcings[[i]][,2])
    imat[i+1]<-imat[i]+nrow(forcings[[i]])
  }
  storage.mode(tmat) <- storage.mode(fmat) <- "double"
  storage.mode(imat) <- "integer"

  return(list(tmat=tmat,fmat=fmat,imat=imat,ModelForc=ModelForc))
}
