diagnostics <- function(obj) {
  Attr <- attributes(obj)
  type <- Attr$type
  if (is.null(type))
    stop("cannot print ODE characteristics; output not of correct type")
  istate <- Attr$istate
  rstate <- Attr$rstate
    cat("--------------------\n")
  if (type == "vode")
    cat("vode return code")
  else if (type == "lsoda")
    cat("lsoda return code")
  else   if (type == "lsodes")
    cat("lsodes return code")
  else if (type == "lsode")
    cat("lsode return code")
  else if (type=="lsodar")
    cat("lsodar return code")
  else if (type == "daspk")
    cat("daspk return code")
  else if (type == "rk")
    cat("rk return code")
    
    cat("\n--------------------\n")

  idid <- istate[1]
  cat(paste("idid = ",idid), "\n")

  if (type %in% c("vode", "lsoda", "lsodes", "lsode", "lsodar")) {
    if (type == "lsodar" && idid ==2 )
      cat(" integration was successful but no root was found\n") else
    if (idid == 2)  cat("integration was successful\n") else
    if (idid == 3) print(" integration was successful and a root was found before reaching the end\n") else
    if (idid == -1) cat(" excess work done on this call. (Perhaps wrong jacobian type MF.)\n") else
    if (idid == -2) cat(" excess accuracy requested. (Tolerances too small.)\n") else
    if (idid == -3) cat(" illegal input detected. (See printed message.)\n") else
    if (idid == -4) cat(" repeated error test failures. (Check all input.)\n") else
    if (idid == -5) cat(" repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.)\n") else
    if (idid == -6) cat(" error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)\n") else
    if (type == lsodes && idid == -7)
       cat(" a fatal error came from sparse solver CDRV by way of DPRJS or DSOLSS\n") else
    if (idid == -7) cat(" work space insufficient to finish (see messages)\n")
  } else if (type == "daspk\n") {
    if (idid >0)   {
      cat (" *** TASK COMPLETED *** \n")
      if (idid == 1) cat("a step was successfully taken in the intermediate-output mode.  The code has not yet reached TOUT\n")
      if (idid == 2) cat("the integration to TSTOP was successfully completed (T = TSTOP) by stepping exactly to TSTOP\n")
      if (idid == 3) cat("the integration to TOUT was successfullycompleted (T = TOUT) by stepping past TOUT. Y(*) and YPRIME(*) are obtained by interpolation\n")
      if (idid == 4) cat("the initial condition calculation, with INFO(11) > 0, was successful, and INFO(14) = 1. No integration steps were taken, and the solution is not considered to have been started\n")
    } else if (idid < 0 & idid > -33)  {
      cat (" *** TASK INTERRUPTED *** \n")
      if (idid == -1) cat("a large amount of work has been expended (about 500 steps)\n") else
      if (idid == -2) cat("the error tolerances are too stringent\n") else
      if (idid == -3) cat("the local error test cannot be satisfied because a zero component in ATOL was specified and the corresponding computed solution component is zero.  Thus, a pure relative error test is impossible for this component\n") else
      if (idid == -5) cat("there were repeated failures in the evaluation or processing of the preconditioner (in jacfunc)\n") else
      if (idid == -6) cat("DDASPK had repeated error test failures on the last attempted step\n") else
      if (idid == -7) cat("the nonlinear system solver in the time integration could not converge\n") else
      if (idid == -8) cat("the matrix of partial derivatives appears to be singular (direct method)\n") else
      if (idid == -9) cat("the nonlinear system solver in the time integration failed to achieve convergence, and there were repeated error test failures in this step\n") else
      if (idid == -10) cat("the nonlinear system solver in the time integration failed to achieve convergence because IRES was equal to -1\n") else
      if (idid == -11) cat("IRES = -2 was encountered and control is being returned to the calling program\n") else
      if (idid == -12) cat("DDASPK failed to compute the initial Y, YPRIME\n") else
      if (idid == -13) cat("unrecoverable error encountered inside user's PSOL routine, and control is being returned to the calling program\n") else
      if (idid == -14) cat("the Krylov linear system solver could not achieve convergence\n")
    } else if (idid ==-33)  {
      cat (" *** TASK TERMINATED *** \n")
      cat("the code has encountered trouble from which it cannot recover.  A message is printed explaining the trouble and control is returned to the calling program.\n")
    }
  } else if (type == "rk") {
      if (idid == 0)  cat("integration was successful\n") else
      if (idid == -1) cat(" a large amount of work has been expended. (Increase maxsteps.)\n") else
      if (idid == -2) cat(" excess accuracy requested. (Tolerances too small.)\n") else
      cat("rk returned with unknown return code\n")
  } else {
    warning("Unknown return type!")
  }

#### istate

  cat("\n--------------------\n")
  cat("ISTATE values\n")
  cat("--------------------\n")
  df <- c( " istate, the return code",
           " The number of steps taken for the problem so far.",
           " The number of function evaluations for the problem so far.",
           " The number of Jacobian evaluations so far.",
           " The method order last used (successfully).",
           " The order to be attempted on the next step.",
           " if istate=-4,-5: the largest comp in error vector",
           " The length of rwork actually required.",
           " The length of iwork actually required.",
           " The number of matrix LU decompositions so far.",
           " The number of nonlinear (Newton) iterations so far.",
           " The number of convergence failures of the solver so far ",
           " The number of error test failures of the integrator so far.")
  if (type == "vode")  {
    ii <- c(1,12:23)
  } else if (type %in% c("lsoda","lsodar")) {
    df[4] <- " The number of Jacobian evaluations and LU decompositions so far."
    df[10]<- " The method indicator for the last succesful step, 1=adams (nonstiff), 2= bdf (stiff)"
    df[11]<- " The current method indicator to be attempted on th next step, 1=adams (nonstiff), 2= bdf (stiff)"
    ii <- c(1,12:21)

  } else if (type == "lsodes") {
    df[4] <- " The number of Jacobian evaluations and LU decompositions so far."
    df[10]<- "  The number of nonzero elements in the sparse jacobian"
    ii <- c(1,12:20)

  } else if (type == "lsode") {
    df[4] <- " The number of Jacobian evaluations and LU decompositions so far."
    ii <- c(1,12:19)

  } else if (type == "daspk") {
    df <- c( " idid, the return code",
             " The order of the method to be attempted on th next step",
             " The order of the method used on the last step",
             " The number of steps taken for the problem so far.",
             " The number of res evaluations for the problem so far.",
             " The number of jacobian evaluations for the problem so far.",
             " The total number of error test failures so far.",
             " The total number of nonlinear convergence failures so far.",
             " The number of convergence failures of the linear iteration so far.",
             " The length of iwork actually required.",
             " The length of rwork actually required.",
             " The number of nonlinear iterations so far.",
             " The number of linear (Krylov) iterations so far ",
             " The number of psol calls so far.")
    ii <- c(1,8:9,12:22)
  } else if (type == "rk") {
    df <- c(" The return code.",
            " The number of steps taken for the problem so far.",
            " The number of function evaluations for the problem so far.",
            " The order of the method.")
    ii <- c(1, 12, 13, 15)
  }

### rstate
  print(data.frame(mess=df[1:length(ii)], val=istate[ii]))

  if (type != "rk") {
    cat("\n--------------------\n")
    cat("RSTATE values\n")
    cat("--------------------\n")
    ii <- 1:4
    df <- c( " The step size in t last used (successfully).",
      " The step size to be attempted on the next step.",
      " The current value of the independent variable which the solver has actually reached",
      " Tolerance scale factor > 1.0 computed when a request for too much accuracy was detected")

    if (type%in% c("lsoda","lsodar")) {
      df <- c(df," the value of t at the time of the last method switch, if any.")
      ii <- 1:5
    }
    print(data.frame(mess=df, val=rstate[ii]))
  }

}
