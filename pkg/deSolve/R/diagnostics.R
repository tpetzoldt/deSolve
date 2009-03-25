printODE <- function(Attr) {
  type <- Attr$type
  if (is.null(type))
    stop("cannot print ODE characteristics; output not of correct type")
  istate <- Attr$istate
  rstate <- Attr$rstate
    print("--------------------")
  if (type == "vode")
    print("vode return code")
  else if (type == "lsoda")
    print("lsoda return code")
  else   if (type == "lsodes")
    print("lsodes return code")
  else if (type == "lsode")
    print("lsode return code")
  else if (type=="lsodar")
    print("lsodar return code")
  else if (type == "daspk")
    print("daspk return code")
    print("--------------------")

  idid <- istate[1]
  print(paste("idid = ",idid))

  if (type != "daspk") {
    if (type == "lsodar" && idid ==2 )
      print(" integration was successful but no root was found") else
    if (idid == 2)  print("integration was successful") else
    if (idid == 3) print(" integration was successful and a root was found before reaching the end") else
    if (idid == -1) print(" excess work done on this call. (Perhaps wrong jacobian type MF.)") else
    if (idid == -2) print(" excess accuracy requested. (Tolerances too small.)") else
    if (idid == -3) print(" illegal input detected. (See printed message.)") else
    if (idid == -4) print(" repeated error test failures. (Check all input.)") else
    if (idid == -5) print(" repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.)") else
    if (idid == -6) print(" error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)") else
    if (type == lsodes && idid == -7)
       print(" a fatal error came from sparse solver CDRV by way of DPRJS or DSOLSS") else
    if (idid == -7) print(" work space insufficient to finish (see messages)")
  } else {
    if (idid >0)   {
      print (" *** TASK COMPLETED *** ")
      if (idid == 1) print("a step was successfully taken in the intermediate-output mode.  The code has not yet reached TOUT")
      if (idid == 2) print("the integration to TSTOP was successfully completed (T = TSTOP) by stepping exactly to TSTOP")
      if (idid == 3) print("the integration to TOUT was successfullycompleted (T = TOUT) by stepping past TOUT. Y(*) and YPRIME(*) are obtained by interpolation")
      if (idid == 4) print("the initial condition calculation, with INFO(11) > 0, was successful, and INFO(14) = 1. No integration steps were taken, and the solution is not considered to have been started")
    } else if (idid<0& idid>-33)  {
      print (" *** TASK INTERRUPTED *** ")
      if (idid == -1) print("a large amount of work has been expended (about 500 steps)") else
      if (idid == -2) print("the error tolerances are too stringent") else
      if (idid == -3) print("the local error test cannot be satisfied because a zero component in ATOL was specified and the corresponding computed solution component is zero.  Thus, a pure relative error test is impossible for this component") else
      if (idid == -5) print("there were repeated failures in the evaluation or processing of the preconditioner (in jacfunc)") else
      if (idid == -6) print("DDASPK had repeated error test failures on the last attempted step") else
      if (idid == -7) print("the nonlinear system solver in the time integration could not converge") else
      if (idid == -8) print("the matrix of partial derivatives appears to be singular (direct method)") else
      if (idid == -9) print("the nonlinear system solver in the time integration failed to achieve convergence, and there were repeated error test failures in this step") else
      if (idid == -10) print("the nonlinear system solver in the time integration failed to achieve convergence because IRES was equal to -1") else
      if (idid == -11) print("IRES = -2 was encountered and control is being returned to the calling program") else
      if (idid == -12) print("DDASPK failed to compute the initial Y, YPRIME") else
      if (idid == -13) print("unrecoverable error encountered inside user's PSOL routine, and control is being returned to the calling program") else
      if (idid == -14) print("the Krylov linear system solver could not achieve convergence")
    } else if (idid ==-33)  {
      print (" *** TASK TERMINATED *** ")
      print("the code has encountered trouble from which it cannot recover.  A message is printed explaining the trouble and control is returned to the calling program.")
    }
  }

#### istate

  print("")
  print("--------------------")
  print("ISTATE values")
  print("--------------------")
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
  }
    
### rstate
  print(data.frame(mess=df[1:length(ii)], val=istate[ii]))
  print("")
  print("--------------------")
  print("RSTATE values")
  print("--------------------")
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
