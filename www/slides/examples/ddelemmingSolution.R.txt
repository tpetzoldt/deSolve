
lemming <- function(t, y, parms) {
  tlag <- t - 0.74

	if (tlag <= 0)
		ylag <- 19
	else 
		ylag <- lagvalue(tlag)
  
  dy <- r * y* (1 - ylag/m)
	list(dy, ylag = ylag)
}
r = 3.5; m = 19
yinit <- c(y = 19.001)
times <- seq(from = 0, to = 40, by = 0.01)
yout <- dede(y = yinit, times = times, func = lemming, parms = NULL)
par(mfrow=c(1,2))
plot(yout[,1:2], type = "l", lwd = 2, main = "Lemming model")
plot(yout[,2:3], xlab = "y", ylab = "y(t-0.74)", type = "l", lwd = 2)
