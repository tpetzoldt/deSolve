
plotstates <- function(out1, out2) {
  matplot(out1[,1], out1[,-1], type="l", ylim=c(0,2), lwd=4)
  matplot(out2[,1], out2[,-1], type="l", ylim=c(0,2), add=T, col="black")
}

## error plot
ploterrors <- function(out1, out2) {
  matplot(out1[,1], out1[,2:4] - out2[,2:4], type="l")
}

plotres <- function(out1, out2) {
  par(mfrow=c(2,1))
  plotstates(out1, out2)
  ploterrors(out1, out2)
}

