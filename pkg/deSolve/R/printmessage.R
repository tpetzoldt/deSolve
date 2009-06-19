## internal helper function for printing solver return code messages
## this function is not exported

printmessage <-function(df, state) {
  cat("\n", paste(formatC(1:length(df), "##",width=2), df,
              signif(state, digits = getOption("digits")), "\n"), "\n")
}
