## internal helper functions for printing solver return code messages
## these functions are not exported

## print combined messages (message and numeric output)
printmessage <-function(message, state) {
  cat("\n", paste(formatC(1:length(message), "##",width=2), message,
              signif(state, digits = getOption("digits")), "\n"), "\n")
}

## print short messages
printM <- function(message) cat(message, "\n")
