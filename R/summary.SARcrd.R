#' @export

summary.SARcrd <- function(x) {
  cat("      Summary of CRD","\n","\n")
  cat("Parameters tested:","\n")
  print(x$Par)
  cat("\n")
  cat("Selected parameters:","\n")
  print(x$Par[which.min(x$Par[,3]),])
}
