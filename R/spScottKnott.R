#' @export

spScottKnott <- function(x, sig.level = 0.05) {
  UseMethod("spScottKnott", x)
}


#' @importFrom ScottKnott SK
#' @importFrom utils capture.output
#' @export


spScottKnott.SARanova <- function(x, sig.level = 0.05) {
  invisible(capture.output(out <- summary(SK(x=x$modelAdj$model,
                                             model = 'Y_ajus ~ treat', which = 'treat',
                                             dispersion = 's',
                                             sig.level = sig.level))))
  colnames(out) <- c("Treatment", "Mean", "Groups")
  cat("Scott-Knott Test","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different","\n")
  cat("\n")
  return(out)
}
