#' @importFrom agricolae HSD.test
#' @export

spTukey.SARcrd <- function(x, sig.level = 0.05){
  out <- HSD.test(x$model$model[,1], x$model$model[,2],
                  MSerror = x$MS[3], DFerror = x$DF[3], group=TRUE, alpha = sig.level)
  result <- cbind(rownames(out$group), out$group)
  colnames(result) <- c("Treatments", "Mean", "Groups")
  cat("Tukey's Test","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different","\n")
  cat("\n")
  return(result)
}
