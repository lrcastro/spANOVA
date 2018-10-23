#' @importFrom multcomp glht cld
#' @export

spMVT.SARcrd <-function(x, sig.level = 0.05) {
  comp <-glht(x$model, linfct = mcp(treat = "Tukey"))
  let <- cld(comp,decreasing = TRUE, level = sig.level)
  mbg <- tapply(x$model$model[,1], x$model$model[,2], mean)
  result <- data.frame(Treatment = order(mbg,decreasing = TRUE),
                       Mean = mbg[order(mbg,decreasing = TRUE)] ,
                       Groups = let$mcletters$Letters[order(mbg,decreasing = TRUE)])
  cat("Test based on multivariate t-student distribution","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different","\n")
  cat("\n")
  return(result)
}
