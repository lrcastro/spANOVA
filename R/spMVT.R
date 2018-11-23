#' @export

spMVT <- function(x, sig.level = 0.05){
  UseMethod("spMVT", x)
}


#' @importFrom multcomp glht cld mcp
#' @export

spMVT.SARanova <-function(x, sig.level = 0.05) {
  comp <-glht(x$modelAdj, linfct = mcp(treat = "Tukey"))
  let <- cld(comp,decreasing = TRUE, level = sig.level)
  mbg <- tapply(x$modelAdj$model[,1], x$modelAdj$model[,2], mean)
  result <- data.frame(Treatment = order(mbg,decreasing = TRUE),
                       Mean = mbg[order(mbg,decreasing = TRUE)] ,
                       Groups = let$mcletters$Letters[order(mbg,decreasing = TRUE)])
  cat("Test based on multivariate t-student distribution","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different","\n")
  cat("\n")
  return(result)
}
