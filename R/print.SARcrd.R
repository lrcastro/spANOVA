print.SARcrd <- function(x) {
  cat("Terms:","\n")
  trm <- data.frame(treat = c(as.character(round(x$SS[2],3)),as.character(x$DF[2])),
                    Residuals = c(as.character(round(x$SS[3],3)),as.character(x$DF[3])))
  rownames(trm) <- c("Sum of Squares","Deg. of Freedom")
  print(trm)
  rse <- sqrt(x$MS[3])
  cat("\n")
  cat("Residual standard error:",rse)
  cat("\n")
  cat("Spatial autoregressive parameter:", x$rho,"\n")
  cat("Samples considered neighbor within a",x$Par[which.min(x$Par[,3]),1],"units radius")
}
