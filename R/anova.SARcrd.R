#' @export

anova.SARcrd <- function(x) {
  cat("Analysis of Variance With Spatially Correlated Errors","\n")
  anova.p1 <- data.frame("DF" = x$DF[1:3],
                         SS = x$SS[1:3],
                         MS = x$MS[1:3],
                         Fc = c(NA,x$Fc,NA),
                         "Pv" = c(NA,cv$p.value,NA))

  colnames(anova.p1) <- c("DF", "SS", "MS", "Fc" ,"Pr(>Fc)")
  rownames(anova.p1) <- c("Rho","Treatment","Residuals")
  printCoefmat(anova.p1, P.values=TRUE, na.print="", check.names=FALSE)
}
