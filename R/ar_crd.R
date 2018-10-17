ar_crd <- function(resp, treat, coords, radius.min, radius.max, by) {

  seq.radius <- seq(radius.min, radius.max, by=by)
  params <- data.frame(radius = 0, rho = 0, AIC = 0)
  anova.list <- list()
  n <- length(resp)
  p.radius <- length(seq.radius)
  Y_ajus <- NULL

  for (i in 1:p.radius) {
    nb <- dnearneigh(coords, 0, seq.radius[i])
    w <- nb2mat(nb, style = "W")
    listw <- nb2listw(nb, glist = NULL, style = "W")

    # SAR model
    SAR <- lagsarlm(resp ~ treat, listw = listw,
                    method = "eigen", tol.solve = 1e-15)
    ajuste <- summary(SAR)
    rho <- as.numeric(ajuste["rho"]$rho)
    params[i, ] <- c(raio = seq.radius[i], rho = rho, AIC = AIC(SAR))
  }

  # Adjusting the data and constructing the ANOVA table
  best.par <- which.min(params$AIC)
  beta <- mean(resp)
  nb <- dnearneigh(coords, 0, seq.radius[best.par])
  w <- nb2mat(nb, style = "W")
  Y_ajus <- resp - (params[best.par,"rho"] * w%*%resp - params[best.par,"rho"] * beta)
  aov.cl <- anova(aov(resp ~ factor(treat)))
  aov.ar <- anova(aov(Y_ajus ~ factor(treat)))
  outpt <- list(DF = c(aov.ar[1][[1]],sum(aov.ar[1][[1]])), SS = c(aov.ar[2][[1]],sum(aov.ar[2][[1]])),
                MS = c(aov.ar[3][[1]],""), Fc = c(aov.ar[4][[1]],""),
                p.value = c(aov.ar[5][[1]],""),rho = params[best.par,"rho"],Par = params,
                Sqt.nadj = sum(aov.cl[,2]))
  class(outpt)<-c("arAnovaCRD",class(aov.ar))
  return(outpt)
}

# Definindo a funcao print
# print <- function(x) {
#   UseMethod("print",x)
# }

# Metodo print para a classe arAnovaCRD
print.arAnovaCRD <- function(x) {
  cat("Terms:","\n")
  trm <- data.frame(treat = c(as.character(round(x$SS[1],3)),as.character(x$DF[1])),
                    Residuals = c(as.character(round(x$SS[2],3)),as.character(x$DF[2])))
  rownames(trm) <- c("Sum of Squares","Deg. of Freedom")
  print(trm)
  rse <- sqrt(x$SS[2]/x$DF[2])
  cat("\n")
  cat("Residual standard error:",rse)
  cat("\n")
  cat("Spatial autoregressive parameter:", x$rho,"\n")
  cat("Samples considered neighbor within a",x$Par[which.min(x$Par[,3]),1],"units radius")
}
