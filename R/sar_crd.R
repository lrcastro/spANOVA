#' @name sar_crd
#'
#' @title Using a SAR model to handle spatial dependence in a Completely Randomized Design
#' @description Fit a completely randomized design when the experimental units have some degree of
#' spatial dependence using a Spatial Lag Model (SAR).
#' @usage sar_crd(resp, treat, coord, radius.min, radius.max, by)
#'
#' @param resp Numeric or complex vector containing the values of response variable.
#' @param treat Numeric or complex vector containing the treatment applied to each experimental unit.
#' @param coord Matrix of point coordinates or a SpatialPoints Object.
#' @param ray.min Numeric value specifyng the minimum ray to be considered as a neighbor.
#' @param ray.max Numeric value specifying the maximum ray to be considered as a neighbor.
#' @param by Numeric value specifying the increment of the ray sequence.
#'
#' @return \code{sar_crd} returns an object of \code{\link[base]{class}} "SARanova".
#' The functions summary and anova are used to obtain and print a summary and analysis of variance
#' table of the results.
#' An object of class "SARanova" is a list containing the following components:
#'
#' \item{DF}{degrees of freedom of rho, treatments, residual and total.}
#' \item{SS}{sum of squares of rho, treatments and residual.}
#' \item{Fc}{F statistic calculated for treatment.}
#' \item{p.value}{p-value associated to F statistic for treatment.}
#' \item{rho}{the autoregressive parameter.}
#' \item{Par}{data.frame with the radius tested and its AIC.}
#' \item{y_orig}{vector of response.}
#' \item{treat}{vector of treatment applied to each experimental unit.}
#' \item{fitted}{model of class \code{\link[stats]{aov}} using the adjusted response.}
#'
#' @references Long, D. S. "Spatial statistics for analysis of variance of agronomic field trials."
#' Practical handbook of spatial statistics. CRC Press, Boca Raton, FL (1996): 251-278.
#'
#' @examples
#' \dontrun{
#' data("carrancas")
#' resp <- carrancas$DAP16
#' treat <- carrancas$T
#' coord <- cbind(carrancas$X, carrancas$Y)
#' radius.min <- 50
#' radius.max <- 80
#' by <- 10
#' cv<-sar_crd(resp, treat, coord, radius.min, radius.max, by)
#' cv
#'
#' #Summary for class SARanova
#' summary(cv)
#'
#' #Anova for class SARanova
#' anova(cv)
#'
#' #Test based on multivariate t-student distribution
#' spMVT(cv)
#'
#' #Tukey's test
#' spTukey(cv)
#'
#' #Scott-Knott test
#' spScottKnott(cv)
#'
#' }
#'
#' @import spdep
#' @importFrom gtools stars.pval
#' @export sar_crd print.SARcrd summary.SARcrd anova.SARcrd


sar_crd <- function(resp, treat, coord, radius.seq) {

  # Defensive programming
  if(!(is.vector(resp) | is.numeric(resp))) {
    stop("'resp' must be a vector or numeric")
  }

  if(!(is.vector(treat) | is.numeric(treat))) {
    stop("'treat' must be a vector or numeric")
  }

  if(!(is.matrix(coord) | class(coord)=="SpatialPoints")) {
    stop("'coord' must be a matrix or SpatialPoints object")
  }

  if(ncol(coord) < 2){
    stop("'coord' must have at least two columns")
  }

  if(missing(radius.seq)){
    max.dist <- max(dist(coord))
    seq.radius <- seq(0, 0.5*max.dist, l = 11)[-1]
  }

  params <- data.frame(radius = 0, rho = 0, AIC = 0)
  anova.list <- list()
  n <- length(resp)
  p.radius <- length(seq.radius)
  Y_ajus <- NULL

  for (i in 1:p.radius) {
    nb <- dnearneigh(coord, 0, seq.radius[i])
    w <- try(nb2mat(nb, style = "W"), silent = TRUE)
    test <- grepl("Error", w)

    # Se caso nao forem encontradas amostras dentro do raio especificado
    k <- 0.1 # incremento
    while(test[1] == TRUE){
      seq.radius <- seq(0, (0.5+k)*max.dist, l = 11)[-1]
      nb <- dnearneigh(coord, 0, seq.radius[i])
      w <- try(nb2mat(nb, style = "W"), silent = TRUE)
      test <- grepl("Error", w)
      k <- k + 0.1
    }

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
  nb <- dnearneigh(coord, 0, seq.radius[best.par])
  w <- nb2mat(nb, style = "W")
  Y_ajus <- resp - (params[best.par,"rho"] * w%*%resp - params[best.par,"rho"] * beta)
  treat <- factor(treat)
  aov.cl <- anova(aov(resp ~ treat))
  model <- aov(Y_ajus ~ treat)
  aov.adj <- anova(model)
  Sqt.nadj <- sum(aov.cl[,2])

  #Degres of freedom
  gltrat <- aov.adj[1][[1]][1]
  glerror <- aov.adj[1][[1]][2]
  gltot <- sum(gltrat,glerror)

  #Sum of squares
  sqtrat <- aov.adj[2][[1]][1]
  sqerror <- aov.adj[2][[1]][2]
  sqtot <- sum(sqtrat, sqerror)
  sqtotcor <- Sqt.nadj - sqtot

  #Mean Squares
  mstrat <- sqtrat/gltrat
  mserror <- sqerror/glerror

  #F statistics
  ftrat <- mstrat/mserror
  pvalue <- pf(ftrat,gltrat,glerror,lower.tail = FALSE)
  name.y <- paste(deparse(substitute(resp)))
  name.x <- paste(deparse(substitute(treat)))

  outpt <- list(DF = round(c(gltrat, glerror, gltot),0),
                SS = c(sqtrat, sqerror, sqtotcor),
                MS = c(mstrat, mserror),
                Fc = c(ftrat),
                p.value = c(pvalue), rho = params[best.par,"rho"], Par = params,
                y_orig = resp, treat = treat , model = model, namey = name.y,
                namex = name.x)
  class(outpt)<-c("SARanova","SARcrd",class(aov.adj))
  return(outpt)
}


# Print method for this class

print.SARcrd <- function(x) {
  cat("Response: ", x$namey, "\n")
  cat("Terms:","\n")
  trm <- data.frame(treat = c(as.character(round(x$SS[2],3)),as.character(x$DF[2])),
                    Residuals = c(as.character(round(x$SS[3],3)),as.character(x$DF[3])))
  rownames(trm) <- c("Sum of Squares","Deg. of Freedom")
  print(trm)
  rse <- sqrt(x$MS[2])
  cat("\n")
  cat("Residual standard error:",rse)
  cat("\n")
  cat("Spatial autoregressive parameter:", x$rho,"\n")
  cat("Samples considered neighbor within a",x$Par[which.min(x$Par[,3]),1],"units radius")
}


# Summary method for this class

summary.SARcrd <- function(x) {
  cat("      Summary of CRD","\n","\n")
  cat("Parameters tested:","\n")
  print(x$Par)
  cat("\n")
  cat("Selected parameters:","\n")
  print(x$Par[which.min(x$Par[,3]),])
}

# Anova method for this class

anova.SARcrd <- function(x) {
  cat("Analysis of Variance With Spatially Correlated Errors","\n")
  cat("\n")
  cat("Response:", ifelse(length(x$namey)>1,"resp",x$namey), "\n")
  star <- stars.pval(x$p.value)
  anova.p1 <- data.frame("DF" = round(x$DF[1:3],0),
                         "SS" = round(x$SS[1:3],4),
                         "MS" = c(round(x$MS[1:2],4), ""),
                         "Fc" = c(round(x$Fc,4), "", ""),
                         "Pv" = c(format.pval(x$p.value), "", ""),
                         "St" = c(star, "", "")
  )
  colnames(anova.p1) <- c("Df", "Sum Sq", "Mean Sq", "F Value" ,"Pr(>Fc)", "")
  rownames(anova.p1) <- c("Treatment","Residuals","Corrected Total")
  print(anova.p1)
  cat("---","\n")
  cat("Signif. codes: ",attr(star, "legend"))
}


#exportar a funcao stars.pval do pacote gtools

