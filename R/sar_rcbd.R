#' @name sar_rcbd
#'
#' @title Using a SAR model to handle spatial dependence in a Randomized Complete Block Design
#' @description Fit a randomized complete block design when the experimental units have some degree of
#' spatial dependence using a Spatial Lag Model (SAR).
#' @usage sar_rcbd(resp, treat, coord, seq.radius = NULL)
#'
#' @param resp Numeric or complex vector containing the values of response variable.
#' @param treat Numeric or complex vector containing the treatment applied to each experimental unit.
#' @param block Numeric or complex vector specifying the blocks.
#' @param coord Matrix of point coordinates or a SpatialPoints Object.
#' @param seq.radius Complex vector containing a radii sequence used to set the neighborhood pattern.
#' The default sequence has ten numbers from 0 to half of the maximum distance between the samples.
#'
#' @return \code{sar_rcbd} returns an object of \code{\link[base]{class}} "SARanova".
#' The functions summary and anova are used to obtain and print a summary and analysis of variance
#' table of the results.
#' An object of class "SARanova" is a list containing the following components:
#'
#' \item{DF}{degrees of freedom of rho, treatments, residual and total.}
#' \item{SS}{sum of squares of rho, treatments and residual.}
#' \item{Fc}{F statistic calculated for treatment.}
#' \item{p.value}{p-value associated to F statistic for treatment.}
#' \item{rho}{the autoregressive parameter.}
#' \item{Par}{data.frame with the radii tested and its AIC.}
#' \item{y_orig}{vector of response.}
#' \item{treat}{vector of treatment applied to each experimental unit.}
#' \item{modelAdj}{model of class \code{\link[stats]{aov}} using the adjusted response.}
#' \item{namey}{response variable name.}
#' \item{namex}{treatment variable name.}
#' \item{modelstd}{data frame containing the ANOVA table using non-adjusted response.}
#'
#' @references Scolforo, Henrique Ferraço, et al. "Autoregressive spatial analysis and individual
#' tree modeling as strategies for the management of Eremanthus erythropappus." Journal of
#' forestry research 27.3 (2016): 595-603.
#'
#' @examples
#' \dontrun{
#' data("carrancas")
#' resp <- carrancas$DAP16
#' treat <- carrancas$T
#' block <- carrancas$Bloco
#' coord <- cbind(carrancas$X, carrancas$Y)
#' cv<-sar_rcbd(resp, treat, block, coord)
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
#' @export sar_rcbd print.SARrcbd summary.SARrcbd anova.SARrcbd

sar_rcbd <- function(resp, treat, block, coord, seq.radius) {

  # Defensive programming
  if(!(is.vector(resp) | is.numeric(resp))) {
    stop("'resp' must be a vector or numeric")
  }

  if(!(is.vector(treat) | is.numeric(treat))) {
    stop("'treat' must be a vector or numeric")
  }

  if(!(is.vector(block) | is.numeric(block))) {
    stop("'block' must be a vector or numeric")
  }

  if(!(is.matrix(coord) | class(coord)=="SpatialPoints")) {
    stop("'coord' must be a matrix or SpatialPoints object")
  }

  if(ncol(coord) < 2){
    stop("'coord' must have at least two columns")
  }

  if(missing(seq.radius)){
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
    SAR <- lagsarlm(resp ~ treat + block, listw = listw,
                    method = "eigen", tol.solve = 1e-15)
    ajuste <- summary(SAR)
    rho <- as.numeric(ajuste["rho"]$rho)
    params[i, ] <- c(raio = seq.radius[i], rho = rho, AIC = AIC(SAR))
  }

  treat <- factor(treat)
  block <- factor(block)
  # Adjusting the data and constructing the ANOVA table
  best.par <- which.min(params$AIC)
  beta <- mean(resp)
  nb <- dnearneigh(coord, 0, seq.radius[best.par])
  w <- nb2mat(nb, style = "W")
  Y_ajus <- resp - (params[best.par,"rho"] * w%*%resp - params[best.par,"rho"] * beta)
  aov.cl <- anova(aov(resp ~ treat + block))
  model.adj <- aov(Y_ajus ~ treat + block)
  aov.adj <- anova(model.adj)
  Sqt.nadj <- sum(aov.cl[,2])

  #Degres of freedom
  gltrat <- aov.adj[1][[1]][1]
  glblock <- aov.adj[1][[1]][2]
  glerror <- aov.adj[1][[1]][3]

  gltot <- sum(gltrat, glblock, glerror)

  #Sum of squares
  sqtrat <- aov.adj[2][[1]][1]
  sqblock <- aov.adj[2][[1]][2]
  sqerror <- aov.adj[2][[1]][3]
  sqtot <- sum(sqtrat, sqblock, sqerror)
  sqtotcor <- Sqt.nadj - sqtot

  #Mean Squares
  mstrat <- sqtrat/gltrat
  msblock <- sqblock/glblock
  mserror <- sqerror/glerror

  #F statistics
  ftrat <- mstrat/mserror
  fblock <- msblock/mserror
  pvalue.trat <- pf(ftrat, gltrat, glerror, lower.tail = FALSE)
  pvalue.block <- pf(fblock, glblock, glerror, lower.tail = FALSE)
  name.y <- paste(deparse(substitute(resp)))
  name.x <- paste(deparse(substitute(treat)))

  outpt <- list(DF = round(c(gltrat, glblock, glerror, gltot),0),
                SS = c(sqtrat, sqblock, sqerror, sqtotcor),
                MS = c(mstrat, sqblock, mserror),
                Fc = c(ftrat, fblock),
                p.value = c(pvalue.trat, pvalue.block), rho = params[best.par,"rho"],
                Par = params, y_orig = resp, treat = treat , modelAdj = model.adj,
                namey = name.y, namex = name.x, modelstd = aov.cl)
  class(outpt)<-c("SARanova","SARrcbd",class(aov.adj))
  return(outpt)
}


# Print method for this class

print.SARrcbd <- function(x) {
  cat("Response: ", x$namey, "\n")
  cat("Terms:","\n")
  trm <- data.frame(treat = c(as.character(round(x$SS[1],3)),as.character(x$DF[1])),
                    block = c(as.character(round(x$SS[2],3)),as.character(x$DF[2])),
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


# Summary method for this class

summary.SARrcbd <- function(x) {
  cat("      Summary of RCBD","\n","\n")
  cat("Parameters tested:","\n")
  print(x$Par)
  cat("\n")
  cat("Selected parameters:","\n")
  print(x$Par[which.min(x$Par[,3]),])
}

# Anova method for this class

anova.SARrcbd <- function(x, compare = FALSE) {

  if(is.logical(compare) == FALSE){
    warning("'compare' must be logical. Assuming compare == FALSE")
    compare = FALSE
  }

  cat("Analysis of Variance With Spatially Correlated Errors","\n")
  cat("\n")
  cat("Response:", ifelse(length(x$namey)>1,"resp",x$namey), "\n")
  star <- stars.pval(x$p.value[1])
  star2 <- stars.pval(x$p.value[2])
  anova.p1 <- data.frame("DF" = round(x$DF[1:4],0),
                         "SS" = round(x$SS[1:4],4),
                         "MS" = c(round(x$MS[1:3],4), ""),
                         "Fc" = c(round(x$Fc[1],4), round(x$Fc[2],4), "", ""),
                         "Pv" = c(format.pval(x$p.value[1]), format.pval(x$p.value[2]), "", ""),
                         "St" = c(star, star2, "", "")
  )
  colnames(anova.p1) <- c("Df", "Sum Sq", "Mean Sq", "F Value", "Pr(>Fc)", "")
  rownames(anova.p1) <- c("Treatment", "Block", "Residuals", "Corrected Total")
  print(anova.p1)
  cat("---","\n")
  cat("Signif. codes: ",attr(star, "legend"))

  if(compare){
    cat("\n", "\n")
    cat("---------------------------------------------------------------","\n")
    cat("Standard Analysis of Variance", "\n")
    cat("---------------------------------------------------------------")
    cat("\n")
    print(x$modelstd)
  }

}


#exportar a funcao stars.pval do pacote gtools

