#' @name sar_crd
#'
#' @title Using a SAR model to handle spatial dependence in a Completely Randomized Design
#' @description Fit a completely randomized design when the experimental units have some degree of
#' spatial dependence using a Spatial Lag Model (SAR).
#' @usage sar_crd(resp, treat, coords, radius.min, radius.max, by)
#'
#' @param resp Numeric or complex vector containing the values of response variable.
#' @param treat Numeric or complex vector containing the treatment applied to each experimental unit.
#' @param coords Matrix of point coordinates or a SpatialPoints Object.
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
#' coords <- cbind(carrancas$X, carrancas$Y)
#' radius.min <- 50
#' radius.max <- 80
#' by <- 10
#' cv<-sar_crd(resp, treat, coords, radius.min, radius.max, by)
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
#'
#' @export


sar_crd <- function(resp, treat, coords, radius.min, radius.max, by) {

  # Defensive programming
  if(!(is.vector(resp) | is.numeric(resp))) {
    stop("'resp' must be a vector or numeric")
  }

  if(!(is.vector(treat) | is.numeric(treat))) {
    stop("'treat' must be a vector or numeric")
  }

  if(!(is.matrix(coords) | class(coords)=="SpatialPoints")) {
    stop("'coords' must be a matrix or SpatialPoints object")
  }

  if(!(is.numeric(radius.min))) {
    stop("'radius.min' must be a numeric")
  }

  if(!(is.numeric(radius.max))) {
    stop("'radius.max' must be a numeric")
  }

  if(!(is.numeric(by))) {
    stop("'radius.by' must be a numeric")
  }

  if(length(radius.min)!=1) {
    stop("'radius.min' must be of length 1")
  }

  if(length(radius.max)!=1) {
    stop("'radius.max' must be of length 1")
  }

  if(length(by)!=1) {
    stop("'by' must be of length 1")
  }


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
  treat <- factor(treat)
  aov.cl <- anova(aov(resp ~ treat))
  model <- aov(Y_ajus ~ treat)
  aov.adj <- anova(model)
  Sqt.nadj <- sum(aov.cl[,2])

  #Degres of freedom
  glrho <- 1
  gltrat <- aov.adj[1][[1]][1]
  glerror <- aov.adj[1][[1]][2]-1
  gltot <- sum(glrho,gltrat,glerror)

  #Sum of squares
  sqrho <- Sqt.nadj - sum(aov.adj[2][[1]])
  sqtrat <- aov.adj[2][[1]][1]
  sqerror <- aov.adj[2][[1]][2]
  sqtot <- sum(sqrho, sqtrat, sqerror)

  #Mean Squares
  msrho <- sqrho/glrho
  mstrat <- sqtrat/gltrat
  mserror <- sqerror/glerror

  #F statistics
  ftrat <- mstrat/mserror
  pvalue <- pf(ftrat,gltrat,glerror,lower.tail = FALSE)

  outpt <- list(DF = c(glrho, gltrat, glerror, gltot),
                SS = c(sqrho, sqtrat, sqerror),
                MS = c(msrho, mstrat, mserror),
                Fc = c(ftrat),
                p.value = c(pvalue),rho = params[best.par,"rho"],Par = params,
                y_orig = resp, treat = treat , model = model)
  class(outpt)<-c("SARanova","SARcrd",class(aov.adj))
  return(outpt)
}
