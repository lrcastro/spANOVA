#' @name spVariofit
#'
#' @title Fit a variogram model
#'
#' @description Fit a parametric model to a variogram created by the function spVariog.
#'
#' @usage spVariofit(x, ...)
#'
#' @param x an object of class \code{\link[spANOVA]{spVariog}}.
#' @param ... further arguments to be passed to \code{\link[geoR]{variofit}} function.
#'
#' @details
#' This function is a wrapper to \code{\link[geoR]{variofit}} and can be used to fit a
#' parametric model to a variogram using either ordinary least squares or weighted least squares.
#' It takes as the main argument a spVariog object and others arguments should be passed
#' to \code{...} such as "cov.model" and so on.
#'
#' @return an object of class \code{SpVariofit} which is a list containing the following
#' components:
#'
#' \item{mod}{an object of class \code{\link[geoR]{variofit}}}
#' \item{data.geo}{an object of class geodata}
#' \item{des.mat}{the design matrix}
#' \item{trend}{a character specifying the type of spatial trend}
#'
#' @examples
#' data("crd_simulated")
#' dados <- crd_simulated
#'
#' #Geodata object
#' geodados <- as.geodata(dados, coords.col = 1:2, data.col = 3,
#'                       covar.col = 4)
#' h_max <- summary(geodados)[[3]][[2]]
#' dist <- 0.6*h_max
#'
#' # Computing the variogram
#' variograma <- spVariog(geodata = geodados,
#'                       trend = "cte", max.dist = dist, design = "crd",
#'                       scale = FALSE)
#'
#' plot(variograma, ylab = "Semivariance", xlab = "Distance")
#'
#' # Spherical Model
#' ols1 <- spVariofit(variograma, cov.model = "spherical", weights = "equal",
#'                   max.dist = dist)
#' lines(ols1)
#'
#'
#' @seealso
#' \code{\link[geoR]{variofit}}
#' @export
spVariofit <- function(x, ...) {
  UseMethod("spVariofit", x)
}

#' @export
#' @method spVariofit spVariogcrd
spVariofit.spVariogcrd <- function(x, ...){
  mod <- variofit(x$vario.res, ...)
  result <- list(mod = mod, data.geo = x$data.geo, des.mat = x$des.mat, trend = x$trend)
  class(result) <- c("spVariofit", "spVariofitCRD")
  return(result)
}

#' @export
#' @method spVariofit spVariogrcbd
spVariofit.spVariogrcbd <- function(x, ...){
  mod <- variofit(x$vario.res, ...)
  result <- list(mod = mod, data.geo = x$data.geo, des.mat = x$des.mat, trend = x$trend)
  class(result) <- c("spVariofit", "spVariofitRCBD")
  return(result)
}

#' @export
#' @method lines spVariofit
lines.spVariofit <- function(x, ...){
  lines(x$mod, ...)
}
