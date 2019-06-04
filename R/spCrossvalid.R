#' @name spCrossvalid
#' @title Cross-validation by kriging
#'
#' @description
#' Compute cross-validation for an object of class \code{\link[spANOVA]{spVariofit}}.
#'
#' @usage spCrossvalid(x, ...)
#'
#' @param x an object of class \code{\link[spANOVA]{spVariofit}}.
#' @param ... further arguments to be passed to \code{\link[geoR]{xvalid}} function.
#'
#' @details This function is a wrapper to \code{\link[geoR]{xvalid}} function of the package geoR.
#' Please check its documentation for additional information.
#'
#' @examples
#' data("crd_simulated")
#' dados <- crd_simulated
#'
#' #Geodata object
#' geodados <- as.geodata(dados, coords.col = 1:2, data.col = 3,
#'                             covar.col = 4)
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
#'
#' #Using crossvalidation to assess the error
#' ols1.cv <- spCrossvalid(ols1)
#'
#' @export
spCrossvalid <- function(x, ...) {
  UseMethod("spCrossvalid", x)
}

#' @export
#' @method spCrossvalid spVariofit
spCrossvalid.spVariofit <- function(x, ...){
  result <- xvalid(coords = x$data.geo$coords, data = x$data.geo$data, model = x$mod, ...)
  return(result)
}
