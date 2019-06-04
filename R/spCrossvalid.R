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
