#' @name spANOVAapp
#'
#' @title Shiny app for spANOVA
#'
#' @description  Shiny app for analysis of variance with spatially correlated errors
#'
#' @usage spANOVAapp(external = TRUE)
#'
#' @param external logical. If true, the system's default web browser will be
#' launched automatically after the app is started.
#'
#' @examples
#' \dontrun{
#' library(spANOVA)
#  spANOVAapp()
#' }
#'
#' @export
#' @importFrom DT datatable
#' @importFrom stats dist
#' @importFrom ape Moran.I
#' @import shiny shinyBS xtable shinythemes rmarkdown knitr
#'
spANOVAapp <- function(external = TRUE) {
  appDir <- system.file("shiny", "spANOVAapp", package = "spANOVA")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `spANOVA`.", call. = FALSE)
  }

  if(is.logical(external) == FALSE){
    stop("'external' should be logical")
  }


  shiny::runApp(appDir, display.mode = "normal", launch.browser = external)
}
