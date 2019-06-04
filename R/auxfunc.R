#' @name contr.tuk
#' @title Tukey's contrast matrix
#'
#' @description Compute Tukey's contrast matrix
#'
#' @usage contr.tuk(x)
#'
#' @param x a vector of means.
#'
#' @details
#' Computes the matrix of contrasts for comparisons of mean levels.
#'
#' @return
#' The matrix of contrasts with treatments in row names is returned
#'
#' @export
contr.tuk <- function(x){
  k <- length(x)
  kindx <- 1:k
  CM <- c()
  rnames <- c()
  if (!is.null(names(x)))
    varnames <- names(x) else varnames <- 1:length(x)
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      CM <- rbind(CM, as.numeric(kindx == i) - as.numeric(kindx == j))
      rnames <- c(rnames, paste(varnames[i], varnames[j], sep="-"))
    } }
  rownames(CM) <- rnames
  return(CM)
}

