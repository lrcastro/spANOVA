#' @name spVariog
#'
#' @title Compute empirical residual variogram for CRD or RCBD.
#'
#' @description Compute empirical residual variogram for a Completely Randomized
#' Design (CRD) or a Randomized Complete Block Design (RCBD) by a call to variog function of the
#' package geoR.
#'
#' @usage spVariog(geodata, resp = NULL, treat = NULL, block = NULL, coords = NULL,
#' data = NULL, trend = c("cte", "1st"), scale = FALSE, max.dist,
#' design = c("crd", "rcbd"), ...)
#'
#' @param geodata an object of class \code{geodata} in which the response variable
#' should be given in 'data.col', the coordinates in 'coords.col', the treatment vector
#' should be given as the first column of 'covar.col' and block as the second one.
#' @param resp either a vector of response variables or a character giving
#' the column name where it can be found in 'data'. Optional argument, just required if
#' geodata is not provided.
#' @param treat either a vector of treatment factors or a character giving
#' the column name where it can be found in 'data'. Optional argument, just required if
#' geodata is not provided.
#' @param block either a vector of block factors or a character giving the column name where
#' it can be found in 'data'. Optional argument, just required if geodata is not provided.
#' @param coords either a 2-column matrix containing the spatial coordinates or a
#' character vector giving the columns name where the coordinates can be found in 'data'.
#' Optional argument, just required if geodata is not provided.
#' @param data a data frame in which the variables specified as characters will be found.
#' Optional argument, just required if geodata is not provided.
#' @param trend type of spatial trend considered.
#' @param scale logical argument. Should the coordinates be scaled? We recommend this
#' argument to be set as TRUE if your spatial coordinates have high values as in
#' UTM coordinate system otherwise, you could get errors in the calculations.
#' See ‘Details’.
#' @param max.dist numerical value defining the maximum distance for the variogram.
#' See \code{\link[geoR]{variog}} documentation for additional information.
#' @param design type of experimental design. "crd" corresponds to Completely Randomized
#' Design and "rcbd" corresponds to Randomized Complete Block Design.
#' @param ... further arguments to be passed to \code{\link[geoR]{variog}} function.
#'
#' @details
#' This function provides a wrapper to variog to compute residual variogram for experimental designs.
#' The residuals are obtained by
#'
#' \deqn{\varepsilon = Y-X\beta,}
#'
#' where Y is the vector of response, X is the design matrix built according to the experimental design
#' chosen, and \eqn{\beta} is the vector of coefficients estimated by the linear model.
#'
#' If scale = TRUE the spatial coordinates will be scaled for numerical reasons. The scale is made by
#' subtracting the minimum spatial coordinate value from all others.
#'
#' @return
#' An object of class spVariog which is a list with the following components:
#'
#' \item{vario.res}{an object of class variogram}
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
#' @seealso \code{\link[geoR]{variog}}
#' @export

spVariog <- function(geodata, resp = NULL, treat = NULL,
                     block = NULL, coords = NULL, data = NULL,
                     trend = c("cte", "1st"), scale = FALSE,
                     max.dist, design = c("crd", "rcbd"), ...){


  trend <- as.character(trend)
  trend <- match.arg(trend)

  dln <- as.character(design)
  dln <- match.arg(design)

  if(is.logical(scale)==FALSE){
    stop("'scale' must be logical")
  }

  if(is.numeric(max.dist)==FALSE){
    stop("'max.dist' must be numeric")
  }

  if(length(max.dist) != 1){
    stop("a single numerical value must be provided in the argument max.dist")
  }


  if(missing(design)){
    stop("'design' must be specified")
  }

  if(missing(trend)){
    trend = "cte"
    warning("trend was not specified. Assuming 'trend = cte'")
  }

  classes <- c(class(resp)[1], class(treat)[1])
  numerics <- c(is.numeric(resp), is.numeric(treat))
  characters <- c(is.character(resp), is.character(treat))

  # Se geodata nao for fornecido
  if(missing(geodata)){

    if(is.null(resp)){
      stop("'resp' must be provided")
    }


    if(is.null(treat)){
      stop("'treat' must be provided")
    }

    if(is.null(coords)){
      stop("'coords' must be provided")
    }

    if(all(classes != "character")){

      # Se coords nao for matrix ou data.frame
      if(all((is.matrix(coords) == FALSE) & (is.data.frame(coords) == FALSE))){
        stop("coords must be a matrix or a data.frame")
      }

      # Se as coords nao tiverem 2 colunas
      if(dim(coords)[2] != 2){
        stop("coords must have two columns")
      }

      # se a resposta nao for numerica
      if(classes[1] != "numeric"){
        stop("resp must be numeric")
      }

    }


    # Se o objeto geodata nao for fornecido, mas sim nomes para resp, treat e coords, e
    # o data.frame data esta ausente
    if(all(classes == "character") & is.null(data)){
      stop("data must be provided")
    }


    #verificar se os argumentos principais sao do mesmo tipo
    if(any(numerics == FALSE)){
      if(any(characters == FALSE)){
        stop("Arguments should have a common class")
      }
    }

    # Se o objeto geodata nao for fornecido e  data nao for data.frame
    if(all(classes == "character") & is.data.frame(data) == FALSE){
      stop("data must be a data.frame")
    }


    # Se o objeto geodata nao for fornecido, mas resp e treat forem vetores e coords mat ou df
    if(is.vector(resp) & is.vector(treat)){
      geodata <- list(coords = coords, data = resp, covariate = as.matrix(treat))
    } else {
      stop("resp and treat must be vectors")
    }

    # Se o objeto geodata nao for fornecido, mas sim os nomes para resp, treat e coords juntamente
    # com o data.frame data
    if(all(classes == "character") & (is.null(data) == FALSE)){
      geodata <- list(coords = data[,coords], data = data[ ,resp],
                      covariate = as.matrix(data[ ,treat]))
    }

    # Se o objeto geodata nao for fornecido, mas resp e treat forem vetores e coords mat ou df
    if(design == "rcbd"){
      if(is.null(block)){
        stop("'block' must be provided")
      }

      blockn <- is.numeric(block)
      blockc <- is.character(block)

      #verificar se os argumentos principais sao do mesmo tipo
      if(all(numerics == blockn) == FALSE){
        if(all(characters == blockc) == FALSE){
          stop("Arguments should have a common class")
        }
      }

      if(is.vector(resp) & is.vector(treat) & is.vector(block)){
        geodata <- list(coords = coords, data = resp, covariate = cbind(treat, block))
      } else {
        stop("resp and treat must be vectors")
      }
      if(all(classes == "character") & (is.null(data) == FALSE)){
        geodata <- list(coords = data[,coords], data = data[ ,resp],
                        covariate = cbind(data[ ,treat], data[ ,block]))
      }
    } else {
      if(is.null(block) == FALSE) stop("'block' is just expected in an RCBD")
    }

  } else {

    if(design == "rcbd"){
      block <- geodata$covariate[,2]
    }

    if(is.geodata(geodata) == FALSE){
      stop("'geodata' must be of class geodata.")
    }

    # Se for fornecido o geodata e as variaveis principais
    if(missing(geodata) == FALSE){
      if(is.null(resp) == FALSE){
        stop("only geodata or resp must be provided")
      }

      if(is.null(treat) == FALSE){
        stop("only geodata or treat must be provided")
      }

      if(is.null(coords) == FALSE){
        stop("only geodata or coords must be provided")
      }

      # Se for RCBD e geodata estiver presente
      if((ncol(geodata$covariate) != 2) & design == "rcbd"){
        stop("Treat and Block must be specified as covariates in geodata object")
      }

      # Se for CRD e o tratamento estiver faltando
      if(is.null(geodata$covariate)){
        stop("Treat must be specified as covariate in geodata object")
      }
    }

    resp <- geodata$data
    treat <- geodata$covariate[,1]
    coords <- geodata$coords

  }

  if(scale == TRUE){
    geodata$coords[ ,1] <- geodata$coords[,1] - min(geodata$coords[,1])
    geodata$coords[ ,2] <- geodata$coords[,2] - min(geodata$coords[,2])
  }

  if(design == "crd"){
    trt <- length(unique(geodata$covariate[,1])) #numero de tratamentos
    r <- tapply(geodata$data, geodata$covariate[,1], length) # O nº de obs por categoria
    n <- sum(r)

    #primeira coluna da matriz X (vetor de 1's)
    X1 <- matrix(rep(1,n))

    #X2 <- model.matrix(~as.factor(geodata$covariate[,1])-1)

    #matriz de incidencia (sem a primeira coluna)
    row <- c(0)
    for (i in 1:trt){
      row[i+1] = sum(r[1:i])
    }

    X2 <- matrix(0,n,trt)
    for(i in 1:n){
      for(j in 1:trt){
        if(i > row[j] & i<= row[j+1]) X2[i,j] <- 1
      }
    }

    Y <- geodata$data

    if(trend == "cte"){
      X <- as.matrix(cbind(X1,X2))
      theta <- ginv(t(X)%*%X)%*%t(X)%*%Y
    }else{
      X <- as.matrix(cbind(X1, X2, geodata$coords))
      theta <- ginv(t(X)%*%X)%*%t(X)%*%Y
    }
  } else {
    if(trend != "cte"){
      trend <- "cte"
      # Se a analise RCBD for executada exibir o warning
      warning("trend is not supported in RCBD. Assuming 'trend = cte'")
    }

    Trat <- as.factor(geodata$covariate[,1])
    Bloco <- as.factor(geodata$covariate[,2])
    Y <- geodata$data

    X <- model.matrix(~ Trat + Bloco,
                 contrasts.arg=list(Trat=diag(nlevels(Trat)),
                                    Bloco=diag(nlevels(Bloco))))

    theta <- ginv(t(X)%*%X)%*%t(X)%*%Y

  }

  #Calculando o erro
  erro <- Y - X%*%theta

  #Converter os dados para a classe geodata
  geodados <- as.geodata(data.frame(geodata$coords, erro),
                         coords.col = 1:2, data.col = 3)

  #Obter o semiovariograma empírico
  vario.res <- variog(geodados, max.dist = max.dist, trend = "cte", ...)
  output <- list(vario.res = vario.res, data.geo = geodata, des.mat = X, trend = trend)
  class(output) <- c("spVariog", paste0("spVariog", design))
  return(output)
}

#' @export
plot.spVariog <- function(x, ...) {
  vari <- x$vario.res
  plot(vari, ...)
}


