#' @name spMVT
#' @title
#' Multiple comparison test based on multivariate t student distribution
#'
#' @description
#' Use a multivariate t student distribution to assess the equality of means.
#'
#' @usage spMVT(x, sig.level = 0.05)
#'
#' @param x a fitted model object of class SARcrd, SARrcbd or GEOanova.
#' @param sig.level a numeric value between zero and one giving the
#' significance level to use.
#'
#' @details
#' For objects of class SARcrd or SARrcbd this function performs the general
#' linear hypothesis method provided by the function \code{\link[multcomp]{glht}} on the adjusted
#' response.
#'
#' For objects of class GEOanova, the test is modified to accommodate the
#' spatial dependence among the observations as pointed out by Nogueira (2017)
#'
#' @examples
#'
#' \donttest{
#' data("crd_simulated")
#'
#' #Geodata object
#' geodados <- as.geodata(crd_simulated, coords.col = 1:2, data.col = 3,
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
#' # Gaussian Model
#' ols <- spVariofit(variograma, cov.model = "gaussian", weights = "equal",
#'                  max.dist = dist)
#'
#' lines(ols, col = 1)
#'
#' # Compute the model and get the analysis of variance table
#' mod <- aovGeo(ols, cutoff = 0.6)
#'
#' # Multivariate T test
#' spMVT(mod)
#'}
#' @return
#' a data frame containing the original mean, the spatially filtered mean and its group.
#' For the class GEOanova, the spatial dependence is filtered out using geostatistics,
#' while for the class SARanova the adjusted response based on SAR model is employed.
#'
#' @references
#' Nogueira, C. H. Testes para comparações múltiplas de
#' médias em experimentos com tendência e dependência espacial.
#' 142 f. Tese (Doutorado em Estatística e Experimentação
#' Agropecuária) | Universidade Federal de Lavras, Lavras, 2017
#'
#' @export
spMVT <- function(x, sig.level = 0.05){
  UseMethod("spMVT", x)
}

#' @export
#' @importFrom multcomp glht cld mcp
#' @method spMVT SARanova
#' @rdname spMVT

spMVT.SARanova <-function(x, sig.level = 0.05) {
  comp <-glht(x$modelAdj, linfct = mcp(treat = "Tukey"))
  let <- cld(comp, decreasing = TRUE, level = sig.level)
  mbg <- tapply(x$modelAdj$model[,1], x$modelAdj$model[,2], mean)
  mgb.orig <- tapply(x$y_orig, x$modelAdj$model[,2], mean)
  result <- data.frame(mean =  round(mgb.orig[order(mbg, decreasing = TRUE)],3),
                       filtered.mean = round(mbg[order(mbg, decreasing = TRUE)],3) ,
                       groups = let$mcletters$Letters[order(mbg,decreasing = TRUE)])
  #row.names = x$modelAdj$model[,2][order(mbg,decreasing = TRUE)]
  cat("Test based on multivariate t-student distribution","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  return(result)
}

#' @export
#' @importFrom multcompView multcompLetters
#' @importFrom mvtnorm pmvt
#' @rdname spMVT
#' @method spMVT GEOanova

spMVT.GEOanova <- function(x, sig.level = 0.05){
  # variaveis
  dados <- x$data
  resp <- dados$data
  trat <- dados$covariate[ ,1]
  D <- dados$coords

  # Atributos do modelo geoestatistico
  covMod <- x$mod
  nugget <- x$params[3]
  psill <- x$params[1]
  phi <- x$params[2]
  trend <- x$type

  Y <- c(resp)
  k <- nlevels(as.factor(trat)) ## numero de tratamentos
  r <- tapply(Y, trat, "length") ## vetor com n de repeticoes
  n <- sum(r) ## n de observacoes

  # Matriz de incidencia dos tratamentos
  #if(trend == "cte"){
  #  X2 <- x$des.mat[ ,-1]
  #}else{
  X2 <- x$des.mat[ ,2:(k+1)]
  #}

  ## sigma: a matriz de covariancia espacial
  sigma <- varcov.spatial(coords = D, cov.model = covMod, nugget = nugget,
                          kappa = 0.5, cov.pars = c(psill,phi))$varcov
  V<-sigma
  i.sig <- chol2inv(chol(V))

  # Fase 1 - Estimacao dos parametros s/ trend--------------------------------
  if(trend == "cte"){
    # 1- Obter as medias espaciais dos tratamentos

    ## medias espaciais dos tratamentos
    media <- solve(t(X2)%*%i.sig%*%X2)%*%t(X2)%*%i.sig%*%Y

    # 2- Fazer a matriz de contrastes
    C <- contr.tuk(media) ## matriz dos contrastes de tukey

    # 3 - Obter a matriz de covariancia das medias
    var.mu <- chol2inv(chol(t(X2)%*%i.sig%*%X2)) ## matriz de covariancia das medias

    # 4 - Obter a matriz W que representa a matriz de cov dos contrastes
    W <- C%*%var.mu%*%t(C) ## matriz de covariancia dos contrastes
  } else {
    # Fase 1 - Estimacao dos parametros c/ trend--------------------------------
    P <- i.sig - i.sig%*%X2%*%solve(t(X2)%*%i.sig%*%X2)%*%t(X2)%*%i.sig
    theta <- solve(t(D)%*%P%*%D)%*%t(D)%*%P%*%Y
    media <- solve(t(X2)%*%i.sig%*%X2)%*%t(X2)%*%i.sig%*%Y -
      solve(t(X2)%*%i.sig%*%X2)%*%t(X2)%*%i.sig%*%D%*%theta
    ## medias ajustadas
    s <- rankMatrix(D)[1]
    D.m <- c()
    for (i in 1:s) {
      D.m[i] <- mean(D[,i])
    }
    media.aj <- media + as.numeric(t(D.m)%*%theta)
    ## variancia das medias
    A <- D%*%solve(t(D)%*%P%*%D)%*%t(D)
    var.mu <- (solve(t(X2)%*%i.sig%*%X2) + solve(t(X2)%*%i.sig%*%X2)%*%
              t(X2)%*%i.sig%*%A%*%i.sig%*%X2%*%solve(t(X2)%*%i.sig%*%X2))
    C <- contr.tuk(media) ## matriz dos contrastes de tukey
    W <- C%*%var.mu%*%t(C)
  }

  #  Fase 2 - Aplicacao do teste---------------------------------------------------------

  S.i <- diag(1/sqrt(diag(W)))
  R <- S.i%*%W%*%S.i;R ## matriz de correlacao dos contrastes
  est <- C%*%media ## estimativa dos contrastes
  ep <- sqrt(diag(W)) ## erro padrao dos contrastes
  tc <- est/ep ## estatistica t
  glres <- sum(r) - k #ou sum(r) - k - s
  p <- nrow(C) ## numero de contrastes
  r <- with(dados,tapply(Y, trat, "length")) ## n de repeticoes
  val.p <- as.vector(0)
  for (i in 1:p){
    val.p[i] <- 1-pmvt(-rep(abs(tc[i]), p), rep(abs(tc[i]), p),
                       corr = R, df = glres)
  }
  teste <- data.frame(Estimativa=round(est,4), Erro_Padrao = round(ep,4),
                   t_calc=round(tc,3), Valor_p=round(val.p,4))

  names(val.p) <- rownames(C)


  tukeyset <- data.frame(datacol = media, faccol = unique(trat))

  tukey <- multcompLetters2(datacol ~ faccol, val.p, data = tukeyset)

  letras <- tukey$Letters

  orig.mean <- tapply(resp, trat, mean)

  print.tuk <- data.frame(mean = round(orig.mean[order(media, decreasing = TRUE)],3),
                          filtered.mean = round(unname(sort(media, decreasing = TRUE)),3),
                          groups = unname(letras),
                          row.names =  unique(trat)[order(media,decreasing = TRUE)])


  cat("Test based on multivariate t-student distribution","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  print(print.tuk)
  return(invisible(print.tuk))
}




