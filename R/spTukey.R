#' @name spTukey
#' @title Compute Tukey Honest Significant Differences for a Spatially Correlated Model
#'
#' @description
#' Perform multiple comparisons of means treatments based on the Studentized range
#' statistic when the errors are spatially correlated.
#'
#' @usage
#' spTukey(x, sig.level = 0.05)
#'
#' @param x A fitted model object of class SARcrd, SARrcbd or GEOanova.
#' @param sig.level A numeric value between zero and one giving the significance
#' level to use.
#'
#' @details
#' For objects of class SARcrd or SARrcbd this function performs the standard
#' Tukey's ‘Honest Significant Difference’ method provided by the function
#' \code{\link[stats]{TukeyHSD}} on the adjusted response.
#'
#' For objects of class GEOanova, the method is modified to take into account
#' the spatial dependence among the observations. First, we estimate a
#' contrast matrix (\eqn{C}) using cont.tuk function and then after estimate the
#' spatial mean of each treatment (\eqn{\mu_i}) we can assess the significance of
#' the contrast by
#'
#' \deqn{|c_i \mu_i| > {HSD}_i}
#'
#' where \eqn{HSD_i = q(\alpha, k, \nu) * sqrt(0.5*{w}_ii)} and
#' \eqn{k} is the number of treatments, \eqn{\alpha} is the level of significance,
#' \eqn{\nu} is the degree of freedom of the model, \eqn{{w}_ii} is the variance of
#' the i-th contrast.
#'
#' @return
#' a data frame containing the original mean, the spatially filtered mean and its group.
#' For the class GEOanova, the spatial dependence is filtered out using geostatistics,
#' while for the class SARanova the adjusted response based on SAR model is employed.
#'
#' @references
#' NOGUEIRA, C. H. Testes para comparações múltiplas de
#' médias em experimentos com tendência e dependência espacial.
#' 142 f. Tese (Doutorado em Estatística e Experimentação
#' Agropecuária) | Universidade Federal de Lavras, Lavras, 2017
#'
#' @examples
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
#' # Tukey's HSD
#' spTukey(mod)
#'
#'
#' @export
spTukey <- function(x, sig.level = 0.05) {
  UseMethod("spTukey", x)
}

#' @export
#' @importFrom multcompView multcompLetters2 multcompLetters3
#' @importFrom stats TukeyHSD
#' @rdname spTukey
#' @method spTukey SARcrd

spTukey.SARcrd <- function(x, sig.level = 0.05){
  exp_tukey <- TukeyHSD(x$modelAdj, conf.level = 1 - sig.level)
  let <- multcompLetters3(2, 1, exp_tukey$treat[,"p adj"], x$modelAdj$model)
  mbg <- tapply(x$modelAdj$model[,1], x$modelAdj$model[,2], mean)
  mgb.orig <- round(tapply(x$y_orig, x$modelAdj$model[,2], mean),3)
  result <- data.frame(mean = round(mgb.orig[order(mbg, decreasing = TRUE)],3),
                       filtered.mean = round(mbg[order(mbg,decreasing = TRUE)],3) ,
                       groups = let$Letters)
  cat("Tukey's Test","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  return(result)
}

#' @export
#' @importFrom multcompView multcompLetters2
#' @rdname spTukey
#' @method spTukey SARrcbd

spTukey.SARrcbd <- function(x, sig.level = 0.05){
  exp_tukey <- TukeyHSD(x$modelAdj, conf.level = 1 - sig.level)
  let <- multcompLetters3(2, 1, exp_tukey$treat[,"p adj"], x$modelAdj$model)
  mbg <- tapply(x$modelAdj$model[,1], x$modelAdj$model[,2], mean)
  result <- data.frame(means = round(mbg[order(mbg,decreasing = TRUE)],3) ,
                       groups = let$Letters)
  cat("Tukey's Test","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  return(result)
}

#' @export
#' @rdname spTukey
#' @method spTukey GEOanova

spTukey.GEOanova <- function(x, sig.level = 0.05){

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
  n <- sum(r) ## n de observações

  # # Matriz de incidencia dos tratamentos
  # row <- c(0)
  # for (i in 1:k){
  #   row[i+1] = sum(r[1:i])
  # }
  #
  # X2 <- matrix(0,n,k)
  # for(i in 1:n){
  #   for(j in 1:k){
  #     if(i > row[j] & i<= row[j+1]) X2[i,j] <- 1
  #   }
  # }
  #
  # if(trend != "cte"){
  #   X1 <- rep(1, n)
  #   X2 <- as.matrix(cbind(X1, X2, D))
  # }

  # if(trend == "cte"){
  #
  #   #X2 <- x$des.mat[ ,-1]
  # }else{
  X2 <- x$des.mat[ ,2:(k+1)]
  # }

  ## sigma: a matriz de covariancia espacial
  sigma <- varcov.spatial(coords = D, cov.model = covMod, nugget = nugget,
                          kappa = 0.5, cov.pars = c(psill,phi))$varcov
  V <- sigma
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

  #Aplicacao do teste
  glres <- sum(r) - k ## ou glres= sum(r) - k - s
  q <- qtukey(sig.level, k, glres, lower.tail = F)
  dms <- q*sqrt(0.5*diag(W))
  est <- C%*%media

  result <- c()
  p <- choose(k,2) ### numero de contrastes
  for (i in 1:p){
    if ((abs(est[i])-dms[i]) >= 0 ){
      result[i] <- "Sig"
    } else {
      result[i] <- "No Sig"
    }
  }

  teste <- data.frame(Estimativa = round(abs(est), 3), DMS=round(dms, 3),
                      Resultado = result)


  dif <- c(abs(est)-dms >= 0)
  names(dif) <- rownames(C)
  #tukey <- multcompLetters(dif, rev = TRUE)

  tukeyset <- data.frame(datacol = media, faccol = unique(trat))
  tukey <- multcompLetters2(datacol ~ faccol, dif, data = tukeyset)

  letras <- tukey$Letters


  orig.mean <- tapply(resp, trat, mean)

  print.tuk <- data.frame(mean = round(orig.mean[order(media, decreasing = TRUE)],3),
                          filtered.mean = round(unname(sort(media, decreasing = TRUE)),3),
                          groups = unname(letras),
                          row.names =  unique(trat)[order(media,decreasing = TRUE)])
 # print.tuk

  #print.tuk
  cat("Tukey's Test","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  print(print.tuk)
  return(invisible(print.tuk))
}

#@rdname spTukey
#@method spTukey GEOrcbd
#@export

# spTukey.GEOrcbd <- function(x, sig.level = 0.05){
#
#   # variaveis
#   dados <- x$data
#   resp <- dados$data
#   D <- dados$coords
#   trat <- dados$covariate[ ,1]
#   block <- dados$covariate[,2]
#
#   # Atributos do modelo geoestatistico
#   covMod <- x$mod
#   nugget <- x$params[3]
#   psill <- x$params[1]
#   phi <- x$params[2]
#   trend <- "cte"
#
#   Y <- c(resp)
#   k <- nlevels(as.factor(trat)) ## numero de tratamentos
#   r <- tapply(Y, trat, "length") ## vetor com n de repeticoes
#   n <- sum(r) ## n de observações
#
#   # Matriz de incidencia dos tratamentos
#   blk <- length(unique(block)) #numero de blocos
#   X1 <- as.matrix(x$des.mat[ ,1], ncol = 1)
#   X2 <- x$des.mat[ ,-c(1,(k+2):(blk+k+1))]
#
#   ## sigma: a matriz de covariancia espacial
#   sigma <- varcov.spatial(coords = D, cov.model = covMod, nugget = nugget,
#                           kappa = 0.5, cov.pars = c(psill,phi))$varcov
#   V <- sigma
#   i.sig <- chol2inv(chol(V))
#
#   # Fase 1 - Estimacao dos parametros c/ trend--------------------------------
#   P <- i.sig - i.sig%*%X2%*%solve(t(X2)%*%i.sig%*%X2)%*%t(X2)%*%i.sig
#   theta <- solve(t(D)%*%P%*%D)%*%t(D)%*%P%*%Y
#   media <- solve(t(X2)%*%i.sig%*%X2)%*%t(X2)%*%i.sig%*%Y -
#     solve(t(X2)%*%i.sig%*%X2)%*%t(X2)%*%i.sig%*%D%*%theta
#   ## medias ajustadas
#   s <- rankMatrix(D)[1]
#   D.m <- c()
#   for (i in 1:s) {
#     D.m[i] <- mean(D[,i])
#   }
#
#   media.aj <- media + as.numeric(t(D.m)%*%theta)
#
#   ## variancia das medias
#   A <- D%*%solve(t(D)%*%P%*%D)%*%t(D)
#   var.mu <- (solve(t(X2)%*%i.sig%*%X2) + solve(t(X2)%*%i.sig%*%X2)%*%
#                t(X2)%*%i.sig%*%A%*%i.sig%*%X2%*%solve(t(X2)%*%i.sig%*%X2))
#   C <- contr.tuk(media) ## matriz dos contrastes de tukey
#   W <- C%*%var.mu%*%t(C)
#
#   #Aplicacao do teste
#   glres <- x$DF[3] ## ou glres= sum(r) - k - s
#   q <- qtukey(sig.level, k, glres, lower.tail = F)
#   dms <- q*sqrt(0.5*diag(W))
#   est <- C%*%media
#
#   result <- c()
#   p <- choose(k,2) ### numero de contrastes
#   for (i in 1:p){
#     if ((abs(est[i])-dms[i]) >= 0 ){
#       result[i] <- "Sig"
#     } else {
#       result[i] <- "No Sig"
#     }
#   }
#
#   teste <- data.frame(Estimativa = round(abs(est), 3), DMS=round(dms, 3),
#                       Resultado = result)
#
#
#   dif <- c(abs(est)-dms >= 0)
#   names(dif) <- rownames(C)
#   #tukey <- multcompLetters(dif, rev = TRUE)
#
#   tukeyset <- data.frame(datacol = media, faccol = unique(trat))
#   tukey <- multcompLetters2(datacol ~ faccol, dif, data = tukeyset)
#
#   letras <- tukey$Letters
#
#
#
#   print.tuk <- data.frame(means = round(unname(sort(media, decreasing = TRUE)),3),
#                           groups = unname(letras),
#                           row.names =  unique(trat)[order(media,decreasing = TRUE)])
#   #print.tuk
#
#   #print.tuk
#   cat("Tukey's Test","\n")
#   cat("\n")
#   cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
#   cat("\n")
#   print(print.tuk)
#   return(invisible(print.tuk))
# }

