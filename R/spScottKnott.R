#' @name spScottKnott
#' @title The Scott-Knott Clustering Algorithm
#'
#' @description
#' This function implements the Scott-Knott Clustering
#' Algorithm for objects of class SARcrd, SARrcbd, and GEOanova.
#'
#' @usage spScottKnott(x, sig.level = 0.05)
#'
#' @param x a fitted model object of class SARcrd, SARrcbd or GEOanova.
#' @param sig.level a numeric value between zero and one giving the significance
#' level to use.
#'
#' @details
#' For objects of class SARcrd or SARrcbd this function performs the standard Scott-Knott
#' Clustering Algorithm provided by the function \code{\link[ScottKnott]{SK}} on the
#' adjusted response.
#'
#' For objects of class GEOanova, the method is modified to take into account the spatial
#' dependence among the observations. The method is described in Nogueira (2017).
#'
#' @return
#' a data frame containing the means and its group
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
#'
#' # Compute the model and get the analysis of variance table
#' mod <- aovGeo(ols, cutoff = 0.6)
#'
#' # Scott-Knott clustering algorithm
#' spScottKnott(mod)
#'
#' @references
#' Nogueira, C. H. Testes para comparações múltiplas de
#' médias em experimentos com tendência e dependência espacial.
#' 142 f. Tese (Doutorado em Estatística e Experimentação
#' Agropecuária) | Universidade Federal de Lavras, Lavras, 2017
#'
#' @export
spScottKnott <- function(x, sig.level = 0.05) {
  UseMethod("spScottKnott", x)
}

#' @export
#' @importFrom ScottKnott SK
#' @rdname spScottKnott
#' @method spScottKnott SARanova

spScottKnott.SARanova <- function(x, sig.level = 0.05) {
  invisible(capture.output(out <- summary(SK(x = x$modelAdj,
                                             which = 'treat',
                                             dispersion = 's',
                                             sig.level = sig.level))))
  mgb.orig <- round(tapply(x$y_orig, x$modelAdj$model[,2], mean),3)
  out <- cbind(sort(mgb.orig, decreasing = T), out[,-1])
  colnames(out) <- c("mean","filtered.mean", "groups")
  cat("Scott-Knott Test","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  return(out)
}

#' @export
#' @rdname spScottKnott
#' @method spScottKnott GEOanova

spScottKnott.GEOanova <- function(x, sig.level = 0.05){
  # variaveis
  dados <- x$data
  resp <- dados$data
  trat <- dados$covariate[ ,1]
  D <- dados$coords
  #design <- x$des

  # Atributos do modelo geoestatistico
  covMod <- x$mod
  nugget <- x$params[3]
  psill <- x$params[1]
  phi <- x$params[2]
  trend <- x$type

  y <- c(resp)
  k <- nlevels(as.factor(trat)) ## numero de tratamentos
  r <- tapply(y, trat, "length") ## vetor com n de repeticoes
  n <- sum(r) ## n de observacoes

  # Matriz de incidencia dos tratamentos
  #if(trend == "cte"){
    #X2 <- x$des.mat[ ,-1]
  #}else{
   X2 <- x$des.mat[ ,2:(k+1)]
  #}


  ## sigma: a matriz de covariancia espacial
  sigma <- varcov.spatial(coords = D, cov.model = covMod, nugget = nugget,
                          kappa = 0.5, cov.pars = c(psill,phi))$varcov
  V <- sigma/sigma[1,1]
  i.V <- chol2inv(chol(V))

  if (trend != "cte") {
    P <- i.V- i.V%*%X2%*%solve(t(X2)%*%i.V%*%X2)%*%t(X2)%*%i.V
    theta <- solve(t(D)%*%P%*%D)%*%t(D)%*%P%*%y
    beta <- solve(t(X2)%*%i.V%*%X2)%*%t(X2)%*%i.V%*%y -
      solve(t(X2)%*%i.V%*%X2)%*%t(X2)%*%i.V%*%D%*%theta
  }else{
    beta <- solve(t(X2)%*%i.V%*%X2)%*%t(X2)%*%i.V%*%y
  }

  trt <- c()
  for (i in 1:k){
    trt <- as.factor(c(trt, rep(i, r[i])))
  }

  dados.na <- data.frame(obs = 1:length(y), y = y, trat = trt,
                      D1 = D[,1], D2 = D[,2])
  obs <- c()
  for (i in 1:k){
    obs <- c(obs, subset(dados.na, (trat==order(beta)[i]))$obs)
  }

  dados.ord <- dados.na[obs,]
  beta.ord <- sort(beta)
  sigma <- varcov.spatial(dados.ord[,4:5], cov.model = covMod, nugget = nugget,
                       cov.pars = c(psill,phi))$varcov
  V <- sigma/sigma[1,1]
  i.V <- chol2inv(chol(V))

  B <- c()
  for (i in 1:(k-1)){
    n1 <- sum(r[1:i])
    n2 <- sum(r[(i+1):k])
    n <- n1+n2
    K1 <- matrix(c(rep(1,n1),rep(0,n2)),ncol=1)
    K2 <- matrix(c(rep(0,n1),rep(1,n2)),ncol=1)
    K <- cbind(K1,K2)
    X1 <- matrix(c(rep(1,n)),ncol=1)
    Y <- c(dados.ord$y)
    sqd <- t(Y)%*%i.V%*%Y- t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%
      t(X1)%*%i.V%*%Y
    sqn <- t(Y)%*%i.V%*%K%*%solve(t(K)%*%i.V%*%K)%*%t(K)%*%i.V%*%Y -
      t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%t(X1)%*%i.V%*%Y
    B <- c(B, sqn/sqd)}
  posmax <- grep(max(B),B)
  lamb <- as.numeric((pi/(2*(pi-2)))*n*B[posmax])
  pval <- pchisq(lamb, k/(pi-2), lower.tail = FALSE)
  groups <- rep(0, times = k)
  ng <- 1
  if (pval > sig.level) {
    groups[1:k] <- 1
  } else {
    groups[1:posmax] <- 1
    groups[(posmax+1):k] <- 2
    ng <- ng + 1
  }
  if ((any(groups!=1))) {
    continua1 = TRUE
  } else {
    continua1 = FALSE
  }
  if (continua1 == TRUE) {
    posI <- 1
    fim <- FALSE
    while(posI <= k) {
      ini <- groups[posI]
      posF <- max(which(groups == ini))
      if (((posF - posI) > 0) & (ini > 0)) {
        B <- c()
        for (i in 1:(posF-posI)){
          n1 <- sum(r[posI:(posI+i-1)])
          n2 <- sum(r[(i+posI):posF])
          n <- n1+n2
          K1 <- matrix(c(rep(1,n1),rep(0,n2)),ncol=1)
          K2 <- matrix(c(rep(0,n1),rep(1,n2)),ncol=1)
          K <- cbind(K1,K2)
          X1 <- matrix(c(rep(1,n)),ncol=1)
          if (posI==1) { g=0 } else { g=sum(r[1:(posI-1)]) }
          dados.ord1 <- dados.ord[(g+1):(g+n),]
          Y <- c(dados.ord1$y)
          sigma <- varcov.spatial(dados.ord1[,4:5],cov.model = covMod, nugget = nugget,
                               cov.pars = c(psill,phi))$varcov
          V <- sigma/sigma[1,1]
          i.V <- chol2inv(chol(V))
          sqd <- t(Y)%*%i.V%*%Y- t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%
            t(X1)%*%i.V%*%Y
          sqn <- t(Y)%*%i.V%*%K%*%solve(t(K)%*%i.V%*%K)%*%t(K)%*%i.V%*%Y -
            t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%t(X1)%*%i.V%*%Y
          B <- c(B, sqn/sqd) }
        posmax <- grep(max(B),B)
        lamb <- as.numeric((pi/(2*(pi-2)))*n*B[posmax])
        pval <- pchisq(lamb, (posF-posI+1)/(pi-2), lower.tail = FALSE)
        aux <- groups[posI]
        if (pval > sig.level) {
          groups[posI:posF] <- - aux
          posI <- posF + 1
          if (posI > k) {
            posI <- 1
          }
        }
        else {
          groups[(posI+posmax):posF] <- ng + 1
          ng <- ng + 1
          posI <- 1
        }
      }
      else posI <- posF + 1
      #if (posI >= k) fim <- TRUE
      #if (fim == TRUE) break
    }
    groups <- abs(groups)
  }
  Teste <- data.frame(trat=levels(factor(trat))[order(beta)],
                     medias=beta.ord, grupos=groups)
  k <- nrow(Teste)
  ordertest <- Teste[order(Teste[, 2], decreasing = TRUE),]
  M <- rep("", k)
  M[1] <- "a"
  j <- 2
  for (i in 2:k) {
    if (ordertest[i, 3] == ordertest[i-1, 3]) {
      M[i] <- M[i-1] }
    else { M[i] <- letters[j]
    j <- j + 1 }
  }

  orig.mean <- tapply(resp, trat, mean)
  saida <- data.frame(mean = round(orig.mean[order(Teste[, 2], decreasing = TRUE)], 3) ,
                      filtered.mean = unname(round(ordertest[, 2],3)),
                      groups = M,
                      row.names = ordertest[, 1])
  cat("Scott-Knott Test","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  print(saida)
  return(invisible(saida))
}

#@rdname spScottKnott
#@method spScottKnott GEOrcbd
#@export

# spScottKnott.GEOrcbd <- function(x, sig.level = 0.05){
#   # variaveis
#   dados <- x$data
#   resp <- dados$data
#   trat <- dados$covariate[ ,1]
#   block <- dados$covariate[,2]
#   D <- dados$coords
#
#   # Atributos do modelo geoestatistico
#   covMod <- x$mod
#   nugget <- x$params[3]
#   psill <- x$params[1]
#   phi <- x$params[2]
#
#   y <- c(resp)
#   k <- nlevels(as.factor(trat)) ## numero de tratamentos
#   r <- tapply(y, trat, "length") ## vetor com n de repeticoes
#   n <- sum(r) ## n de observacoes
#
#   # Matriz de incidencia dos tratamentos
#   blk <- length(unique(block)) #numero de blocos
#   X1 <- as.matrix(x$des.mat[ ,1], ncol = 1)
#   X2 <- x$des.mat[ ,-c(1,(k+2):(blk+k+1))]
#
#   ## sigma: a matriz de covariancia espacial
#   sigma <- varcov.spatial(coords = D, cov.model = covMod, nugget = nugget,
#                           kappa = 0.5, cov.pars = c(psill,phi))$varcov
#   V <- sigma/sigma[1,1]
#   i.V <- chol2inv(chol(V))
#
#   #Estimacao considerando tendencia inerente ao modelo dbc
#
#   P <- i.V- i.V%*%X2%*%solve(t(X2)%*%i.V%*%X2)%*%t(X2)%*%i.V
#   theta <- solve(t(D)%*%P%*%D)%*%t(D)%*%P%*%y
#   beta <- solve(t(X2)%*%i.V%*%X2)%*%t(X2)%*%i.V%*%y -
#     solve(t(X2)%*%i.V%*%X2)%*%t(X2)%*%i.V%*%D%*%theta
#
#
#   trt <- c()
#   for (i in 1:k){
#     trt <- as.factor(c(trt, rep(i, r[i])))
#   }
#
#   dados.na <- data.frame(obs = 1:length(y), y = y, trat = trt,
#                          D1 = D[,1], D2 = D[,2])
#   obs <- c()
#   for (i in 1:k){
#     obs <- c(obs, subset(dados.na, (trat == order(beta)[i]))$obs)
#   }
#
#   dados.ord <- dados.na[obs,]
#   beta.ord <- sort(beta)
#   sigma <- varcov.spatial(dados.ord[,4:5], cov.model = covMod, nugget = nugget,
#                           cov.pars = c(psill, phi))$varcov
#   V <- sigma/sigma[1,1]
#   i.V <- chol2inv(chol(V))
#
#   B <- c()
#   for (i in 1:(k-1)){
#     n1 <- sum(r[1:i])
#     n2 <- sum(r[(i+1):k])
#     n <- n1+n2
#     K1 <- matrix(c(rep(1,n1),rep(0,n2)),ncol=1)
#     K2 <- matrix(c(rep(0,n1),rep(1,n2)),ncol=1)
#     K <- cbind(K1,K2)
#     X1 <- matrix(c(rep(1,n)),ncol=1)
#     Y <- c(dados.ord$y)
#     sqd <- t(Y)%*%i.V%*%Y- t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%
#       t(X1)%*%i.V%*%Y
#     sqn <- t(Y)%*%i.V%*%K%*%solve(t(K)%*%i.V%*%K)%*%t(K)%*%i.V%*%Y -
#       t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%t(X1)%*%i.V%*%Y
#     B <- c(B, sqn/sqd)}
#   posmax <- grep(max(B),B)
#   lamb <- as.numeric((pi/(2*(pi-2)))*n*B[posmax])
#   pval <- pchisq(lamb, k/(pi-2), lower.tail = FALSE)
#   groups <- rep(0, times = k)
#   ng <- 1
#   if (pval > sig.level) {
#     groups[1:k] <- 1
#   } else {
#     groups[1:posmax] <- 1
#     groups[(posmax+1):k] <- 2
#     ng <- ng + 1
#   }
#   if ((any(groups!=1))) {
#     continua1 = TRUE
#   } else {
#     continua1 = FALSE
#   }
#   if (continua1 == TRUE) {
#     posI <- 1
#     fim <- FALSE
#     while(posI <= k) {
#       ini <- groups[posI]
#       posF <- max(which(groups == ini))
#       if (((posF - posI) > 0) & (ini > 0)) {
#         B <- c()
#         for (i in 1:(posF-posI)){
#           n1 <- sum(r[posI:(posI+i-1)])
#           n2 <- sum(r[(i+posI):posF])
#           n <- n1+n2
#           K1 <- matrix(c(rep(1,n1),rep(0,n2)),ncol=1)
#           K2 <- matrix(c(rep(0,n1),rep(1,n2)),ncol=1)
#           K <- cbind(K1,K2)
#           X1 <- matrix(c(rep(1,n)),ncol=1)
#           if (posI==1) { g=0 } else { g=sum(r[1:(posI-1)]) }
#           dados.ord1 <- dados.ord[(g+1):(g+n),]
#           Y <- c(dados.ord1$y)
#           sigma <- varcov.spatial(dados.ord1[,4:5],cov.model = covMod, nugget = nugget,
#                                   cov.pars = c(psill,phi))$varcov
#           V <- sigma/sigma[1,1]
#           i.V <- chol2inv(chol(V))
#           sqd <- t(Y)%*%i.V%*%Y- t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%
#             t(X1)%*%i.V%*%Y
#           sqn <- t(Y)%*%i.V%*%K%*%solve(t(K)%*%i.V%*%K)%*%t(K)%*%i.V%*%Y -
#             t(Y)%*%i.V%*%X1%*%solve(t(X1)%*%i.V%*%X1)%*%t(X1)%*%i.V%*%Y
#           B <- c(B, sqn/sqd) }
#         posmax <- grep(max(B),B)
#         lamb <- as.numeric((pi/(2*(pi-2)))*n*B[posmax])
#         pval <- pchisq(lamb, (posF-posI+1)/(pi-2), lower.tail = FALSE)
#         aux <- groups[posI]
#         if (pval > sig.level) {
#           groups[posI:posF] <- - aux
#           posI <- posF + 1
#           if (posI > k) {
#             posI <- 1
#           }
#         }
#         else {
#           groups[(posI+posmax):posF] <- ng + 1
#           ng <- ng + 1
#           posI <- 1
#         }
#       }
#       else posI <- posF + 1
#       #if (posI >= k) fim <- TRUE
#       #if (fim == TRUE) break
#     }
#     groups <- abs(groups)
#   }
#   Teste <- data.frame(trat=levels(factor(trat))[order(beta)],
#                       medias=beta.ord, grupos=groups)
#   k <- nrow(Teste)
#   ordertest <- Teste[order(Teste[, 2], decreasing = TRUE),]
#   M <- rep("", k)
#   M[1] <- "a"
#   j <- 2
#   for (i in 2:k) {
#     if (ordertest[i, 3] == ordertest[i-1, 3]) {
#       M[i] <- M[i-1] }
#     else { M[i] <- letters[j]
#     j <- j + 1 }
#   }
#   saida <- data.frame(means = unname(round(ordertest[, 2],3)), groups = M,
#                       row.names = ordertest[, 1])
#   cat("Scott-Knott Test","\n")
#   cat("\n")
#   cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
#   cat("\n")
#   print(saida)
#   return(invisible(saida))
# }
