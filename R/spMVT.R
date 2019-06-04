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
#' @return
#' a data frame containing the means and its group.
#'
#' @references
#' NOGUEIRA, C. H. Testes para comparações múltiplas de
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
  let <- cld(comp,decreasing = TRUE, level = sig.level)
  mbg <- tapply(x$modelAdj$model[,1], x$modelAdj$model[,2], mean)
  result <- data.frame(means = mbg[order(mbg,decreasing = TRUE)] ,
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
  if(trend == "cte"){
    X2 <- x$des.mat[ ,-1]
  }else{
    X2 <- x$des.mat[ ,2:(k+1)]
  }

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

  print.tuk <- data.frame(means = round(unname(sort(media, decreasing = TRUE)),3),
                          groups = unname(letras),
                          row.names =  unique(trat)[order(media,decreasing = TRUE)])

  # names(val.p) <- rownames(C)
  # tukey <- multcompLetters(val.p, rev=TRUE)
  #
  # if (sum(val.p <= sig.level) == 0){
  #   letras=tukey$Letters
  # }
  # else{
  #   letras = tukey$monospacedLetters
  # }
  #
  # print.tuk <- data.frame(Treatment = order(media, decreasing=TRUE),
  #                      Mean = round(sort(media, decreasing=TRUE),3),
  #                      Group = letras[order(media, decreasing=TRUE)], row.names = NULL)

  cat("Test based on multivariate t-student distribution","\n")
  cat("\n")
  cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
  cat("\n")
  print(print.tuk)
  return(invisible(print.tuk))
}


#@rdname spMVT
#@method spMVT GEOrcbd
#@export

# spMVT.GEOrcbd<- function(x, sig.level = 0.05){
#
#   # variaveis
#   dados <- x$data
#   resp <- dados$data
#   trat <- dados$covariate[ ,1]
#   block <- dados$covariate[ ,2]
#   D <- dados$coords
#
#   # Atributos do modelo geoestatistico
#   covMod <- x$mod
#   nugget <- x$params[3]
#   psill <- x$params[1]
#   phi <- x$params[2]
#
#   Y <- c(resp)
#   k <- nlevels(as.factor(trat)) ## numero de tratamentos
#   r <- tapply(Y, trat, "length") ## vetor com n de repeticoes
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
#   V<-sigma
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
#   media.aj <- media + as.numeric(t(D.m)%*%theta)
#   ## variancia das medias
#   A <- D%*%solve(t(D)%*%P%*%D)%*%t(D)
#   var.mu <- (solve(t(X2)%*%i.sig%*%X2) + solve(t(X2)%*%i.sig%*%X2)%*%
#                t(X2)%*%i.sig%*%A%*%i.sig%*%X2%*%solve(t(X2)%*%i.sig%*%X2))
#   C <- contr.tuk(media) ## matriz dos contrastes de tukey
#   W <- C%*%var.mu%*%t(C)
#
#
#   #  Fase 2 - Aplicacao do teste---------------------------------------------------------
#
#   S.i <- diag(1/sqrt(diag(W)))
#   R <- S.i%*%W%*%S.i;R ## matriz de correlacao dos contrastes
#   est <- C%*%media ## estimativa dos contrastes
#   ep <- sqrt(diag(W)) ## erro padrao dos contrastes
#   tc <- est/ep ## estatistica t
#   glres <- sum(r) - k #ou sum(r) - k - s
#   p <- nrow(C) ## numero de contrastes
#   r <- with(dados,tapply(Y, trat, "length")) ## n de repeticoes
#   val.p <- as.vector(0)
#   for (i in 1:p){
#     val.p[i] <- 1-pmvt(-rep(abs(tc[i]), p), rep(abs(tc[i]), p),
#                        corr = R, df = glres)
#   }
#   teste <- data.frame(Estimativa=round(est,4), Erro_Padrao = round(ep,4),
#                       t_calc=round(tc,3), Valor_p=round(val.p,4))
#
#   #Visualizacao
#   Q <- matrix(1, ncol = k, nrow = k)
#   p <- val.p
#   idx <-0
#   for(i in 1:(k-1)){
#     for(j in (i+1):k){
#       idx <- idx + 1
#       Q[i,j] <- val.p[idx]
#       Q[j,i] <- val.p[idx]
#     }
#   }
#   print.tuk <- orderPvalue(unique(trat), round(as.numeric(media),3),
#                            sig.level, Q, console=FALSE)
#
#   # names(val.p) <- rownames(C)
#   # tukey <- multcompLetters(val.p, rev=TRUE)
#   #
#   # if (sum(val.p <= sig.level) == 0){
#   #   letras=tukey$Letters
#   # }
#   # else{
#   #   letras = tukey$monospacedLetters
#   # }
#   #
#   # print.tuk <- data.frame(Treatment = order(media, decreasing=TRUE),
#   #                      Mean = round(sort(media, decreasing=TRUE),3),
#   #                      Group = letras[order(media, decreasing=TRUE)], row.names = NULL)
#
#   cat("Test based on multivariate t-student distribution","\n")
#   cat("\n")
#   cat("Treatments with the same letter are not significantly different at a", sig.level * 100,"%" ,"significance level.", "\n")
#   cat("\n")
#   print(print.tuk)
#   return(invisible(print.tuk))
# }


