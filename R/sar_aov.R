#' @name sar_aov
#'
#' @title Using a SAR model to handle spatial dependence in a aov model.
#' @description Fit a completely randomized design when the experimental units have some degree of
#' spatial dependence using a Spatial Lag Model (SAR).
#' @usage sar_aov(formula, coord, data, seq.radius = NULL)
#'
#' @param formula A formula specifying the model.
#' @param coord A matrix or data.frame of point coordinates or a SpatialPoints Object.
#' @param data A data frame in which the variables specified in the formula will be found.
#' @param seq.radius A complex vector containing a radii sequence used to set the neighborhood pattern.
#' The default sequence has ten numbers from 0 to half of the maximum distance between the samples.
#'
#' @return \code{sar_aov} returns an object of \code{\link[base]{class}} "SARaov".
#' The functions summary and anova are used to obtain and print a summary and analysis of variance
#' table of the results.
#' An object of class "SARaov" is a list containing the following components:
#'
#' \item{DF}{degrees of freedom of rho, treatments, residual and total.}
#' \item{SS}{sum of squares of rho, treatments and residual.}
#' \item{rho}{the autoregressive parameter.}
#' \item{Par}{data.frame with the radius tested and its AIC.}
#' \item{modelAdj}{model of class \code{\link[stats]{aov}} using the adjusted response.}
#' \item{modelstd}{data frame containing the ANOVA table using non-adjusted response.}
#' \item{namey}{response variable name.}
#'
#'
#' @references Long, D. S. "Spatial statistics for analysis of variance of agronomic field trials."
#' Practical handbook of spatial statistics. CRC Press, Boca Raton, FL (1996): 251-278.
#'
#' @examples
#' \dontrun{
#' data("carrancas")
#' resp <- carrancas$DAP16
#' treat <- carrancas$T
#' coord <- cbind(carrancas$X, carrancas$Y)
#' mdl <- aov(resp ~ treat, data = carrancas)
#' cv<-sar_aov(mdl, coord)
#' cv
#'
#' #Summary for class SARanova
#' summary(cv)
#'
#' #Anova for class SARanova
#' anova(cv)
#' }
#'
#' @import spdep
#' @importFrom gtools stars.pval
#' @importFrom car Anova
#' @export sar_aov print.SARaov summary.SARaov anova.SARaov

sar_aov <- function(formula, coord, data, seq.radius) {

  # Defensive programming
  if(!(is.matrix(coord) | class(coord)=="SpatialPoints")) {
    stop("'coord' must be a matrix or SpatialPoints object")
  }

  if(ncol(coord) < 2){
    stop("'coord' must have at least two columns")
  }

  if(missing(data)){
    stop("'data' must be provided")
  }

  if(class(data) != "data.frame"){
    stop("'data' must be a data.frame")
  }

  if(is.data.frame(coord)){
    coord <- as.matrix(coord, ncol = ncol(coord))
  }

  if(missing(seq.radius)){
    max.dist <- max(dist(coord))
    seq.radius <- seq(0, 0.5*max.dist, l = 11)[-1]
  }

  params <- data.frame(radius = 0, rho = 0, AIC = 0)
  anova.list <- list()
  p.radius <- length(seq.radius)
  Y_ajus <- NULL
  form.char <- invisible(capture.output(print(formula, showEnv = FALSE)))
  formula.sar <- as.formula(gsub("as.factor|factor", "", form.char))


  for (i in 1:p.radius) {
    nb <- dnearneigh(coord, 0, seq.radius[i])
    w <- try(nb2mat(nb, style = "W"), silent = TRUE)
    test <- grepl("Error", w)

    # Se caso nao forem encontradas amostras dentro do raio especificado
    k <- 0.1 # incremento
    while(test[1] == TRUE){
      seq.radius <- seq(0, (0.5+k)*max.dist, l = 11)[-1]
      nb <- dnearneigh(coord, 0, seq.radius[i])
      w <- try(nb2mat(nb, style = "W"), silent = TRUE)
      test <- grepl("Error", w)
      k <- k + 0.1
    }

    listw <- nb2listw(nb, glist = NULL, style = "W")

    # SAR model
    SAR <- lagsarlm(formula.sar, data = data, listw = listw,
                    method = "eigen", tol.solve = 1e-15)
    ajuste <- summary(SAR)
    rho <- as.numeric(ajuste["rho"]$rho)
    params[i, ] <- c(raio = seq.radius[i], rho = rho, AIC = AIC(SAR))
  }

  # Separar o nome da variavel resposta na formula para substitui-la por Y_ajus
  form.split <- strsplit(form.char, "~")[[1]]
  new.formula <- as.formula(paste("Y_ajus","~",form.split[2]))
  resp.name <- sub(" ","",form.split[1])
  resp <- data[ ,resp.name]

  # Adjusting the data and constructing the ANOVA table
  best.par <- which.min(params$AIC)
  beta <- mean(resp)
  nb <- dnearneigh(coord, 0, seq.radius[best.par])
  w <- nb2mat(nb, style = "W")
  Y_ajus <- resp - (params[best.par,"rho"] * w%*%resp - params[best.par,"rho"] * beta)
  model.cl <- aov(formula, data = data)
  aov.cl <- anova(model.cl)

  # Fazer a analise com os dados ajustados
  new.data <- data
  new.data[ ,resp.name] <- Y_ajus
  model.adj <- aov(new.formula, data = new.data)
  aov.adj <- anova(model.adj)
  Sqt.nadj <- sum(aov.cl[,2])

  #Degres of freedom
  glerror <- model.cl$df.residual

  #Sum of squares
  sqerror <- sum(model.cl$residuals^2)
  sqtot <- sum(aov.adj[,2])
  sqtotcor <- Sqt.nadj - sqtot

  #Mean Squares
  mserror <- sqerror/glerror

  name.y <- names(attr(model.cl$terms,"dataClasses")[1])

  outpt <- list(DF = glerror,
                SS = c(sqerror, sqtotcor),
                MS = c(mserror),
                rho = params[best.par,"rho"], Par = params,
                modelAdj = model.adj, modelstd = aov.cl, namey = name.y)

  class(outpt)<-c("SARaov",class(aov.adj))
  return(outpt)
}

# Print method for this class

print.SARaov <- function(x) {
  cat("Response: ", x$namey, "\n")
  rse <- sqrt(x$MS)
  cat("\n")
  cat("Residual standard error:",rse)
  cat("\n")
  cat("Spatial autoregressive parameter:", x$rho,"\n")
  cat("Samples considered neighbor within a",x$Par[which.min(x$Par[,3]),1],"units radius")
}

# Summary method for this class

summary.SARaov <- function(x) {
  cat("      Summary of SARaov","\n","\n")
  cat("Parameters tested:","\n")
  print(x$Par)
  cat("\n")
  cat("Selected parameters:","\n")
  print(x$Par[which.min(x$Par[,3]),])
}

# Anova method for this class

anova.SARaov <- function(x, type = c("II","III", 2, 3), compare = FALSE) {
  type <- as.character(type)
  type <- match.arg(type)

  if(missing(type)){
    type = "II"
  }

  if(is.logical(compare) == FALSE){
    warning("'compare' must be logical. Assuming compare == FALSE")
    compare = FALSE
  }

  cat("Analysis of Variance With Spatially Correlated Errors","\n")
  cat("\n")
  print(Anova(x$modelAdj, type = type))


  if(compare){
    cat("\n", "\n")
    cat("---------------------------------------------------------------","\n")
    cat("Standard Analysis of Variance", "\n")
    cat("---------------------------------------------------------------")
    cat("\n")
    print(x$modelstd)
  }

}


#exportar a funcao stars.pval do pacote gtools

