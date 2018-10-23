sar_crd <- function(resp, treat, coords, radius.min, radius.max, by) {

  # Defensive programming
  if(!(is.vector(resp) | is.numeric(resp))) {
    stop("'resp' must be a vector or numeric")
  }

  if(!(is.vector(treat) | is.numeric(treat))) {
    stop("'treat' must be a vector or numeric")
  }

  if(!(is.matrix(coords) | class(coords)=="SpatialPoints")) {
    stop("'coords' must be a matrix or SpatialPoints object")
  }

  if(!(is.numeric(radius.min))) {
    stop("'radius.min' must be a numeric")
  }

  if(!(is.numeric(radius.max))) {
    stop("'radius.max' must be a numeric")
  }

  if(!(is.numeric(by))) {
    stop("'radius.by' must be a numeric")
  }

  if(length(radius.min)!=1) {
    stop("'radius.min' must be of length 1")
  }

  if(length(radius.max)!=1) {
    stop("'radius.max' must be of length 1")
  }

  if(length(by)!=1) {
    stop("'by' must be of length 1")
  }


  seq.radius <- seq(radius.min, radius.max, by=by)
  params <- data.frame(radius = 0, rho = 0, AIC = 0)
  anova.list <- list()
  n <- length(resp)
  p.radius <- length(seq.radius)
  Y_ajus <- NULL

  for (i in 1:p.radius) {
    nb <- dnearneigh(coords, 0, seq.radius[i])
    w <- nb2mat(nb, style = "W")
    listw <- nb2listw(nb, glist = NULL, style = "W")

    # SAR model
    SAR <- lagsarlm(resp ~ treat, listw = listw,
                    method = "eigen", tol.solve = 1e-15)
    ajuste <- summary(SAR)
    rho <- as.numeric(ajuste["rho"]$rho)
    params[i, ] <- c(raio = seq.radius[i], rho = rho, AIC = AIC(SAR))
  }

  # Adjusting the data and constructing the ANOVA table
  best.par <- which.min(params$AIC)
  beta <- mean(resp)
  nb <- dnearneigh(coords, 0, seq.radius[best.par])
  w <- nb2mat(nb, style = "W")
  Y_ajus <- resp - (params[best.par,"rho"] * w%*%resp - params[best.par,"rho"] * beta)
  treat <- factor(treat)
  aov.cl <- anova(aov(resp ~ treat))
  model <- aov(Y_ajus ~ treat)
  aov.adj <- anova(model)
  Sqt.nadj <- sum(aov.cl[,2])

  #Degres of freedom
  glrho <- 1
  gltrat <- aov.adj[1][[1]][1]
  glerror <- aov.adj[1][[1]][2]-1
  gltot <- sum(glrho,gltrat,glerror)

  #Sum of squares
  sqrho <- Sqt.nadj - sum(aov.adj[2][[1]])
  sqtrat <- aov.adj[2][[1]][1]
  sqerror <- aov.adj[2][[1]][2]
  sqtot <- sum(sqrho, sqtrat, sqerror)

  #Mean Squares
  msrho <- sqrho/glrho
  mstrat <- sqtrat/gltrat
  mserror <- sqerror/glerror

  #F statistics
  ftrat <- mstrat/mserror
  pvalue <- pf(ftrat,gltrat,glerror,lower.tail = FALSE)

  outpt <- list(DF = c(glrho, gltrat, glerror, gltot),
                SS = c(sqrho, sqtrat, sqerror),
                MS = c(msrho, mstrat, mserror),
                Fc = c(ftrat),
                p.value = c(pvalue),rho = params[best.par,"rho"],Par = params,
                y_orig = resp, treat = treat , model = model)
  class(outpt)<-c("SARanova","SARcrd",class(aov.adj))
  return(outpt)
}
