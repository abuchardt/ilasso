#' Identify interactions with lasso for multivariate outcome
#'
#' Identifies main effects and pairwise interactions in a multivariate linear or logistic regression model via hierarchical lasso penalisation. Can deal with all shapes of data, including very large sparse data matrices.
#'
#' Step 1 assumes a multivariate pure main effects model and returns the main effects selected by lasso corresponding to the set of non-zero estimates.
#'
#' Step 2 assumes a multivariate pairwise interaction model subject to a hierarchical restriction and returns the main effects and interactions selected by the MSGLasso corresponding to the set of non-zero estimates.
#'
#' @param X input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.
#' @param Y response matrix, of dimension nobs x nouts. Quantitative for family="gaussian".
#' @param step character. Which step should be run. For \code{step = "both"} or \code{step = "second"}, a value for \code{hierarchy} must also be provided. For \code{step = "second"}, values for \code{X1} must also be provided.
#' @param hierarchy character. The hierarchy of the model.
#' @param X1 numeric vector. This is for the \code{step = "second"}, and allows the user to provide the main effects for which interactions are included (or not, depending on the \code{hierarchy}).
#' @param Y1 numeric vector. This is for the \code{step = "second"}, and allows the user to provide the responses for which interactions are included (or not, depending on the \code{hierarchy}).
#' @param standardize logical. Should X be standardised prior to fitting the model. Defaults to \code{TRUE}.
#' @param standardize.response logical. Should Y be standardised prior to fitting the model. Defaults to \code{TRUE}.
#' @param fold a positive integer for the corss validation fold. Default is \code{fold = 5}.
#' @param seed a numeric scaler, specifying the seed of the random number generator in R for generating cross validation subset for each fold. Default=1.
#' @param reruns numeric. The number of reruns of the procedure, the default is 10.
#' @param ... Other arguments passed to the MSGLasso.cv function.
#'
#' @importFrom stats coef sd
#' @importFrom utils combn
#' @importFrom wrapr match_order
#' @importFrom dplyr arrange
#' @importFrom MASS mvrnorm
#'
#' @import MSGLasso
#'
#' @return An object with S3 class "milasso". Selected main effects and interactions. \item{X1 }{For \code{step = "both"} or \code{step = "first"} only. A vector consisting of the indecies of the main effects selected in step 1.} \item{X1 }{For \code{step = "both"} or \code{step = "second"} only. A vector consisting of the indecies of the main effects selected in step 2.} \item{interactions2 }{For \code{step = "both"} or \code{step = "second"} only. A matrix consisting of the indecies of the main effects (columns) corresponding to the pairwise interactions (rows) selected in step 2.}
#'
#'
#' @seealso \code{\link{MSGLasso}} and \code{\link{MSGLasso.cv}}
#'
#' @examples
#' # Gaussian
#' p <- 20
#' q <- 10
#' n <- 150
#' # Working variables
#' nz <- dplyr::arrange(expand.grid(a=1:p,b=1:p),a)
#' # Features
#' set.seed(1)
#' X <- matrix(sample(0:2, n*p, replace=TRUE), nrow=n, ncol=p)
#' colnames(X) <- paste0("X[, ", 1:p, "]")
#' # Interactions
#' ind <- t(combn(ncol(X),2))
#' out <- apply(ind, 1, function(x) X[,x[1]] * X[,x[2]])
#' colnames(out) <- apply(ind, 1, function(x) paste0("X[, ",x[1],"]:X[, ",x[2],"]"))
#' # Design
#' design <- cbind(X, out)
#' # X1 and X2 affect Y1
#' B <- matrix(0, nrow = p, ncol = q)
#' B[1, 1] <- 1 # X1 påvirker Y1
#' B[2, 1] <- 1 # X2 påvirker Y1
#' BB <- matrix(0, nrow = choose(p, 2), ncol = q)
#' BB[1,1] <- 1 # VV mellem X1 og X2 på Y1
#' # Parameters
#' sigmaE <- c(1,rep(5, q-1))
#' SigmaE <- diag(x = sigmaE, nrow = q, ncol = q)
#' Y <- design %*%  rbind(B, BB) + MASS::mvrnorm(n = n, rep(0, q), SigmaE)
#' # Step 1
#' fit1 <- milasso(X = X, Y = Y, step = "first")
#' me1x <- fit1$X1
#' me1y <- fit1$Y1
#' print(fit1)
#' # Step 2
#' fit2 <- milasso(X = X, Y = Y, step = "second", hierarchy = "strong", X1 = me1x, Y1 = me1y)
#' print(fit2)
#'
#' @export milasso
#'
milasso <- function(X, Y, step = c("both", "first", "second"),
                   hierarchy = c("strong", "weak"),
                   preserve = TRUE,
                   X1, Y1,
                   #family = c("gaussian"), #' @param family Response type (see above).
                   standardize = TRUE, standardize.response = TRUE,
                   fold = 5,
                   seed = NULL,
                   reruns = 1, ...) {

  this.call <- match.call()

  #family <- match.arg(family)
  step <- match.arg(step)
  hierarchy <- match.arg(hierarchy)

  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if(step == "second") {
    X1 <- as.numeric(X1)
    Y1 <- as.numeric(Y1)
    if(hierarchy == "weak" && length(X1) < 1) {
      stop("Number of elements in X1 is less than one")
    }
    if(hierarchy == "strong" && length(X1) < 2) {
      stop("Number of elements in X1 is less than two")
    }
  }

  if(step %in% c("both", "first")) {
    step1 <- .milasso1(X = X, Y = Y, standardize = standardize,
                       fold = fold, seed = seed,
                       reruns = reruns, ...)
    step2 <- NULL
  }

  if (step %in% c("both", "second")) {
    if(step == "both") {
      if(length(step1$X1) < 1) {
        stop("No main effects where selected in step 1")
      } else if (length(step1$X1) == 1 && hierarchy == "strong") {
        stop("Strong hierarchy needs two or more main effects from step 1; only one was selected")
      }
      X1 <- step1$X1
      Y1 <- step1$Y1
      coef1 <- step1$coef1
    } else {
      step1 <- NULL
    }
    if(missing(X1) && step == "second") {
      stop("Step 2 has no selected main effects from step 1")
    }
    step2 <- .milasso2(X = X, Y = Y,
                       hierarchy = hierarchy, preserve = preserve,
                       X1 = X1, Y1 = Y1, coef1 = coef1,
                       standardize = standardize,
                       fold = fold, seed = seed,
                       reruns = reruns, ...) #.milasso2(X, Y, hierarchy, X1, Y1, standardize, reruns)
  }

  if(step == "first"){

    fit <- list(call = this.call,
                hierarchy = hierarchy,
                coef1 = step1$coef1,
                X1 = step1$X1,
                Y1 = step1$Y1)

    class(fit) <- "milasso"
    fit

  } else if (step == "second") {


    fit <- list(call = this.call,
                hierarchy = hierarchy,
                coef2 = step2$coef2,
                ME2 = step2$ME2,
                I2 = step2$I2,
                X2 = step2$X2,
                Y2 = step2$Y2,
                iX2 = step2$iX2,
                iY2 = step2$iY2)

    class(fit) <- "milasso"
    fit

  } else if (step == "both") {

    fit <- list(call = this.call,
                hierarchy = hierarchy,
                coef1 = step1$coef1,
                coef2 = step2$coef2,
                ME2 = step2$ME2,
                I2 = step2$I2,
                X1 = step1$X1,
                Y1 = step1$Y1,
                X2 = step2$X2,
                Y2 = step2$Y2,
                iX2 = step2$iX2,
                iY2 = step2$iY2)

    class(fit) <- "milasso"
    fit

  }

}
#'
#' Internal milasso function
#'
#' @inheritParams milasso
#'
.milasso1 <- function(X, Y, standardize, fold, seed, reruns, ...) {
  if (reruns < 1) {
    warning("reruns<1; set to 1")
    reruns = 1
  }
  reruns <- as.double(reruns)

  np <- dim(X)
  if (is.null(np) | (np[2] <= 1)) {
    stop("X should be a matrix with 2 or more columns")
  }
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  dimy <- dim(Y)
  nrowy <- ifelse(is.null(dimy), length(Y), dimy[1])
  nouts <- as.integer(dimy[2])
  if (nrowy != nobs) {
    stop(paste("number of observations in Y (", nrowy, ") not equal to the number of rows of X (",
               nobs, ")", sep = ""))
  }
  vnames <- colnames(X)
  if (is.null(vnames)) {
    vnames <- paste("V", seq(nvars), sep = "")
  }

  # Lasso regression - wothout repetitions
  if (reruns == 1) {

    # Fit an MSGLasso
    G <- nvars #a positive interger indicating number of predictor groups
    R <- nouts #a positive interger indicating number of response groups
    gmax <- nouts #a positive interger indicating the max number of different groups a single variable (either a predictor or response variable) belongs to.
    cmax <- nvars #a positive interger indicating the max number of variables a single group (either a predictor or response group) contains.
    GarrStarts <- 1:nvars #1 #a vector of starting coordinates for the predictor groups.
    GarrEnds <- 1:nvars #p #a vector of ending coordinates for the predictor groups.
    RarrStarts <- 1:nouts #1 #a vector of starting coordinates for the response groups.
    RarrEnds <- 1:nouts #q #a vector of ending coordinates for the response groups.

    tmp <- try(MSGLasso::FindingPQGrps(nvars, nouts, G, R, gmax,
                                       GarrStarts, GarrEnds, RarrStarts, RarrEnds),
               silent=TRUE)
    if(inherits(tmp, "try-error")) stop("MSGLasso error")

    PQgrps <- tmp$PQgrps

    #tmp1 <- MSGLasso::Cal_grpWTs(nvars, nouts, G, R, gmax, PQgrps)
    #grpWTs <- tmp1$grpWTs
    grpWTs <- matrix(1, nrow = nvars, ncol = nouts)

    tmp2 <- try(MSGLasso::FindingGRGrps(nvars, nouts, G, R, cmax,
                                        GarrStarts, GarrEnds, RarrStarts, RarrEnds),
                silent=TRUE)
    if(inherits(tmp, "try-error")) stop("MSGLasso error")

    GRgrps <- tmp2$GRgrps

    Pen.L <- matrix(rep(1,nvars*nouts), nvars, nouts, byrow=TRUE)
    ## Due to MSGLasso bug (?)
    #Pen_L <- Pen.L
    ## End bug fix
    Pen.G <- matrix(rep(1,G*R),G,R, byrow=TRUE)
    grp_Norm0 <- matrix(rep(0, G*R), nrow=G, byrow=TRUE)

    # No CV
    #MSGLassolam1 <- 1.6
    #MSGLassolamG <- 0.26
    #MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)
    #msglassofit <- MSGLasso::MSGLasso(X, Y, grpWTs, Pen.L, Pen.G, PQgrps, GRgrps, grp_Norm0,
    #                                  MSGLassolam1, MSGLassolamG.m, Beta0=NULL)

    # With CV
    lam1.v <- seq(1.0, 1.5, length=6)
    lamG.v <- seq(0.19, 0.25, length=7)

    msglassofitcv <- MSGLasso::MSGLasso.cv(#.MSGLasso.cv(#
      X, Y, grpWTs, Pen.L, Pen.G,
      PQgrps, GRgrps,
      lam1.v, lamG.v, fold, seed = seed)
    MSGLassolam1 <- msglassofitcv$lams.c[which.min(as.vector(msglassofitcv$rss.cv))][[1]]$lam1
    MSGLassolamG  <- msglassofitcv$lams.c[which.min(as.vector(msglassofitcv$rss.cv))][[1]]$lam3
    MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)

    msglassofit <- MSGLasso::MSGLasso(#.MSGLasso(#
      X, Y, grpWTs, Pen.L, Pen.G, PQgrps, GRgrps,
      grp_Norm0, MSGLassolam1, MSGLassolamG.m, Beta0=NULL)

    #which(msglassofit$Beta != 0, arr.ind = TRUE)

    # Extract coefficients
    vcoef <- msglassofit$Beta

    # Lasso regression - with repetitions
  } else {

    stop("only reruns = 1 implemented")

  }

  # Step-1-selected main effects (Non-zero coefficients excluding intercept)
  meXY <- which(vcoef != 0, arr.ind = TRUE)
  X1 <- meXY[,1]
  Y1 <- meXY[,2]

  list(coef1 = vcoef,
       X1 = X1, Y1 = Y1)
}
#'
#' Internal milasso function
#'
#' @inheritParams milasso
#'
.milasso2 <- function(X, Y, hierarchy, preserve, X1, Y1, coef1, standardize, fold, seed, reruns, ...) {

  if (reruns < 1) {
    warning("reruns<1; set to 1")
    reruns = 1
  }
  reruns <- as.double(reruns)
  np <- dim(X)
  if (is.null(np) | (np[2] <= 1))
    stop("X should be a matrix with 2 or more columns")
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  dimy <- dim(Y)
  nrowy <- ifelse(is.null(dimy), length(Y), dimy[1])
  nouts <- as.integer(dimy[2])
  if (nrowy != nobs)
    stop(paste("number of observations in Y (", nrowy, ") not equal to the number of rows of X (",
               nobs, ")", sep = ""))
  colnames(X) <- paste0("X[, ", 1:nvars, "]")

  # Features:
  # Include all main effects
  # There should be no penalty on main effects selected in step 1
  # Include only interactions between main effects selected in step 1 and another main effect
  if(hierarchy == "weak"){

    # Interactions
    ind <- expand.grid(1:ncol(X), X1)
    out <- apply(ind, 1, function(i) X[,i[1]] * X[,i[2]])
    if (requireNamespace("Matrix", quietly = TRUE)) {
      out <- Matrix::Matrix(out, sparse = TRUE)
    }
    colnames(out) <- apply(ind, 1, function(x) paste0("X[, ",x[1],"]:X[, ",x[2],"]"))

    newX <- cbind(X, out)

    # Include only interactions between step-1-selected main effects
    # Include all main effects
    # To honour hierarchy: no penalty on step-1-selected main effects.
  } else if(hierarchy == "strong") {

    #mesx <- X1 #which(msglassofit$Beta != 0, arr.ind = TRUE)[,1]
    #mesy <- Y1 #which(msglassofit$Beta != 0, arr.ind = TRUE)[,2]

    # Interactions to include
    ind <- t(combn(unique(X1),2))
    out <- apply(ind, 1, function(i) X[,i[1]] * X[,i[2]])
    colnames(out) <- apply(ind, 1, function(x) paste0("X[, ",x[1],"]:X[, ",x[2],"]"))

    newX <- cbind(X, out)
  }

  P <- ncol(newX)
  G <- P #a positive interger indicating number of predictor groups
  R <- nouts #a positive interger indicating number of response groups
  gmax <- nouts #a positive interger indicating the max number of different groups a single variable (either a predictor or response variable) belongs to.
  cmax <- P #a positive interger indicating the max number of variables a single group (either a predictor or response group) contains.
  GarrStarts <- 1:P #a vector of starting coordinates for the predictor groups.
  GarrEnds <- 1:P   #a vector of ending coordinates for the predictor groups.
  RarrStarts <- 1:nouts #a vector of starting coordinates for the response groups.
  RarrEnds <- 1:nouts #a vector of ending coordinates for the response groups.

  tmp <- MSGLasso::FindingPQGrps(P, nouts, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  PQgrps <- tmp$PQgrps

  #tmp1 <- MSGLasso::Cal_grpWTs(P, nouts, G, R, gmax, PQgrps)
  #grpWTs <- tmp1$grpWTs
  grpWTs <- matrix(1, nrow = P, ncol = nouts)
  # Penalty factor that multiplies lambda to allow
  # no shrinkage of selected main effects
  if (preserve) {
    grpWTs[X1, ] <- 0
  }
  if (!preserve) {
    grpWTs[abs(coef1) > 1.5e-8] <- 0
  }

  tmp2 <- MSGLasso::FindingGRGrps(P, nouts, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
  GRgrps <- tmp2$GRgrps

  Pen.L <- matrix(rep(1,P*nouts),P,nouts, byrow=TRUE)
  if (preserve) {
    Pen.L[X1, ] <- 0
  }
  if (!preserve) {
    Pen.L[abs(coef1) > 1.5e-8] <- 0
  }
  ## Due to MSGLasso bug (?)
  #Pen_L <- Pen.L
  ## End bug fix
  Pen.G <- matrix(rep(1,G*R),G,R, byrow=TRUE)
  if (preserve) {
    Pen.G[X1, ] <- 0
  }
  if (!preserve) {
    #Pen.G[X1, Y1] <- 0
    Pen.G[abs(coef1) > 1.5e-8] <- 0
  }

  grp_Norm0 <- matrix(rep(0, G*R), nrow=G, byrow=TRUE)

  # No CV
  #MSGLassolam1 <- 1.6
  #MSGLassolamG <- 0.26
  #MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)

  # With CV
  lam1.v <- seq(1.0, 1.5, length=6)
  lamG.v <- seq(0.19, 0.25, length=7)


  # MSGLasso regression
  if(reruns==1){

    # No CV
    # msglassofit2 <- try(MSGLasso::MSGLasso(newX, Y, grpWTs, Pen.L, Pen.G,
    #                                        PQgrps, GRgrps, grp_Norm0,
    #                                        MSGLassolam1, MSGLassolamG.m, Beta0=NULL),
    #                     silent=TRUE)

    # With CV
    msglassofitcv2 <- try(MSGLasso::MSGLasso.cv(#.MSGLasso.cv(#
      newX, Y, grpWTs, Pen.L, Pen.G, PQgrps, GRgrps,
      lam1.v, lamG.v, fold, seed = seed), silent=TRUE)

    if(inherits(msglassofitcv2, "try-error")) stop("MSGLasso error")

    MSGLassolam1 <- msglassofitcv2$lams.c[which.min(as.vector(msglassofitcv2$rss.cv))][[1]]$lam1
    MSGLassolamG  <- msglassofitcv2$lams.c[which.min(as.vector(msglassofitcv2$rss.cv))][[1]]$lam3
    MSGLassolamG.m <- matrix(rep(MSGLassolamG, G*R),G,R,byrow=TRUE)

    msglassofit2 <- try(MSGLasso::MSGLasso(#.MSGLasso(#
      newX, Y, grpWTs, Pen.L, Pen.G,
      PQgrps, GRgrps, grp_Norm0,
      MSGLassolam1, MSGLassolamG.m, Beta0=NULL),
      silent=TRUE)

    if(inherits(msglassofit2, "try-error")) stop("MSGLasso error")

    # Zero and non-zero coefficients
    #c2c <- coef(msglassofit2, s = "lambda.min")
    c2c <- msglassofit2$Beta


  } else {

    stop("only reruns = 1 implemented")

  }

  # Non-zero coefficients (excluding intercept)
  nzc <- which(c2c != 0, arr.ind = TRUE)

  # Step-2-selected main effects(Non-zero coefficients excluding intercept)
  X2 <- nzc[nzc[,1] <= ncol(X),1]
  Y2 <- nzc[nzc[,1] <= ncol(X),2]
  ME2 <- matrix(0, nrow = p, ncol = nouts)
  ME2[abs(c2c[1:p, ]) > 1.5e-8] <- 1

  # Step-2-selected interactions
  iX2 <- ind[nzc[nzc[,1] > ncol(X),1]-ncol(X),1]
  names(iX2) <- colnames(out)[iX2]
  iY2 <- nzc[nzc[,1] > ncol(X),2]
  I2 <- matrix(0, nrow = ncol(out), ncol = nouts)
  I2[abs(c2c[(p+1):(p+ncol(out)),]) > 1.5e-8] <- 1
  rownames(I2) <- colnames(out)

  # Return

  list(coef2 = c2c, newX = newX,
       ME2, I2,
       X2 = X2, Y2 = Y2,
       iX2 = iX2, iY2 = iY2)

}

