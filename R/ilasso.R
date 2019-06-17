#' Identify interactions with lasso
#'
#' Identifies main effects and pairwise interactions in a linear or logistic regression model via hierarchical lasso penalisation. Can deal with all shapes of data, including very large sparse data matrices.
#'
#' Step 1 assumes the the pure main effects model and returns the main effects selected by lasso corresponding to the set of non-zero estimates.
#'
#' Step 2 assumes a pairwise interaction model subject to a hierarchical restriction and returns the main effects and interactions selected by the adaptive lasso corresponding to the set of non-zero estimates.
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.
#' @param y response variable. Quantitative for family="gaussian". For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For "binomial" if y is presented as a vector, it will be coerced into a factor.
#' @param step character. Which step should be run. For \code{step = "both"} or \code{step = "second"}, a value for \code{hierarchy} must also be provided. For \code{step = "second"}, values for \code{maineffects1} must also be provided.
#' @param hierarchy character. The hierarchy of the model.
#' @param maineffects1 numeric vector. This is for the \code{step = "second"}, and allows the user to provide the main effects for which interactions are included (or not, depending on the \code{hierarchy}).
#' @param family Response type (see above).
#' @param standardize logical. Should x be standardised prior to fitting the model. Defaults to \code{TRUE}.
#' @param reruns numeric. The number of reruns of the procedure, the default is 10.
#' @param ... Other arguments passed to the cv.glmnet function.
#'
#' @return An object with S3 class "ilasso". Selected main effects and interactions. \item{maineffects1 }{For \code{step = "both"} or \code{step = "first"} only. A vector consisting of the indecies of the main effects selected in step 1.} \item{maineffects1 }{For \code{step = "both"} or \code{step = "second"} only. A vector consisting of the indecies of the main effects selected in step 2.} \item{interactions2 }{For \code{step = "both"} or \code{step = "second"} only. A matrix consisting of the indecies of the main effects (columns) corresponding to the pairwise interactions (rows) selected in step 2.}
#'
#'
#' @seealso \code{\link{glmnet}} and \code{\link{cv.glmnet}}
#'
#' @examples
#' # Gaussian
#' p <- 200
#' n <- 100
#' nz <- dplyr::arrange(expand.grid(a = 1:p, b = 1:p), a)
#' X <- matrix(sample(c(TRUE, FALSE), n * p, replace = TRUE), nrow = n, ncol = p)
#' ind <- t(combn(ncol(X),2))
#' out <- apply(ind, 1, function(x) X[,x[1]] * X[,x[2]])
#' XX <- cbind(X, out)
#' beta <- c(3,3, rep(0, ncol(X)-2))
#' beta12 <- c(3, rep(0, choose(p, 2)-1))
#' y <- XX %*% c(beta,beta12) + rnorm(n)
#' # Step 1
#' fit1 <- ilasso(x = X, y = y, step = "first")
#' me1 <- fit1$maineffects1
#' print(fit1)
#' # Step 2
#' fit2 <- ilasso(x = X, y = y, step = "second", maineffects1 = me1)
#' print(fit2)
#'
#' @export ilasso
#'
ilasso <- function(x, y, step = c("both", "first", "second"),
                           hierarchy = c("strong", "weak"), maineffects1,
                           family = c("gaussian", "binomial"), standardize = TRUE, reruns = 10, ...) {
  family <- match.arg(family)
  step <- match.arg(step)
  hierarchy <- match.arg(hierarchy)

  x <- as.matrix(x)
  y <- as.numeric(y)

  if(step == "second") {
    maineffects1 <- as.numeric(maineffects1)
    if(hierarchy == "weak" && length(maineffects1) < 1) {
      stop("Number of elements in maineffects1 is less than one")
    }
    if(hierarchy == "strong" && length(maineffects1) < 2) {
      stop("Number of elements in maineffects1 is less than two")
    }
  }

  if(step %in% c("both", "first")) {
    step1 <- .ilasso1(x, y, family, standardize, reruns)
    step2 <- NULL
  }

  if (step %in% c("both", "second")) {
    if(step == "both") {
      if(length(step1$maineffects1) < 1) {
        stop("No main effects where selected in step 1")
      } else if (length(step1$maineffects1) == 1 && hierachy == "strong") {
        stop("Strong hierarchy needs two or more main effects from step 1; only one was selected")
      }
        maineffects1 <- step1$maineffects1
    } else {
      step1 <- NULL
    }
    if(missing(maineffects1) && step == "second") {
      stop("Step 2 has no selected main effects from step 1")
    }
    step2 <- .ilasso2(x, y, hierarchy, maineffects1, family, standardize, reruns)
  }

  if(step == "first"){
    structure(list(maineffects1 = step1$maineffects1), class = "ilasso")
  } else if (step == "second") {
    structure(list(maineffects2 = step2$maineffects2,
                   interactions2 = step2$interactions2), class = "ilasso")
  } else if (step == "both") {
    structure(list(maineffects1 = step1$maineffects1,
                   maineffects2 = step2$maineffects2,
                   interactions2 = step2$interactions2), class = "ilasso")
  }

}
#'
#' Internal ilasso function
#'
#' @inheritParams ilasso
#'
.ilasso1 <- function(x, y, family, standardize, reruns, ...) {
  if (reruns < 1) {
    warning("reruns<1; set to 1")
    reruns = 1
  }
  reruns <- as.double(reruns)
  #this.call <- match.call()
  y <- drop(y)
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  if (nrowy != nobs)
    stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (",
               nobs, ")", sep = ""))
  vnames <- colnames(x)
  if (is.null(vnames))
    vnames <- paste("V", seq(nvars), sep = "")

  # Lasso regression - wothout repetitions
  if(reruns==1){

    # Fit a GLM with lasso and 10-fold cross-validation
    fit <- glmnet::cv.glmnet(x, y, family, standardize)
    # Extract coefficients
    vcoef <- coef(fit, s = "lambda.min")

    # Lasso regression - with repetitions
  } else {

    mmse <- matrix(ncol=reruns, nrow=100)
    mlambda <- matrix(ncol=reruns, nrow=100)

    for(w in 1:reruns){

      # Fit a GLM with lasso and 10-fold cross-validation
      fit <- glmnet::cv.glmnet(x, y, family=family, standardize = standardize)

      # Extract the values of lambda used in the fit
      if(w == 1) {
        matchLambdas <- 1:length(fit$lambda)
      } else {
        matchLambdas <- wrapr::match_order(mlambda[,1], fit$lambda)
      }
      mlambda[matchLambdas,w] <- fit$lambda[matchLambdas]
      # and the corresponding mean cross-validated error
      mmse[matchLambdas,w] <- fit$cvm[matchLambdas]

      # Output the progress
      #cat(".")
    }

    # Select lambda corresponding to the smallest average MSE
    if(sum(apply(mlambda, 1, function(x) sd(x)), na.rm = TRUE) == 0) {
      lambda.min <- mlambda[which.min(rowMeans(mmse, na.rm=TRUE)), 1]
    } else {
      warning("Lambdas do not match!")
    }

    # Extract coefficients
    vcoef <- coef(fit, s = lambda.min)

  }

  # Step-1-selected main effects (Non-zero coefficients excluding intercept)
  maineffects1 <- which(as.numeric(vcoef)[-1] != 0)

  list(maineffects1 = maineffects1)
}
#'
#' Internal ilasso function
#'
#' @inheritParams ilasso
#'
.ilasso2 <- function(x, y, hierarchy, maineffects1, family, standardize, reruns, ...) {

  if (reruns < 1) {
    warning("reruns<1; set to 1")
    reruns = 1
  }
  reruns <- as.double(reruns)
  #this.call <- match.call()
  y <- drop(y)
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  if (nrowy != nobs)
    stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (",
               nobs, ")", sep = ""))

  # Features:
  # Include all main effects
  # There should be no penalty on main effects selected in step 1
  # Include only interactions between main effects selected in step 1 and another main effect
  if(hierarchy == "weak"){

    # Interactions
    ind <- expand.grid(1:ncol(x), maineffects1)
    out <- apply(ind, 1, function(i) x[,i[1]] * x[,i[2]])
    if (requireNamespace("Matrix", quietly = TRUE)) {
      out <- Matrix::Matrix(out, sparse = TRUE)
    }

  newX <- cbind(x, out)

    # Include only interactions between step-1-selected main effects
    # Include all main effects
    # To honour hierarchy: no penalty on step-1-selected main effects.
  } else if(hierarchy == "strong") {

    # Interactions to include
    ind <- t(combn(maineffects1,2))
    out <- apply(ind, 1, function(i) x[,i[1]] * x[,i[2]])

    newX <- cbind(x, out)
  }

  # Penalty factor that multiplies lambda to allow
  # no shrinkage of selected main effects
  pf <- rep(1, ncol(newX))
  pf[maineffects1] <- 0

  # Lasso regression
  if(reruns==1){

    fit2c <- try(glmnet::cv.glmnet(newX, y, family=family, standardize = standardize,
                                   penalty.factor = pf), silent=TRUE)
    if(inherits(fit2c, "try-error")) stop("cv.glmnet error")

    # Zero and non-zero coefficients
    c2c <- coef(fit2c, s = "lambda.min")

  } else {

    MSE2s <- matrix(ncol=reruns, nrow=100)
    LAMBDA2s <- matrix(ncol=reruns, nrow=100)

    for(w in 1:reruns){
      fit2c <- glmnet::cv.glmnet(newX, y, family=family, standardize = standardize,
                                 penalty.factor = pf)

      # Extract the values of lambda used in the fits
      if(w == 1) {
        matchLambdas <- 1:length(fit2c$lambda)
      } else {
        matchLambdas <- match_order(LAMBDA2s[,1], fit2c$lambda)
      }
      LAMBDA2s[matchLambdas,w] <- fit2c$lambda[matchLambdas]
      # and the corresponding mean cross-validated error
      MSE2s[matchLambdas,w] <- fit2c$cvm[matchLambdas]

      if(w < reruns) rm(fit2c)
      # Output the progress
      #cat(".")
    }

    # Select lambda corresponding to the smallest average MSE
    if(sum(apply(LAMBDA2s, 1, function(i) sd(i)), na.rm = TRUE) == 0) {
      lambda.min <- LAMBDA2s[which.min(rowMeans(MSE2s, na.rm=TRUE)), 1]
    } else {
      warning("Lambdas do not match!")
    }

    # Zero and non-zero coefficients
    c2c <- coef(fit2c, s = lambda.min)
  }

  # Non-zero coefficients (excluding intercept)
  nzc <- which(as.numeric(c2c)[-1] != 0)

  # Step-2-selected main effects
  maineffects2 <- nzc[nzc %in% 1:ncol(x)]

  # Step-2-selected interactions
  interactions2 <- ind[nzc[!(nzc %in% 1:ncol(x))]-ncol(x),]

  # Return
  list(interactions2 = interactions2, maineffects2 = maineffects2)

}
