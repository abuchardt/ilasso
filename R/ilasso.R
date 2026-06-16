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
#' @importFrom stats coef sd
#' @importFrom utils combn
#' @importFrom wrapr match_order
#'
#' @import glmnet
#' @import dplyr
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
#' X <- matrix(sample(0:2, n*p, replace=TRUE), nrow=n, ncol=p)
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
                   family = c("gaussian", "binomial"),
                   standardize = TRUE, reruns = 10, ...) {

  this.call <- match.call()

  family <- match.arg(family)
  step <- match.arg(step)
  hierarchy <- match.arg(hierarchy)

  x <- .prepare_x(x)
  y <- .prepare_y(y, family = family)

  if (step == "second") {
    if (missing(maineffects1)) {
      stop("Step 2 requires maineffects1.", call. = FALSE)
    }
    maineffects1 <- as.integer(maineffects1)
    maineffects1 <- maineffects1[!is.na(maineffects1)]
    maineffects1 <- sort(unique(maineffects1))

    if (hierarchy == "weak" && length(maineffects1) < 1) {
      stop("Weak hierarchy needs at least one selected main effect.", call. = FALSE)
    }
    if (hierarchy == "strong" && length(maineffects1) < 2) {
      stop("Strong hierarchy needs at least two selected main effects.", call. = FALSE)
    }
  }

  if (step %in% c("both", "first")) {
    step1 <- .ilasso1(x, y, family, standardize, reruns, ...)
    step2 <- NULL
  }

  if (step %in% c("both", "second")) {
    if (step == "both") {
      if (length(step1$maineffects1) < 1) {
        stop("No main effects were selected in step 1.", call. = FALSE)
      } else if (length(step1$maineffects1) == 1 && hierarchy == "strong") {
        stop("Strong hierarchy needs two or more main effects from step 1; only one was selected.", call. = FALSE)
      }
      maineffects1 <- step1$maineffects1
    } else {
      step1 <- NULL
    }

    step2 <- .ilasso2(x, y, hierarchy, maineffects1, family, standardize, reruns, ...)
  }

  if (step == "first") {
    fit <- list(call = this.call,
                hierarchy = hierarchy,
                coef1 = step1$coef1,
                maineffects1 = step1$maineffects1)
  } else if (step == "second") {
    fit <- list(call = this.call,
                hierarchy = hierarchy,
                coef2 = step2$coef2,
                maineffects2 = step2$maineffects2,
                interactions2 = step2$interactions2,
                interactions_all = step2$interactions_all,
                newX = step2$newX)
  } else {
    fit <- list(call = this.call,
                hierarchy = hierarchy,
                coef1 = step1$coef1,
                coef2 = step2$coef2,
                maineffects1 = step1$maineffects1,
                maineffects2 = step2$maineffects2,
                interactions2 = step2$interactions2,
                interactions_all = step2$interactions_all,
                newX = step2$newX)
  }

  class(fit) <- "ilasso"
  fit
}

.prepare_x <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    return(x)
  }
  as.matrix(x)
}

.prepare_y <- function(y, family) {
  y <- drop(y)

  if (family == "gaussian") {
    y_num <- suppressWarnings(as.numeric(y))
    if (anyNA(y_num)) {
      stop("For family='gaussian', y must be coercible to numeric without NAs.", call. = FALSE)
    }
    return(y_num)
  }

  # glmnet accepts a two-column matrix for binomial responses; keep it as-is.
  if (is.matrix(y) && ncol(y) == 2) {
    return(y)
  }

  if (is.logical(y)) {
    return(as.integer(y))
  }

  if (is.factor(y) || is.character(y)) {
    y_fac <- as.factor(y)
    if (nlevels(y_fac) != 2) {
      stop("For family='binomial', factor/character y must have exactly two levels.", call. = FALSE)
    }
    return(y_fac)
  }

  y_num <- suppressWarnings(as.numeric(y))
  if (anyNA(y_num)) {
    stop("For family='binomial', y must be 0/1, logical, two-level factor/character, or a two-column matrix.", call. = FALSE)
  }
  if (!all(y_num %in% c(0, 1))) {
    stop("For family='binomial', numeric y must contain only 0 and 1.", call. = FALSE)
  }
  y_num
}

.check_xy <- function(x, y) {
  np <- dim(x)
  if (is.null(np) || np[2] <= 1) {
    stop("x should be a matrix with 2 or more columns.", call. = FALSE)
  }
  nobs <- as.integer(np[1])
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  if (nrowy != nobs) {
    stop(paste0("number of observations in y (", nrowy,
                ") not equal to the number of rows of x (", nobs, ")"),
         call. = FALSE)
  }
  invisible(TRUE)
}

.cv_glmnet_repeated <- function(x, y, family, standardize, reruns, penalty.factor = NULL, ...) {
  if (reruns < 1) {
    warning("reruns < 1; set to 1")
    reruns <- 1
  }
  reruns <- as.integer(reruns)

  cv_args <- list(
    x = x,
    y = y,
    family = family,
    standardize = standardize,
    ...
  )
  if (!is.null(penalty.factor)) {
    cv_args$penalty.factor <- penalty.factor
  }

  if (reruns == 1) {
    fit <- do.call(glmnet::cv.glmnet, cv_args)
    return(list(fit = fit, lambda = "lambda.min"))
  }

  fits <- vector("list", reruns)
  for (w in seq_len(reruns)) {
    fits[[w]] <- do.call(glmnet::cv.glmnet, cv_args)
  }

  # Robust repeated-CV summary. Rather than relying on exact lambda matching across
  # runs, use the median lambda.min from the repeated fits.
  lambda.min <- stats::median(vapply(fits, function(z) z$lambda.min, numeric(1)))

  list(fit = fits[[reruns]], lambda = lambda.min)
}

.ilasso1 <- function(x, y, family, standardize, reruns, ...) {
  .check_xy(x, y)

  cv <- .cv_glmnet_repeated(
    x = x,
    y = y,
    family = family,
    standardize = standardize,
    reruns = reruns,
    ...
  )

  vcoef <- coef(cv$fit, s = cv$lambda)
  maineffects1 <- which(as.numeric(vcoef)[-1] != 0)

  list(coef1 = vcoef,
       maineffects1 = maineffects1)
}

.ilasso2 <- function(x, y, hierarchy, maineffects1, family, standardize, reruns, ...) {
  .check_xy(x, y)

  nvars <- ncol(x)
  maineffects1 <- sort(unique(as.integer(maineffects1)))

  if (any(maineffects1 < 1 | maineffects1 > nvars)) {
    stop("maineffects1 contains indices outside the columns of x.", call. = FALSE)
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X", seq_len(nvars))
  }

  if (hierarchy == "weak") {
    # All pairwise interactions with at least one step-1-selected main effect.
    # This avoids squares and duplicate pairs.
    all_pairs <- utils::combn(seq_len(nvars), 2)
    keep <- all_pairs[1, ] %in% maineffects1 | all_pairs[2, ] %in% maineffects1
    ind <- t(all_pairs[, keep, drop = FALSE])
  } else {
    # Only interactions between step-1-selected main effects.
    ind <- t(utils::combn(maineffects1, 2))
  }

  if (nrow(ind) == 0) {
    stop("No pairwise interactions are available under the requested hierarchy.", call. = FALSE)
  }

  out <- .make_interactions(x, ind)
  colnames(out) <- paste0(colnames(x)[ind[, 1]], ":", colnames(x)[ind[, 2]])

  if (inherits(x, "sparseMatrix") || inherits(out, "sparseMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The Matrix package is required for sparse input.", call. = FALSE)
    }
    newX <- cbind(x, out)
  } else {
    newX <- cbind(x, out)
  }

  # Penalty factor: selected step-1 main effects are unpenalised to honour hierarchy.
  pf <- rep(1, ncol(newX))
  pf[maineffects1] <- 0

  cv <- .cv_glmnet_repeated(
    x = newX,
    y = y,
    family = family,
    standardize = standardize,
    reruns = reruns,
    penalty.factor = pf,
    ...
  )

  c2c <- coef(cv$fit, s = cv$lambda)
  nzc <- which(as.numeric(c2c)[-1] != 0)

  maineffects2 <- nzc[nzc <= nvars]

  interaction_pos <- nzc[nzc > nvars] - nvars
  if (length(interaction_pos) == 0) {
    interactions2 <- matrix(integer(0), ncol = 2,
                            dimnames = list(NULL, c("var1", "var2")))
  } else {
    interactions2 <- ind[interaction_pos, , drop = FALSE]
    colnames(interactions2) <- c("var1", "var2")
  }

  # list(coef2 = c2c,
  #      newX = newX,
  #      interactions2 = interactions2,
  #      maineffects2 = maineffects2)
  list(
    coef2 = c2c,
    newX = newX,
    interactions_all = ind,
    interactions2 = interactions2,
    maineffects2 = maineffects2
  )
}

.make_interactions <- function(x, ind) {
  if (inherits(x, "sparseMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The Matrix package is required for sparse input.", call. = FALSE)
    }
    out <- lapply(seq_len(nrow(ind)), function(k) x[, ind[k, 1]] * x[, ind[k, 2]])
    out <- do.call(cbind, out)
    return(out)
  }

  out <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(ind))
  for (k in seq_len(nrow(ind))) {
    out[, k] <- x[, ind[k, 1]] * x[, ind[k, 2]]
  }
  out
}
