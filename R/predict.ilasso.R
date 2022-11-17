#' Extract coefficients from an ilasso object
#'
#' Similar to other predict methods, this function predicts fitted values and extract coefficients from a fitted "edgwas" object.
#'
#' @param object Fitted "ilasso" object.
#' @param newx Matrix of new values for x of dimension nobs x nouts.
#' @param ... additional arguments affecting the predictions produced.
#'
#' @return The object returned is a vector of length nobs of the fitted values.
#'
#' @examples
#' # Gaussian
#' p <- 20
#' n <- 100
#' nz <- dplyr::arrange(expand.grid(a = 1:p, b = 1:p), a)
#' X <- matrix(sample(0:2, n*p, replace=TRUE), nrow=n, ncol=p)
#' ind <- t(combn(ncol(X),2))
#' out <- apply(ind, 1, function(x) X[,x[1]] * X[,x[2]])
#' XX <- cbind(X, out)
#' beta <- c(3,3, rep(0, ncol(X)-2))
#' beta12 <- c(3, rep(0, choose(p, 2)-1))
#' y <- XX %*% c(beta,beta12) + rnorm(n)
#' # Fit weak hierarchy
#' fit <- ilasso(x = X, y = y, hierarchy = "weak")
#' # Predict
#' newx <- matrix(sample(0:2, n*p, replace=TRUE), nrow=n, ncol=p)
#' yHat <- predict(fit, newx = newx)
#'
#' @export
#'
predict.ilasso <- function(object, newx, ...
){
  if (missing(newx)){
    stop("Value for 'newx' missing")
  }

  nouts <- ncol(newx)
  coefs <- object$coef2
  hierarchy <- object$hierarchy

  colnames(newx) <- paste0("X[, ", 1:nouts, "]")

  if(length(object$interactions2) == 0) {

    design <- cbind(1, newx)

  } else {

    if(hierarchy == "strong") {

      ind <- t(combn(object$maineffects1,2))
      out <- apply(ind, 1, function(i) newx[,i[1]] * newx[,i[2]])
      colnames(out) <- apply(ind, 1, function(x)
        paste0("X[, ",x[1],"]:X[, ",x[2],"]")
      )
      design <- cbind(1, newx, out)

    } else if(hierarchy == "weak") {

      ind <- expand.grid(1:nouts, object$maineffects1)
      out <- apply(ind, 1, function(i) newx[,i[1]] * newx[,i[2]])
      colnames(out) <- apply(ind, 1, function(x)
        paste0("X[, ",x[1],"]:X[, ",x[2],"]")
      )
      design <- cbind(1, newx, out)

    }

  }

  nfit <- design %*% coefs

  nfit
}
