% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/milasso.R
\name{milasso}
\alias{milasso}
\title{Identify interactions with lasso for multivariate outcome}
\usage{
milasso(
  X,
  Y,
  step = c("both", "first", "second"),
  hierarchy = c("strong", "weak"),
  preserve = TRUE,
  X1,
  Y1,
  standardize = TRUE,
  standardize.response = TRUE,
  fold = 5,
  seed = NULL,
  reruns = 1,
  ...
)
}
\arguments{
\item{X}{input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.}

\item{Y}{response matrix, of dimension nobs x nouts. Quantitative for family="gaussian".}

\item{step}{character. Which step should be run. For \code{step = "both"} or \code{step = "second"}, a value for \code{hierarchy} must also be provided. For \code{step = "second"}, values for \code{X1} must also be provided.}

\item{hierarchy}{character. The hierarchy of the model.}

\item{X1}{numeric vector. This is for the \code{step = "second"}, and allows the user to provide the main effects for which interactions are included (or not, depending on the \code{hierarchy}).}

\item{Y1}{numeric vector. This is for the \code{step = "second"}, and allows the user to provide the responses for which interactions are included (or not, depending on the \code{hierarchy}).}

\item{standardize}{logical. Should X be standardised prior to fitting the model. Defaults to \code{TRUE}.}

\item{standardize.response}{logical. Should Y be standardised prior to fitting the model. Defaults to \code{TRUE}.}

\item{fold}{a positive integer for the corss validation fold. Default is \code{fold = 5}.}

\item{seed}{a numeric scaler, specifying the seed of the random number generator in R for generating cross validation subset for each fold. Default=1.}

\item{reruns}{numeric. The number of reruns of the procedure, the default is 10.}

\item{...}{Other arguments passed to the MSGLasso.cv function.}
}
\value{
An object with S3 class "milasso". Selected main effects and interactions. \item{X1 }{For \code{step = "both"} or \code{step = "first"} only. A vector consisting of the indecies of the main effects selected in step 1.} \item{X1 }{For \code{step = "both"} or \code{step = "second"} only. A vector consisting of the indecies of the main effects selected in step 2.} \item{interactions2 }{For \code{step = "both"} or \code{step = "second"} only. A matrix consisting of the indecies of the main effects (columns) corresponding to the pairwise interactions (rows) selected in step 2.}
}
\description{
Identifies main effects and pairwise interactions in a multivariate linear or logistic regression model via hierarchical lasso penalisation. Can deal with all shapes of data, including very large sparse data matrices.
}
\details{
Step 1 assumes a multivariate pure main effects model and returns the main effects selected by lasso corresponding to the set of non-zero estimates.

Step 2 assumes a multivariate pairwise interaction model subject to a hierarchical restriction and returns the main effects and interactions selected by the MSGLasso corresponding to the set of non-zero estimates.
}
\examples{
# Gaussian
p <- 20
q <- 10
n <- 150
# Working variables
nz <- dplyr::arrange(expand.grid(a=1:p,b=1:p),a)
# Features
set.seed(1)
X <- matrix(sample(0:2, n*p, replace=TRUE), nrow=n, ncol=p)
colnames(X) <- paste0("X[, ", 1:p, "]")
# Interactions
ind <- t(combn(ncol(X),2))
out <- apply(ind, 1, function(x) X[,x[1]] * X[,x[2]])
colnames(out) <- apply(ind, 1, function(x) paste0("X[, ",x[1],"]:X[, ",x[2],"]"))
# Design
design <- cbind(X, out)
# X1 and X2 affect Y1
B <- matrix(0, nrow = p, ncol = q)
B[1, 1] <- 1 # X1 påvirker Y1
B[2, 1] <- 1 # X2 påvirker Y1
BB <- matrix(0, nrow = choose(p, 2), ncol = q)
BB[1,1] <- 1 # VV mellem X1 og X2 på Y1
# Parameters
sigmaE <- c(1,rep(5, q-1))
SigmaE <- diag(x = sigmaE, nrow = q, ncol = q)
Y <- design \%*\%  rbind(B, BB) + MASS::mvrnorm(n = n, rep(0, q), SigmaE)
# Step 1
fit1 <- milasso(X = X, Y = Y, step = "first")
me1x <- fit1$X1
me1y <- fit1$Y1
print(fit1)
# Step 2
fit2 <- milasso(X = X, Y = Y, step = "second", hierarchy = "strong", X1 = me1x, Y1 = me1y)
print(fit2)

}
\seealso{
\code{\link{MSGLasso}} and \code{\link{MSGLasso.cv}}
}
