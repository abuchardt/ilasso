% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/milasso.R
\name{.milasso1}
\alias{.milasso1}
\title{Internal milasso function}
\usage{
.milasso1(X, Y, standardize, fold, seed, reruns, ...)
}
\arguments{
\item{X}{input matrix, of dimension nobs x nvars; each row is an observation vector. Can be in sparse matrix format.}

\item{Y}{response matrix, of dimension nobs x nouts. Quantitative for family="gaussian".}

\item{standardize}{logical. Should X be standardised prior to fitting the model. Defaults to \code{TRUE}.}

\item{fold}{a positive integer for the corss validation fold. Default is \code{fold = 5}.}

\item{seed}{a numeric scaler, specifying the seed of the random number generator in R for generating cross validation subset for each fold. Default=1.}

\item{reruns}{numeric. The number of reruns of the procedure, the default is 10.}

\item{...}{Other arguments passed to the MSGLasso.cv function.}
}
\description{
Internal milasso function
}
