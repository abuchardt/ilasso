
<!-- README.md is generated from README.Rmd. Please edit that file -->
ilasso
======

<!-- badges: start -->
<!-- badges: end -->
The goal of ilasso is to help identifying interactions via hierarchical lasso regularisation.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("abuchardt/ilasso")
```

Example
-------

This is a basic example on simulated data:

``` r
library(ilasso)

# Gaussian
p <- 200
n <- 100
nz <- dplyr::arrange(expand.grid(a = 1:p, b = 1:p), a)
X <- matrix(sample(c(TRUE, FALSE), n * p, replace = TRUE), nrow = n, ncol = p)
ind <- t(combn(ncol(X),2))
out <- apply(ind, 1, function(x) X[,x[1]] * X[,x[2]])
XX <- cbind(X, out)
beta <- c(3,3, rep(0, ncol(X)-2))
beta12 <- c(3, rep(0, choose(p, 2)-1))
y <- XX %*% c(beta,beta12) + rnorm(n)
```

Run step 1

``` r
# Step 1
fit1 <- ilasso(x = X, y = y, step = "first")
me1 <- fit1$maineffects1
print(fit1)
#> $maineffects1
#> [1]   1   2  99 102 122
#> 
#> attr(,"class")
#> [1] "ilasso"
```

Run step 2 under strong hierarchy (default)

``` r
# Step 2
fit2 <- ilasso(x = X, y = y, step = "second", maineffects1 = me1)
print(fit2)
#> $maineffects2
#> [1]   1   2  99 102 122
#> 
#> $interactions2
#>      [,1] [,2]
#> 
#> attr(,"class")
#> [1] "ilasso"
```

Run both steps under strong hierarchy (default)

``` r
# Both steps
fit12 <- ilasso(x = X, y = y, step = "both")
print(fit2)
#> $maineffects2
#> [1]   1   2  99 102 122
#> 
#> $interactions2
#>      [,1] [,2]
#> 
#> attr(,"class")
#> [1] "ilasso"
```
