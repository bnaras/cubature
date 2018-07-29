library(cubature)
context("Test known Cuba results")

test_that("Test Factorial Function", {
    expected <- 0.1666667
    relTol <- 1e-3
    absTol <- 0
    ## Factorial Function
    j <- 3
    nDim <- 3
    gTilde2 <- function(x) {
	x[1]^(j - 1) / factorial(j - 1)
    }
    result <- cubature::cuhre(f = gTilde2,
                              lowerLimit = rep(0, nDim),
                              upperLimit = rep(1, nDim),
                              relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    gTilde2_v <- function(x) {
        r <- apply(x, 2, function(z) z[1]^(j-1) / factorial(j - 1))
        matrix(r, ncol = ncol(x))
    }

    result <- cubature::cuhre(f = gTilde2_v,
                              lowerLimit = rep(0, nDim),
                              upperLimit = rep(1, nDim),
                              nVec = 128L,
                              relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

})

## test_that("Test based on an example of Pierre de Villemereuil", {
##     expected <- 15.00001
##     relTol <- 1e-4
##     absTol <- 0
##     requireNamespace("mvtnorm")
##     ## Minimal working example
##     mu <- c(0, 10)			        #Mean
##     G <- matrix(c(0.5, 0, 0, 1), nrow = 2)	#Some variance-covariance matrix
##     P <- matrix(c(1, 0, 0, 2), nrow = 2) 	#Some other VCV matrix

##     ## Arbitrary function yielding a scalar
##     arb.func<-function(x) x[1] + 0.5 * x[2]
##     ## We want to compute the covariance between a vector v and an arbitrary
##     ## function of another vector which depends on v
##     ## A way to do that is first to compute the expectation of the function
##     ## given a value of the vector v
##     exp_func<- function(v) {
##         integrand <- function(x) arb.func(x) * mvtnorm::dmvnorm(x, mu + v, P)
##         cuhre(f = integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
##               relTol = relTol, absTol = absTol)$integral
##     }
##     integrand <- function(v) {
##         v * exp_func(v) * mvtnorm::dmvnorm(v, mean =  c(0, 0), sigma = G)
##     }
##     a <- cuhre(f = integrand, nComp = 2, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
##                relTol = relTol, absTol = absTol)
##     arb.fun_v <- function(x) c(1, 0.5) %*% x
##     exp_func_v <- function(v) {
##         integrand <- function(x) apply(x, 2, function(z) arb.fun_v(z) * mvtnorm::dmvnorm(z, mu + v, P))
##         cuhre(f = integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
##               relTol = relTol, absTol = absTol, nVec = 2048)$integral
##     }
##     integrand_v <- function(z) {
##         apply(z, 2, function(v) v * exp_func_v(v) * mvtnorm::dmvnorm(v, mean =  c(0, 0), sigma = G))
##     }
##     a <- cuhre(f = integrand_v, nComp = 2, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
##                relTol = relTol, absTol = absTol, nVec = 2048)

##     ## Then we average the expectation above over all values of v
## })
