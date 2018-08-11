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
    gTilde2_v <- function(x) {
        r <- apply(x, 2, function(z) z[1]^(j-1) / factorial(j - 1))
        matrix(r, ncol = ncol(x))
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

    result <- cubature::vegas(f = gTilde2,
                              lowerLimit = rep(0, nDim),
                              upperLimit = rep(1, nDim),
                              relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- cubature::suave(f = gTilde2,
                              lowerLimit = rep(0, nDim),
                              upperLimit = rep(1, nDim),
                              relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- cubature::divonne(f = gTilde2,
                                lowerLimit = rep(0, nDim),
                                upperLimit = rep(1, nDim),
                                relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

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

    result <- cubature::vegas(f = gTilde2_v,
                              lowerLimit = rep(0, nDim),
                              upperLimit = rep(1, nDim),
                              nVec = 128L,
                              relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- cubature::suave(f = gTilde2_v,
                              lowerLimit = rep(0, nDim),
                              upperLimit = rep(1, nDim),
                              nVec = 128L,
                              relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- cubature::divonne(f = gTilde2_v,
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

test_that("Test based on an example of Pierre de Villemereuil", {
    expected <- 15.00001
    relTol <- 1e-4
    absTol <- 0
    requireNamespace("mvtnorm")
    ## Minimal working example
    mu <- c(0, 10)			        #Mean
    G <- matrix(c(0.5, 0, 0, 1), nrow = 2)	#Some variance-covariance matrix
    P <- matrix(c(1, 0, 0, 2), nrow = 2) 	#Some other VCV matrix

    ## Arbitrary function yielding a scalar
    arb.func<-function(x) x[1] + 0.5 * x[2]
    ## We want to compute the covariance between a vector v and an arbitrary
    ## function of another vector which depends on v
    ## A way to do that is first to compute the expectation of the function
    ## given a value of the vector v
    v <- c(10, 0)
    integrand <- function(x) arb.func(x) * mvtnorm::dmvnorm(x, mu + v, P)
    arb.fun_v <- function(x) c(1, 0.5) %*% x
    integrand_v <- function(x) apply(x, 2, function(z) arb.fun_v(z) * mvtnorm::dmvnorm(z, mu + v, P))

    result <- cuhre(f = integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
                    relTol = relTol, absTol = absTol, maxEval = 10000L)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- cuhre(f = integrand_v, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
                    relTol = relTol, absTol = absTol, nVec = 512L, maxEval = 10000L)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    ## result <- vegas(f = integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
    ##                 relTol = relTol, absTol = absTol, maxEval = 10000L)

    ## testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
    ##                        info = "Relative error not reached")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
    ##                        info = "Absolute error not reached")

    ## result <- vegas(f = integrand_v, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
    ##                 relTol = relTol, absTol = absTol, nVec = 512L, maxEval = 10000L)

    ## testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
    ##                        info = "Relative error not reached")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
    ##                        info = "Absolute error not reached")

    ## result <- suave(f = integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
    ##                 relTol = relTol, absTol = absTol, maxEval = 10000L)

    ## testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")
    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
    ##                        info = "Relative error not reached")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
    ##                        info = "Absolute error not reached")

    ## result <- suave(f = integrand_v, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
    ##                 relTol = relTol, absTol = absTol, nVec = 512L, maxEval = 10000L)

    ## testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
    ##                        info = "Relative error not reached")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
    ##                        info = "Absolute error not reached")

    ## result <- divonne(f = integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
    ##                 relTol = relTol, absTol = absTol, maxEval = 10000L)

    ## testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
    ##                        info = "Relative error not reached")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
    ##                        info = "Absolute error not reached")

    ## result <- divonne(f = integrand_v, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf),
    ##                 relTol = relTol, absTol = absTol, nVec = 512L, maxEval = 10000L)

    ## testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
    ##                        info = "Relative error not reached")

    ## testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
    ##                        info = "Absolute error not reached")

})


test_that("Test displacement of origin", {
    expected <- -0.02460334
    relTol <- 1e-3
    absTol <- 0
    ## Factorial Function
    j <- 3
    nDim <- 3
    integrand2 <- function(x) {
        sin(x[1] - 3) * cos(x[2] - 2) * exp(x[3] - 1)
    } # End integrand2

    integrand2_v <- function(x) {
        apply(x, 2, function(z) sin(z[1] - 3) * cos(z[2] - 2) * exp(z[3] - 1))
    } # End integrand2

    result <- cuhre(f = integrand2,
                    lowerLimit = rep(0, nDim),
                    upperLimit = rep(1, nDim),
                    relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- cuhre(f = integrand2_v,
                    lowerLimit = rep(0, nDim),
                    upperLimit = rep(1, nDim),
                    relTol = relTol, absTol = absTol, nVec = 512L)

    result <- vegas(f = integrand2_v,
                    lowerLimit = rep(0, nDim),
                    upperLimit = rep(1, nDim),
                    relTol = relTol, absTol = absTol,
                    maxEval = 10L^4,
                    nVec = 512L)

    testthat::expect_equal(1, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- suave(f = integrand2_v,
                    lowerLimit = rep(0, nDim),
                    upperLimit = rep(1, nDim),
                    relTol = relTol, absTol = absTol,
                    maxEval = 10L^4,
                    nVec = 512L)

    testthat::expect_equal(1, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- divonne(f = integrand2_v,
                      lowerLimit = rep(0, nDim),
                      upperLimit = rep(1, nDim),
                      relTol = relTol, absTol = absTol,
                      nVec = 512L)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")


})

test_that("Test phase and shift", {
    expected <- c(0.6646297, 0.3078155)
    relTol <- 1e-3
    absTol <- 0
    ## Factorial Function

    integrand3 <- function(arg) {
        x <- arg[1]; y <- arg[2]; z <- arg[3];
        ff <- sin(x) * cos(y) * exp(z);
        gg <-  1 / (3.75 - cos(pi * x) - cos(pi * y) - cos(pi * z));
        return(c(ff, gg))
    } # End integrand3

    integrand3_v <- function(mat) {
        apply(mat, 2,
              function(arg) {
                  x <- arg[1]; y <- arg[2]; z <- arg[3];
                  ff <- sin(x) * cos(y) * exp(z);
                  gg <-  1 / (3.75 - cos(pi * x) - cos(pi * y) - cos(pi * z));
                  return(c(ff, gg))
              })
    }

    result <- divonne(f = integrand3, nComp = 2, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
                      relTol = relTol, absTol = absTol)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

    result <- divonne(f = integrand3_v, nComp = 2, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
                      relTol = relTol, absTol = absTol, nVec = 512L)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")

}) # End integrand3

test_that("Test Divonne cuba_phase arg", {
    expected <- 0.6646098
    relTol <- 1e-3
    absTol <- 0

    NDIM <- 3
    NCOMP <- 1
    NMAX <- 4
    integrand <- function(arg, cuba_phase) {
        x <- arg[1]
        y <- arg[2]
        z <- arg[3]
        ##cat("PHASE", cuba_phase)
        ff <- sin(x)*cos(y)*exp(z);
        return(ff)
    } # fin integrand

    peakf <- function(bounds, nMax) {
        ##  print(bounds)
        nDim <- ncol(bounds)
        x <- matrix(0, ncol = nMax, nrow = nDim)
        pas <- 1 / (nMax - 1)
        ## 1ier point
        x[, 1] <- rep(0, nDim)
        ## Les autres points
        for (i in 2:nMaX) {
            x[, i] <- x[, (i - 1)] + pas
        }
        return(x)
    } #end peakf

    result <- divonne(f = integrand, nComp = NCOMP, lowerLimit = rep(0, 3), upperLimit = rep(1, 3),
                      relTol = relTol, absTol = absTol, peakFinder = peakf)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = relTol, scale = 1,
                           info = "Absolute error not reached")
})

