library(cubature)
context("Test known integration results")

test_that("Test a product of Cosine functions", {
    expected <- 0.708073
    tol <- 1e-4
    ## Product of cosines
    testFn0 <- function(x) prod(cos(x))
    result <- cubature::adaptIntegrate(f = testFn0,
                                       lowerLimit = rep(0, 2),
                                       upperLimit = rep(1, 2),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                           info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")
})

test_that("Test a Gaussian integral remapped to [0, infinity] limits", {
    expected <- 1.00001
    tol <- 1e-4
    ## Gaussian function
    testFn1 <- function(x) {
        val <- sum(((1 - x) / x)^2)
        scale <- prod((2 / sqrt(pi)) / x^2)
        exp(-val) * scale
    }
    result <- cubature::adaptIntegrate(f = testFn1,
                                       lowerLimit = rep(0, 3),
                                       upperLimit = rep(1, 3),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test volume of a hypersphere (integrating a discountinuous function)", {
    expected <- 0.19728
    tol <- 1e-4
    ## discontinuous objective: volume of hypersphere
    testFn2 <- function(x) {
        radius <- 0.50124145262344534123412
        ifelse(sum(x * x) < radius * radius, 1, 0)
    }
    result <- cubature::adaptIntegrate(f = testFn2,
                                       lowerLimit = rep(0, 2),
                                       upperLimit = rep(1, 2),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test a simple polynomial (product of coordinates)", {
    expected <- 1
    tol <- 1e-4
    ## product of coordinates
    testFn3 <- function(x) prod(2 * x)
    result <- cubature::adaptIntegrate(f = testFn3,
                                       lowerLimit = rep(0, 3),
                                       upperLimit = rep(1, 3),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test Gaussian centered at 1/2", {
    expected <- 1
    tol <- 1e-4
    ## guassian centered at 1/2
    testFn4 <- function(x) {
        a <- 0.1
        s <- sum((x - 0.5)^2)
        ((2 / sqrt(pi)) / (2. * a))^length(x) * exp (-s / (a * a))
    }
    result <- cubature::adaptIntegrate(f = testFn4,
                                       lowerLimit = rep(0, 2),
                                       upperLimit = rep(1, 2),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test double Gaussian", {
    expected <- 0.999994
    tol <- 1e-4
    ## double guassian
    testFn5 <- function(x) {
        a = 0.1
        s1 = sum((x - 1 / 3)^2)
        s2 = sum((x - 2 / 3)^2)
        0.5 * ((2 / sqrt(pi)) / (2. * a))^length(x) * (exp(-s1 / (a * a)) + exp(-s2 / (a * a)))
    }
    result <- cubature::adaptIntegrate(f = testFn5,
                                       lowerLimit = rep(0, 3),
                                       upperLimit = rep(1, 3),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test Tsuda's example", {
    expected <- 0.999998
    tol <- 1e-4
    ## Tsuda's example
    testFn6 <- function(x) {
        a <- (1 + sqrt(10.0)) / 9.0
        prod( a / (a + 1) * ((a + 1) / (a + x))^2)
    }

    result <- cubature::adaptIntegrate(f = testFn6,
                                       lowerLimit = rep(0, 3),
                                       upperLimit = rep(1, 3),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test Morokoff & Calflish example", {
    expected <- 1.00001
    tol <- 1e-4
    ## Morokoff & Calflish function
    testFn7 <- function(x) {
        n <- length(x)
        p <- 1/n
        (1 + p)^n * prod(x^p)
    }
    result <- cubature::adaptIntegrate(f = testFn7,
                                       lowerLimit = rep(0, 3),
                                       upperLimit = rep(1, 3),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test cubature web example", {
    expected <- 13.69609
    tol <- 1e-4
    result <- cubature::adaptIntegrate(f = function(x) exp(-0.5 * sum(x^2)),
                                       lowerLimit = rep(-2, 3),
                                       upperLimit = rep(2, 3),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")

})

test_that("Test Wang-Landau sampling 1d example", {
    expected <- 1.63564436296
    tol <- 1e-5
    ## Numerical integration using Wang-Landau sampling
    ## Y. W. Li, T. Wust, D. P. Landau, H. Q. Lin
    ## Computer Physics Communications, 2007, 524-529
    ## Compare with exact answer: 1.63564436296
    ##
    I.1d <- function(x) {
        sin(4 * x) *
            x * ((x * ( x * (x * x - 4) + 1) - 1))
    }

    result <- cubature::adaptIntegrate(f = I.1d,
                                       lowerLimit = -2,
                                       upperLimit = 2,
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")
})

test_that("Test Wang-Landau sampling 2d example", {
    expected <- -0.01797992646
    tol <- 1e-5
    ## Numerical integration using Wang-Landau sampling
    ## Y. W. Li, T. Wust, D. P. Landau, H. Q. Lin
    ## Computer Physics Communications, 2007, 524-529
    ## Compare with exact answer: 1.63564436296
    ##
    I.2d <- function(x) {
        x1 <- x[1]; x2 <- x[2]
        sin(4 * x1 + 1) * cos(4 * x2) * x1 * (x1 * (x1 * x1)^2 - x2 * (x2 * x2 - x1) +2)
    }


    result <- cubature::adaptIntegrate(f = I.2d,
                                       lowerLimit = rep(-1, 2),
                                       upperLimit = rep(1, 2),
                                       tol = tol)
    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")
})

test_that("Test Multivariate Normal", {
    expected <- 0.3341125
    tol <- 1e-5

    dmvnorm <- function (x, mean, sigma, log = FALSE) {
        if (is.vector(x)) {
            x <- matrix(x, ncol = length(x))
        }
        if (missing(mean)) {
            mean <- rep(0, length = ncol(x))
        }
        if (missing(sigma)) {
            sigma <- diag(ncol(x))
        }
        if (NCOL(x) != NCOL(sigma)) {
            stop("x and sigma have non-conforming size")
        }
        if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                         check.attributes = FALSE)) {
            stop("sigma must be a symmetric matrix")
        }
        if (length(mean) != NROW(sigma)) {
            stop("mean and sigma have non-conforming size")
        }
        distval <- mahalanobis(x, center = mean, cov = sigma)
        logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
        logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
        if (log)
            return(logretval)
        exp(logretval)
    }

    m <- 3
    sigma <- diag(3)
    sigma[2,1] <- sigma[1, 2] <- 3/5 ; sigma[3,1] <- sigma[1, 3] <- 1/3
    sigma[3,2] <- sigma[2, 3] <- 11/15

    result <- cubature::adaptIntegrate(f = dmvnorm,
                                       lowerLimit = rep(-0.5, 3),
                                       upperLimit = c(1, 4, 2),
                                       tol = tol,
                                       mean=rep(0, m), sigma=sigma, log=FALSE)

    testthat::expect_equal(0, result$returnCode, info = "Integration unsuccessful!")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = result$integral,
                 info = "Relative error not reached")

    testthat::expect_equal(expected, result$integral, tolerance = tol, scale = 1,
                           info = "Absolute error not reached")
})
