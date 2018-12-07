library(cubature)

## Test harness
do_test <- function(info_string,
                    expected, ## A list of lists of two named items, code and value, one per method below
                    scalar_f, vector_f,
                    relTol = 1e-3, absTol = 0, nDim = 3, nVec = 256L, nComp = 1L,
                    lowerLimit = rep(0, nDim), upperLimit = rep(1, nDim),
                    methods = list(Cuhre = cubature::cuhre, Divonne = cubature::divonne,
                                   Suave = cubature::suave, Vegas = cubature::vegas)) {
    test_that(info_string, {
        method_names <- names(methods)
        names(expected) <- method_names
        for (name in method_names) {
            method <- methods[[name]]
            expect <- expected[[name]]
            result <- do.call(what = method,
                              args = list(f = scalar_f,
                                          lowerLimit = lowerLimit, upperLimit = upperLimit,
                                          relTol = relTol, absTol = absTol, nComp = nComp))
            ##print(name)
            ##print(result)
            testthat::expect_equal(expect$code, result$returnCode, info = paste(name, " :Integration unsuccessful!"))
            testthat::expect_equal(expect$value, result$integral, tolerance = relTol,
                                   scale = abs(result$integral),
                                   info = paste(name, ": Relative error not reached"))
            ## testthat::expect_equal(expect$value, result$integral, tolerance = relTol, scale = 1,
            ##                        info = paste(name, ": Absolute error not reached"))
            result <- do.call(what = method,
                              args = list(f = vector_f,
                                          lowerLimit = lowerLimit, upperLimit = upperLimit,
                                          relTol = relTol, absTol = absTol, nVec = nVec, nComp = nComp))
            ##print(result)
            testthat::expect_equal(expect$code, result$returnCode, info = paste(name, "Vectorized :Integration unsuccessful!"))
            testthat::expect_equal(expect$value, result$integral, tolerance = relTol,
                                   scale = abs(result$integral),
                                   info = paste(name, "Vectorized : Relative error not reached"))
            ## testthat::expect_equal(expect$value, result$integral, tolerance = relTol, scale = 1,
            ##                        info = paste(name, "Vectorized : Absolute error not reached"))
        }
    })
}

context("Test known Cuba results")

## Numerical integration using Wang-Landau sampling
## Y. W. Li, T. Wust, D. P. Landau, H. Q. Lin
## Computer Physics Communications, 2007, 524-529
## Compare with exact answer: 1.63564436296
##
do_test(
    info_string = "Wang-Landau sampling 1d example",
    expected = rep(list(list(code = 0, value = 1.63564436296)), 4),
    scalar_f = function(x)  sin(4 * x) * x * ((x * ( x * (x * x - 4) + 1) - 1)),
    vector_f = function(x) {
        matrix(apply(x, 2, function(z)
            sin(4 * z) *
            z * ((z * ( z * (z * z - 4) + 1) - 1))),
            ncol = ncol(x))
    },
    relTol = 1e-5,
    lowerLimit = -2, upperLimit = 2, nDim = 1,
    methods = list(Cuhre = cubature::cuhre)
)

do_test(
    info_string = "Factorial",
    expected = rep(list(list(code = 0, value = 0.1666667)), 4),
    scalar_f = function(x) x[1]^(3 - 1) / factorial(3 - 1),
    vector_f = function(x) {
        r <- apply(x, 2, function(z) z[1]^(3-1) / factorial(3 - 1))
        matrix(r, ncol = ncol(x))
    }
)

## do_test(
##     info_string = "Displacement of Origin",
##     ## Suave and Vegas don't return code of 0!
##     expected = list(cuhre = list(code = 0, value = -0.02460334),
##                     divonne = list(code = 0, value = -0.02460334),
##                     suave = list(code = 0, value = -0.02460334),
##                     vegas = list(code = 1, value = -0.02460334)),
##     scalar_f = function(x) sin(x[1] - 3) * cos(x[2] - 2) * exp(x[3] - 1),
##     vector_f = function(x) {
##         apply(x, 2, function(z) sin(z[1] - 3) * cos(z[2] - 2) * exp(z[3] - 1))
##     }
## )

do_test(
    info_string = "Phase and Shift",
    ## Vegas fails
    expected = rep(list(list(code = 0, value = c(0.6646297, 0.3078155))), 3),
    scalar_f = function(arg) {
        x <- arg[1]; y <- arg[2]; z <- arg[3];
        ff <- sin(x) * cos(y) * exp(z);
        gg <-  1 / (3.75 - cos(pi * x) - cos(pi * y) - cos(pi * z));
        return(c(ff, gg))
    },
    vector_f = function(mat) {
        apply(mat, 2,
              function(arg) {
                  x <- arg[1]; y <- arg[2]; z <- arg[3];
                  ff <- sin(x) * cos(y) * exp(z);
                  gg <-  1 / (3.75 - cos(pi * x) - cos(pi * y) - cos(pi * z));
                  return(c(ff, gg))
              })
    },
    nComp = 2L,
    methods = list(Cuhre = cubature::cuhre, Divonne = cubature::divonne,
                   Suave = cubature::suave)
)

test_that("Divonne cuba_phase arg", {
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

## do_test(
##     ## This only passes for certain settings, even in the original version
##     info_string = "Cuba test function 2 (FUN = 2)",
##     expected = rep(list(list(code = 0, value = 5.26851507)), 4),
##     scalar_f = function(x) 1 / ( (x[1] + x[2])^2 + 0.003) * cos(x[2]) * exp(x[3]),
##     vector_f = function(xmat) {
##         matrix(apply(xmat, 2,
##                      function(x) 1 / ( (x[1] + x[2])^2 + 0.003) * cos(x[2]) * exp(x[3])))
##     }

## )

