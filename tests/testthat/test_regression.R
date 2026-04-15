library(cubature)

# These tests lock in known-good values to catch future regressions.
# The expected values are computed with high precision and should not change
# unless there's a fundamental change in the algorithms.

test_that("product of cosines on [0,1]^2 returns expected value", {
    f <- function(x) prod(cos(x))
    result <- hcubature(f, rep(0, 2), rep(1, 2), tol = 1e-8)

    # sin(1)^2 = 0.7080734182735712...
    expect_equal(result$integral, 0.7080734183, tolerance = 1e-7)
    expect_equal(result$returnCode, 0)
})

test_that("simple polynomial on [0,1]^3 returns exact value", {
    f <- function(x) prod(2 * x)
    result <- hcubature(f, rep(0, 3), rep(1, 3), tol = 1e-8)

    # Integral of 8*x*y*z over [0,1]^3 = 8 * (1/2)^3 = 1
    expect_equal(result$integral, 1.0, tolerance = 1e-10)
    expect_equal(result$returnCode, 0)
})

test_that("Gaussian centered at 1/2 on [0,1]^2 integrates to ~1", {
    f <- function(x) {
        a <- 0.1
        s <- sum((x - 0.5)^2)
        ((2 / sqrt(pi)) / (2 * a))^length(x) * exp(-s / (a * a))
    }
    result <- hcubature(f, rep(0, 2), rep(1, 2), tol = 1e-6)

    expect_equal(result$integral, 1.0, tolerance = 1e-4)
    expect_equal(result$returnCode, 0)
})

test_that("pcubature gives same result for product of cosines", {
    f <- function(x) prod(cos(x))
    result <- pcubature(f, rep(0, 2), rep(1, 2), tol = 1e-8)

    expect_equal(result$integral, 0.7080734183, tolerance = 1e-7)
    expect_equal(result$returnCode, 0)
})

test_that("1D Wang-Landau example returns known value", {
    # From: Y. W. Li, T. Wust, D. P. Landau, H. Q. Lin
    # Computer Physics Communications, 2007, 524-529
    f <- function(x) {
        sin(4 * x) * x * ((x * (x * (x * x - 4) + 1) - 1))
    }

    result <- hcubature(f, -2, 2, tol = 1e-7)
    expect_equal(result$integral, 1.63564436296, tolerance = 1e-6)
})

test_that("2D Wang-Landau example returns known value", {
    f <- function(x) {
        x1 <- x[1]; x2 <- x[2]
        sin(4 * x1 + 1) * cos(4 * x2) * x1 *
            (x1 * (x1 * x1)^2 - x2 * (x2 * x2 - x1) + 2)
    }

    result <- hcubature(f, rep(-1, 2), rep(1, 2), tol = 1e-6)
    expect_equal(result$integral, -0.01797992646, tolerance = 1e-5)
})

test_that("cuhre returns expected value for cosine product", {
    f <- function(x) prod(cos(x))
    result <- cuhre(f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2), relTol = 1e-6)

    expect_equal(result$integral, 0.7080734183, tolerance = 1e-5)
    expect_equal(result$returnCode, 0)
})

test_that("default_args returns expected structure", {
    args <- default_args()

    # Check all methods are present
    expect_named(args, c("hcubature", "pcubature", "cuhre",
                         "divonne", "suave", "vegas"))

    # Check cubature methods have norm
    expect_true("norm" %in% names(args$hcubature))
    expect_true("norm" %in% names(args$pcubature))

    # Check Cuba methods have expected parameters
    expect_true("key" %in% names(args$cuhre))
    expect_true("key1" %in% names(args$divonne))
    expect_true("nNew" %in% names(args$suave))
    expect_true("nStart" %in% names(args$vegas))

    # Check common Cuba parameters
    expect_true("minEval" %in% names(args$cuhre))
    expect_true("stateFile" %in% names(args$divonne))
    expect_true("flags" %in% names(args$suave))
})

test_that("norm options are all present", {
    args <- default_args()
    expected_norms <- c("INDIVIDUAL", "PAIRED", "L2", "L1", "LINF")

    expect_equal(args$hcubature$norm, expected_norms)
    expect_equal(args$pcubature$norm, expected_norms)
})

test_that("vector-valued integral returns expected components", {
    # Integrate (x^2, y^2) over [0,1]^2
    # Each component should integrate to 1/3
    f <- function(x) c(x[1]^2, x[2]^2)

    result <- hcubature(f, c(0, 0), c(1, 1), fDim = 2, tol = 1e-6)

    expect_equal(result$integral[1], 1/3, tolerance = 1e-5)
    expect_equal(result$integral[2], 1/3, tolerance = 1e-5)
    expect_length(result$integral, 2)
    expect_length(result$error, 2)
})

test_that("Tsuda's example returns expected value", {
    f <- function(x) {
        a <- (1 + sqrt(10.0)) / 9.0
        prod(a / (a + 1) * ((a + 1) / (a + x))^2)
    }

    result <- hcubature(f, rep(0, 3), rep(1, 3), tol = 1e-5)
    expect_equal(result$integral, 1.0, tolerance = 1e-4)
})

test_that("Morokoff & Caflisch example returns expected value", {
    f <- function(x) {
        n <- length(x)
        p <- 1/n
        (1 + p)^n * prod(x^p)
    }

    result <- hcubature(f, rep(0, 3), rep(1, 3), tol = 1e-5)
    expect_equal(result$integral, 1.0, tolerance = 1e-4)
})
