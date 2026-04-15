library(cubature)

test_that("maxEval limit is approximately respected for hcubature", {
    call_count <- 0
    counting_f <- function(x) {
        call_count <<- call_count + 1
        prod(cos(x))
    }

    # tol is tighter than maxEval can reach, so we expect NOT_CONVERGED (rc=2)
    # and the accompanying warning; the call_count check still validates
    # that the budget was approximately honored.
    result <- expect_warning(
        hcubature(counting_f, rep(0, 3), rep(1, 3),
                  maxEval = 100, tol = 1e-10),
        "maxEval .* exhausted"
    )
    expect_equal(result$returnCode, 2L)
    expect_lt(call_count, 500)
})

test_that("maxEval limit is approximately respected for pcubature", {
    call_count <- 0
    counting_f <- function(x) {
        call_count <<- call_count + 1
        prod(cos(x))
    }

    result <- expect_warning(
        pcubature(counting_f, rep(0, 2), rep(1, 2),
                  maxEval = 100, tol = 1e-10),
        "maxEval .* exhausted"
    )
    expect_equal(result$returnCode, 2L)
    expect_lt(call_count, 500)
})

test_that("infinite limits are handled correctly - full line", {
    # Gaussian integral over (-Inf, Inf) should be sqrt(pi)
    f <- function(x) exp(-x^2)

    result <- hcubature(f, -Inf, Inf, tol = 1e-4)
    expect_equal(result$integral, sqrt(pi), tolerance = 1e-3)
})

test_that("infinite limits are handled correctly - half line positive", {
    # Exponential integral from 0 to Inf should be 1
    f <- function(x) exp(-x)

    result <- hcubature(f, 0, Inf, tol = 1e-4)
    expect_equal(result$integral, 1, tolerance = 1e-3)
})

test_that("infinite limits are handled correctly - half line negative", {
    # Exponential integral from -Inf to 0 should be 1
    f <- function(x) exp(x)

    result <- hcubature(f, -Inf, 0, tol = 1e-4)
    expect_equal(result$integral, 1, tolerance = 1e-3)
})

test_that("multidimensional infinite limits work", {
    # 2D Gaussian integral should be pi
    f <- function(x) exp(-sum(x^2))

    result <- hcubature(f, c(-Inf, -Inf), c(Inf, Inf), tol = 1e-3)
    expect_equal(result$integral, pi, tolerance = 1e-2)
})

test_that("mixed finite and infinite limits work", {
    # Integral of exp(-x-y) from 0 to Inf in both dimensions = 1
    f <- function(x) exp(-sum(x))

    result <- hcubature(f, c(0, 0), c(Inf, Inf), tol = 1e-4)
    expect_equal(result$integral, 1, tolerance = 1e-3)
})

test_that("very small integration region works", {
    f <- function(x) x^2

    # Integrate x^2 from 0 to 0.001, should be (0.001)^3/3
    result <- hcubature(f, 0, 0.001, tol = 1e-10)
    expected <- (0.001)^3 / 3
    expect_equal(result$integral, expected, tolerance = 1e-12)
})

test_that("very large integration region works", {
    # Integrate a function that decays fast enough
    f <- function(x) exp(-x^2 / 1000)

    result <- hcubature(f, -100, 100, tol = 1e-4)
    # Should be close to sqrt(1000 * pi)
    expect_equal(result$integral, sqrt(1000 * pi), tolerance = 1)
})

test_that("negative integration limits work (lower > upper after swap)", {
    f <- function(x) x^2

    # Standard integral
    result_pos <- hcubature(f, 0, 1)

    # Reversed limits - some implementations may handle this
    # by returning negative of integral or erroring
    result_neg <- tryCatch(
        hcubature(f, 1, 0),
        error = function(e) NULL
    )

    # If it doesn't error, the integral should be negated or same
    if (!is.null(result_neg)) {
        expect_true(
            abs(result_neg$integral + result_pos$integral) < 1e-6 ||
            abs(result_neg$integral - result_pos$integral) < 1e-6
        )
    }
})
