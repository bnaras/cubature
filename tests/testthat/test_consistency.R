library(cubature)

test_that("hcubature and pcubature agree on smooth functions in 2D", {
    f <- function(x) prod(cos(x))

    h_result <- hcubature(f, rep(0, 2), rep(1, 2), tol = 1e-6)
    p_result <- pcubature(f, rep(0, 2), rep(1, 2), tol = 1e-6)

    expect_equal(h_result$integral, p_result$integral, tolerance = 1e-5)
})

test_that("hcubature and pcubature agree on polynomials", {
    f <- function(x) prod(x)^2

    h_result <- hcubature(f, rep(0, 2), rep(1, 2), tol = 1e-8)
    p_result <- pcubature(f, rep(0, 2), rep(1, 2), tol = 1e-8)

    expect_equal(h_result$integral, p_result$integral, tolerance = 1e-6)
})

test_that("cubintegrate dispatches correctly to hcubature", {
    f <- function(x) prod(cos(x))

    direct <- hcubature(f, rep(0, 2), rep(1, 2), tol = 1e-5)
    unified <- cubintegrate(f, rep(0, 2), rep(1, 2), method = "hcubature",
                            relTol = 1e-5)

    expect_equal(direct$integral, unified$integral, tolerance = 1e-8)
})

test_that("cubintegrate dispatches correctly to pcubature", {
    f <- function(x) prod(cos(x))

    direct <- pcubature(f, rep(0, 2), rep(1, 2), tol = 1e-5)
    unified <- cubintegrate(f, rep(0, 2), rep(1, 2), method = "pcubature",
                            relTol = 1e-5)

    expect_equal(direct$integral, unified$integral, tolerance = 1e-8)
})

test_that("cubintegrate dispatches correctly to cuhre", {
    f <- function(x) prod(cos(x))

    direct <- cuhre(f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2), relTol = 1e-5)
    unified <- cubintegrate(f, rep(0, 2), rep(1, 2), method = "cuhre",
                            relTol = 1e-5)

    expect_equal(direct$integral, unified$integral, tolerance = 1e-8)
})

test_that("scalar and vector interfaces give same results for hcubature", {
    scalar_f <- function(x) prod(cos(x))
    vector_f <- function(x) {
        matrix(apply(x, 2, function(z) prod(cos(z))), ncol = ncol(x))
    }

    scalar_result <- hcubature(scalar_f, rep(0, 2), rep(1, 2), tol = 1e-6)
    vector_result <- hcubature(vector_f, rep(0, 2), rep(1, 2),
                               vectorInterface = TRUE, tol = 1e-6)

    expect_equal(scalar_result$integral, vector_result$integral,
                 tolerance = 1e-5)
})

test_that("scalar and vector interfaces give same results for pcubature", {
    scalar_f <- function(x) prod(cos(x))
    vector_f <- function(x) {
        matrix(apply(x, 2, function(z) prod(cos(z))), ncol = ncol(x))
    }

    scalar_result <- pcubature(scalar_f, rep(0, 2), rep(1, 2), tol = 1e-6)
    vector_result <- pcubature(vector_f, rep(0, 2), rep(1, 2),
                               vectorInterface = TRUE, tol = 1e-6)

    expect_equal(scalar_result$integral, vector_result$integral,
                 tolerance = 1e-5)
})

test_that("scalar and vector interfaces give same results for cuhre", {
    scalar_f <- function(x) prod(cos(x))
    vector_f <- function(x) {
        apply(x, 2, function(z) prod(cos(z)))
    }

    scalar_result <- cuhre(scalar_f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2), relTol = 1e-6)
    vector_result <- cuhre(vector_f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2),
                           nVec = 128, relTol = 1e-6)

    expect_equal(scalar_result$integral, vector_result$integral,
                 tolerance = 1e-5)
})

test_that("adaptIntegrate is alias for hcubature", {
    f <- function(x) prod(cos(x))

    h_result <- hcubature(f, rep(0, 2), rep(1, 2))
    a_result <- adaptIntegrate(f, rep(0, 2), rep(1, 2))

    expect_equal(h_result$integral, a_result$integral)
    expect_equal(h_result$error, a_result$error)
    expect_equal(h_result$functionEvaluations, a_result$functionEvaluations)
})

test_that("Cuba methods agree on simple integrals", {
    f <- function(x) prod(cos(x))
    tol <- 1e-4
    expected <- 0.7080734  # Known value

    cuhre_result <- cuhre(f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2), relTol = tol)
    divonne_result <- divonne(f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2), relTol = tol)
    suave_result <- suave(f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2), relTol = tol)
    vegas_result <- vegas(f, lowerLimit = rep(0, 2), upperLimit = rep(1, 2), relTol = tol)

    expect_equal(cuhre_result$integral, expected, tolerance = tol * 10)
    expect_equal(divonne_result$integral, expected, tolerance = tol * 10)
    expect_equal(suave_result$integral, expected, tolerance = tol * 10)
    expect_equal(vegas_result$integral, expected, tolerance = tol * 10)
})

test_that("vector-valued integrands work consistently", {
    # Returns 2 components
    f <- function(x) c(x[1]^2, x[2]^2)

    h_result <- hcubature(f, c(0, 0), c(1, 1), fDim = 2)
    p_result <- pcubature(f, c(0, 0), c(1, 1), fDim = 2)
    c_result <- cuhre(f, lowerLimit = c(0, 0), upperLimit = c(1, 1), nComp = 2)

    expected <- c(1/3, 1/3)

    expect_equal(h_result$integral, expected, tolerance = 1e-4)
    expect_equal(p_result$integral, expected, tolerance = 1e-4)
    expect_equal(c_result$integral, expected, tolerance = 1e-4)
})
