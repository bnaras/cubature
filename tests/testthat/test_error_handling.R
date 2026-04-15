library(cubature)

test_that("hcubature propagates integrand errors", {
    bad_f <- function(x) {
        if (x[1] > 0.5) stop("Intentional error")
        sum(x)
    }

    expect_error(hcubature(bad_f, c(0, 0), c(1, 1)))
})

test_that("pcubature propagates integrand errors", {
    bad_f <- function(x) {
        if (x[1] > 0.5) stop("Intentional error")
        sum(x)
    }

    expect_error(pcubature(bad_f, c(0, 0), c(1, 1)))
})

## Note: Tests for non-numeric integrand returns are skipped because
## the Rcpp layer doesn't safely handle type mismatches and can cause
## crashes. This is documented as a known limitation - users must ensure
## their integrand returns numeric values.

test_that("cuhre propagates integrand errors", {
    bad_f <- function(x) {
        if (x[1] > 0.5) stop("Intentional error")
        sum(x)
    }

    expect_error(cuhre(bad_f, lowerLimit = c(0, 0), upperLimit = c(1, 1)))
})

test_that("vector interface errors on scalar function", {
    scalar_f <- function(x) sum(x)

    # A scalar function doesn't return a matrix, so this should error
    expect_error(
        hcubature(scalar_f, c(0, 0), c(1, 1), vectorInterface = TRUE)
    )
})

## Note: Tests for wrong matrix dimensions in vector interface are skipped
## because the C layer doesn't safely validate dimensions and may produce
## undefined behavior. Users must ensure their vectorized integrand returns
## the correct matrix dimensions (fDim x nPoints).
