library(cubature)

test_that("hcubature rejects mismatched limit lengths", {
    f <- function(x) sum(x)

    expect_error(hcubature(f, c(0, 0), c(1)),
                 "lowerLimit and upperLimit must have same length")
})

test_that("hcubature rejects zero-length limits", {
    f <- function(x) sum(x)

    expect_error(hcubature(f, numeric(0), numeric(0)),
                 "dimension >= 1")
})

test_that("hcubature rejects fDim <= 0", {
    f <- function(x) x

    expect_error(hcubature(f, 0, 1, fDim = 0),
                 "dimension >= 1")

    expect_error(hcubature(f, 0, 1, fDim = -1),
                 "dimension >= 1")
})

test_that("pcubature rejects mismatched limit lengths", {
    f <- function(x) sum(x)

    expect_error(pcubature(f, c(0, 0), c(1)),
                 "lowerLimit and upperLimit must have same length")
})

test_that("pcubature rejects zero-length limits", {
    f <- function(x) sum(x)

    expect_error(pcubature(f, numeric(0), numeric(0)),
                 "dimension >= 1")
})

test_that("pcubature warns for dimensions > 3", {
    f <- function(x) sum(x)

    expect_warning(pcubature(f, rep(0, 4), rep(1, 4)),
                   "not recommended for dimensions > 3")
})

test_that("hcubature rejects non-positive tolerance", {
    f <- function(x) x

    expect_error(hcubature(f, 0, 1, tol = 0),
                 "tol should be positive")

    expect_error(hcubature(f, 0, 1, tol = -1),
                 "tol should be positive")
})

test_that("pcubature rejects non-positive tolerance", {
    f <- function(x) x

    expect_error(pcubature(f, 0, 1, tol = 0),
                 "tol should be positive")

    expect_error(pcubature(f, 0, 1, tol = -1),
                 "tol should be positive")
})

test_that("cuhre rejects mismatched limit lengths", {
    f <- function(x) sum(x)

    expect_error(cuhre(f, lowerLimit = c(0, 0), upperLimit = c(1)))
})

test_that("divonne rejects mismatched limit lengths", {
    f <- function(x) sum(x)

    expect_error(divonne(f, lowerLimit = c(0, 0), upperLimit = c(1)))
})

test_that("suave rejects mismatched limit lengths", {
    f <- function(x) sum(x)

    expect_error(suave(f, lowerLimit = c(0, 0), upperLimit = c(1)))
})

test_that("vegas rejects mismatched limit lengths", {
    f <- function(x) sum(x)

    expect_error(vegas(f, lowerLimit = c(0, 0), upperLimit = c(1)))
})

test_that("cubintegrate validates method argument", {
    f <- function(x) x

    expect_error(cubintegrate(f, 0, 1, method = "invalid_method"),
                 "should be one of")
})
