# Cubature Vectorization Results

## Introduction

This R `cubature` package exposes both the `hcubature` and `pcubature`
routines of the underlying C `cubature` library, including the
vectorized interfaces.

Per the documentation, use of `pcubature` is advisable only for smooth
integrands in dimensions up to three at most. In fact, the `pcubature`
routines perform significantly worse than the vectorized `hcubature` in
inappropriate cases. So when in doubt, you are better off using
`hcubature`.

Version 2.0 of this package integrates the `Cuba` library as well, once
again providing vectorized interfaces.

The main point of this note is to examine the difference vectorization
makes. My recommendations are below in the summary section.

## A Timing Harness

Our harness will provide timing results for `hcubature`, `pcubature`
(where appropriate) and Cuba `cuhre` calls. We begin by creating a
harness for these calls.

``` r
library(bench)
library(cubature)

harness <- function(which = NULL,
                    f, fv, lowerLimit, upperLimit, tol = 1e-3, iterations = 20, ...) {
  
  fns <- c(hc = "Non-vectorized Hcubature",
           hc.v = "Vectorized Hcubature",
           pc = "Non-vectorized Pcubature",
           pc.v = "Vectorized Pcubature",
           cc = "Non-vectorized cubature::cuhre",
           cc_v = "Vectorized cubature::cuhre")
  cc <- function() cubature::cuhre(f = f,
                                   lowerLimit = lowerLimit, upperLimit = upperLimit,
                                   relTol = tol,
                                   ...)
  cc_v <- function() cubature::cuhre(f = fv,
                                     lowerLimit = lowerLimit, upperLimit = upperLimit,
                                     relTol = tol,
                                     nVec = 1024L,
                                     ...)
  hc <- function() cubature::hcubature(f = f,
                                       lowerLimit = lowerLimit,
                                       upperLimit = upperLimit,
                                       tol = tol,
                                       ...)
  hc.v <- function() cubature::hcubature(f = fv,
                                         lowerLimit = lowerLimit,
                                         upperLimit = upperLimit,
                                         tol = tol,
                                         vectorInterface = TRUE,
                                         ...)
  pc <- function() cubature::pcubature(f = f,
                                       lowerLimit = lowerLimit,
                                       upperLimit = upperLimit,
                                       tol = tol,
                                       ...)
  pc.v <- function() cubature::pcubature(f = fv,
                                         lowerLimit = lowerLimit,
                                         upperLimit = upperLimit,
                                         tol = tol,
                                         vectorInterface = TRUE,
                                         ...)

  ndim <- length(lowerLimit)
  
  if (is.null(which)) {
    fnIndices <- seq_along(fns)
  } else {
    fnIndices <- na.omit(match(which, names(fns)))
  }
  fnList <- lapply(names(fns)[fnIndices], function(x) call(x))
  
  argList <- c(fnList, iterations = iterations, check = FALSE)
  result <- do.call(bench::mark, args = argList)
  d <- result[seq_along(fnIndices), 1:9]
  d$expression <- fns[fnIndices]
  d
}
```

We reel off the timing runs.

## Example 1.

``` r
func <- function(x) sin(x[1]) * cos(x[2]) * exp(x[3])
func.v <- function(x) {
    matrix(apply(x, 2, function(z) sin(z[1]) * cos(z[2]) * exp(z[3])), ncol = ncol(x))
}

d <- harness(f = func, fv = func.v,
             lowerLimit = rep(0, 3),
             upperLimit = rep(1, 3),
             tol = 1e-5,
             iterations = 100)[, 1:9]

knitr::kable(d, digits = 3, row.names = FALSE)
```

| expression                     |      min |   median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 271.24µs | 284.01µs | 3491.207 |    59.3KB | 35.265 |    99 |    1 |     28.4ms |
| Vectorized Hcubature           |  390.3µs | 412.85µs | 2388.755 |    54.5KB | 24.129 |    99 |    1 |     41.4ms |
| Non-vectorized Pcubature       |  867.5µs | 894.27µs | 1102.125 |    31.2KB | 34.086 |    97 |    3 |       88ms |
| Vectorized Pcubature           |   1.18ms |   1.22ms |  811.554 |    58.4KB | 16.562 |    98 |    2 |    120.8ms |
| Non-vectorized cubature::cuhre | 589.32µs | 609.42µs | 1624.346 |    40.5KB | 16.408 |    99 |    1 |     60.9ms |
| Vectorized cubature::cuhre     | 596.45µs |  630.2µs | 1583.485 |    39.8KB | 32.316 |    98 |    2 |     61.9ms |

## Multivariate Normal

Using `cubature`, we evaluate $$\int_{R}\phi(x)dx$$ where $\phi(x)$ is
the three-dimensional multivariate normal density with mean 0, and
variance $$\Sigma = \left( \begin{array}{rrr}
1 & \frac{3}{5} & \frac{1}{3} \\
\frac{3}{5} & 1 & \frac{11}{15} \\
\frac{1}{3} & \frac{11}{15} & 1
\end{array} \right)$$ and $R$ is
$\left\lbrack - \frac{1}{2},1 \right\rbrack \times \left\lbrack - \frac{1}{2},4 \right\rbrack \times \left\lbrack - \frac{1}{2},2 \right\rbrack.$

We construct a scalar function (`my_dmvnorm`) and a vector analog
(`my_dmvnorm_v`). First the functions.

``` r
m <- 3
sigma <- diag(3)
sigma[2,1] <- sigma[1, 2] <- 3/5 ; sigma[3,1] <- sigma[1, 3] <- 1/3
sigma[3,2] <- sigma[2, 3] <- 11/15
logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
my_dmvnorm <- function (x, mean, sigma, logdet) {
    x <- matrix(x, ncol = length(x))
    distval <- stats::mahalanobis(x, center = mean, cov = sigma)
    exp(-(3 * log(2 * pi) + logdet + distval)/2)
}

my_dmvnorm_v <- function (x, mean, sigma, logdet) {
    distval <- stats::mahalanobis(t(x), center = mean, cov = sigma)
    exp(matrix(-(3 * log(2 * pi) + logdet + distval)/2, ncol = ncol(x)))
}
```

Now the timing.

``` r
d <- harness(f = my_dmvnorm, fv = my_dmvnorm_v,
             lowerLimit = rep(-0.5, 3),
             upperLimit = c(1, 4, 2),
             tol = 1e-5,
             iterations = 10,
             mean = rep(0, m), sigma = sigma, logdet = logdet)
knitr::kable(d, digits = 3)
```

| expression                     |      min |   median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 815.22ms | 829.15ms |   1.189 |  139.68KB | 18.429 |    10 |  155 |      8.41s |
| Vectorized Hcubature           |   1.77ms |   2.45ms | 438.862 |    1.89MB |  0.000 |    10 |    0 |    22.79ms |
| Non-vectorized Pcubature       | 353.16ms |  359.1ms |   2.744 |        0B | 17.836 |    10 |   65 |      3.64s |
| Vectorized Pcubature           |   1.24ms |   1.27ms | 661.570 |  810.26KB | 66.157 |    10 |    1 |    15.12ms |
| Non-vectorized cubature::cuhre | 342.23ms | 344.47ms |   2.900 |        0B | 18.273 |    10 |   63 |      3.45s |
| Vectorized cubature::cuhre     |   3.22ms |   3.28ms | 303.423 |  898.41KB |  0.000 |    10 |    0 |    32.96ms |

The effect of vectorization is huge. So it makes sense for users to
vectorize the integrands as much as possible for efficiency.

Furthermore, for this particular example, we expect
[`mvtnorm::pmvnorm`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html) to do
pretty well since it is specialized for the multivariate normal. The
vectorized versions of `hcubature` and `pcubature` seem competitive and
in some cases better, if you compare the table above to the one below.

``` r
library(mvtnorm)
g1 <- function() pmvnorm(lower = rep(-0.5, m),
                                  upper = c(1, 4, 2), mean = rep(0, m), corr = sigma,
                                  alg = Miwa(), abseps = 1e-5, releps = 1e-5)
g2 <- function() pmvnorm(lower = rep(-0.5, m),
                         upper = c(1, 4, 2), mean = rep(0, m), corr = sigma,
                         alg = GenzBretz(), abseps = 1e-5, releps = 1e-5)
g3 <- function() pmvnorm(lower = rep(-0.5, m),
                         upper = c(1, 4, 2), mean = rep(0, m), corr = sigma,
                         alg = TVPACK(), abseps = 1e-5, releps = 1e-5)

knitr::kable(bench::mark(g1(), g2(), g3(), iterations = 20, check = FALSE)[, 1:9],
             digits = 3, row.names = FALSE)
```

| expression |   min | median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-----------|------:|-------:|--------:|----------:|-------:|------:|-----:|-----------:|
| g1()       | 765µs | 1.37ms | 737.572 |  355.13KB |      0 |    20 |    0 |     27.1ms |
| g2()       | 750µs | 1.37ms | 760.937 |    2.49KB |      0 |    20 |    0 |     26.3ms |
| g3()       | 754µs | 1.37ms | 725.544 |    2.49KB |      0 |    20 |    0 |     27.6ms |

## Product of cosines

``` r
testFn0 <- function(x) prod(cos(x))
testFn0_v <- function(x) matrix(apply(x, 2, function(z) prod(cos(z))), ncol = ncol(x))

d <- harness(f = testFn0, fv = testFn0_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 1000)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  36.7µs |  39.2µs | 24593.047 |    7.66KB | 24.618 |   999 |    1 |     40.6ms |
| Vectorized Hcubature           |  64.4µs |  67.8µs | 14338.266 |    1.16KB | 28.734 |   998 |    2 |     69.6ms |
| Non-vectorized Pcubature       |  49.2µs |  50.9µs | 19137.458 |        0B | 19.157 |   999 |    1 |     52.2ms |
| Vectorized Pcubature           | 112.7µs |   120µs |  8120.587 |   18.68KB | 24.435 |   997 |    3 |    122.8ms |
| Non-vectorized cubature::cuhre | 333.3µs | 350.5µs |  2836.890 |        0B | 25.764 |   991 |    9 |    349.3ms |
| Vectorized cubature::cuhre     | 365.9µs | 384.1µs |  2592.187 |   16.38KB | 26.184 |   990 |   10 |    381.9ms |

## Gaussian function

``` r
testFn1 <- function(x) {
    val <- sum(((1 - x) / x)^2)
    scale <- prod((2 / sqrt(pi)) / x^2)
    exp(-val) * scale
}

testFn1_v <- function(x) {
    val <- matrix(apply(x, 2, function(z) sum(((1 - z) / z)^2)), ncol(x))
    scale <- matrix(apply(x, 2, function(z) prod((2 / sqrt(pi)) / z^2)), ncol(x))
    exp(-val) * scale
}

d <- harness(f = testFn1, fv = testFn1_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 10)

knitr::kable(d, digits = 3)
```

| expression                     |      min |   median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |    2.9ms |   2.92ms |   336.492 |    67.3KB |  0.000 |    10 |    0 |    29.72ms |
| Vectorized Hcubature           |   4.95ms |   5.04ms |   198.542 |   290.5KB | 49.636 |     8 |    2 |    40.29ms |
| Non-vectorized Pcubature       |  73.25µs |  77.39µs | 12325.295 |        0B |  0.000 |    10 |    0 |   811.34µs |
| Vectorized Pcubature           | 154.27µs | 162.05µs |  6023.480 |     4.1KB |  0.000 |    10 |    0 |     1.66ms |
| Non-vectorized cubature::cuhre |  13.66ms |  13.83ms |    72.333 |        0B | 31.000 |     7 |    3 |    96.78ms |
| Vectorized cubature::cuhre     |  20.63ms |  20.67ms |    48.225 |   971.5KB | 48.225 |     5 |    5 |   103.68ms |

## Discontinuous function

``` r
testFn2 <- function(x) {
    radius <- 0.50124145262344534123412
    ifelse(sum(x * x) < radius * radius, 1, 0)
}

testFn2_v <- function(x) {
    radius <- 0.50124145262344534123412
    matrix(apply(x, 2, function(z) ifelse(sum(z * z) < radius * radius, 1, 0)), ncol = ncol(x))
}

d <- harness(which = c("hc", "hc.v", "cc", "cc_v"),
             f = testFn2, fv = testFn2_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 10)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  47.8ms |  48.5ms |  17.513 |    17.7KB | 35.026 |    10 |   20 |   571.01ms |
| Vectorized Hcubature           |  48.9ms |  49.6ms |  19.989 |  1011.8KB | 25.985 |    10 |   13 |   500.28ms |
| Non-vectorized cubature::cuhre | 814.1ms | 827.8ms |   1.199 |        0B | 28.168 |    10 |  235 |      8.34s |
| Vectorized cubature::cuhre     | 910.4ms | 924.4ms |   1.074 |    21.2MB | 24.390 |    10 |  227 |      9.31s |

## A Simple Polynomial (product of coordinates)

``` r
testFn3 <- function(x) prod(2 * x)
testFn3_v <- function(x) matrix(apply(x, 2, function(z) prod(2 * z)), ncol = ncol(x))

d <- harness(f = testFn3, fv = testFn3_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 50)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc |  gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|--------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  59.2µs |  60.9µs | 15714.003 |    6.75KB |   0.000 |    50 |    0 |     3.18ms |
| Vectorized Hcubature           |  91.5µs |  96.7µs |  9996.213 |    2.91KB |   0.000 |    50 |    0 |        5ms |
| Non-vectorized Pcubature       |  52.3µs |  54.8µs | 17614.605 |        0B |   0.000 |    50 |    0 |     2.84ms |
| Vectorized Pcubature           |  80.5µs |  83.7µs | 11108.906 |   19.12KB | 226.712 |    49 |    1 |     4.41ms |
| Non-vectorized cubature::cuhre | 650.7µs | 674.2µs |  1478.406 |        0B |  30.172 |    49 |    1 |    33.14ms |
| Vectorized cubature::cuhre     | 658.6µs | 694.3µs |  1439.396 |   39.84KB |   0.000 |    50 |    0 |    34.74ms |

## Gaussian centered at $\frac{1}{2}$

``` r
testFn4 <- function(x) {
    a <- 0.1
    s <- sum((x - 0.5)^2)
    ((2 / sqrt(pi)) / (2. * a))^length(x) * exp (-s / (a * a))
}

testFn4_v <- function(x) {
    a <- 0.1
    r <- apply(x, 2, function(z) {
        s <- sum((z - 0.5)^2)
        ((2 / sqrt(pi)) / (2. * a))^length(z) * exp (-s / (a * a))
    })
    matrix(r, ncol = ncol(x))
}

d <- harness(f = testFn4, fv = testFn4_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |    min | median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|-------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 1.37ms |  1.4ms | 713.642 |    85.2KB | 37.560 |    19 |    1 |     26.6ms |
| Vectorized Hcubature           | 1.79ms | 1.83ms | 541.683 |   147.2KB | 28.510 |    19 |    1 |     35.1ms |
| Non-vectorized Pcubature       | 2.03ms | 2.08ms | 480.350 |        0B | 25.282 |    19 |    1 |     39.6ms |
| Vectorized Pcubature           | 2.61ms | 2.68ms | 372.897 |    68.9KB | 41.433 |    18 |    2 |     48.3ms |
| Non-vectorized cubature::cuhre |  3.3ms | 3.36ms | 297.015 |        0B | 33.002 |    18 |    2 |     60.6ms |
| Vectorized cubature::cuhre     | 3.77ms | 3.87ms | 258.418 |   125.6KB | 13.601 |    19 |    1 |     73.5ms |

## Double Gaussian

``` r
testFn5 <- function(x) {
    a <- 0.1
    s1 <- sum((x - 1 / 3)^2)
    s2 <- sum((x - 2 / 3)^2)
    0.5 * ((2 / sqrt(pi)) / (2. * a))^length(x) * (exp(-s1 / (a * a)) + exp(-s2 / (a * a)))
}
testFn5_v <- function(x) {
    a <- 0.1
    r <- apply(x, 2, function(z) {
        s1 <- sum((z - 1 / 3)^2)
        s2 <- sum((z - 2 / 3)^2)
        0.5 * ((2 / sqrt(pi)) / (2. * a))^length(z) * (exp(-s1 / (a * a)) + exp(-s2 / (a * a)))
    })
    matrix(r, ncol = ncol(x))
}

d <- harness(f = testFn5, fv = testFn5_v,
             lowerLimit = rep(0, 2), upperLimit = rep(1, 2), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |    min | median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|-------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 3.66ms | 3.71ms | 267.753 |     133KB | 47.250 |    17 |    3 |     63.5ms |
| Vectorized Hcubature           | 4.72ms | 4.79ms | 208.361 |     249KB | 23.151 |    18 |    2 |     86.4ms |
| Non-vectorized Pcubature       | 2.54ms | 2.59ms | 384.867 |        0B | 42.763 |    18 |    2 |     46.8ms |
| Vectorized Pcubature           | 3.33ms | 3.42ms | 290.594 |      69KB | 15.294 |    19 |    1 |     65.4ms |
| Non-vectorized cubature::cuhre |  7.2ms | 7.26ms | 137.650 |        0B | 45.883 |    15 |    5 |      109ms |
| Vectorized cubature::cuhre     | 8.39ms | 8.46ms | 117.819 |     224KB | 20.792 |    17 |    3 |    144.3ms |

## Tsuda’s Example

``` r
testFn6 <- function(x) {
    a <- (1 + sqrt(10.0)) / 9.0
    prod( a / (a + 1) * ((a + 1) / (a + x))^2)
}

testFn6_v <- function(x) {
    a <- (1 + sqrt(10.0)) / 9.0
    r <- apply(x, 2, function(z) prod( a / (a + 1) * ((a + 1) / (a + z))^2))
    matrix(r, ncol = ncol(x))
}

d <- harness(f = testFn6, fv = testFn6_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  1.58ms |  1.61ms | 618.345 |    64.6KB | 32.544 |    19 |    1 |     30.7ms |
| Vectorized Hcubature           |  2.13ms |  2.16ms | 461.510 |     156KB | 24.290 |    19 |    1 |     41.2ms |
| Non-vectorized Pcubature       |  8.22ms |  8.35ms | 120.199 |        0B | 40.066 |    15 |    5 |    124.8ms |
| Vectorized Pcubature           | 10.39ms | 10.57ms |  94.089 |   386.2KB | 31.363 |    15 |    5 |    159.4ms |
| Non-vectorized cubature::cuhre |  4.33ms |  4.42ms | 226.170 |        0B | 25.130 |    18 |    2 |     79.6ms |
| Vectorized cubature::cuhre     |  4.83ms |  4.92ms | 202.758 |   225.8KB | 22.529 |    18 |    2 |     88.8ms |

## Morokoff & Calflish Example

``` r
testFn7 <- function(x) {
    n <- length(x)
    p <- 1/n
    (1 + p)^n * prod(x^p)
}
testFn7_v <- function(x) {
    matrix(apply(x, 2, function(z) {
        n <- length(z)
        p <- 1/n
        (1 + p)^n * prod(z^p)
    }), ncol = ncol(x))
}

d <- harness(f = testFn7, fv = testFn7_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 20)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  3.36ms |  3.43ms | 272.514 |   32.89KB | 27.251 |    20 |    2 |     73.4ms |
| Vectorized Hcubature           |  4.07ms |  4.21ms | 224.062 |  205.61KB | 22.406 |    20 |    2 |     89.3ms |
| Non-vectorized Pcubature       |  8.41ms |  8.68ms | 108.359 |        0B | 27.090 |    20 |    5 |    184.6ms |
| Vectorized Pcubature           |  9.74ms |  9.96ms |  94.818 |  386.24KB | 23.704 |    20 |    5 |    210.9ms |
| Non-vectorized cubature::cuhre | 44.48ms | 45.83ms |  21.790 |        0B | 23.969 |    20 |   22 |    917.8ms |
| Vectorized cubature::cuhre     | 43.99ms | 45.22ms |  22.085 |    2.04MB | 23.189 |    20 |   21 |    905.6ms |

## Wang-Landau Sampling 1d, 2d Examples

``` r
I.1d <- function(x) {
    sin(4 * x) *
        x * ((x * ( x * (x * x - 4) + 1) - 1))
}
I.1d_v <- function(x) {
    matrix(apply(x, 2, function(z)
        sin(4 * z) *
        z * ((z * ( z * (z * z - 4) + 1) - 1))),
        ncol = ncol(x))
}
d <- harness(f = I.1d, fv = I.1d_v,
             lowerLimit = -2, upperLimit = 2, iterations = 100)
knitr::kable(d, digits = 3)
```

| expression                     |     min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |   126µs | 137.1µs |  7007.805 |    53.7KB | 70.786 |    99 |    1 |    14.13ms |
| Vectorized Hcubature           | 211.1µs | 220.6µs |  4335.328 |    68.5KB |  0.000 |   100 |    0 |    23.07ms |
| Non-vectorized Pcubature       |  51.6µs |  54.4µs | 17695.283 |        0B |  0.000 |   100 |    0 |     5.65ms |
| Vectorized Pcubature           | 154.9µs | 164.3µs |  5727.471 |        0B | 57.853 |    99 |    1 |    17.29ms |
| Non-vectorized cubature::cuhre | 220.1µs | 228.4µs |  4318.719 |        0B |  0.000 |   100 |    0 |    23.16ms |
| Vectorized cubature::cuhre     | 502.1µs | 527.1µs |  1869.066 |        0B | 18.879 |    99 |    1 |    52.97ms |

``` r
I.2d <- function(x) {
    x1 <- x[1]; x2 <- x[2]
    sin(4 * x1 + 1) * cos(4 * x2) * x1 * (x1 * (x1 * x1)^2 - x2 * (x2 * x2 - x1) +2)
}
I.2d_v <- function(x) {
    matrix(apply(x, 2,
                 function(z) {
                     x1 <- z[1]; x2 <- z[2]
                     sin(4 * x1 + 1) * cos(4 * x2) * x1 * (x1 * (x1 * x1)^2 - x2 * (x2 * x2 - x1) +2)
                 }),
           ncol = ncol(x))
}
d <- harness(f = I.2d, fv = I.2d_v,
             lowerLimit = rep(-1, 2), upperLimit = rep(1, 2), iterations = 100)
knitr::kable(d, digits = 3)
```

| expression                     |      min |   median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |   4.48ms |   4.56ms |  216.767 |    78.4KB | 35.288 |    86 |   14 |    396.7ms |
| Vectorized Hcubature           |   5.47ms |    5.6ms |  177.918 |   304.7KB | 26.585 |    87 |   13 |      489ms |
| Non-vectorized Pcubature       | 408.42µs |  424.5µs | 2317.877 |        0B | 23.413 |    99 |    1 |     42.7ms |
| Vectorized Pcubature           | 614.52µs | 640.81µs | 1533.359 |    18.3KB | 31.293 |    98 |    2 |     63.9ms |
| Non-vectorized cubature::cuhre |   1.24ms |   1.27ms |  777.161 |        0B | 24.036 |    97 |    3 |    124.8ms |
| Vectorized cubature::cuhre     |   1.35ms |    1.4ms |  708.661 |    60.1KB | 21.917 |    97 |    3 |    136.9ms |

## Implementation Notes

About the only real modification we have made to the underlying
`cubature` library is that we use `M = 16` rather than the default
`M = 19` suggested by the original author for `pcubature`. This allows
us to comply with CRAN package size limits and seems to work reasonably
well for the above tests. Future versions will allow for such
customization on demand.

The changes made to the `Cuba` library are managed in a [Github repo
branch](https://github.com/bnaras/Cuba/tree/R_pkg): each time a new
release is made, we update the main branch, and keep all changes for
Unix platforms in a branch named `R_pkg` against the current main
branch. Customization for windows is done in the package itself using
the `Makevars.win` script.

## Summary

The recommendations are:

1.  *Vectorize* your function. The time spent in so doing pays back
    enormously. This is easy to do and the examples above show how.

2.  Vectorized `hcubature` seems to be a good starting point.

3.  For smooth integrands in low dimensions ($\leq 3$), `pcubature`
    might be worth trying out. Experiment before using in a production
    package.

## Session Info

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] mvtnorm_1.3-3    cubature_2.1.4-1 bench_1.1.4     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.5       cli_3.6.5         knitr_1.50        rlang_1.1.6      
    ##  [5] xfun_0.54         textshaping_1.0.4 jsonlite_2.0.0    glue_1.8.0       
    ##  [9] htmltools_0.5.8.1 ragg_1.5.0        sass_0.4.10       rmarkdown_2.30   
    ## [13] tibble_3.3.0      evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0    
    ## [17] profmem_0.7.0     yaml_2.3.10       lifecycle_1.0.4   compiler_4.5.2   
    ## [21] fs_1.6.6          pkgconfig_2.0.3   Rcpp_1.1.0        systemfonts_1.3.1
    ## [25] digest_0.6.39     R6_2.6.1          pillar_1.11.1     magrittr_2.0.4   
    ## [29] bslib_0.9.0       tools_4.5.2       pkgdown_2.2.0     cachem_1.1.0     
    ## [33] desc_1.4.3
