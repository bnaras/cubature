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
| Non-vectorized Hcubature       | 252.74µs | 281.93µs | 3516.802 |    59.3KB | 35.523 |    99 |    1 |     28.2ms |
| Vectorized Hcubature           | 371.74µs | 392.53µs | 2511.624 |    54.5KB | 25.370 |    99 |    1 |     39.4ms |
| Non-vectorized Pcubature       | 836.09µs | 899.53µs | 1107.845 |    31.2KB | 34.263 |    97 |    3 |     87.6ms |
| Vectorized Pcubature           |   1.14ms |   1.19ms |  832.897 |    58.4KB | 16.998 |    98 |    2 |    117.7ms |
| Non-vectorized cubature::cuhre | 540.38µs | 588.11µs | 1692.355 |    40.5KB | 17.094 |    99 |    1 |     58.5ms |
| Vectorized cubature::cuhre     | 575.57µs | 610.75µs | 1627.377 |    39.8KB | 33.212 |    98 |    2 |     60.2ms |

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
| Non-vectorized Hcubature       | 673.62ms | 683.49ms |   1.422 |  139.68KB | 22.036 |    10 |  155 |      7.03s |
| Vectorized Hcubature           |   1.68ms |   2.32ms | 463.807 |    1.89MB |  0.000 |    10 |    0 |    21.56ms |
| Non-vectorized Pcubature       | 287.38ms |  295.2ms |   3.332 |        0B | 21.661 |    10 |   65 |         3s |
| Vectorized Pcubature           |   1.11ms |   1.22ms | 692.083 |  810.26KB | 69.208 |    10 |    1 |    14.45ms |
| Non-vectorized cubature::cuhre | 282.43ms |  284.3ms |   3.514 |        0B | 22.137 |    10 |   63 |      2.85s |
| Vectorized cubature::cuhre     |   2.65ms |   2.79ms | 350.634 |  898.41KB |  0.000 |    10 |    0 |    28.52ms |

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
| g1()       | 675µs | 1.24ms | 809.138 |  355.13KB |      0 |    20 |    0 |     24.7ms |
| g2()       | 670µs | 1.24ms | 831.680 |    2.49KB |      0 |    20 |    0 |       24ms |
| g3()       | 680µs | 1.24ms | 811.326 |    2.49KB |      0 |    20 |    0 |     24.7ms |

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
| Non-vectorized Hcubature       |    32µs |  35.8µs | 27305.257 |    7.66KB | 27.333 |   999 |    1 |     36.6ms |
| Vectorized Hcubature           |  54.2µs |  58.4µs | 16702.135 |    1.16KB | 33.471 |   998 |    2 |     59.8ms |
| Non-vectorized Pcubature       |  42.2µs |  46.4µs | 21025.150 |        0B | 21.046 |   999 |    1 |     47.5ms |
| Vectorized Pcubature           |  96.2µs | 103.6µs |  9481.220 |   18.68KB | 28.529 |   997 |    3 |    105.2ms |
| Non-vectorized cubature::cuhre | 303.4µs | 348.4µs |  2863.060 |        0B | 26.002 |   991 |    9 |    346.1ms |
| Vectorized cubature::cuhre     | 332.6µs | 358.2µs |  2770.079 |   16.38KB | 27.981 |   990 |   10 |    357.4ms |

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
| Non-vectorized Hcubature       |   2.63ms |   2.67ms |   372.475 |    67.3KB | 41.386 |     9 |    1 |    24.16ms |
| Vectorized Hcubature           |   4.62ms |   4.71ms |   212.755 |   290.5KB | 23.639 |     9 |    1 |     42.3ms |
| Non-vectorized Pcubature       |  79.37µs |  81.16µs | 11662.588 |        0B |  0.000 |    10 |    0 |   857.44µs |
| Vectorized Pcubature           | 146.13µs | 151.91µs |  6288.140 |     4.1KB |  0.000 |    10 |    0 |     1.59ms |
| Non-vectorized cubature::cuhre |  12.21ms |  12.44ms |    80.191 |        0B | 34.368 |     7 |    3 |    87.29ms |
| Vectorized cubature::cuhre     |   18.9ms |     19ms |    52.502 |   971.5KB | 52.502 |     5 |    5 |    95.23ms |

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
| Non-vectorized Hcubature       |  46.5ms |  47.1ms |  18.079 |    17.7KB | 36.157 |    10 |   20 |   553.14ms |
| Vectorized Hcubature           |  47.3ms |  47.9ms |  20.651 |  1011.8KB | 26.846 |    10 |   13 |   484.25ms |
| Non-vectorized cubature::cuhre | 805.7ms | 818.3ms |   1.205 |        0B | 28.315 |    10 |  235 |       8.3s |
| Vectorized cubature::cuhre     | 872.4ms | 882.8ms |   1.120 |    21.2MB | 25.433 |    10 |  227 |      8.93s |

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
| Non-vectorized Hcubature       |  53.3µs |  57.3µs | 16759.695 |    6.75KB |   0.000 |    50 |    0 |     2.98ms |
| Vectorized Hcubature           |  79.7µs |  85.4µs | 10851.973 |    2.91KB | 221.469 |    49 |    1 |     4.51ms |
| Non-vectorized Pcubature       |    44µs |  45.1µs | 21084.959 |        0B |   0.000 |    50 |    0 |     2.37ms |
| Vectorized Pcubature           |  67.9µs |  71.4µs | 13515.659 |   19.12KB |   0.000 |    50 |    0 |      3.7ms |
| Non-vectorized cubature::cuhre | 596.4µs | 641.5µs |  1538.272 |        0B |  31.393 |    49 |    1 |    31.85ms |
| Vectorized cubature::cuhre     | 612.9µs | 652.7µs |  1531.361 |   39.84KB |  31.252 |    49 |    1 |       32ms |

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
| Non-vectorized Hcubature       | 1.23ms | 1.31ms | 757.544 |    85.2KB | 39.871 |    19 |    1 |     25.1ms |
| Vectorized Hcubature           |  1.7ms | 1.76ms | 568.723 |   147.2KB |  0.000 |    20 |    0 |     35.2ms |
| Non-vectorized Pcubature       | 1.85ms | 1.99ms | 496.789 |        0B | 55.199 |    18 |    2 |     36.2ms |
| Vectorized Pcubature           | 2.45ms | 2.55ms | 390.605 |    68.9KB | 20.558 |    19 |    1 |     48.6ms |
| Non-vectorized cubature::cuhre | 2.95ms | 3.12ms | 318.540 |        0B | 35.393 |    18 |    2 |     56.5ms |
| Vectorized cubature::cuhre     | 3.52ms | 3.62ms | 275.377 |   125.6KB | 30.597 |    18 |    2 |     65.4ms |

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
| Non-vectorized Hcubature       | 3.32ms | 3.59ms | 281.564 |     133KB | 31.285 |    18 |    2 |     63.9ms |
| Vectorized Hcubature           | 4.58ms | 4.78ms | 210.121 |     249KB | 23.347 |    18 |    2 |     85.7ms |
| Non-vectorized Pcubature       | 2.31ms | 2.39ms | 412.413 |        0B | 45.824 |    18 |    2 |     43.6ms |
| Vectorized Pcubature           | 3.13ms | 3.28ms | 304.179 |      69KB | 33.798 |    18 |    2 |     59.2ms |
| Non-vectorized cubature::cuhre | 6.64ms | 6.75ms | 147.482 |        0B | 36.870 |    16 |    4 |    108.5ms |
| Vectorized cubature::cuhre     |  7.9ms | 7.99ms | 124.664 |     224KB | 31.166 |    16 |    4 |    128.3ms |

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

| expression                     |    min |  median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|--------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 1.51ms |  1.59ms | 625.555 |    64.6KB | 32.924 |    19 |    1 |     30.4ms |
| Vectorized Hcubature           | 2.03ms |  2.05ms | 483.449 |     156KB | 25.445 |    19 |    1 |     39.3ms |
| Non-vectorized Pcubature       | 7.63ms |  7.74ms | 127.032 |        0B | 42.344 |    15 |    5 |    118.1ms |
| Vectorized Pcubature           | 9.96ms | 10.05ms |  99.263 |   386.2KB | 24.816 |    16 |    4 |    161.2ms |
| Non-vectorized cubature::cuhre |    4ms |  4.18ms | 239.557 |        0B | 42.275 |    17 |    3 |       71ms |
| Vectorized cubature::cuhre     |  4.5ms |  4.57ms | 219.081 |   225.8KB | 24.342 |    18 |    2 |     82.2ms |

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
| Non-vectorized Hcubature       |  3.11ms |  3.35ms | 271.504 |   32.89KB | 27.150 |    20 |    2 |     73.7ms |
| Vectorized Hcubature           |  3.85ms |  3.96ms | 233.513 |  205.61KB | 23.351 |    20 |    2 |     85.6ms |
| Non-vectorized Pcubature       |   7.8ms |  8.64ms | 106.683 |        0B | 32.005 |    20 |    6 |    187.5ms |
| Vectorized Pcubature           |  9.08ms |  9.43ms |  97.782 |  386.24KB | 24.446 |    20 |    5 |    204.5ms |
| Non-vectorized cubature::cuhre | 42.05ms |  43.1ms |  22.738 |        0B | 23.875 |    20 |   21 |    879.6ms |
| Vectorized cubature::cuhre     | 42.43ms | 43.13ms |  23.024 |    2.04MB | 25.327 |    20 |   22 |    868.6ms |

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

| expression                     |     min |  median |   itr/sec | mem_alloc |  gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|--------:|----------:|----------:|--------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 112.9µs | 119.4µs |  7949.969 |    53.7KB |   0.000 |   100 |    0 |    12.58ms |
| Vectorized Hcubature           | 194.3µs | 206.6µs |  4761.105 |    68.5KB |   0.000 |   100 |    0 |       21ms |
| Non-vectorized Pcubature       |  49.2µs |  56.8µs | 16143.375 |        0B | 163.064 |    99 |    1 |     6.13ms |
| Vectorized Pcubature           | 131.8µs | 136.7µs |  6932.497 |        0B |   0.000 |   100 |    0 |    14.43ms |
| Non-vectorized cubature::cuhre | 200.8µs | 224.8µs |  4303.243 |        0B |  43.467 |    99 |    1 |    23.01ms |
| Vectorized cubature::cuhre     | 457.3µs | 487.8µs |  1988.945 |        0B |  20.090 |    99 |    1 |    49.77ms |

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
| Non-vectorized Hcubature       |   4.32ms |   4.48ms |  218.946 |    78.4KB | 35.642 |    86 |   14 |    392.8ms |
| Vectorized Hcubature           |   5.46ms |    5.6ms |  178.056 |   304.7KB | 26.606 |    87 |   13 |    488.6ms |
| Non-vectorized Pcubature       | 373.51µs | 420.07µs | 2308.003 |        0B | 47.102 |    98 |    2 |     42.5ms |
| Vectorized Pcubature           | 580.31µs | 611.19µs | 1619.748 |    18.3KB | 16.361 |    99 |    1 |     61.1ms |
| Non-vectorized cubature::cuhre |   1.15ms |   1.19ms |  825.542 |        0B | 25.532 |    97 |    3 |    117.5ms |
| Vectorized cubature::cuhre     |   1.29ms |   1.33ms |  745.039 |    60.1KB | 23.042 |    97 |    3 |    130.2ms |

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
    ## [25] digest_0.6.38     R6_2.6.1          pillar_1.11.1     magrittr_2.0.4   
    ## [29] bslib_0.9.0       tools_4.5.2       pkgdown_2.2.0     cachem_1.1.0     
    ## [33] desc_1.4.3
