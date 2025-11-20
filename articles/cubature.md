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
| Non-vectorized Hcubature       | 269.88µs |    282µs | 3531.961 |    59.3KB | 35.676 |    99 |    1 |       28ms |
| Vectorized Hcubature           |  391.4µs | 408.08µs | 2418.316 |    54.5KB | 24.427 |    99 |    1 |     40.9ms |
| Non-vectorized Pcubature       | 853.34µs | 877.23µs | 1122.984 |    31.2KB | 34.731 |    97 |    3 |     86.4ms |
| Vectorized Pcubature           |   1.18ms |   1.21ms |  819.137 |    58.4KB | 16.717 |    98 |    2 |    119.6ms |
| Non-vectorized cubature::cuhre | 578.79µs | 599.58µs | 1657.366 |    40.5KB | 16.741 |    99 |    1 |     59.7ms |
| Vectorized cubature::cuhre     | 589.67µs | 630.46µs | 1590.675 |    39.8KB | 32.463 |    98 |    2 |     61.6ms |

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
| Non-vectorized Hcubature       | 806.67ms | 814.05ms |   1.213 |  139.68KB | 18.795 |    10 |  155 |      8.25s |
| Vectorized Hcubature           |   1.77ms |   2.45ms | 441.859 |    1.89MB |  0.000 |    10 |    0 |    22.63ms |
| Non-vectorized Pcubature       | 345.66ms | 349.81ms |   2.815 |        0B | 18.299 |    10 |   65 |      3.55s |
| Vectorized Pcubature           |   1.25ms |   1.28ms | 660.090 |  810.26KB | 66.009 |    10 |    1 |    15.15ms |
| Non-vectorized cubature::cuhre | 334.58ms | 338.71ms |   2.940 |        0B | 18.522 |    10 |   63 |       3.4s |
| Vectorized cubature::cuhre     |   3.18ms |   3.22ms | 308.492 |  898.41KB |  0.000 |    10 |    0 |    32.42ms |

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
| g1()       | 771µs | 1.38ms | 730.011 |  355.13KB |      0 |    20 |    0 |     27.4ms |
| g2()       | 760µs | 1.37ms | 762.057 |    2.49KB |      0 |    20 |    0 |     26.2ms |
| g3()       | 760µs | 1.37ms | 741.999 |    2.49KB |      0 |    20 |    0 |       27ms |

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
| Non-vectorized Hcubature       |  36.7µs |  38.6µs | 25216.742 |    7.66KB | 25.242 |   999 |    1 |     39.6ms |
| Vectorized Hcubature           |  63.7µs |  66.7µs | 14566.850 |    1.16KB | 29.192 |   998 |    2 |     68.5ms |
| Non-vectorized Pcubature       |  48.5µs |  50.1µs | 19428.731 |        0B | 19.448 |   999 |    1 |     51.4ms |
| Vectorized Pcubature           | 114.2µs | 118.9µs |  8200.614 |   18.68KB | 24.676 |   997 |    3 |    121.6ms |
| Non-vectorized cubature::cuhre | 337.4µs | 349.4µs |  2836.483 |        0B | 25.760 |   991 |    9 |    349.4ms |
| Vectorized cubature::cuhre     | 362.8µs |   381µs |  2601.642 |   16.38KB | 26.279 |   990 |   10 |    380.5ms |

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
| Non-vectorized Hcubature       |    2.9ms |   2.93ms |   340.402 |    67.3KB |  0.000 |    10 |    0 |     29.4ms |
| Vectorized Hcubature           |   4.96ms |   5.04ms |   198.613 |   290.5KB | 49.653 |     8 |    2 |     40.3ms |
| Non-vectorized Pcubature       |  73.29µs |  74.26µs | 12236.013 |        0B |  0.000 |    10 |    0 |    817.3µs |
| Vectorized Pcubature           | 151.08µs | 156.08µs |  6241.746 |     4.1KB |  0.000 |    10 |    0 |      1.6ms |
| Non-vectorized cubature::cuhre |  13.49ms |  13.57ms |    73.359 |        0B | 31.440 |     7 |    3 |     95.4ms |
| Vectorized cubature::cuhre     |  20.34ms |  20.56ms |    47.055 |   971.5KB | 47.055 |     5 |    5 |    106.3ms |

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
| Non-vectorized Hcubature       |  44.8ms |  45.9ms |  18.731 |    17.7KB | 37.463 |    10 |   20 |   533.87ms |
| Vectorized Hcubature           |  48.7ms |  49.2ms |  20.085 |  1011.8KB | 26.110 |    10 |   13 |    497.9ms |
| Non-vectorized cubature::cuhre | 804.5ms | 824.5ms |   1.213 |        0B | 28.498 |    10 |  235 |      8.25s |
| Vectorized cubature::cuhre     | 897.7ms | 910.1ms |   1.086 |    21.2MB | 24.643 |    10 |  227 |      9.21s |

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
| Non-vectorized Hcubature       |  57.8µs |  59.1µs | 16446.168 |    6.75KB |   0.000 |    50 |    0 |     3.04ms |
| Vectorized Hcubature           |  90.6µs |  96.5µs | 10128.527 |    2.91KB |   0.000 |    50 |    0 |     4.94ms |
| Non-vectorized Pcubature       |  49.7µs |  53.5µs | 17360.574 |        0B | 354.297 |    49 |    1 |     2.82ms |
| Vectorized Pcubature           |  80.7µs |  83.7µs | 11495.409 |   19.12KB |   0.000 |    50 |    0 |     4.35ms |
| Non-vectorized cubature::cuhre | 651.2µs | 671.8µs |  1467.516 |        0B |  29.949 |    49 |    1 |    33.39ms |
| Vectorized cubature::cuhre     | 657.7µs | 684.5µs |  1451.110 |   39.84KB |   0.000 |    50 |    0 |    34.46ms |

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
| Non-vectorized Hcubature       | 1.37ms | 1.39ms | 715.486 |    85.2KB | 37.657 |    19 |    1 |     26.6ms |
| Vectorized Hcubature           |  1.8ms | 1.83ms | 541.295 |   147.2KB | 28.489 |    19 |    1 |     35.1ms |
| Non-vectorized Pcubature       | 2.05ms | 2.08ms | 479.085 |        0B | 25.215 |    19 |    1 |     39.7ms |
| Vectorized Pcubature           | 2.62ms | 2.67ms | 374.333 |    68.9KB | 41.593 |    18 |    2 |     48.1ms |
| Non-vectorized cubature::cuhre | 3.23ms | 3.31ms | 303.096 |        0B | 33.677 |    18 |    2 |     59.4ms |
| Vectorized cubature::cuhre     | 3.77ms | 3.86ms | 259.156 |   125.6KB | 13.640 |    19 |    1 |     73.3ms |

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
| Non-vectorized Hcubature       | 3.58ms |  3.6ms | 275.429 |     133KB | 30.603 |    18 |    2 |     65.4ms |
| Vectorized Hcubature           | 4.61ms | 4.72ms | 211.189 |     249KB | 37.269 |    17 |    3 |     80.5ms |
| Non-vectorized Pcubature       | 2.51ms | 2.56ms | 390.582 |        0B | 43.398 |    18 |    2 |     46.1ms |
| Vectorized Pcubature           | 3.29ms | 3.35ms | 298.252 |      69KB | 15.697 |    19 |    1 |     63.7ms |
| Non-vectorized cubature::cuhre | 7.57ms | 7.67ms | 130.448 |        0B | 32.612 |    16 |    4 |    122.7ms |
| Vectorized cubature::cuhre     | 8.21ms | 8.28ms | 120.104 |     224KB | 30.026 |    16 |    4 |    133.2ms |

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
| Non-vectorized Hcubature       |   1.6ms |  1.62ms | 616.365 |    64.6KB | 32.440 |    19 |    1 |     30.8ms |
| Vectorized Hcubature           |  2.16ms |  2.18ms | 457.225 |     156KB | 24.064 |    19 |    1 |     41.6ms |
| Non-vectorized Pcubature       |  8.27ms |  8.38ms | 119.182 |        0B | 39.727 |    15 |    5 |    125.9ms |
| Vectorized Pcubature           | 10.47ms | 10.66ms |  93.949 |   386.2KB | 31.316 |    15 |    5 |    159.7ms |
| Non-vectorized cubature::cuhre |  4.34ms |  4.43ms | 226.439 |        0B | 25.160 |    18 |    2 |     79.5ms |
| Vectorized cubature::cuhre     |  4.75ms |  4.84ms | 205.730 |   225.8KB | 22.859 |    18 |    2 |     87.5ms |

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
| Non-vectorized Hcubature       |   3.4ms |  3.48ms | 271.922 |   32.89KB | 27.192 |    20 |    2 |     73.5ms |
| Vectorized Hcubature           |  4.05ms |  4.17ms | 229.240 |  205.61KB | 22.924 |    20 |    2 |     87.2ms |
| Non-vectorized Pcubature       |  8.51ms |  8.71ms | 109.699 |        0B | 27.425 |    20 |    5 |    182.3ms |
| Vectorized Pcubature           |  9.56ms |  9.82ms |  97.716 |  386.24KB | 24.429 |    20 |    5 |    204.7ms |
| Non-vectorized cubature::cuhre | 43.18ms | 43.76ms |  22.697 |        0B | 24.967 |    20 |   22 |    881.2ms |
| Vectorized cubature::cuhre     | 43.52ms | 44.03ms |  22.633 |    2.04MB | 23.765 |    20 |   21 |    883.7ms |

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
| Non-vectorized Hcubature       | 124.5µs | 128.6µs |  7608.038 |    53.7KB |  0.000 |   100 |    0 |    13.14ms |
| Vectorized Hcubature           | 209.1µs | 222.6µs |  4384.940 |    68.5KB | 44.292 |    99 |    1 |    22.58ms |
| Non-vectorized Pcubature       |  50.5µs |  51.9µs | 18556.261 |        0B |  0.000 |   100 |    0 |     5.39ms |
| Vectorized Pcubature           | 153.1µs |   158µs |  6111.214 |        0B |  0.000 |   100 |    0 |    16.36ms |
| Non-vectorized cubature::cuhre | 220.4µs |   227µs |  4346.414 |        0B | 43.903 |    99 |    1 |    22.78ms |
| Vectorized cubature::cuhre     | 494.6µs | 515.9µs |  1922.005 |        0B | 19.414 |    99 |    1 |    51.51ms |

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
| Non-vectorized Hcubature       |   4.52ms |   4.57ms |  218.470 |    78.4KB | 35.565 |    86 |   14 |    393.6ms |
| Vectorized Hcubature           |   5.41ms |    5.5ms |  181.037 |   304.7KB | 27.052 |    87 |   13 |    480.6ms |
| Non-vectorized Pcubature       | 408.26µs | 425.56µs | 2335.780 |        0B | 23.594 |    99 |    1 |     42.4ms |
| Vectorized Pcubature           | 607.78µs | 633.71µs | 1567.523 |    18.3KB | 15.834 |    99 |    1 |     63.2ms |
| Non-vectorized cubature::cuhre |   1.23ms |   1.25ms |  794.719 |        0B | 33.113 |    96 |    4 |    120.8ms |
| Vectorized cubature::cuhre     |   1.33ms |   1.38ms |  709.631 |    60.1KB | 21.947 |    97 |    3 |    136.7ms |

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
    ## [1] mvtnorm_1.3-3    cubature_2.1.4-5 bench_1.1.4     
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
