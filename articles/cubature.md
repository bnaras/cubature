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
| Non-vectorized Hcubature       | 272.62µs | 284.25µs | 3467.632 |    59.3KB | 35.027 |    99 |    1 |     28.6ms |
| Vectorized Hcubature           |    395µs | 418.68µs | 2357.743 |    54.5KB | 23.816 |    99 |    1 |       42ms |
| Non-vectorized Pcubature       | 872.97µs | 896.86µs | 1099.629 |    31.2KB | 34.009 |    97 |    3 |     88.2ms |
| Vectorized Pcubature           |   1.21ms |   1.24ms |  795.401 |    58.4KB | 16.233 |    98 |    2 |    123.2ms |
| Non-vectorized cubature::cuhre | 591.12µs | 607.86µs | 1633.375 |    40.5KB | 16.499 |    99 |    1 |     60.6ms |
| Vectorized cubature::cuhre     |  601.8µs | 637.24µs | 1566.886 |    39.8KB | 31.977 |    98 |    2 |     62.5ms |

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
| Non-vectorized Hcubature       | 845.58ms | 851.77ms |   1.151 |  139.68KB | 17.844 |    10 |  155 |      8.69s |
| Vectorized Hcubature           |   1.82ms |   2.55ms | 423.708 |    1.89MB |  0.000 |    10 |    0 |     23.6ms |
| Non-vectorized Pcubature       | 363.38ms | 366.59ms |   2.679 |        0B | 17.416 |    10 |   65 |      3.73s |
| Vectorized Pcubature           |   1.28ms |   1.32ms | 626.560 |  810.26KB | 62.656 |    10 |    1 |    15.96ms |
| Non-vectorized cubature::cuhre | 352.33ms | 353.85ms |   2.822 |        0B | 17.778 |    10 |   63 |      3.54s |
| Vectorized cubature::cuhre     |   3.34ms |   3.39ms | 292.965 |  898.41KB |  0.000 |    10 |    0 |    34.13ms |

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
| g1()       | 767µs | 1.38ms | 731.020 |  355.13KB |      0 |    20 |    0 |     27.4ms |
| g2()       | 748µs | 1.37ms | 757.928 |    2.49KB |      0 |    20 |    0 |     26.4ms |
| g3()       | 762µs | 1.38ms | 734.542 |    2.49KB |      0 |    20 |    0 |     27.2ms |

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
| Non-vectorized Hcubature       |  36.6µs |  39.4µs | 24555.168 |    7.66KB | 24.580 |   999 |    1 |     40.7ms |
| Vectorized Hcubature           |  65.3µs |  68.8µs | 13930.744 |    1.16KB | 27.917 |   998 |    2 |     71.6ms |
| Non-vectorized Pcubature       |    49µs |  51.2µs | 18906.070 |        0B | 18.925 |   999 |    1 |     52.8ms |
| Vectorized Pcubature           | 114.4µs | 121.8µs |  7967.540 |   18.68KB | 23.975 |   997 |    3 |    125.1ms |
| Non-vectorized cubature::cuhre | 341.3µs |   356µs |  2787.078 |        0B | 25.312 |   991 |    9 |    355.6ms |
| Vectorized cubature::cuhre     |   372µs | 395.4µs |  2514.324 |   16.38KB | 25.397 |   990 |   10 |    393.7ms |

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

| expression                     |    min |   median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|---------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  2.9ms |   2.95ms |   336.323 |    67.3KB |  0.000 |    10 |    0 |    29.73ms |
| Vectorized Hcubature           |  5.1ms |   5.17ms |   192.985 |   290.5KB | 48.246 |     8 |    2 |    41.45ms |
| Non-vectorized Pcubature       | 74.9µs |  77.09µs | 11727.996 |        0B |  0.000 |    10 |    0 |   852.66µs |
| Vectorized Pcubature           |  159µs | 165.87µs |  5736.449 |     4.1KB |  0.000 |    10 |    0 |     1.74ms |
| Non-vectorized cubature::cuhre | 13.8ms |  13.96ms |    71.547 |        0B | 30.663 |     7 |    3 |    97.84ms |
| Vectorized cubature::cuhre     | 21.2ms |  21.24ms |    46.996 |   971.5KB | 46.996 |     5 |    5 |   106.39ms |

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
| Non-vectorized Hcubature       |  48.9ms |  49.2ms |  17.149 |    17.7KB | 34.298 |    10 |   20 |   583.13ms |
| Vectorized Hcubature           |  49.9ms |  50.6ms |  19.628 |  1011.8KB | 23.553 |    10 |   12 |   509.49ms |
| Non-vectorized cubature::cuhre | 828.8ms | 839.2ms |   1.183 |        0B | 27.923 |    10 |  236 |      8.45s |
| Vectorized cubature::cuhre     | 940.1ms | 955.3ms |   1.037 |    21.2MB | 23.544 |    10 |  227 |      9.64s |

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
| Non-vectorized Hcubature       |  57.7µs |  60.9µs | 15522.842 |    6.75KB |   0.000 |    50 |    0 |     3.22ms |
| Vectorized Hcubature           |  91.8µs |  95.8µs | 10003.803 |    2.91KB |   0.000 |    50 |    0 |        5ms |
| Non-vectorized Pcubature       |    51µs |  53.9µs | 17386.257 |        0B |   0.000 |    50 |    0 |     2.88ms |
| Vectorized Pcubature           |  82.7µs |  89.4µs | 10289.507 |   19.12KB | 209.990 |    49 |    1 |     4.76ms |
| Non-vectorized cubature::cuhre | 651.4µs | 678.3µs |  1473.807 |        0B |  30.078 |    49 |    1 |    33.25ms |
| Vectorized cubature::cuhre     | 666.1µs | 703.1µs |  1422.627 |   39.84KB |   0.000 |    50 |    0 |    35.15ms |

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
| Non-vectorized Hcubature       | 1.37ms |  1.4ms | 708.239 |    85.2KB | 37.276 |    19 |    1 |     26.8ms |
| Vectorized Hcubature           | 1.84ms |  1.9ms | 527.843 |   147.2KB | 27.781 |    19 |    1 |       36ms |
| Non-vectorized Pcubature       | 2.04ms | 2.08ms | 477.612 |        0B | 25.137 |    19 |    1 |     39.8ms |
| Vectorized Pcubature           | 2.67ms | 2.78ms | 361.781 |    68.9KB | 40.198 |    18 |    2 |     49.8ms |
| Non-vectorized cubature::cuhre | 3.26ms | 3.33ms | 298.684 |        0B | 15.720 |    19 |    1 |     63.6ms |
| Vectorized cubature::cuhre     | 3.88ms | 3.93ms | 253.759 |   125.6KB | 28.195 |    18 |    2 |     70.9ms |

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
| Non-vectorized Hcubature       | 3.69ms | 3.74ms | 265.373 |     133KB | 46.830 |    17 |    3 |     64.1ms |
| Vectorized Hcubature           | 4.82ms | 4.92ms | 203.239 |     249KB | 22.582 |    18 |    2 |     88.6ms |
| Non-vectorized Pcubature       | 2.55ms | 2.59ms | 385.540 |        0B | 42.838 |    18 |    2 |     46.7ms |
| Vectorized Pcubature           | 3.37ms | 3.46ms | 288.877 |      69KB | 15.204 |    19 |    1 |     65.8ms |
| Non-vectorized cubature::cuhre | 7.15ms | 7.21ms | 138.804 |        0B | 46.268 |    15 |    5 |    108.1ms |
| Vectorized cubature::cuhre     | 8.46ms | 8.61ms | 116.132 |     224KB | 20.494 |    17 |    3 |    146.4ms |

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
| Non-vectorized Hcubature       |  1.61ms |  1.63ms | 613.061 |    64.6KB | 32.266 |    19 |    1 |       31ms |
| Vectorized Hcubature           |  2.18ms |  2.24ms | 446.971 |     156KB |  0.000 |    20 |    0 |     44.7ms |
| Non-vectorized Pcubature       |  8.23ms |  8.37ms | 119.499 |        0B | 39.833 |    15 |    5 |    125.5ms |
| Vectorized Pcubature           | 10.45ms | 10.72ms |  93.443 |   386.2KB | 31.148 |    15 |    5 |    160.5ms |
| Non-vectorized cubature::cuhre |  4.33ms |  4.42ms | 226.651 |        0B | 25.183 |    18 |    2 |     79.4ms |
| Vectorized cubature::cuhre     |  4.87ms |  5.01ms | 199.387 |   225.8KB | 22.154 |    18 |    2 |     90.3ms |

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
| Non-vectorized Hcubature       |  3.39ms |  3.48ms | 267.142 |   32.89KB | 26.714 |    20 |    2 |     74.9ms |
| Vectorized Hcubature           |   4.1ms |  4.24ms | 220.062 |  205.61KB | 22.006 |    20 |    2 |     90.9ms |
| Non-vectorized Pcubature       |  8.48ms |  8.74ms | 106.842 |        0B | 26.711 |    20 |    5 |    187.2ms |
| Vectorized Pcubature           |  9.87ms | 10.04ms |  93.003 |  386.24KB | 23.251 |    20 |    5 |      215ms |
| Non-vectorized cubature::cuhre | 44.53ms | 45.16ms |  21.974 |        0B | 24.172 |    20 |   22 |    910.2ms |
| Vectorized cubature::cuhre     | 45.36ms | 46.02ms |  21.593 |    2.04MB | 22.673 |    20 |   21 |    926.2ms |

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
| Non-vectorized Hcubature       |   127µs | 138.1µs |  6996.101 |    53.7KB | 70.668 |    99 |    1 |    14.15ms |
| Vectorized Hcubature           | 212.9µs | 225.3µs |  4290.781 |    68.5KB |  0.000 |   100 |    0 |    23.31ms |
| Non-vectorized Pcubature       |  51.8µs |  54.1µs | 17704.325 |        0B |  0.000 |   100 |    0 |     5.65ms |
| Vectorized Pcubature           | 154.8µs | 165.3µs |  5643.709 |        0B | 57.007 |    99 |    1 |    17.54ms |
| Non-vectorized cubature::cuhre | 219.8µs | 229.7µs |  4280.917 |        0B |  0.000 |   100 |    0 |    23.36ms |
| Vectorized cubature::cuhre     | 502.7µs | 530.4µs |  1848.474 |        0B | 18.671 |    99 |    1 |    53.56ms |

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
| Non-vectorized Hcubature       |   4.56ms |   4.66ms |  210.996 |    78.4KB | 34.348 |    86 |   14 |    407.6ms |
| Vectorized Hcubature           |   5.54ms |   5.72ms |  173.807 |   304.7KB | 25.971 |    87 |   13 |    500.6ms |
| Non-vectorized Pcubature       | 408.58µs | 432.05µs | 2272.399 |        0B | 22.954 |    99 |    1 |     43.6ms |
| Vectorized Pcubature           | 617.04µs | 650.51µs | 1496.325 |    18.3KB | 30.537 |    98 |    2 |     65.5ms |
| Non-vectorized cubature::cuhre |   1.25ms |    1.3ms |  752.656 |        0B | 23.278 |    97 |    3 |    128.9ms |
| Vectorized cubature::cuhre     |   1.36ms |    1.4ms |  698.721 |    60.1KB | 21.610 |    97 |    3 |    138.8ms |

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
    ## [1] mvtnorm_1.3-3    cubature_2.1.4-2 bench_1.1.4     
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
