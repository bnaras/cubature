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
| Non-vectorized Hcubature       | 270.25µs | 284.07µs | 3476.469 |    59.3KB | 35.116 |    99 |    1 |     28.5ms |
| Vectorized Hcubature           | 387.82µs | 412.92µs | 2407.100 |    54.5KB | 24.314 |    99 |    1 |     41.1ms |
| Non-vectorized Pcubature       | 844.24µs | 874.87µs | 1127.413 |    31.2KB | 34.868 |    97 |    3 |       86ms |
| Vectorized Pcubature           |   1.18ms |   1.21ms |  818.736 |    58.4KB | 16.709 |    98 |    2 |    119.7ms |
| Non-vectorized cubature::cuhre | 582.65µs | 600.95µs | 1649.178 |    40.5KB | 16.658 |    99 |    1 |       60ms |
| Vectorized cubature::cuhre     | 588.46µs | 631.46µs | 1584.577 |    39.8KB | 32.338 |    98 |    2 |     61.8ms |

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
| Non-vectorized Hcubature       | 827.53ms | 839.27ms |   1.177 |  139.68KB | 18.245 |    10 |  155 |       8.5s |
| Vectorized Hcubature           |   1.79ms |   2.54ms | 428.078 |    1.89MB |  0.000 |    10 |    0 |    23.36ms |
| Non-vectorized Pcubature       | 356.51ms | 361.05ms |   2.721 |        0B | 17.684 |    10 |   65 |      3.68s |
| Vectorized Pcubature           |   1.28ms |   1.31ms | 752.851 |  810.26KB |  0.000 |    10 |    0 |    13.28ms |
| Non-vectorized cubature::cuhre | 342.37ms | 346.32ms |   2.884 |        0B | 18.169 |    10 |   63 |      3.47s |
| Vectorized cubature::cuhre     |   3.27ms |   3.33ms | 298.967 |  898.41KB |  0.000 |    10 |    0 |    33.45ms |

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
| g1()       | 771µs | 1.38ms | 731.187 |  355.13KB |      0 |    20 |    0 |     27.4ms |
| g2()       | 748µs | 1.37ms | 757.342 |    2.49KB |      0 |    20 |    0 |     26.4ms |
| g3()       | 755µs | 1.38ms | 733.515 |    2.49KB |      0 |    20 |    0 |     27.3ms |

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
| Non-vectorized Hcubature       |  36.5µs |  38.9µs | 24818.854 |    7.66KB | 24.844 |   999 |    1 |     40.3ms |
| Vectorized Hcubature           |  64.4µs |  67.9µs | 14271.906 |    1.16KB | 28.601 |   998 |    2 |     69.9ms |
| Non-vectorized Pcubature       |  48.3µs |  50.2µs | 19362.142 |        0B | 19.382 |   999 |    1 |     51.6ms |
| Vectorized Pcubature           | 113.2µs | 119.5µs |  8153.994 |   18.68KB | 24.536 |   997 |    3 |    122.3ms |
| Non-vectorized cubature::cuhre | 348.4µs | 362.9µs |  2736.571 |        0B | 24.853 |   991 |    9 |    362.1ms |
| Vectorized cubature::cuhre     | 367.2µs |   386µs |  2577.510 |   16.38KB | 26.035 |   990 |   10 |    384.1ms |

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

| expression                     |      min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |    2.9ms |  2.91ms |   342.310 |    67.3KB | 38.034 |     9 |    1 |    26.29ms |
| Vectorized Hcubature           |   5.07ms |  5.12ms |   195.188 |   290.5KB | 21.688 |     9 |    1 |    46.11ms |
| Non-vectorized Pcubature       |  75.32µs | 80.41µs | 11953.006 |        0B |  0.000 |    10 |    0 |   836.61µs |
| Vectorized Pcubature           | 158.15µs | 162.7µs |  5920.752 |     4.1KB |  0.000 |    10 |    0 |     1.69ms |
| Non-vectorized cubature::cuhre |  13.66ms | 13.81ms |    72.322 |        0B | 30.995 |     7 |    3 |    96.79ms |
| Vectorized cubature::cuhre     |  20.83ms | 20.99ms |    47.631 |   971.5KB | 47.631 |     5 |    5 |   104.97ms |

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
| Non-vectorized Hcubature       |  45.3ms |  46.7ms |  18.068 |    17.7KB | 36.136 |    10 |   20 |   553.47ms |
| Vectorized Hcubature           |  48.4ms |  48.9ms |  20.241 |  1011.8KB | 26.313 |    10 |   13 |   494.05ms |
| Non-vectorized cubature::cuhre | 795.5ms | 806.8ms |   1.223 |        0B | 28.752 |    10 |  235 |      8.17s |
| Vectorized cubature::cuhre     | 900.6ms | 909.4ms |   1.086 |    21.2MB | 24.652 |    10 |  227 |      9.21s |

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
| Non-vectorized Hcubature       |  59.2µs |  61.1µs | 15720.451 |    6.75KB |   0.000 |    50 |    0 |     3.18ms |
| Vectorized Hcubature           |  93.3µs |  98.8µs |  9436.436 |    2.91KB | 192.580 |    49 |    1 |     5.19ms |
| Non-vectorized Pcubature       |  49.8µs |  51.6µs | 18848.520 |        0B |   0.000 |    50 |    0 |     2.65ms |
| Vectorized Pcubature           |  82.1µs |  84.7µs | 11393.948 |   19.12KB |   0.000 |    50 |    0 |     4.39ms |
| Non-vectorized cubature::cuhre | 672.2µs | 696.6µs |  1434.700 |        0B |  29.280 |    49 |    1 |    34.15ms |
| Vectorized cubature::cuhre     | 659.4µs | 689.6µs |  1446.316 |   39.84KB |  29.517 |    49 |    1 |    33.88ms |

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
| Non-vectorized Hcubature       | 1.36ms |  1.4ms | 709.732 |    85.2KB | 37.354 |    19 |    1 |     26.8ms |
| Vectorized Hcubature           | 1.82ms | 1.89ms | 532.534 |   147.2KB |  0.000 |    20 |    0 |     37.6ms |
| Non-vectorized Pcubature       | 2.05ms |  2.1ms | 474.668 |        0B | 52.741 |    18 |    2 |     37.9ms |
| Vectorized Pcubature           | 2.63ms |  2.7ms | 369.651 |    68.9KB | 19.455 |    19 |    1 |     51.4ms |
| Non-vectorized cubature::cuhre | 3.23ms | 3.29ms | 302.620 |        0B | 33.624 |    18 |    2 |     59.5ms |
| Vectorized cubature::cuhre     |  3.8ms | 3.88ms | 257.109 |   125.6KB | 28.568 |    18 |    2 |       70ms |

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
| Non-vectorized Hcubature       | 3.65ms |  3.7ms | 269.276 |     133KB | 29.920 |    18 |    2 |     66.8ms |
| Vectorized Hcubature           | 4.71ms | 4.79ms | 207.838 |     249KB | 23.093 |    18 |    2 |     86.6ms |
| Non-vectorized Pcubature       | 2.52ms | 2.56ms | 388.661 |        0B | 43.185 |    18 |    2 |     46.3ms |
| Vectorized Pcubature           | 3.35ms |  3.4ms | 293.714 |      69KB | 32.635 |    18 |    2 |     61.3ms |
| Non-vectorized cubature::cuhre | 6.95ms |  7.1ms | 141.114 |        0B | 35.279 |    16 |    4 |    113.4ms |
| Vectorized cubature::cuhre     | 8.28ms | 8.35ms | 119.303 |     224KB | 29.826 |    16 |    4 |    134.1ms |

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
| Non-vectorized Hcubature       |  1.56ms |   1.6ms | 622.257 |    64.6KB | 32.750 |    19 |    1 |     30.5ms |
| Vectorized Hcubature           |  2.14ms |  2.17ms | 458.951 |     156KB | 24.155 |    19 |    1 |     41.4ms |
| Non-vectorized Pcubature       |  8.16ms |  8.32ms | 120.287 |        0B | 40.096 |    15 |    5 |    124.7ms |
| Vectorized Pcubature           | 10.39ms | 10.56ms |  94.407 |   386.2KB | 23.602 |    16 |    4 |    169.5ms |
| Non-vectorized cubature::cuhre |  4.43ms |  4.49ms | 222.457 |        0B | 39.257 |    17 |    3 |     76.4ms |
| Vectorized cubature::cuhre     |  4.78ms |  4.87ms | 205.430 |   225.8KB | 22.826 |    18 |    2 |     87.6ms |

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
| Non-vectorized Hcubature       |  3.38ms |  3.49ms | 271.443 |   32.89KB | 27.144 |    20 |    2 |     73.7ms |
| Vectorized Hcubature           |  4.07ms |  4.23ms | 224.096 |  205.61KB | 22.410 |    20 |    2 |     89.2ms |
| Non-vectorized Pcubature       |  8.41ms |  8.62ms | 108.451 |        0B | 32.535 |    20 |    6 |    184.4ms |
| Vectorized Pcubature           |  9.71ms |  9.91ms |  95.325 |  386.24KB | 23.831 |    20 |    5 |    209.8ms |
| Non-vectorized cubature::cuhre | 43.28ms | 43.83ms |  22.626 |        0B | 24.888 |    20 |   22 |      884ms |
| Vectorized cubature::cuhre     |  43.9ms | 44.24ms |  22.466 |    2.04MB | 23.589 |    20 |   21 |    890.2ms |

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

| expression                     |   min |  median |   itr/sec | mem_alloc |  gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|------:|--------:|----------:|----------:|--------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 125µs | 128.2µs |  7553.943 |    53.7KB |   0.000 |   100 |    0 |    13.24ms |
| Vectorized Hcubature           | 212µs | 219.8µs |  4455.221 |    68.5KB |   0.000 |   100 |    0 |    22.45ms |
| Non-vectorized Pcubature       |  52µs |  55.3µs | 16929.692 |        0B | 171.007 |    99 |    1 |     5.85ms |
| Vectorized Pcubature           | 155µs |   161µs |  6022.636 |        0B |   0.000 |   100 |    0 |     16.6ms |
| Non-vectorized cubature::cuhre | 224µs | 232.8µs |  4218.364 |        0B |  42.610 |    99 |    1 |    23.47ms |
| Vectorized cubature::cuhre     | 502µs |   533µs |  1839.397 |        0B |  18.580 |    99 |    1 |    53.82ms |

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
| Non-vectorized Hcubature       |   4.54ms |   4.62ms |  215.793 |    78.4KB | 35.129 |    86 |   14 |    398.5ms |
| Vectorized Hcubature           |   5.48ms |   5.57ms |  179.211 |   304.7KB | 26.779 |    87 |   13 |    485.5ms |
| Non-vectorized Pcubature       | 407.01µs | 424.84µs | 2339.966 |        0B | 47.754 |    98 |    2 |     41.9ms |
| Vectorized Pcubature           | 601.66µs | 629.56µs | 1577.380 |    18.3KB | 15.933 |    99 |    1 |     62.8ms |
| Non-vectorized cubature::cuhre |   1.22ms |   1.25ms |  791.760 |        0B | 24.487 |    97 |    3 |    122.5ms |
| Vectorized cubature::cuhre     |   1.34ms |   1.39ms |  717.139 |    60.1KB | 22.180 |    97 |    3 |    135.3ms |

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
