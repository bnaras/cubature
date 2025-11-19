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
| Non-vectorized Hcubature       | 270.83µs | 283.31µs | 3501.728 |    59.3KB | 35.371 |    99 |    1 |     28.3ms |
| Vectorized Hcubature           | 383.76µs | 406.01µs | 2425.728 |    54.5KB | 24.502 |    99 |    1 |     40.8ms |
| Non-vectorized Pcubature       | 851.15µs | 874.99µs | 1123.736 |    31.2KB | 34.755 |    97 |    3 |     86.3ms |
| Vectorized Pcubature           |   1.16ms |   1.19ms |  833.164 |    58.4KB | 17.003 |    98 |    2 |    117.6ms |
| Non-vectorized cubature::cuhre | 584.89µs | 603.94µs | 1642.542 |    40.5KB | 16.591 |    99 |    1 |     60.3ms |
| Vectorized cubature::cuhre     | 584.58µs | 622.43µs | 1600.641 |    39.8KB | 32.666 |    98 |    2 |     61.2ms |

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
| Non-vectorized Hcubature       | 818.18ms | 831.06ms |   1.186 |  139.68KB | 18.377 |    10 |  155 |      8.43s |
| Vectorized Hcubature           |    1.8ms |   2.52ms | 437.474 |    1.89MB |  0.000 |    10 |    0 |    22.86ms |
| Non-vectorized Pcubature       | 353.97ms | 356.67ms |   2.756 |        0B | 17.914 |    10 |   65 |      3.63s |
| Vectorized Pcubature           |   1.25ms |   1.29ms | 652.714 |  810.26KB | 65.271 |    10 |    1 |    15.32ms |
| Non-vectorized cubature::cuhre | 341.58ms | 347.64ms |   2.882 |        0B | 18.157 |    10 |   63 |      3.47s |
| Vectorized cubature::cuhre     |   3.25ms |   3.33ms | 298.864 |  898.41KB |  0.000 |    10 |    0 |    33.46ms |

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
| g1()       | 770µs | 1.38ms | 729.132 |  355.13KB |      0 |    20 |    0 |     27.4ms |
| g2()       | 750µs | 1.37ms | 754.272 |    2.49KB |      0 |    20 |    0 |     26.5ms |
| g3()       | 760µs | 1.37ms | 736.666 |    2.49KB |      0 |    20 |    0 |     27.1ms |

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
| Non-vectorized Hcubature       |  36.6µs |    39µs | 24958.746 |    7.66KB | 24.984 |   999 |    1 |       40ms |
| Vectorized Hcubature           |  64.8µs |    68µs | 14252.966 |    1.16KB | 28.563 |   998 |    2 |       70ms |
| Non-vectorized Pcubature       |  48.2µs |  50.6µs | 19120.831 |        0B | 19.140 |   999 |    1 |     52.2ms |
| Vectorized Pcubature           | 112.1µs | 117.9µs |  8302.082 |   18.68KB | 24.981 |   997 |    3 |    120.1ms |
| Non-vectorized cubature::cuhre | 339.6µs | 352.3µs |  2821.473 |        0B | 25.624 |   991 |    9 |    351.2ms |
| Vectorized cubature::cuhre     | 359.4µs | 376.3µs |  2640.433 |   16.38KB | 26.671 |   990 |   10 |    374.9ms |

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
| Non-vectorized Hcubature       |   2.85ms |   2.91ms |   344.234 |    67.3KB | 38.248 |     9 |    1 |    26.14ms |
| Vectorized Hcubature           |   4.89ms |   4.97ms |   201.400 |   290.5KB | 22.378 |     9 |    1 |    44.69ms |
| Non-vectorized Pcubature       |  73.21µs |  74.58µs | 12820.170 |        0B |  0.000 |    10 |    0 |   780.02µs |
| Vectorized Pcubature           | 151.93µs | 159.45µs |  6115.988 |     4.1KB |  0.000 |    10 |    0 |     1.64ms |
| Non-vectorized cubature::cuhre |   13.7ms |  13.76ms |    72.645 |        0B | 31.134 |     7 |    3 |    96.36ms |
| Vectorized cubature::cuhre     |  19.83ms |  20.28ms |    49.527 |   971.5KB | 49.527 |     5 |    5 |   100.95ms |

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
| Non-vectorized Hcubature       |  45.2ms |  45.6ms |  18.564 |    17.7KB | 37.128 |    10 |   20 |   538.68ms |
| Vectorized Hcubature           |  48.2ms |  48.7ms |  20.360 |  1011.8KB | 24.432 |    10 |   12 |   491.16ms |
| Non-vectorized cubature::cuhre | 798.8ms | 806.9ms |   1.230 |        0B | 29.024 |    10 |  236 |      8.13s |
| Vectorized cubature::cuhre     | 893.1ms | 908.6ms |   1.090 |    21.2MB | 24.734 |    10 |  227 |      9.18s |

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
| Non-vectorized Hcubature       |  57.6µs |  60.2µs | 15875.925 |    6.75KB |   0.000 |    50 |    0 |     3.15ms |
| Vectorized Hcubature           |  90.9µs |  95.2µs | 10140.928 |    2.91KB |   0.000 |    50 |    0 |     4.93ms |
| Non-vectorized Pcubature       |  51.4µs |  54.2µs | 17911.336 |        0B |   0.000 |    50 |    0 |     2.79ms |
| Vectorized Pcubature           |  81.6µs |  87.3µs | 10692.427 |   19.12KB | 218.213 |    49 |    1 |     4.58ms |
| Non-vectorized cubature::cuhre | 645.4µs | 665.9µs |  1492.657 |        0B |  30.462 |    49 |    1 |    32.83ms |
| Vectorized cubature::cuhre     | 646.5µs | 677.1µs |  1471.558 |   39.84KB |   0.000 |    50 |    0 |    33.98ms |

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
| Non-vectorized Hcubature       | 1.35ms | 1.38ms | 720.731 |    85.2KB | 37.933 |    19 |    1 |     26.4ms |
| Vectorized Hcubature           | 1.79ms | 1.82ms | 546.245 |   147.2KB | 28.750 |    19 |    1 |     34.8ms |
| Non-vectorized Pcubature       | 2.04ms | 2.08ms | 476.361 |        0B | 25.072 |    19 |    1 |     39.9ms |
| Vectorized Pcubature           | 2.57ms | 2.66ms | 375.420 |    68.9KB | 41.713 |    18 |    2 |     47.9ms |
| Non-vectorized cubature::cuhre | 3.29ms | 3.37ms | 296.932 |        0B | 15.628 |    19 |    1 |       64ms |
| Vectorized cubature::cuhre     | 3.71ms | 3.77ms | 265.063 |   125.6KB | 29.451 |    18 |    2 |     67.9ms |

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
| Non-vectorized Hcubature       | 3.63ms | 3.67ms | 270.953 |     133KB | 47.815 |    17 |    3 |     62.7ms |
| Vectorized Hcubature           | 4.65ms | 4.76ms | 210.297 |     249KB | 23.366 |    18 |    2 |     85.6ms |
| Non-vectorized Pcubature       | 2.48ms | 2.54ms | 393.179 |        0B | 43.687 |    18 |    2 |     45.8ms |
| Vectorized Pcubature           | 3.26ms | 3.34ms | 299.380 |      69KB | 15.757 |    19 |    1 |     63.5ms |
| Non-vectorized cubature::cuhre | 7.08ms | 7.17ms | 139.613 |        0B | 46.538 |    15 |    5 |    107.4ms |
| Vectorized cubature::cuhre     |  8.3ms | 8.39ms | 119.107 |     224KB | 21.019 |    17 |    3 |    142.7ms |

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
| Non-vectorized Hcubature       |  1.56ms |  1.59ms | 628.962 |    64.6KB | 33.103 |    19 |    1 |     30.2ms |
| Vectorized Hcubature           |  2.15ms |   2.2ms | 454.154 |     156KB |  0.000 |    20 |    0 |       44ms |
| Non-vectorized Pcubature       |  8.11ms |  8.28ms | 121.058 |        0B | 40.353 |    15 |    5 |    123.9ms |
| Vectorized Pcubature           | 10.27ms | 10.47ms |  95.705 |   386.2KB | 31.902 |    15 |    5 |    156.7ms |
| Non-vectorized cubature::cuhre |  4.54ms |   4.6ms | 217.376 |        0B | 24.153 |    18 |    2 |     82.8ms |
| Vectorized cubature::cuhre     |  4.67ms |  4.77ms | 209.784 |   225.8KB | 23.309 |    18 |    2 |     85.8ms |

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
| Non-vectorized Hcubature       |  3.35ms |  3.46ms | 271.442 |   32.89KB | 27.144 |    20 |    2 |     73.7ms |
| Vectorized Hcubature           |  4.04ms |  4.19ms | 227.392 |  205.61KB | 22.739 |    20 |    2 |       88ms |
| Non-vectorized Pcubature       |  8.27ms |   8.4ms | 111.545 |        0B | 27.886 |    20 |    5 |    179.3ms |
| Vectorized Pcubature           |  9.54ms |  9.73ms |  97.667 |  386.24KB | 24.417 |    20 |    5 |    204.8ms |
| Non-vectorized cubature::cuhre | 46.33ms |  46.8ms |  21.181 |        0B | 23.299 |    20 |   22 |    944.2ms |
| Vectorized cubature::cuhre     | 43.54ms | 44.02ms |  22.581 |    2.04MB | 23.710 |    20 |   21 |    885.7ms |

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
| Non-vectorized Hcubature       | 126.3µs | 137.7µs |  7056.159 |    53.7KB | 71.274 |    99 |    1 |       14ms |
| Vectorized Hcubature           | 209.1µs | 216.9µs |  4481.481 |    68.5KB |  0.000 |   100 |    0 |     22.3ms |
| Non-vectorized Pcubature       |  52.1µs |  53.8µs | 17852.223 |        0B |  0.000 |   100 |    0 |      5.6ms |
| Vectorized Pcubature           | 156.2µs | 166.2µs |  5808.610 |        0B | 58.673 |    99 |    1 |       17ms |
| Non-vectorized cubature::cuhre | 253.3µs | 260.8µs |  3791.155 |        0B |  0.000 |   100 |    0 |     26.4ms |
| Vectorized cubature::cuhre     | 498.9µs | 525.5µs |  1884.318 |        0B | 19.034 |    99 |    1 |     52.5ms |

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
| Non-vectorized Hcubature       |   4.54ms |   4.63ms |  213.886 |    78.4KB | 34.819 |    86 |   14 |    402.1ms |
| Vectorized Hcubature           |   5.42ms |   5.58ms |  178.468 |   304.7KB | 26.668 |    87 |   13 |    487.5ms |
| Non-vectorized Pcubature       | 405.72µs |  425.3µs | 2336.200 |        0B | 23.598 |    99 |    1 |     42.4ms |
| Vectorized Pcubature           | 602.37µs | 632.55µs | 1561.982 |    18.3KB | 31.877 |    98 |    2 |     62.7ms |
| Non-vectorized cubature::cuhre |   1.23ms |   1.27ms |  785.515 |        0B | 24.294 |    97 |    3 |    123.5ms |
| Vectorized cubature::cuhre     |   1.33ms |   1.38ms |  719.601 |    60.1KB | 22.256 |    97 |    3 |    134.8ms |

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
