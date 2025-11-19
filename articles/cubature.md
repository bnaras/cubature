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

| expression                     |      min |  median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|--------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 270.88µs | 282.6µs | 3512.711 |    59.3KB | 35.482 |    99 |    1 |     28.2ms |
| Vectorized Hcubature           | 389.52µs | 412.1µs | 2395.089 |    54.5KB | 24.193 |    99 |    1 |     41.3ms |
| Non-vectorized Pcubature       | 853.12µs | 886.1µs | 1114.087 |    31.2KB | 34.456 |    97 |    3 |     87.1ms |
| Vectorized Pcubature           |   1.16ms |   1.2ms |  826.148 |    58.4KB | 16.860 |    98 |    2 |    118.6ms |
| Non-vectorized cubature::cuhre | 586.12µs | 602.8µs | 1643.600 |    40.5KB | 16.602 |    99 |    1 |     60.2ms |
| Vectorized cubature::cuhre     | 588.08µs | 621.9µs | 1606.397 |    39.8KB | 32.784 |    98 |    2 |       61ms |

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
| Non-vectorized Hcubature       | 819.85ms | 827.25ms |   1.186 |  139.68KB | 18.389 |    10 |  155 |      8.43s |
| Vectorized Hcubature           |   1.77ms |   2.49ms | 440.758 |    1.89MB |  0.000 |    10 |    0 |    22.69ms |
| Non-vectorized Pcubature       | 343.79ms | 348.71ms |   2.834 |        0B | 18.421 |    10 |   65 |      3.53s |
| Vectorized Pcubature           |   1.25ms |   1.27ms | 682.377 |  810.26KB | 68.238 |    10 |    1 |    14.65ms |
| Non-vectorized cubature::cuhre | 330.49ms | 332.85ms |   3.008 |        0B | 18.948 |    10 |   63 |      3.33s |
| Vectorized cubature::cuhre     |   3.19ms |   3.25ms | 306.119 |  898.41KB |  0.000 |    10 |    0 |    32.67ms |

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
| g1()       | 760µs | 1.38ms | 735.101 |  355.13KB |      0 |    20 |    0 |     27.2ms |
| g2()       | 757µs | 1.36ms | 760.078 |    2.49KB |      0 |    20 |    0 |     26.3ms |
| g3()       | 747µs | 1.37ms | 738.318 |    2.49KB |      0 |    20 |    0 |     27.1ms |

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
| Non-vectorized Hcubature       |  37.7µs |  39.8µs | 24409.837 |    7.66KB | 24.434 |   999 |    1 |     40.9ms |
| Vectorized Hcubature           |  64.4µs |  67.6µs | 14315.604 |    1.16KB | 28.689 |   998 |    2 |     69.7ms |
| Non-vectorized Pcubature       |    49µs |  50.6µs | 19208.439 |        0B | 19.228 |   999 |    1 |       52ms |
| Vectorized Pcubature           | 113.5µs | 118.4µs |  8212.793 |   18.68KB | 24.713 |   997 |    3 |    121.4ms |
| Non-vectorized cubature::cuhre | 347.7µs | 359.9µs |  2760.175 |        0B | 25.067 |   991 |    9 |      359ms |
| Vectorized cubature::cuhre     | 360.8µs | 377.6µs |  2636.488 |   16.38KB | 26.631 |   990 |   10 |    375.5ms |

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
| Non-vectorized Hcubature       |   3.03ms |   3.06ms |   326.147 |    67.3KB | 36.239 |     9 |    1 |     27.6ms |
| Vectorized Hcubature           |   4.94ms |   5.02ms |   197.433 |   290.5KB | 21.937 |     9 |    1 |     45.6ms |
| Non-vectorized Pcubature       |  75.18µs |  76.66µs | 12524.679 |        0B |  0.000 |    10 |    0 |    798.4µs |
| Vectorized Pcubature           | 151.19µs | 154.62µs |  6261.958 |     4.1KB |  0.000 |    10 |    0 |      1.6ms |
| Non-vectorized cubature::cuhre |  13.38ms |  13.46ms |    74.142 |        0B | 31.775 |     7 |    3 |     94.4ms |
| Vectorized cubature::cuhre     |  20.03ms |  20.05ms |    49.443 |   971.5KB | 49.443 |     5 |    5 |    101.1ms |

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
| Non-vectorized Hcubature       |  46.8ms |  47.3ms |  18.439 |    17.7KB | 36.878 |    10 |   20 |   542.34ms |
| Vectorized Hcubature           |    47ms |  47.7ms |  20.900 |  1011.8KB | 27.170 |    10 |   13 |   478.47ms |
| Non-vectorized cubature::cuhre | 785.9ms |   821ms |   1.228 |        0B | 28.866 |    10 |  235 |      8.14s |
| Vectorized cubature::cuhre     | 884.9ms | 899.5ms |   1.103 |    21.2MB | 25.048 |    10 |  227 |      9.06s |

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
| Non-vectorized Hcubature       |    60µs |  61.7µs | 15783.593 |    6.75KB |   0.000 |    50 |    0 |     3.17ms |
| Vectorized Hcubature           |  91.1µs |    97µs |  9942.231 |    2.91KB |   0.000 |    50 |    0 |     5.03ms |
| Non-vectorized Pcubature       |  51.3µs |  55.3µs | 17228.781 |        0B | 351.608 |    49 |    1 |     2.84ms |
| Vectorized Pcubature           |  81.2µs |  83.7µs | 11607.434 |   19.12KB |   0.000 |    50 |    0 |     4.31ms |
| Non-vectorized cubature::cuhre | 640.6µs | 664.3µs |  1503.319 |        0B |  30.680 |    49 |    1 |    32.59ms |
| Vectorized cubature::cuhre     | 654.6µs | 679.3µs |  1465.846 |   39.84KB |   0.000 |    50 |    0 |    34.11ms |

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
| Non-vectorized Hcubature       | 1.41ms | 1.45ms | 690.288 |    85.2KB | 36.331 |    19 |    1 |     27.5ms |
| Vectorized Hcubature           | 1.76ms | 1.81ms | 550.014 |   147.2KB | 28.948 |    19 |    1 |     34.5ms |
| Non-vectorized Pcubature       | 2.06ms | 2.07ms | 479.560 |        0B | 25.240 |    19 |    1 |     39.6ms |
| Vectorized Pcubature           | 2.58ms | 2.64ms | 377.627 |    68.9KB | 19.875 |    19 |    1 |     50.3ms |
| Non-vectorized cubature::cuhre | 3.27ms | 3.35ms | 298.744 |        0B | 33.194 |    18 |    2 |     60.3ms |
| Vectorized cubature::cuhre     | 3.75ms | 3.79ms | 263.372 |   125.6KB | 13.862 |    19 |    1 |     72.1ms |

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
| Non-vectorized Hcubature       | 3.76ms | 3.83ms | 260.061 |     133KB | 28.896 |    18 |    2 |     69.2ms |
| Vectorized Hcubature           | 4.66ms | 4.76ms | 209.691 |     249KB | 37.004 |    17 |    3 |     81.1ms |
| Non-vectorized Pcubature       | 2.53ms | 2.57ms | 388.714 |        0B | 20.459 |    19 |    1 |     48.9ms |
| Vectorized Pcubature           | 3.28ms | 3.34ms | 298.501 |      69KB | 33.167 |    18 |    2 |     60.3ms |
| Non-vectorized cubature::cuhre | 7.02ms |  7.1ms | 140.698 |        0B | 35.174 |    16 |    4 |    113.7ms |
| Vectorized cubature::cuhre     | 8.27ms | 8.32ms | 118.237 |     224KB | 29.559 |    16 |    4 |    135.3ms |

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
| Non-vectorized Hcubature       |  1.64ms |  1.66ms | 597.769 |    64.6KB | 31.462 |    19 |    1 |     31.8ms |
| Vectorized Hcubature           |  2.14ms |  2.17ms | 458.678 |     156KB | 24.141 |    19 |    1 |     41.4ms |
| Non-vectorized Pcubature       |  8.29ms |  8.45ms | 118.735 |        0B | 39.578 |    15 |    5 |    126.3ms |
| Vectorized Pcubature           | 10.44ms | 10.69ms |  93.258 |   386.2KB | 31.086 |    15 |    5 |    160.8ms |
| Non-vectorized cubature::cuhre |  4.45ms |  4.51ms | 221.317 |        0B | 24.591 |    18 |    2 |     81.3ms |
| Vectorized cubature::cuhre     |   4.8ms |  4.88ms | 204.856 |   225.8KB | 22.762 |    18 |    2 |     87.9ms |

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
| Non-vectorized Hcubature       |  3.58ms |  3.66ms | 259.614 |   32.89KB | 25.961 |    20 |    2 |       77ms |
| Vectorized Hcubature           |  4.05ms |  4.15ms | 230.495 |  205.61KB | 23.050 |    20 |    2 |     86.8ms |
| Non-vectorized Pcubature       |   8.4ms |  8.56ms | 110.484 |        0B | 33.145 |    20 |    6 |      181ms |
| Vectorized Pcubature           |  9.57ms |  9.76ms |  97.284 |  386.24KB | 24.321 |    20 |    5 |    205.6ms |
| Non-vectorized cubature::cuhre |  43.9ms | 45.56ms |  21.902 |        0B | 24.092 |    20 |   22 |    913.2ms |
| Vectorized cubature::cuhre     | 44.52ms | 45.14ms |  22.086 |    2.04MB | 23.190 |    20 |   21 |    905.6ms |

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
| Non-vectorized Hcubature       | 126.2µs | 130.5µs |  7464.944 |    53.7KB |  0.000 |   100 |    0 |     13.4ms |
| Vectorized Hcubature           | 212.1µs | 223.4µs |  4302.687 |    68.5KB | 43.461 |    99 |    1 |    23.01ms |
| Non-vectorized Pcubature       |  51.4µs |  53.3µs | 18206.392 |        0B |  0.000 |   100 |    0 |     5.49ms |
| Vectorized Pcubature           | 154.5µs |   162µs |  6039.379 |        0B |  0.000 |   100 |    0 |    16.56ms |
| Non-vectorized cubature::cuhre | 223.1µs | 229.8µs |  4273.118 |        0B | 43.163 |    99 |    1 |    23.17ms |
| Vectorized cubature::cuhre     | 497.9µs | 523.2µs |  1893.390 |        0B | 19.125 |    99 |    1 |    52.29ms |

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
| Non-vectorized Hcubature       |   4.75ms |   4.82ms |  206.691 |    78.4KB | 30.885 |    87 |   13 |    420.9ms |
| Vectorized Hcubature           |   5.42ms |   5.53ms |  180.407 |   304.7KB | 29.369 |    86 |   14 |    476.7ms |
| Non-vectorized Pcubature       | 411.64µs | 426.05µs | 2334.373 |        0B | 23.580 |    99 |    1 |     42.4ms |
| Vectorized Pcubature           | 599.31µs | 625.27µs | 1590.094 |    18.3KB | 16.062 |    99 |    1 |     62.3ms |
| Non-vectorized cubature::cuhre |   1.23ms |   1.26ms |  789.506 |        0B | 24.418 |    97 |    3 |    122.9ms |
| Vectorized cubature::cuhre     |   1.33ms |   1.36ms |  728.257 |    60.1KB | 30.344 |    96 |    4 |    131.8ms |

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
