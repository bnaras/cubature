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
| Non-vectorized Hcubature       | 273.03µs | 285.87µs | 3458.035 |    59.3KB | 34.930 |    99 |    1 |     28.6ms |
| Vectorized Hcubature           | 394.41µs | 416.51µs | 2373.099 |    54.5KB | 23.971 |    99 |    1 |     41.7ms |
| Non-vectorized Pcubature       | 835.62µs | 875.92µs | 1126.957 |    31.2KB | 34.854 |    97 |    3 |     86.1ms |
| Vectorized Pcubature           |   1.17ms |   1.21ms |  819.833 |    58.4KB | 16.731 |    98 |    2 |    119.5ms |
| Non-vectorized cubature::cuhre |    583µs | 606.57µs | 1638.584 |    40.5KB | 16.551 |    99 |    1 |     60.4ms |
| Vectorized cubature::cuhre     | 593.49µs | 630.93µs | 1580.865 |    39.8KB | 32.263 |    98 |    2 |       62ms |

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
| Non-vectorized Hcubature       | 826.74ms | 831.51ms |   1.182 |  139.68KB | 18.317 |    10 |  155 |      8.46s |
| Vectorized Hcubature           |   1.96ms |   2.59ms | 409.841 |    1.89MB |  0.000 |    10 |    0 |     24.4ms |
| Non-vectorized Pcubature       | 352.55ms | 357.11ms |   2.759 |        0B | 17.930 |    10 |   65 |      3.62s |
| Vectorized Pcubature           |    1.3ms |   1.32ms | 644.505 |  810.26KB | 64.450 |    10 |    1 |    15.52ms |
| Non-vectorized cubature::cuhre |  339.7ms | 346.11ms |   2.885 |        0B | 18.173 |    10 |   63 |      3.47s |
| Vectorized cubature::cuhre     |   3.28ms |   3.34ms | 296.827 |  898.41KB |  0.000 |    10 |    0 |    33.69ms |

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
| g1()       | 759µs | 1.37ms | 734.081 |  355.13KB |      0 |    20 |    0 |     27.2ms |
| g2()       | 756µs | 1.37ms | 757.925 |    2.49KB |      0 |    20 |    0 |     26.4ms |
| g3()       | 748µs | 1.37ms | 740.777 |    2.49KB |      0 |    20 |    0 |       27ms |

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
| Non-vectorized Hcubature       |  36.8µs |  39.2µs | 24783.435 |    7.66KB | 24.808 |   999 |    1 |     40.3ms |
| Vectorized Hcubature           |  64.5µs |  67.8µs | 14291.686 |    1.16KB | 28.641 |   998 |    2 |     69.8ms |
| Non-vectorized Pcubature       |  49.1µs |  51.8µs | 18551.706 |        0B | 18.570 |   999 |    1 |     53.8ms |
| Vectorized Pcubature           | 114.8µs | 121.9µs |  7994.687 |   18.68KB | 24.056 |   997 |    3 |    124.7ms |
| Non-vectorized cubature::cuhre | 350.3µs | 365.6µs |  2718.281 |        0B | 24.687 |   991 |    9 |    364.6ms |
| Vectorized cubature::cuhre     | 363.8µs | 381.9µs |  2602.766 |   16.38KB | 26.291 |   990 |   10 |    380.4ms |

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
| Non-vectorized Hcubature       |   2.93ms |   2.99ms |   332.879 |    67.3KB |  0.000 |    10 |    0 |    30.04ms |
| Vectorized Hcubature           |      5ms |   5.12ms |   196.367 |   290.5KB | 49.092 |     8 |    2 |    40.74ms |
| Non-vectorized Pcubature       |  74.81µs |  76.77µs | 12327.287 |        0B |  0.000 |    10 |    0 |   811.21µs |
| Vectorized Pcubature           | 155.57µs | 168.14µs |  5951.872 |     4.1KB |  0.000 |    10 |    0 |     1.68ms |
| Non-vectorized cubature::cuhre |  13.64ms |  13.73ms |    72.580 |        0B | 31.106 |     7 |    3 |    96.44ms |
| Vectorized cubature::cuhre     |  20.32ms |  20.42ms |    48.912 |   971.5KB | 48.912 |     5 |    5 |   102.22ms |

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
| Non-vectorized Hcubature       |  44.8ms |  45.4ms |  18.402 |    17.7KB | 36.804 |    10 |   20 |   543.42ms |
| Vectorized Hcubature           |  47.9ms |  48.7ms |  20.394 |  1011.8KB | 26.512 |    10 |   13 |   490.35ms |
| Non-vectorized cubature::cuhre | 792.3ms | 807.9ms |   1.232 |        0B | 28.944 |    10 |  235 |      8.12s |
| Vectorized cubature::cuhre     | 907.4ms | 918.7ms |   1.069 |    21.2MB | 24.256 |    10 |  227 |      9.36s |

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
| Non-vectorized Hcubature       |  58.8µs |  60.6µs | 16011.875 |    6.75KB |   0.000 |    50 |    0 |     3.12ms |
| Vectorized Hcubature           |  94.3µs |   100µs |  9709.057 |    2.91KB |   0.000 |    50 |    0 |     5.15ms |
| Non-vectorized Pcubature       |  50.9µs |  54.9µs | 17285.306 |        0B | 352.761 |    49 |    1 |     2.83ms |
| Vectorized Pcubature           |    81µs |  84.6µs | 11299.722 |   19.12KB |   0.000 |    50 |    0 |     4.42ms |
| Non-vectorized cubature::cuhre | 644.5µs | 678.4µs |  1468.276 |        0B |  29.965 |    49 |    1 |    33.37ms |
| Vectorized cubature::cuhre     | 653.7µs | 685.6µs |  1454.797 |   39.84KB |   0.000 |    50 |    0 |    34.37ms |

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
| Non-vectorized Hcubature       | 1.35ms |  1.4ms | 714.250 |    85.2KB | 37.592 |    19 |    1 |     26.6ms |
| Vectorized Hcubature           | 1.79ms | 1.83ms | 545.258 |   147.2KB | 28.698 |    19 |    1 |     34.8ms |
| Non-vectorized Pcubature       | 2.03ms | 2.06ms | 480.389 |        0B | 53.377 |    18 |    2 |     37.5ms |
| Vectorized Pcubature           | 2.63ms |  2.7ms | 369.269 |    68.9KB | 19.435 |    19 |    1 |     51.5ms |
| Non-vectorized cubature::cuhre | 3.27ms | 3.35ms | 298.337 |        0B | 33.149 |    18 |    2 |     60.3ms |
| Vectorized cubature::cuhre     | 3.79ms | 3.87ms | 258.439 |   125.6KB | 13.602 |    19 |    1 |     73.5ms |

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
| Non-vectorized Hcubature       | 3.66ms | 3.71ms | 268.961 |     133KB | 29.885 |    18 |    2 |     66.9ms |
| Vectorized Hcubature           | 4.64ms | 4.74ms | 210.723 |     249KB | 37.186 |    17 |    3 |     80.7ms |
| Non-vectorized Pcubature       | 2.51ms | 2.57ms | 389.023 |        0B | 20.475 |    19 |    1 |     48.8ms |
| Vectorized Pcubature           | 3.26ms | 3.31ms | 300.823 |      69KB | 33.425 |    18 |    2 |     59.8ms |
| Non-vectorized cubature::cuhre |    7ms | 7.09ms | 141.076 |        0B | 35.269 |    16 |    4 |    113.4ms |
| Vectorized cubature::cuhre     | 8.26ms | 8.35ms | 119.641 |     224KB | 29.910 |    16 |    4 |    133.7ms |

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

| expression                     |     min | median | itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|-------:|--------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |  1.57ms | 1.61ms | 597.743 |    64.6KB | 31.460 |    19 |    1 |     31.8ms |
| Vectorized Hcubature           |  2.15ms | 2.17ms | 455.106 |     156KB | 23.953 |    19 |    1 |     41.7ms |
| Non-vectorized Pcubature       |  8.12ms | 8.28ms | 121.012 |        0B | 40.337 |    15 |    5 |      124ms |
| Vectorized Pcubature           | 10.48ms | 10.7ms |  93.390 |   386.2KB | 31.130 |    15 |    5 |    160.6ms |
| Non-vectorized cubature::cuhre |   4.4ms | 4.51ms | 222.274 |        0B | 24.697 |    18 |    2 |       81ms |
| Vectorized cubature::cuhre     |  4.73ms | 4.81ms | 207.640 |   225.8KB | 23.071 |    18 |    2 |     86.7ms |

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
| Non-vectorized Hcubature       |  3.37ms |  3.47ms | 273.478 |   32.89KB | 27.348 |    20 |    2 |     73.1ms |
| Vectorized Hcubature           |  4.03ms |  4.21ms | 228.125 |  205.61KB | 22.812 |    20 |    2 |     87.7ms |
| Non-vectorized Pcubature       |  8.23ms |  8.49ms | 112.600 |        0B | 28.150 |    20 |    5 |    177.6ms |
| Vectorized Pcubature           |  9.56ms |  9.75ms |  97.605 |  386.24KB | 24.401 |    20 |    5 |    204.9ms |
| Non-vectorized cubature::cuhre |  42.9ms | 43.57ms |  22.424 |        0B | 24.667 |    20 |   22 |    891.9ms |
| Vectorized cubature::cuhre     | 43.75ms | 44.48ms |  22.055 |    2.04MB | 23.158 |    20 |   21 |    906.8ms |

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
| Non-vectorized Hcubature       | 125.8µs | 134.9µs |  7164.550 |    53.7KB | 72.369 |    99 |    1 |    13.82ms |
| Vectorized Hcubature           | 211.6µs | 222.2µs |  4369.893 |    68.5KB |  0.000 |   100 |    0 |    22.88ms |
| Non-vectorized Pcubature       |  52.7µs |  54.7µs | 17461.127 |        0B |  0.000 |   100 |    0 |     5.73ms |
| Vectorized Pcubature           | 158.1µs | 169.4µs |  5577.978 |        0B | 56.343 |    99 |    1 |    17.75ms |
| Non-vectorized cubature::cuhre | 222.2µs | 232.9µs |  4219.745 |        0B |  0.000 |   100 |    0 |     23.7ms |
| Vectorized cubature::cuhre     | 502.6µs | 530.3µs |  1856.510 |        0B | 18.753 |    99 |    1 |    53.33ms |

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
| Non-vectorized Hcubature       |   4.48ms |   4.61ms |  216.544 |    78.4KB | 35.251 |    86 |   14 |    397.1ms |
| Vectorized Hcubature           |   5.44ms |   5.57ms |  178.670 |   304.7KB | 26.698 |    87 |   13 |    486.9ms |
| Non-vectorized Pcubature       | 408.69µs | 425.71µs | 2315.690 |        0B | 23.391 |    99 |    1 |     42.8ms |
| Vectorized Pcubature           | 606.83µs | 638.73µs | 1542.423 |    18.3KB | 31.478 |    98 |    2 |     63.5ms |
| Non-vectorized cubature::cuhre |   1.23ms |   1.26ms |  786.757 |        0B | 24.333 |    97 |    3 |    123.3ms |
| Vectorized cubature::cuhre     |   1.34ms |   1.38ms |  720.052 |    60.1KB | 22.270 |    97 |    3 |    134.7ms |

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
    ## [1] mvtnorm_1.3-3    cubature_2.1.4-4 bench_1.1.4     
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
