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
| Non-vectorized Hcubature       | 276.89µs | 290.28µs | 3415.139 |    59.3KB | 34.496 |    99 |    1 |       29ms |
| Vectorized Hcubature           | 398.28µs | 431.64µs | 2168.452 |    54.5KB | 21.904 |    99 |    1 |     45.7ms |
| Non-vectorized Pcubature       | 891.36µs | 925.32µs | 1033.260 |    31.2KB | 31.956 |    97 |    3 |     93.9ms |
| Vectorized Pcubature           |   1.18ms |   1.36ms |  670.072 |    58.4KB | 13.675 |    98 |    2 |    146.3ms |
| Non-vectorized cubature::cuhre | 588.62µs | 607.82µs | 1635.423 |    40.5KB | 16.519 |    99 |    1 |     60.5ms |
| Vectorized cubature::cuhre     | 589.36µs | 623.38µs | 1601.297 |    39.8KB | 32.680 |    98 |    2 |     61.2ms |

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
| Non-vectorized Hcubature       | 829.86ms | 874.92ms |   1.133 |  139.68KB | 17.566 |    10 |  155 |      8.82s |
| Vectorized Hcubature           |   1.85ms |   2.54ms | 425.762 |    1.89MB |  0.000 |    10 |    0 |    23.49ms |
| Non-vectorized Pcubature       | 360.51ms | 376.35ms |   2.580 |        0B | 17.029 |    10 |   66 |      3.88s |
| Vectorized Pcubature           |    1.3ms |   1.36ms | 720.689 |  810.26KB |  0.000 |    10 |    0 |    13.88ms |
| Non-vectorized cubature::cuhre | 347.28ms |  351.8ms |   2.813 |        0B | 17.720 |    10 |   63 |      3.56s |
| Vectorized cubature::cuhre     |   3.22ms |   3.27ms | 304.121 |  898.41KB |  0.000 |    10 |    0 |    32.88ms |

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
| g1()       | 774µs | 1.38ms | 731.863 |  355.13KB |      0 |    20 |    0 |     27.3ms |
| g2()       | 763µs | 1.37ms | 757.563 |    2.49KB |      0 |    20 |    0 |     26.4ms |
| g3()       | 758µs | 1.37ms | 740.041 |    2.49KB |      0 |    20 |    0 |       27ms |

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
| Non-vectorized Hcubature       |    37µs |  39.3µs | 24585.260 |    7.66KB | 24.610 |   999 |    1 |     40.6ms |
| Vectorized Hcubature           |  64.4µs |  68.1µs | 14247.830 |    1.16KB | 28.553 |   998 |    2 |       70ms |
| Non-vectorized Pcubature       |  48.8µs |  50.7µs | 19197.970 |        0B | 19.217 |   999 |    1 |       52ms |
| Vectorized Pcubature           | 113.3µs | 119.5µs |  8160.193 |   18.68KB | 24.554 |   997 |    3 |    122.2ms |
| Non-vectorized cubature::cuhre |   338µs | 350.7µs |  2828.758 |        0B | 25.690 |   991 |    9 |    350.3ms |
| Vectorized cubature::cuhre     | 365.9µs | 384.8µs |  2574.473 |   16.38KB | 26.005 |   990 |   10 |    384.5ms |

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

| expression                     |      min |   median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|---------:|---------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       |   2.94ms |   3.15ms |  309.208 |    67.3KB | 34.356 |     9 |    1 |    29.11ms |
| Vectorized Hcubature           |   5.29ms |   5.54ms |  174.544 |   290.5KB | 19.394 |     9 |    1 |    51.56ms |
| Non-vectorized Pcubature       |  92.04µs | 101.04µs | 9184.563 |        0B |  0.000 |    10 |    0 |     1.09ms |
| Vectorized Pcubature           | 201.21µs | 219.94µs | 4265.123 |     4.1KB |  0.000 |    10 |    0 |     2.35ms |
| Non-vectorized cubature::cuhre |  13.67ms |  14.95ms |   62.895 |        0B | 26.955 |     7 |    3 |    111.3ms |
| Vectorized cubature::cuhre     |  20.48ms |  20.59ms |   44.202 |   971.5KB | 44.202 |     5 |    5 |   113.12ms |

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
| Non-vectorized Hcubature       |  47.9ms |  56.8ms |  15.473 |    17.7KB | 30.947 |    10 |   20 |   646.27ms |
| Vectorized Hcubature           |    49ms |  54.7ms |  17.799 |  1011.8KB | 23.139 |    10 |   13 |   561.82ms |
| Non-vectorized cubature::cuhre | 810.7ms | 837.6ms |   1.180 |        0B | 27.729 |    10 |  235 |      8.47s |
| Vectorized cubature::cuhre     | 921.9ms | 980.6ms |   1.013 |    21.2MB | 23.000 |    10 |  227 |      9.87s |

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
| Non-vectorized Hcubature       |  59.2µs |  62.3µs | 15520.410 |    6.75KB |   0.000 |    50 |    0 |     3.22ms |
| Vectorized Hcubature           |  94.3µs |  99.5µs |  9723.500 |    2.91KB | 198.439 |    49 |    1 |     5.04ms |
| Non-vectorized Pcubature       |    50µs |  52.1µs | 17579.928 |        0B |   0.000 |    50 |    0 |     2.84ms |
| Vectorized Pcubature           |  82.4µs |  85.5µs | 11325.295 |   19.12KB |   0.000 |    50 |    0 |     4.42ms |
| Non-vectorized cubature::cuhre | 656.3µs | 784.9µs |  1299.259 |        0B |  26.515 |    49 |    1 |    37.71ms |
| Vectorized cubature::cuhre     | 677.4µs | 889.9µs |  1145.959 |   39.84KB |  23.387 |    49 |    1 |    42.76ms |

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
| Non-vectorized Hcubature       | 1.36ms | 1.42ms | 702.491 |    85.2KB | 36.973 |    19 |    1 |       27ms |
| Vectorized Hcubature           | 1.82ms | 1.87ms | 534.314 |   147.2KB |  0.000 |    20 |    0 |     37.4ms |
| Non-vectorized Pcubature       | 2.02ms | 2.08ms | 478.056 |        0B | 53.117 |    18 |    2 |     37.7ms |
| Vectorized Pcubature           | 2.65ms | 2.71ms | 369.083 |    68.9KB | 19.425 |    19 |    1 |     51.5ms |
| Non-vectorized cubature::cuhre | 3.27ms | 3.35ms | 298.468 |        0B | 33.163 |    18 |    2 |     60.3ms |
| Vectorized cubature::cuhre     | 3.84ms |  3.9ms | 255.221 |   125.6KB | 28.358 |    18 |    2 |     70.5ms |

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
| Non-vectorized Hcubature       | 3.62ms | 3.71ms | 269.164 |     133KB | 29.907 |    18 |    2 |     66.9ms |
| Vectorized Hcubature           | 4.67ms | 4.79ms | 208.373 |     249KB | 23.153 |    18 |    2 |     86.4ms |
| Non-vectorized Pcubature       | 2.52ms | 2.56ms | 390.300 |        0B | 43.367 |    18 |    2 |     46.1ms |
| Vectorized Pcubature           |  3.4ms | 3.46ms | 287.341 |      69KB | 31.927 |    18 |    2 |     62.6ms |
| Non-vectorized cubature::cuhre |    7ms | 7.22ms | 137.702 |        0B | 34.426 |    16 |    4 |    116.2ms |
| Vectorized cubature::cuhre     | 8.33ms | 8.48ms | 116.317 |     224KB | 29.079 |    16 |    4 |    137.6ms |

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
| Non-vectorized Hcubature       |   1.6ms |  1.65ms | 603.422 |    64.6KB | 31.759 |    19 |    1 |     31.5ms |
| Vectorized Hcubature           |  2.11ms |  2.19ms | 455.684 |     156KB | 23.983 |    19 |    1 |     41.7ms |
| Non-vectorized Pcubature       |  8.19ms |  8.43ms | 116.220 |        0B | 38.740 |    15 |    5 |    129.1ms |
| Vectorized Pcubature           | 10.49ms | 11.33ms |  86.018 |   386.2KB | 21.504 |    16 |    4 |      186ms |
| Non-vectorized cubature::cuhre |  4.33ms |  4.97ms | 203.336 |        0B | 35.883 |    17 |    3 |     83.6ms |
| Vectorized cubature::cuhre     |  4.82ms |  4.95ms | 201.839 |   225.8KB | 22.427 |    18 |    2 |     89.2ms |

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
| Non-vectorized Hcubature       |   3.4ms |  3.48ms | 267.658 |   32.89KB | 26.766 |    20 |    2 |     74.7ms |
| Vectorized Hcubature           |  4.03ms |  4.21ms | 221.814 |  205.61KB | 22.181 |    20 |    2 |     90.2ms |
| Non-vectorized Pcubature       |  8.48ms |  8.73ms | 107.298 |        0B | 32.189 |    20 |    6 |    186.4ms |
| Vectorized Pcubature           |  9.84ms | 10.21ms |  89.874 |  386.24KB | 22.468 |    20 |    5 |    222.5ms |
| Non-vectorized cubature::cuhre | 43.36ms | 45.79ms |  21.057 |        0B | 22.110 |    20 |   21 |    949.8ms |
| Vectorized cubature::cuhre     | 45.36ms |  46.9ms |  20.932 |    2.04MB | 23.025 |    20 |   22 |    955.5ms |

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
| Non-vectorized Hcubature       | 128.9µs | 139.3µs |  6965.644 |    53.7KB |   0.000 |   100 |    0 |    14.36ms |
| Vectorized Hcubature           | 213.4µs | 223.9µs |  4375.313 |    68.5KB |   0.000 |   100 |    0 |    22.86ms |
| Non-vectorized Pcubature       |  52.5µs |  59.2µs | 15803.141 |        0B | 159.628 |    99 |    1 |     6.26ms |
| Vectorized Pcubature           | 180.2µs | 202.9µs |  4828.139 |        0B |   0.000 |   100 |    0 |    20.71ms |
| Non-vectorized cubature::cuhre | 236.8µs | 260.9µs |  3644.289 |        0B |  36.811 |    99 |    1 |    27.17ms |
| Vectorized cubature::cuhre     | 518.8µs | 647.7µs |  1503.560 |        0B |  15.187 |    99 |    1 |    65.84ms |

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
| Non-vectorized Hcubature       |   4.64ms |   4.83ms |  203.735 |    78.4KB | 33.166 |    86 |   14 |    422.1ms |
| Vectorized Hcubature           |   5.45ms |   5.62ms |  176.360 |   304.7KB | 26.353 |    87 |   13 |    493.3ms |
| Non-vectorized Pcubature       | 415.93µs | 433.35µs | 2271.984 |        0B | 46.367 |    98 |    2 |     43.1ms |
| Vectorized Pcubature           | 605.47µs |  644.5µs | 1505.084 |    18.3KB | 15.203 |    99 |    1 |     65.8ms |
| Non-vectorized cubature::cuhre |   1.26ms |   1.29ms |  767.877 |        0B | 23.749 |    97 |    3 |    126.3ms |
| Vectorized cubature::cuhre     |   1.34ms |    1.4ms |  709.948 |    60.1KB | 21.957 |    97 |    3 |    136.6ms |

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
