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
| Non-vectorized Hcubature       | 276.92µs | 292.33µs | 3353.993 |    59.3KB | 33.879 |    99 |    1 |     29.5ms |
| Vectorized Hcubature           | 397.63µs | 423.24µs | 2322.795 |    54.5KB | 23.463 |    99 |    1 |     42.6ms |
| Non-vectorized Pcubature       | 888.68µs | 915.84µs | 1077.131 |    31.2KB | 33.313 |    97 |    3 |     90.1ms |
| Vectorized Pcubature           |   1.22ms |   1.36ms |  709.154 |    58.4KB | 14.473 |    98 |    2 |    138.2ms |
| Non-vectorized cubature::cuhre | 603.03µs | 634.87µs | 1526.842 |    40.5KB | 15.423 |    99 |    1 |     64.8ms |
| Vectorized cubature::cuhre     | 609.71µs | 648.38µs | 1524.448 |    39.8KB | 31.111 |    98 |    2 |     64.3ms |

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
| Non-vectorized Hcubature       | 850.98ms | 860.41ms |   1.135 |  139.68KB | 17.600 |    10 |  155 |      8.81s |
| Vectorized Hcubature           |   1.84ms |   2.58ms | 414.241 |    1.89MB |  0.000 |    10 |    0 |    24.14ms |
| Non-vectorized Pcubature       |    359ms | 366.25ms |   2.685 |        0B | 17.723 |    10 |   66 |      3.72s |
| Vectorized Pcubature           |   1.27ms |   1.31ms | 750.199 |  810.26KB |  0.000 |    10 |    0 |    13.33ms |
| Non-vectorized cubature::cuhre | 350.82ms | 359.44ms |   2.781 |        0B | 17.520 |    10 |   63 |       3.6s |
| Vectorized cubature::cuhre     |    3.3ms |   3.35ms | 294.462 |  898.41KB |  0.000 |    10 |    0 |    33.96ms |

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
| g1()       | 762µs | 1.38ms | 726.222 |  355.13KB |      0 |    20 |    0 |     27.5ms |
| g2()       | 751µs | 1.38ms | 744.797 |    2.49KB |      0 |    20 |    0 |     26.9ms |
| g3()       | 758µs | 1.38ms | 729.606 |    2.49KB |      0 |    20 |    0 |     27.4ms |

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
| Non-vectorized Hcubature       |  37.2µs |  39.3µs | 24542.831 |    7.66KB | 24.567 |   999 |    1 |     40.7ms |
| Vectorized Hcubature           |  64.7µs |  68.6µs | 13920.306 |    1.16KB | 27.896 |   998 |    2 |     71.7ms |
| Non-vectorized Pcubature       |  49.1µs |  52.1µs | 18269.983 |        0B | 18.288 |   999 |    1 |     54.7ms |
| Vectorized Pcubature           | 115.5µs | 123.5µs |  7669.989 |   18.68KB | 23.079 |   997 |    3 |      130ms |
| Non-vectorized cubature::cuhre | 338.3µs | 353.8µs |  2809.813 |        0B | 25.518 |   991 |    9 |    352.7ms |
| Vectorized cubature::cuhre     | 373.1µs | 394.6µs |  2500.307 |   16.38KB | 25.256 |   990 |   10 |      396ms |

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
| Non-vectorized Hcubature       |   2.94ms |     3ms |   333.665 |    67.3KB | 37.074 |     9 |    1 |    26.97ms |
| Vectorized Hcubature           |   5.05ms |  5.18ms |   193.427 |   290.5KB | 21.492 |     9 |    1 |    46.53ms |
| Non-vectorized Pcubature       |  74.67µs | 78.73µs | 11565.752 |        0B |  0.000 |    10 |    0 |   864.62µs |
| Vectorized Pcubature           | 155.95µs | 160.2µs |  6031.400 |     4.1KB |  0.000 |    10 |    0 |     1.66ms |
| Non-vectorized cubature::cuhre |  13.76ms | 13.99ms |    71.690 |        0B | 47.793 |     6 |    4 |    83.69ms |
| Vectorized cubature::cuhre     |  20.79ms | 21.17ms |    47.450 |   971.5KB | 31.633 |     6 |    4 |   126.45ms |

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
| Non-vectorized Hcubature       |  46.7ms |  49.3ms |  17.532 |    17.7KB | 33.312 |    10 |   19 |   570.37ms |
| Vectorized Hcubature           |  48.9ms |  50.4ms |  19.763 |  1011.8KB | 25.692 |    10 |   13 |      506ms |
| Non-vectorized cubature::cuhre | 828.8ms | 839.9ms |   1.177 |        0B | 27.788 |    10 |  236 |      8.49s |
| Vectorized cubature::cuhre     | 936.7ms | 953.1ms |   1.039 |    21.2MB | 23.485 |    10 |  226 |      9.62s |

## A Simple Polynomial (product of coordinates)

``` r
testFn3 <- function(x) prod(2 * x)
testFn3_v <- function(x) matrix(apply(x, 2, function(z) prod(2 * z)), ncol = ncol(x))

d <- harness(f = testFn3, fv = testFn3_v,
             lowerLimit = rep(0, 3), upperLimit = rep(1, 3), iterations = 50)
knitr::kable(d, digits = 3)
```

| expression                     |    min |  median |   itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|-------:|--------:|----------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 57.7µs |  59.5µs | 15921.904 |    6.75KB |  0.000 |    50 |    0 |     3.14ms |
| Vectorized Hcubature           |   92µs |  95.3µs | 10025.112 |    2.91KB |  0.000 |    50 |    0 |     4.99ms |
| Non-vectorized Pcubature       | 50.6µs |  52.3µs | 18118.257 |        0B |  0.000 |    50 |    0 |     2.76ms |
| Vectorized Pcubature           | 82.3µs |  86.5µs | 11151.040 |   19.12KB |  0.000 |    50 |    0 |     4.48ms |
| Non-vectorized cubature::cuhre |  660µs | 680.3µs |  1450.728 |        0B | 29.607 |    49 |    1 |    33.78ms |
| Vectorized cubature::cuhre     |  661µs |   690µs |  1433.483 |   39.84KB | 29.255 |    49 |    1 |    34.18ms |

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
| Non-vectorized Hcubature       | 1.36ms | 1.41ms | 698.601 |    85.2KB | 36.768 |    19 |    1 |     27.2ms |
| Vectorized Hcubature           | 1.82ms | 1.87ms | 532.639 |   147.2KB | 28.034 |    19 |    1 |     35.7ms |
| Non-vectorized Pcubature       | 2.07ms |  2.1ms | 473.751 |        0B | 24.934 |    19 |    1 |     40.1ms |
| Vectorized Pcubature           | 2.63ms | 2.72ms | 367.396 |    68.9KB | 19.337 |    19 |    1 |     51.7ms |
| Non-vectorized cubature::cuhre | 3.26ms | 3.32ms | 299.060 |        0B | 33.229 |    18 |    2 |     60.2ms |
| Vectorized cubature::cuhre     | 3.81ms |  3.9ms | 255.631 |   125.6KB | 28.403 |    18 |    2 |     70.4ms |

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
| Non-vectorized Hcubature       | 3.61ms | 3.69ms | 269.869 |     133KB | 29.985 |    18 |    2 |     66.7ms |
| Vectorized Hcubature           |  4.8ms | 4.91ms | 203.435 |     249KB | 35.900 |    17 |    3 |     83.6ms |
| Non-vectorized Pcubature       | 2.56ms |  2.6ms | 384.240 |        0B | 20.223 |    19 |    1 |     49.4ms |
| Vectorized Pcubature           | 3.34ms | 3.46ms | 290.663 |      69KB | 32.296 |    18 |    2 |     61.9ms |
| Non-vectorized cubature::cuhre | 6.99ms | 7.17ms | 139.493 |        0B | 34.873 |    16 |    4 |    114.7ms |
| Vectorized cubature::cuhre     | 8.52ms | 8.64ms | 115.238 |     224KB | 28.809 |    16 |    4 |    138.8ms |

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
| Non-vectorized Hcubature       | 1.61ms |  1.65ms | 601.684 |    64.6KB | 31.668 |    19 |    1 |     31.6ms |
| Vectorized Hcubature           | 2.18ms |  2.22ms | 446.708 |     156KB | 23.511 |    19 |    1 |     42.5ms |
| Non-vectorized Pcubature       | 8.31ms |  8.45ms | 118.032 |        0B | 39.344 |    15 |    5 |    127.1ms |
| Vectorized Pcubature           | 10.5ms | 10.68ms |  93.491 |   386.2KB | 23.373 |    16 |    4 |    171.1ms |
| Non-vectorized cubature::cuhre |  4.3ms |  4.42ms | 225.278 |        0B | 39.755 |    17 |    3 |     75.5ms |
| Vectorized cubature::cuhre     | 4.87ms |  4.96ms | 199.766 |   225.8KB | 22.196 |    18 |    2 |     90.1ms |

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
| Non-vectorized Hcubature       |   3.4ms |   3.5ms | 266.871 |   32.89KB | 26.687 |    20 |    2 |     74.9ms |
| Vectorized Hcubature           |  4.11ms |  4.25ms | 221.491 |  205.61KB | 22.149 |    20 |    2 |     90.3ms |
| Non-vectorized Pcubature       |  8.55ms |  8.76ms | 106.642 |        0B | 26.661 |    20 |    5 |    187.5ms |
| Vectorized Pcubature           |   9.8ms | 10.27ms |  90.023 |  386.24KB | 22.506 |    20 |    5 |    222.2ms |
| Non-vectorized cubature::cuhre | 44.62ms | 45.51ms |  21.600 |        0B | 23.760 |    20 |   22 |    925.9ms |
| Vectorized cubature::cuhre     |  45.2ms | 46.49ms |  21.341 |    2.04MB | 22.408 |    20 |   21 |    937.2ms |

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
| Non-vectorized Hcubature       | 129.8µs |   136µs |  7173.363 |    53.7KB |  0.000 |   100 |    0 |    13.94ms |
| Vectorized Hcubature           | 215.3µs | 226.8µs |  4109.128 |    68.5KB | 41.506 |    99 |    1 |    24.09ms |
| Non-vectorized Pcubature       |  52.2µs |  55.2µs | 17406.580 |        0B |  0.000 |   100 |    0 |     5.75ms |
| Vectorized Pcubature           | 159.3µs | 166.4µs |  5687.127 |        0B |  0.000 |   100 |    0 |    17.58ms |
| Non-vectorized cubature::cuhre | 222.2µs |   237µs |  4039.016 |        0B | 40.798 |    99 |    1 |    24.51ms |
| Vectorized cubature::cuhre     | 508.3µs | 544.3µs |  1795.310 |        0B | 18.134 |    99 |    1 |    55.14ms |

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
| Non-vectorized Hcubature       |   4.65ms |   4.75ms |  208.162 |    78.4KB | 33.887 |    86 |   14 |    413.1ms |
| Vectorized Hcubature           |   5.69ms |    5.9ms |  168.298 |   304.7KB | 25.148 |    87 |   13 |    516.9ms |
| Non-vectorized Pcubature       | 417.17µs | 435.76µs | 2244.896 |        0B | 22.676 |    99 |    1 |     44.1ms |
| Vectorized Pcubature           | 628.56µs | 658.09µs | 1488.488 |    18.3KB | 30.377 |    98 |    2 |     65.8ms |
| Non-vectorized cubature::cuhre |   1.25ms |   1.29ms |  766.829 |        0B | 23.716 |    97 |    3 |    126.5ms |
| Vectorized cubature::cuhre     |   1.38ms |   1.43ms |  689.789 |    60.1KB | 21.334 |    97 |    3 |    140.6ms |

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
    ## [1] mvtnorm_1.3-3    cubature_2.1.4-6 bench_1.1.4     
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
