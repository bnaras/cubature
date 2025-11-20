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

| expression                     |     min |   median |  itr/sec | mem_alloc | gc/sec | n_itr | n_gc | total_time |
|:-------------------------------|--------:|---------:|---------:|----------:|-------:|------:|-----:|-----------:|
| Non-vectorized Hcubature       | 273.4µs | 283.98µs | 3476.985 |    59.3KB | 35.121 |    99 |    1 |     28.5ms |
| Vectorized Hcubature           | 392.9µs | 413.93µs | 2387.530 |    54.5KB | 24.116 |    99 |    1 |     41.5ms |
| Non-vectorized Pcubature       | 865.4µs | 889.83µs | 1108.375 |    31.2KB | 34.280 |    97 |    3 |     87.5ms |
| Vectorized Pcubature           |   1.2ms |   1.23ms |  803.568 |    58.4KB | 16.399 |    98 |    2 |      122ms |
| Non-vectorized cubature::cuhre | 589.7µs | 608.27µs | 1630.477 |    40.5KB | 16.469 |    99 |    1 |     60.7ms |
| Vectorized cubature::cuhre     | 599.4µs |  634.2µs | 1573.101 |    39.8KB | 32.104 |    98 |    2 |     62.3ms |

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
| Non-vectorized Hcubature       | 824.16ms | 835.03ms |   1.171 |  139.68KB | 18.152 |    10 |  155 |      8.54s |
| Vectorized Hcubature           |   1.81ms |   2.46ms | 435.490 |    1.89MB |  0.000 |    10 |    0 |    22.96ms |
| Non-vectorized Pcubature       | 362.99ms | 364.93ms |   2.686 |        0B | 17.460 |    10 |   65 |      3.72s |
| Vectorized Pcubature           |   1.27ms |   1.29ms | 652.764 |  810.26KB | 65.276 |    10 |    1 |    15.32ms |
| Non-vectorized cubature::cuhre | 344.81ms | 350.35ms |   2.843 |        0B | 17.912 |    10 |   63 |      3.52s |
| Vectorized cubature::cuhre     |   3.23ms |    3.3ms | 300.433 |  898.41KB |  0.000 |    10 |    0 |    33.28ms |

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
| g1()       | 767µs | 1.37ms | 737.155 |  355.13KB |      0 |    20 |    0 |     27.1ms |
| g2()       | 760µs | 1.37ms | 759.279 |    2.49KB |      0 |    20 |    0 |     26.3ms |
| g3()       | 761µs | 1.38ms | 736.959 |    2.49KB |      0 |    20 |    0 |     27.1ms |

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
| Non-vectorized Hcubature       |  36.9µs |  39.1µs | 24825.795 |    7.66KB | 24.851 |   999 |    1 |     40.2ms |
| Vectorized Hcubature           |  64.6µs |  67.9µs | 14084.502 |    1.16KB | 28.225 |   998 |    2 |     70.9ms |
| Non-vectorized Pcubature       |  48.8µs |  50.7µs | 19088.426 |        0B | 19.108 |   999 |    1 |     52.3ms |
| Vectorized Pcubature           | 113.6µs | 120.1µs |  8119.501 |   18.68KB | 24.432 |   997 |    3 |    122.8ms |
| Non-vectorized cubature::cuhre | 339.8µs | 352.7µs |  2827.734 |        0B | 25.681 |   991 |    9 |    350.5ms |
| Vectorized cubature::cuhre     | 366.3µs |   386µs |  2574.422 |   16.38KB | 26.004 |   990 |   10 |    384.6ms |

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
| Non-vectorized Hcubature       |   2.86ms |  2.92ms |   340.453 |    67.3KB |  0.000 |    10 |    0 |     29.4ms |
| Vectorized Hcubature           |   4.97ms |  5.08ms |   196.921 |   290.5KB | 49.230 |     8 |    2 |     40.6ms |
| Non-vectorized Pcubature       |  74.62µs | 75.54µs | 12330.201 |        0B |  0.000 |    10 |    0 |      811µs |
| Vectorized Pcubature           | 155.25µs | 167.2µs |  5900.646 |     4.1KB |  0.000 |    10 |    0 |      1.7ms |
| Non-vectorized cubature::cuhre |  13.63ms | 13.79ms |    72.467 |        0B | 31.057 |     7 |    3 |     96.6ms |
| Vectorized cubature::cuhre     |  20.68ms | 20.75ms |    48.061 |   971.5KB | 48.061 |     5 |    5 |      104ms |

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
| Non-vectorized Hcubature       |  47.4ms |  48.5ms |  17.574 |    17.7KB | 35.147 |    10 |   20 |   569.03ms |
| Vectorized Hcubature           |  49.9ms |  50.5ms |  19.443 |  1011.8KB | 25.275 |    10 |   13 |   514.33ms |
| Non-vectorized cubature::cuhre | 822.7ms |   839ms |   1.189 |        0B | 27.936 |    10 |  235 |      8.41s |
| Vectorized cubature::cuhre     |   913ms | 933.9ms |   1.063 |    21.2MB | 24.121 |    10 |  227 |      9.41s |

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
| Non-vectorized Hcubature       |  58.6µs |  60.6µs | 15798.461 |    6.75KB |   0.000 |    50 |    0 |     3.17ms |
| Vectorized Hcubature           |  92.5µs |  96.9µs | 10045.412 |    2.91KB |   0.000 |    50 |    0 |     4.98ms |
| Non-vectorized Pcubature       |  52.8µs |  55.1µs | 17102.217 |        0B | 349.025 |    49 |    1 |     2.87ms |
| Vectorized Pcubature           |  81.4µs |  83.9µs | 11421.537 |   19.12KB |   0.000 |    50 |    0 |     4.38ms |
| Non-vectorized cubature::cuhre | 648.5µs | 688.9µs |  1225.759 |        0B |  25.015 |    49 |    1 |    39.98ms |
| Vectorized cubature::cuhre     | 660.4µs | 692.9µs |  1437.311 |   39.84KB |   0.000 |    50 |    0 |    34.79ms |

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
| Non-vectorized Hcubature       | 1.37ms |  1.4ms | 708.533 |    85.2KB | 37.291 |    19 |    1 |     26.8ms |
| Vectorized Hcubature           |  1.8ms | 1.84ms | 539.464 |   147.2KB | 28.393 |    19 |    1 |     35.2ms |
| Non-vectorized Pcubature       | 2.06ms | 2.12ms | 473.174 |        0B | 24.904 |    19 |    1 |     40.2ms |
| Vectorized Pcubature           | 2.59ms | 2.67ms | 372.258 |    68.9KB | 41.362 |    18 |    2 |     48.4ms |
| Non-vectorized cubature::cuhre | 3.26ms | 3.36ms | 297.002 |        0B | 33.000 |    18 |    2 |     60.6ms |
| Vectorized cubature::cuhre     | 3.78ms | 3.87ms | 258.543 |   125.6KB | 13.608 |    19 |    1 |     73.5ms |

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
| Non-vectorized Hcubature       | 3.58ms | 3.65ms | 272.232 |     133KB | 48.041 |    17 |    3 |     62.4ms |
| Vectorized Hcubature           |  4.7ms |  4.8ms | 207.708 |     249KB | 23.079 |    18 |    2 |     86.7ms |
| Non-vectorized Pcubature       | 2.51ms | 2.57ms | 389.034 |        0B | 43.226 |    18 |    2 |     46.3ms |
| Vectorized Pcubature           | 3.27ms | 3.37ms | 296.798 |      69KB | 15.621 |    19 |    1 |       64ms |
| Non-vectorized cubature::cuhre | 6.98ms | 7.07ms | 141.129 |        0B | 47.043 |    15 |    5 |    106.3ms |
| Vectorized cubature::cuhre     | 8.21ms |  8.3ms | 120.195 |     224KB | 21.211 |    17 |    3 |    141.4ms |

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
| Non-vectorized Hcubature       |  1.62ms |  1.65ms | 601.654 |    64.6KB | 31.666 |    19 |    1 |     31.6ms |
| Vectorized Hcubature           |  2.17ms |  2.23ms | 446.607 |     156KB | 23.506 |    19 |    1 |     42.5ms |
| Non-vectorized Pcubature       |  8.34ms |  8.51ms | 117.186 |        0B | 39.062 |    15 |    5 |      128ms |
| Vectorized Pcubature           | 10.39ms | 10.64ms |  93.990 |   386.2KB | 31.330 |    15 |    5 |    159.6ms |
| Non-vectorized cubature::cuhre |  4.36ms |  4.47ms | 223.863 |        0B | 24.874 |    18 |    2 |     80.4ms |
| Vectorized cubature::cuhre     |  4.84ms |  4.99ms | 200.041 |   225.8KB | 22.227 |    18 |    2 |       90ms |

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
| Non-vectorized Hcubature       |  3.44ms |   3.6ms | 257.498 |   32.89KB | 25.750 |    20 |    2 |     77.7ms |
| Vectorized Hcubature           |  4.11ms |  4.27ms | 219.427 |  205.61KB | 21.943 |    20 |    2 |     91.1ms |
| Non-vectorized Pcubature       |  8.51ms |  8.76ms | 106.016 |        0B | 26.504 |    20 |    5 |    188.7ms |
| Vectorized Pcubature           |  9.73ms | 10.11ms |  93.695 |  386.24KB | 23.424 |    20 |    5 |    213.5ms |
| Non-vectorized cubature::cuhre | 44.05ms | 45.88ms |  21.731 |        0B | 23.904 |    20 |   22 |    920.3ms |
| Vectorized cubature::cuhre     | 45.23ms | 46.06ms |  21.669 |    2.04MB | 22.752 |    20 |   21 |      923ms |

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
| Non-vectorized Hcubature       | 124.6µs | 128.4µs |  7578.763 |    53.7KB |  0.000 |   100 |    0 |     13.2ms |
| Vectorized Hcubature           | 210.1µs | 225.2µs |  4255.681 |    68.5KB | 42.987 |    99 |    1 |    23.26ms |
| Non-vectorized Pcubature       |  50.9µs |  52.7µs | 18466.789 |        0B |  0.000 |   100 |    0 |     5.42ms |
| Vectorized Pcubature           | 155.8µs | 162.9µs |  5928.579 |        0B |  0.000 |   100 |    0 |    16.87ms |
| Non-vectorized cubature::cuhre | 221.2µs | 233.8µs |  4156.490 |        0B | 41.985 |    99 |    1 |    23.82ms |
| Vectorized cubature::cuhre     | 505.3µs |   542µs |  1812.692 |        0B | 18.310 |    99 |    1 |    54.62ms |

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
| Non-vectorized Hcubature       |   4.54ms |   4.63ms |  214.133 |    78.4KB | 34.859 |    86 |   14 |    401.6ms |
| Vectorized Hcubature           |   5.48ms |   5.62ms |  176.553 |   304.7KB | 26.381 |    87 |   13 |    492.8ms |
| Non-vectorized Pcubature       | 411.52µs | 433.27µs | 2238.560 |        0B | 22.612 |    99 |    1 |     44.2ms |
| Vectorized Pcubature           | 612.62µs | 642.55µs | 1527.907 |    18.3KB | 15.433 |    99 |    1 |     64.8ms |
| Non-vectorized cubature::cuhre |   1.25ms |   1.28ms |  770.809 |        0B | 32.117 |    96 |    4 |    124.5ms |
| Vectorized cubature::cuhre     |   1.34ms |    1.4ms |  703.928 |    60.1KB | 21.771 |    97 |    3 |    137.8ms |

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
